#include "TimeStepPF.h"

#include "SimulationDataPF.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/Simulation.h"

#include <atomic>
#include <iostream>

using namespace SPH;
using namespace std;
using namespace GenParam;

int TimeStepPF::STIFFNESS = -1;

#define Vec3Block(i) template segment<3>(3 * (i))

// helper functions
namespace
{
	template <typename T>
	struct atomic_wrapper
	{
		std::atomic<T> _a;

		atomic_wrapper() : _a(0) {}
		atomic_wrapper(const std::atomic<T> &a) :_a(a.load()) {}
		atomic_wrapper(const atomic_wrapper &other) :_a(other._a.load()) {}
		atomic_wrapper &operator=(const atomic_wrapper &other)
		{
			_a.store(other._a.load());
			return *this;
		}
	};

	inline void addToAtomicReal(std::atomic<Real> & a, const Real & r)
	{
		Real current = a;
		while (!a.compare_exchange_weak(current, current + r))
			;
	}
}

using AtomicRealVec = std::vector < atomic_wrapper<Real> >;

TimeStepPF::TimeStepPF() :
	TimeStep()
{
	m_simulationData.init();
	m_counter = 0;
	m_stiffness = 50000.0;
}

TimeStepPF::~TimeStepPF(void)
{
}

void TimeStepPF::initParameters()
{
	TimeStep::initParameters();

	STIFFNESS = createNumericParameter("stiffness", "Stiffness", &m_stiffness);
	setGroup(STIFFNESS, "PF");
	setDescription(STIFFNESS, "Stiffness coefficient.");
	static_cast<RealParameter*>(getParameter(STIFFNESS))->setMinValue(1e-6);
}

void TimeStepPF::step()
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getModel();
	TimeManager *tm = TimeManager::getCurrent ();

	clearAccelerations();
	initialGuessForPositions();
	performNeighborhoodSearch();

	START_TIMING("solvePDConstraints");
	solvePDConstraints();
	STOP_TIMING_AVG;

	computeDensities();
	sim->computeNonPressureForces();
	addAccellerationToVelocity();

	// Compute new time	
	sim->updateTimeStepSize();
	sim->emitParticles();
	tm->setTime(tm->getTime () + tm->getTimeStepSize());
}

void TimeStepPF::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
}

void TimeStepPF::initialGuessForPositions()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const auto numParticles = model->numActiveParticles();
	const auto h = TimeManager::getCurrent()->getTimeStepSize();

	#pragma omp parallel for
	for (int i = 0; i < (int)numParticles; i++)
	{
		m_simulationData.setOldPosition(i, model->getPosition(0, i));
		const auto newPos = (model->getPosition(0, i) + h * model->getVelocity(0, i) + (h * h) * model->getAcceleration(i)).eval();
		model->setPosition(0, i, newPos);
		m_simulationData.setS(i, newPos);
	}
}

void TimeStepPF::solvePDConstraints()
{
	const auto model = Simulation::getCurrent()->getModel();
	const auto numParticles = model->numActiveParticles();

	//////////////////////////////////////////////////////////////////////////
	// initialize positions
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel for schedule(static)  
	for (int i = 0; i < (int)numParticles; i++)
	{
		m_simulationData.getX()[i] = m_simulationData.getS(i);
	}

	//////////////////////////////////////////////////////////////////////////
	// count number of fluid neighbors for relaxation
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel for schedule(static)  
	for (int i = 0; i < (int)numParticles; i++)
	{
		m_simulationData.setNumFluidNeighbors(i, model->numberOfNeighbors(0, i) + 1u);
	}

	//////////////////////////////////////////////////////////////////////////
	// Init linear system solver and preconditioner
	//////////////////////////////////////////////////////////////////////////
	MatrixReplacement A(3 * numParticles, matrixVecProd, (void*) this);
#ifdef PD_USE_DIAGONAL_PRECONDITIONER
	preparePreconditioner();
	m_solver.preconditioner().init(numParticles, diagonalMatrixElement, (void*)this);
#endif
	m_solver.setMaxIterations(m_maxIterations);
	m_solver.compute(A);

	auto x = Eigen::Map<VectorXr>(m_simulationData.getX().data()->data(), 3 * numParticles,	1);
	VectorXr b(3 * numParticles);
	for (m_iterations = 0u; m_iterations < m_maxIterations; m_iterations++)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute RHS
		//////////////////////////////////////////////////////////////////////////
		matrixFreeRHS(b);

		//////////////////////////////////////////////////////////////////////////
		// Solve linear system 
		//////////////////////////////////////////////////////////////////////////
#ifdef PD_USE_DIAGONAL_PRECONDITIONER
		// hack to make the solver perform at least min_iter CG iterations for stability
		constexpr unsigned int min_iter = 1;
		m_solver.setTolerance(m_iterations < min_iter ? 1e-32 : 1e-10);
#endif
		x = m_solver.solveWithGuess(b, x);
		if (m_solver.iterations() == 0)
			break;
	}
	updatePositionsAndVelocity();
}

void TimeStepPF::updatePositionsAndVelocity()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const auto  numParticles = model->numActiveParticles();
	const auto  h            = TimeManager::getCurrent()->getTimeStepSize();
	const auto& x            = m_simulationData.getX();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			auto vel = (x[i] - m_simulationData.getOldPosition(i)) / h;
			model->setPosition(0, i, x[i]);
			model->setVelocity(0, i, vel);
		}
	}
}

#ifdef PD_USE_DIAGONAL_PRECONDITIONER
FORCE_INLINE void SPH::TimeStepPF::diagonalMatrixElement(const unsigned int row, Vector3r &result, void *userData)
{
	TimeStepPF *timeStepPF = static_cast<TimeStepPF*>(userData);
	result = timeStepPF->m_simulationData.getDiag(row);
}

void SPH::TimeStepPF::preparePreconditioner()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const auto h = TimeManager::getCurrent()->getTimeStepSize();
	const auto system_scale = h * h * m_stiffness;
	const auto numParticles = model->numActiveParticles();
	#pragma omp parallel for schedule(static)  
	for (int i = 0; i < (int)numParticles; i++)
	{
		const auto numFluidNeighbors = model->numberOfNeighbors(0, i);
		m_simulationData.setDiag(i, Vector3r::Constant(system_scale * (numFluidNeighbors + 1) + model->getMass(i)));
	}
}

#endif

void SPH::TimeStepPF::addAccellerationToVelocity()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	#pragma omp parallel default(shared)
	{
		const auto  numParticles = model->numActiveParticles();
		const auto  h = TimeManager::getCurrent()->getTimeStepSize();
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			model->setVelocity(0, i, model->getVelocity(0, i) + h * model->getAcceleration(i));
		}
	}
}

/** \brief compute the right hand side of the system in a matrix-free fashion and store the result in result*/
void SPH::TimeStepPF::matrixFreeRHS(VectorXr & result)
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const auto numParticles = model->numActiveParticles();
	const auto h            = TimeManager::getCurrent()->getTimeStepSize();
	const Real density0 = model->getValue<Real>(FluidModel::DENSITY0);

	AtomicRealVec accumulator(3 * numParticles);

	// influence of pressure
	#pragma omp parallel default(shared)
	{
		const auto & x = m_simulationData.getX();
		const auto density0_inv = 1 / density0;

		// local step for fluid constraints
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			unsigned int numNeighbors = 0;
			for (unsigned int pid = 0; pid < model->numberOfPointSets(); pid++)
			{
				numNeighbors += model->numberOfNeighbors(pid, i);
			}
			//////////////////////////////////////////////////////////////////////////
			// particle positions in current constraint, will be projected
			//////////////////////////////////////////////////////////////////////////
			std::vector<Vector3r> p(numNeighbors + 1);
			// the i'th particle itself
			p[0] = x[i];
			unsigned int counter = 1;
			// fluid neighbors
			for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = model->getNeighbor(0, i, j);
				p[counter++] = x[neighborIndex];
			}
			// boundary neighbors
			for (unsigned int pid = 1; pid < model->numberOfPointSets(); pid++)
			{
				for (unsigned int j = 0; j < model->numberOfNeighbors(pid, i); j++)
				{
					const unsigned int neighborIndex = model->getNeighbor(pid, i, j);
					p[counter++] = model->getPosition(pid, neighborIndex);
				}
			}
			//////////////////////////////////////////////////////////////////////////
			// helper functions
			//////////////////////////////////////////////////////////////////////////
			// constraint value
			auto calculateC = [&]() -> Real
			{
				// Compute current density for particle i
				Real density = model->getMass(i) * model->W_zero();
				const Vector3r &xi = p[0];
				unsigned int counter = 1;
				for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = model->getNeighbor(0, i, j);
					const auto& xj = p[counter++];
					density += model->getMass(neighborIndex) * model->W(xi - xj);
				}
				for (unsigned int pid = 1; pid < model->numberOfPointSets(); pid++)
				{
					for (unsigned int j = 0; j < model->numberOfNeighbors(pid, i); j++)
					{
						const unsigned int neighborIndex = model->getNeighbor(pid, i, j);
						const auto& xj = p[counter++];
						// Boundary: Akinci2012
						density += model->getBoundaryPsi(pid, neighborIndex) * model->W(xi - xj);
					}
				}
				// constraint value = density / density0 - 1
				const auto C = density * density0_inv - 1;
				// pressure clamping
				return (C < 0) ? 0 : C;
			};
			// constraint gradient
			auto calculateNablaC = [&]() -> std::vector<Vector3r>
			{
				std::vector<Vector3r> nablaC(p.size());
				nablaC[0].setZero();
				const Vector3r &xi = p[0];

				unsigned int counter = 1;
				for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = model->getNeighbor(0, i, j);
					const auto& xj = p[counter];
					nablaC[counter] = (-density0_inv * model->getMass(neighborIndex)) * model->gradW(xi - xj);
					nablaC[0] -= nablaC[counter];
					counter++;
				}
				for (unsigned int pid = 1; pid < model->numberOfPointSets(); pid++)
				{
					for (unsigned int j = 0; j < model->numberOfNeighbors(pid, i); j++)
					{
						const unsigned int neighborIndex = model->getNeighbor(pid, i, j);
						const auto& xj = p[counter];
						// Boundary: Akinci2012
						nablaC[counter] = (-density0_inv * model->getBoundaryPsi(pid, neighborIndex)) * model->gradW(xi - xj);
						nablaC[0] -= nablaC[counter];
						counter++;
					}
				}
				return nablaC;
			};
			//////////////////////////////////////////////////////////////////////////
			// constraint projection
			//////////////////////////////////////////////////////////////////////////
			const auto C_goal    = Real(1e-14);
			const auto max_steps = 100u;
			      auto it        = 0u;
			      auto C	     = calculateC();
			while ((std::abs(C) > C_goal) && it++ < max_steps)
			{
				const auto g = calculateNablaC();
				const auto dg = [&g]() { Real s = 0;  for (const auto & x : g) { s += x.squaredNorm(); } return s; }();
				if (dg == 0) break;	// found a minimum
				const auto cdg = -C / (dg + 1e-6f); // add regularization factor

				// move fluid particles along constraint gradient
				p[0] += (cdg * m_simulationData.getNumFluidNeighbors(i)) * g[0];

				// only fluid particles are projected
				for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
				{
					const auto neighborIndex = model->getNeighbor(0, i, j);
					const auto nfn = m_simulationData.getNumFluidNeighbors(neighborIndex);
					p[j + 1] += (cdg * nfn) * g[j + 1];
				}

				if (it + 1 < max_steps)
				{
					C = calculateC();
				}
			}
			// update RHS
			for (auto c = 0u; c < 3; c++)
				addToAtomicReal(accumulator[3 * i + c]._a, p[0][c]);
			for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
			{
				const auto neighborIndex = model->getNeighbor(0, i, j);
				for (auto c = 0u; c < 3; c++)
					addToAtomicReal(accumulator[3 * neighborIndex + c]._a, p[j+1][c]);
			}
		}
	}

	// influence of momentum
	#pragma omp parallel default(shared)
	{
		const auto system_scale = h * h * m_stiffness;
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const auto  m = model->getMass(i);
			const auto& s = m_simulationData.getS(i);
			for (auto c = 0u; c < 3; c++)
			{
				const auto id = 3 * i + c;
				result[id] = system_scale * accumulator[id]._a + m * s[c];
			}
		}
	}
}

void SPH::TimeStepPF::matrixVecProd(const Real* vec, Real *result, void *userData)
{
	TimeStepPF *timeStepPF = static_cast<TimeStepPF*>(userData);

	FluidModel *model = Simulation::getCurrent()->getModel();
	const auto numParticles = model->numActiveParticles();
	const auto h = TimeManager::getCurrent()->getTimeStepSize();

	AtomicRealVec accumulator(3 * numParticles);

	// influence of pressure
	#pragma omp parallel for schedule(static)  
	for (int i = 0; i < (int)numParticles; i++)
	{
		const auto numNeighbors = model->numberOfNeighbors(0, i);
		for (auto c = 0u; c < 3; c++)
			addToAtomicReal(accumulator[3 * i + c]._a, vec[3 * i + c]);
		for (auto j = 0u; j < numNeighbors; j++)
		{
			const auto neighborIndex = model->getNeighbor(0, i, j);
			for (auto c = 0u; c < 3; c++)
				addToAtomicReal(accumulator[3 * neighborIndex + c]._a, vec[3 * neighborIndex + c]);
		}
	}

	// influence of momentum
	const auto system_scale = h * h * timeStepPF->m_stiffness;
	#pragma omp parallel for schedule(static)  
	for (int i = 0; i < (int)numParticles; i++)
		for (auto c = 0u; c < 3; c++)
			result[3 * i + c] = system_scale * accumulator[3 * i + c]._a + model->getMass(i) * vec[3 * i + c];
}

void TimeStepPF::performNeighborhoodSearch()
{
	if (m_counter % 500 == 0)
	{
		FluidModel *model = Simulation::getCurrent()->getModel();
		model->performNeighborhoodSearchSort();
		m_simulationData.performNeighborhoodSearchSort();
		Simulation::getCurrent()->performNeighborhoodSearchSort();
	}
	m_counter++;

	Simulation::getCurrent()->performNeighborhoodSearch();
}

void TimeStepPF::emittedParticles(const unsigned int startIndex)
{
	m_simulationData.emittedParticles(startIndex);
	Simulation::getCurrent()->emittedParticles(startIndex);
}

void TimeStepPF::resize()
{
	m_simulationData.init();
}
