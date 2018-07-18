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
	TimeStep(),
	m_stiffness(50000.0),
	m_counter(0),
	m_numActiveParticlesTotal(0)
{
	m_simulationData.init();
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
	const unsigned int nModels = sim->numberOfFluidModels();
	TimeManager *tm = TimeManager::getCurrent ();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		clearAccelerations(fluidModelIndex);
		initialGuessForPositions(fluidModelIndex);
	}
	performNeighborhoodSearch();

	START_TIMING("solvePDConstraints");
	solvePDConstraints();
	STOP_TIMING_AVG;

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		computeDensities(fluidModelIndex);
	sim->computeNonPressureForces();
	addAccellerationToVelocity();

	// update emitters
	sim->emitParticles();
	// Compute new time	
	sim->updateTimeStepSize();
	tm->setTime(tm->getTime () + tm->getTimeStepSize());
}

void TimeStepPF::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
}

void TimeStepPF::initialGuessForPositions(const unsigned int fluidModelIndex)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const auto numParticles = model->numActiveParticles();
	const auto h = TimeManager::getCurrent()->getTimeStepSize();

	#pragma omp parallel for
	for (int i = 0; i < (int)numParticles; i++)
	{
		m_simulationData.setOldPosition(fluidModelIndex, i, model->getPosition(i));
		const auto newPos = (model->getPosition(i) + h * model->getVelocity(i) + (h * h) * model->getAcceleration(i)).eval();
		model->setPosition(i, newPos);
		m_simulationData.setS(fluidModelIndex, i, newPos);
	}
}

void TimeStepPF::solvePDConstraints()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	if (nModels == 0)
		return;

	// total number of active fluid particles
	m_numActiveParticlesTotal = m_simulationData.getParticleOffset(nModels - 1) + sim->getFluidModel(nModels - 1)->numActiveParticles();
	
	VectorXr x(3 * m_numActiveParticlesTotal);
	VectorXr b(3 * m_numActiveParticlesTotal);

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const unsigned int offset = m_simulationData.getParticleOffset(fluidModelIndex);
		
		#pragma omp parallel for schedule(static)  
		for (int i = 0; i < (int)model->numActiveParticles(); i++)
		{
			//////////////////////////////////////////////////////////////////////////
			// initialize positions
			//////////////////////////////////////////////////////////////////////////
			x.Vec3Block(offset + i) = m_simulationData.getS(fluidModelIndex, i);

			//////////////////////////////////////////////////////////////////////////
			// count number of fluid neighbors for relaxation
			//////////////////////////////////////////////////////////////////////////
			unsigned int nNeighbors = 0;
			for (unsigned int pid = 0; pid < nModels; pid++)
				nNeighbors += sim->numberOfNeighbors(fluidModelIndex, pid, i);
			m_simulationData.setNumFluidNeighbors(fluidModelIndex, i, nNeighbors + 1u);
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Init linear system solver and preconditioner
	//////////////////////////////////////////////////////////////////////////
	MatrixReplacement A(3 * m_numActiveParticlesTotal, matrixVecProd, (void*) this);
#ifdef PD_USE_DIAGONAL_PRECONDITIONER
	preparePreconditioner();
	m_solver.preconditioner().init(m_numActiveParticlesTotal, diagonalMatrixElement, (void*)this);
#endif
	m_solver.setMaxIterations(m_maxIterations);
	m_solver.compute(A);
	
	for (m_iterations = 0u; m_iterations < m_maxIterations; m_iterations++)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute RHS
		//////////////////////////////////////////////////////////////////////////
		matrixFreeRHS(x,b);

		//////////////////////////////////////////////////////////////////////////
		// Solve linear system 
		//////////////////////////////////////////////////////////////////////////
#ifdef PD_USE_DIAGONAL_PRECONDITIONER
		// hack to make the solver perform at least min_iter CG iterations for stability
		constexpr unsigned int min_iter = 1;
		m_solver.setTolerance(m_iterations < min_iter ? static_cast<Real>(1e-32) : static_cast<Real>(1e-10));
#endif
		x = m_solver.solveWithGuess(b, x);
		if (m_solver.iterations() == 0)
			break;
	}

	updatePositionsAndVelocity(x);
}

void SPH::TimeStepPF::updatePositionsAndVelocity(const VectorXr & x)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	const Real h_inv = 1 / TimeManager::getCurrent()->getTimeStepSize();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();
		const unsigned int offset = m_simulationData.getParticleOffset(fluidModelIndex);
		
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				const Vector3r vel = h_inv * (x.Vec3Block(offset + i) - m_simulationData.getOldPosition(fluidModelIndex, i));
				model->setPosition(i, x.Vec3Block(offset + i));
				model->setVelocity(i, vel);
			}
		}
	}
}

#ifdef PD_USE_DIAGONAL_PRECONDITIONER
FORCE_INLINE void SPH::TimeStepPF::diagonalMatrixElement(const unsigned int row, Vector3r &result, void *userData)
{
	TimeStepPF *timeStepPF = static_cast<TimeStepPF*>(userData);
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	// TODO: use m_simulationData->getParticleOffset instead of accessing FluidModels
	// find corresponding fluid model and particle index for the current row
	unsigned int fluidRow = row;
	unsigned int fluidModelIndex;
	for (fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);		
		if (fluidRow >= model->numActiveParticles())
			fluidRow -= model->numActiveParticles();
		else
			break;
	}

	result = timeStepPF->m_simulationData.getDiag(fluidModelIndex, fluidRow);
}

void SPH::TimeStepPF::preparePreconditioner()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	const auto h = TimeManager::getCurrent()->getTimeStepSize();
	const auto system_scale = h * h * m_stiffness;

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const auto numParticles = model->numActiveParticles();
		const Real density0_inv = 1 / model->getDensity0();

		#pragma omp parallel for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Real density_scale = 1;
			for (unsigned int pid = 0; pid < nModels; pid++)
			{
				density_scale += density0_inv * sim->getFluidModelFromPointSet(pid)->getDensity0() * sim->numberOfNeighbors(fluidModelIndex, pid, i);
			}
			m_simulationData.setDiag(fluidModelIndex, i, Vector3r::Constant(system_scale * density_scale + model->getMass(i)));
		}
	}
}

#endif

void SPH::TimeStepPF::addAccellerationToVelocity()
{
	Simulation *sim = Simulation::getCurrent();
	const auto  h = TimeManager::getCurrent()->getTimeStepSize();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const auto  numParticles = model->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				model->setVelocity(i, model->getVelocity(i) + h * model->getAcceleration(i));
			}
		}
	}
}

/** \brief compute the right hand side of the system in a matrix-free fashion and store the result in result*/
void SPH::TimeStepPF::matrixFreeRHS(const VectorXr & x, VectorXr & result)
{
	Simulation *sim = Simulation::getCurrent();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();;
	const unsigned int nModels = sim->numberOfFluidModels();

	AtomicRealVec accumulator(3 * m_numActiveParticlesTotal);

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();
		const unsigned int offset = m_simulationData.getParticleOffset(fluidModelIndex);
		const Real density0 = model->getDensity0();

		// influence of pressure
#pragma omp parallel default(shared)
		{
			// local step for fluid constraints
#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// total number of neighbors in all point sets (fluids and boundaries)
				//////////////////////////////////////////////////////////////////////////
				unsigned int numParticlesInConstraint = 1;
				for (unsigned int pid = 0; pid < sim->numberOfPointSets(); pid++)
				{
					numParticlesInConstraint += sim->numberOfNeighbors(fluidModelIndex, pid, i);
				}
				//////////////////////////////////////////////////////////////////////////
				// particle positions in current constraint, will be projected
				//////////////////////////////////////////////////////////////////////////
				std::vector<Vector3r> p(numParticlesInConstraint);
				// the i'th particle itself
				p[0] = x.Vec3Block(offset + i);
				unsigned int counter = 1;
				// fluid neighbors
				for (unsigned int pid = 0; pid < nModels; pid++)
				{
					const unsigned int neighborOffset = m_simulationData.getParticleOffset(pid);
					for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
					{
						const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
						p[counter++] = x.Vec3Block(neighborOffset + neighborIndex);
					}
				}
				// boundary neighbors
				for (unsigned int pid = nModels; pid < sim->numberOfPointSets(); pid++)
				{
					const BoundaryModel *bm_neighbor = sim->getBoundaryModelFromPointSet(pid);
					for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
					{
						const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
						p[counter++] = bm_neighbor->getPosition(neighborIndex);
					}
				}
				//////////////////////////////////////////////////////////////////////////
				// helper functions
				//////////////////////////////////////////////////////////////////////////
				// constraint value
				auto calculateC = [&]() -> Real
				{
					// Compute current density for particle i
					Real density = model->getVolume(i) * sim->W_zero();
					const Vector3r &xi = p[0];
					unsigned int counter = 1;
					for (unsigned int pid = 0; pid < nModels; pid++)
					{
						const FluidModel *fm_neighbor = sim->getFluidModelFromPointSet(pid);
						for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
						{
							const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
							const Vector3r & xj = p[counter++];
							density += fm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
						}
					}
					for (unsigned int pid = nModels; pid < sim->numberOfPointSets(); pid++)
					{
						const BoundaryModel *bm_neighbor = sim->getBoundaryModelFromPointSet(pid);
						for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
						{
							const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
							const Vector3r & xj = p[counter++];
							// Boundary: Akinci2012
							density += bm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
						}
					}
					// constraint value = density / density0 - 1
					const auto C = density - 1;
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
					for (unsigned int pid = 0; pid < nModels; pid++)
					{
						const FluidModel *fm_neighbor = sim->getFluidModelFromPointSet(pid);
						for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
						{
							const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
							const Vector3r & xj = p[counter];
							nablaC[counter] = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
							nablaC[0] -= nablaC[counter];
							counter++;
						}
					}
					for (unsigned int pid = nModels; pid < sim->numberOfPointSets(); pid++)
					{
						const BoundaryModel *bm_neighbor = sim->getBoundaryModelFromPointSet(pid);
						for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
						{
							const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
							const Vector3r & xj = p[counter];
							// Boundary: Akinci2012
							nablaC[counter] = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
							nablaC[0] -= nablaC[counter];
							counter++;
						}
					}
					return nablaC;
				};
				//////////////////////////////////////////////////////////////////////////
				// constraint projection
				//////////////////////////////////////////////////////////////////////////
				const auto C_goal = Real(1e-14);
				const auto max_steps = 100u;
				auto it = 0u;
				auto C = calculateC();
				while ((std::abs(C) > C_goal) && it++ < max_steps)
				{
					const auto g = calculateNablaC();
					const Real dg = [&g]() { Real s = 0;  for (const auto & x : g) { s += x.squaredNorm(); } return s; }();
					if (dg == 0) break;	// found a minimum
					const Real cdg = -C / (dg + 1e-6f); // add regularization factor

					// move fluid particles along constraint gradient
					p[0] += (cdg * m_simulationData.getNumFluidNeighbors(fluidModelIndex, i)) * g[0];

					// only fluid particles are projected
					unsigned int counter = 1;
					for (unsigned int pid = 0; pid < nModels; pid++)
					{
						for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
						{
							const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
							const unsigned int nfn = m_simulationData.getNumFluidNeighbors(pid, neighborIndex);
							p[counter] += (cdg * nfn) * g[counter];
							counter++;
						}
					}

					if (it + 1 < max_steps)
					{
						C = calculateC();
					}
				}
				// update RHS
				const unsigned int index_i = offset + i;
				for (auto c = 0u; c < 3; c++)
					addToAtomicReal(accumulator[3 * index_i + c]._a,  density0 * p[0][c]);
				counter = 1;
				for (unsigned int pid = 0; pid < nModels; pid++)
				{
					const unsigned int neighborOffset = m_simulationData.getParticleOffset(pid);
					for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
					{
						const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
						const unsigned int index_j = neighborOffset + neighborIndex;
						for (auto c = 0u; c < 3; c++)
							addToAtomicReal(accumulator[3 * index_j + c]._a, density0 * p[counter][c]);
						counter++;
					}
				}
			}
		}
	}

	const Real system_scale = h * h * m_stiffness;
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();
		const unsigned int offset = m_simulationData.getParticleOffset(fluidModelIndex);
		const Real density_scaled_system = system_scale / model->getDensity0();

		// influence of momentum
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				const auto  m = model->getMass(i);
				const auto& s = m_simulationData.getS(fluidModelIndex, i);
				for (auto c = 0u; c < 3; c++)
				{
					const auto id = 3 * (offset + i) + c;
					result[id] = density_scaled_system * accumulator[id]._a + m * s[c];
				}
			}
		}
	}
}

void SPH::TimeStepPF::matrixVecProd(const Real* vec, Real *result, void *userData)
{
	TimeStepPF *timeStepPF = static_cast<TimeStepPF*>(userData);
	
	Simulation *sim = Simulation::getCurrent();
	const Real  h = TimeManager::getCurrent()->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();

	AtomicRealVec accumulator(3 * timeStepPF->m_numActiveParticlesTotal);
	const Real system_scale = h * h * timeStepPF->m_stiffness;

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();
		const unsigned int offset = timeStepPF->m_simulationData.getParticleOffset(fluidModelIndex);
		const Real density0 = model->getDensity0();
		// influence of pressure
		#pragma omp parallel for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			// current particle
			const unsigned int index_i = i + offset;
			for (auto c = 0u; c < 3; c++)
				addToAtomicReal(accumulator[3 * index_i + c]._a, density0 * vec[3 * index_i + c]);
			// fluid neighbors in all fluid models
			for (unsigned int pid = 0; pid < nFluids; pid++)
			{
				const unsigned int neighborOffset = timeStepPF->m_simulationData.getParticleOffset(pid);
				for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
				{
					const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
					const unsigned int index_j = neighborIndex + neighborOffset;
					for (auto c = 0u; c < 3; c++)
						addToAtomicReal(accumulator[3 * index_j + c]._a, density0 * vec[3 * index_j + c]);
				}
			}
		}
	}

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();
		const unsigned int offset = timeStepPF->m_simulationData.getParticleOffset(fluidModelIndex);
		const Real density_scaled_system = system_scale / model->getDensity0();

		// influence of momentum
		#pragma omp parallel for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int index_i = i + offset;
			for (auto c = 0u; c < 3; c++)
			{
				const unsigned int idx = 3 * index_i + c;
				result[idx] = density_scaled_system * accumulator[idx]._a + model->getMass(i) * vec[idx];
			}
		}
	}
}

void TimeStepPF::performNeighborhoodSearch()
{
	if (m_counter % 500 == 0)
	{
		Simulation::getCurrent()->performNeighborhoodSearchSort();
		m_simulationData.performNeighborhoodSearchSort();
	}
	m_counter++;

	Simulation::getCurrent()->performNeighborhoodSearch();
}

void TimeStepPF::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepPF::resize()
{
	m_simulationData.init();
}
