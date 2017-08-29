#include "TimeStepPF.h"

#include "SimulationDataPF.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/Utilities/Timing.h"

#include <atomic>
#include <iostream>

using namespace SPH;
using namespace std;

#define Vec3Block(i) template block<3, 1> (3 * (i), 0)

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

TimeStepPF::TimeStepPF(FluidModel *model) :
	TimeStep(model)
{
	m_simulationData.init(model);
	m_counter = 0;
	m_iterationsV = 0;
	m_stiffness = 50000.0;
}

TimeStepPF::~TimeStepPF(void)
{
}

void TimeStepPF::step()
{
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	const unsigned int numParticles = m_model->numActiveParticles();

	clearAccelerations();
	initialGuessForPositions();
	performNeighborhoodSearch();

	START_TIMING("solvePDConstraints");
	solvePDConstraints();
	STOP_TIMING_AVG;

	computeDensities();
	computeNonPressureForces();
	addAccellerationToVelocity();

	updateTimeStepSize();

	// Compute new time	
	tm->setTime (tm->getTime () + h);
}

void TimeStepPF::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
	m_iterationsV = 0;
}

void TimeStepPF::initialGuessForPositions()
{
	const auto numParticles = m_model->numActiveParticles();
	const auto h = TimeManager::getCurrent()->getTimeStepSize();

#pragma omp parallel for
	for (int i = 0; i < (int)numParticles; i++)
	{
		m_simulationData.setOldPosition(i, m_model->getPosition(0, i));
		const auto newPos = (m_model->getPosition(0, i) + h * m_model->getVelocity(0, i) + (h * h) * m_model->getAcceleration(i)).eval();
		m_model->setPosition(0, i, newPos);
		m_simulationData.setS(i, newPos);
	}
}

void TimeStepPF::prepareSolve()
{
	const auto numParticles = m_model->numActiveParticles();
	auto      x             = VectorXrMap(m_simulationData.getX().data(), m_simulationData.getX().size(), 1);

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			x.Vec3Block(i) = m_model->getPosition(0, i);
			auto nfn = m_model->numberOfNeighbors(0, i) + 1u;

			m_simulationData.setNumFluidNeighbors(i, nfn);
		}
	}
}

void TimeStepPF::solvePDConstraints()
{
	const auto numParticles = m_model->numActiveParticles();

	prepareSolve();

	for (m_iterations = 0u; m_iterations < m_maxIterations; m_iterations++)
	{
		const auto s = cgSolve();
		if (s == CGSolveState::ALREADY_SOLVED) break;
	}

	updatePositionsAndVelocity();
}

void TimeStepPF::updatePositionsAndVelocity()
{
	const auto  numParticles = m_model->numActiveParticles();
	const auto  h            = TimeManager::getCurrent()->getTimeStepSize();
	const auto& x            = VectorXrMap(m_simulationData.getX().data(), m_simulationData.getX().size(), 1);

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			m_model->setPosition(0, i, x.Vec3Block(i));
			auto vel = (m_model->getPosition(0, i) - m_simulationData.getOldPosition(i)) / h;
			m_model->setVelocity(0, i, vel);
		}
	}
}

void SPH::TimeStepPF::addAccellerationToVelocity()
{
	#pragma omp parallel default(shared)
	{
		const auto  numParticles = m_model->numActiveParticles();
		const auto  h = TimeManager::getCurrent()->getTimeStepSize();
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			m_model->setVelocity(0, i, m_model->getVelocity(0, i) + h * m_model->getAcceleration(i));
		}
	}
}

SPH::TimeStepPF::CGSolveState SPH::TimeStepPF::cgSolve()
{
	const auto numVariables = 3u * m_model->numActiveParticles();
	const auto restart_iterations = 50u;
	
	auto x = VectorXrMap(m_simulationData.getX().data(), m_simulationData.getX().size(), 1);
	// initialization of CG
	VectorXr d(numVariables);
	VectorXr r(numVariables);
	VectorXr q(numVariables);
	VectorXr b(numVariables);
	calculateNegativeGradient(r, b, true);
	d = r;

	const auto tol_abs   = 1e-10;
	const auto tol_rel   = 1e-8;
	      auto delta_new = r.squaredNorm();
	const auto delta_0   = delta_new;
	      auto delta_old = std::numeric_limits<Real>::max();

	if ((delta_new < tol_abs) || (delta_new < tol_rel * delta_0))
		return CGSolveState::ALREADY_SOLVED;

	// CG iterations
	for (auto cg_it = 0u; cg_it < numVariables; cg_it++)
	{
		matrixFreeLHS(d, q);
		const auto alpha = delta_new / d.dot(q);
		x += alpha * d;
		
		if ((cg_it + 1) % restart_iterations == 0)
			calculateNegativeGradient(r, b, false);
		else
			r -= alpha * q;

		// test for convergence
		delta_old = delta_new;
		delta_new = r.squaredNorm();
		if ((delta_new < tol_abs) || (delta_new < tol_rel * delta_0))
		{
			return CGSolveState::CONVERGED;
		}

		const auto beta = delta_new / delta_old;
		d *= beta;
		d += r;
	}
	return CGSolveState::MAX_ITER_REACHED;
}

/** \brief calculate the negative gradient for CG*/
void SPH::TimeStepPF::calculateNegativeGradient(VectorXr & r, VectorXr & b, const bool updateRhs)
{
	const auto numVariables = 3u * m_model->numActiveParticles();
	auto x = VectorXrMap(m_simulationData.getX().data(), m_simulationData.getX().size(), 1);
	// use r as temporary buffer for matrix-vector product
	matrixFreeLHS(x, r);
	if (updateRhs)
		matrixFreeRHS(b);
	// -grad_f = b - A*x
	r = b - r;
}

/** \brief compute product of system matrix with x in a matrix-free way and store the result in result*/
void SPH::TimeStepPF::matrixFreeLHS(const VectorXr & x, VectorXr & result)
{
	const auto numParticles = m_model->numActiveParticles();
	const auto numVariables = 3u * numParticles;
	const auto h = TimeManager::getCurrent()->getTimeStepSize();
	
	AtomicRealVec accumulator(numVariables);

	// influence of pressure
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const auto numNeighbors = m_model->numberOfNeighbors(0, i);
			const Vector3r& xi = x.Vec3Block(i);
			for (auto c = 0u; c < 3; c++)
				addToAtomicReal(accumulator[3 * i + c]._a, xi[c]);
			for (auto j = 0u; j < numNeighbors; j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r& xj = x.Vec3Block(neighborIndex);
				for (auto c = 0u; c < 3; c++)
					addToAtomicReal(accumulator[3 * neighborIndex + c]._a, xj[c]);
			}
		}
	}

	// influence of momentum
	#pragma omp parallel default(shared)
	{
		const auto system_scale = h * h * getStiffness();
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const auto m = m_model->getMass(i);
			for (auto c = 0u; c < 3; c++)
			{
				const auto id = 3 * i + c;
				result[id] = system_scale * accumulator[id]._a + m * x[id];
			}
		}
	}
}

/** \brief compute the right hand side of the system in a matrix-free fashion and store the result in result*/
void SPH::TimeStepPF::matrixFreeRHS(VectorXr & result)
{
	const auto numParticles = m_model->numActiveParticles();
	const auto numVariables = 3u * numParticles;
	const auto h            = TimeManager::getCurrent()->getTimeStepSize();

	AtomicRealVec accumulator(numVariables);

	// influence of pressure
	#pragma omp parallel default(shared)
	{
		auto x = VectorXrMap(m_simulationData.getX().data(), m_simulationData.getX().size(), 1);
		const auto density0_inv = 1 / m_model->getDensity0();

		// local step for fluid constraints
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			unsigned int numNeighbors = 0;
			for (unsigned int pid = 0; pid < m_model->numberOfPointSets(); pid++)
			{
				numNeighbors += m_model->numberOfNeighbors(pid, i);
			}
			// particle positions in current constraint, will be projected
			std::vector<Vector3r> p(numNeighbors + 1);
			p[0] = x.Vec3Block(i);
			int counter = 1;
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				p[counter++] = x.Vec3Block(neighborIndex);
			}
			for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
			{
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(pid, i, j);
					p[counter++] = m_model->getPosition(pid, neighborIndex);
				}
			}
			// helper functions
			auto calculateC = [&]() -> Real
			{
				// Compute current density for particle i
				Real density = m_model->getMass(i) * m_model->W_zero();
				const Vector3r &xi = p[0];
				int counter = 1;
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
					const auto& xj = p[counter++];
					density += m_model->getMass(neighborIndex) * m_model->W(xi - xj);
				}
				for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
				{
					for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, i); j++)
					{
						const unsigned int neighborIndex = m_model->getNeighbor(pid, i, j);
						const auto& xj = p[counter++];
						// Boundary: Akinci2012
						density += m_model->getBoundaryPsi(pid, neighborIndex) * m_model->W(xi - xj);
					}
				}
				// constraint value = density / density0 - 1
				const auto C = density * density0_inv - 1;
				// pressure clamping
				return (C < 0) ? 0 : C;
			};
			auto calculateNablaC = [&]() -> VectorXr
			{
				VectorXr nablaC(3 * p.size());
				nablaC.Vec3Block(0).setZero();
				const Vector3r &xi = p[0];



				int counter = 1;
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
					const auto& xj = p[counter];
					nablaC.Vec3Block(counter) = (-density0_inv * m_model->getMass(neighborIndex)) * m_model->gradW(xi - xj);
					nablaC.Vec3Block(0) -= nablaC.Vec3Block(counter);
					counter++;
				}
				for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
				{
					for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, i); j++)
					{
						const unsigned int neighborIndex = m_model->getNeighbor(pid, i, j);
						const auto& xj = p[counter];
						// Boundary: Akinci2012
						nablaC.Vec3Block(counter) = (-density0_inv * m_model->getBoundaryPsi(pid, neighborIndex)) * m_model->gradW(xi - xj);
						nablaC.Vec3Block(0) -= nablaC.Vec3Block(counter);
						counter++;
					}
				}
				return nablaC;
			};
			// projection
			const auto C_goal    = Real(1e-14);
			const auto max_steps = 100u;
			      auto it        = 0u;
			      auto C	     = calculateC();
			while ((std::abs(C) > C_goal) && it++ < max_steps)
			{
				const auto g = calculateNablaC();
				const auto dg = g.squaredNorm();
				if (dg == 0) break;	// found a minimum
				const auto cdg = -C / (dg + 1e-6f); // add regularization factor

				// move fluid particles along constraint gradient
				p[0] += (cdg * m_simulationData.getNumFluidNeighbors(i)) * g.Vec3Block(0);

				// only fluid particles are projected
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
					const auto nfn = m_simulationData.getNumFluidNeighbors(neighborIndex);
					p[j + 1] += (cdg * nfn) * g.Vec3Block(j + 1);
				}

				if (it + 1 < max_steps)
				{
					C = calculateC();
				}
			}
			// update RHS
			for (auto c = 0u; c < 3; c++)
				addToAtomicReal(accumulator[3 * i + c]._a, p[0][c]);

			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				for (auto c = 0u; c < 3; c++)
					addToAtomicReal(accumulator[3 * neighborIndex + c]._a, p[j+1][c]);
			}
		}
	}

	// influence of momentum
	#pragma omp parallel default(shared)
	{
		const auto system_scale = h * h * getStiffness();
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const auto  m = m_model->getMass(i);
			const auto& s = m_simulationData.getS(i);
			for (auto c = 0u; c < 3; c++)
			{
				const auto id = 3 * i + c;
				result[id] = system_scale * accumulator[id]._a + m * s[c];
			}
		}
	}
}

void TimeStepPF::performNeighborhoodSearch()
{
	if (m_counter % 500 == 0)
	{
		m_model->performNeighborhoodSearchSort();
		m_simulationData.performNeighborhoodSearchSort();
		TimeStep::performNeighborhoodSearchSort();
	}
	m_counter++;

	TimeStep::performNeighborhoodSearch();
}

void TimeStepPF::emittedParticles(const unsigned int startIndex)
{

	m_simulationData.emittedParticles(startIndex);
	TimeStep::emittedParticles(startIndex);
}

