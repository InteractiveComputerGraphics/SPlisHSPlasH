#include "TimeStepPF.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataPF.h"
#include <iostream>
#include "SPlisHSPlasH/Utilities/Timing.h"

using namespace SPH;
using namespace std;

TimeStepPF::TimeStepPF(FluidModel *model) :
	TimeStep(model)
{
	m_simulationData.init(model);
	model->updateBoundaryPsi();
	m_counter = 0;
	m_iterationsV = 0;
}

TimeStepPF::~TimeStepPF(void)
{
}

void TimeStepPF::step()
{
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	const unsigned int numParticles = m_model->numParticles();

	clearAccelerations();
	initialGuessForPositions();
	performNeighborhoodSearch();

	START_TIMING("solvePDConstraints");
	solvePDConstraints();
	STOP_TIMING_AVG;

	computeSurfaceTension();
	computeViscosity();

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
	const auto numParticles = m_model->numParticles();
	const auto h = TimeManager::getCurrent()->getTimeStepSize();

#pragma omp parallel for
	for (int i = 0; i < numParticles; i++)
	{
		m_simulationData.setOldPosition(i, m_model->getPosition(0, i));
		const auto newPos = m_model->getPosition(0, i) + h * m_model->getVelocity(0, i) + (h * h) * m_model->getAcceleration(i);
		m_model->setPosition(0, i, newPos);
	}
}

void TimeStepPF::prepareSolve()
{
	const auto numParticles = m_model->numParticles();
	auto&      x            = m_simulationData.getX();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			x[3 * i + 0] = m_model->getPosition(0, i)[0];
			x[3 * i + 1] = m_model->getPosition(0, i)[1];
			x[3 * i + 2] = m_model->getPosition(0, i)[2];

			auto nfn = 0u;
			for (auto j = 0u; j < m_model->numberOfNeighbors(i); j++)
				if (m_model->getNeighbor(i, j).point_set_id == 0u)
					nfn++;
			m_simulationData.setNumFluidNeighbors(i, nfn);
		}
	}
}

void TimeStepPF::solvePDConstraints()
{
	const auto numParticles = m_model->numParticles();

	prepareSolve();

	for (auto it = 0u; it < m_maxIterations; it++)
	{
		const auto s = cgSolve();
		if (s == CGSolveState::ALREADY_SOLVED) break;
	}

	updatePositionsAndVelocity();
}

void TimeStepPF::updatePositionsAndVelocity()
{
	const auto  numParticles = m_model->numParticles();
	const auto  h            = TimeManager::getCurrent()->getTimeStepSize();
	const auto& x            = m_simulationData.getX();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			
			m_model->setPosition(0, i,
				{
					x[3 * i + 0],
					x[3 * i + 1],
					x[3 * i + 2]
				}
			);
			Vector3r &vel = m_model->getVelocity(0, i);
			vel = (m_model->getPosition(0, i) - m_simulationData.getOldPosition(i)) /  h;
		}
	}
}

SPH::TimeStepPF::CGSolveState SPH::TimeStepPF::cgSolve()
{
	const auto numVariables = 3u * m_model->numParticles();
	const auto restart_iterations = 50u;
	
	auto x = VectorXrMap(m_simulationData.getX().data(), m_simulationData.getX().size(), 1);
	// initialization of CG
	VectorXr d(numVariables);
	VectorXr r(numVariables);
	VectorXr q(numVariables);
	VectorXr b(numVariables);
	calculateNegativeGradient(r, b, true);
	d = r;

	auto delta_new = r.squaredNorm();
	auto delta_0   = delta_new;
	auto delta_old = std::numeric_limits<Real>::max();

	if ((delta_new < 1e-12) || (delta_new < 1e-10 * delta_0))
		return CGSolveState::ALREADY_SOLVED;

	// CG iterations
	for (auto cg_it = 0u; cg_it < numVariables; cg_it++)
	{
		matrixFreeLHS(d, q);
		const auto alpha = delta_new / d.dot(q);
		x = alpha * d;
		
		if ((cg_it + 1) % restart_iterations == 0)
			calculateNegativeGradient(r, b, false);
		else
			r = -alpha * q;

		// test for convergence
		delta_old = delta_new;
		delta_new = r.squaredNorm();
		if ((delta_new < 1e-12) || (delta_new < 1e-10 * delta_0))
		{
			return CGSolveState::CONVERGED;
		}

		const auto beta = delta_new / delta_old;
		d *= beta;
		d += r;
	}
	return CGSolveState::MAX_ITER_REACHED;
}

void SPH::TimeStepPF::calculateNegativeGradient(VectorXr & r, VectorXr & b, const bool updateRhs)
{
	const auto numVariables = 3u * m_model->numParticles();
	auto x = VectorXrMap(m_simulationData.getX().data(), m_simulationData.getX().size(), 1);

	matrixFreeLHS(x, r);
	if (updateRhs)
		matrixFreeRHS(b);
	r = b - r;
}

void SPH::TimeStepPF::matrixFreeLHS(const VectorXr & v, VectorXr & result)
{

}

void SPH::TimeStepPF::matrixFreeRHS(VectorXr & result)
{

}

void TimeStepPF::performNeighborhoodSearch()
{
	const unsigned int numParticles = m_model->numParticles();
	const Real supportRadius = m_model->getSupportRadius();

	if (m_counter % 100 == 0)
	{
		m_model->performNeighborhoodSearchSort();
		m_simulationData.performNeighborhoodSearchSort();
	}
	m_counter++;

	TimeStep::performNeighborhoodSearch();
}
