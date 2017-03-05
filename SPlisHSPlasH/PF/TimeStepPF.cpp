#include "TimeStepPF.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataPF.h"
#include <iostream>
#include "SPlisHSPlasH/Utilities/Timing.h"

using namespace SPH;
using namespace std;

#define USE_WARMSTART
#define USE_WARMSTART_V

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

	updateVelocity();

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

void TimeStepPF::countFluidNeighbors()
{
	const auto numParticles = m_model->numParticles();
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
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

	countFluidNeighbors();

	for (auto it = 0u; it < m_maxIterations; it++)
	{

	}
}

void TimeStepPF::updateVelocity()
{
	const auto numParticles = m_model->numParticles();
	const auto h = TimeManager::getCurrent()->getTimeStepSize();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &vel = m_model->getVelocity(0, i);
			vel = (m_model->getPosition(0, i) - m_simulationData.getOldPosition(i)) /  h;
		}
	}
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
