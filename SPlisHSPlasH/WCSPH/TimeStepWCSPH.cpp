#include "TimeStepWCSPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataWCSPH.h"
#include <iostream>
#include "SPlisHSPlasH/Utilities/Timing.h"

using namespace SPH;
using namespace std;

TimeStepWCSPH::TimeStepWCSPH(FluidModel *model) :
	TimeStep(model)
{
	m_simulationData.init(model);
	m_counter = 0;

	m_stiffness = 50000.0;
	m_exponent = 7.0;
}

TimeStepWCSPH::~TimeStepWCSPH(void)
{
}

void TimeStepWCSPH::step()
{
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	performNeighborhoodSearch();

	// Compute accelerations: a(t)
	clearAccelerations();
	computeDensities();
	computeNonPressureForces();

	const Real stiffness = getStiffness();	
	const Real density0 = m_model->getDensity0();

	const Real exponent = getExponent();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) m_model->numActiveParticles(); i++)
		{
			Real &density = m_model->getDensity(i);
			density = max(density, density0);
			m_simulationData.getPressure(i) = stiffness * (pow(density/density0, exponent) - 1.0);
		}
	}

	computePressureAccels();

	updateTimeStepSize();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)m_model->numActiveParticles(); i++)
		{
			Vector3r &pos = m_model->getPosition(0, i);
			Vector3r &vel = m_model->getVelocity(0, i);
			Vector3r &accel = m_model->getAcceleration(i);
			accel += m_simulationData.getPressureAccel(i);
			vel += accel * h;
			pos += vel * h;
		}
	}

	emitParticles();

	// Compute new time	
	tm->setTime (tm->getTime () + h);
}


void TimeStepWCSPH::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
}

void TimeStepWCSPH::computePressureAccels()
{
	const unsigned int numParticles = m_model->numActiveParticles();

	// Compute pressure forces
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			const Real &density_i = m_model->getDensity(i);

			Vector3r &ai = m_simulationData.getPressureAccel(i);
			ai.setZero();

			const Real dpi = m_simulationData.getPressure(i) / (density_i*density_i);
			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
#
				// Pressure 
				const Real &density_j = m_model->getDensity(neighborIndex);
				const Real dpj = m_simulationData.getPressure(neighborIndex) / (density_j*density_j);
				ai -= m_model->getMass(neighborIndex) * (dpi + dpj) * m_model->gradW(xi - xj);
			}

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
			{
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(pid, i, j);
					const Vector3r &xj = m_model->getPosition(pid, neighborIndex);

					const Vector3r a = m_model->getBoundaryPsi(pid, neighborIndex) * (dpi)* m_model->gradW(xi - xj);
					ai -= a;

					m_model->getForce(pid, neighborIndex) += m_model->getMass(i) * a;
				}
			}
		}
	}
}

void TimeStepWCSPH::performNeighborhoodSearch()
{
	const unsigned int numParticles = m_model->numActiveParticles();

	if (m_counter % 500 == 0)
	{
		m_model->performNeighborhoodSearchSort();
		m_simulationData.performNeighborhoodSearchSort();
		TimeStep::performNeighborhoodSearchSort();
	}
	m_counter++;

	TimeStep::performNeighborhoodSearch();
}

void TimeStepWCSPH::emittedParticles(const unsigned int startIndex)
{
	m_simulationData.emittedParticles(startIndex);
	TimeStep::emittedParticles(startIndex);
}
