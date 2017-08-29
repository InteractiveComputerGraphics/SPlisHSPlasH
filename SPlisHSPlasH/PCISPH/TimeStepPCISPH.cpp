#include "TimeStepPCISPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataPCISPH.h"
#include <iostream>
#include "SPlisHSPlasH/Utilities/Timing.h"

using namespace SPH;
using namespace std;

TimeStepPCISPH::TimeStepPCISPH(FluidModel *model) :
	TimeStep(model)
{
	m_simulationData.init(model);
	m_counter = 0;
}

TimeStepPCISPH::~TimeStepPCISPH(void)
{
}

void TimeStepPCISPH::step()
{
	TimeManager *tm = TimeManager::getCurrent();
	const int numParticles = (int)m_model->numActiveParticles();
	const Real h = tm->getTimeStepSize();

	performNeighborhoodSearch();

	// Compute accelerations: a(t)
	clearAccelerations();
	computeDensities();
	computeNonPressureForces();

	updateTimeStepSize();

	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			const Vector3r accel = m_model->getAcceleration(i) + m_simulationData.getPressureAccel(i);
			const Vector3r &lastX = m_simulationData.getLastPosition(i);
			const Vector3r &lastV = m_simulationData.getLastVelocity(i);
			Vector3r &x = m_model->getPosition(0, i);
			Vector3r &v = m_model->getVelocity(0, i);
			v = lastV + h*accel;
			x = lastX + h*v;
		}
	}

	emitParticles();

	// Compute new time		
	tm->setTime(tm->getTime() + h);
}

void TimeStepPCISPH::pressureSolve()
{
	const int numParticles = (int)m_model->numActiveParticles();
	const Real density0 = m_model->getDensity0();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real h2 = h*h;
	const Real invH2 = 1.0 / h2;

	// Initialization
	for (int i = 0; i <numParticles; i++)
	{
		Vector3r &lastX = m_simulationData.getLastPosition(i);
		Vector3r &lastV = m_simulationData.getLastVelocity(i);
		lastX = m_model->getPosition(0, i);
		lastV = m_model->getVelocity(0, i);
		m_simulationData.getPressure(i) = 0.0;
		m_simulationData.getPressureAccel(i).setZero();
	}

	Real avg_density_err = 0;
	m_iterations = 0;

	// Maximal allowed density fluctuation
	const Real eta = m_maxError * 0.01 * density0;  // maxError is given in percent

	while (((avg_density_err > eta) || (m_iterations < 3)) && (m_iterations < m_maxIterations))
	{
		avg_density_err = 0.0;

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				const Vector3r accel = m_model->getAcceleration(i) + m_simulationData.getPressureAccel(i);
				const Vector3r &lastX = m_simulationData.getLastPosition(i);
				const Vector3r &lastV = m_simulationData.getLastVelocity(i);
				Vector3r &x = m_model->getPosition(0, i);
				Vector3r &v = m_model->getVelocity(0, i);
				v = lastV + h*accel;
				x = lastX + h*v;
			}
		}


		// Predict density 
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				const Vector3r &xi = m_model->getPosition(0, i);
				Real &densityAdv = m_simulationData.getDensityAdv(i);
				densityAdv = m_model->getMass(i) * m_model->W_zero();
				
				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
					const Vector3r &xj = m_model->getPosition(0, neighborIndex);
					densityAdv += m_model->getMass(neighborIndex) * m_model->W(xi - xj);
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
						densityAdv += m_model->getBoundaryPsi(pid, neighborIndex) * m_model->W(xi - xj);
					}
				}

				densityAdv = max(densityAdv, density0);
				const Real density_err = densityAdv - density0;
				Real &pressure = m_simulationData.getPressure(i);
				pressure += invH2 * m_simulationData.getPCISPH_ScalingFactor() * density_err;

				#pragma omp atomic
				avg_density_err += density_err;
			}
		}

		avg_density_err /= numParticles;

		// Compute pressure forces
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				const Vector3r &xi = m_simulationData.getLastPosition(i);

				Vector3r &ai = m_simulationData.getPressureAccel(i);
				ai.setZero();

				const Real dpi = m_simulationData.getPressure(i) / (density0*density0);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
					// Pressure 
					const Real dpj = m_simulationData.getPressure(neighborIndex) / (density0*density0);
					const Vector3r &xj = m_simulationData.getLastPosition(neighborIndex);
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
						// Pressure 
						const Vector3r &xj = m_model->getPosition(pid, neighborIndex);
						const Vector3r a = m_model->getBoundaryPsi(pid, neighborIndex) * (dpi)* m_model->gradW(xi - xj);
						ai -= a;

						m_model->getForce(pid, neighborIndex) += m_model->getMass(i) * a;
					}
				}
			}
		}


		if (m_iterations > m_maxIterations)
			break;

		m_iterations++;
	}
}


void TimeStepPCISPH::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
}

void TimeStepPCISPH::performNeighborhoodSearch()
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

void TimeStepPCISPH::emittedParticles(const unsigned int startIndex)
{

	m_simulationData.emittedParticles(startIndex);
	TimeStep::emittedParticles(startIndex);
}
