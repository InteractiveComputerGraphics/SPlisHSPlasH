#include "TimeStepPBF.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "TimeIntegration.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataPBF.h"
#include <iostream>
#include "SPlisHSPlasH/Utilities/Timing.h"

using namespace SPH;
using namespace std;

TimeStepPBF::TimeStepPBF(FluidModel *model) :
	TimeStep(model)
{
	m_simulationData.init(model);
	model->updateBoundaryPsi();
	m_counter = 0;
}

TimeStepPBF::~TimeStepPBF(void)
{
}

void TimeStepPBF::step()
{
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	clearAccelerations();

	// Time integration
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) m_model->numParticles(); i++)
		{
			m_simulationData.getDeltaX(i).setZero();
			m_simulationData.getLastPosition(i) = m_simulationData.getOldPosition(i);
			m_simulationData.getOldPosition(i) = m_model->getPosition(0, i);
			TimeIntegration::semiImplicitEuler(h, m_model->getMass(i), m_model->getPosition(0, i), m_model->getVelocity(0, i), m_model->getAcceleration(i));
		}
	}

	// Perform neighborhood search
	performNeighborhoodSearch();

	// Solve density constraint
	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;

	// Update velocities	
	if (m_model->getVelocityUpdateMethod() == 0)
	{
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)m_model->numParticles(); i++)
			{
				TimeIntegration::velocityUpdateFirstOrder(h, m_model->getMass(i), m_model->getPosition(0, i), m_simulationData.getOldPosition(i), m_model->getVelocity(0, i));
				m_model->getAcceleration(i).setZero();
			}
		}
	}
	else
	{
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)m_model->numParticles(); i++)
			{
				TimeIntegration::velocityUpdateSecondOrder(h, m_model->getMass(i), m_model->getPosition(0, i), m_simulationData.getOldPosition(i), m_simulationData.getLastPosition(i), m_model->getVelocity(0, i));
				m_model->getAcceleration(i).setZero();
			}
		}
	}

	// Compute viscosity 
	computeViscosity();
	computeSurfaceTension();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)m_model->numParticles(); i++)
		{
			m_model->getVelocity(0, i) += h * m_model->getAcceleration(i);
		}
	}

	updateTimeStepSize();

	// Compute new time	
	tm->setTime (tm->getTime () + h);
}


void TimeStepPBF::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
}


void TimeStepPBF::pressureSolve()
{
	m_iterations = 0;
	const Real eps = 1.0e-6;

	const unsigned int numParticles = m_model->numParticles();
	const Real invH = 1.0 / TimeManager::getCurrent()->getTimeStepSize();
	const Real invH2 = invH*invH;

	const Real density0 = m_model->getDensity0();

	const Real eta = m_maxError * 0.01 * density0;  // maxError is given in percent

	Real avg_density_err = 0.0;
	while (((avg_density_err > eta) || (m_iterations < 2)) && (m_iterations < m_maxIterations))
	{
		avg_density_err = 0.0;

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int) numParticles; i++)
			{
				Real &density = m_model->getDensity(i);				

				// Compute current density for particle i
				density = m_model->getMass(i) * m_model->W_zero();
				const Vector3r &xi = m_model->getPosition(0, i);
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
				{
					const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
					const unsigned int &neighborIndex = particleId.point_id;
					const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);

					if (particleId.point_set_id == 0)		// Test if fluid particle
					{
						density += m_model->getMass(neighborIndex) * m_model->W(xi - xj);
					}
					else 
					{
						// Boundary: Akinci2012
						density += m_model->getBoundaryPsi(particleId.point_set_id, neighborIndex) * m_model->W(xi - xj);
					}
				}

				const Real density_err = max(density, density0) - density0;
				#pragma omp atomic
				avg_density_err += density_err / numParticles;

				// Evaluate constraint function
				const Real C = std::max(density / density0 - 1.0, 0.0);			// clamp to prevent particle clumping at surface

				if (C != 0.0)
				{
					// Compute gradients dC/dx_j 
					Real sum_grad_C2 = 0.0;
					Vector3r gradC_i(0.0, 0.0, 0.0);

					for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
					{
						const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
						const unsigned int &neighborIndex = particleId.point_id;
						const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);

						if (particleId.point_set_id == 0)		// Test if fluid particle
						{
							const Vector3r gradC_j = -m_model->getMass(neighborIndex) / density0 * m_model->gradW(xi - xj);
							sum_grad_C2 += gradC_j.squaredNorm();
							gradC_i -= gradC_j;
						}
						else
						{
							// Boundary: Akinci2012
							const Vector3r gradC_j = -m_model->getBoundaryPsi(particleId.point_set_id, neighborIndex) / density0 * m_model->gradW(xi - xj);
							sum_grad_C2 += gradC_j.squaredNorm();
							gradC_i -= gradC_j;
						}
					}

					sum_grad_C2 += gradC_i.squaredNorm();

					// Compute lambda
					Real &lambda = m_simulationData.getLambda(i);
					lambda = -C / (sum_grad_C2 + eps);
				}
				else
					m_simulationData.getLambda(i) = 0.0;
			}


			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				Vector3r &corr = m_simulationData.getDeltaX(i);
				
				// Compute position correction
				corr.setZero();
				const Vector3r &xi = m_model->getPosition(0, i);
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
				{
					const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
					const unsigned int &neighborIndex = particleId.point_id;
					const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);

					if (particleId.point_set_id == 0)		// Test if fluid particle
					{
						const Vector3r gradC_j = -m_model->getMass(neighborIndex) / density0 * m_model->gradW(xi - xj);
						corr -= (m_simulationData.getLambda(i) + m_simulationData.getLambda(neighborIndex)) * gradC_j;
					}
					else 
					{
						// Boundary: Akinci2012
						const Vector3r gradC_j = -m_model->getBoundaryPsi(particleId.point_set_id, neighborIndex) / density0 * m_model->gradW(xi - xj);
						const Vector3r dx = 2.0 * m_simulationData.getLambda(i) * gradC_j;
						corr -= dx;

						m_model->getForce(particleId.point_set_id, neighborIndex) += m_model->getMass(i) * dx * invH2;
					}
				}
			}

			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				m_model->getPosition(0, i) += m_simulationData.getDeltaX(i);
			}
		}

		m_iterations++;
	}
}

void TimeStepPBF::performNeighborhoodSearch()
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
