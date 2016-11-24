#include "TimeStepIISPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataIISPH.h"
#include <iostream>
#include "SPlisHSPlasH/Utilities/Timing.h"

using namespace SPH;
using namespace std;

TimeStepIISPH::TimeStepIISPH(FluidModel *model) :
	TimeStep(model)
{
	m_simulationData.init(model);
	model->updateBoundaryPsi();
	m_counter = 0;
}

TimeStepIISPH::~TimeStepIISPH(void)
{
}

void TimeStepIISPH::step()
{
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	clearAccelerations();

	performNeighborhoodSearch();

	computeDensities();

	// Compute viscosity 
	computeViscosity();
	computeSurfaceTension();

	updateTimeStepSize();

	// Solve density constraint	
	START_TIMING("predictAdvection");
	predictAdvection();
	STOP_TIMING_AVG;

	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;

	integration();

	// Compute new time	
	tm->setTime (tm->getTime () + h);
}


void TimeStepIISPH::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
}

void TimeStepIISPH::predictAdvection()
{
	const unsigned int numParticles = m_model->numParticles();

	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	// Predict v_adv
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &vel = m_model->getVelocity(0, i);
			const Vector3r &accel = m_model->getAcceleration(i);
			Vector3r &dii = m_simulationData.getDii(i);

			vel += h * accel;

			// Compute d_ii
			dii.setZero();
			const Real density2 = m_model->getDensity(i)*m_model->getDensity(i);
			const Vector3r &xi = m_model->getPosition(0, i);
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
			{
				const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
				const unsigned int &neighborIndex = particleId.point_id;
				const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);

				if (particleId.point_set_id == 0)		// Test if fluid particle
				{
					dii -= m_model->getMass(neighborIndex) / density2 * m_model->gradW(xi - xj);
				}
				else
				{
					dii -= m_model->getBoundaryPsi(particleId.point_set_id, neighborIndex) / density2 * m_model->gradW(xi - xj);
				}
			}
		}
	}

	// Compute rho_adv
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Real &density = m_model->getDensity(i);
			Real &densityAdv = m_simulationData.getDensityAdv(i);
			densityAdv = density;
			const Vector3r &xi = m_model->getPosition(0, i);
			const Vector3r &vi = m_model->getVelocity(0, i);
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
			{
				const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
				const unsigned int &neighborIndex = particleId.point_id;
				const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);
				const Vector3r &vj = m_model->getVelocity(particleId.point_set_id, neighborIndex);

				if (particleId.point_set_id == 0)		// Test if fluid particle
				{
					densityAdv += h*m_model->getMass(neighborIndex) * (vi - vj).dot(m_model->gradW(xi - xj));
				}
				else
				{
					densityAdv += h*m_model->getBoundaryPsi(particleId.point_set_id, neighborIndex) * (vi - vj).dot(m_model->gradW(xi - xj));
				}
			}

			const Real &pressure = m_simulationData.getPressure(i);
			Real &lastPressure = m_simulationData.getLastPressure(i);
			lastPressure = 0.5*pressure;

			// Compute a_ii
			Real &aii = m_simulationData.getAii(i);
			aii = 0.0;
			const Vector3r &dii = m_simulationData.getDii(i);

			const Real dpi = m_model->getMass(i) / (density*density);
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
			{
				const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
				const unsigned int &neighborIndex = particleId.point_id;
				const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);

				if (particleId.point_set_id == 0)		// Test if fluid particle
				{
					// Compute d_ji
					const Vector3r kernel = m_model->gradW(xi - xj);
					const Vector3r dji = dpi * kernel;			

					aii += m_model->getMass(neighborIndex) * (dii - dji).dot(kernel);
				}
				else
				{
					const Vector3r kernel = m_model->gradW(xi - xj);
					const Vector3r dji = dpi * kernel;			
					aii += m_model->getBoundaryPsi(particleId.point_set_id, neighborIndex) * (dii - dji).dot(kernel);
				}
			}
		}
	}
}

void TimeStepIISPH::pressureSolve()
{
	const unsigned int numParticles = m_model->numParticles();

	const Real density0 = m_model->getDensity0();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real h2 = h*h;

	m_iterations = 0;
	const Real omega = 0.5;

	const Real eta = m_maxError * 0.01 * density0;  // maxError is given in percent

	Real avg_density = 0.0;
	while ((((avg_density - density0) > eta) || (m_iterations < 2)) && (m_iterations < m_maxIterations))
	{
		avg_density = 0.0;

		// Compute dij_pj
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int) numParticles; i++)
			{
				Vector3r &dij_pj = m_simulationData.getDij_pj(i);
				dij_pj.setZero();
				const Vector3r &xi = m_model->getPosition(0, i);
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
				{
					const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
					const unsigned int &neighborIndex = particleId.point_id;
					if (particleId.point_set_id == 0)		// Test if fluid particle
					{
						const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);
						const Real &densityj = m_model->getDensity(neighborIndex);
						dij_pj -= m_model->getMass(neighborIndex) / (densityj*densityj) * m_simulationData.getLastPressure(neighborIndex) * m_model->gradW(xi - xj);
					}
				}
			}
		}

		// Compute new pressure
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int) numParticles; i++)
			{
				const Real &aii = m_simulationData.getAii(i);
				const Real density = m_model->getDensity(i);
				const Vector3r &xi = m_model->getPosition(0, i);
				const Real dpi = m_model->getMass(i) / (density*density);
				Real sum = 0.0;
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
				{
					const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
					const unsigned int &neighborIndex = particleId.point_id;
					const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);

					if (particleId.point_set_id == 0)		// Test if fluid particle
					{
						const Vector3r &d_jk_pk = m_simulationData.getDij_pj(neighborIndex);

						// Compute \sum_{k \neq i} djk*pk
						// Compute d_ji
						const Vector3r kernel = m_model->gradW(xi - xj);
						const Vector3r dji = dpi * kernel;			
						const Vector3r d_ji_pi = dji * m_simulationData.getLastPressure(i);

						// \sum ( mj * (\sum dij*pj - djj*pj - \sum_{k \neq i} djk*pk) * m_model->gradW)
						sum += m_model->getMass(neighborIndex) * (m_simulationData.getDij_pj(i) - m_simulationData.getDii(neighborIndex)*m_simulationData.getLastPressure(neighborIndex) - (d_jk_pk - d_ji_pi)).dot(kernel);
					}
					else
					{
						sum += m_model->getBoundaryPsi(particleId.point_set_id, neighborIndex) * m_simulationData.getDij_pj(i).dot(m_model->gradW(xi - xj));
					}
				}

				const Real b = density0 - m_simulationData.getDensityAdv(i);
			
				Real &pi = m_simulationData.getPressure(i);
				const Real &lastPi = m_simulationData.getLastPressure(i);
				const Real denom = aii*h2;
				if (fabs(denom) > 1.0e-9)
					pi = max((1.0 - omega)*lastPi + omega / denom * (b - h2*sum), 0.0);
				else
					pi = 0.0;
 
				if (pi != 0.0)
				{					
					const Real newDensity = (aii*pi + sum)*h2 - b + density0;

					#pragma omp atomic
					avg_density += newDensity;
				}
				else
				{
					#pragma omp atomic
					avg_density += density0;
				}
			}
		}

		for (int i = 0; i < (int)numParticles; i++)
		{
			const Real &pi = m_simulationData.getPressure(i);
			Real &lastPi = m_simulationData.getLastPressure(i);
			lastPi = pi;
		}
		
		avg_density /= numParticles;

		m_iterations++;
	}
}

void TimeStepIISPH::integration()
{
	const unsigned int numParticles = m_model->numParticles();

	// Compute pressure forces
	computePressureAccels();

	Real h = TimeManager::getCurrent()->getTimeStepSize();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numParticles; i++)
		{
			Vector3r &pos = m_model->getPosition(0, i);
			Vector3r &vel = m_model->getVelocity(0, i);
			vel += m_simulationData.getPressureAccel(i) * h;
			pos += vel * h;
		}
	}
}

void TimeStepIISPH::computePressureAccels()
{
	const unsigned int numParticles = m_model->numParticles();

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
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
			{
				const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
				const unsigned int &neighborIndex = particleId.point_id;
				const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);

				if (particleId.point_set_id == 0)		// Test if fluid particle
				{					
					// Pressure 
					const Real &density_j = m_model->getDensity(neighborIndex);

					const Real dpj = m_simulationData.getPressure(neighborIndex) / (density_j*density_j);

					ai -= m_model->getMass(neighborIndex) * (dpi + dpj) * m_model->gradW(xi - xj);
				}
				else
				{
					const Vector3r a = m_model->getBoundaryPsi(particleId.point_set_id, neighborIndex) * (dpi)* m_model->gradW(xi - xj);
					ai -= a;

					m_model->getForce(particleId.point_set_id, neighborIndex) += m_model->getMass(i) * a;
				}
			}
		}
	}
}

void TimeStepIISPH::performNeighborhoodSearch()
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
