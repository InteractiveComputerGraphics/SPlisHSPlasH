#include "TimeStepIISPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataIISPH.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;
using namespace std;

TimeStepIISPH::TimeStepIISPH() :
	TimeStep()
{
	m_simulationData.init();
	m_counter = 0;
}

TimeStepIISPH::~TimeStepIISPH(void)
{
}

void TimeStepIISPH::step()
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getModel();
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	clearAccelerations();

	performNeighborhoodSearch();

	computeDensities();

	sim->computeNonPressureForces();

	sim->updateTimeStepSize();

	// Solve density constraint	
	START_TIMING("predictAdvection");
	predictAdvection();
	STOP_TIMING_AVG;

	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;

	integration();

	sim->emitParticles();

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
	FluidModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numParticles = model->numActiveParticles();

	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	// Predict v_adv
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &vel = model->getVelocity(0, i);
			const Vector3r &accel = model->getAcceleration(i);
			Vector3r &dii = m_simulationData.getDii(i);

			vel += h * accel;

			// Compute d_ii
			dii.setZero();
			const Real density2 = model->getDensity(i)*model->getDensity(i);
			const Vector3r &xi = model->getPosition(0, i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = model->getNeighbor(0, i, j);
				const Vector3r &xj = model->getPosition(0, neighborIndex);
				dii -= model->getMass(neighborIndex) / density2 * model->gradW(xi - xj);
			}

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int pid = 1; pid < model->numberOfPointSets(); pid++)
			{
				for (unsigned int j = 0; j < model->numberOfNeighbors(pid, i); j++)
				{
					const unsigned int neighborIndex = model->getNeighbor(pid, i, j);
					const Vector3r &xj = model->getPosition(pid, neighborIndex);
					dii -= model->getBoundaryPsi(pid, neighborIndex) / density2 * model->gradW(xi - xj);
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
			const Real &density = model->getDensity(i);
			Real &densityAdv = m_simulationData.getDensityAdv(i);
			densityAdv = density;
			const Vector3r &xi = model->getPosition(0, i);
			const Vector3r &vi = model->getVelocity(0, i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = model->getNeighbor(0, i, j);
				const Vector3r &xj = model->getPosition(0, neighborIndex);
				const Vector3r &vj = model->getVelocity(0, neighborIndex);
				densityAdv += h*model->getMass(neighborIndex) * (vi - vj).dot(model->gradW(xi - xj));
			}

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int pid = 1; pid < model->numberOfPointSets(); pid++)
			{
				for (unsigned int j = 0; j < model->numberOfNeighbors(pid, i); j++)
				{
					const unsigned int neighborIndex = model->getNeighbor(pid, i, j);
					const Vector3r &xj = model->getPosition(pid, neighborIndex);
					const Vector3r &vj = model->getVelocity(pid, neighborIndex);
					densityAdv += h*model->getBoundaryPsi(pid, neighborIndex) * (vi - vj).dot(model->gradW(xi - xj));
				}
			}

			const Real &pressure = m_simulationData.getPressure(i);
			Real &lastPressure = m_simulationData.getLastPressure(i);
			lastPressure = 0.5*pressure;

			// Compute a_ii
			Real &aii = m_simulationData.getAii(i);
			aii = 0.0;
			const Vector3r &dii = m_simulationData.getDii(i);

			const Real dpi = model->getMass(i) / (density*density);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = model->getNeighbor(0, i, j);
				const Vector3r &xj = model->getPosition(0, neighborIndex);

				// Compute d_ji
				const Vector3r kernel = model->gradW(xi - xj);
				const Vector3r dji = dpi * kernel;			
				aii += model->getMass(neighborIndex) * (dii - dji).dot(kernel);
			}

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int pid = 1; pid < model->numberOfPointSets(); pid++)
			{
				for (unsigned int j = 0; j < model->numberOfNeighbors(pid, i); j++)
				{
					const unsigned int neighborIndex = model->getNeighbor(pid, i, j);
					const Vector3r &xj = model->getPosition(pid, neighborIndex); 
					const Vector3r kernel = model->gradW(xi - xj);
					const Vector3r dji = dpi * kernel;			
					aii += model->getBoundaryPsi(pid, neighborIndex) * (dii - dji).dot(kernel);
				}
			}
		}
	}
}

void TimeStepIISPH::pressureSolve()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numParticles = model->numActiveParticles();

	const Real density0 = model->getValue<Real>(FluidModel::DENSITY0);
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
				const Vector3r &xi = model->getPosition(0, i);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = model->getNeighbor(0, i, j);
					const Vector3r &xj = model->getPosition(0, neighborIndex);
					const Real &densityj = model->getDensity(neighborIndex);
					dij_pj -= model->getMass(neighborIndex) / (densityj*densityj) * m_simulationData.getLastPressure(neighborIndex) * model->gradW(xi - xj);
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
				const Real density = model->getDensity(i);
				const Vector3r &xi = model->getPosition(0, i);
				const Real dpi = model->getMass(i) / (density*density);
				Real sum = 0.0;

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = model->getNeighbor(0, i, j);
					const Vector3r &xj = model->getPosition(0, neighborIndex);

					const Vector3r &d_jk_pk = m_simulationData.getDij_pj(neighborIndex);

					// Compute \sum_{k \neq i} djk*pk
					// Compute d_ji
					const Vector3r kernel = model->gradW(xi - xj);
					const Vector3r dji = dpi * kernel;
					const Vector3r d_ji_pi = dji * m_simulationData.getLastPressure(i);

					// \sum ( mj * (\sum dij*pj - djj*pj - \sum_{k \neq i} djk*pk) * model->gradW)
					sum += model->getMass(neighborIndex) * (m_simulationData.getDij_pj(i) - m_simulationData.getDii(neighborIndex)*m_simulationData.getLastPressure(neighborIndex) - (d_jk_pk - d_ji_pi)).dot(kernel);
				}

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int pid = 1; pid < model->numberOfPointSets(); pid++)
				{
					for (unsigned int j = 0; j < model->numberOfNeighbors(pid, i); j++)
					{
						const unsigned int neighborIndex = model->getNeighbor(pid, i, j);
						const Vector3r &xj = model->getPosition(pid, neighborIndex);
						sum += model->getBoundaryPsi(pid, neighborIndex) * m_simulationData.getDij_pj(i).dot(model->gradW(xi - xj));
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
	FluidModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numParticles = model->numActiveParticles();

	// Compute pressure forces
	computePressureAccels();

	Real h = TimeManager::getCurrent()->getTimeStepSize();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numParticles; i++)
		{
			Vector3r &pos = model->getPosition(0, i);
			Vector3r &vel = model->getVelocity(0, i);
			vel += m_simulationData.getPressureAccel(i) * h;
			pos += vel * h;
		}
	}
}

void TimeStepIISPH::computePressureAccels()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numParticles = model->numActiveParticles();

	// Compute pressure forces
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = model->getPosition(0, i);
			const Real &density_i = model->getDensity(i);

			Vector3r &ai = m_simulationData.getPressureAccel(i);
			ai.setZero();

			const Real dpi = m_simulationData.getPressure(i) / (density_i*density_i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = model->getNeighbor(0, i, j);
				const Vector3r &xj = model->getPosition(0, neighborIndex);

				// Pressure 
				const Real &density_j = model->getDensity(neighborIndex);
				const Real dpj = m_simulationData.getPressure(neighborIndex) / (density_j*density_j);
				ai -= model->getMass(neighborIndex) * (dpi + dpj) * model->gradW(xi - xj);
			}

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int pid = 1; pid < model->numberOfPointSets(); pid++)
			{
				for (unsigned int j = 0; j < model->numberOfNeighbors(pid, i); j++)
				{
					const unsigned int neighborIndex = model->getNeighbor(pid, i, j);
					const Vector3r &xj = model->getPosition(pid, neighborIndex);
					const Vector3r a = model->getBoundaryPsi(pid, neighborIndex) * (dpi)* model->gradW(xi - xj);
					ai -= a;
					model->getForce(pid, neighborIndex) += model->getMass(i) * a;
				}
			}
		}
	}
}

void TimeStepIISPH::performNeighborhoodSearch()
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

void TimeStepIISPH::emittedParticles(const unsigned int startIndex)
{

	m_simulationData.emittedParticles(startIndex);
	Simulation::getCurrent()->emittedParticles(startIndex);
}

void TimeStepIISPH::resize()
{
	m_simulationData.init();
}
