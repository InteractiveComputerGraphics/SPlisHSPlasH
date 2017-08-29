#include "TimeStepDFSPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataDFSPH.h"
#include <iostream>
#include "SPlisHSPlasH/Utilities/Timing.h"

using namespace SPH;
using namespace std;

#define USE_WARMSTART
#define USE_WARMSTART_V

TimeStepDFSPH::TimeStepDFSPH(FluidModel *model) :
	TimeStep(model),
	m_simulationData()
{
	m_simulationData.init(model);
	m_counter = 0;
	m_iterationsV = 0;
	m_enableDivergenceSolver = true;
}

TimeStepDFSPH::~TimeStepDFSPH(void)
{
}

void TimeStepDFSPH::step()
{
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	const unsigned int numParticles = m_model->numActiveParticles();

	performNeighborhoodSearch();

	computeDensities();

	START_TIMING("computeDFSPHFactor");
	computeDFSPHFactor();
	STOP_TIMING_AVG;

	if (m_enableDivergenceSolver)
	{
		START_TIMING("divergenceSolve");
		divergenceSolve();
		STOP_TIMING_AVG
	}
	else
		m_iterationsV = 0;

	// Compute accelerations: a(t)
	clearAccelerations();

	computeNonPressureForces();

	updateTimeStepSize();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &vel = m_model->getVelocity(0, i);
			vel += h * m_model->getAcceleration(i);
		}
	}

	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &xi = m_model->getPosition(0, i);
			const Vector3r &vi = m_model->getVelocity(0, i);
			xi += h * vi;
		}
	}

	emitParticles();

	// Compute new time	
	tm->setTime (tm->getTime () + h);
}

void TimeStepDFSPH::computeDFSPHFactor()
{
	//////////////////////////////////////////////////////////////////////////
	// Init parameters
	//////////////////////////////////////////////////////////////////////////

	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const int numParticles = (int) m_model->numActiveParticles();

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure stiffness denominator
		//////////////////////////////////////////////////////////////////////////

		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			//////////////////////////////////////////////////////////////////////////
			// Compute gradient dp_i/dx_j * (1/k)  and dp_j/dx_j * (1/k)
			//////////////////////////////////////////////////////////////////////////
			const Vector3r &xi = m_model->getPosition(0, i);
			Real sum_grad_p_k = 0.0;
			Vector3r grad_p_i;
			grad_p_i.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Vector3r grad_p_j = -m_model->getMass(neighborIndex) * m_model->gradW(xi - xj);
				sum_grad_p_k += grad_p_j.squaredNorm();
				grad_p_i -= grad_p_j;
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
					const Vector3r grad_p_j = -m_model->getBoundaryPsi(pid, neighborIndex) * m_model->gradW(xi - xj);
					sum_grad_p_k += grad_p_j.squaredNorm();
					grad_p_i -= grad_p_j;
				}
			}

			sum_grad_p_k += grad_p_i.squaredNorm();

			//////////////////////////////////////////////////////////////////////////
			// Compute pressure stiffness denominator
			//////////////////////////////////////////////////////////////////////////
			Real &factor = m_simulationData.getFactor(i);

			sum_grad_p_k = max(sum_grad_p_k, m_eps);
			factor = -1.0 / (sum_grad_p_k);
		}
	}
}

void TimeStepDFSPH::pressureSolve()
{
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real h2 = h*h;
	const Real invH = 1.0 / h;
	const Real invH2 = 1.0/h2;
	const Real density0 = m_model->getDensity0();
	const int numParticles = (int)m_model->numActiveParticles();
	Real avg_density_err = 0.0;

#ifdef USE_WARMSTART			
	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Divide by h^2, the time step size has been removed in 
		// the last step to make the stiffness value independent 
		// of the time step size
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			m_simulationData.getKappa(i) = max(m_simulationData.getKappa(i)*invH2, -0.5);
			//computeDensityAdv(i, numParticles, h, density0);
		}

		//////////////////////////////////////////////////////////////////////////
		// Predict v_adv with external velocities
		////////////////////////////////////////////////////////////////////////// 

		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			//if (m_simulationData.getDensityAdv(i) > density0)
			{
				Vector3r &vel = m_model->getVelocity(0, i);
				const Real ki = m_simulationData.getKappa(i);
				const Vector3r &xi = m_model->getPosition(0, i);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
					const Real kj = m_simulationData.getKappa(neighborIndex);

					const Real kSum = (ki + kj);
					if (fabs(kSum) > m_eps)
					{
						const Vector3r &xj = m_model->getPosition(0, neighborIndex);
						const Vector3r grad_p_j = -m_model->getMass(neighborIndex) * m_model->gradW(xi - xj);
						vel -= h * kSum * grad_p_j;					// ki, kj already contain inverse density
					}
				}

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				if (fabs(ki) > m_eps)
				{
					for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
					{
						for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, i); j++)
						{
							const unsigned int neighborIndex = m_model->getNeighbor(pid, i, j);
							const Vector3r &xj = m_model->getPosition(pid, neighborIndex);
							const Vector3r grad_p_j = -m_model->getBoundaryPsi(pid, neighborIndex) * m_model->gradW(xi - xj);
							const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
							vel += velChange;

							m_model->getForce(pid, neighborIndex) -= m_model->getMass(i) * velChange * invH;
						}
					}
				}
			}
		}
	}
#endif

	//////////////////////////////////////////////////////////////////////////
	// Compute rho_adv
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			computeDensityAdv(i, numParticles, h, density0);
			m_simulationData.getFactor(i) *= invH2;
#ifdef USE_WARMSTART
			m_simulationData.getKappa(i) = 0.0;
#endif
		}
	}

	m_iterations = 0;

	//////////////////////////////////////////////////////////////////////////
	// Start solver
	//////////////////////////////////////////////////////////////////////////
	
	// Maximal allowed density fluctuation
	const Real eta = m_maxError * 0.01 * density0;  // maxError is given in percent
	
	while (((avg_density_err > eta) || (m_iterations < 2)) && (m_iterations < m_maxIterations))
	{
		avg_density_err = 0.0;

		#pragma omp parallel default(shared)
		{
			//////////////////////////////////////////////////////////////////////////
			// Compute pressure forces
			//////////////////////////////////////////////////////////////////////////
			#pragma omp for schedule(static) 
			for (int i = 0; i < numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Evaluate rhs
				//////////////////////////////////////////////////////////////////////////
				const Real b_i = m_simulationData.getDensityAdv(i) - density0;
				const Real ki = b_i*m_simulationData.getFactor(i);
#ifdef USE_WARMSTART
				m_simulationData.getKappa(i) += ki;
#endif

				Vector3r &v_i = m_model->getVelocity(0, i);
				const Vector3r &xi = m_model->getPosition(0, i);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
					const Real b_j = m_simulationData.getDensityAdv(neighborIndex) - density0;
					const Real kj = b_j*m_simulationData.getFactor(neighborIndex);
					const Real kSum = (ki + kj);
					if (fabs(kSum) > m_eps)
					{
						const Vector3r &xj = m_model->getPosition(0, neighborIndex);
						const Vector3r grad_p_j = -m_model->getMass(neighborIndex) * m_model->gradW(xi - xj);

						// Directly update velocities instead of storing pressure accelerations
						v_i -= h * kSum * grad_p_j;			// ki, kj already contain inverse density						
					}
				}

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				if (fabs(ki) > m_eps)
				{
					for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
					{
						for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, i); j++)
						{
							const unsigned int neighborIndex = m_model->getNeighbor(pid, i, j);
							const Vector3r &xj = m_model->getPosition(pid, neighborIndex);
							const Vector3r grad_p_j = -m_model->getBoundaryPsi(pid, neighborIndex) * m_model->gradW(xi - xj);

							// Directly update velocities instead of storing pressure accelerations
							const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
							v_i += velChange;

							m_model->getForce(pid, neighborIndex) -= m_model->getMass(i) * velChange * invH;
						}
					}
				}
			}

		
			//////////////////////////////////////////////////////////////////////////
			// Update rho_adv and density error
			//////////////////////////////////////////////////////////////////////////
			#pragma omp for reduction(+:avg_density_err) schedule(static) 
			for (int i = 0; i < numParticles; i++)
			{
				computeDensityAdv(i, numParticles, h, density0);

				const Real density_err = m_simulationData.getDensityAdv(i) - density0;
				avg_density_err += density_err;
			}
		}

		avg_density_err /= numParticles;

		m_iterations++;
	}


#ifdef USE_WARMSTART
	//////////////////////////////////////////////////////////////////////////
	// Multiply by h^2, the time step size has to be removed 
	// to make the stiffness value independent 
	// of the time step size
	//////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < numParticles; i++)
		m_simulationData.getKappa(i) *= h2;
#endif
}

void TimeStepDFSPH::divergenceSolve()
{
	//////////////////////////////////////////////////////////////////////////
	// Init parameters
	//////////////////////////////////////////////////////////////////////////

	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = 1.0 / h;
	const int numParticles = (int)m_model->numActiveParticles();
	const unsigned int maxIter = m_maxIterationsV;
	const Real maxError = m_maxErrorV;
	const Real density0 = m_model->getDensity0();


#ifdef USE_WARMSTART_V
	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Divide by h^2, the time step size has been removed in 
		// the last step to make the stiffness value independent 
		// of the time step size
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			m_simulationData.getKappaV(i) = 0.5*max(m_simulationData.getKappaV(i)*invH, -0.5);
			computeDensityChange(i, h, density0);
		}

		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (m_simulationData.getDensityAdv(i) > 0.0)
			{
				Vector3r &vel = m_model->getVelocity(0, i);
				const Real ki = m_simulationData.getKappaV(i);
				const Vector3r &xi = m_model->getPosition(0, i);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
					const Real kj = m_simulationData.getKappaV(neighborIndex);

					const Real kSum = (ki + kj);
					if (fabs(kSum) > m_eps)
					{
						const Vector3r &xj = m_model->getPosition(0, neighborIndex);
						const Vector3r grad_p_j = -m_model->getMass(neighborIndex) * m_model->gradW(xi - xj);
						vel -= h * kSum * grad_p_j;					// ki, kj already contain inverse density
					}
				}

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				if (fabs(ki) > m_eps)
				{
					for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
					{
						for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, i); j++)
						{
							const unsigned int neighborIndex = m_model->getNeighbor(pid, i, j);
							const Vector3r &xj = m_model->getPosition(pid, neighborIndex);
							const Vector3r grad_p_j = -m_model->getBoundaryPsi(pid, neighborIndex) * m_model->gradW(xi - xj);

							const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
							vel += velChange;

							m_model->getForce(pid, neighborIndex) -= m_model->getMass(i) * velChange * invH;
						}
					}
				}
			}
		}
	}
#endif

	//////////////////////////////////////////////////////////////////////////
	// Compute velocity of density change
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			computeDensityChange(i, h, density0);
			m_simulationData.getFactor(i) *= invH;

#ifdef USE_WARMSTART_V
			m_simulationData.getKappaV(i) = 0.0;
#endif
		}
	}

	m_iterationsV = 0;

	//////////////////////////////////////////////////////////////////////////
	// Start solver
	//////////////////////////////////////////////////////////////////////////
	
	// Maximal allowed density fluctuation
	// use maximal density error divided by time step size
	const Real eta = (1.0/h) * maxError * 0.01 * density0;  // maxError is given in percent
	
	Real avg_density_err = 0.0;
	while (((avg_density_err > eta) || (m_iterationsV < 1)) && (m_iterationsV < maxIter))
	{
		avg_density_err = 0.0;
		
		//////////////////////////////////////////////////////////////////////////
		// Perform Jacobi iteration over all blocks
		//////////////////////////////////////////////////////////////////////////	
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static) 
			for (int i = 0; i < (int)numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Evaluate rhs
				//////////////////////////////////////////////////////////////////////////
				const Real b_i = m_simulationData.getDensityAdv(i);
				const Real ki = b_i*m_simulationData.getFactor(i);
#ifdef USE_WARMSTART_V
				m_simulationData.getKappaV(i) += ki;
#endif

				Vector3r &v_i = m_model->getVelocity(0, i);

				const Vector3r &xi = m_model->getPosition(0, i);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);					
					const Real b_j = m_simulationData.getDensityAdv(neighborIndex);
					const Real kj = b_j*m_simulationData.getFactor(neighborIndex);

					const Real kSum = (ki + kj);
					if (fabs(kSum) > m_eps)
					{
						const Vector3r &xj = m_model->getPosition(0, neighborIndex);
						const Vector3r grad_p_j = -m_model->getMass(neighborIndex) * m_model->gradW(xi - xj);
						v_i -= h * kSum * grad_p_j;			// ki, kj already contain inverse density
					}
				}

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				if (fabs(ki) > m_eps)
				{
					for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
					{
						for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, i); j++)
						{
							const unsigned int neighborIndex = m_model->getNeighbor(pid, i, j);
							const Vector3r &xj = m_model->getPosition(pid, neighborIndex);
							const Vector3r grad_p_j = -m_model->getBoundaryPsi(pid, neighborIndex) * m_model->gradW(xi - xj);

							const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
							v_i += velChange;

							m_model->getForce(pid, neighborIndex) -= m_model->getMass(i) * velChange * invH;
						}
					}
				}
			}

			//////////////////////////////////////////////////////////////////////////
			// Update rho_adv and density error
			//////////////////////////////////////////////////////////////////////////
			#pragma omp for reduction(+:avg_density_err) schedule(static) 
			for (int i = 0; i < (int)numParticles; i++)
			{
				computeDensityChange(i, h, density0);
				avg_density_err += m_simulationData.getDensityAdv(i);
			}
		}	
	
		avg_density_err /= numParticles;
		m_iterationsV++;
	}

#ifdef USE_WARMSTART_V
	//////////////////////////////////////////////////////////////////////////
	// Multiply by h, the time step size has to be removed 
	// to make the stiffness value independent 
	// of the time step size
	//////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < numParticles; i++)
		m_simulationData.getKappaV(i) *= h;
 #endif

	for (int i = 0; i < numParticles; i++)
	{
		m_simulationData.getFactor(i) *= h;
	}
}


void TimeStepDFSPH::computeDensityAdv(const unsigned int index, const int numParticles, const Real h, const Real density0)
{
	const Real &density = m_model->getDensity(index);
	Real &densityAdv = m_simulationData.getDensityAdv(index);
	const Vector3r &xi = m_model->getPosition(0, index);
	const Vector3r &vi = m_model->getVelocity(0, index);
	Real delta = 0.0;

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, index); j++)
	{
		const unsigned int neighborIndex = m_model->getNeighbor(0, index, j);
		const Vector3r &xj = m_model->getPosition(0, neighborIndex);
		const Vector3r &vj = m_model->getVelocity(0, neighborIndex);
		delta += m_model->getMass(neighborIndex) * (vi - vj).dot(m_model->gradW(xi - xj));
	}

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
	{
		for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, index); j++)
		{
			const unsigned int neighborIndex = m_model->getNeighbor(pid, index, j);
			const Vector3r &xj = m_model->getPosition(pid, neighborIndex);
			const Vector3r &vj = m_model->getVelocity(pid, neighborIndex);
			delta += m_model->getBoundaryPsi(pid, neighborIndex) * (vi - vj).dot(m_model->gradW(xi - xj));
		}
	}

	densityAdv = density + h*delta;
	densityAdv = max(densityAdv, density0);
}

void TimeStepDFSPH::computeDensityChange(const unsigned int index, const Real h, const Real density0)
{
	Real &densityAdv = m_simulationData.getDensityAdv(index);
	const Vector3r &xi = m_model->getPosition(0, index);
	const Vector3r &vi = m_model->getVelocity(0, index);
	densityAdv = 0.0;
	unsigned int numNeighbors = m_model->numberOfNeighbors(0, index);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int j = 0; j < numNeighbors; j++)
	{
		const unsigned int neighborIndex = m_model->getNeighbor(0, index, j);
		const Vector3r &xj = m_model->getPosition(0, neighborIndex);
		const Vector3r &vj = m_model->getVelocity(0, neighborIndex);
		densityAdv += m_model->getMass(neighborIndex) * (vi - vj).dot(m_model->gradW(xi - xj));
	}

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
	{
		numNeighbors += m_model->numberOfNeighbors(pid, index);
		for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, index); j++)
		{
			const unsigned int neighborIndex = m_model->getNeighbor(pid, index, j);
			const Vector3r &xj = m_model->getPosition(pid, neighborIndex);
			const Vector3r &vj = m_model->getVelocity(pid, neighborIndex);
			densityAdv += m_model->getBoundaryPsi(pid, neighborIndex) * (vi - vj).dot(m_model->gradW(xi - xj));
		}
	}

	// only correct positive divergence
	densityAdv = max(densityAdv, 0.0);

	// in case of particle deficiency do not perform a divergence solve
	if (numNeighbors < 20)
		densityAdv = 0.0;
}

void TimeStepDFSPH::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
	m_iterationsV = 0;
}

void TimeStepDFSPH::performNeighborhoodSearch()
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

void TimeStepDFSPH::emittedParticles(const unsigned int startIndex)
{
	m_simulationData.emittedParticles(startIndex);
	TimeStep::emittedParticles(startIndex);
}
