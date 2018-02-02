#include "TimeStepPBF.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "TimeIntegration.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataPBF.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;
using namespace std;
using namespace GenParam;

int TimeStepPBF::VELOCITY_UPDATE_METHOD = -1;
int TimeStepPBF::ENUM_PBF_FIRST_ORDER = -1;
int TimeStepPBF::ENUM_PBF_SECOND_ORDER = -1;


TimeStepPBF::TimeStepPBF() :
	TimeStep()
{
	m_simulationData.init();
	m_counter = 0;
	m_velocityUpdateMethod = 0;
}

TimeStepPBF::~TimeStepPBF(void)
{
}

void TimeStepPBF::initParameters()
{
	TimeStep::initParameters();

	VELOCITY_UPDATE_METHOD = createEnumParameter("velocityUpdateMethod", "Velocity update method", &m_velocityUpdateMethod);
	setGroup(VELOCITY_UPDATE_METHOD, "PBF");
	setDescription(VELOCITY_UPDATE_METHOD, "Method for the velocity integration.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(VELOCITY_UPDATE_METHOD));
	enumParam->addEnumValue("First Order Update", ENUM_PBF_FIRST_ORDER);
	enumParam->addEnumValue("Second Order Update", ENUM_PBF_SECOND_ORDER);
}

void TimeStepPBF::step()
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getModel();
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	clearAccelerations();

	// Time integration
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) model->numActiveParticles(); i++)
		{
			m_simulationData.getLastPosition(i) = m_simulationData.getOldPosition(i);
			m_simulationData.getOldPosition(i) = model->getPosition(0, i);
			TimeIntegration::semiImplicitEuler(h, model->getMass(i), model->getPosition(0, i), model->getVelocity(0, i), model->getAcceleration(i));
		}
	}

	// Perform neighborhood search
	performNeighborhoodSearch();

	// Solve density constraint
	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;

	// Update velocities	
	if (m_velocityUpdateMethod == ENUM_PBF_FIRST_ORDER)
	{
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int) model->numActiveParticles(); i++)
			{
				TimeIntegration::velocityUpdateFirstOrder(h, model->getMass(i), model->getPosition(0, i), m_simulationData.getOldPosition(i), model->getVelocity(0, i));
				model->getAcceleration(i).setZero();
			}
		}
	}
	else
	{
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				TimeIntegration::velocityUpdateSecondOrder(h, model->getMass(i), model->getPosition(0, i), m_simulationData.getOldPosition(i), m_simulationData.getLastPosition(i), model->getVelocity(0, i));
				model->getAcceleration(i).setZero();
			}
		}
	}

	sim->computeNonPressureForces();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)model->numActiveParticles(); i++)
		{
			model->getVelocity(0, i) += h * model->getAcceleration(i);
		}
	}

	sim->updateTimeStepSize();

	sim->emitParticles();

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

	FluidModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numParticles = model->numActiveParticles();
	const Real invH = 1.0 / TimeManager::getCurrent()->getTimeStepSize();
	const Real invH2 = invH*invH;

	const Real density0 = model->getValue<Real>(FluidModel::DENSITY0);

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
				Real &density = model->getDensity(i);				

				// Compute current density for particle i
				density = model->getMass(i) * model->W_zero();
				const Vector3r &xi = model->getPosition(0, i);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = model->getNeighbor(0, i, j);
					const Vector3r &xj = model->getPosition(0, neighborIndex);
					density += model->getMass(neighborIndex) * model->W(xi - xj);
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
						// Boundary: Akinci2012
						density += model->getBoundaryPsi(pid, neighborIndex) * model->W(xi - xj);
					}
				}

				const Real density_err = max(density, density0) - density0;
				#pragma omp atomic
				avg_density_err += density_err;

				// Evaluate constraint function
				const Real C = std::max(density / density0 - 1.0, 0.0);			// clamp to prevent particle clumping at surface

				if (C != 0.0)
				{
					// Compute gradients dC/dx_j 
					Real sum_grad_C2 = 0.0;
					Vector3r gradC_i(0.0, 0.0, 0.0);

					//////////////////////////////////////////////////////////////////////////
					// Fluid
					//////////////////////////////////////////////////////////////////////////
					for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
					{
						const unsigned int neighborIndex = model->getNeighbor(0, i, j);
						const Vector3r &xj = model->getPosition(0, neighborIndex);
						const Vector3r gradC_j = -model->getMass(neighborIndex) / density0 * model->gradW(xi - xj);
						sum_grad_C2 += gradC_j.squaredNorm();
						gradC_i -= gradC_j;
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
							
							// Boundary: Akinci2012
							const Vector3r gradC_j = -model->getBoundaryPsi(pid, neighborIndex) / density0 * model->gradW(xi - xj);
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
		}

		avg_density_err /= numParticles;

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				Vector3r &corr = m_simulationData.getDeltaX(i);

				// Compute position correction
				corr.setZero();
				const Vector3r &xi = model->getPosition(0, i);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
				{
					const unsigned int neighborIndex = model->getNeighbor(0, i, j);
					const Vector3r &xj = model->getPosition(0, neighborIndex);
					const Vector3r gradC_j = -model->getMass(neighborIndex) / density0 * model->gradW(xi - xj);
					corr -= (m_simulationData.getLambda(i) + m_simulationData.getLambda(neighborIndex)) * gradC_j;
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
						// Boundary: Akinci2012
						const Vector3r gradC_j = -model->getBoundaryPsi(pid, neighborIndex) / density0 * model->gradW(xi - xj);
						const Vector3r dx = 2.0 * m_simulationData.getLambda(i) * gradC_j;
						corr -= dx;

						model->getForce(pid, neighborIndex) += model->getMass(i) * dx * invH2;
					}
				}
			}
		}
		
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				model->getPosition(0, i) += m_simulationData.getDeltaX(i);
			}
		}

		m_iterations++;
	}
}

void TimeStepPBF::performNeighborhoodSearch()
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

void TimeStepPBF::emittedParticles(const unsigned int startIndex)
{

	m_simulationData.emittedParticles(startIndex);
	Simulation::getCurrent()->emittedParticles(startIndex);
}

void TimeStepPBF::resize()
{
	m_simulationData.init();
}
