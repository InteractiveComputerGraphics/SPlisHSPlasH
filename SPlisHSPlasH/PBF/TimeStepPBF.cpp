#include "TimeStepPBF.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "TimeIntegration.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataPBF.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

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

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "lambda", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getLambda(fluidModelIndex, i); } });
		model->addField({ "deltaX", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDeltaX(fluidModelIndex, i)[0]; } });
	}
}

TimeStepPBF::~TimeStepPBF(void)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("lambda");
		model->removeFieldByName("deltaX");
	}
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
	const unsigned int nModels = sim->numberOfFluidModels();
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		clearAccelerations(fluidModelIndex);

		// Time integration
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int) model->numActiveParticles(); i++)
			{
				m_simulationData.getLastPosition(fluidModelIndex, i) = m_simulationData.getOldPosition(fluidModelIndex, i);
				m_simulationData.getOldPosition(fluidModelIndex, i) = model->getPosition(i);
				if (model->getParticleState(i) == ParticleState::Active)
					TimeIntegration::semiImplicitEuler(h, model->getMass(i), model->getPosition(i), model->getVelocity(i), model->getAcceleration(i));
			}
		}
	}

	// Perform neighborhood search
	performNeighborhoodSearch();

	// Solve density constraint
	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;

	// Update velocities	
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		if (m_velocityUpdateMethod == ENUM_PBF_FIRST_ORDER)
		{
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static)  
				for (int i = 0; i < (int)model->numActiveParticles(); i++)
				{
					if (model->getParticleState(i) == ParticleState::Active)
						TimeIntegration::velocityUpdateFirstOrder(h, model->getMass(i), model->getPosition(i), m_simulationData.getOldPosition(fluidModelIndex, i), model->getVelocity(i));
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
					if (model->getParticleState(i) == ParticleState::Active)
						TimeIntegration::velocityUpdateSecondOrder(h, model->getMass(i), model->getPosition(i), m_simulationData.getOldPosition(fluidModelIndex, i), m_simulationData.getLastPosition(fluidModelIndex, i), model->getVelocity(i));
					model->getAcceleration(i).setZero();
				}
			}
		}
	}

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		computeDensities(fluidModelIndex);
	sim->computeNonPressureForces();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				if (model->getParticleState(i) == ParticleState::Active)
					model->getVelocity(i) += h * model->getAcceleration(i);
			}
		}
	}

	sim->updateTimeStepSize();

	sim->emitParticles();
	sim->animateParticles();

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
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

	Real avg_density_err = 0;
	m_iterations = 0;
	bool chk = false;

	while ((!chk || (m_iterations < m_minIterations)) && (m_iterations < m_maxIterations))
	{
		chk = true;
		for (unsigned int i = 0; i < nFluids; i++)
		{
			FluidModel *model = sim->getFluidModel(i);
			const Real density0 = model->getDensity0();

			avg_density_err = 0.0;
			pressureSolveIteration(i, avg_density_err);

			// Maximal allowed density fluctuation
			const Real eta = m_maxError * static_cast<Real>(0.01) * density0;  // maxError is given in percent
			chk = chk && (avg_density_err <= eta);
		}
		m_iterations++;
	}
	INCREASE_COUNTER("PBF - iterations", m_iterations);
}

void TimeStepPBF::pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	if (numParticles == 0)
		return;

	const Real invH = static_cast<Real>(1.0) / TimeManager::getCurrent()->getTimeStepSize();
	const Real invH2 = invH*invH;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real eps = 1.0e-6;

	const Real density0 = model->getDensity0();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numParticles; i++)
		{
			Real &density = model->getDensity(i);				
			m_simulationData.getLambda(fluidModelIndex, i) = 0.0;

			// Compute current density for particle i
			density = model->getVolume(i) * sim->W_zero();
			const Vector3r &xi = model->getPosition(i);

			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				computeVolumeAndBoundaryX(fluidModelIndex, i, xi);
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				computeDensityAndGradient(fluidModelIndex, i, xi);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				density += fm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					density += bm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					density += rho;
				);
			}	
			else // Bender2019
			{
				forall_volume_maps(
					density += Vj * sim->W(xi - xj);
				);
			}

			const Real density_err = density0 * (max(density, static_cast<Real>(1.0)) - static_cast<Real>(1.0));
			#pragma omp atomic
			avg_density_err += density_err;

			// Evaluate constraint function
			const Real C = std::max(density - static_cast<Real>(1.0), static_cast<Real>(0.0));			// clamp to prevent particle clumping at surface

			if (C != 0.0)
			{
				// Compute gradients dC/dx_j 
				Real sum_grad_C2 = 0.0;
				Vector3r gradC_i(0.0, 0.0, 0.0);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				forall_fluid_neighbors(
					const Vector3r gradC_j = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
					sum_grad_C2 += gradC_j.squaredNorm();
					gradC_i -= gradC_j;
				);

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					forall_boundary_neighbors(
						// Boundary: Akinci2012
						const Vector3r gradC_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
						gradC_i -= gradC_j;
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						gradC_i -= gradRho;
					);
				}
				else   // Bender2019
				{
					forall_volume_maps(
						const Vector3r gradC_j = -Vj * sim->gradW(xi - xj);
						gradC_i -= gradC_j;
					);
				}

				sum_grad_C2 += gradC_i.squaredNorm();

				// Compute lambda
				Real &lambda = m_simulationData.getLambda(fluidModelIndex, i);
				lambda = -C / (sum_grad_C2 + eps);
			}
			else
				m_simulationData.getLambda(fluidModelIndex, i) = 0.0;
		}
	}

	avg_density_err /= numParticles;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &corr = m_simulationData.getDeltaX(fluidModelIndex, i);

			// Compute position correction
			corr.setZero();

			const Vector3r &xi = model->getPosition(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Vector3r gradC_j = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
				corr -= (m_simulationData.getLambda(fluidModelIndex, i) + (fm_neighbor->getDensity0() / density0) * m_simulationData.getLambda(pid, neighborIndex)) * gradC_j;
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					// Boundary: Akinci2012
					const Vector3r gradC_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
					const Vector3r dx = m_simulationData.getLambda(fluidModelIndex, i) * gradC_j;
					corr -= dx;

					bm_neighbor->addForce(xj, model->getMass(i) * dx * invH2);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					const Vector3r dx = m_simulationData.getLambda(fluidModelIndex, i) * gradRho;
					corr -= dx;
					bm_neighbor->addForce(xj, model->getMass(i) * dx * invH2);	
				);
			}
			else   // Bender2019
			{
				forall_volume_maps(
					const Vector3r gradC_j = -Vj * sim->gradW(xi - xj);	
					const Vector3r dx = m_simulationData.getLambda(fluidModelIndex, i) * gradC_j;
					corr -= dx;
					bm_neighbor->addForce(xj, model->getMass(i) * dx * invH2);
				);
			}
		}
	}
		
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Active)
				model->getPosition(i) += m_simulationData.getDeltaX(fluidModelIndex, i);
		}
	}
}

void TimeStepPBF::performNeighborhoodSearch()
{
	if (Simulation::getCurrent()->zSortEnabled())
	{
		if (m_counter % 500 == 0)
		{
			Simulation::getCurrent()->performNeighborhoodSearchSort();
			m_simulationData.performNeighborhoodSearchSort();
		}
		m_counter++;
	}

	Simulation::getCurrent()->performNeighborhoodSearch();
}

void TimeStepPBF::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepPBF::resize()
{
	m_simulationData.init();
}
