#include "TimeStepWCSPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataWCSPH.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "../Simulation.h"

using namespace SPH;
using namespace std;
using namespace GenParam;

int TimeStepWCSPH::STIFFNESS = -1;
int TimeStepWCSPH::EXPONENT = -1;


TimeStepWCSPH::TimeStepWCSPH() :
	TimeStep()
{
	m_simulationData.init();
	m_counter = 0;

	m_stiffness = 50.0;
	m_exponent = 7.0;
}

TimeStepWCSPH::~TimeStepWCSPH(void)
{
}

void TimeStepWCSPH::initParameters()
{
	TimeStep::initParameters();

	STIFFNESS = createNumericParameter("stiffness", "Stiffness", &m_stiffness);
	setGroup(STIFFNESS, "WCSPH");
	setDescription(STIFFNESS, "Stiffness coefficient of EOS.");
	static_cast<RealParameter*>(getParameter(STIFFNESS))->setMinValue(1e-6);

	EXPONENT = createNumericParameter("exponent", "Exponent (gamma)", &m_exponent);
	setGroup(EXPONENT, "WCSPH");
	setDescription(EXPONENT, "Exponent of EOS.");
	static_cast<RealParameter*>(getParameter(EXPONENT))->setMinValue(1e-6);

}

void TimeStepWCSPH::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	performNeighborhoodSearch();

	// Compute accelerations: a(t)
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		clearAccelerations(fluidModelIndex);
		computeDensities(fluidModelIndex);
	}
	sim->computeNonPressureForces();


	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const Real density0 = model->getValue<Real>(FluidModel::DENSITY0);
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				Real &density = model->getDensity(i);
				density = max(density, density0) / density0;
				m_simulationData.getPressure(fluidModelIndex, i) = m_stiffness * (pow(density / 1.0, m_exponent) - 1.0);
			}
		}

		computePressureAccels(fluidModelIndex);
	}

	sim->updateTimeStepSize();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static) 
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				Vector3r &pos = model->getPosition(i);
				Vector3r &vel = model->getVelocity(i);
				Vector3r &accel = model->getAcceleration(i);
				accel += m_simulationData.getPressureAccel(fluidModelIndex, i);
				vel += accel * h;
				pos += vel * h;
			}
		}
	}

	sim->emitParticles();

	// Compute new time	
	tm->setTime (tm->getTime () + h);
}


void TimeStepWCSPH::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
}

void TimeStepWCSPH::computePressureAccels(const unsigned int fluidModelIndex)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();

	// Compute pressure forces
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = model->getPosition(i);
			const Real &density_i = model->getDensity(i);

			Vector3r &ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
			ai.setZero();

			const Real dpi = m_simulationData.getPressure(fluidModelIndex, i) / (density_i*density_i);
			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				// Pressure 
				const Real &density_j = fm_neighbor->getDensity(neighborIndex);
				const Real dpj = m_simulationData.getPressure(pid, neighborIndex) / (density_j*density_j);
				ai -= fm_neighbor->getVolume(neighborIndex) * (dpi + (fm_neighbor->getDensity0() / density0) * dpj) * sim->gradW(xi - xj);
			)

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			forall_boundary_neighbors(
				const Vector3r a = bm_neighbor->getVolume(neighborIndex) * (dpi)* sim->gradW(xi - xj);
				ai -= a;
				bm_neighbor->getForce(neighborIndex) += model->getMass(i) * a;
			)
		}
	}
}

void TimeStepWCSPH::performNeighborhoodSearch()
{
	if (m_counter % 500 == 0)
	{
		Simulation::getCurrent()->performNeighborhoodSearchSort();
		m_simulationData.performNeighborhoodSearchSort();
	}
	m_counter++;

	Simulation::getCurrent()->performNeighborhoodSearch();
}

void TimeStepWCSPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepWCSPH::resize()
{
	m_simulationData.init();
}

