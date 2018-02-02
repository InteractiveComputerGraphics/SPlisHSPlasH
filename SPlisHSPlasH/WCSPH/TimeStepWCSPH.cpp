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

	m_stiffness = 50000.0;
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
	FluidModel *model = sim->getModel();
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	performNeighborhoodSearch();

	// Compute accelerations: a(t)
	clearAccelerations();
	computeDensities();
	sim->computeNonPressureForces();

	const Real density0 = model->getValue<Real>(FluidModel::DENSITY0);

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) model->numActiveParticles(); i++)
		{
			Real &density = model->getDensity(i);
			density = max(density, density0);
			m_simulationData.getPressure(i) = m_stiffness * (pow(density/density0, m_exponent) - 1.0);
		}
	}

	computePressureAccels();

	sim->updateTimeStepSize();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)model->numActiveParticles(); i++)
		{
			Vector3r &pos = model->getPosition(0, i);
			Vector3r &vel = model->getVelocity(0, i);
			Vector3r &accel = model->getAcceleration(i);
			accel += m_simulationData.getPressureAccel(i);
			vel += accel * h;
			pos += vel * h;
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

void TimeStepWCSPH::computePressureAccels()
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getModel();
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
#
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

void TimeStepWCSPH::performNeighborhoodSearch()
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getModel();
	const unsigned int numParticles = model->numActiveParticles();

	if (m_counter % 500 == 0)
	{
		model->performNeighborhoodSearchSort();
		m_simulationData.performNeighborhoodSearchSort();
		sim->performNeighborhoodSearchSort();
	}
	m_counter++;

	sim->performNeighborhoodSearch();
}

void TimeStepWCSPH::emittedParticles(const unsigned int startIndex)
{
	m_simulationData.emittedParticles(startIndex);
	Simulation::getCurrent()->emittedParticles(startIndex);
}

void TimeStepWCSPH::resize()
{
	m_simulationData.init();
}

