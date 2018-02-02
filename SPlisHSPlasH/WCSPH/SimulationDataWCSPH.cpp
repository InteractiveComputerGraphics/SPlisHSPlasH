#include "SimulationDataWCSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "../Simulation.h"

using namespace SPH;

SimulationDataWCSPH::SimulationDataWCSPH()
{
}

SimulationDataWCSPH::~SimulationDataWCSPH(void)
{
	cleanup();
}

void SimulationDataWCSPH::init()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	m_pressure.resize(model->numParticles(), 0.0);
	m_pressureAccel.resize(model->numParticles(), Vector3r::Zero());
}

void SimulationDataWCSPH::cleanup()
{
	m_pressure.clear();
	m_pressureAccel.clear();
}

void SimulationDataWCSPH::reset()
{
}

void SimulationDataWCSPH::performNeighborhoodSearchSort()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numPart = model->numActiveParticles();
	if (numPart == 0)
		return;

	auto const& d = model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_pressure[0]);
	d.sort_field(&m_pressureAccel[0]);
}


void SimulationDataWCSPH::emittedParticles(const unsigned int startIndex)
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	for (unsigned int i = startIndex; i < model->numActiveParticles(); i++)
	{
		m_pressure[i] = 0.0;
		m_pressureAccel[i].setZero();
	}
}