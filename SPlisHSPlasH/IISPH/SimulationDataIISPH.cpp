#include "SimulationDataIISPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;

SimulationDataIISPH::SimulationDataIISPH()
{
}

SimulationDataIISPH::~SimulationDataIISPH(void)
{
	cleanup();
}

void SimulationDataIISPH::init()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	m_aii.resize(model->numParticles(), 0.0);
	m_dii.resize(model->numParticles(), Vector3r::Zero());
	m_dij_pj.resize(model->numParticles(), Vector3r::Zero());
	m_density_adv.resize(model->numParticles(), 0.0);
	m_pressure.resize(model->numParticles(), 0.0);
	m_lastPressure.resize(model->numParticles(), 0.0);
	m_pressureAccel.resize(model->numParticles(), Vector3r::Zero());
}

void SimulationDataIISPH::cleanup()
{
	m_aii.clear();
	m_dii.clear();
	m_dij_pj.clear();
	m_density_adv.clear();
	m_pressure.clear();
	m_lastPressure.clear();
	m_pressureAccel.clear();
}

void SimulationDataIISPH::reset()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	for(unsigned int i=0; i < model->numActiveParticles(); i++)
	{
		m_lastPressure[i] = 0.0;
	}
}

void SimulationDataIISPH::performNeighborhoodSearchSort()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numPart = model->numActiveParticles();
	if (numPart == 0)
		return;

	auto const& d = model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_aii[0]);
	d.sort_field(&m_dii[0]);
	d.sort_field(&m_dij_pj[0]);
	d.sort_field(&m_density_adv[0]);
	d.sort_field(&m_pressure[0]);
	d.sort_field(&m_lastPressure[0]);
	d.sort_field(&m_pressureAccel[0]);
}

void SimulationDataIISPH::emittedParticles(const unsigned int startIndex)
{
	// initialize last pressure values for new particles
	FluidModel *model = Simulation::getCurrent()->getModel();
	for (unsigned int i = startIndex; i < model->numActiveParticles(); i++)
	{
		m_lastPressure[i] = 0.0;
	}
}
