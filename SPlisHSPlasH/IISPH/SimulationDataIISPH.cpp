#include "SimulationDataIISPH.h"
#include "SPlisHSPlasH/SPHKernels.h"

using namespace SPH;

SimulationDataIISPH::SimulationDataIISPH()
{
	m_model = NULL;
}

SimulationDataIISPH::~SimulationDataIISPH(void)
{
	cleanup();
}

void SimulationDataIISPH::init(FluidModel *model)
{
	m_model = model;

	m_aii.resize(model->numParticles(), 0.0);
	m_dii.resize(model->numParticles(), SPH::Vector3r::Zero());
	m_dij_pj.resize(model->numParticles(), SPH::Vector3r::Zero());
	m_density_adv.resize(model->numParticles(), 0.0);
	m_pressure.resize(model->numParticles(), 0.0);
	m_lastPressure.resize(model->numParticles(), 0.0);
	m_pressureAccel.resize(model->numParticles(), SPH::Vector3r::Zero());
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
	for(unsigned int i=0; i < m_model->numActiveParticles(); i++)
	{
		m_lastPressure[i] = 0.0;
	}
}

void SimulationDataIISPH::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	auto const& d = m_model->getNeighborhoodSearch()->point_set(0);
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
	for (unsigned int i = startIndex; i < m_model->numActiveParticles(); i++)
	{
		m_lastPressure[i] = 0.0;
	}
}
