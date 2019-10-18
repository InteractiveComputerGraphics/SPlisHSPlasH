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
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_aii.resize(nModels);
	m_dii.resize(nModels);
	m_dij_pj.resize(nModels);
	m_density_adv.resize(nModels);
	m_pressure.resize(nModels);
	m_lastPressure.resize(nModels);
	m_pressureAccel.resize(nModels);
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		m_aii[i].resize(fm->numParticles(), 0.0);
		m_dii[i].resize(fm->numParticles(), Vector3r::Zero());
		m_dij_pj[i].resize(fm->numParticles(), Vector3r::Zero());
		m_density_adv[i].resize(fm->numParticles(), 0.0);
		m_pressure[i].resize(fm->numParticles(), 0.0);
		m_lastPressure[i].resize(fm->numParticles(), 0.0);
		m_pressureAccel[i].resize(fm->numParticles(), Vector3r::Zero());
	}
}

void SimulationDataIISPH::cleanup()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_aii[i].clear();
		m_dii[i].clear();
		m_dij_pj[i].clear();
		m_density_adv[i].clear();
		m_pressure[i].clear();
		m_lastPressure[i].clear();
		m_pressureAccel[i].clear();
	}
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
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numActiveParticles(); j++)
		{
			m_lastPressure[i][j] = 0.0;
		}
	}
}

void SimulationDataIISPH::performNeighborhoodSearchSort()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		const unsigned int numPart = fm->numActiveParticles();
		if (numPart != 0)
		{
			auto const& d = sim->getNeighborhoodSearch()->point_set(fm->getPointSetIndex());
			d.sort_field(&m_aii[i][0]);
			d.sort_field(&m_dii[i][0]);
			d.sort_field(&m_dij_pj[i][0]);
			d.sort_field(&m_density_adv[i][0]);
			d.sort_field(&m_pressure[i][0]);
			d.sort_field(&m_lastPressure[i][0]);
			d.sort_field(&m_pressureAccel[i][0]);
		}
	}
}

void SimulationDataIISPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	// initialize last pressure values for new particles
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_lastPressure[fluidModelIndex][j] = 0.0;
	}
}
