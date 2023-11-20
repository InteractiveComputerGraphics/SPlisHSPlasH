#include "SimulationDataICSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;

SimulationDataICSPH::SimulationDataICSPH()
{
}

SimulationDataICSPH::~SimulationDataICSPH(void)
{
	cleanup();
}

void SimulationDataICSPH::init()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_aii.resize(nModels);
	m_density_adv.resize(nModels);
	m_pressure.resize(nModels);
	m_pressureAccel.resize(nModels);
	m_pressureGradient.resize(nModels);
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		m_aii[i].resize(fm->numParticles(), 0.0);
		m_density_adv[i].resize(fm->numParticles(), 0.0);
		m_pressure[i].resize(fm->numParticles(), 0.0);
		m_pressureAccel[i].resize(fm->numParticles(), Vector3r::Zero());
		m_pressureGradient[i].resize(fm->numParticles(), Vector3r::Zero());
	}
}

void SimulationDataICSPH::cleanup()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_aii[i].clear();
		m_density_adv[i].clear();
		m_pressure[i].clear();
		m_pressureAccel[i].clear();
		m_pressureGradient[i].clear();
	}
	m_aii.clear();
	m_density_adv.clear();
	m_pressure.clear();
	m_pressureAccel.clear();
	m_pressureGradient.clear();
}

void SimulationDataICSPH::reset()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numParticles(); j++)
		{
			m_pressure[i][j] = 0.0;
		}
	}
}

void SimulationDataICSPH::performNeighborhoodSearchSort()
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
			d.sort_field(&m_pressure[i][0]);
		}
	}
}

void SimulationDataICSPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	// initialize last pressure values for new particles
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_pressure[fluidModelIndex][j] = 0.0;
	}
}
