#include "SimulationDataDFSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;

SimulationDataDFSPH::SimulationDataDFSPH() :
	m_factor(),
	m_kappa(),
	m_kappaV(),
	m_density_adv()
{
}

SimulationDataDFSPH::~SimulationDataDFSPH(void)
{
	cleanup();
}


void SimulationDataDFSPH::init()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_factor.resize(nModels);
	m_kappa.resize(nModels);
	m_kappaV.resize(nModels);
	m_density_adv.resize(nModels);
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		m_factor[i].resize(fm->numParticles(), 0.0);
		m_kappa[i].resize(fm->numParticles(), 0.0);
		m_kappaV[i].resize(fm->numParticles(), 0.0);
		m_density_adv[i].resize(fm->numParticles(), 0.0);
	}
}

void SimulationDataDFSPH::cleanup()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_factor[i].clear();
		m_kappa[i].clear();
		m_kappaV[i].clear();
		m_density_adv[i].clear();
	}
	m_factor.clear();
	m_kappa.clear();
	m_kappaV.clear();
	m_density_adv.clear();
}

void SimulationDataDFSPH::reset()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numActiveParticles(); j++)
		{
			m_kappa[i][j] = 0.0;
			m_kappaV[i][j] = 0.0;
		}
	}
}

void SimulationDataDFSPH::performNeighborhoodSearchSort()
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
			d.sort_field(&m_factor[i][0]);
			d.sort_field(&m_kappa[i][0]);
			d.sort_field(&m_kappaV[i][0]);
			d.sort_field(&m_density_adv[i][0]);
		}
	}
}

void SimulationDataDFSPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	// initialize kappa values for new particles
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_kappa[fluidModelIndex][j] = 0.0;
		m_kappaV[fluidModelIndex][j] = 0.0;
	}
}
