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
	FluidModel *model = Simulation::getCurrent()->getModel();

	m_factor.resize(model->numParticles(), 0.0);
	m_kappa.resize(model->numParticles(), 0.0);
	m_kappaV.resize(model->numParticles(), 0.0);
	m_density_adv.resize(model->numParticles(), 0.0);
}

void SimulationDataDFSPH::cleanup()
{
	m_factor.clear();
	m_kappa.clear();
	m_kappaV.clear();
	m_density_adv.clear();
}

void SimulationDataDFSPH::reset()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	for (unsigned int i = 0; i < model->numActiveParticles(); i++)
	{
		m_kappa[i] = 0.0;
		m_kappaV[i] = 0.0;
	}
}

void SimulationDataDFSPH::performNeighborhoodSearchSort()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numPart = model->numActiveParticles();
	if (numPart == 0)
		return;

	auto const& d = model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_factor[0]);
	d.sort_field(&m_kappa[0]);
	d.sort_field(&m_kappaV[0]);
	d.sort_field(&m_density_adv[0]);
}

void SimulationDataDFSPH::emittedParticles(const unsigned int startIndex)
{
	// initialize kappa values for new particles
	FluidModel *model = Simulation::getCurrent()->getModel();
	for (unsigned int i = startIndex; i < model->numActiveParticles(); i++)
	{
		m_kappa[i] = 0.0;
		m_kappaV[i] = 0.0;
	}
}
