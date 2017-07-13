#include "SimulationDataDFSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"

using namespace SPH;

SimulationDataDFSPH::SimulationDataDFSPH() :
	m_factor(),
	m_kappa(),
	m_kappaV(),
	m_density_adv()
{
	m_model = NULL;
}

SimulationDataDFSPH::~SimulationDataDFSPH(void)
{
	cleanup();
}


void SimulationDataDFSPH::init(FluidModel *model)
{
	m_model = model;

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
	for (unsigned int i = 0; i < m_model->numActiveParticles(); i++)
	{
		m_kappa[i] = 0.0;
		m_kappaV[i] = 0.0;
	}
}

void SimulationDataDFSPH::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	auto const& d = m_model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_factor[0]);
	d.sort_field(&m_kappa[0]);
	d.sort_field(&m_kappaV[0]);
	d.sort_field(&m_density_adv[0]);
}

void SimulationDataDFSPH::emittedParticles(const unsigned int startIndex)
{
	// initialize kappa values for new particles
	for (unsigned int i = startIndex; i < m_model->numActiveParticles(); i++)
	{
		m_kappa[i] = 0.0;
		m_kappaV[i] = 0.0;
	}
}
