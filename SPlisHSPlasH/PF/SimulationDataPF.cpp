#include "SimulationDataPF.h"
#include "SPlisHSPlasH/SPHKernels.h"

using namespace SPH;

SimulationDataPF::SimulationDataPF()
{
	m_model = NULL;
}

SimulationDataPF::~SimulationDataPF(void)
{
	cleanup();
}


void SimulationDataPF::init(FluidModel *model)
{
	m_model = model;

	m_factor.resize(model->numParticles(), 0.0);
	m_kappa.resize(model->numParticles(), 0.0);
	m_kappaV.resize(model->numParticles(), 0.0);
	m_density_adv.resize(model->numParticles(), 0.0);
}

void SimulationDataPF::cleanup()
{
	m_factor.clear();
	m_kappa.clear();
	m_kappaV.clear();
	m_density_adv.clear();
}

void SimulationDataPF::reset()
{
	for (unsigned int i = 0; i < m_kappa.size(); i++)
	{
		m_kappa[i] = 0.0;
		m_kappaV[i] = 0.0;
	}
}

void SimulationDataPF::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numParticles();
	if (numPart == 0)
		return;

	auto const& d = m_model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_factor[0]);
	d.sort_field(&m_kappa[0]);
	d.sort_field(&m_kappaV[0]);
	d.sort_field(&m_density_adv[0]);
}
