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

	m_old_position.resize(model->numParticles(), Vector3r::Zero());
	m_num_fluid_neighbors.resize(model->numParticles(), 0);
	m_x.resize(3 * model->numParticles(), 0);
}

void SimulationDataPF::cleanup()
{
	m_old_position.clear();
	m_num_fluid_neighbors.clear();
}

void SimulationDataPF::reset()
{
	for (unsigned int i = 0; i < m_old_position.size(); i++)
	{
		m_old_position[i].setZero();
		m_num_fluid_neighbors[i] = 0;
	}
	for (unsigned int i = 0; i < m_x.size(); i++)
	{
		m_x[i] = 0;
	}
}

void SimulationDataPF::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numParticles();
	if (numPart == 0)
		return;

	auto const& d = m_model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_old_position[0]);
	d.sort_field(&m_num_fluid_neighbors[0]);
}
