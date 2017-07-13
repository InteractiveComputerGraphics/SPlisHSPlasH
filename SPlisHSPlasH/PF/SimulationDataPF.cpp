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

	m_old_position       .resize(    model->numParticles(), Vector3r::Zero());
	m_num_fluid_neighbors.resize(    model->numParticles(), 0);
	m_x                  .resize(3 * model->numParticles(), 0);
	m_s                  .resize(    model->numParticles(), Vector3r::Zero());
}

void SimulationDataPF::cleanup()
{
	m_old_position       .clear();
	m_num_fluid_neighbors.clear();
	m_x                  .clear();
	m_s                  .clear();
}

void SimulationDataPF::reset()
{
	for (unsigned int i = 0; i < m_old_position.size(); i++)
	{
		m_old_position[i].setZero();
		m_num_fluid_neighbors[i] = 0;
		m_s[i].setZero();
	}
	for (unsigned int i = 0; i < m_x.size(); i++)
	{
		m_x[i] = 0;
	}
}

void SimulationDataPF::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	auto const& d = m_model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_old_position[0]);
	d.sort_field(&m_num_fluid_neighbors[0]);
	d.sort_field((Vector3r*)&m_x[0]);
	d.sort_field(&m_s[0]);
}

void SimulationDataPF::emittedParticles(const unsigned int startIndex)
{
	// initialize values for new particles
	for (unsigned int i = startIndex; i < m_model->numActiveParticles(); i++)
	{
		m_old_position[i] = m_model->getPosition(0, i);
		m_num_fluid_neighbors[i] = 0;
		m_s[i].setZero();
	}
}
