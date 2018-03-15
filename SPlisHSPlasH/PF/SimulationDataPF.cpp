#include "SimulationDataPF.h"

#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;

SimulationDataPF::SimulationDataPF()
{
}

SimulationDataPF::~SimulationDataPF(void)
{
	cleanup();
}


void SimulationDataPF::init()
{
	FluidModel *model = Simulation::getCurrent()->getModel();

	m_old_position       .resize(model->numParticles(), Vector3r::Zero());
	m_num_fluid_neighbors.resize(model->numParticles(), 0);
	m_x                  .resize(model->numParticles(), Vector3r::Zero());
	m_s                  .resize(model->numParticles(), Vector3r::Zero());
	m_mat_diag           .resize(model->numParticles(), Vector3r::Zero());
}

void SimulationDataPF::cleanup()
{
	m_old_position       .clear();
	m_num_fluid_neighbors.clear();
	m_x                  .clear();
	m_s                  .clear();
	m_mat_diag           .clear();
}

void SimulationDataPF::reset()
{
	#pragma omp parallel for
	for (int i = 0; i < (int) m_x.size(); i++)
	{
		m_num_fluid_neighbors[i] = 0;
		m_x[i].setZero();
		m_old_position[i].setZero();
		m_s[i].setZero();
		m_mat_diag[i].setZero();
	}
}

void SimulationDataPF::performNeighborhoodSearchSort()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numPart = model->numActiveParticles();
	if (numPart == 0)
		return;

	auto const& d = model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_old_position[0]);
	d.sort_field(&m_num_fluid_neighbors[0]);
	d.sort_field((Vector3r*)&m_x[0]);
	d.sort_field(&m_s[0]);
	d.sort_field(&m_mat_diag[0]);
}

void SimulationDataPF::emittedParticles(const unsigned int startIndex)
{
	// initialize values for new particles
	FluidModel *model = Simulation::getCurrent()->getModel();
	for (unsigned int i = startIndex; i < model->numActiveParticles(); i++)
	{
		m_old_position[i] = model->getPosition(0, i);
		m_num_fluid_neighbors[i] = 0;
		m_s[i].setZero();
		m_mat_diag[i].setZero();
	}
}
