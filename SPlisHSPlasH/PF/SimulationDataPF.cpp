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
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_old_position.resize(nModels);
	m_num_fluid_neighbors.resize(nModels);
	m_s.resize(nModels);
	m_mat_diag.resize(nModels);
	m_particleOffset.resize(nModels);

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);

		m_old_position[i].resize(fm->numParticles(), Vector3r::Zero());
		m_num_fluid_neighbors[i].resize(fm->numParticles(), 0);
		m_s[i].resize(fm->numParticles(), Vector3r::Zero());
		m_mat_diag[i].resize(fm->numParticles(), Vector3r::Zero());
		
		m_particleOffset[i] = (i==0) ? 0 : sim->getFluidModel(i - 1)->numActiveParticles() + m_particleOffset[i - 1];
	}
}

void SimulationDataPF::cleanup()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_old_position[i].clear();
		m_num_fluid_neighbors[i].clear();
		m_s[i].clear();
		m_mat_diag[i].clear();
	}
	m_old_position.clear();
	m_num_fluid_neighbors.clear();
	m_s.clear();
	m_mat_diag.clear();
	m_particleOffset.clear();
}

void SimulationDataPF::reset()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numParticles(); j++)
		{
			m_num_fluid_neighbors[i][j] = 0;
			m_old_position[i][j].setZero();
			m_s[i][j].setZero();
			m_mat_diag[i][j].setZero();
		}
		m_particleOffset[i] = (i == 0) ? 0 : sim->getFluidModel(i - 1)->numActiveParticles() + m_particleOffset[i - 1];
	}
}

void SimulationDataPF::performNeighborhoodSearchSort()
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
			d.sort_field(&m_old_position[i][0]);
			d.sort_field(&m_num_fluid_neighbors[i][0]);
			d.sort_field(&m_s[i][0]);
			d.sort_field(&m_mat_diag[i][0]);
		}
	}
}

void SimulationDataPF::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	// initialize lastX values for new particles
	Simulation *sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_old_position[fluidModelIndex][j] = model->getPosition(j);
		m_num_fluid_neighbors[fluidModelIndex][j] = 0;
		m_s[fluidModelIndex][j].setZero();
		m_mat_diag[fluidModelIndex][j].setZero();
	}
	// increase particle offsets of following fluid models
	const unsigned int numEmittedParticles = model->numActiveParticles() - startIndex;
	for (unsigned int j = fluidModelIndex + 1; j < sim->numberOfFluidModels(); j++)
	{
		m_particleOffset[j] += numEmittedParticles;
	}
}
