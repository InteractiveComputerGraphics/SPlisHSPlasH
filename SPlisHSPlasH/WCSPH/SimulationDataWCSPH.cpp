#include "SimulationDataWCSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "../Simulation.h"

using namespace SPH;

SimulationDataWCSPH::SimulationDataWCSPH() : 
	m_pressure(), 
	m_pressureAccel()
{
}

SimulationDataWCSPH::~SimulationDataWCSPH(void)
{
	cleanup();
}

void SimulationDataWCSPH::init()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_pressure.resize(nModels);
	m_pressureAccel.resize(nModels);
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		m_pressure[i].resize(fm->numParticles(), 0.0);
		m_pressureAccel[i].resize(fm->numParticles(), Vector3r::Zero());
	}
}

void SimulationDataWCSPH::cleanup()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_pressure[i].clear();
		m_pressureAccel[i].clear();
	}
	m_pressure.clear();
	m_pressureAccel.clear();
}

void SimulationDataWCSPH::reset()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numParticles(); j++)
		{
			m_pressure[i][j] = 0.0;
			m_pressureAccel[i][j].setZero();
		}
	}
}

void SimulationDataWCSPH::performNeighborhoodSearchSort()
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
			d.sort_field(&m_pressureAccel[i][0]);
		}
	}
}


void SimulationDataWCSPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	// initialize kappa values for new particles
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_pressure[fluidModelIndex][j] = 0.0;
		m_pressureAccel[fluidModelIndex][j].setZero();
	}
}