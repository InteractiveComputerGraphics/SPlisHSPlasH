#include "SimulationDataPBF.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;

SimulationDataPBF::SimulationDataPBF() :
	m_deltaX(), 
	m_lambda(), 
	m_lastX(),
	m_oldX()
{
}

SimulationDataPBF::~SimulationDataPBF(void)
{
	cleanup();
}

void SimulationDataPBF::init()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_lambda.resize(nModels);
	m_deltaX.resize(nModels);
	m_oldX.resize(nModels);
	m_lastX.resize(nModels);
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		m_lambda[i].resize(fm->numParticles(), 0.0);
		m_deltaX[i].resize(fm->numParticles(), Vector3r::Zero());
		m_oldX[i].resize(fm->numParticles(), Vector3r::Zero());
		m_lastX[i].resize(fm->numParticles(), Vector3r::Zero());
	}
	reset();
}

void SimulationDataPBF::cleanup()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_lambda[i].clear();
		m_deltaX[i].clear();
		m_oldX[i].clear();
		m_lastX[i].clear();
	}
	m_lambda.clear();
	m_deltaX.clear();
	m_oldX.clear();
	m_lastX.clear();
}

void SimulationDataPBF::reset()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numActiveParticles(); j++)
		{
			m_deltaX[i][j].setZero();
			m_lambda[i][j] = 0.0;
			getLastPosition(i, j) = fm->getPosition(j);
			getOldPosition(i, j) = fm->getPosition(j);
		}
	}
}

void SimulationDataPBF::performNeighborhoodSearchSort()
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
			d.sort_field(&m_lambda[i][0]);
			d.sort_field(&m_deltaX[i][0]);
			d.sort_field(&m_oldX[i][0]);
			d.sort_field(&m_lastX[i][0]);
		}
	}
}

void SimulationDataPBF::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	// initialize lastX values for new particles
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_lastX[fluidModelIndex][j] = model->getPosition(j);
	}
}

