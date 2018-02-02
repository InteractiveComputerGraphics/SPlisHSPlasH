#include "SimulationDataPBF.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;

SimulationDataPBF::SimulationDataPBF()
{
}

SimulationDataPBF::~SimulationDataPBF(void)
{
	cleanup();
}

void SimulationDataPBF::init()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	m_lambda.resize(model->numParticles(), 0.0);
	m_deltaX.resize(model->numParticles(), Vector3r::Zero());
	m_oldX.resize(model->numParticles(), Vector3r::Zero());
	m_lastX.resize(model->numParticles(), Vector3r::Zero());
	reset();
}

void SimulationDataPBF::cleanup()
{
	m_lambda.clear();
	m_deltaX.clear();
	m_oldX.clear();
	m_lastX.clear();
}

void SimulationDataPBF::reset()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	for(unsigned int i=0; i < model->numActiveParticles(); i++)
	{
		m_deltaX[i].setZero();
		m_lambda[i] = 0.0;
		getLastPosition(i) = model->getPosition(0, i);
		getOldPosition(i) = model->getPosition(0, i);
	}
}

void SimulationDataPBF::performNeighborhoodSearchSort()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numPart = model->numActiveParticles();
	if (numPart == 0)
		return;

	auto const& d = model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_lambda[0]);
	d.sort_field(&m_deltaX[0]);
	d.sort_field(&m_oldX[0]);
	d.sort_field(&m_lastX[0]);
}

void SimulationDataPBF::emittedParticles(const unsigned int startIndex)
{
	// initialize lastX values for new particles
	FluidModel *model = Simulation::getCurrent()->getModel();
	for (unsigned int i = startIndex; i < model->numActiveParticles(); i++)
	{
		m_lastX[i] = model->getPosition(0, i);
	}
}

