#include "SimulationDataPCISPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"
#include <iostream>
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Logger.h"

using namespace SPH;

SimulationDataPCISPH::SimulationDataPCISPH()
{
}

SimulationDataPCISPH::~SimulationDataPCISPH(void)
{
	cleanup();
}

void SimulationDataPCISPH::init()
{
	FluidModel *model = Simulation::getCurrent()->getModel();

	m_lastX.resize(model->numParticles(), Vector3r::Zero());
	m_lastV.resize(model->numParticles(), Vector3r::Zero());
	m_pressureAccel.resize(model->numParticles(), Vector3r::Zero());
	m_densityAdv.resize(model->numParticles(), 0.0);
	m_pressure.resize(model->numParticles(), 0.0);
	m_pressureAccel.resize(model->numParticles(), Vector3r::Zero());

	LOG_INFO << "Initialize PCISPH scaling factor";
	m_pcisph_factor = 0.0;
	model->getNeighborhoodSearch()->find_neighbors();

	// Find prototype particle
	// => particle with max. fluid neighbors
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real h2 = h*h;
	const Real density0 = model->getValue<Real>(FluidModel::DENSITY0);
	unsigned int index = 0;
	unsigned int maxNeighbors = 0;

	for (int i = 0; i < (int)model->numActiveParticles(); i++)
	{
		if (model->numberOfNeighbors(0, i) > maxNeighbors)
		{
			maxNeighbors = model->numberOfNeighbors(0, i);
			index = i;
		}
	}

	Vector3r sumGradW = Vector3r::Zero();
	Real sumGradW2 = 0.0;
	const Vector3r &xi = model->getPosition(0, index);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int j = 0; j < model->numberOfNeighbors(0, index); j++)
	{
		const unsigned int neighborIndex = model->getNeighbor(0, index, j);
		const Vector3r &xj = model->getPosition(0, neighborIndex);
		const Vector3r gradW = model->gradW(xi - xj);
		sumGradW += gradW;
		sumGradW2 += gradW.squaredNorm();
	}

	const Real beta = 2.0 * model->getMass(index)*model->getMass(index) / (density0*density0);	// h^2 is multiplied in each iteration
																								// to make the factor independent of h
	m_pcisph_factor = 1.0 / (beta * (sumGradW.squaredNorm() + sumGradW2));
}

void SimulationDataPCISPH::cleanup()
{
	m_lastX.clear();
	m_lastV.clear();
	m_densityAdv.clear();
	m_pressure.clear();
	m_pressureAccel.clear();
}

void SimulationDataPCISPH::reset()
{
}

void SimulationDataPCISPH::performNeighborhoodSearchSort()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numPart = model->numActiveParticles();
	if (numPart == 0)
		return;

	auto const& d = model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_lastX[0]);
	d.sort_field(&m_lastV[0]);
	d.sort_field(&m_densityAdv[0]);
	d.sort_field(&m_pressure[0]);
	d.sort_field(&m_pressureAccel[0]);
}

void SimulationDataPCISPH::emittedParticles(const unsigned int startIndex)
{
	// initialize values for new particles
	FluidModel *model = Simulation::getCurrent()->getModel();
	for (unsigned int i = startIndex; i < model->numActiveParticles(); i++)
	{
		m_lastX[i] = model->getPosition(0, i);
		m_lastV[i] = model->getVelocity(0, i);
	}
}