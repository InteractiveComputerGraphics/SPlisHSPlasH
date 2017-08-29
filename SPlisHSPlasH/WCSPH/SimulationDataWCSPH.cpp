#include "SimulationDataWCSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"

using namespace SPH;

SimulationDataWCSPH::SimulationDataWCSPH()
{
	m_model = NULL;
}

SimulationDataWCSPH::~SimulationDataWCSPH(void)
{
	cleanup();
}

void SimulationDataWCSPH::init(FluidModel *model)
{
	m_model = model;

	m_pressure.resize(model->numParticles(), 0.0);
	m_pressureAccel.resize(model->numParticles(), SPH::Vector3r::Zero());
}

void SimulationDataWCSPH::cleanup()
{
	m_pressure.clear();
	m_pressureAccel.clear();
}

void SimulationDataWCSPH::reset()
{
}

void SimulationDataWCSPH::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	auto const& d = m_model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_pressure[0]);
	d.sort_field(&m_pressureAccel[0]);
}


void SimulationDataWCSPH::emittedParticles(const unsigned int startIndex)
{
	for (unsigned int i = startIndex; i < m_model->numActiveParticles(); i++)
	{
		m_pressure[i] = 0.0;
		m_pressureAccel[i].setZero();
	}
}