#include "SimulationDataDFSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;

SimulationDataDFSPH::SimulationDataDFSPH() :
	m_factor(),
	m_pressure_rho2(),
	m_pressure_rho2_V(),
	m_pressureAccel(),
	m_density_adv()
{
}

SimulationDataDFSPH::~SimulationDataDFSPH(void)
{
	cleanup();
}


void SimulationDataDFSPH::init()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_factor.resize(nModels);
	m_density_adv.resize(nModels);
	m_pressure_rho2.resize(nModels);
	m_pressure_rho2_V.resize(nModels);
	m_pressureAccel.resize(nModels);

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		m_factor[i].resize(fm->numParticles(), 0.0);
		m_density_adv[i].resize(fm->numParticles(), 0.0);
		m_pressure_rho2[i].resize(fm->numParticles(), 0.0);
		m_pressure_rho2_V[i].resize(fm->numParticles(), 0.0);
		m_pressureAccel[i].resize(fm->numParticles(), Vector3r::Zero());
	}
}

void SimulationDataDFSPH::cleanup()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_factor[i].clear();
		m_density_adv[i].clear();
		m_pressure_rho2[i].clear();
		m_pressure_rho2_V[i].clear();
		m_pressureAccel[i].clear();
	}
	m_factor.clear();
	m_density_adv.clear();
	m_pressure_rho2.clear();
	m_pressure_rho2_V.clear();
	m_pressureAccel.clear();
}

void SimulationDataDFSPH::reset()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numParticles(); j++)
		{
			m_density_adv[i][j] = 0.0;
			m_pressure_rho2[i][j] = 0.0;
			m_pressure_rho2_V[i][j] = 0.0;
			m_factor[i][j] = 0.0;
			m_pressureAccel[i][j].setZero();
		}
	}
}

void SimulationDataDFSPH::performNeighborhoodSearchSort()
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
			//d.sort_field(&m_factor[i][0]);
			//d.sort_field(&m_density_adv[i][0]);
			d.sort_field(&m_pressure_rho2[i][0]);
			d.sort_field(&m_pressure_rho2_V[i][0]);
		}
	}
}

void SimulationDataDFSPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	// initialize kappa values for new particles
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_pressure_rho2[fluidModelIndex][j] = 0.0;
		m_pressure_rho2_V[fluidModelIndex][j] = 0.0;
	}
}
