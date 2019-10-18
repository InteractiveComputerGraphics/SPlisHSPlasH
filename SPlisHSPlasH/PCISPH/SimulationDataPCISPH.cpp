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
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_lastX.resize(nModels);
	m_lastV.resize(nModels);
	m_densityAdv.resize(nModels);
	m_pressure.resize(nModels);
	m_pressureAccel.resize(nModels);
	m_pcisph_factor.resize(nModels);
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		m_lastX[i].resize(fm->numParticles(), Vector3r::Zero());
		m_lastV[i].resize(fm->numParticles(), Vector3r::Zero());
		m_densityAdv[i].resize(fm->numParticles(), 0.0);
		m_pressure[i].resize(fm->numParticles(), 0.0);
		m_pressureAccel[i].resize(fm->numParticles(), Vector3r::Zero());
	}

	LOG_INFO << "Initialize PCISPH scaling factor";
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		m_pcisph_factor[fluidModelIndex] = 0.0;

		// Find prototype particle
		// => particle with max. fluid neighbors
		const Real density0 = model->getDensity0();
 
 		Vector3r sumGradW = Vector3r::Zero();
 		Real sumGradW2 = 0.0;
 		const Real supportRadius = sim->getSupportRadius();
 		const Real particleRadius = sim->getParticleRadius();
 		const Real diam = 2.0 * particleRadius;
 		const Vector3r xi(0,0,0);
 
 		// use a regular sampling around (0,0,0)
		if (sim->is2DSimulation())
		{
			Vector3r xj = { -supportRadius, -supportRadius, 0.0 };
			while (xj[0] <= supportRadius)
			{
				while (xj[1] <= supportRadius)
				{
					// check if xj is in the support of xi
					if ((xi - xj).squaredNorm() < supportRadius*supportRadius)
					{
						const Vector3r gradW = sim->gradW(xi - xj);
						sumGradW += gradW;
						sumGradW2 += gradW.squaredNorm();
					}
					xj[1] += diam;
				}
				xj[0] += diam;
				xj[1] = -supportRadius;
			}
		}
		else
		{
			Vector3r xj = { -supportRadius, -supportRadius, -supportRadius };
			while (xj[0] <= supportRadius)
			{
				while (xj[1] <= supportRadius)
				{
					while (xj[2] <= supportRadius)
					{
						// check if xj is in the support of xi
						if ((xi - xj).squaredNorm() < supportRadius*supportRadius)
						{
							const Vector3r gradW = sim->gradW(xi - xj);
							sumGradW += gradW;
							sumGradW2 += gradW.squaredNorm();
						}
						xj[2] += diam;
					}
					xj[1] += diam;
					xj[2] = -supportRadius;
				}
				xj[0] += diam;
				xj[1] = -supportRadius;
				xj[2] = -supportRadius;
			}
		}
 
		const Real beta = static_cast<Real>(2.0) * model->getVolume(0)*model->getVolume(0);
		m_pcisph_factor[fluidModelIndex] = static_cast<Real>(1.0) / (beta * (sumGradW.squaredNorm() + sumGradW2));
	}
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
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		const unsigned int numPart = fm->numActiveParticles();
		if (numPart != 0)
		{
			auto const& d = sim->getNeighborhoodSearch()->point_set(fm->getPointSetIndex());
			d.sort_field(&m_lastX[i][0]);
			d.sort_field(&m_lastV[i][0]);
			d.sort_field(&m_densityAdv[i][0]);
			d.sort_field(&m_pressure[i][0]);
			d.sort_field(&m_pressureAccel[i][0]);
		}
	}
}

void SimulationDataPCISPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	// initialize values for new particles
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_lastX[fluidModelIndex][j] = model->getPosition(j);
		m_lastV[fluidModelIndex][j] = model->getVelocity(j);
	}
}