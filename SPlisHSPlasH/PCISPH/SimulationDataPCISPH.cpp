#include "SimulationDataPCISPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include <iostream>
#include "SPlisHSPlasH/TimeManager.h"

using namespace SPH;

SimulationDataPCISPH::SimulationDataPCISPH()
{
	m_model = NULL;
}

SimulationDataPCISPH::~SimulationDataPCISPH(void)
{
	cleanup();
}

void SimulationDataPCISPH::init(FluidModel *model)
{
	m_model = model;

	m_lastX.resize(model->numParticles(), SPH::Vector3r::Zero());
	m_lastV.resize(model->numParticles(), SPH::Vector3r::Zero());
	m_pressureAccel.resize(model->numParticles(), SPH::Vector3r::Zero());
	m_densityAdv.resize(model->numParticles(), 0.0);
	m_pressure.resize(model->numParticles(), 0.0);
	m_pressureAccel.resize(model->numParticles(), SPH::Vector3r::Zero());

	std::cout << "Initialize PCISPH scaling factor\n";
	m_pcisph_factor = 0.0;
	model->getNeighborhoodSearch()->find_neighbors();

	// Find prototype particle
	// => particle with max. fluid neighbors
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real h2 = h*h;
	const Real density0 = m_model->getDensity0();
	unsigned int index = 0;
	unsigned int maxNeighbors = 0;

	#pragma omp parallel default(shared)
	{
		
		#pragma omp for schedule(static) nowait 
		for (int i = 0; i < (int)model->numParticles(); i++)
		{
			unsigned int counter = 0;
			for (unsigned int j = 0; j < model->numberOfNeighbors(i); j++)
			{
				const CompactNSearch::PointID &particleId = model->getNeighbor(i, j);
				if (particleId.point_set_id == 0)
					counter++;
			}
			if (counter > maxNeighbors)
			{
				maxNeighbors = counter;
				index = i;
			}
		}
	}

	Vector3r sumGradW = Vector3r::Zero();
	Real sumGradW2 = 0.0;
	const Vector3r &xi = model->getPosition(0, index);
	for (unsigned int j = 0; j < model->numberOfNeighbors(index); j++)
	{
		const CompactNSearch::PointID &particleId = model->getNeighbor(index, j);
		if (particleId.point_set_id == 0)
		{
			const unsigned int &neighborIndex = particleId.point_id;
			const Vector3r &xj = model->getPosition(particleId.point_set_id, neighborIndex);
			const Vector3r gradW = m_model->gradW(xi - xj);
			sumGradW += gradW;
			sumGradW2 += gradW.squaredNorm();
		}
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
	const unsigned int numPart = m_model->numParticles();
	if (numPart == 0)
		return;

	auto const& d = m_model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_lastX[0]);
	d.sort_field(&m_lastV[0]);
	d.sort_field(&m_densityAdv[0]);
	d.sort_field(&m_pressure[0]);
	d.sort_field(&m_pressureAccel[0]);
}
