#include "BoundaryModel_Koschier2017.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "Utilities/Logger.h"
#include "NeighborhoodSearch.h"
#include "Simulation.h"

using namespace SPH;


BoundaryModel_Koschier2017::BoundaryModel_Koschier2017() :
	m_boundaryDensity(),
	m_boundaryDensityGradient(),
	m_boundaryXj()
{		
	m_map = nullptr;
	m_maxDist = 0.0;
	m_maxVel = 0.0;
}

BoundaryModel_Koschier2017::~BoundaryModel_Koschier2017(void)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int i = 0; i < nModels; i++)
	{
		m_boundaryDensity[i].clear();
		m_boundaryDensityGradient[i].clear();
		m_boundaryXj[i].clear();
	}
	m_boundaryDensity.clear();
	m_boundaryDensityGradient.clear();
	m_boundaryXj.clear();
	delete m_map;
}


void BoundaryModel_Koschier2017::initModel(RigidBodyObject *rbo)
{
	Simulation *sim = Simulation::getCurrent();

	const unsigned int nModels = sim->numberOfFluidModels();
	m_boundaryDensity.resize(nModels);
	m_boundaryDensityGradient.resize(nModels);
	m_boundaryXj.resize(nModels);
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		m_boundaryDensity[i].resize(fm->numParticles(), 0.0);
		m_boundaryDensityGradient[i].resize(fm->numParticles(), Vector3r::Zero());
		m_boundaryXj[i].resize(fm->numParticles(), Vector3r::Zero());
	}

	#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
	#else
	const int maxThreads = 1;
	#endif
	m_forcePerThread.resize(maxThreads, Vector3r::Zero());
	m_torquePerThread.resize(maxThreads, Vector3r::Zero());

	m_rigidBody = rbo;
}

void BoundaryModel_Koschier2017::reset()
{
	BoundaryModel::reset();

	m_maxDist = 0.0;
	m_maxVel = 0.0;
}