#include "BoundaryModel.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "Utilities/Logger.h"
#include "CompactNSearch.h"
#include "Simulation.h"

using namespace SPH;


BoundaryModel::BoundaryModel() :
	m_x0(),
	m_x(),
	m_v(),
	m_V(),
	m_boundaryPsi(),
	m_forcePerThread(),
	m_torquePerThread()
{		
	m_sorted = false;
	m_pointSetIndex = 0;
}

BoundaryModel::~BoundaryModel(void)
{
	m_x0.clear();
	m_x.clear();
	m_v.clear();
	m_V.clear();
	m_boundaryPsi.clear();
	m_forcePerThread.clear();
	m_torquePerThread.clear();

	delete m_rigidBody;
}

void BoundaryModel::reset()
{
	// reset velocities and accelerations
	for (int j = 0; j < (int)numberOfParticles(); j++)
	{
		m_x[j] = m_x0[j];
		m_v[j].setZero();
	}

	#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
	#else
	const int maxThreads = 1;
	#endif

	for (int j = 0; j < maxThreads; j++)
	{
		m_forcePerThread[j].setZero();
		m_torquePerThread[j].setZero();
	}
}

void BoundaryModel::computeBoundaryPsi(const Real density0)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	CompactNSearch::NeighborhoodSearch *neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();
	
	const unsigned int numBoundaryParticles = numberOfParticles();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numBoundaryParticles; i++)
		{
			Real delta = sim->W_zero();
			for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++)
			{
				BoundaryModel *bm_neighbor = sim->getBoundaryModelFromPointSet(pid);
				for (unsigned int j = 0; j < neighborhoodSearch->point_set(m_pointSetIndex).n_neighbors(pid, i); j++)
				{
					const unsigned int neighborIndex = neighborhoodSearch->point_set(m_pointSetIndex).neighbor(pid, i, j);
					delta += sim->W(getPosition(i) - bm_neighbor->getPosition(neighborIndex));
				}
			}
			const Real volume = static_cast<Real>(1.0) / delta;
			m_V[i] = volume;
			m_boundaryPsi[i] = density0 * volume; 
		}
	}
}

void BoundaryModel::initModel(RigidBodyObject *rbo, const unsigned int numBoundaryParticles, Vector3r *boundaryParticles)
{
	m_x0.resize(numBoundaryParticles);
	m_x.resize(numBoundaryParticles);
	m_v.resize(numBoundaryParticles);
	m_V.resize(numBoundaryParticles);
	m_boundaryPsi.resize(numBoundaryParticles);
	
	#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
	#else
	const int maxThreads = 1;
	#endif

	m_forcePerThread.resize(maxThreads, Vector3r::Zero());
	m_torquePerThread.resize(maxThreads, Vector3r::Zero());

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numBoundaryParticles; i++)
		{
			m_x0[i] = boundaryParticles[i];
			m_x[i] = boundaryParticles[i];
			m_v[i].setZero();
			m_V[i] = 0.0;
			m_boundaryPsi[i] = 0.0;
		}
	}
	m_rigidBody = rbo;

	CompactNSearch::NeighborhoodSearch *neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();
	m_pointSetIndex = neighborhoodSearch->add_point_set(&m_x[0][0], m_x.size(), m_rigidBody->isDynamic(), false, true, this);
}

void BoundaryModel::performNeighborhoodSearchSort()
{
	const unsigned int numPart = numberOfParticles();

	// sort static boundaries only once
	if ((numPart == 0) || (!m_rigidBody->isDynamic() && !m_sorted))
		return;

	CompactNSearch::NeighborhoodSearch *neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();

	auto const& d = neighborhoodSearch->point_set(m_pointSetIndex);  
	d.sort_field(&m_x[0]);
	d.sort_field(&m_v[0]);
	d.sort_field(&m_V[0]);
	d.sort_field(&m_boundaryPsi[0]);
	m_sorted = true;
}

void BoundaryModel::getForceAndTorque(Vector3r &force, Vector3r &torque)
{
	#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
	#else
	const int maxThreads = 1;
	#endif

	force.setZero();
	for (int j = 0; j < maxThreads; j++)
	{
		force += m_forcePerThread[j];
		torque += m_torquePerThread[j];
	}
}

void BoundaryModel::clearForceAndTorque()
{
	#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
	#else
	const int maxThreads = 1;
	#endif

	for (int j = 0; j < maxThreads; j++)
	{
		m_forcePerThread[j].setZero();
		m_torquePerThread[j].setZero();
	}
}
