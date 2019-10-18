#include "BoundaryModel_Akinci2012.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "Utilities/Logger.h"
#include "NeighborhoodSearch.h"
#include "Simulation.h"

using namespace SPH;


BoundaryModel_Akinci2012::BoundaryModel_Akinci2012() :
	m_x0(),
	m_x(),
	m_v(),
	m_V()
{		
	m_sorted = false;
	m_pointSetIndex = 0;
}

BoundaryModel_Akinci2012::~BoundaryModel_Akinci2012(void)
{
	m_x0.clear();
	m_x.clear();
	m_v.clear();
	m_V.clear();
}

void BoundaryModel_Akinci2012::reset()
{
	BoundaryModel::reset();

	// reset velocities and accelerations
	for (int j = 0; j < (int)numberOfParticles(); j++)
	{
		m_x[j] = m_x0[j];
		m_v[j].setZero();
	}
}

void BoundaryModel_Akinci2012::computeBoundaryVolume()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	NeighborhoodSearch *neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();

	const unsigned int numBoundaryParticles = numberOfParticles();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numBoundaryParticles; i++)
		{
			Real delta = sim->W_zero();
			for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++)
			{
				BoundaryModel_Akinci2012 *bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
				for (unsigned int j = 0; j < neighborhoodSearch->point_set(m_pointSetIndex).n_neighbors(pid, i); j++)
				{
					const unsigned int neighborIndex = neighborhoodSearch->point_set(m_pointSetIndex).neighbor(pid, i, j);
					delta += sim->W(getPosition(i) - bm_neighbor->getPosition(neighborIndex));
				}
			}
			const Real volume = static_cast<Real>(1.0) / delta;
			m_V[i] = volume;
		}
	}
}

void BoundaryModel_Akinci2012::initModel(RigidBodyObject *rbo, const unsigned int numBoundaryParticles, Vector3r *boundaryParticles)
{
	m_x0.resize(numBoundaryParticles);
	m_x.resize(numBoundaryParticles);
	m_v.resize(numBoundaryParticles);
	m_V.resize(numBoundaryParticles);

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
		}
	}
	m_rigidBody = rbo;

	NeighborhoodSearch *neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();
	m_pointSetIndex = neighborhoodSearch->add_point_set(&m_x[0][0], m_x.size(), m_rigidBody->isDynamic(), false, true, this);
}

void BoundaryModel_Akinci2012::performNeighborhoodSearchSort()
{
	const unsigned int numPart = numberOfParticles();

	// sort static boundaries only once
	if ((numPart == 0) || (!m_rigidBody->isDynamic() && m_sorted))
		return;

	NeighborhoodSearch *neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();

	auto const& d = neighborhoodSearch->point_set(m_pointSetIndex);  
	d.sort_field(&m_x0[0]);
	d.sort_field(&m_x[0]);
	d.sort_field(&m_v[0]);
	d.sort_field(&m_V[0]);
	m_sorted = true;
}

void SPH::BoundaryModel_Akinci2012::saveState(BinaryFileWriter &binWriter)
{
	binWriter.write(m_sorted);
	binWriter.write(m_pointSetIndex);
}

void SPH::BoundaryModel_Akinci2012::loadState(BinaryFileReader &binReader)
{
	binReader.read(m_sorted);
	binReader.read(m_pointSetIndex);
}

void SPH::BoundaryModel_Akinci2012::resize(const unsigned int numBoundaryParticles)
{
	m_x0.resize(numBoundaryParticles);
	m_x.resize(numBoundaryParticles);
	m_v.resize(numBoundaryParticles);
	m_V.resize(numBoundaryParticles);
}
