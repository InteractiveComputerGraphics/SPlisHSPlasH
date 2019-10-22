#include "BoundaryModel.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "Utilities/Logger.h"
#include "NeighborhoodSearch.h"
#include "Simulation.h"

using namespace SPH;


BoundaryModel::BoundaryModel() :
	m_forcePerThread(),
	m_torquePerThread()
{		
}

BoundaryModel::~BoundaryModel(void)
{
	m_forcePerThread.clear();
	m_torquePerThread.clear();

	delete m_rigidBody;
}

void BoundaryModel::reset()
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

void BoundaryModel::getForceAndTorque(Vector3r &force, Vector3r &torque)
{
	#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
	#else
	const int maxThreads = 1;
	#endif

	force.setZero();
	torque.setZero();
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
