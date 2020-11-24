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
	for (int j = 0; j < m_forcePerThread.size(); j++)
	{
		m_forcePerThread[j].setZero();
		m_torquePerThread[j].setZero();
	}
}

void BoundaryModel::getForceAndTorque(Vector3r &force, Vector3r &torque)
{
	force.setZero();
	torque.setZero();
	for (int j = 0; j < m_forcePerThread.size(); j++)
	{
		force += m_forcePerThread[j];
		torque += m_torquePerThread[j];
	}
}

void BoundaryModel::clearForceAndTorque()
{
	for (int j = 0; j < m_forcePerThread.size(); j++)
	{
		m_forcePerThread[j].setZero();
		m_torquePerThread[j].setZero();
	}
}
