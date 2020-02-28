#include "BoundarySimulator.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/BoundaryModel.h"

using namespace SPH;

void BoundarySimulator::updateBoundaryForces()
{
	Real h = TimeManager::getCurrent()->getTimeStepSize();
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nObjects = sim->numberOfBoundaryModels();
	for (unsigned int i = 0; i < nObjects; i++)
	{
		BoundaryModel *bm = sim->getBoundaryModel(i);
		RigidBodyObject *rbo = bm->getRigidBodyObject();
		if (rbo->isDynamic())
		{
			Vector3r force, torque;
			bm->getForceAndTorque(force, torque);
			rbo->addForce(force);
			rbo->addTorque(torque);
			bm->clearForceAndTorque();
		}
	}
}
