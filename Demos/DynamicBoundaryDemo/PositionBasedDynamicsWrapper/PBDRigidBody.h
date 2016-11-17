#ifndef __PBDRigidBody_h__
#define __PBDRigidBody_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/RigidBodyObject.h"
#include "Demos/Simulation/RigidBody.h"
#include "SPlisHSPlasH/TimeManager.h"

namespace SPH 
{	
	class PBDRigidBody : public RigidBodyObject 
	{
	protected: 
		Real m_h;
		PBD::RigidBody *m_rigidBody;	

	public:
		PBDRigidBody(PBD::RigidBody *rigidBody) : m_rigidBody(rigidBody), m_h(0.0) {}

		void updateTimeStepSize() { m_h = TimeManager::getCurrent()->getTimeStepSize(); }

		virtual bool isDynamic() const { return m_rigidBody->getMass() != 0.0; }

		virtual Real const getMass() const { return m_rigidBody->getMass(); }
		virtual Vector3r const& getPosition() { return m_rigidBody->getPosition(); }
		virtual Vector3r const& getVelocity() const { return m_rigidBody->getVelocity(); }
		virtual Matrix3r const& getRotation() const { return m_rigidBody->getRotationMatrix(); }
		virtual Vector3r const& getAngularVelocity() const { return m_rigidBody->getAngularVelocity(); }
		virtual void addForce(const Vector3r &f) { m_rigidBody->getVelocity() += (1.0/ m_rigidBody->getMass()) * f * m_h; }
		virtual void addTorque(const Vector3r &t) { m_rigidBody->getAngularVelocity() += m_rigidBody->getInertiaTensorInverseW() * t * m_h; }
	};
}

#endif 