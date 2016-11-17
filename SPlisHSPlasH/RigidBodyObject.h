#ifndef __RigidBodyObject_h__
#define __RigidBodyObject_h__

#include "Common.h"

namespace SPH 
{	
	/** \brief Base class for rigid body objects. 
	*/
	class RigidBodyObject 
	{
	public:
		virtual ~RigidBodyObject() {};

		virtual bool isDynamic() const = 0;

		virtual Real const getMass() const = 0;
		virtual Vector3r const& getPosition() = 0;
		virtual Vector3r const& getVelocity() const = 0;
		virtual Matrix3r const& getRotation() const = 0;
		virtual Vector3r const& getAngularVelocity() const = 0;
		virtual void addForce(const Vector3r &f) = 0;
		virtual void addTorque(const Vector3r &t) = 0;
	};
}

#endif 