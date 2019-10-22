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
		virtual Vector3r const& getPosition() const = 0;
		virtual void setPosition(const Vector3r &x) = 0;
		virtual Vector3r getWorldSpacePosition() const = 0;
		virtual Vector3r const& getVelocity() const = 0;
		virtual void setVelocity(const Vector3r &v) = 0;
		virtual Matrix3r const& getRotation() const = 0;
		virtual void setRotation(const Matrix3r &r) = 0;
		virtual Matrix3r getWorldSpaceRotation() const = 0;
		virtual Vector3r const& getAngularVelocity() const = 0;
		virtual void setAngularVelocity(const Vector3r &v) = 0;
		virtual void addForce(const Vector3r &f) = 0;
		virtual void addTorque(const Vector3r &t) = 0;
	};
}

#endif 