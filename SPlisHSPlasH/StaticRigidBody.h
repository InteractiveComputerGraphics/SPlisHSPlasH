#ifndef __StaticRigidBody_h__
#define __StaticRigidBody_h__

#include "Common.h"
#include "RigidBodyObject.h"
#include "TriangleMesh.h"

namespace SPH 
{	
	/** \brief This class stores the information of a static rigid body which 
	* is not part of a rigid body simulation. 
	*/
	class StaticRigidBody : public RigidBodyObject 
	{
	protected: 
		Vector3r m_x;
		Vector3r m_zero;
		Matrix3r m_R;
		TriangleMesh m_geometry;

	public:
		StaticRigidBody() { m_zero = Vector3r::Zero(); }

		virtual bool isDynamic() const { return false; }

		virtual Real const getMass() const { return 0.0; }
		virtual Vector3r const& getPosition() { return m_x; }
		virtual Vector3r const& getVelocity() const { return m_zero; }
		virtual Matrix3r const& getRotation() const { return m_R; }
		virtual Vector3r const& getAngularVelocity() const { return m_zero; }
		virtual void addForce(const Vector3r &f) {}
		virtual void addTorque(const Vector3r &t) {}

		void setPosition(const Vector3r &x) { m_x = x; }
		void setRotation(const Matrix3r &r) { m_R = r; }
		TriangleMesh& getGeometry() { return m_geometry; }
	};
}

#endif 