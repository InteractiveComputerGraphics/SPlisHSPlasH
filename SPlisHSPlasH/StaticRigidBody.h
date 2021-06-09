#ifndef __StaticRigidBody_h__
#define __StaticRigidBody_h__

#include "Common.h"
#include "RigidBodyObject.h"
#include "TriangleMesh.h"
#include "SPlisHSPlasH/TimeManager.h"

namespace SPH 
{	
	/** \brief This class stores the information of a static rigid body which 
	* is not part of a rigid body simulation. 
	*/
	class StaticRigidBody : public RigidBodyObject 
	{
	protected: 
		Vector3r m_x0;
		Vector3r m_x;
		Quaternionr m_q;
		Quaternionr m_q0;
		Vector3r m_velocity;
		Vector3r m_angularVelocity;
		TriangleMesh m_geometry;

	public:
		StaticRigidBody() 
		{ 
			m_isAnimated = false; 
			m_velocity = Vector3r::Zero(); 
			m_angularVelocity = Vector3r::Zero();
		}

		virtual bool isDynamic() const { return false; }

		virtual Real const getMass() const { return 0.0; }
		virtual Vector3r const& getPosition() const  { return m_x; }
		virtual void setPosition(const Vector3r &x) { m_x = x; }
		Vector3r const& getPosition0() const { return m_x0; }
		void setPosition0(const Vector3r& x) { m_x0 = x; }
		virtual Vector3r getWorldSpacePosition() const { return m_x; }
		virtual Vector3r const& getVelocity() const { return m_velocity; }
		virtual void setVelocity(const Vector3r& v) { if (m_isAnimated) m_velocity = v; }
		virtual Quaternionr const& getRotation() const { return m_q; }
		virtual void setRotation(const Quaternionr& q) { m_q = q; }
		Quaternionr const& getRotation0() const { return m_q0; }
		void setRotation0(const Quaternionr& q) { m_q0 = q; }
		virtual Matrix3r getWorldSpaceRotation() const { return m_q.toRotationMatrix(); }
		virtual Vector3r const& getAngularVelocity() const { return m_angularVelocity; }
		virtual void setAngularVelocity(const Vector3r &v) { if (m_isAnimated) m_angularVelocity = v; }
		virtual void addForce(const Vector3r &f) {}
		virtual void addTorque(const Vector3r &t) {}
		void animate() 
		{
			const Real dt = TimeManager::getCurrent()->getTimeStepSize();
			m_x += m_velocity * dt; 
			Quaternionr angVelQ(0.0, m_angularVelocity[0], m_angularVelocity[1], m_angularVelocity[2]);
			m_q.coeffs() += dt * 0.5 * (angVelQ * m_q).coeffs();
			m_q.normalize();
			updateMeshTransformation();
		}
		
		virtual const std::vector<Vector3r> &getVertices() const { return m_geometry.getVertices(); };
		virtual const std::vector<Vector3r> &getVertexNormals() const { return m_geometry.getVertexNormals(); };
		virtual const std::vector<unsigned int> &getFaces() const { return m_geometry.getFaces(); };

		void setWorldSpacePosition(const Vector3r &x) { m_x = x; }
		void setWorldSpaceRotation(const Matrix3r &r) { m_q = Quaternionr(r); }
		TriangleMesh& getGeometry() { return m_geometry; }

		virtual void updateMeshTransformation()
		{
			m_geometry.updateMeshTransformation(m_x, m_q.toRotationMatrix());
			m_geometry.updateNormals();
			m_geometry.updateVertexNormals();
		}

		void reset()
		{
			m_x = m_x0;
			m_q = m_q0;
			updateMeshTransformation();
		}
	};
}

#endif 