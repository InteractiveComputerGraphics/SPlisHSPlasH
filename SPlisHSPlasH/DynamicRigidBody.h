#ifndef __DynamicRigidBody_h__
#define __DynamicRigidBody_h__

#include "Common.h"
#include "RigidBodyObject.h"
#include "TriangleMesh.h"
#include "TimeManager.h"
#include "Simulation.h"

namespace SPH {
	/** \brief This class stores the information of a dynamic rigid body which
	* is used for the strong coupling method introduced in 
	* Interlinked SPH Pressure Solvers for Strong Fluid-Rigid Coupling. Gissler et al. https://doi.org/10.1145/3284980
	*/
	class DynamicRigidBody : public RigidBodyObject {
		// Some fields are from PBD::RigidBody
	private:
		bool m_isDynamic;

		Real m_density;

		/** mass */
		Real m_mass;
		/** inverse mass */
		Real m_invMass;
		/** center of mass */
		Vector3r m_x;
		Vector3r m_lastX;
		Vector3r m_oldX;
		Vector3r m_x0;
		/** center of mass velocity */
		Vector3r m_v;
		Vector3r m_v0;
		/** acceleration (by external forces) */
		Vector3r m_a;
		/* external forces*/
		Vector3r m_force;

		/** Inertia tensor in the principal axis system: \n
		* After the main axis transformation the inertia tensor is a diagonal matrix.
		* So only three values are required to store the inertia tensor. These values
		* are constant over time.
		*/
		Vector3r m_inertiaTensor;
		/** 3x3 matrix, inertia tensor in world space */
		Matrix3r m_inertiaTensorW;
		/** Inverse inertia tensor in body space */
		Vector3r m_inertiaTensorInverse;
		/** 3x3 matrix, inverse of the inertia tensor in world space */
		Matrix3r m_inertiaTensorInverseW;
		/** Quaternion that describes the rotation of the body in world space */
		Quaternionr m_q;
		Quaternionr m_lastQ;
		Quaternionr m_oldQ;
		Quaternionr m_q0;
		/** Quaternion representing the rotation of the main axis transformation
		that is performed to get a diagonal inertia tensor */
		Quaternionr m_q_mat;
		/** Quaternion representing the initial rotation of the geometry */
		Quaternionr m_q_initial;
		/** difference of the initial translation and the translation of the main axis transformation */
		Vector3r m_x0_mat;
		/** rotationMatrix = 3x3 matrix.
		* Important for the transformation from world in body space and vice versa.
		* When using quaternions the rotation matrix is computed out of the quaternion.
		*/
		Matrix3r m_rot;
		/** Angular velocity, defines rotation axis and velocity (magnitude of the vector) */
		Vector3r m_omega;
		Vector3r m_omega0;
		/* Angular accelaration*/
		Vector3r m_a_omega;
		/** external torque */
		Vector3r m_torque;

		// used to recompute mass properties, may dont need after implement the volume integrateion
		Vector3r m_scale;

		// Real m_restitutionCoeff;
		Real m_frictionCoeff;

		// transformation required to transform a point to local space or vice vera
		Matrix3r m_transformation_R;
		Vector3r m_transformation_v1;
		Vector3r m_transformation_v2;
		Vector3r m_transformation_R_X_v1;
	protected:

		TriangleMesh m_geometry;

	public:  
		DynamicRigidBody() {
			m_isAnimated = true;
		}


		~DynamicRigidBody(void)
		{
		}

		void initBody(const Real density, const bool isDynamic, const Vector3r &position, const Quaternionr &rotation,
			 const Vector3r &scale)
		{
			m_density = density;
			m_scale = scale;
			determineMassProperties(density, isDynamic, scale);
			m_x = position;
			m_x0 = position;
			m_lastX = position;
			m_oldX = position;
			m_v.setZero();
			m_v0.setZero();
			m_a.setZero();
			m_force.setZero();

			m_q = rotation;
			m_q0 = rotation;
			m_lastQ = rotation;
			m_oldQ = rotation;
			m_rot = m_q.matrix();
			rotationUpdated();
			m_omega.setZero();
			m_omega0.setZero();
			m_torque.setZero();

			//m_restitutionCoeff = static_cast<Real>(0.6);
			m_frictionCoeff = static_cast<Real>(0.2);

			updateMeshTransformation();
		}

		void reset()
		{
			m_x = m_x0;
			getOldPosition() = getPosition0();
			getLastPosition() = getPosition0();

			m_q = m_q0;
			getOldRotation() = getRotation0();
			getLastRotation() = getRotation0();

			m_v = m_v0;
			m_omega = m_omega0;

			getAcceleration().setZero();
			getForce().setZero();
			getTorque().setZero();

			rotationUpdated();

			updateMeshTransformation();

			clearForceAndTorque();
		}

		void updateInverseTransformation()
		{
			// remove the rotation of the main axis transformation that is performed
			// to get a diagonal inertia tensor since the distance function is 
			// evaluated in local coordinates
			//
			// transformation world to local:
			// p_local = R_initial^T ( R_MAT R^T (p_world - x) - x_initial + x_MAT)
			// 
			// transformation local to world:
			// p_world = R R_MAT^T (R_initial p_local + x_initial - x_MAT) + x
			//
			m_transformation_R = (getRotationInitial().inverse() * getRotationMAT() * getRotation().inverse()).matrix();
			m_transformation_v1 = -getRotationInitial().inverse().matrix() * getPositionInitial_MAT();
			m_transformation_v2 = (getRotation()*getRotationMAT().inverse()).matrix() * getPositionInitial_MAT() + getPosition();
			m_transformation_R_X_v1 = -m_transformation_R * getPosition() + m_transformation_v1;
		}

		// Determine mass and inertia tensor
		void determineMassProperties(const Real density, bool isDynamic, const Vector3r scale) {
			m_isDynamic = isDynamic;
			// for now only consider cubiod which is scaled from a unit cube
			setMass(density * scale.x() * scale.y() * scale.z());
			Vector3r value = m_mass * Vector3r((scale.y() * scale.y() + scale.z() * scale.z()) / 12, (scale.x() * scale.x() + scale.z() * scale.z()) / 12, (scale.x() * scale.x() + scale.z() * scale.z()) / 12);
			m_inertiaTensor = value;
			m_inertiaTensorInverse = Vector3r(static_cast<Real>(1.0) / value[0], static_cast<Real>(1.0) / value[1], static_cast<Real>(1.0) / value[2]);
		}

		void rotationUpdated()
		{
			if (m_mass != 0.0)
			{
				m_rot = m_q.matrix();
				updateInertiaW();
				updateInverseTransformation();
			}
		}

		void updateInertiaW()
		{
			if (m_mass != 0.0)
			{
				m_inertiaTensorW = m_rot * m_inertiaTensor.asDiagonal() * m_rot.transpose();
				m_inertiaTensorInverseW = m_rot * m_inertiaTensorInverse.asDiagonal() * m_rot.transpose();
			}
		}

		/** Determine mass and inertia tensor of the given geometry.
			*/
		//void determineMassProperties(const Real density)
		//{
		//	// apply initial rotation
		//	VertexData &vd = m_geometry.getVertexDataLocal();
		//		
		//	Utilities::VolumeIntegration vi(m_geometry.getVertexDataLocal().size(), m_geometry.getMesh().numFaces(), &m_geometry.getVertexDataLocal().getPosition(0), m_geometry.getMesh().getFaces().data());
		//	vi.compute_inertia_tensor(density);
		//
		//	// Diagonalize Inertia Tensor
		//	Eigen::SelfAdjointEigenSolver<Matrix3r> es(vi.getInertia());
		//	Vector3r inertiaTensor = es.eigenvalues();
		//	Matrix3r R = es.eigenvectors();
		//
		//	setMass(1.0);
		//	setInertiaTensor(inertiaTensor);
		//
		//	if (R.determinant() < 0.0)
		//		R = -R;
		//
		//	for (unsigned int i = 0; i < vd.size(); i++)
		//		vd.getPosition(i) = m_rot * vd.getPosition(i) + m_x0;
		//
		//	Vector3r x_MAT = vi.getCenterOfMass();
		//	R = m_rot * R;
		//	x_MAT = m_rot * x_MAT + m_x0;
		//
		//	// rotate vertices back				
		//	for (unsigned int i = 0; i < vd.size(); i++)
		//		vd.getPosition(i) = R.transpose() * (vd.getPosition(i) - x_MAT);
		//
		//	// set rotation
		//	Quaternionr qR = Quaternionr(R);
		//	qR.normalize();
		//	m_q_mat = qR;
		//	m_q_initial = m_q0;
		//	m_x0_mat = m_x0 - x_MAT;
		//
		//	m_q0 = qR;
		//	m_q = m_q0;
		//	m_lastQ = m_q0;
		//	m_oldQ = m_q0;
		//	rotationUpdated();
		//
		//	// set translation
		//	m_x0 = x_MAT;
		//	m_x = m_x0;
		//	m_lastX = m_x0;
		//	m_oldX = m_x0;
		//	updateInverseTransformation();
		//}

		const Matrix3r &getTransformationR() { return m_transformation_R;  }
		const Vector3r &getTransformationV1() { return m_transformation_v1; }
		const Vector3r &getTransformationV2() { return m_transformation_v2; }
		const Vector3r &getTransformationRXV1() { return m_transformation_R_X_v1; }

		FORCE_INLINE const Vector3r& getScale() {
			return m_scale;
		}

		FORCE_INLINE void setMass(const Real &value)
		{
			m_mass = value;
			if (m_mass != 0.0)
				m_invMass = static_cast<Real>(1.0) / m_mass;
			else
				m_invMass = 0.0;
		}

		FORCE_INLINE const Real &getInvMass() const
		{
			return m_invMass;
		}

		FORCE_INLINE const Real& getDensity() const {
			return m_density;
		}

		FORCE_INLINE Real& getDensity() {
			return m_density;
		}

		FORCE_INLINE void setDensity(const Real& density) {
			m_density = density;
			determineMassProperties(density, isDynamic(), m_scale);
		}

		FORCE_INLINE Vector3r &getLastPosition()
		{
			return m_lastX;
		}

		FORCE_INLINE const Vector3r &getLastPosition() const
		{
			return m_lastX;
		}

		FORCE_INLINE void setLastPosition(const Vector3r &pos)
		{
			m_lastX = pos;
		}

		FORCE_INLINE Vector3r &getOldPosition()
		{
			return m_oldX;
		}

		FORCE_INLINE const Vector3r &getOldPosition() const
		{
			return m_oldX;
		}

		FORCE_INLINE void setOldPosition(const Vector3r &pos)
		{
			m_oldX = pos;
		}

		FORCE_INLINE Vector3r &getForce()
		{
			return m_force;
		}

		FORCE_INLINE const Vector3r &getForce() const
		{
			return m_force;
		}

		FORCE_INLINE void setForce(const Vector3r &force)
		{
			m_force = force;
		}

		FORCE_INLINE Vector3r &getPositionInitial_MAT()
		{
			return m_x0_mat;
		}

		FORCE_INLINE const Vector3r &getPositionInitial_MAT() const
		{
			return m_x0_mat;
		}

		FORCE_INLINE void setPositionInitial_MAT(const Vector3r &pos)
		{
			m_x0_mat = pos;
		}

		FORCE_INLINE Vector3r &getAcceleration()
		{
			return m_a;
		}

		FORCE_INLINE const Vector3r &getAcceleration() const 
		{
			return m_a;
		}

		FORCE_INLINE void setAcceleration(const Vector3r &accel)
		{
			m_a = accel;
		}

		FORCE_INLINE const Vector3r &getInertiaTensor() const
		{
			return m_inertiaTensor;
		}

		FORCE_INLINE void setInertiaTensor(const Vector3r &value)
		{
			m_inertiaTensor = value;
			m_inertiaTensorInverse = Vector3r(static_cast<Real>(1.0) / value[0], static_cast<Real>(1.0) / value[1], static_cast<Real>(1.0) / value[2]);
		}

		FORCE_INLINE Matrix3r& getInertiaTensorW()
		{
			return m_inertiaTensorW;
		}

		FORCE_INLINE const Matrix3r& getInertiaTensorW() const
		{
			return m_inertiaTensorW;
		}

		FORCE_INLINE const Vector3r &getInertiaTensorInverse() const
		{
			return m_inertiaTensorInverse;
		}

		FORCE_INLINE Matrix3r &getInertiaTensorInverseW()
		{
			return m_inertiaTensorInverseW;
		}

		FORCE_INLINE const Matrix3r &getInertiaTensorInverseW() const
		{
			return m_inertiaTensorInverseW;
		}

		FORCE_INLINE void setInertiaTensorInverseW(const Matrix3r &value)
		{
			m_inertiaTensorInverseW = value;
		}

		FORCE_INLINE Quaternionr &getLastRotation()
		{
			return m_lastQ;
		}

		FORCE_INLINE const Quaternionr &getLastRotation() const
		{
			return m_lastQ;
		}

		FORCE_INLINE void setLastRotation(const Quaternionr &value)
		{
			m_lastQ = value;
		}

		FORCE_INLINE Quaternionr &getOldRotation()
		{
			return m_oldQ;
		}

		FORCE_INLINE const Quaternionr &getOldRotation() const
		{
			return m_oldQ;
		}

		FORCE_INLINE void setOldRotation(const Quaternionr &value)
		{
			m_oldQ = value;
		}

		FORCE_INLINE Quaternionr &getRotationMAT()
		{
			return m_q_mat;
		}

		FORCE_INLINE const Quaternionr &getRotationMAT() const
		{
			return m_q_mat;
		}

		FORCE_INLINE void setRotationMAT(const Quaternionr &value)
		{
			m_q_mat = value;
		}

		FORCE_INLINE Quaternionr &getRotationInitial()
		{
			return m_q_initial;
		}

		FORCE_INLINE const Quaternionr &getRotationInitial() const
		{
			return m_q_initial;
		}

		FORCE_INLINE void setRotationInitial(const Quaternionr &value)
		{
			m_q_initial = value;
		}

		FORCE_INLINE Matrix3r &getRotationMatrix()
		{
			return m_rot;
		}

		FORCE_INLINE const Matrix3r &getRotationMatrix() const
		{
			return m_rot;
		}

		FORCE_INLINE void setRotationMatrix(const Matrix3r &value)
		{
			m_rot = value;
		}

		FORCE_INLINE Vector3r &getAngularAccelaration()
		{
			return m_a_omega;
		}

		FORCE_INLINE const Vector3r &getAngularAccelaration() const
		{
			return m_a_omega;
		}

		FORCE_INLINE void setAngularAccelaration(const Vector3r &value)
		{
			m_a_omega = value;
		}

		FORCE_INLINE Vector3r &getAngularVelocity0()
		{
			return m_omega0;
		}

		FORCE_INLINE const Vector3r &getAngularVelocity0() const
		{
			return m_omega0;
		}

		FORCE_INLINE void setAngularVelocity0(const Vector3r &value)
		{
			m_omega0 = value;
		}

		FORCE_INLINE Vector3r &getTorque()
		{
			return m_torque;
		}

		FORCE_INLINE const Vector3r &getTorque() const
		{
			return m_torque;
		}

		FORCE_INLINE void setTorque(const Vector3r &value)
		{
			m_torque = value;
		}

		FORCE_INLINE Real getFrictionCoeff() const 
		{ 
			return m_frictionCoeff; 
		}

		FORCE_INLINE void setFrictionCoeff(Real val) 
		{ 
			m_frictionCoeff = val; 
		}


		virtual bool isDynamic() const {
			return m_isDynamic;
		}

		virtual Real const getMass() const {
			return m_mass;
		}
		virtual Vector3r const& getPosition() const {
			return m_x;
		}
		virtual void setPosition(const Vector3r& x) {
			m_x = x;
		}
		Vector3r const& getPosition0() const {
			return m_x0;
		}
		void setPosition0(const Vector3r& x) {
			m_x0 = x;
		}
		virtual Vector3r getWorldSpacePosition() const {
			return m_x;
		}
		virtual Vector3r const& getVelocity() const {
			return m_v;
		}

		FORCE_INLINE Vector3r& getVelocity() {
			return m_v;
		}

		virtual void setVelocity(const Vector3r& v) {
			if (m_isAnimated) m_v = v;
		}
		virtual Quaternionr const& getRotation() const {
			return m_q;
		}
		virtual void setRotation(const Quaternionr& q) {
			m_q = q;
		}
		Quaternionr const& getRotation0() const {
			return m_q0;
		}
		void setRotation0(const Quaternionr& q) {
			m_q0 = q;
		}
		virtual Matrix3r getWorldSpaceRotation() const {
			return m_q.toRotationMatrix();
		}
		virtual Vector3r const& getAngularVelocity() const {
			return m_omega;
		}
		virtual void setAngularVelocity(const Vector3r& v) {
			if (m_isAnimated) m_omega = v;
		}
		virtual void addForce(const Vector3r& f) {
			m_force += f;
		}

		virtual void addTorque(const Vector3r& t) {
			m_torque += t;
		}

		void clearForceAndTorque() {
			m_force.setZero();
			m_torque.setZero();
		}

		void timeStep() {
			Simulation* sim = Simulation::getCurrent();
			const Vector3r gravAccel(sim->getVecValue<Real>(Simulation::GRAVITATION));
			const Real dt = TimeManager::getCurrent()->getTimeStepSize();
			// Save old values
			m_lastX = m_oldX;
			m_oldX = m_x;
			m_lastQ = m_oldQ;
			m_oldQ = m_q;
			// Semi-implicit Euler
			m_a = m_invMass * m_force + gravAccel;
			m_v += m_a * dt;
			m_x += m_v * dt;
			m_a_omega = m_inertiaTensorInverseW * (m_torque - (m_omega.cross(m_inertiaTensorW * m_omega)));
			m_omega += m_a_omega * dt;
			Quaternionr omegaTilde(0.0, m_omega[0], m_omega[1], m_omega[2]);
			m_q.coeffs() += 0.5 * (omegaTilde * m_q).coeffs() * dt;
			m_q.normalize();
			rotationUpdated();
			updateMeshTransformation();
			clearForceAndTorque();
		}

		virtual const std::vector<Vector3r>& getVertices() const {
			return m_geometry.getVertices();
		};
		virtual const std::vector<Vector3r>& getVertexNormals() const {
			return m_geometry.getVertexNormals();
		};
		virtual const std::vector<unsigned int>& getFaces() const {
			return m_geometry.getFaces();
		};

		void setWorldSpacePosition(const Vector3r& x) {
			m_x = x;
		}
		void setWorldSpaceRotation(const Matrix3r& r) {
			m_q = Quaternionr(r);
		}
		TriangleMesh& getGeometry() {
			return m_geometry;
		}

		virtual void updateMeshTransformation() {
			m_geometry.updateMeshTransformation(m_x, m_q.toRotationMatrix());
			m_geometry.updateNormals();
			m_geometry.updateVertexNormals();
		}
	};
}

#endif 