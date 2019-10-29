#ifndef __PBDWrapper_h__
#define __PBDWrapper_h__

#include "SPlisHSPlasH/Common.h"
#include <string>
#include "Visualization/Shader.h"
#include "extern/AntTweakBar/include/AntTweakBar.h"
#include "Simulation/SimulationModel.h"
#include "Visualization/MiniGL.h"
#include "Simulation/CubicSDFCollisionDetection.h"


namespace PBD
{
	class TimeStepController;
}

class PBDWrapper
{
protected:
	SPH::Shader *m_shader;
	PBD::SimulationModel m_model;
	PBD::CubicSDFCollisionDetection m_cd;
	PBD::TimeStepController *m_timeStep;

	short m_clothSimulationMethod = 2;
	short m_solidSimulationMethod = 2;
	short m_bendingMethod = 2;
	bool m_drawAABB;
	bool m_drawSDF;
	bool m_drawStaticBodies;
	int m_drawBVHDepth;
	std::string m_sceneName;
	std::string m_sceneFileName;
	bool m_enableMayaExport = false;
	Real m_dampingCoeff = 0.0;

	SPH::Shader *createShader(const std::string &vertexShader, const std::string &geometryShader, const std::string &fragmentShader);

public:
	struct RBData 
	{
		Vector3r x;
		Matrix3r R;
		Vector3r scale;
		std::string objFile;
		int collisionType;
		Real restitution;
		Real friction;
	};

	PBDWrapper();
	~PBDWrapper();

	void reset();
	void initGUI();
	void initModel(const Real timeStepSize);

	/** Read rigid body scene and create the rigid body model
	*/
	void readScene(const std::string &sceneFileName, const std::vector< RBData> &additionalRigidBodies);
	void initTriangleModelConstraints();
	void initTetModelConstraints();

	void timeStep();
	void updateVisModels();

	void initShader();

	void loadObj(const std::string &filename, PBD::VertexData &vd, Utilities::IndexedFaceMesh &mesh, const Vector3r &scale);

	PBD::SimulationModel &getSimulationModel() { return m_model; }
	PBD::DistanceFieldCollisionDetection &getCollisionDetection() { return m_cd; }
	PBD::TimeStepController &getTimeStepController();

	void shaderBegin(const float *col);
	void shaderEnd();
	void renderTriangleModels();
	void renderTetModels();
	void renderConstraints();
	void renderAABB(PBD::AABB &aabb);
	void renderBVH();
	void renderSDF();
	void renderSDF(PBD::CollisionDetection::CollisionObject* co);
	void renderScene();
	void renderBallJoint(PBD::BallJoint &bj);
	void renderRigidBodyParticleBallJoint(PBD::RigidBodyParticleBallJoint &bj);
	void renderBallOnLineJoint(PBD::BallOnLineJoint &bj);
	void renderHingeJoint(PBD::HingeJoint &hj);
	void renderUniversalJoint(PBD::UniversalJoint &uj);
	void renderSliderJoint(PBD::SliderJoint &joint);
	void renderTargetPositionMotorSliderJoint(PBD::TargetPositionMotorSliderJoint &joint);
	void renderTargetVelocityMotorSliderJoint(PBD::TargetVelocityMotorSliderJoint &joint);
	void renderTargetAngleMotorHingeJoint(PBD::TargetAngleMotorHingeJoint &hj);
	void renderTargetVelocityMotorHingeJoint(PBD::TargetVelocityMotorHingeJoint &hj);
	void renderRigidBodyContact(PBD::RigidBodyContactConstraint &cc);
	void renderParticleRigidBodyContact(PBD::ParticleRigidBodyContactConstraint &cc);
	void renderSpring(PBD::RigidBodySpring &s);
	void renderDistanceJoint(PBD::DistanceJoint &j);
	void renderDamperJoint(PBD::DamperJoint &joint);

	static void TW_CALL setVelocityUpdateMethod(const void *value, void *clientData);
	static void TW_CALL getVelocityUpdateMethod(void *value, void *clientData);
	static void TW_CALL setMaxIterations(const void *value, void *clientData);
	static void TW_CALL getMaxIterations(void *value, void *clientData);
	static void TW_CALL setMaxIterationsV(const void *value, void *clientData);
	static void TW_CALL getMaxIterationsV(void *value, void *clientData);
	static void TW_CALL setContactTolerance(const void *value, void *clientData);
	static void TW_CALL getContactTolerance(void *value, void *clientData);
	static void TW_CALL setContactStiffnessRigidBody(const void *value, void *clientData);
	static void TW_CALL getContactStiffnessRigidBody(void *value, void *clientData);
	static void TW_CALL setContactStiffnessParticleRigidBody(const void *value, void *clientData);
	static void TW_CALL getContactStiffnessParticleRigidBody(void *value, void *clientData);
	static void TW_CALL setStiffness(const void *value, void *clientData);
	static void TW_CALL getStiffness(void *value, void *clientData);
	static void TW_CALL setXXStiffness(const void *value, void *clientData);
	static void TW_CALL getXXStiffness(void *value, void *clientData);
	static void TW_CALL setYYStiffness(const void *value, void *clientData);
	static void TW_CALL getYYStiffness(void *value, void *clientData);
	static void TW_CALL setXYStiffness(const void *value, void *clientData);
	static void TW_CALL getXYStiffness(void *value, void *clientData);
	static void TW_CALL setXYPoissonRatio(const void *value, void *clientData);
	static void TW_CALL getXYPoissonRatio(void *value, void *clientData);
	static void TW_CALL setYXPoissonRatio(const void *value, void *clientData);
	static void TW_CALL getYXPoissonRatio(void *value, void *clientData);
	static void TW_CALL setNormalizeStretch(const void *value, void *clientData);
	static void TW_CALL getNormalizeStretch(void *value, void *clientData);
	static void TW_CALL setNormalizeShear(const void *value, void *clientData);
	static void TW_CALL getNormalizeShear(void *value, void *clientData);
	static void TW_CALL setBendingStiffness(const void *value, void *clientData);
	static void TW_CALL getBendingStiffness(void *value, void *clientData);
	static void TW_CALL setBendingMethod(const void *value, void *clientData);
	static void TW_CALL getBendingMethod(void *value, void *clientData);
	static void TW_CALL setClothSimulationMethod(const void *value, void *clientData);
	static void TW_CALL getClothSimulationMethod(void *value, void *clientData);
	static void TW_CALL setSolidSimulationMethod(const void *value, void *clientData);
	static void TW_CALL getSolidSimulationMethod(void *value, void *clientData);

	
public:
	template<class PositionData>
	static void drawMesh(const PositionData &pd, const Utilities::IndexedFaceMesh &mesh, const unsigned int offset, const float * const color);
};


template<class PositionData>
void PBDWrapper::drawMesh(const PositionData &pd, const Utilities::IndexedFaceMesh &mesh, const unsigned int offset, const float * const color)
{
	// draw mesh 
	const unsigned int *faces = mesh.getFaces().data();
	const unsigned int nFaces = mesh.numFaces();
	const Vector3r *vertexNormals = mesh.getVertexNormals().data();

	if (SPH::MiniGL::checkOpenGLVersion(3, 3))
	{
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &pd.getPosition(offset)[0]);
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 3, GL_REAL, GL_FALSE, 0, &vertexNormals[0][0]);
	}
	else
	{
		float speccolor[4] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, speccolor);
		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0f);
		glColor3fv(color);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glVertexPointer(3, GL_REAL, 0, &pd.getPosition(0)[0]);
		glNormalPointer(GL_REAL, 0, &vertexNormals[0][0]);
	}

	glDrawElements(GL_TRIANGLES, (GLsizei)3 * mesh.numFaces(), GL_UNSIGNED_INT, mesh.getFaces().data());

	if (SPH::MiniGL::checkOpenGLVersion(3, 3))
	{
		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(2);
	}
	else
	{
		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
	}
}



#endif
