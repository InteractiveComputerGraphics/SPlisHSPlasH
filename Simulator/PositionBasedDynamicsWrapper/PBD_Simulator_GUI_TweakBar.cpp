#include "PBD_Simulator_GUI_TweakBar.h"
#include "Utilities/FileSystem.h"
#include "Simulation/Simulation.h"
#include "Simulation/TimeStepController.h"
#include "GUI/TweakBar/TweakBarParameters.h"


using namespace SPH;
using namespace Utilities;

PBD_Simulator_GUI_TweakBar::PBD_Simulator_GUI_TweakBar(SimulatorBase *base, PBDWrapper *pbdWrapper) :
	Simulator_GUI_TweakBar(base)
{
	m_pbdWrapper = pbdWrapper;
	m_drawAABB = false;
	m_drawBVHDepth = -1;
	m_drawSDF = false;
	m_shader = nullptr;
	m_jointColor[0] = 0.0f;
	m_jointColor[1] = 0.6f;
	m_jointColor[2] = 0.2f;
	m_jointColor[3] = 1.0f;
}

PBD_Simulator_GUI_TweakBar::~PBD_Simulator_GUI_TweakBar(void)
{	
	delete m_shader;
}

void PBD_Simulator_GUI_TweakBar::init(int argc, char **argv, const char *name)
{
	Simulator_GUI_TweakBar::init(argc, argv, name);
	initShader();
}

void PBD_Simulator_GUI_TweakBar::render()
{
	Simulator_GUI_TweakBar::render();

	renderTriangleModels();
	renderTetModels();
	renderConstraints();
	renderBVH();
	renderSDF();
}


SPH::Shader *PBD_Simulator_GUI_TweakBar::createShader(const std::string &vertexShader, const std::string &geometryShader, const std::string &fragmentShader)
{
	if (SPH::MiniGL::checkOpenGLVersion(3, 3))
	{
		SPH::Shader *shader = new SPH::Shader();

		if (vertexShader != "")
			shader->compileShaderFile(GL_VERTEX_SHADER, vertexShader);
		if (geometryShader != "")
			shader->compileShaderFile(GL_GEOMETRY_SHADER, geometryShader);
		if (fragmentShader != "")
			shader->compileShaderFile(GL_FRAGMENT_SHADER, fragmentShader);
		shader->createAndLinkProgram();
		return shader;
	}
	return NULL;
}

void PBD_Simulator_GUI_TweakBar::initShader()
{
	std::string exePath = Utilities::FileSystem::getProgramPath();
	std::string vertFile = Utilities::FileSystem::normalizePath(exePath + "/resources/pbd_shaders/vs_smooth.glsl");
	std::string fragFile = Utilities::FileSystem::normalizePath(exePath + "/resources/pbd_shaders/fs_smooth.glsl");
	m_shader = createShader(vertFile, "", fragFile);

	if (m_shader == NULL)
		return;

	m_shader->begin();
	m_shader->addUniform("modelview_matrix");
	m_shader->addUniform("projection_matrix");
	m_shader->addUniform("surface_color");
	m_shader->addUniform("shininess");
	m_shader->addUniform("specular_factor");
	m_shader->end();
}


void PBD_Simulator_GUI_TweakBar::shaderBegin(const float *col)
{
	if (m_shader)
	{
		m_shader->begin();
		glUniform1f(m_shader->getUniform("shininess"), 5.0f);
		glUniform1f(m_shader->getUniform("specular_factor"), 0.2f);

		GLfloat matrix[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
		glUniformMatrix4fv(m_shader->getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
		GLfloat pmatrix[16];
		glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
		glUniformMatrix4fv(m_shader->getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);
		glUniform3fv(m_shader->getUniform("surface_color"), 1, col);
	}
}

void PBD_Simulator_GUI_TweakBar::shaderEnd()
{
	if (m_shader)
		m_shader->end();
}

void PBD_Simulator_GUI_TweakBar::renderAABB(PBD::AABB &aabb)
{
	Vector3r p1, p2;
	glBegin(GL_LINES);
	for (unsigned char i = 0; i < 12; i++)
	{
		PBD::AABB::getEdge(aabb, i, p1, p2);
		glVertex3d(p1[0], p1[1], p1[2]);
		glVertex3d(p2[0], p2[1], p2[2]);
	}
	glEnd();
}

void PBD_Simulator_GUI_TweakBar::renderBVH()
{
	if (m_drawAABB || (m_drawBVHDepth >= 0))
	{
		float staticColor[4] = { 0.5f, 0.5f, 0.5f, 0.3f };

		PBD::DistanceFieldCollisionDetection &cd = m_pbdWrapper->getCollisionDetection();
		PBD::SimulationModel::RigidBodyVector &rb = m_pbdWrapper->getSimulationModel().getRigidBodies();
		std::vector<PBD::CollisionDetection::CollisionObject*> &collisionObjects = cd.getCollisionObjects();
		for (unsigned int k = 0; k < collisionObjects.size(); k++)
		{
			if (m_drawAABB)
				renderAABB(collisionObjects[k]->m_aabb);

			if (m_drawBVHDepth >= 0)
			{
				if (cd.isDistanceFieldCollisionObject(collisionObjects[k]))
				{
					const PBD::PointCloudBSH &bvh = ((PBD::DistanceFieldCollisionDetection::DistanceFieldCollisionObject*) collisionObjects[k])->m_bvh;

					std::function<bool(unsigned int, unsigned int)> predicate = [&](unsigned int node_index, unsigned int depth) { return (int)depth <= m_drawBVHDepth; };
					std::function<void(unsigned int, unsigned int)> cb = [&](unsigned int node_index, unsigned int depth)
					{
						if (depth == m_drawBVHDepth)
						{
							const PBD::BoundingSphere &bs = bvh.hull(node_index);
							if (collisionObjects[k]->m_bodyType == PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType)
							{
								PBD::RigidBody *body = rb[collisionObjects[k]->m_bodyIndex];
								const Vector3r &sphere_x = bs.x();
								const Vector3r sphere_x_w = body->getRotation() * sphere_x + body->getPosition();
								SPH::MiniGL::drawSphere(sphere_x_w, std::max((float)bs.r(), 0.05f), staticColor);
							}
							else
								SPH::MiniGL::drawSphere(bs.x(), std::max((float)bs.r(), 0.05f), staticColor);
						}
					};

					bvh.traverse_depth_first(predicate, cb);
				}
			}
		}
	}
}


void PBD_Simulator_GUI_TweakBar::renderSDF()
{
	if (m_drawSDF)
	{
		PBD::DistanceFieldCollisionDetection &cd = m_pbdWrapper->getCollisionDetection();
		std::vector<PBD::CollisionDetection::CollisionObject*> &collisionObjects = cd.getCollisionObjects();
		for (unsigned int k = 0; k < collisionObjects.size(); k++)
		{
			renderSDF(collisionObjects[k]);
		}
	}
}

void PBD_Simulator_GUI_TweakBar::renderSDF(PBD::CollisionDetection::CollisionObject* co)
{
	PBD::DistanceFieldCollisionDetection &cd = m_pbdWrapper->getCollisionDetection();
	if ((!cd.isDistanceFieldCollisionObject(co)) || (co->m_bodyType != PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType))
		return;

	PBD::SimulationModel &model = m_pbdWrapper->getSimulationModel();
	const PBD::SimulationModel::RigidBodyVector &rigidBodies = model.getRigidBodies();
	PBD::RigidBody *rb = rigidBodies[co->m_bodyIndex];

	const Vector3r &com = rb->getPosition();
	const Matrix3r &R = rb->getTransformationR();
	const Vector3r &v1 = rb->getTransformationV1();
	const Vector3r &v2 = rb->getTransformationV2();

	PBD::DistanceFieldCollisionDetection::DistanceFieldCollisionObject *dfco = (PBD::DistanceFieldCollisionDetection::DistanceFieldCollisionObject *)co;
	const Vector3r &startX = co->m_aabb.m_p[0] - 0.1 * Vector3r::Ones();
	const Vector3r &endX = co->m_aabb.m_p[1] + 0.1 * Vector3r::Ones();
	Vector3r diff = endX - startX;
	const unsigned int steps = 20;
	Vector3r stepSize = (1.0 / steps) * diff;
	for (Real x = startX[0]; x < endX[0]; x += stepSize[0])
	{
		for (Real y = startX[1]; y < endX[1]; y += stepSize[1])
		{
			for (Real z = startX[2]; z < endX[2]; z += stepSize[2])
			{
				Vector3r pos_w(x, y, z);
				const Vector3r pos = R * (pos_w - com) + v1;
				double dist = dfco->distance(pos.template cast<double>(), 0.0);

				if (dist < 0.0)
				{
					if (dist < -1.0)
						dist = -1.0;
					float col[4] = { 1.0f + (float)dist, 0.0f, 0.0f, 1.0f };
					SPH::MiniGL::drawPoint(pos_w, 3.0f, col);
				}
				else if (dist < 1.0)
				{
					float col[4] = { 0.0f, 1.0f - (float)dist, 0.0f, 1.0f };
					SPH::MiniGL::drawPoint(pos_w, 3.0f, col);
				}
			}
		}
	}
}


void PBD_Simulator_GUI_TweakBar::renderTriangleModels()
{
	PBD::SimulationModel &model = m_pbdWrapper->getSimulationModel();
	const PBD::ParticleData &pd = model.getParticles();
	float surfaceColor[4] = { 0.8f, 0.9f, 0.2f, 1 };

	shaderBegin(surfaceColor);

	for (unsigned int i = 0; i < model.getTriangleModels().size(); i++)
	{
		// mesh 
		const Utilities::IndexedFaceMesh &mesh = model.getTriangleModels()[i]->getParticleMesh();
		const unsigned int offset = model.getTriangleModels()[i]->getIndexOffset();
		drawMesh(pd, mesh, offset, surfaceColor);
	}

	shaderEnd();
}

void PBD_Simulator_GUI_TweakBar::renderTetModels()
{
	PBD::SimulationModel &model = m_pbdWrapper->getSimulationModel();
	const PBD::ParticleData &pd = model.getParticles();
	float surfaceColor[4] = { 0.8f, 0.4f, 0.7f, 1 };

	shaderBegin(surfaceColor);

	for (unsigned int i = 0; i < model.getTetModels().size(); i++)
	{
		const PBD::VertexData &vdVis = model.getTetModels()[i]->getVisVertices();
		if (vdVis.size() > 0)
		{
			const Utilities::IndexedFaceMesh &visMesh = model.getTetModels()[i]->getVisMesh();
			drawMesh(vdVis, visMesh, 0, surfaceColor);
		}
		else
		{
			const Utilities::IndexedFaceMesh &surfaceMesh = model.getTetModels()[i]->getSurfaceMesh();
			const unsigned int offset = model.getTetModels()[i]->getIndexOffset();
			drawMesh(pd, surfaceMesh, offset, surfaceColor);
		}
	}

	shaderEnd();
}


void PBD_Simulator_GUI_TweakBar::renderBallJoint(PBD::BallJoint &bj)
{
	SPH::MiniGL::drawSphere(bj.m_jointInfo.col(2), 0.15f, m_jointColor);
}

void PBD_Simulator_GUI_TweakBar::renderRigidBodyParticleBallJoint(PBD::RigidBodyParticleBallJoint &bj)
{
	SPH::MiniGL::drawSphere(bj.m_jointInfo.col(1), 0.1f, m_jointColor);
}

void PBD_Simulator_GUI_TweakBar::renderBallOnLineJoint(PBD::BallOnLineJoint &bj)
{
	SPH::MiniGL::drawSphere(bj.m_jointInfo.col(5), 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(bj.m_jointInfo.col(5) - bj.m_jointInfo.col(7), bj.m_jointInfo.col(5) + bj.m_jointInfo.col(7), m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_TweakBar::renderHingeJoint(PBD::HingeJoint &joint)
{
	PBD::SimulationModel *model = PBD::Simulation::getCurrent()->getModel();
	const PBD::SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
	PBD::RigidBody *rb = rigidBodies[joint.m_bodies[0]];

	const Vector3r &c = joint.m_jointInfo.block<3, 1>(0, 4);
	const Vector3r &axis_local = joint.m_jointInfo.block<3, 1>(0, 6);
	const Vector3r axis = rb->getRotation().matrix() * axis_local;

	SPH::MiniGL::drawSphere(c - 0.5*axis, 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(c + 0.5*axis, 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(c - 0.5*axis, c + 0.5*axis, m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_TweakBar::renderUniversalJoint(PBD::UniversalJoint &uj)
{
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(4) - 0.5*uj.m_jointInfo.col(6), 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(4) + 0.5*uj.m_jointInfo.col(6), 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(5) - 0.5*uj.m_jointInfo.col(7), 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(5) + 0.5*uj.m_jointInfo.col(7), 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(uj.m_jointInfo.col(4) - 0.5*uj.m_jointInfo.col(6), uj.m_jointInfo.col(4) + 0.5*uj.m_jointInfo.col(6), m_jointColor, 0.05f);
	SPH::MiniGL::drawCylinder(uj.m_jointInfo.col(5) - 0.5*uj.m_jointInfo.col(7), uj.m_jointInfo.col(5) + 0.5*uj.m_jointInfo.col(7), m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_TweakBar::renderSliderJoint(PBD::SliderJoint &joint)
{
	PBD::SimulationModel *model = PBD::Simulation::getCurrent()->getModel();
	const PBD::SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
	PBD::RigidBody *rb = rigidBodies[joint.m_bodies[0]];

	Quaternionr qR0;
	qR0.coeffs() = joint.m_jointInfo.col(1);
	const Vector3r &c = rb->getPosition();
	Vector3r axis = qR0.matrix().col(0);
	SPH::MiniGL::drawSphere(c, 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(c - axis, c + axis, m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_TweakBar::renderTargetPositionMotorSliderJoint(PBD::TargetPositionMotorSliderJoint &joint)
{
	PBD::SimulationModel *model = PBD::Simulation::getCurrent()->getModel();
	const PBD::SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
	PBD::RigidBody *rb = rigidBodies[joint.m_bodies[0]];

	const Vector3r &c = rb->getPosition();
	Vector3r axis = joint.m_jointInfo.block<3, 1>(0, 1);
	SPH::MiniGL::drawSphere(c, 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(c - axis, c + axis, m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_TweakBar::renderTargetVelocityMotorSliderJoint(PBD::TargetVelocityMotorSliderJoint &joint)
{
	PBD::SimulationModel *model = PBD::Simulation::getCurrent()->getModel();
	const PBD::SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
	PBD::RigidBody *rb = rigidBodies[joint.m_bodies[0]];

	Quaternionr qR0;
	qR0.coeffs() = joint.m_jointInfo.col(1);
	const Vector3r &c = rb->getPosition();
	Vector3r axis = qR0.matrix().col(0);
	SPH::MiniGL::drawSphere(c, 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(c - axis, c + axis, m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_TweakBar::renderTargetAngleMotorHingeJoint(PBD::TargetAngleMotorHingeJoint &joint)
{
	PBD::SimulationModel *model = PBD::Simulation::getCurrent()->getModel();
	const PBD::SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
	PBD::RigidBody *rb = rigidBodies[joint.m_bodies[0]];

	const Vector3r &c = joint.m_jointInfo.block<3, 1>(0, 5);
	const Vector3r &axis_local = joint.m_jointInfo.block<3, 1>(0, 7);
	const Vector3r axis = rb->getRotation().matrix() * axis_local;

	SPH::MiniGL::drawSphere(c - 0.5*axis, 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(c + 0.5*axis, 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(c - 0.5*axis, c + 0.5*axis, m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_TweakBar::renderTargetVelocityMotorHingeJoint(PBD::TargetVelocityMotorHingeJoint &joint)
{
	PBD::SimulationModel *model = PBD::Simulation::getCurrent()->getModel();
	const PBD::SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
	PBD::RigidBody *rb = rigidBodies[joint.m_bodies[0]];

	const Vector3r &c = joint.m_jointInfo.block<3, 1>(0, 5);
	const Vector3r axis = joint.m_jointInfo.block<3, 1>(0, 7);

	SPH::MiniGL::drawSphere(c - 0.5*axis, 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(c + 0.5*axis, 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(c - 0.5*axis, c + 0.5*axis, m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_TweakBar::renderRigidBodyContact(PBD::RigidBodyContactConstraint &cc)
{
	float col1[4] = { 0.0f, 0.6f, 0.2f, 1 };
	float col2[4] = { 0.6f, 0.0f, 0.2f, 1 };
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(0), 5.0f, col1);
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(1), 5.0f, col2);
	SPH::MiniGL::drawVector(cc.m_constraintInfo.col(1), cc.m_constraintInfo.col(1) + cc.m_constraintInfo.col(2), 1.0f, col2);
}

void PBD_Simulator_GUI_TweakBar::renderParticleRigidBodyContact(PBD::ParticleRigidBodyContactConstraint &cc)
{
	float col1[4] = { 0.0f, 0.6f, 0.2f, 1 };
	float col2[4] = { 0.6f, 0.0f, 0.2f, 1 };
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(0), 5.0f, col1);
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(1), 5.0f, col2);
	SPH::MiniGL::drawVector(cc.m_constraintInfo.col(1), cc.m_constraintInfo.col(1) + cc.m_constraintInfo.col(2), 1.0f, col2);
}

void PBD_Simulator_GUI_TweakBar::renderSpring(PBD::RigidBodySpring &s)
{
	SPH::MiniGL::drawSphere(s.m_jointInfo.col(2), 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(s.m_jointInfo.col(3), 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(s.m_jointInfo.col(2), s.m_jointInfo.col(3), m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_TweakBar::renderDistanceJoint(PBD::DistanceJoint &j)
{
	SPH::MiniGL::drawSphere(j.m_jointInfo.col(2), 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(j.m_jointInfo.col(3), 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(j.m_jointInfo.col(2), j.m_jointInfo.col(3), m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_TweakBar::renderDamperJoint(PBD::DamperJoint &joint)
{
	PBD::SimulationModel *model = PBD::Simulation::getCurrent()->getModel();
	const PBD::SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
	PBD::RigidBody *rb = rigidBodies[joint.m_bodies[0]];

	Quaternionr qR0;
	qR0.coeffs() = joint.m_jointInfo.col(1);
	const Vector3r &c = rb->getPosition();
	Vector3r axis = qR0.matrix().col(0);
	SPH::MiniGL::drawSphere(c, 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(c - axis, c + axis, m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_TweakBar::renderConstraints()
{
	PBD::SimulationModel &model = m_pbdWrapper->getSimulationModel();
	PBD::SimulationModel::ConstraintVector &constraints = model.getConstraints();
	PBD::SimulationModel::RigidBodyContactConstraintVector &rigidBodyContacts = model.getRigidBodyContactConstraints();
	PBD::SimulationModel::ParticleRigidBodyContactConstraintVector &particleRigidBodyContacts = model.getParticleRigidBodyContactConstraints();

	for (unsigned int i = 0; i < rigidBodyContacts.size(); i++)
		renderRigidBodyContact(rigidBodyContacts[i]);
	for (unsigned int i = 0; i < particleRigidBodyContacts.size(); i++)
		renderParticleRigidBodyContact(particleRigidBodyContacts[i]);

	for (size_t i = 0; i < constraints.size(); i++)
	{
		if (constraints[i]->getTypeId() == PBD::BallJoint::TYPE_ID)
		{
			renderBallJoint(*(PBD::BallJoint*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::BallOnLineJoint::TYPE_ID)
		{
			renderBallOnLineJoint(*(PBD::BallOnLineJoint*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::HingeJoint::TYPE_ID)
		{
			renderHingeJoint(*(PBD::HingeJoint*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::UniversalJoint::TYPE_ID)
		{
			renderUniversalJoint(*(PBD::UniversalJoint*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::SliderJoint::TYPE_ID)
		{
			renderSliderJoint(*(PBD::SliderJoint*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::TargetAngleMotorHingeJoint::TYPE_ID)
		{
			renderTargetAngleMotorHingeJoint(*(PBD::TargetAngleMotorHingeJoint*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::TargetVelocityMotorHingeJoint::TYPE_ID)
		{
			renderTargetVelocityMotorHingeJoint(*(PBD::TargetVelocityMotorHingeJoint*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::TargetPositionMotorSliderJoint::TYPE_ID)
		{
			renderTargetPositionMotorSliderJoint(*(PBD::TargetPositionMotorSliderJoint*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::TargetVelocityMotorSliderJoint::TYPE_ID)
		{
			renderTargetVelocityMotorSliderJoint(*(PBD::TargetVelocityMotorSliderJoint*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::RigidBodyParticleBallJoint::TYPE_ID)
		{
			renderRigidBodyParticleBallJoint(*(PBD::RigidBodyParticleBallJoint*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::RigidBodySpring::TYPE_ID)
		{
			renderSpring(*(PBD::RigidBodySpring*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::DistanceJoint::TYPE_ID)
		{
			renderDistanceJoint(*(PBD::DistanceJoint*)constraints[i]);
		}
		else if (constraints[i]->getTypeId() == PBD::DamperJoint::TYPE_ID)
		{
			renderDamperJoint(*(PBD::DamperJoint*)constraints[i]);
		}
	}
}


void PBD_Simulator_GUI_TweakBar::initSimulationParameterGUI()
{
	Simulator_GUI_TweakBar::initSimulationParameterGUI();
	PBD::TimeStepController &timeStep = m_pbdWrapper->getTimeStepController();
	PBD::DistanceFieldCollisionDetection &cd = m_pbdWrapper->getCollisionDetection();
	PBD::SimulationModel &model = m_pbdWrapper->getSimulationModel();

	TwAddVarRW(getTweakBar(), "RenderAABBs", TW_TYPE_BOOLCPP, &m_drawAABB, " label='Render AABBs' group=PBD ");
	TwAddVarRW(getTweakBar(), "RenderBVH", TW_TYPE_INT32, &m_drawBVHDepth, " label='Render BVH depth' group=PBD ");
	TwAddVarRW(getTweakBar(), "RenderDistanceFields", TW_TYPE_BOOLCPP, &m_drawSDF, " label='Render distance fields' group=PBD ");
	TwAddVarCB(getTweakBar(), "DampingCoeff", TW_TYPE_REAL, setDampingCoeff, getDampingCoeff, m_pbdWrapper, " label='Damping' group=PBD ");
	TwType enumType = TwDefineEnum("VelocityUpdateMethodType", NULL, 0);
	TwAddVarCB(getTweakBar(), "VelocityUpdateMethod", enumType, setVelocityUpdateMethod, getVelocityUpdateMethod, &timeStep, " label='Velocity update method' enum='0 {First Order Update}, 1 {Second Order Update}' group=PBD");
	TwAddVarCB(getTweakBar(), "MaxIter", TW_TYPE_UINT32, setMaxIterations, getMaxIterations, &timeStep, " label='Max. iterations'  min=1 step=1 group=PBD ");
	TwAddVarCB(getTweakBar(), "MaxIterV", TW_TYPE_UINT32, setMaxIterationsV, getMaxIterationsV, &timeStep, " label='Max. iterations Vel.'  min=1 step=1 group=PBD ");
	TwAddVarCB(getTweakBar(), "ContactTolerance", TW_TYPE_REAL, setContactTolerance, getContactTolerance, &cd, " label='Contact tolerance'  step=0.001 precision=3 group=PBD ");
	TwAddVarCB(getTweakBar(), "ContactStiffnessRigidBody", TW_TYPE_REAL, setContactStiffnessRigidBody, getContactStiffnessRigidBody, &model, " label='Contact stiffness RB'  min=0.0 step=0.1 precision=2 group=PBD ");
	TwAddVarCB(getTweakBar(), "ContactStiffnessParticleRigidBody", TW_TYPE_REAL, setContactStiffnessParticleRigidBody, getContactStiffnessParticleRigidBody, &model, " label='Contact stiffness Particle-RB'  min=0.0 step=0.1 precision=2 group=PBD ");

	TwType enumType2 = TwDefineEnum("ClothSimulationMethodType", NULL, 0);
	TwAddVarCB(getTweakBar(), "ClothSimulationMethod", enumType2, setClothSimulationMethod, getClothSimulationMethod, m_pbdWrapper, " label='Cloth sim. method' enum='0 {None}, 1 {Distance constraints}, 2 {FEM based PBD}, 3 {Strain based dynamics}' group=PBD");
	TwType enumType3 = TwDefineEnum("SolidSimulationMethodType", NULL, 0);
	TwAddVarCB(getTweakBar(), "SolidSimulationMethod", enumType3, setSolidSimulationMethod, getSolidSimulationMethod, m_pbdWrapper,
		" label='Solid sim. method' enum='0 {None}, 1 {Volume constraints}, 2 {FEM based PBD}, 3 {Strain based dynamics (no inversion handling)}, 4 {Shape matching (no inversion handling)}' group=PBD");
	TwAddVarCB(getTweakBar(), "Stiffness", TW_TYPE_REAL, setStiffness, getStiffness, &model, " label='Stiffness'  min=0.0 step=0.1 precision=4 group='Distance constraints' ");
	TwAddVarCB(getTweakBar(), "XXStiffness", TW_TYPE_REAL, setXXStiffness, getXXStiffness, &model, " label='Stiffness XX'  min=0.0 step=0.1 precision=4 group='Strain based dynamics' ");
	TwAddVarCB(getTweakBar(), "YYStiffness", TW_TYPE_REAL, setYYStiffness, getYYStiffness, &model, " label='Stiffness YY'  min=0.0 step=0.1 precision=4 group='Strain based dynamics' ");
	TwAddVarCB(getTweakBar(), "XYStiffness", TW_TYPE_REAL, setXYStiffness, getXYStiffness, &model, " label='Stiffness XY'  min=0.0 step=0.1 precision=4 group='Strain based dynamics' ");
	TwAddVarCB(getTweakBar(), "XXStiffnessFEM", TW_TYPE_REAL, setXXStiffness, getXXStiffness, &model, " label='Youngs modulus XX'  min=0.0 step=0.1 precision=4 group='FEM based PBD' ");
	TwAddVarCB(getTweakBar(), "YYStiffnessFEM", TW_TYPE_REAL, setYYStiffness, getYYStiffness, &model, " label='Youngs modulus YY'  min=0.0 step=0.1 precision=4 group='FEM based PBD' ");
	TwAddVarCB(getTweakBar(), "XYStiffnessFEM", TW_TYPE_REAL, setXYStiffness, getXYStiffness, &model, " label='Youngs modulus XY'  min=0.0 step=0.1 precision=4 group='FEM based PBD' ");
	TwAddVarCB(getTweakBar(), "XYPoissonRatioFEM", TW_TYPE_REAL, setXYPoissonRatio, getXYPoissonRatio, &model, " label='Poisson ratio XY'  min=0.0 step=0.1 precision=4 group='FEM based PBD' ");
	TwAddVarCB(getTweakBar(), "YXPoissonRatioFEM", TW_TYPE_REAL, setYXPoissonRatio, getYXPoissonRatio, &model, " label='Poisson ratio YX'  min=0.0 step=0.1 precision=4 group='FEM based PBD' ");
	TwAddVarCB(getTweakBar(), "NormalizeStretch", TW_TYPE_BOOL32, setNormalizeStretch, getNormalizeStretch, &model, " label='Normalize stretch' group='Strain based dynamics' ");
	TwAddVarCB(getTweakBar(), "NormalizeShear", TW_TYPE_BOOL32, setNormalizeShear, getNormalizeShear, &model, " label='Normalize shear' group='Strain based dynamics' ");
	TwType enumType4 = TwDefineEnum("BendingMethodType", NULL, 0);
	TwAddVarCB(getTweakBar(), "BendingMethod", enumType4, setBendingMethod, getBendingMethod, m_pbdWrapper, " label='Bending method' enum='0 {None}, 1 {Dihedral angle}, 2 {Isometric bending}' group=Bending");
	TwAddVarCB(getTweakBar(), "BendingStiffness", TW_TYPE_REAL, setBendingStiffness, getBendingStiffness, &model, " label='Bending stiffness'  min=0.0 step=0.01 precision=4 group=Bending ");
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setVelocityUpdateMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	((PBD::TimeStepController*)clientData)->setValue(PBD::TimeStepController::VELOCITY_UPDATE_METHOD, (int)val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getVelocityUpdateMethod(void *value, void *clientData)
{
	*(short *)(value) = (short)((PBD::TimeStepController*)clientData)->getValue<int>(PBD::TimeStepController::VELOCITY_UPDATE_METHOD);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setDampingCoeff(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBDWrapper*)clientData)->setDampingCoeff(val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getDampingCoeff(void *value, void *clientData)
{
	*(Real *)(value) = ((PBDWrapper*)clientData)->getDampingCoeff();
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setMaxIterations(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	((PBD::TimeStepController*)clientData)->setValue(PBD::TimeStepController::MAX_ITERATIONS, val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getMaxIterations(void *value, void *clientData)
{
	*(unsigned int *)(value) = ((PBD::TimeStepController*)clientData)->getValue<unsigned int>(PBD::TimeStepController::MAX_ITERATIONS);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setMaxIterationsV(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	((PBD::TimeStepController*)clientData)->setValue(PBD::TimeStepController::MAX_ITERATIONS_V, val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getMaxIterationsV(void *value, void *clientData)
{
	*(unsigned int *)(value) = ((PBD::TimeStepController*)clientData)->getValue<unsigned int>(PBD::TimeStepController::MAX_ITERATIONS_V);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setContactTolerance(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::DistanceFieldCollisionDetection*)clientData)->setTolerance(val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getContactTolerance(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::DistanceFieldCollisionDetection*)clientData)->getTolerance();
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setContactStiffnessRigidBody(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setContactStiffnessRigidBody(val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getContactStiffnessRigidBody(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getContactStiffnessRigidBody();
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setContactStiffnessParticleRigidBody(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setContactStiffnessParticleRigidBody(val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getContactStiffnessParticleRigidBody(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getContactStiffnessParticleRigidBody();
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setValue(PBD::SimulationModel::CLOTH_STIFFNESS, val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getValue<Real>(PBD::SimulationModel::CLOTH_STIFFNESS);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setXXStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setValue(PBD::SimulationModel::CLOTH_STIFFNESS_XX, val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getXXStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getValue<Real>(PBD::SimulationModel::CLOTH_STIFFNESS_XX);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setYYStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setValue(PBD::SimulationModel::CLOTH_STIFFNESS_YY, val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getYYStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getValue<Real>(PBD::SimulationModel::CLOTH_STIFFNESS_YY);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setXYStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setValue(PBD::SimulationModel::CLOTH_STIFFNESS_XY, val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getXYStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getValue<Real>(PBD::SimulationModel::CLOTH_STIFFNESS_XY);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setYXPoissonRatio(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setValue(PBD::SimulationModel::CLOTH_POISSON_RATIO_YX, val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getYXPoissonRatio(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getValue<Real>(PBD::SimulationModel::CLOTH_POISSON_RATIO_YX);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setXYPoissonRatio(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setValue(PBD::SimulationModel::CLOTH_POISSON_RATIO_XY, val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getXYPoissonRatio(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getValue<Real>(PBD::SimulationModel::CLOTH_POISSON_RATIO_XY);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setNormalizeStretch(const void *value, void *clientData)
{
	const bool val = *(const bool *)(value);
	((PBD::SimulationModel*)clientData)->setValue(PBD::SimulationModel::CLOTH_NORMALIZE_STRETCH, val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getNormalizeStretch(void *value, void *clientData)
{
	*(bool *)(value) = ((PBD::SimulationModel*)clientData)->getValue<unsigned int>(PBD::SimulationModel::CLOTH_NORMALIZE_STRETCH);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setNormalizeShear(const void *value, void *clientData)
{
	const bool val = *(const bool *)(value);
	((PBD::SimulationModel*)clientData)->setValue(PBD::SimulationModel::CLOTH_NORMALIZE_SHEAR, val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getNormalizeShear(void *value, void *clientData)
{
	*(bool *)(value) = ((PBD::SimulationModel*)clientData)->getValue<unsigned int>(PBD::SimulationModel::CLOTH_NORMALIZE_SHEAR);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setBendingStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setValue(PBD::SimulationModel::CLOTH_BENDING_STIFFNESS, val);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getBendingStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getValue<Real>(PBD::SimulationModel::CLOTH_BENDING_STIFFNESS);
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setBendingMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	((PBDWrapper*)clientData)->setBendingMethod(val);
	((PBDWrapper*)clientData)->reset();
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getBendingMethod(void *value, void *clientData)
{
	*(short *)(value) = ((PBDWrapper*)clientData)->getBendingMethod();
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setClothSimulationMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	((PBDWrapper*)clientData)->setClothSimulationMethod(val);
	((PBDWrapper*)clientData)->reset();
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getClothSimulationMethod(void *value, void *clientData)
{
	*(short *)(value) = ((PBDWrapper*)clientData)->getClothSimulationMethod();
}

void TW_CALL PBD_Simulator_GUI_TweakBar::setSolidSimulationMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	((PBDWrapper*)clientData)->setSolidSimulationMethod(val);
	((PBDWrapper*)clientData)->reset();
}

void TW_CALL PBD_Simulator_GUI_TweakBar::getSolidSimulationMethod(void *value, void *clientData)
{
	*(short *)(value) = ((PBDWrapper*)clientData)->getSolidSimulationMethod();
}
