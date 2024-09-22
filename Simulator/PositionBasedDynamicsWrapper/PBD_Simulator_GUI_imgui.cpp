#include "PBD_Simulator_GUI_imgui.h"
#include "Utilities/FileSystem.h"
#include "Simulation/Simulation.h"
#include "Simulation/TimeStepController.h"
#include "GUI/imgui/imguiParameters.h"


using namespace SPH;
using namespace Utilities;

PBD_Simulator_GUI_imgui::PBD_Simulator_GUI_imgui(SimulatorBase *base, PBDWrapper *pbdWrapper) :
	Simulator_GUI_imgui(base)
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

PBD_Simulator_GUI_imgui::~PBD_Simulator_GUI_imgui(void)
{	
	delete m_shader;
}

void PBD_Simulator_GUI_imgui::init(const char *name)
{
	Simulator_GUI_imgui::init(name);
	initShader();
}

void PBD_Simulator_GUI_imgui::render()
{
	Simulator_GUI_imgui::render();

	renderTriangleModels();
	renderTetModels();
	renderConstraints();
	renderBVH();
	renderSDF();
}


SPH::Shader *PBD_Simulator_GUI_imgui::createShader(const std::string &vertexShader, const std::string &geometryShader, const std::string &fragmentShader)
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

void PBD_Simulator_GUI_imgui::initShader()
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


void PBD_Simulator_GUI_imgui::shaderBegin(const float *col)
{
	if (m_shader)
	{
		m_shader->begin();
		glUniform1f(m_shader->getUniform("shininess"), 5.0f);
		glUniform1f(m_shader->getUniform("specular_factor"), 0.2f);

		const Matrix4f matrix(SPH::MiniGL::getModelviewMatrix().cast<float>());
		glUniformMatrix4fv(m_shader->getUniform("modelview_matrix"), 1, GL_FALSE, &matrix(0,0));
		const Matrix4f pmatrix(SPH::MiniGL::getProjectionMatrix().cast<float>());
		glUniformMatrix4fv(m_shader->getUniform("projection_matrix"), 1, GL_FALSE, &pmatrix(0,0));
		glUniform3fv(m_shader->getUniform("surface_color"), 1, col);
	}
}

void PBD_Simulator_GUI_imgui::shaderEnd()
{
	if (m_shader)
		m_shader->end();
}

void PBD_Simulator_GUI_imgui::renderAABB(PBD::AABB &aabb)
{
	float color[4] = { 0.5f, 0.5f, 0.5f, 0.3f };

	Vector3r p1, p2;
	for (unsigned char i = 0; i < 12; i++)
	{
		PBD::AABB::getEdge(aabb, i, p1, p2);
		SPH::MiniGL::drawVector(p1, p2, 1.0f, color);
	}
}

void PBD_Simulator_GUI_imgui::renderBVH()
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


void PBD_Simulator_GUI_imgui::renderSDF()
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

void PBD_Simulator_GUI_imgui::renderSDF(PBD::CollisionDetection::CollisionObject* co)
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


void PBD_Simulator_GUI_imgui::renderTriangleModels()
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

void PBD_Simulator_GUI_imgui::renderTetModels()
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


void PBD_Simulator_GUI_imgui::renderBallJoint(PBD::BallJoint &bj)
{
	SPH::MiniGL::drawSphere(bj.m_jointInfo.col(2), 0.15f, m_jointColor);
}

void PBD_Simulator_GUI_imgui::renderRigidBodyParticleBallJoint(PBD::RigidBodyParticleBallJoint &bj)
{
	SPH::MiniGL::drawSphere(bj.m_jointInfo.col(1), 0.1f, m_jointColor);
}

void PBD_Simulator_GUI_imgui::renderBallOnLineJoint(PBD::BallOnLineJoint &bj)
{
	SPH::MiniGL::drawSphere(bj.m_jointInfo.col(5), 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(bj.m_jointInfo.col(5) - bj.m_jointInfo.col(7), bj.m_jointInfo.col(5) + bj.m_jointInfo.col(7), m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_imgui::renderHingeJoint(PBD::HingeJoint &joint)
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

void PBD_Simulator_GUI_imgui::renderUniversalJoint(PBD::UniversalJoint &uj)
{
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(4) - 0.5*uj.m_jointInfo.col(6), 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(4) + 0.5*uj.m_jointInfo.col(6), 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(5) - 0.5*uj.m_jointInfo.col(7), 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(5) + 0.5*uj.m_jointInfo.col(7), 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(uj.m_jointInfo.col(4) - 0.5*uj.m_jointInfo.col(6), uj.m_jointInfo.col(4) + 0.5*uj.m_jointInfo.col(6), m_jointColor, 0.05f);
	SPH::MiniGL::drawCylinder(uj.m_jointInfo.col(5) - 0.5*uj.m_jointInfo.col(7), uj.m_jointInfo.col(5) + 0.5*uj.m_jointInfo.col(7), m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_imgui::renderSliderJoint(PBD::SliderJoint &joint)
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

void PBD_Simulator_GUI_imgui::renderTargetPositionMotorSliderJoint(PBD::TargetPositionMotorSliderJoint &joint)
{
	PBD::SimulationModel *model = PBD::Simulation::getCurrent()->getModel();
	const PBD::SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
	PBD::RigidBody *rb = rigidBodies[joint.m_bodies[0]];

	const Vector3r &c = rb->getPosition();
	Vector3r axis = joint.m_jointInfo.block<3, 1>(0, 1);
	SPH::MiniGL::drawSphere(c, 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(c - axis, c + axis, m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_imgui::renderTargetVelocityMotorSliderJoint(PBD::TargetVelocityMotorSliderJoint &joint)
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

void PBD_Simulator_GUI_imgui::renderTargetAngleMotorHingeJoint(PBD::TargetAngleMotorHingeJoint &joint)
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

void PBD_Simulator_GUI_imgui::renderTargetVelocityMotorHingeJoint(PBD::TargetVelocityMotorHingeJoint &joint)
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

void PBD_Simulator_GUI_imgui::renderRigidBodyContact(PBD::RigidBodyContactConstraint &cc)
{
	float col1[4] = { 0.0f, 0.6f, 0.2f, 1 };
	float col2[4] = { 0.6f, 0.0f, 0.2f, 1 };
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(0), 5.0f, col1);
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(1), 5.0f, col2);
	SPH::MiniGL::drawVector(cc.m_constraintInfo.col(1), cc.m_constraintInfo.col(1) + cc.m_constraintInfo.col(2), 1.0f, col2);
}

void PBD_Simulator_GUI_imgui::renderParticleRigidBodyContact(PBD::ParticleRigidBodyContactConstraint &cc)
{
	float col1[4] = { 0.0f, 0.6f, 0.2f, 1 };
	float col2[4] = { 0.6f, 0.0f, 0.2f, 1 };
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(0), 5.0f, col1);
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(1), 5.0f, col2);
	SPH::MiniGL::drawVector(cc.m_constraintInfo.col(1), cc.m_constraintInfo.col(1) + cc.m_constraintInfo.col(2), 1.0f, col2);
}

void PBD_Simulator_GUI_imgui::renderSpring(PBD::RigidBodySpring &s)
{
	SPH::MiniGL::drawSphere(s.m_jointInfo.col(2), 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(s.m_jointInfo.col(3), 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(s.m_jointInfo.col(2), s.m_jointInfo.col(3), m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_imgui::renderDistanceJoint(PBD::DistanceJoint &j)
{
	SPH::MiniGL::drawSphere(j.m_jointInfo.col(2), 0.1f, m_jointColor);
	SPH::MiniGL::drawSphere(j.m_jointInfo.col(3), 0.1f, m_jointColor);
	SPH::MiniGL::drawCylinder(j.m_jointInfo.col(2), j.m_jointInfo.col(3), m_jointColor, 0.05f);
}

void PBD_Simulator_GUI_imgui::renderDamperJoint(PBD::DamperJoint &joint)
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

void PBD_Simulator_GUI_imgui::renderConstraints()
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


void PBD_Simulator_GUI_imgui::initSimulationParameterGUI()
{
	Simulator_GUI_imgui::initSimulationParameterGUI();
	PBD::TimeStepController &timeStep = m_pbdWrapper->getTimeStepController();
	PBD::DistanceFieldCollisionDetection &cd = m_pbdWrapper->getCollisionDetection();
	PBD::SimulationModel &model = m_pbdWrapper->getSimulationModel();

	imguiParameters::imguiBoolParameter* bparam = new imguiParameters::imguiBoolParameter();
	bparam->description = "Render AABBs";
	bparam->label = "Render AABBs";
	bparam->getFct = [this]() -> bool { return m_drawAABB; };
	bparam->setFct = [this](bool v) { m_drawAABB = v; };
	imguiParameters::addParam("PBD", "Visualization", bparam);

	imguiParameters::imguiNumericParameter<int> *iparam = new imguiParameters::imguiNumericParameter<int>();
	iparam->description = "Render BVH depth";
	iparam->label = "Render BVH depth";
	iparam->getFct = [this]() -> int { return m_drawBVHDepth; };
	iparam->setFct = [this](int v) { m_drawBVHDepth = v; };
	imguiParameters::addParam("PBD", "Visualization", iparam);

	bparam = new imguiParameters::imguiBoolParameter();
	bparam->description = "Render distance fields";
	bparam->label = "Render distance fields";
	bparam->getFct = [this]() -> bool { return m_drawSDF; };
	bparam->setFct = [this](bool v) { m_drawSDF = v; };
	imguiParameters::addParam("PBD", "Visualization", bparam);

	imguiParameters::imguiNumericParameter<Real>* rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Damping coefficient";
	rparam->label = "Damping";
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getDampingCoeff(); };
	rparam->setFct = [this](Real v) { m_pbdWrapper->setDampingCoeff(v); };
	imguiParameters::addParam("PBD", "Simulation", rparam);

	imguiParameters::imguiEnumParameter* eparam = new imguiParameters::imguiEnumParameter();
	eparam->description = "Velocity update method";
	eparam->label = "Velocity update method";
	eparam->items.push_back("First Order Update");
	eparam->items.push_back("Second Order Update");
	eparam->getFct = [this]() -> int { return m_pbdWrapper->getTimeStepController().getValue<int>(PBD::TimeStepController::VELOCITY_UPDATE_METHOD); };
	eparam->setFct = [this](int v) { m_pbdWrapper->getTimeStepController().setValue(PBD::TimeStepController::VELOCITY_UPDATE_METHOD, v); };
	imguiParameters::addParam("PBD", "Simulation", eparam);

	imguiParameters::imguiNumericParameter<int>* uparam = new imguiParameters::imguiNumericParameter<int>();
	uparam->description = "Number of sub steps of the solver.";
	uparam->label = "# sub steps";
	uparam->minValue = 1;
	uparam->getFct = [this]() -> unsigned int { return m_pbdWrapper->getTimeStepController().getValue<int>(PBD::TimeStepController::NUM_SUB_STEPS); };
	uparam->setFct = [this](unsigned int v) { m_pbdWrapper->getTimeStepController().setValue(PBD::TimeStepController::NUM_SUB_STEPS, v); };
	imguiParameters::addParam("PBD", "Simulation", uparam);

	uparam = new imguiParameters::imguiNumericParameter<int>();
	uparam->description = "Max. iterations";
	uparam->label = "Max. iterations";
	uparam->minValue = 1;
	uparam->getFct = [this]() -> unsigned int { return m_pbdWrapper->getTimeStepController().getValue<int>(PBD::TimeStepController::MAX_ITERATIONS); };
	uparam->setFct = [this](unsigned int v) { m_pbdWrapper->getTimeStepController().setValue(PBD::TimeStepController::MAX_ITERATIONS, v); };
	imguiParameters::addParam("PBD", "Simulation", uparam);

	uparam = new imguiParameters::imguiNumericParameter<int>();
	uparam->description = "Max. iterations Vel.";
	uparam->label = "Max. iterations Vel.";
	uparam->minValue = 1;
	uparam->getFct = [this]() -> unsigned int { return m_pbdWrapper->getTimeStepController().getValue<int>(PBD::TimeStepController::MAX_ITERATIONS_V); };
	uparam->setFct = [this](unsigned int v) { m_pbdWrapper->getTimeStepController().setValue(PBD::TimeStepController::MAX_ITERATIONS_V, v); };
	imguiParameters::addParam("PBD", "Simulation", uparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Contact tolerance";
	rparam->label = "Contact tolerance";
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getCollisionDetection().getTolerance(); };
	rparam->setFct = [this](Real v) { m_pbdWrapper->getCollisionDetection().setTolerance(v); };
	imguiParameters::addParam("PBD", "Simulation", rparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Contact stiffness RB";
	rparam->label = "Contact stiffness RB";
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getSimulationModel().getContactStiffnessRigidBody(); };
	rparam->setFct = [this](Real v) { m_pbdWrapper->getSimulationModel().setContactStiffnessRigidBody(v); };
	imguiParameters::addParam("PBD", "Simulation", rparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Contact stiffness Particle-RB";
	rparam->label = "Contact stiffness Particle-RB";
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getSimulationModel().getContactStiffnessParticleRigidBody(); };
	rparam->setFct = [this](Real v) { m_pbdWrapper->getSimulationModel().setContactStiffnessParticleRigidBody(v); };
	imguiParameters::addParam("PBD", "Simulation", rparam);
	
	eparam = new imguiParameters::imguiEnumParameter();
	eparam->description = "Cloth sim. method";
	eparam->label = "Cloth sim. method";
	eparam->items.push_back("None");
	eparam->items.push_back("Distance constraints");
	eparam->items.push_back("FEM based PBD");
	eparam->items.push_back("Strain based dynamics");
	eparam->getFct = [this]() -> int { return m_pbdWrapper->getClothSimulationMethod(); };
	eparam->setFct = [this](int v) { m_pbdWrapper->setClothSimulationMethod(v); m_pbdWrapper->reset(); };
	imguiParameters::addParam("PBD", "Simulation", eparam);

	eparam = new imguiParameters::imguiEnumParameter();
	eparam->description = "Solid sim. method";
	eparam->label = "Solid sim. method";
	eparam->items.push_back("None");
	eparam->items.push_back("Volume constraints");
	eparam->items.push_back("FEM based PBD");
	eparam->items.push_back("Strain based dynamics (no inversion handling)");
	eparam->items.push_back("Shape matching (no inversion handling)");
	eparam->getFct = [this]() -> int { return m_pbdWrapper->getSolidSimulationMethod(); };
	eparam->setFct = [this](int v) { m_pbdWrapper->setSolidSimulationMethod(v); m_pbdWrapper->reset(); };
	imguiParameters::addParam("PBD", "Simulation", eparam);

	eparam = new imguiParameters::imguiEnumParameter();
	eparam->description = "Bending method";
	eparam->label = "Bending method";
	eparam->items.push_back("None");
	eparam->items.push_back("Dihedral angle");
	eparam->items.push_back("Isometric bending");
	eparam->getFct = [this]() -> int { return m_pbdWrapper->getBendingMethod(); };
	eparam->setFct = [this](int v) { m_pbdWrapper->setBendingMethod(v); m_pbdWrapper->reset(); };
	imguiParameters::addParam("PBD", "Bending", eparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Bending stiffness";
	rparam->label = "Bending stiffness";
	rparam->minValue = 0.0;
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getBendingStiffness(); };
	rparam->setFct = [this](Real v) { 
		m_pbdWrapper->setBendingStiffness(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::DihedralConstraint, Real, &PBD::DihedralConstraint::m_stiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::IsometricBendingConstraint, Real, &PBD::IsometricBendingConstraint::m_stiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::IsometricBendingConstraint_XPBD, Real, &PBD::IsometricBendingConstraint_XPBD::m_stiffness>(v);
	};
	imguiParameters::addParam("PBD", "Bending", rparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Distance constraint stiffness";
	rparam->label = "Distance constraint stiffness";
	rparam->minValue = 0.0;
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getDistanceStiffness(); };
	rparam->setFct = [this](Real v) {
		m_pbdWrapper->setDistanceStiffness(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::DistanceConstraint, Real, &PBD::DistanceConstraint::m_stiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::DistanceConstraint_XPBD, Real, &PBD::DistanceConstraint_XPBD::m_stiffness>(v);
	};
	imguiParameters::addParam("PBD", "Cloth", rparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Youngs modulus XX";
	rparam->label = "Youngs modulus XX";
	rparam->minValue = 0.0;
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getXXStiffness(); };
	rparam->setFct = [this](Real v) {
		m_pbdWrapper->setXXStiffness(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::FEMTriangleConstraint, Real, &PBD::FEMTriangleConstraint::m_xxStiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::StrainTriangleConstraint, Real, &PBD::StrainTriangleConstraint::m_xxStiffness>(v);
	};
	imguiParameters::addParam("PBD", "Cloth", rparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Youngs modulus YY";
	rparam->label = "Youngs modulus YY";
	rparam->minValue = 0.0;
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getYYStiffness(); };
	rparam->setFct = [this](Real v) {
		m_pbdWrapper->setYYStiffness(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::FEMTriangleConstraint, Real, &PBD::FEMTriangleConstraint::m_yyStiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::StrainTriangleConstraint, Real, &PBD::StrainTriangleConstraint::m_yyStiffness>(v);
	};
	imguiParameters::addParam("PBD", "Cloth", rparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Youngs modulus XY";
	rparam->label = "Youngs modulus XY";
	rparam->minValue = 0.0;
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getXYStiffness(); };
	rparam->setFct = [this](Real v) {
		m_pbdWrapper->setXYStiffness(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::FEMTriangleConstraint, Real, &PBD::FEMTriangleConstraint::m_xyStiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::StrainTriangleConstraint, Real, &PBD::StrainTriangleConstraint::m_xyStiffness>(v);
	};
	imguiParameters::addParam("PBD", "Cloth", rparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Poisson ratio XY";
	rparam->label = "Poisson ratio XY";
	rparam->minValue = 0.0;
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getXYPoissonRatio(); };
	rparam->setFct = [this](Real v) {
		m_pbdWrapper->setXYPoissonRatio(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::FEMTriangleConstraint, Real, &PBD::FEMTriangleConstraint::m_xyPoissonRatio>(v);
	};
	imguiParameters::addParam("PBD", "Cloth", rparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Poisson ratio YX";
	rparam->label = "Poisson ratio YX";
	rparam->minValue = 0.0;
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getYXPoissonRatio(); };
	rparam->setFct = [this](Real v) {
		m_pbdWrapper->setYXPoissonRatio(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::FEMTriangleConstraint, Real, &PBD::FEMTriangleConstraint::m_yxPoissonRatio>(v);
	};
	imguiParameters::addParam("PBD", "Cloth", rparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Stiffness";
	rparam->label = "Stiffness";
	rparam->minValue = 0.0;
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getSolidStiffness(); };
	rparam->setFct = [this](Real v) {
		m_pbdWrapper->setSolidStiffness(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::FEMTetConstraint, Real, &PBD::FEMTetConstraint::m_stiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::StrainTetConstraint, Real, &PBD::StrainTetConstraint::m_stretchStiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::StrainTetConstraint, Real, &PBD::StrainTetConstraint::m_shearStiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::DistanceConstraint, Real, &PBD::DistanceConstraint::m_stiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::DistanceConstraint_XPBD, Real, &PBD::DistanceConstraint_XPBD::m_stiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::ShapeMatchingConstraint, Real, &PBD::ShapeMatchingConstraint::m_stiffness>(v);
	};
	imguiParameters::addParam("PBD", "Solid", rparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Poisson ratio";
	rparam->label = "Poisson ratio";
	rparam->minValue = 0.0;
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getSolidPoissonRatio(); };
	rparam->setFct = [this](Real v) {
		m_pbdWrapper->setSolidPoissonRatio(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::FEMTetConstraint, Real, &PBD::FEMTetConstraint::m_poissonRatio>(v);
	};
	imguiParameters::addParam("PBD", "Solid", rparam);

	rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Volume stiffness";
	rparam->label = "Volume stiffness";
	rparam->minValue = 0.0;
	rparam->getFct = [this]() -> Real { return m_pbdWrapper->getVolumeStiffness(); };
	rparam->setFct = [this](Real v) {
		m_pbdWrapper->setVolumeStiffness(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::VolumeConstraint, Real, &PBD::VolumeConstraint::m_stiffness>(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::VolumeConstraint_XPBD, Real, &PBD::VolumeConstraint_XPBD::m_stiffness>(v);
	};
	imguiParameters::addParam("PBD", "Solid", rparam);

	bparam = new imguiParameters::imguiBoolParameter();
	bparam->description = "Normalize stretch";
	bparam->label = "Normalize stretch";
	bparam->getFct = [this]() -> bool { return m_pbdWrapper->getNormalizeStretch(); };
	bparam->setFct = [this](bool v) {
		m_pbdWrapper->setNormalizeStretch(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::StrainTriangleConstraint, bool, &PBD::StrainTriangleConstraint::m_normalizeStretch>(v);
	};
	imguiParameters::addParam("PBD", "Cloth - Strain based dynamics", bparam);

	bparam = new imguiParameters::imguiBoolParameter();
	bparam->description = "Normalize shear";
	bparam->label = "Normalize shear";
	bparam->getFct = [this]() -> bool { return m_pbdWrapper->getNormalizeShear(); };
	bparam->setFct = [this](bool v) {
		m_pbdWrapper->setNormalizeShear(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::StrainTriangleConstraint, bool, &PBD::StrainTriangleConstraint::m_normalizeShear>(v);
	};
	imguiParameters::addParam("PBD", "Cloth - Strain based dynamics", bparam);

	bparam = new imguiParameters::imguiBoolParameter();
	bparam->description = "Normalize stretch";
	bparam->label = "Normalize stretch";
	bparam->getFct = [this]() -> bool { return m_pbdWrapper->getSolidNormalizeStretch(); };
	bparam->setFct = [this](bool v) {
		m_pbdWrapper->setSolidNormalizeStretch(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::StrainTetConstraint, bool, &PBD::StrainTetConstraint::m_normalizeStretch>(v);
	};
	imguiParameters::addParam("PBD", "Solid - Strain based dynamics", bparam);

	bparam = new imguiParameters::imguiBoolParameter();
	bparam->description = "Normalize shear";
	bparam->label = "Normalize shear";
	bparam->getFct = [this]() -> bool { return m_pbdWrapper->getSolidNormalizeShear(); };
	bparam->setFct = [this](bool v) {
		m_pbdWrapper->setSolidNormalizeShear(v);
		m_pbdWrapper->getSimulationModel().setConstraintValue<PBD::StrainTetConstraint, bool, &PBD::StrainTetConstraint::m_normalizeShear>(v);
	};
	imguiParameters::addParam("PBD", "Solid - Strain based dynamics", bparam);
}
