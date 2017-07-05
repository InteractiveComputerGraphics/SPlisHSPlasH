#include "PBDWrapper.h"
#include "GL/glew.h"
#include "GL/glut.h"
#include "Common/Common.h"
#include <iostream>
#include "Demos/Utils/OBJLoader.h"
#include "Demos/Utils/Utilities.h"
#include "Demos/Utils/SceneLoader.h"
#include "Demos/Utils/TetGenLoader.h"
#include "Demos/Utils/Utilities.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/Utilities/Timing.h"
#include "GL/freeglut_ext.h"
#include "Demos/Visualization/Shader.h"
#include "Demos/Simulation/TimeManager.h"
#include "Demos/Visualization/Visualization.h"

#define _USE_MATH_DEFINES
#include "math.h"


using namespace Eigen;
using namespace std;


float jointColor[4] = { 0.0f, 0.6f, 0.2f, 1 };

PBDWrapper::PBDWrapper()
{
	m_dataPath = std::string(PBD_DATA_PATH);

	m_clothSimulationMethod = 2;
	m_solidSimulationMethod = 2;
	m_bendingMethod = 2;
	m_sceneName = "";
	m_sceneFileName = "";
	m_enableMayaExport = false;
	m_dampingCoeff = 0.0;
	m_drawAABB = false;
	m_drawStaticBodies = true;
	m_drawDistanceFields = -1;
	m_drawBVHDepth = -1;
}

PBDWrapper::~PBDWrapper()
{
	delete m_shader;
	delete PBD::TimeManager::getCurrent();
}

void PBDWrapper::initGUI()
{
	TwAddVarRW(SPH::MiniGL::getTweakBar(), "RenderAABBs", TW_TYPE_BOOLCPP, &m_drawAABB, " label='Render AABBs' group=PBD ");
	TwAddVarRW(SPH::MiniGL::getTweakBar(), "RenderBVH", TW_TYPE_INT32, &m_drawBVHDepth, " label='Render BVH depth' group=PBD ");
	TwAddVarRW(SPH::MiniGL::getTweakBar(), "RenderStaticBodies", TW_TYPE_BOOLCPP, &m_drawStaticBodies, " label='Render static bodies' group=PBD ");
	TwAddVarRW(SPH::MiniGL::getTweakBar(), "RenderDistanceFields", TW_TYPE_INT32, &m_drawDistanceFields, " label='Render distance fields' step=1 group=PBD ");
	TwAddVarRW(SPH::MiniGL::getTweakBar(), "DampingCoeff", TW_TYPE_REAL, &m_dampingCoeff, " label='Damping' group=PBD ");
	TwType enumType = TwDefineEnum("VelocityUpdateMethodType", NULL, 0);
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "VelocityUpdateMethod", enumType, setVelocityUpdateMethod, getVelocityUpdateMethod, &m_sim, " label='Velocity update method' enum='0 {First Order Update}, 1 {Second Order Update}' group=PBD");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "MaxIter", TW_TYPE_UINT32, setMaxIterations, getMaxIterations, &m_sim, " label='Max. iterations'  min=1 step=1 group=PBD ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "MaxIterV", TW_TYPE_UINT32, setMaxIterationsV, getMaxIterationsV, &m_sim, " label='Max. iterations Vel.'  min=1 step=1 group=PBD ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "ContactTolerance", TW_TYPE_REAL, setContactTolerance, getContactTolerance, &m_cd, " label='Contact tolerance'  step=0.001 precision=3 group=PBD ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "ContactStiffnessRigidBody", TW_TYPE_REAL, setContactStiffnessRigidBody, getContactStiffnessRigidBody, &m_model, " label='Contact stiffness RB'  min=0.0 step=0.1 precision=2 group=PBD ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "ContactStiffnessParticleRigidBody", TW_TYPE_REAL, setContactStiffnessParticleRigidBody, getContactStiffnessParticleRigidBody, &m_model, " label='Contact stiffness Particle-RB'  min=0.0 step=0.1 precision=2 group=PBD ");
	TwType enumType2 = TwDefineEnum("ClothSimulationMethodType", NULL, 0);
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "ClothSimulationMethod", enumType2, setClothSimulationMethod, getClothSimulationMethod, this, " label='Cloth sim. method' enum='0 {None}, 1 {Distance constraints}, 2 {FEM based PBD}, 3 {Strain based dynamics}' group=PBD");
	TwType enumType3 = TwDefineEnum("SolidSimulationMethodType", NULL, 0);
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "SolidSimulationMethod", enumType3, setSolidSimulationMethod, getSolidSimulationMethod, this,
		" label='Solid sim. method' enum='0 {None}, 1 {Volume constraints}, 2 {FEM based PBD}, 3 {Strain based dynamics (no inversion handling)}, 4 {Shape matching (no inversion handling)}' group=PBD");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "Stiffness", TW_TYPE_REAL, setStiffness, getStiffness, &m_model, " label='Stiffness'  min=0.0 step=0.1 precision=4 group='Distance constraints' ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "XXStiffness", TW_TYPE_REAL, setXXStiffness, getXXStiffness, &m_model, " label='Stiffness XX'  min=0.0 step=0.1 precision=4 group='Strain based dynamics' ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "YYStiffness", TW_TYPE_REAL, setYYStiffness, getYYStiffness, &m_model, " label='Stiffness YY'  min=0.0 step=0.1 precision=4 group='Strain based dynamics' ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "XYStiffness", TW_TYPE_REAL, setXYStiffness, getXYStiffness, &m_model, " label='Stiffness XY'  min=0.0 step=0.1 precision=4 group='Strain based dynamics' ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "XXStiffnessFEM", TW_TYPE_REAL, setXXStiffness, getXXStiffness, &m_model, " label='Youngs modulus XX'  min=0.0 step=0.1 precision=4 group='FEM based PBD' ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "YYStiffnessFEM", TW_TYPE_REAL, setYYStiffness, getYYStiffness, &m_model, " label='Youngs modulus YY'  min=0.0 step=0.1 precision=4 group='FEM based PBD' ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "XYStiffnessFEM", TW_TYPE_REAL, setXYStiffness, getXYStiffness, &m_model, " label='Youngs modulus XY'  min=0.0 step=0.1 precision=4 group='FEM based PBD' ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "XYPoissonRatioFEM", TW_TYPE_REAL, setXYPoissonRatio, getXYPoissonRatio, &m_model, " label='Poisson ratio XY'  min=0.0 step=0.1 precision=4 group='FEM based PBD' ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "YXPoissonRatioFEM", TW_TYPE_REAL, setYXPoissonRatio, getYXPoissonRatio, &m_model, " label='Poisson ratio YX'  min=0.0 step=0.1 precision=4 group='FEM based PBD' ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "NormalizeStretch", TW_TYPE_BOOL32, setNormalizeStretch, getNormalizeStretch, &m_model, " label='Normalize stretch' group='Strain based dynamics' ");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "NormalizeShear", TW_TYPE_BOOL32, setNormalizeShear, getNormalizeShear, &m_model, " label='Normalize shear' group='Strain based dynamics' ");
	TwType enumType4 = TwDefineEnum("BendingMethodType", NULL, 0);
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "BendingMethod", enumType4, setBendingMethod, getBendingMethod, this, " label='Bending method' enum='0 {None}, 1 {Dihedral angle}, 2 {Isometric bending}' group=Bending");
	TwAddVarCB(SPH::MiniGL::getTweakBar(), "BendingStiffness", TW_TYPE_REAL, setBendingStiffness, getBendingStiffness, &m_model, " label='Bending stiffness'  min=0.0 step=0.01 precision=4 group=Bending ");
}

SPH::Shader *PBDWrapper::createShader(const std::string &vertexShader, const std::string &geometryShader, const std::string &fragmentShader)
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

void PBDWrapper::initShader()
{
	std::string vertFile = PBD::Utilities::normalizePath(m_dataPath + "/shaders/vs_smooth.glsl");
	std::string fragFile = PBD::Utilities::normalizePath(m_dataPath + "/shaders/fs_smooth.glsl");
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

 void PBDWrapper::reset()
 {
	m_model.reset();
	m_sim.reset();
 }
 

 
 void PBDWrapper::timeStep()
 {
	PBD::ParticleData &pd = m_model.getParticles();
	PBD::SimulationModel::RigidBodyVector &rb = m_model.getRigidBodies();
	PBD::TimeManager::getCurrent()->setTimeStepSize(SPH::TimeManager::getCurrent()->getTimeStepSize());

	m_sim.step(m_model);

	for (unsigned int i = 0; i < pd.size(); i++)
	{
		pd.getVelocity(i) *= (1.0 - m_dampingCoeff);
	}
	for (unsigned int i = 0; i < rb.size(); i++)
	{
		rb[i]->getVelocity() *= (1.0 - m_dampingCoeff);
		rb[i]->getAngularVelocity() *= (1.0 - m_dampingCoeff);
	}
}
 
void PBDWrapper::updateVisModels()
{
	PBD::ParticleData &pd = m_model.getParticles();

	// Update visualization models	
	for (unsigned int i = 0; i < m_model.getTetModels().size(); i++)
	{
		m_model.getTetModels()[i]->updateMeshNormals(pd);
		m_model.getTetModels()[i]->updateVisMesh(pd);
	}
	for (unsigned int i = 0; i < m_model.getTriangleModels().size(); i++)
	{
		m_model.getTriangleModels()[i]->updateMeshNormals(pd);
	}
}

void PBDWrapper::readScene(const std::string &sceneFileName)
{
	m_sceneFileName = sceneFileName;

	PBD::SimulationModel::RigidBodyVector &rb = m_model.getRigidBodies();
	PBD::SimulationModel::TriangleModelVector &triModels = m_model.getTriangleModels();
	PBD::SimulationModel::TetModelVector &tetModels = m_model.getTetModels();
	PBD::SimulationModel::ConstraintVector &constraints = m_model.getConstraints();

	PBD::SceneLoader::SceneData data;
	PBD::SceneLoader loader;
	loader.readScene(sceneFileName, data);
	std::cout << "Scene: " << sceneFileName << "\n";

	std::string basePath = PBD::Utilities::getFilePath(sceneFileName);
	
	m_sceneName = data.m_sceneName;

	m_sim.setGravity(data.m_gravity);
	m_sim.setMaxIterations(data.m_maxIter);
	m_sim.setMaxIterationsV(data.m_maxIterVel);
	m_sim.setVelocityUpdateMethod(data.m_velocityUpdateMethod);
	if (data.m_triangleModelSimulationMethod != -1)
		m_clothSimulationMethod = data.m_triangleModelSimulationMethod;
	if (data.m_tetModelSimulationMethod != -1)
		m_solidSimulationMethod = data.m_tetModelSimulationMethod;
	if (data.m_triangleModelBendingMethod != -1)
		m_bendingMethod = data.m_triangleModelBendingMethod;
	m_cd.setTolerance(data.m_contactTolerance);
	m_model.setContactStiffnessRigidBody(data.m_contactStiffnessRigidBody);
	m_model.setContactStiffnessParticleRigidBody(data.m_contactStiffnessParticleRigidBody);

	m_model.setClothStiffness(data.m_cloth_stiffness);
	m_model.setClothBendingStiffness(data.m_cloth_bendingStiffness);
	m_model.setClothXXStiffness(data.m_cloth_xxStiffness);
	m_model.setClothYYStiffness(data.m_cloth_yyStiffness);
	m_model.setClothXYStiffness(data.m_cloth_xyStiffness);
	m_model.setClothXYPoissonRatio(data.m_cloth_xyPoissonRatio);
	m_model.setClothYXPoissonRatio(data.m_cloth_yxPoissonRatio);
	m_model.setClothNormalizeStretch(data.m_cloth_normalizeStretch);
	m_model.setClothNormalizeShear(data.m_cloth_normalizeShear);

	//////////////////////////////////////////////////////////////////////////
	// rigid bodies
	//////////////////////////////////////////////////////////////////////////

	// map file names to loaded geometry to prevent multiple imports of same files
	std::map<std::string, pair<PBD::VertexData, PBD::IndexedFaceMesh>> objFiles;
	for (unsigned int i = 0; i < data.m_rigidBodyData.size(); i++)
	{
		const PBD::SceneLoader::RigidBodyData &rbd = data.m_rigidBodyData[i];

		// Check if already loaded
		if (objFiles.find(rbd.m_modelFile) == objFiles.end())
		{
			PBD::IndexedFaceMesh mesh;
			PBD::VertexData vd;
			PBD::OBJLoader::loadObj(PBD::Utilities::normalizePath(rbd.m_modelFile), vd, mesh);
			objFiles[rbd.m_modelFile] = { vd, mesh };
		}
	}

	for (unsigned int i = 0; i < data.m_tetModelData.size(); i++)
	{
		const PBD::SceneLoader::TetModelData &tmd = data.m_tetModelData[i];

		// Check if already loaded
		if ((tmd.m_modelFileVis != "") &&
			(objFiles.find(tmd.m_modelFileVis) == objFiles.end()))
		{
			PBD::IndexedFaceMesh mesh;
			PBD::VertexData vd;
			PBD::OBJLoader::loadObj(PBD::Utilities::normalizePath(tmd.m_modelFileVis), vd, mesh);
			objFiles[tmd.m_modelFileVis] = { vd, mesh };
		}
	}


	rb.resize(data.m_rigidBodyData.size());
	std::map<unsigned int, unsigned int> id_index;
	for (unsigned int i = 0; i < data.m_rigidBodyData.size(); i++)
	{
		const PBD::SceneLoader::RigidBodyData &rbd = data.m_rigidBodyData[i];

		if (objFiles.find(rbd.m_modelFile) == objFiles.end())
			continue;

		id_index[rbd.m_id] = i;

		PBD::VertexData &vd = objFiles[rbd.m_modelFile].first;
		PBD::IndexedFaceMesh &mesh = objFiles[rbd.m_modelFile].second;

		rb[i] = new PBD::RigidBody();

		rb[i]->initBody(rbd.m_density,
			rbd.m_x,
			rbd.m_q,
			vd, mesh,
			rbd.m_scale);

		if (!rbd.m_isDynamic)
			rb[i]->setMass(0.0);
		else
		{
			rb[i]->setVelocity(rbd.m_v);
			rb[i]->setAngularVelocity(rbd.m_omega);
		}
		rb[i]->setRestitutionCoeff(rbd.m_restitutionCoeff);
		rb[i]->setFrictionCoeff(rbd.m_frictionCoeff);

		const std::vector<PBD::Vector3r> *vertices = rb[i]->getGeometry().getVertexDataLocal().getVertices();
		const unsigned int nVert = static_cast<unsigned int>(vertices->size());

		switch (rbd.m_collisionObjectType)
		{
		case PBD::SceneLoader::RigidBodyData::No_Collision_Object: break;
		case PBD::SceneLoader::RigidBodyData::Sphere:
			m_cd.addCollisionSphere(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, rbd.m_collisionObjectScale[0], rbd.m_testMesh, rbd.m_invertSDF);
			break;
		case PBD::SceneLoader::RigidBodyData::Box:
			m_cd.addCollisionBox(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, rbd.m_collisionObjectScale, rbd.m_testMesh, rbd.m_invertSDF);
			break;
		case PBD::SceneLoader::RigidBodyData::Cylinder:
			m_cd.addCollisionCylinder(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, rbd.m_collisionObjectScale.head<2>(), rbd.m_testMesh, rbd.m_invertSDF);
			break;
		case PBD::SceneLoader::RigidBodyData::Torus:
			m_cd.addCollisionTorus(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, rbd.m_collisionObjectScale.head<2>(), rbd.m_testMesh, rbd.m_invertSDF);
			break;
		case PBD::SceneLoader::RigidBodyData::HollowSphere:
			m_cd.addCollisionHollowSphere(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, rbd.m_collisionObjectScale[0], rbd.m_thicknessSDF, rbd.m_testMesh, rbd.m_invertSDF);
			break;
		case PBD::SceneLoader::RigidBodyData::HollowBox:
			m_cd.addCollisionHollowBox(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, rbd.m_collisionObjectScale, rbd.m_thicknessSDF, rbd.m_testMesh, rbd.m_invertSDF);
			break;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// triangle models
	//////////////////////////////////////////////////////////////////////////

	// map file names to loaded geometry to prevent multiple imports of same files
	std::map<std::string, pair<PBD::VertexData, PBD::IndexedFaceMesh>> triFiles;
	for (unsigned int i = 0; i < data.m_triangleModelData.size(); i++)
	{
		const PBD::SceneLoader::TriangleModelData &tmd = data.m_triangleModelData[i];

		// Check if already loaded
		if (triFiles.find(tmd.m_modelFile) == triFiles.end())
		{
			PBD::IndexedFaceMesh mesh;
			PBD::VertexData vd;
			PBD::OBJLoader::loadObj(PBD::Utilities::normalizePath(tmd.m_modelFile), vd, mesh);
			triFiles[tmd.m_modelFile] = { vd, mesh };
		}
	}

	triModels.reserve(data.m_triangleModelData.size());
	std::map<unsigned int, unsigned int> tm_id_index;
	for (unsigned int i = 0; i < data.m_triangleModelData.size(); i++)
	{
		const PBD::SceneLoader::TriangleModelData &tmd = data.m_triangleModelData[i];

		if (triFiles.find(tmd.m_modelFile) == triFiles.end())
			continue;

		tm_id_index[tmd.m_id] = i;

		PBD::VertexData vd = triFiles[tmd.m_modelFile].first;
		PBD::IndexedFaceMesh &mesh = triFiles[tmd.m_modelFile].second;

		const PBD::Matrix3r R = tmd.m_q.matrix();
		for (unsigned int j = 0; j < vd.size(); j++)
		{
			vd.getPosition(j) = R * (vd.getPosition(j).cwiseProduct(tmd.m_scale)) + tmd.m_x;
		}

		m_model.addTriangleModel(vd.size(), mesh.numFaces(), &vd.getPosition(0), mesh.getFaces().data(), mesh.getUVIndices(), mesh.getUVs());

		PBD::TriangleModel *tm = triModels[triModels.size() - 1];
		PBD::ParticleData &pd = m_model.getParticles();
		unsigned int offset = tm->getIndexOffset();

		for (unsigned int j = 0; j < tmd.m_staticParticles.size(); j++)
		{
			const unsigned int index = tmd.m_staticParticles[j] + offset;
			pd.setMass(index, 0.0);
		}

		tm->setRestitutionCoeff(tmd.m_restitutionCoeff);
		tm->setFrictionCoeff(tmd.m_frictionCoeff);
	}

	initTriangleModelConstraints();

	//////////////////////////////////////////////////////////////////////////
	// tet models
	//////////////////////////////////////////////////////////////////////////

	// map file names to loaded geometry to prevent multiple imports of same files
	std::map<pair<string, string>, pair<vector<PBD::Vector3r>, vector<unsigned int>>> tetFiles;
	for (unsigned int i = 0; i < data.m_tetModelData.size(); i++)
	{
		const PBD::SceneLoader::TetModelData &tmd = data.m_tetModelData[i];

		// Check if already loaded
		pair<string, string> fileNames = { tmd.m_modelFileNodes, tmd.m_modelFileElements };
		if (tetFiles.find(fileNames) == tetFiles.end())
		{
			vector<PBD::Vector3r> vertices;
			vector<unsigned int> tets;
			PBD::TetGenLoader::loadTetgenModel(PBD::Utilities::normalizePath(tmd.m_modelFileNodes), PBD::Utilities::normalizePath(tmd.m_modelFileElements), vertices, tets);
			tetFiles[fileNames] = { vertices, tets };
		}
	}

	tetModels.reserve(data.m_tetModelData.size());
	std::map<unsigned int, unsigned int> tm_id_index2;
	for (unsigned int i = 0; i < data.m_tetModelData.size(); i++)
	{
		const PBD::SceneLoader::TetModelData &tmd = data.m_tetModelData[i];

		pair<string, string> fileNames = { tmd.m_modelFileNodes, tmd.m_modelFileElements };
		auto geo = tetFiles.find(fileNames);
		if (geo == tetFiles.end())
			continue;

		tm_id_index2[tmd.m_id] = i;

		vector<PBD::Vector3r> vertices = geo->second.first;
		vector<unsigned int> &tets = geo->second.second;

		const PBD::Matrix3r R = tmd.m_q.matrix();
		for (unsigned int j = 0; j < vertices.size(); j++)
		{
			vertices[j] = R * (vertices[j].cwiseProduct(tmd.m_scale)) + tmd.m_x;
		}

		m_model.addTetModel((unsigned int)vertices.size(), (unsigned int)tets.size() / 4, vertices.data(), tets.data());

		PBD::TetModel *tm = tetModels[tetModels.size() - 1];
		PBD::ParticleData &pd = m_model.getParticles();
		unsigned int offset = tm->getIndexOffset();

		for (unsigned int j = 0; j < tmd.m_staticParticles.size(); j++)
		{
			const unsigned int index = tmd.m_staticParticles[j] + offset;
			pd.setMass(index, 0.0);
		}

		// read visualization mesh
		if (tmd.m_modelFileVis != "")
		{
			if (objFiles.find(tmd.m_modelFileVis) != objFiles.end())
			{
				PBD::IndexedFaceMesh &visMesh = tm->getVisMesh();
				PBD::VertexData &vdVis = tm->getVisVertices();
				vdVis = objFiles[tmd.m_modelFileVis].first;
				visMesh = objFiles[tmd.m_modelFileVis].second;

				for (unsigned int j = 0; j < vdVis.size(); j++)
					vdVis.getPosition(j) = R * (vdVis.getPosition(j).cwiseProduct(tmd.m_scale)) + tmd.m_x;

				tm->updateMeshNormals(pd);
				tm->attachVisMesh(pd);
				tm->updateVisMesh(pd);
			}
		}

		tm->setRestitutionCoeff(tmd.m_restitutionCoeff);
		tm->setFrictionCoeff(tmd.m_frictionCoeff);

		tm->updateMeshNormals(pd);
	}

	initTetModelConstraints();

	// init collision objects for deformable models
	PBD::ParticleData &pd = m_model.getParticles();
	for (unsigned int i = 0; i < data.m_triangleModelData.size(); i++)
	{
		PBD::TriangleModel *tm = triModels[i];
		unsigned int offset = tm->getIndexOffset();
		const unsigned int nVert = tm->getParticleMesh().numVertices();
		m_cd.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert);

	}
	for (unsigned int i = 0; i < data.m_tetModelData.size(); i++)
	{
		PBD::TetModel *tm = tetModels[i];
		unsigned int offset = tm->getIndexOffset();
		const unsigned int nVert = tm->getParticleMesh().numVertices();
		m_cd.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TetModelCollisionObjectType, &pd.getPosition(offset), nVert);
	}

	//////////////////////////////////////////////////////////////////////////
	// joints
	//////////////////////////////////////////////////////////////////////////

	for (unsigned int i = 0; i < data.m_ballJointData.size(); i++)
	{
		const PBD::SceneLoader::BallJointData &jd = data.m_ballJointData[i];
		m_model.addBallJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position);
	}

	for (unsigned int i = 0; i < data.m_ballOnLineJointData.size(); i++)
	{
		const PBD::SceneLoader::BallOnLineJointData &jd = data.m_ballOnLineJointData[i];
		m_model.addBallOnLineJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
	}

	for (unsigned int i = 0; i < data.m_hingeJointData.size(); i++)
	{
		const PBD::SceneLoader::HingeJointData &jd = data.m_hingeJointData[i];
		m_model.addHingeJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
	}

	for (unsigned int i = 0; i < data.m_universalJointData.size(); i++)
	{
		const PBD::SceneLoader::UniversalJointData &jd = data.m_universalJointData[i];
		m_model.addUniversalJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis[0], jd.m_axis[1]);
	}

	for (unsigned int i = 0; i < data.m_sliderJointData.size(); i++)
	{
		const PBD::SceneLoader::SliderJointData &jd = data.m_sliderJointData[i];
		m_model.addSliderJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
	}

	for (unsigned int i = 0; i < data.m_rigidBodyParticleBallJointData.size(); i++)
	{
		const PBD::SceneLoader::RigidBodyParticleBallJointData &jd = data.m_rigidBodyParticleBallJointData[i];
		m_model.addRigidBodyParticleBallJoint(id_index[jd.m_bodyID[0]], jd.m_bodyID[1]);
	}

	for (unsigned int i = 0; i < data.m_targetAngleMotorHingeJointData.size(); i++)
	{
		const PBD::SceneLoader::TargetAngleMotorHingeJointData &jd = data.m_targetAngleMotorHingeJointData[i];
		m_model.addTargetAngleMotorHingeJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
		((PBD::TargetAngleMotorHingeJoint*)constraints[constraints.size() - 1])->setTarget(jd.m_target);
	}

	for (unsigned int i = 0; i < data.m_targetVelocityMotorHingeJointData.size(); i++)
	{
		const PBD::SceneLoader::TargetVelocityMotorHingeJointData &jd = data.m_targetVelocityMotorHingeJointData[i];
		m_model.addTargetVelocityMotorHingeJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
		((PBD::TargetVelocityMotorHingeJoint*)constraints[constraints.size() - 1])->setTarget(jd.m_target);
	}

	for (unsigned int i = 0; i < data.m_targetPositionMotorSliderJointData.size(); i++)
	{
		const PBD::SceneLoader::TargetPositionMotorSliderJointData &jd = data.m_targetPositionMotorSliderJointData[i];
		m_model.addTargetPositionMotorSliderJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
		((PBD::TargetPositionMotorSliderJoint*)constraints[constraints.size() - 1])->setTarget(jd.m_target);
	}

	for (unsigned int i = 0; i < data.m_targetVelocityMotorSliderJointData.size(); i++)
	{
		const PBD::SceneLoader::TargetVelocityMotorSliderJointData &jd = data.m_targetVelocityMotorSliderJointData[i];
		m_model.addTargetVelocityMotorSliderJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
		((PBD::TargetVelocityMotorSliderJoint*)constraints[constraints.size() - 1])->setTarget(jd.m_target);
	}

	m_cd.updateAABBs(m_model);
}

void PBDWrapper::initModel (const Real timeStepSize)
{
	PBD::TimeManager::getCurrent ()->setTimeStepSize(timeStepSize);
 
	m_sim.setCollisionDetection(m_model, &m_cd);
}

void PBDWrapper::shaderBegin(const float *col)
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

void PBDWrapper::shaderEnd()
{
	if (m_shader)
		m_shader->end();
}
 
void PBDWrapper::renderTriangleModels()
{
	const PBD::ParticleData &pd = m_model.getParticles();
	float surfaceColor[4] = { 0.8f, 0.9f, 0.2f, 1 };

	shaderBegin(surfaceColor);

	for (unsigned int i = 0; i < m_model.getTriangleModels().size(); i++)
	{
		// mesh 
		const PBD::IndexedFaceMesh &mesh = m_model.getTriangleModels()[i]->getParticleMesh();
		const unsigned int offset = m_model.getTriangleModels()[i]->getIndexOffset();
		drawMesh(pd, mesh, offset, surfaceColor);
	}

	shaderEnd();
}
 
void PBDWrapper::renderTetModels()
{
	const PBD::ParticleData &pd = m_model.getParticles();
	float surfaceColor[4] = { 0.8f, 0.4f, 0.7f, 1 };
 
	shaderBegin(surfaceColor);
 
	for (unsigned int i = 0; i < m_model.getTetModels().size(); i++)
	{
		const PBD::VertexData &vdVis = m_model.getTetModels()[i]->getVisVertices();
		if (vdVis.size() > 0)
		{
			const PBD::IndexedFaceMesh &visMesh = m_model.getTetModels()[i]->getVisMesh();
			drawMesh(vdVis, visMesh, 0, surfaceColor);
		}
		else
		{
			const PBD::IndexedFaceMesh &surfaceMesh = m_model.getTetModels()[i]->getSurfaceMesh();
			const unsigned int offset = m_model.getTetModels()[i]->getIndexOffset();
			drawMesh(pd, surfaceMesh, offset, surfaceColor);
		}
	}
 
	shaderEnd();
}
 
void PBDWrapper::renderAABB(PBD::AABB &aabb)
{
	PBD::Vector3r p1, p2;
	glBegin(GL_LINES);
	for (unsigned char i = 0; i < 12; i++)
	{
		PBD::AABB::getEdge(aabb, i, p1, p2);
		glVertex3d(p1[0], p1[1], p1[2]);
		glVertex3d(p2[0], p2[1], p2[2]);
	}
	glEnd();
}
 
void PBDWrapper::renderBallJoint(PBD::BallJoint &bj)
{	
	SPH::MiniGL::drawSphere(bj.m_jointInfo.col(2), 0.15f, jointColor);
}
 
void PBDWrapper::renderRigidBodyParticleBallJoint(PBD::RigidBodyParticleBallJoint &bj)
{
	SPH::MiniGL::drawSphere(bj.m_jointInfo.col(1), 0.1f, jointColor);
}
 
void PBDWrapper::renderBallOnLineJoint(PBD::BallOnLineJoint &bj)
{
	SPH::MiniGL::drawSphere(bj.m_jointInfo.col(5), 0.1f, jointColor);
	SPH::MiniGL::drawCylinder(bj.m_jointInfo.col(5) - bj.m_jointInfo.col(7), bj.m_jointInfo.col(5) + bj.m_jointInfo.col(7), jointColor, 0.05f);
}
 
void PBDWrapper::renderHingeJoint(PBD::HingeJoint &hj)
{
	SPH::MiniGL::drawSphere(hj.m_jointInfo.col(6) - 0.5*hj.m_jointInfo.col(8), 0.1f, jointColor);
	SPH::MiniGL::drawSphere(hj.m_jointInfo.col(6) + 0.5*hj.m_jointInfo.col(8), 0.1f, jointColor);
	SPH::MiniGL::drawCylinder(hj.m_jointInfo.col(6) - 0.5*hj.m_jointInfo.col(8), hj.m_jointInfo.col(6) + 0.5*hj.m_jointInfo.col(8), jointColor, 0.05f);
}

void PBDWrapper::renderUniversalJoint(PBD::UniversalJoint &uj)
{
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(4) - 0.5*uj.m_jointInfo.col(6), 0.1f, jointColor);
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(4) + 0.5*uj.m_jointInfo.col(6), 0.1f, jointColor);
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(5) - 0.5*uj.m_jointInfo.col(7), 0.1f, jointColor);
	SPH::MiniGL::drawSphere(uj.m_jointInfo.col(5) + 0.5*uj.m_jointInfo.col(7), 0.1f, jointColor);
	SPH::MiniGL::drawCylinder(uj.m_jointInfo.col(4) - 0.5*uj.m_jointInfo.col(6), uj.m_jointInfo.col(4) + 0.5*uj.m_jointInfo.col(6), jointColor, 0.05f);
	SPH::MiniGL::drawCylinder(uj.m_jointInfo.col(5) - 0.5*uj.m_jointInfo.col(7), uj.m_jointInfo.col(5) + 0.5*uj.m_jointInfo.col(7), jointColor, 0.05f);
}
 
void PBDWrapper::renderSliderJoint(PBD::SliderJoint &joint)
{
	SPH::MiniGL::drawSphere(joint.m_jointInfo.col(6), 0.1f, jointColor);
	SPH::MiniGL::drawCylinder(joint.m_jointInfo.col(7) - joint.m_jointInfo.col(8), joint.m_jointInfo.col(7) + joint.m_jointInfo.col(8), jointColor, 0.05f);
}
 
void PBDWrapper::renderTargetPositionMotorSliderJoint(PBD::TargetPositionMotorSliderJoint &joint)
{
	SPH::MiniGL::drawSphere(joint.m_jointInfo.col(6), 0.1f, jointColor);
	SPH::MiniGL::drawCylinder(joint.m_jointInfo.col(7) - joint.m_jointInfo.col(8), joint.m_jointInfo.col(7) + joint.m_jointInfo.col(8), jointColor, 0.05f);
}
 
void PBDWrapper::renderTargetVelocityMotorSliderJoint(PBD::TargetVelocityMotorSliderJoint &joint)
{
	SPH::MiniGL::drawSphere(joint.m_jointInfo.col(6), 0.1f, jointColor);
	SPH::MiniGL::drawCylinder(joint.m_jointInfo.col(7) - joint.m_jointInfo.col(8), joint.m_jointInfo.col(7) + joint.m_jointInfo.col(8), jointColor, 0.05f);
}
 
void PBDWrapper::renderTargetAngleMotorHingeJoint(PBD::TargetAngleMotorHingeJoint &hj)
{
	SPH::MiniGL::drawSphere(hj.m_jointInfo.col(6) - 0.5*hj.m_jointInfo.col(8), 0.1f, jointColor);
	SPH::MiniGL::drawSphere(hj.m_jointInfo.col(6) + 0.5*hj.m_jointInfo.col(8), 0.1f, jointColor);
	SPH::MiniGL::drawCylinder(hj.m_jointInfo.col(6) - 0.5*hj.m_jointInfo.col(8), hj.m_jointInfo.col(6) + 0.5*hj.m_jointInfo.col(8), jointColor, 0.05f);
}
 
void PBDWrapper::renderTargetVelocityMotorHingeJoint(PBD::TargetVelocityMotorHingeJoint &hj)
{
	SPH::MiniGL::drawSphere(hj.m_jointInfo.col(6) - 0.5*hj.m_jointInfo.col(8), 0.1f, jointColor);
	SPH::MiniGL::drawSphere(hj.m_jointInfo.col(6) + 0.5*hj.m_jointInfo.col(8), 0.1f, jointColor);
	SPH::MiniGL::drawCylinder(hj.m_jointInfo.col(6) - 0.5*hj.m_jointInfo.col(8), hj.m_jointInfo.col(6) + 0.5*hj.m_jointInfo.col(8), jointColor, 0.05f);
}
 
void PBDWrapper::renderRigidBodyContact(PBD::RigidBodyContactConstraint &cc)
{
	float col1[4] = { 0.0f, 0.6f, 0.2f, 1 };
	float col2[4] = { 0.6f, 0.0f, 0.2f, 1 };
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(0), 5.0f, col1);
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(1), 5.0f, col2);
	SPH::MiniGL::drawVector(cc.m_constraintInfo.col(1), cc.m_constraintInfo.col(1) + cc.m_constraintInfo.col(2), 1.0f, col2);
}
 
void PBDWrapper::renderParticleRigidBodyContact(PBD::ParticleRigidBodyContactConstraint &cc)
{
	float col1[4] = { 0.0f, 0.6f, 0.2f, 1 };
	float col2[4] = { 0.6f, 0.0f, 0.2f, 1 };
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(0), 5.0f, col1);
	SPH::MiniGL::drawPoint(cc.m_constraintInfo.col(1), 5.0f, col2);
	SPH::MiniGL::drawVector(cc.m_constraintInfo.col(1), cc.m_constraintInfo.col(1) + cc.m_constraintInfo.col(2), 1.0f, col2);
}

void PBDWrapper::renderScene ()
{
	// Draw sim model	
	PBD::SimulationModel::RigidBodyVector &rb = m_model.getRigidBodies();

	//float selectionColor[4] = { 0.8f, 0.0f, 0.0f, 1 };
	float surfaceColor[4] = { 0.3f, 0.5f, 0.8f, 1 };
	float staticColor[4] = { 0.5f, 0.5f, 0.5f, 0.3f };

	for (size_t i = 0; i < rb.size(); i++)
	{
		const PBD::VertexData &vd = rb[i]->getGeometry().getVertexData();
		const PBD::IndexedFaceMesh &mesh = rb[i]->getGeometry().getMesh();
		if (rb[i]->getMass() != 0.0)
		{
			shaderBegin(surfaceColor);
			drawMesh(vd, mesh, 0, surfaceColor);
			shaderEnd();
		}
		else if (m_drawStaticBodies)
		{
			shaderBegin(staticColor);
			drawMesh(vd, mesh, 0, staticColor);
			shaderEnd();
		}
	}

	renderTriangleModels();
	renderTetModels();
	renderConstraints(); 
	renderBVH();
}

void PBDWrapper::renderConstraints()
{
	PBD::SimulationModel::ConstraintVector &constraints = m_model.getConstraints();
	PBD::SimulationModel::RigidBodyContactConstraintVector &rigidBodyContacts = m_model.getRigidBodyContactConstraints();
	PBD::SimulationModel::ParticleRigidBodyContactConstraintVector &particleRigidBodyContacts = m_model.getParticleRigidBodyContactConstraints();

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
	}
}

void PBDWrapper::renderBVH()
{
	if (m_drawAABB || (m_drawBVHDepth >= 0))
	{
		float staticColor[4] = { 0.5f, 0.5f, 0.5f, 0.3f };

		PBD::SimulationModel::RigidBodyVector &rb = m_model.getRigidBodies();
		PBD::ObjectArray<PBD::CollisionDetection::CollisionObject*> &collisionObjects = m_cd.getCollisionObjects();
		for (unsigned int k = 0; k < collisionObjects.size(); k++)
		{
			if (m_drawAABB)
				renderAABB(collisionObjects[k]->m_aabb);

			if (m_drawBVHDepth >= 0)
			{
				if (m_cd.isDistanceFieldCollisionObject(collisionObjects[k]))
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
								const PBD::Vector3r &sphere_x = bs.x();
								const PBD::Vector3r sphere_x_w = body->getRotation() * sphere_x + body->getPosition();
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

 
 void PBDWrapper::initTriangleModelConstraints()
 {
	// init constraints
	for (unsigned int cm = 0; cm < m_model.getTriangleModels().size(); cm++)
	{
		const unsigned int offset = m_model.getTriangleModels()[cm]->getIndexOffset();
		if (m_clothSimulationMethod == 1)
		{
			const unsigned int nEdges = m_model.getTriangleModels()[cm]->getParticleMesh().numEdges();
			const PBD::IndexedFaceMesh::Edge *edges = m_model.getTriangleModels()[cm]->getParticleMesh().getEdges().data();
			for (unsigned int i = 0; i < nEdges; i++)
			{
				const unsigned int v1 = edges[i].m_vert[0] + offset;
				const unsigned int v2 = edges[i].m_vert[1] + offset;
 
				m_model.addDistanceConstraint(v1, v2);
			}
		}
		else if (m_clothSimulationMethod == 2)
		{
			
			PBD::TriangleModel::ParticleMesh &mesh = m_model.getTriangleModels()[cm]->getParticleMesh();
			const unsigned int *tris = mesh.getFaces().data();
			const unsigned int nFaces = mesh.numFaces();
			for (unsigned int i = 0; i < nFaces; i++)
			{
				const unsigned int v1 = tris[3 * i] + offset;
				const unsigned int v2 = tris[3 * i + 1] + offset;
				const unsigned int v3 = tris[3 * i + 2] + offset;
				m_model.addFEMTriangleConstraint(v1, v2, v3);
			}
		}
		else if (m_clothSimulationMethod == 3)
		{
			PBD::TriangleModel::ParticleMesh &mesh = m_model.getTriangleModels()[cm]->getParticleMesh();
			const unsigned int *tris = mesh.getFaces().data();
			const unsigned int nFaces = mesh.numFaces();
			for (unsigned int i = 0; i < nFaces; i++)
			{
				const unsigned int v1 = tris[3 * i] + offset;
				const unsigned int v2 = tris[3 * i + 1] + offset;
				const unsigned int v3 = tris[3 * i + 2] + offset;
				m_model.addStrainTriangleConstraint(v1, v2, v3);
			}
		}
		if (m_bendingMethod != 0)
		{
			PBD::TriangleModel::ParticleMesh &mesh = m_model.getTriangleModels()[cm]->getParticleMesh();
			unsigned int nEdges = mesh.numEdges();
			const PBD::TriangleModel::ParticleMesh::Edge *edges = mesh.getEdges().data();
			const unsigned int *tris = mesh.getFaces().data();
			for (unsigned int i = 0; i < nEdges; i++)
			{
				const int tri1 = edges[i].m_face[0];
				const int tri2 = edges[i].m_face[1];
				if ((tri1 != 0xffffffff) && (tri2 != 0xffffffff))
				{
					// Find the triangle points which do not lie on the axis
					const int axisPoint1 = edges[i].m_vert[0];
					const int axisPoint2 = edges[i].m_vert[1];
					int point1 = -1;
					int point2 = -1;
					for (int j = 0; j < 3; j++)
					{
						if ((tris[3 * tri1 + j] != axisPoint1) && (tris[3 * tri1 + j] != axisPoint2))
						{
							point1 = tris[3 * tri1 + j];
							break;
						}
					}
					for (int j = 0; j < 3; j++)
					{
						if ((tris[3 * tri2 + j] != axisPoint1) && (tris[3 * tri2 + j] != axisPoint2))
						{
							point2 = tris[3 * tri2 + j];
							break;
						}
					}
					if ((point1 != -1) && (point2 != -1))
					{
						const unsigned int vertex1 = point1 + offset;
						const unsigned int vertex2 = point2 + offset;
						const unsigned int vertex3 = edges[i].m_vert[0] + offset;
						const unsigned int vertex4 = edges[i].m_vert[1] + offset;
						if (m_bendingMethod == 1)
							m_model.addDihedralConstraint(vertex1, vertex2, vertex3, vertex4);
						else if (m_bendingMethod == 2)
							m_model.addIsometricBendingConstraint(vertex1, vertex2, vertex3, vertex4);
					}
				}
			}
		}
	}
 }
 
 void PBDWrapper::initTetModelConstraints()
 {
	// init constraints
	for (unsigned int cm = 0; cm < m_model.getTetModels().size(); cm++)
	{
		const unsigned int offset = m_model.getTetModels()[cm]->getIndexOffset();
		const unsigned int nTets = m_model.getTetModels()[cm]->getParticleMesh().numTets();
		const unsigned int *tets = m_model.getTetModels()[cm]->getParticleMesh().getTets().data();
		const PBD::IndexedTetMesh::VertexTets *vTets = m_model.getTetModels()[cm]->getParticleMesh().getVertexTets().data();
		if (m_solidSimulationMethod == 1)
		{
			const unsigned int offset = m_model.getTetModels()[cm]->getIndexOffset();
			const unsigned int nEdges = m_model.getTetModels()[cm]->getParticleMesh().numEdges();
			const PBD::IndexedTetMesh::Edge *edges = m_model.getTetModels()[cm]->getParticleMesh().getEdges().data();
			for (unsigned int i = 0; i < nEdges; i++)
			{
				const unsigned int v1 = edges[i].m_vert[0] + offset;
				const unsigned int v2 = edges[i].m_vert[1] + offset;
 
				m_model.addDistanceConstraint(v1, v2);
			}
 
			for (unsigned int i = 0; i < nTets; i++)
			{
				const unsigned int v1 = tets[4 * i] + offset;
				const unsigned int v2 = tets[4 * i + 1] + offset;
				const unsigned int v3 = tets[4 * i + 2] + offset;
				const unsigned int v4 = tets[4 * i + 3] + offset;
 
				m_model.addVolumeConstraint(v1, v2, v3, v4);
			}
		}
		else if (m_solidSimulationMethod == 2)
		{
			PBD::TetModel::ParticleMesh &mesh = m_model.getTetModels()[cm]->getParticleMesh();
			for (unsigned int i = 0; i < nTets; i++)
			{
				const unsigned int v1 = tets[4 * i] + offset;
				const unsigned int v2 = tets[4 * i + 1] + offset;
				const unsigned int v3 = tets[4 * i + 2] + offset;
				const unsigned int v4 = tets[4 * i + 3] + offset;
 
				m_model.addFEMTetConstraint(v1, v2, v3, v4);
			}
		}
		else if (m_solidSimulationMethod == 3)
		{
			PBD::TetModel::ParticleMesh &mesh = m_model.getTetModels()[cm]->getParticleMesh();
			for (unsigned int i = 0; i < nTets; i++)
			{
				const unsigned int v1 = tets[4 * i] + offset;
				const unsigned int v2 = tets[4 * i + 1] + offset;
				const unsigned int v3 = tets[4 * i + 2] + offset;
				const unsigned int v4 = tets[4 * i + 3] + offset;
 
				m_model.addStrainTetConstraint(v1, v2, v3, v4);
			}
		}
		else if (m_solidSimulationMethod == 4)
		{
			PBD::TetModel::ParticleMesh &mesh = m_model.getTetModels()[cm]->getParticleMesh();
			for (unsigned int i = 0; i < nTets; i++)
			{
				const unsigned int v[4] = { tets[4 * i] + offset,
											tets[4 * i + 1] + offset, 
											tets[4 * i + 2] + offset, 
											tets[4 * i + 3] + offset };
				// Important: Divide position correction by the number of clusters 
				// which contain the vertex.
				const unsigned int nc[4] = { vTets[v[0]].m_numTets, vTets[v[1]].m_numTets, vTets[v[2]].m_numTets, vTets[v[3]].m_numTets };
				m_model.addShapeMatchingConstraint(4, v, nc);
			}
		}
	}
 }
 
 
void TW_CALL PBDWrapper::setVelocityUpdateMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	((PBD::TimeStepController*)clientData)->setVelocityUpdateMethod((unsigned int)val);
}

void TW_CALL PBDWrapper::getVelocityUpdateMethod(void *value, void *clientData)
{
	*(short *)(value) = (short)((PBD::TimeStepController*)clientData)->getVelocityUpdateMethod();
}

void TW_CALL PBDWrapper::setMaxIterations(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	((PBD::TimeStepController*)clientData)->setMaxIterations(val);
}

void TW_CALL PBDWrapper::getMaxIterations(void *value, void *clientData)
{
	*(unsigned int *)(value) = ((PBD::TimeStepController*)clientData)->getMaxIterations();
}

void TW_CALL PBDWrapper::setMaxIterationsV(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	((PBD::TimeStepController*)clientData)->setMaxIterationsV(val);
}

void TW_CALL PBDWrapper::getMaxIterationsV(void *value, void *clientData)
{
	*(unsigned int *)(value) = ((PBD::TimeStepController*)clientData)->getMaxIterationsV();
}

void TW_CALL PBDWrapper::setContactTolerance(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::DistanceFieldCollisionDetection*)clientData)->setTolerance(val);
}

void TW_CALL PBDWrapper::getContactTolerance(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::DistanceFieldCollisionDetection*)clientData)->getTolerance();
}

void TW_CALL PBDWrapper::setContactStiffnessRigidBody(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setContactStiffnessRigidBody(val);
}

void TW_CALL PBDWrapper::getContactStiffnessRigidBody(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getContactStiffnessRigidBody();
}

void TW_CALL PBDWrapper::setContactStiffnessParticleRigidBody(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setContactStiffnessParticleRigidBody(val);
}

void TW_CALL PBDWrapper::getContactStiffnessParticleRigidBody(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getContactStiffnessParticleRigidBody();
}

void TW_CALL PBDWrapper::setStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setClothStiffness(val);
}

void TW_CALL PBDWrapper::getStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getClothStiffness();
}

void TW_CALL PBDWrapper::setXXStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setClothXXStiffness(val);
}

void TW_CALL PBDWrapper::getXXStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getClothXXStiffness();
}

void TW_CALL PBDWrapper::setYYStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setClothYYStiffness(val);
}

void TW_CALL PBDWrapper::getYYStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getClothYYStiffness();
}

void TW_CALL PBDWrapper::setXYStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setClothXYStiffness(val);
}

void TW_CALL PBDWrapper::getXYStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getClothXYStiffness();
}

void TW_CALL PBDWrapper::setYXPoissonRatio(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setClothYXPoissonRatio(val);
}

void TW_CALL PBDWrapper::getYXPoissonRatio(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getClothYXPoissonRatio();
}

void TW_CALL PBDWrapper::setXYPoissonRatio(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setClothXYPoissonRatio(val);
}

void TW_CALL PBDWrapper::getXYPoissonRatio(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getClothXYPoissonRatio();
}

void TW_CALL PBDWrapper::setNormalizeStretch(const void *value, void *clientData)
{
	const bool val = *(const bool *)(value);
	((PBD::SimulationModel*)clientData)->setClothNormalizeStretch(val);
}

void TW_CALL PBDWrapper::getNormalizeStretch(void *value, void *clientData)
{
	*(bool *)(value) = ((PBD::SimulationModel*)clientData)->getClothNormalizeStretch();
}

void TW_CALL PBDWrapper::setNormalizeShear(const void *value, void *clientData)
{
	const bool val = *(const bool *)(value);
	((PBD::SimulationModel*)clientData)->setClothNormalizeShear(val);
}

void TW_CALL PBDWrapper::getNormalizeShear(void *value, void *clientData)
{
	*(bool *)(value) = ((PBD::SimulationModel*)clientData)->getClothNormalizeShear();
}

void TW_CALL PBDWrapper::setBendingStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((PBD::SimulationModel*)clientData)->setClothBendingStiffness(val);
}

void TW_CALL PBDWrapper::getBendingStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((PBD::SimulationModel*)clientData)->getClothBendingStiffness();
}

void TW_CALL PBDWrapper::setBendingMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	((PBDWrapper*)clientData)->m_bendingMethod = val;
	((PBDWrapper*)clientData)->reset();
}

void TW_CALL PBDWrapper::getBendingMethod(void *value, void *clientData)
{
	*(short *)(value) = ((PBDWrapper*)clientData)->m_bendingMethod;
}

void TW_CALL PBDWrapper::setClothSimulationMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	((PBDWrapper*)clientData)->m_clothSimulationMethod = val;	
	((PBDWrapper*)clientData)->reset();
}

void TW_CALL PBDWrapper::getClothSimulationMethod(void *value, void *clientData)
{
	*(short *)(value) = ((PBDWrapper*)clientData)->m_clothSimulationMethod;
}

void TW_CALL PBDWrapper::setSolidSimulationMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	((PBDWrapper*)clientData)->m_solidSimulationMethod = val;	
	((PBDWrapper*)clientData)->reset();
}

void TW_CALL PBDWrapper::getSolidSimulationMethod(void *value, void *clientData)
{
	*(short *)(value) = ((PBDWrapper*)clientData)->m_solidSimulationMethod;
}
