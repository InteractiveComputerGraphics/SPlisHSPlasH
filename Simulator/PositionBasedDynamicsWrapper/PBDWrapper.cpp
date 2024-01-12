#include "PBDWrapper.h"
#include "Common/Common.h"
#include <iostream>
#include "Utilities/OBJLoader.h"
#include "Utils/SceneLoader.h"
#include "Utils/TetGenLoader.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"
#include "Utilities/FileSystem.h"
#include "Simulation/TimeManager.h"
#include "Simulation/Simulation.h"
#include "Simulation/TimeStepController.h"

#define _USE_MATH_DEFINES
#include "math.h"


using namespace Eigen;
using namespace std;


PBDWrapper::PBDWrapper()
{
	m_clothSimulationMethod = 2;
	m_solidSimulationMethod = 2;
	m_bendingMethod = 2;
	m_distanceStiffness = 1.0;
	m_xxStiffness = 1.0;
	m_yyStiffness = 1.0;
	m_xyStiffness = 1.0;
	m_xyPoissonRatio = static_cast<Real>(0.3);
	m_yxPoissonRatio = static_cast<Real>(0.3);
	m_normalizeStretch = false;
	m_normalizeShear = false;
	m_bendingStiffness = static_cast<Real>(0.01);
	m_solidStiffness = 1.0;
	m_solidPoissonRatio = static_cast<Real>(0.3);
	m_solidNormalizeStretch = false;
	m_solidNormalizeShear = false;
	m_volumeStiffness = 1.0;
	m_sceneName = "";
	m_sceneFileName = "";
	m_dampingCoeff = 0.0;
	m_model.init();
	PBD::Simulation::getCurrent()->setModel(&m_model);
}

PBDWrapper::~PBDWrapper()
{
	delete PBD::Simulation::getCurrent();
}

 void PBDWrapper::reset()
 {
	 PBD::Simulation::getCurrent()->reset();
 }
 

 
 void PBDWrapper::timeStep()
 {
	PBD::ParticleData &pd = m_model.getParticles();
	PBD::SimulationModel::RigidBodyVector &rb = m_model.getRigidBodies();
	PBD::TimeManager::getCurrent()->setTimeStepSize(SPH::TimeManager::getCurrent()->getTimeStepSize());
	PBD::TimeManager::getCurrent()->setTime(SPH::TimeManager::getCurrent()->getTime());

	PBD::Simulation::getCurrent()->getTimeStep()->step(m_model);

	for (unsigned int i = 0; i < pd.size(); i++)
	{
		pd.getVelocity(i) *= (static_cast<Real>(1.0) - m_dampingCoeff);
	}
	for (unsigned int i = 0; i < rb.size(); i++)
	{
		rb[i]->getVelocity() *= (static_cast<Real>(1.0) - m_dampingCoeff);
		rb[i]->getAngularVelocity() *= (static_cast<Real>(1.0) - m_dampingCoeff);
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

void PBDWrapper::loadObj(const std::string &filename, PBD::VertexData &vd, Utilities::IndexedFaceMesh &mesh, const Vector3r &scale)
{
	std::vector<Utilities::OBJLoader::Vec3f> x;
	std::vector<Utilities::OBJLoader::Vec3f> normals;
	std::vector<Utilities::OBJLoader::Vec2f> texCoords;
	std::vector<Utilities::MeshFaceIndices> faces;
	Utilities::OBJLoader::Vec3f s = { (float)scale[0], (float)scale[1], (float)scale[2] };
	Utilities::OBJLoader::loadObj(filename, &x, &faces, &normals, &texCoords, s);

	mesh.release();
	const unsigned int nPoints = (unsigned int)x.size();
	const unsigned int nFaces = (unsigned int)faces.size();
	const unsigned int nTexCoords = (unsigned int)texCoords.size();
	mesh.initMesh(nPoints, nFaces * 2, nFaces);
	vd.reserve(nPoints);
	for (unsigned int i = 0; i < nPoints; i++)
	{
		vd.addVertex(Vector3r(x[i][0], x[i][1], x[i][2]));
	}
	for (unsigned int i = 0; i < nTexCoords; i++)
	{
		mesh.addUV(texCoords[i][0], texCoords[i][1]);
	}
	for (unsigned int i = 0; i < nFaces; i++)
	{
		int posIndices[3];
		int texIndices[3];
		for (int j = 0; j < 3; j++)
		{
			posIndices[j] = faces[i].posIndices[j];
			if (nTexCoords > 0)
			{
				texIndices[j] = faces[i].texIndices[j];
				mesh.addUVIndex(texIndices[j]);
			}
		}

		mesh.addFace(&posIndices[0]);
	}
	mesh.buildNeighbors();

	mesh.updateNormals(vd, 0);
	mesh.updateVertexNormals(vd);

	LOG_INFO << "Number of triangles: " << nFaces;
	LOG_INFO << "Number of vertices: " << nPoints;
}

void PBDWrapper::readScene(const std::string &sceneFileName, const std::vector< RBData> &additionalRigidBodies)
{
	m_sceneFileName = sceneFileName;

	PBD::SimulationModel::RigidBodyVector &rb = m_model.getRigidBodies();
	PBD::SimulationModel::TriangleModelVector &triModels = m_model.getTriangleModels();
	PBD::SimulationModel::TetModelVector &tetModels = m_model.getTetModels();
	PBD::SimulationModel::ConstraintVector &constraints = m_model.getConstraints();

	Utilities::SceneLoader::SceneData data;
	Utilities::SceneLoader loader;
	loader.readScene(sceneFileName, data);
	LOG_INFO << "Scene: " << sceneFileName;

	std::string basePath = Utilities::FileSystem::getFilePath(sceneFileName);
	
	m_sceneName = data.m_sceneName;

	PBD::Simulation *sim = PBD::Simulation::getCurrent();
	loader.readParameterObject(sim);
	loader.readParameterObject(sim->getModel());
	loader.readParameterObject(sim->getTimeStep());
	if (sim->getTimeStep()->getCollisionDetection() != nullptr)
		loader.readParameterObject(sim->getTimeStep()->getCollisionDetection());


	sim->setVecValue<Real>(PBD::Simulation::GRAVITATION, &data.m_gravity[0]);

	//////////////////////////////////////////////////////////////////////////
	// rigid bodies
	//////////////////////////////////////////////////////////////////////////

	// add additional bodies
	for (auto i=0; i < additionalRigidBodies.size(); i++)
	{
		Utilities::SceneLoader::RigidBodyData rbd;
		rbd.m_modelFile = additionalRigidBodies[i].objFile;
		rbd.m_collisionObjectFileName = "";
		rbd.m_x = additionalRigidBodies[i].x;
		rbd.m_q = Quaternionr(additionalRigidBodies[i].R);
		rbd.m_scale = additionalRigidBodies[i].scale;
		rbd.m_collisionObjectScale = additionalRigidBodies[i].scale;
		rbd.m_collisionObjectType = additionalRigidBodies[i].collisionType;
		rbd.m_frictionCoeff = additionalRigidBodies[i].friction;
		rbd.m_restitutionCoeff = additionalRigidBodies[i].restitution;
		rbd.m_density = 1000.0;
		rbd.m_invertSDF = false;
		rbd.m_isDynamic = false;
		rbd.m_omega.setZero();
		rbd.m_v.setZero();
		rbd.m_resolutionSDF = Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign>(10, 10, 10);
		rbd.m_thicknessSDF = 0.0;
		rbd.m_testMesh = false;
		data.m_rigidBodyData.push_back(rbd);
	}

	// map file names to loaded geometry to prevent multiple imports of same files
	std::map<std::string, pair<PBD::VertexData, Utilities::IndexedFaceMesh>> objFiles;
	std::map<std::string, PBD::CubicSDFCollisionDetection::GridPtr> distanceFields;
	for (unsigned int i = 0; i < data.m_rigidBodyData.size(); i++)
	{
		Utilities::SceneLoader::RigidBodyData &rbd = data.m_rigidBodyData[i];

		// Check if already loaded
		rbd.m_modelFile = Utilities::FileSystem::normalizePath(rbd.m_modelFile);
		if (objFiles.find(rbd.m_modelFile) == objFiles.end())
		{
			Utilities::IndexedFaceMesh mesh;
			PBD::VertexData vd;
			loadObj(rbd.m_modelFile, vd, mesh, Vector3r::Ones());
			objFiles[rbd.m_modelFile] = { vd, mesh };
		}

		const std::string basePath = Utilities::FileSystem::getFilePath(sceneFileName);
		const string cachePath = basePath + "/Cache";
		const string resStr = to_string(rbd.m_resolutionSDF[0]) + "_" + to_string(rbd.m_resolutionSDF[1]) + "_" + to_string(rbd.m_resolutionSDF[2]);
		const std::string modelFileName = Utilities::FileSystem::getFileNameWithExt(rbd.m_modelFile);
		const string sdfFileName = Utilities::FileSystem::normalizePath(cachePath + "/" + modelFileName + "_" + resStr + ".csdf");

		std::string sdfKey = rbd.m_collisionObjectFileName;
		if (sdfKey == "")
		{
			sdfKey = sdfFileName;
		}
		if (distanceFields.find(sdfKey) == distanceFields.end())
		{
			// Generate SDF
			if (rbd.m_collisionObjectType == Utilities::SceneLoader::SDF)
			{
				if (rbd.m_collisionObjectFileName == "")
				{
					std::string md5FileName = Utilities::FileSystem::normalizePath(cachePath + "/" + modelFileName + ".md5");
					string md5Str = Utilities::FileSystem::getFileMD5(rbd.m_modelFile);
					bool md5 = false;
					if (Utilities::FileSystem::fileExists(md5FileName))
						md5 = Utilities::FileSystem::checkMD5(md5Str, md5FileName);

					// check MD5 if cache file is available
					const string resStr = to_string(rbd.m_resolutionSDF[0]) + "_" + to_string(rbd.m_resolutionSDF[1]) + "_" + to_string(rbd.m_resolutionSDF[2]);
					const string sdfFileName = Utilities::FileSystem::normalizePath(cachePath + "/" + modelFileName + "_" + resStr + ".csdf");
					bool foundCacheFile = Utilities::FileSystem::fileExists(sdfFileName);

					if (foundCacheFile && md5)
					{
						LOG_INFO << "Load cached SDF: " << sdfFileName;
						distanceFields[sdfFileName] = std::make_shared<PBD::CubicSDFCollisionDetection::Grid>(sdfFileName);
					}
					else
					{
						PBD::VertexData &vd = objFiles[rbd.m_modelFile].first;
						Utilities::IndexedFaceMesh &mesh = objFiles[rbd.m_modelFile].second;

						std::vector<unsigned int> &faces = mesh.getFaces();
						const unsigned int nFaces = mesh.numFaces();

#ifdef USE_DOUBLE
						Discregrid::TriangleMesh sdfMesh(&vd.getPosition(0)[0], faces.data(), vd.size(), nFaces);
#else
						// if type is float, copy vector to double vector
						std::vector<double> doubleVec;
						doubleVec.resize(3 * vd.size());
						for (unsigned int i = 0; i < vd.size(); i++)
							for (unsigned int j = 0; j < 3; j++)
								doubleVec[3 * i + j] = vd.getPosition(i)[j];
						Discregrid::TriangleMesh sdfMesh(&doubleVec[0], faces.data(), vd.size(), nFaces);
#endif
						Discregrid::MeshDistance md(sdfMesh);
						Eigen::AlignedBox3d domain;
						for (auto const& x : sdfMesh.vertices())
						{
							domain.extend(x);
						}
						domain.max() += 0.1 * Eigen::Vector3d::Ones();
						domain.min() -= 0.1 * Eigen::Vector3d::Ones();

						LOG_INFO << "Set SDF resolution: " << rbd.m_resolutionSDF[0] << ", " << rbd.m_resolutionSDF[1] << ", " << rbd.m_resolutionSDF[2];
						distanceFields[sdfFileName] = std::make_shared<PBD::CubicSDFCollisionDetection::Grid>(domain, std::array<unsigned int, 3>({ rbd.m_resolutionSDF[0], rbd.m_resolutionSDF[1], rbd.m_resolutionSDF[2] }));
						auto func = Discregrid::DiscreteGrid::ContinuousFunction{};
						func = [&md](Eigen::Vector3d const& xi) {return md.signedDistanceCached(xi); };
						LOG_INFO << "Generate SDF for " << rbd.m_modelFile;
						distanceFields[sdfFileName]->addFunction(func, true);
						if (Utilities::FileSystem::makeDir(cachePath) == 0)
						{
							LOG_INFO << "Save SDF: " << sdfFileName;
							distanceFields[sdfFileName]->save(sdfFileName);
							Utilities::FileSystem::writeMD5File(rbd.m_modelFile, md5FileName);
						}
					}
				}
				else
				{
					std::string fileName = rbd.m_collisionObjectFileName;
					if (Utilities::FileSystem::isRelativePath(fileName))
					{
						fileName = Utilities::FileSystem::normalizePath(basePath + "/" + fileName);
					}
					LOG_INFO << "Load SDF: " << fileName;
					distanceFields[rbd.m_collisionObjectFileName] = std::make_shared<PBD::CubicSDFCollisionDetection::Grid>(fileName);
				}
			}
		}
	}

	for (unsigned int i = 0; i < data.m_tetModelData.size(); i++)
	{
		const Utilities::SceneLoader::TetModelData &tmd = data.m_tetModelData[i];

		// Check if already loaded
		if ((tmd.m_modelFileVis != "") &&
			(objFiles.find(tmd.m_modelFileVis) == objFiles.end()))
		{
			Utilities::IndexedFaceMesh mesh;
			PBD::VertexData vd;
			loadObj(Utilities::FileSystem::normalizePath(tmd.m_modelFileVis), vd, mesh, Vector3r::Ones());
			objFiles[tmd.m_modelFileVis] = { vd, mesh };
		}
	}


	rb.resize(data.m_rigidBodyData.size());
	std::map<unsigned int, unsigned int> id_index;
	for (unsigned int i = 0; i < data.m_rigidBodyData.size(); i++)
	{
		const Utilities::SceneLoader::RigidBodyData &rbd = data.m_rigidBodyData[i];

		if (objFiles.find(rbd.m_modelFile) == objFiles.end())
			continue;

		id_index[rbd.m_id] = i;

		PBD::VertexData &vd = objFiles[rbd.m_modelFile].first;
		Utilities::IndexedFaceMesh &mesh = objFiles[rbd.m_modelFile].second;

		rb[i] = new PBD::RigidBody();

		rb[i]->initBody(rbd.m_density,
			rbd.m_x,
			rbd.m_q,
			vd, mesh,
			rbd.m_scale);

		if (!rbd.m_isDynamic)
		{
			rb[i]->setMass(0.0);
			rb[i]->setVelocity(Vector3r::Zero());
			rb[i]->setVelocity0(Vector3r::Zero());
			rb[i]->setAngularVelocity(Vector3r::Zero());
			rb[i]->setAngularVelocity0(Vector3r::Zero());
		}
		else
		{
			rb[i]->setVelocity(rbd.m_v);
			rb[i]->setVelocity0(rbd.m_v);
			rb[i]->setAngularVelocity(rbd.m_omega);			
			rb[i]->setAngularVelocity0(rbd.m_omega);
		}
		rb[i]->setRestitutionCoeff(rbd.m_restitutionCoeff);
		rb[i]->setFrictionCoeff(rbd.m_frictionCoeff);

		const std::vector<Vector3r>& vertices = rb[i]->getGeometry().getVertexDataLocal().getVertices();
		const unsigned int nVert = static_cast<unsigned int>(vertices.size());

		switch (rbd.m_collisionObjectType)
		{
		case Utilities::SceneLoader::No_Collision_Object: break;
		case Utilities::SceneLoader::Sphere:
			m_cd.addCollisionSphere(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, vertices.data(), nVert, rbd.m_collisionObjectScale[0], rbd.m_testMesh, rbd.m_invertSDF);
			break;
		case Utilities::SceneLoader::Box:
			m_cd.addCollisionBox(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, vertices.data(), nVert, rbd.m_collisionObjectScale, rbd.m_testMesh, rbd.m_invertSDF);
			break;
		case Utilities::SceneLoader::Cylinder:
			m_cd.addCollisionCylinder(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, vertices.data(), nVert, rbd.m_collisionObjectScale.head<2>(), rbd.m_testMesh, rbd.m_invertSDF);
			break;
		case Utilities::SceneLoader::Torus:
			m_cd.addCollisionTorus(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, vertices.data(), nVert, rbd.m_collisionObjectScale.head<2>(), rbd.m_testMesh, rbd.m_invertSDF);
			break;
		case Utilities::SceneLoader::HollowSphere:
			m_cd.addCollisionHollowSphere(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, vertices.data(), nVert, rbd.m_collisionObjectScale[0], rbd.m_thicknessSDF, rbd.m_testMesh, rbd.m_invertSDF);
			break;
		case Utilities::SceneLoader::HollowBox:
			m_cd.addCollisionHollowBox(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, vertices.data(), nVert, rbd.m_collisionObjectScale, rbd.m_thicknessSDF, rbd.m_testMesh, rbd.m_invertSDF);
			break;
		case Utilities::SceneLoader::SDF:
		{
			if (rbd.m_collisionObjectFileName == "")
			{
				const std::string basePath = Utilities::FileSystem::getFilePath(sceneFileName);
				const string cachePath = basePath + "/Cache";
				const string resStr = to_string(rbd.m_resolutionSDF[0]) + "_" + to_string(rbd.m_resolutionSDF[1]) + "_" + to_string(rbd.m_resolutionSDF[2]);
				const std::string modelFileName = Utilities::FileSystem::getFileNameWithExt(rbd.m_modelFile);
				const string sdfFileName = Utilities::FileSystem::normalizePath(cachePath + "/" + modelFileName + "_" + resStr + ".csdf");
				m_cd.addCubicSDFCollisionObject(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, vertices.data(), nVert, distanceFields[sdfFileName], rbd.m_collisionObjectScale, rbd.m_testMesh, rbd.m_invertSDF);
			}
			else
				m_cd.addCubicSDFCollisionObject(i, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, vertices.data(), nVert, distanceFields[rbd.m_collisionObjectFileName], rbd.m_collisionObjectScale, rbd.m_testMesh, rbd.m_invertSDF);
			break;
		}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// triangle models
	//////////////////////////////////////////////////////////////////////////

	// map file names to loaded geometry to prevent multiple imports of same files
	std::map<std::string, pair<PBD::VertexData, Utilities::IndexedFaceMesh>> triFiles;
	for (unsigned int i = 0; i < data.m_triangleModelData.size(); i++)
	{
		const Utilities::SceneLoader::TriangleModelData &tmd = data.m_triangleModelData[i];

		// Check if already loaded
		if (triFiles.find(tmd.m_modelFile) == triFiles.end())
		{
			Utilities::IndexedFaceMesh mesh;
			PBD::VertexData vd;
			loadObj(Utilities::FileSystem::normalizePath(tmd.m_modelFile), vd, mesh, Vector3r::Ones());
			triFiles[tmd.m_modelFile] = { vd, mesh };
		}
	}

	triModels.reserve(data.m_triangleModelData.size());
	std::map<unsigned int, unsigned int> tm_id_index;
	for (unsigned int i = 0; i < data.m_triangleModelData.size(); i++)
	{
		const Utilities::SceneLoader::TriangleModelData &tmd = data.m_triangleModelData[i];

		if (triFiles.find(tmd.m_modelFile) == triFiles.end())
			continue;

		tm_id_index[tmd.m_id] = i;

		PBD::VertexData vd = triFiles[tmd.m_modelFile].first;
		Utilities::IndexedFaceMesh &mesh = triFiles[tmd.m_modelFile].second;

		const Matrix3r R = tmd.m_q.matrix();
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
	std::map<pair<string, string>, pair<vector<Vector3r>, vector<unsigned int>>> tetFiles;
	for (unsigned int i = 0; i < data.m_tetModelData.size(); i++)
	{
		const Utilities::SceneLoader::TetModelData &tmd = data.m_tetModelData[i];

		// Check if already loaded
		pair<string, string> fileNames = { tmd.m_modelFileNodes, tmd.m_modelFileElements };
		if (tetFiles.find(fileNames) == tetFiles.end())
		{
			vector<Vector3r> vertices;
			vector<unsigned int> tets;
			Utilities::TetGenLoader::loadTetgenModel(Utilities::FileSystem::normalizePath(tmd.m_modelFileNodes), Utilities::FileSystem::normalizePath(tmd.m_modelFileElements), vertices, tets);
			tetFiles[fileNames] = { vertices, tets };
		}
	}

	tetModels.reserve(data.m_tetModelData.size());
	std::map<unsigned int, unsigned int> tm_id_index2;
	for (unsigned int i = 0; i < data.m_tetModelData.size(); i++)
	{
		const Utilities::SceneLoader::TetModelData &tmd = data.m_tetModelData[i];

		pair<string, string> fileNames = { tmd.m_modelFileNodes, tmd.m_modelFileElements };
		auto geo = tetFiles.find(fileNames);
		if (geo == tetFiles.end())
			continue;

		tm_id_index2[tmd.m_id] = i;

		vector<Vector3r> vertices = geo->second.first;
		vector<unsigned int> &tets = geo->second.second;

		const Matrix3r R = tmd.m_q.matrix();
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
				Utilities::IndexedFaceMesh &visMesh = tm->getVisMesh();
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
		m_cd.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert, true);

	}
	for (unsigned int i = 0; i < data.m_tetModelData.size(); i++)
	{
		PBD::TetModel *tm = tetModels[i];
		unsigned int offset = tm->getIndexOffset();
		const unsigned int nVert = tm->getParticleMesh().numVertices();
		m_cd.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TetModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
	}

	//////////////////////////////////////////////////////////////////////////
	// joints
	//////////////////////////////////////////////////////////////////////////

	for (unsigned int i = 0; i < data.m_ballJointData.size(); i++)
	{
		const Utilities::SceneLoader::BallJointData &jd = data.m_ballJointData[i];
		m_model.addBallJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position);
	}

	for (unsigned int i = 0; i < data.m_ballOnLineJointData.size(); i++)
	{
		const Utilities::SceneLoader::BallOnLineJointData &jd = data.m_ballOnLineJointData[i];
		m_model.addBallOnLineJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
	}

	for (unsigned int i = 0; i < data.m_hingeJointData.size(); i++)
	{
		const Utilities::SceneLoader::HingeJointData &jd = data.m_hingeJointData[i];
		m_model.addHingeJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
	}

	for (unsigned int i = 0; i < data.m_universalJointData.size(); i++)
	{
		const Utilities::SceneLoader::UniversalJointData &jd = data.m_universalJointData[i];
		m_model.addUniversalJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis[0], jd.m_axis[1]);
	}

	for (unsigned int i = 0; i < data.m_sliderJointData.size(); i++)
	{
		const Utilities::SceneLoader::SliderJointData &jd = data.m_sliderJointData[i];
		m_model.addSliderJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_axis);
	}

	for (unsigned int i = 0; i < data.m_rigidBodyParticleBallJointData.size(); i++)
	{
		const Utilities::SceneLoader::RigidBodyParticleBallJointData &jd = data.m_rigidBodyParticleBallJointData[i];
		m_model.addRigidBodyParticleBallJoint(id_index[jd.m_bodyID[0]], jd.m_bodyID[1]);
	}

	for (unsigned int i = 0; i < data.m_targetAngleMotorHingeJointData.size(); i++)
	{
		const Utilities::SceneLoader::TargetAngleMotorHingeJointData &jd = data.m_targetAngleMotorHingeJointData[i];
		m_model.addTargetAngleMotorHingeJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setTarget(jd.m_target);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setTargetSequence(jd.m_targetSequence);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setRepeatSequence(jd.m_repeat);
	}

	for (unsigned int i = 0; i < data.m_targetVelocityMotorHingeJointData.size(); i++)
	{
		const Utilities::SceneLoader::TargetVelocityMotorHingeJointData &jd = data.m_targetVelocityMotorHingeJointData[i];
		m_model.addTargetVelocityMotorHingeJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setTarget(jd.m_target);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setTargetSequence(jd.m_targetSequence);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setRepeatSequence(jd.m_repeat);
	}

	for (unsigned int i = 0; i < data.m_targetPositionMotorSliderJointData.size(); i++)
	{
		const Utilities::SceneLoader::TargetPositionMotorSliderJointData &jd = data.m_targetPositionMotorSliderJointData[i];
		m_model.addTargetPositionMotorSliderJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_axis);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setTarget(jd.m_target);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setTargetSequence(jd.m_targetSequence);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setRepeatSequence(jd.m_repeat);
	}

	for (unsigned int i = 0; i < data.m_targetVelocityMotorSliderJointData.size(); i++)
	{
		const Utilities::SceneLoader::TargetVelocityMotorSliderJointData &jd = data.m_targetVelocityMotorSliderJointData[i];
		m_model.addTargetVelocityMotorSliderJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_axis);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setTarget(jd.m_target);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setTargetSequence(jd.m_targetSequence);
		((PBD::MotorJoint*)constraints[constraints.size() - 1])->setRepeatSequence(jd.m_repeat);
	}

	for (unsigned int i = 0; i < data.m_damperJointData.size(); i++)
	{
		const Utilities::SceneLoader::DamperJointData &jd = data.m_damperJointData[i];
		m_model.addDamperJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_axis, jd.m_stiffness);
	}

	for (unsigned int i = 0; i < data.m_rigidBodySpringData.size(); i++)
	{
		const Utilities::SceneLoader::RigidBodySpringData &jd = data.m_rigidBodySpringData[i];
		m_model.addRigidBodySpring(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position1, jd.m_position2, jd.m_stiffness);
	}

	for (unsigned int i = 0; i < data.m_distanceJointData.size(); i++)
	{
		const Utilities::SceneLoader::DistanceJointData &jd = data.m_distanceJointData[i];
		m_model.addDistanceJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position1, jd.m_position2);
	}

	m_cd.updateAABBs(m_model);
}

void PBDWrapper::initModel (const Real timeStepSize)
{
	PBD::TimeManager::getCurrent ()->setTimeStepSize(timeStepSize);
 
	PBD::Simulation::getCurrent()->getTimeStep()->setCollisionDetection(m_model, &m_cd);
}

PBD::TimeStepController& PBDWrapper::getTimeStepController()
{
	return *(PBD::TimeStepController*) PBD::Simulation::getCurrent()->getTimeStep();
}
    
 void PBDWrapper::initTriangleModelConstraints()
 {
	// init constraints
	for (unsigned int cm = 0; cm < m_model.getTriangleModels().size(); cm++)
	{
		m_distanceStiffness = 1.0;
		if (m_clothSimulationMethod == 4)
			m_distanceStiffness = 100000;
		m_model.addClothConstraints(m_model.getTriangleModels()[cm], m_clothSimulationMethod, m_distanceStiffness, m_xxStiffness,
			m_yyStiffness, m_xyStiffness, m_xyPoissonRatio, m_yxPoissonRatio, m_normalizeStretch, m_normalizeShear);

		m_bendingStiffness = 0.01;
		if (m_bendingMethod == 3)
			m_bendingStiffness = 100.0;
		m_model.addBendingConstraints(m_model.getTriangleModels()[cm], m_bendingMethod, m_bendingStiffness);
	}
 }
 
 void PBDWrapper::initTetModelConstraints()
 {
	 // init constraints
	 m_solidStiffness = 1.0;
	 if (m_solidSimulationMethod == 5)
		 m_solidStiffness = 100000;

	 m_volumeStiffness = 1.0;
	 if (m_solidSimulationMethod == 5)
		 m_volumeStiffness = 100000;
	 for (unsigned int cm = 0; cm < m_model.getTetModels().size(); cm++)
	 {
		 m_model.addSolidConstraints(m_model.getTetModels()[cm], m_solidSimulationMethod, m_solidStiffness,
			 m_solidPoissonRatio, m_volumeStiffness, m_solidNormalizeStretch, m_solidNormalizeShear);
	 }
 }
 