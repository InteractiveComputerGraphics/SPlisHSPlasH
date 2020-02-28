#include "StaticBoundarySimulator.h"
#include "SimulatorBase.h"
#include "Utilities/FileSystem.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Utilities/PartioReaderWriter.h"
#include "SPlisHSPlasH/StaticRigidBody.h"
#include "Utilities/Logger.h"
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/Utilities/SurfaceSampling.h"
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/TriangleMesh.h"

using namespace std;
using namespace SPH;
using namespace Utilities;


StaticBoundarySimulator::StaticBoundarySimulator(SimulatorBase *base)
{
	m_base = base;
}

StaticBoundarySimulator::~StaticBoundarySimulator()
{

}


void StaticBoundarySimulator::loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale)
{
	std::vector<OBJLoader::Vec3f> x;
	std::vector<OBJLoader::Vec3f> normals;
	std::vector<MeshFaceIndices> faces;
	OBJLoader::Vec3f s = { (float)scale[0], (float)scale[1], (float)scale[2] };
	OBJLoader::loadObj(filename, &x, &faces, &normals, nullptr, s);

	mesh.release();
	const unsigned int nPoints = (unsigned int)x.size();
	const unsigned int nFaces = (unsigned int)faces.size();
	mesh.initMesh(nPoints, nFaces);
	for (unsigned int i = 0; i < nPoints; i++)
	{
		mesh.addVertex(Vector3r(x[i][0], x[i][1], x[i][2]));
	}
	for (unsigned int i = 0; i < nFaces; i++)
	{
		// Reduce the indices by one
		int posIndices[3];
		for (int j = 0; j < 3; j++)
		{
			posIndices[j] = faces[i].posIndices[j] - 1;
		}

		mesh.addFace(&posIndices[0]);
	}

	LOG_INFO << "Number of triangles: " << nFaces;
	LOG_INFO << "Number of vertices: " << nPoints;
}

void StaticBoundarySimulator::initBoundaryData()
{
	std::string scene_path = FileSystem::getFilePath(m_base->getSceneFile());
	std::string scene_file_name = FileSystem::getFileName(m_base->getSceneFile());
	SceneLoader::Scene &scene = m_base->getScene();
	// no cache for 2D scenes
	// 2D sampling is fast, but storing it would require storing the transformation as well
	const bool useCache = m_base->getUseParticleCaching() && !scene.sim2D;
	Simulation *sim = Simulation::getCurrent();

	string cachePath = scene_path + "/Cache";

	for (unsigned int i = 0; i < scene.boundaryModels.size(); i++)
	{
		string meshFileName = scene.boundaryModels[i]->meshFile;
		if (FileSystem::isRelativePath(meshFileName))
			meshFileName = FileSystem::normalizePath(scene_path + "/" + scene.boundaryModels[i]->meshFile);

		// check if mesh file has changed
		std::string md5FileName = FileSystem::normalizePath(cachePath + "/" + FileSystem::getFileNameWithExt(meshFileName) + ".md5");
		bool md5 = false;
		if (useCache)
		{
			string md5Str = FileSystem::getFileMD5(meshFileName);
			if (FileSystem::fileExists(md5FileName))
				md5 = FileSystem::checkMD5(md5Str, md5FileName);
		}

		StaticRigidBody *rb = new StaticRigidBody();
		TriangleMesh &geo = rb->getGeometry();
		loadObj(meshFileName, geo, scene.boundaryModels[i]->scale);

		std::vector<Vector3r> boundaryParticles;
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		{
			// if a samples file is given, use this one
			if (scene.boundaryModels[i]->samplesFile != "")
			{
				string particleFileName = scene_path + "/" + scene.boundaryModels[i]->samplesFile;
				PartioReaderWriter::readParticles(particleFileName, scene.boundaryModels[i]->translation, scene.boundaryModels[i]->rotation, scene.boundaryModels[i]->scale[0], boundaryParticles);
			}
			else		// if no samples file is given, sample the surface model
			{
				// Cache sampling
				std::string mesh_base_path = FileSystem::getFilePath(scene.boundaryModels[i]->meshFile);
				std::string mesh_file_name = FileSystem::getFileName(scene.boundaryModels[i]->meshFile);

				const string resStr = m_base->real2String(scene.boundaryModels[i]->scale[0]) + "_" + m_base->real2String(scene.boundaryModels[i]->scale[1]) + "_" + m_base->real2String(scene.boundaryModels[i]->scale[2]);
				const string modeStr = "_m" + std::to_string(scene.boundaryModels[i]->samplingMode);
				const string particleFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_sb_" + m_base->real2String(scene.particleRadius) + "_" + resStr + modeStr + ".bgeo");

				// check MD5 if cache file is available
				bool foundCacheFile = false;

				if (useCache)
					foundCacheFile = FileSystem::fileExists(particleFileName);

				if (useCache && foundCacheFile && md5)
				{
					PartioReaderWriter::readParticles(particleFileName, scene.boundaryModels[i]->translation, scene.boundaryModels[i]->rotation, 1.0, boundaryParticles);
					LOG_INFO << "Loaded cached boundary sampling: " << particleFileName;
				}

				if (!useCache || !foundCacheFile || !md5)
				{
					if (!scene.sim2D)
					{
						const auto samplePoissonDisk = [&]()
						{
							LOG_INFO << "Poisson disk surface sampling of " << meshFileName;
							START_TIMING("Poisson disk sampling");
							PoissonDiskSampling sampling;
							sampling.sampleMesh(geo.numVertices(), geo.getVertices().data(), geo.numFaces(), geo.getFaces().data(), scene.particleRadius, 10, 1, boundaryParticles);
							STOP_TIMING_AVG;
						};
						const auto sampleRegularTriangle = [&]()
						{
							LOG_INFO << "Regular triangle surface sampling of " << meshFileName;
							START_TIMING("Regular triangle sampling");
							RegularTriangleSampling sampling;
							sampling.sampleMesh(geo.numVertices(), geo.getVertices().data(), geo.numFaces(), geo.getFaces().data(), 1.5f * scene.particleRadius, boundaryParticles);
							STOP_TIMING_AVG;
						};
						if (SurfaceSamplingMode::PoissonDisk == scene.boundaryModels[i]->samplingMode)
							samplePoissonDisk();
						else if (SurfaceSamplingMode::RegularTriangle == scene.boundaryModels[i]->samplingMode)
							sampleRegularTriangle();
						else
						{
							LOG_WARN << "Unknown surface sampling method: " << scene.boundaryModels[i]->samplingMode;
							LOG_WARN << "Falling back to:";
							sampleRegularTriangle();
						}
					}
					else
					{
						LOG_INFO << "2D regular sampling of " << meshFileName;
						START_TIMING("2D regular sampling");
						RegularSampling2D sampling;
						sampling.sampleMesh(scene.boundaryModels[i]->rotation, scene.boundaryModels[i]->translation,
							geo.numVertices(), geo.getVertices().data(), geo.numFaces(),
							geo.getFaces().data(), 1.75f * scene.particleRadius, boundaryParticles);
						STOP_TIMING_AVG;
					}

					// Cache sampling
					if (useCache && (FileSystem::makeDir(cachePath) == 0))
					{
						LOG_INFO << "Save particle sampling: " << particleFileName;
						PartioReaderWriter::writeParticles(particleFileName, (unsigned int)boundaryParticles.size(), boundaryParticles.data(), nullptr, scene.particleRadius);
					}

					// transform particles
					if (!scene.sim2D)
						for (unsigned int j = 0; j < (unsigned int)boundaryParticles.size(); j++)
							boundaryParticles[j] = scene.boundaryModels[i]->rotation * boundaryParticles[j] + scene.boundaryModels[i]->translation;
				}
			}
		}

		rb->setWorldSpacePosition(scene.boundaryModels[i]->translation);
		rb->setWorldSpaceRotation(scene.boundaryModels[i]->rotation);
		rb->setPosition(scene.boundaryModels[i]->translation);
		rb->setRotation(scene.boundaryModels[i]->rotation);

		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		{
			BoundaryModel_Akinci2012 *bm = new BoundaryModel_Akinci2012();
			bm->initModel(rb, static_cast<unsigned int>(boundaryParticles.size()), &boundaryParticles[0]);
			sim->addBoundaryModel(bm);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		{
			BoundaryModel_Koschier2017 *bm = new BoundaryModel_Koschier2017();
			bm->initModel(rb);
			sim->addBoundaryModel(bm);
			SPH::TriangleMesh &mesh = rb->getGeometry();
			m_base->initDensityMap(mesh.getVertices(), mesh.getFaces(), scene.boundaryModels[i], md5, false, bm);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		{
			BoundaryModel_Bender2019 *bm = new BoundaryModel_Bender2019();
			bm->initModel(rb);
			sim->addBoundaryModel(bm);
			SPH::TriangleMesh &mesh = rb->getGeometry();
			m_base->initVolumeMap(mesh.getVertices(), mesh.getFaces(), scene.boundaryModels[i], md5, false, bm);
		}
		if (useCache && !md5)
			FileSystem::writeMD5File(meshFileName, md5FileName);
		for (unsigned int j = 0; j < geo.numVertices(); j++)
			geo.getVertices()[j] = scene.boundaryModels[i]->rotation * geo.getVertices()[j] + scene.boundaryModels[i]->translation;

		geo.updateNormals();
		geo.updateVertexNormals();

	}
	sim->performNeighborhoodSearchSort();
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		sim->updateBoundaryVolume();

#ifdef GPU_NEIGHBORHOOD_SEARCH
	// copy the particle data to the GPU
	sim->getNeighborhoodSearch()->update_point_sets();
#endif 
}