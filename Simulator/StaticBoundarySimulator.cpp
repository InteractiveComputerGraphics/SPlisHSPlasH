#include "StaticBoundarySimulator.h"
#include "SimulatorBase.h"
#include "Utilities/FileSystem.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Utilities/PartioReaderWriter.h"
#include "SPlisHSPlasH/StaticRigidBody.h"
#include "Utilities/Logger.h"
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/Utilities/SurfaceSampling.h"
#include "SPlisHSPlasH/TriangleMesh.h"
#include "Simulator/SceneConfiguration.h"
#include "SPlisHSPlasH/Utilities/MeshImport.h"

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

void StaticBoundarySimulator::initBoundaryData()
{
	const std::string& sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
	std::string scene_path = FileSystem::getFilePath(sceneFile);
	std::string scene_file_name = FileSystem::getFileName(sceneFile);
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
		rb->setIsAnimated(scene.boundaryModels[i]->isAnimated);
		TriangleMesh &geo = rb->getGeometry();
		MeshImport::importMesh(meshFileName, geo, Vector3r::Zero(), Matrix3r::Identity(), scene.boundaryModels[i]->scale);

		std::vector<Vector3r> boundaryParticles;
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		{
			// if a samples file is given, use this one
			if (scene.boundaryModels[i]->samplesFile != "")
			{
				string particleFileName = scene_path + "/" + scene.boundaryModels[i]->samplesFile;
				PartioReaderWriter::readParticles(particleFileName, Vector3r::Zero(), Matrix3r::Identity(), scene.boundaryModels[i]->scale[0], boundaryParticles);
			}
			else		// if no samples file is given, sample the surface model
			{
				// Cache sampling
				std::string mesh_base_path = FileSystem::getFilePath(scene.boundaryModels[i]->meshFile);
				std::string mesh_file_name = FileSystem::getFileName(scene.boundaryModels[i]->meshFile);

				const string resStr = StringTools::real2String(scene.boundaryModels[i]->scale[0]) + "_" + StringTools::real2String(scene.boundaryModels[i]->scale[1]) + "_" + StringTools::real2String(scene.boundaryModels[i]->scale[2]);
				const string modeStr = "_m" + std::to_string(scene.boundaryModels[i]->samplingMode);
				const string particleFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_sb_" + StringTools::real2String(scene.particleRadius) + "_" + resStr + modeStr + ".bgeo");

				// check MD5 if cache file is available
				bool foundCacheFile = false;

				if (useCache)
					foundCacheFile = FileSystem::fileExists(particleFileName);

				if (useCache && foundCacheFile && md5)
				{
					PartioReaderWriter::readParticles(particleFileName, Vector3r::Zero(), Matrix3r::Identity(), 1.0, boundaryParticles);
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
						sampling.sampleMesh(Matrix3r::Identity(), Vector3r::Zero(),
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
				}
			}
		}

		Matrix3r rot = AngleAxisr(scene.boundaryModels[i]->angle, scene.boundaryModels[i]->axis).toRotationMatrix();
		Quaternionr q(rot);
		rb->setPosition0(scene.boundaryModels[i]->translation);
		rb->setPosition(scene.boundaryModels[i]->translation);
		rb->setRotation0(q);
		rb->setRotation(q);

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
		rb->updateMeshTransformation();
	}
}

void StaticBoundarySimulator::deferredInit()
{
	Simulation* sim = Simulation::getCurrent();
	sim->performNeighborhoodSearchSort();
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		m_base->updateBoundaryParticles(true);
		Simulation::getCurrent()->updateBoundaryVolume();
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		m_base->updateDMVelocity();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		m_base->updateVMVelocity();

#ifdef GPU_NEIGHBORHOOD_SEARCH
	// copy the particle data to the GPU
	sim->getNeighborhoodSearch()->update_point_sets();
#endif 
}

void StaticBoundarySimulator::timeStep()
{
	Simulation* sim = Simulation::getCurrent();
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		m_base->updateBoundaryParticles(false);
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		m_base->updateDMVelocity();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		m_base->updateVMVelocity();
}

void StaticBoundarySimulator::reset()
{
	Simulation* sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel* bm = sim->getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic() || bm->getRigidBodyObject()->isAnimated())
		{
			((StaticRigidBody*)bm->getRigidBodyObject())->reset();
		}
	}

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		m_base->updateBoundaryParticles(true);
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		m_base->updateDMVelocity();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		m_base->updateVMVelocity();
}
