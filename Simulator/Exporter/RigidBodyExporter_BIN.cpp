#include "RigidBodyExporter_BIN.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"
#include "Simulator/SceneConfiguration.h"

using namespace SPH;
using namespace Utilities;

RigidBodyExporter_BIN::RigidBodyExporter_BIN(SimulatorBase* base) :
	ExporterBase(base)
{
	m_isFirstFrame = true;
}

RigidBodyExporter_BIN::~RigidBodyExporter_BIN(void)
{
}

void RigidBodyExporter_BIN::init(const std::string& outputPath)
{
	m_exportPath = FileSystem::normalizePath(outputPath + "/rigid_bodies");
}

void RigidBodyExporter_BIN::step(const unsigned int frame)
{
	if (!m_active)
		return;

	std::string fileName = "rb_data_";
	fileName = fileName + std::to_string(frame) + ".bin";
	std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);

	writeRigidBodies(exportFileName);
}

void RigidBodyExporter_BIN::reset()
{
	m_isFirstFrame = true;
}

void RigidBodyExporter_BIN::setActive(const bool active)
{
	ExporterBase::setActive(active);
	if (m_active)
		FileSystem::makeDirs(m_exportPath);
}


void RigidBodyExporter_BIN::writeRigidBodies(const std::string& fileName)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nBoundaryModels = sim->numberOfBoundaryModels();

	const std::string& sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();

	std::string scene_path = FileSystem::getFilePath(sceneFile);

	// check if we have a static model
	bool isStatic = true;
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel* bm = sim->getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic() || bm->getRigidBodyObject()->isAnimated())
		{
			isStatic = false;
			break;
		}
	}

	BinaryFileWriter binWriter;
	if (m_isFirstFrame || !isStatic)
		binWriter.openFile(fileName.c_str());

	if (m_isFirstFrame)
	{
		binWriter.write(nBoundaryModels);

		for (unsigned int i = 0; i < scene.boundaryModels.size(); i++)
		{
			std::string meshFileName = scene.boundaryModels[i]->meshFile;
			if (FileSystem::isRelativePath(meshFileName))
				meshFileName = FileSystem::normalizePath(scene_path + "/" + meshFileName);

			const std::string fileNameMesh = Utilities::FileSystem::getFileNameWithExt(meshFileName);
			binWriter.write(fileNameMesh);
			Eigen::Vector3f s = scene.boundaryModels[i]->scale.template cast<float>();
			binWriter.writeMatrix(s);
			std::string targetFilePath = m_exportPath + "/" + fileNameMesh;
			if (!Utilities::FileSystem::fileExists(targetFilePath))
			{
				Utilities::FileSystem::copyFile(meshFileName, targetFilePath);
			}
			binWriter.write((char)scene.boundaryModels[i]->isWall);
			binWriter.writeMatrix(scene.boundaryModels[i]->color);
		}
	}

	if (m_isFirstFrame || !isStatic)
	{
		for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
		{
			BoundaryModel* bm = sim->getBoundaryModel(i);
			const Vector3r& x = bm->getRigidBodyObject()->getWorldSpacePosition();
			const Eigen::Vector3f x_f = x.template cast<float>();
			binWriter.writeMatrix(x_f);

			const Matrix3r& R = bm->getRigidBodyObject()->getWorldSpaceRotation();
			const Eigen::Matrix3f R_f = R.template cast<float>();
			//const Eigen::Matrix3f RT = R.transpose().template cast<float>();
			binWriter.writeMatrix(R_f);
		}
		binWriter.closeFile();
	}

	m_isFirstFrame = false;
}