#include "SPlisHSPlasH/Common.h"
#include "GL/glew.h"
#include "Visualization/MiniGL.h"
#include "GL/glut.h"
#include "SPlisHSPlasH/TimeManager.h"
#include <Eigen/Dense>
#include <iostream>
#include "Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "SPlisHSPlasH/StaticRigidBody.h"
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/Utilities/PoissonDiskSampling.h"
#include "Demos/Common/DemoBase.h"
#include "Utilities/FileSystem.h"
#include "Utilities/Version.h"
#include <fstream>
#include "Utilities/Logger.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Demos/Common/TweakBarParameters.h"

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
	#define new DEBUG_NEW 
#endif
#endif

using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;
using namespace GenParam;

void timeStep ();
void initBoundaryData();
void render ();
void renderBoundary();
void reset();
void initParameters();
void loadObj(const std::string &filename, TriangleMesh &geo, const Vector3r &scale);
void TW_CALL setCurrentFluidModel(const void *value, void *clientData);
void TW_CALL getCurrentFluidModel(void *value, void *clientData);

DemoBase *base;
unsigned int currentFluidModel = 0;

// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;

	base = new DemoBase();
	base->init(argc, argv, "StaticBoundaryDemo");

	Simulation *sim = Simulation::getCurrent();
	sim->init(base->getScene().particleRadius);
	
	base->buildModel();

	initBoundaryData();
	unsigned int nBoundaryParticles = 0;
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
		nBoundaryParticles += sim->getBoundaryModel(i)->numberOfParticles();

	LOG_INFO << "Number of boundary particles: " << nBoundaryParticles;

	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		const std::string &key = model->getId();
		model->setDragMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getDragBase()); });
		model->setSurfaceMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getSurfaceTensionBase()); });
		model->setViscosityMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getViscosityBase()); });
		model->setVorticityMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getVorticityBase()); });
	}

	initParameters();
	base->readParameters();

	Simulation::getCurrent()->setSimulationMethodChangedCallback([&]() { reset(); initParameters(); base->getSceneLoader()->readParameterObject("Configuration", Simulation::getCurrent()->getTimeStep()); });

	MiniGL::setClientIdleFunc(50, timeStep);
	MiniGL::setKeyFunc(0, 'r', reset);
	MiniGL::setClientSceneFunc(render);

	glutMainLoop ();	

	base->cleanup ();

	Utilities::Timing::printAverageTimes();
	Utilities::Timing::printTimeSums();

	Utilities::Counting::printAverageCounts();
	Utilities::Counting::printCounterSums();

	delete Simulation::getCurrent();
	delete base;
	
	return 0;
}

void initParameters()
{
	TwRemoveAllVars(MiniGL::getTweakBar());
	TweakBarParameters::cleanup();

	MiniGL::initTweakBarParameters();

	Simulation *sim = Simulation::getCurrent();
	TweakBarParameters::createParameterGUI();
	TweakBarParameters::createParameterObjectGUI(base);
	TweakBarParameters::createParameterObjectGUI(sim);
	TweakBarParameters::createParameterObjectGUI(sim->getTimeStep());

	// Enum for all fluid models
	if (sim->numberOfFluidModels() > 0)
	{
		TwType enumType = TwDefineEnum("CurrentFluidModelEnum", NULL, 0);
		std::ostringstream oss;
		oss << 0 << " {" << sim->getFluidModel(0)->getId().c_str() << "}";
		for (unsigned int j = 1; j < sim->numberOfFluidModels(); j++)
		{
			oss << ", " << j << " {" << sim->getFluidModel(j)->getId().c_str() << "}";
		}
		std::string enumStr = " label='Current fluid model' enum='" + oss.str() + "' group='Fluid model'";
		TwAddVarCB(MiniGL::getTweakBar(), "CurrentFluidModel", enumType, setCurrentFluidModel, getCurrentFluidModel, &currentFluidModel, enumStr.c_str());
	}

	// show GUI only for currently selected fluid model
	unsigned int i = currentFluidModel;
	FluidModel *model = sim->getFluidModel(currentFluidModel);
	TweakBarParameters::createParameterObjectGUI(model);
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getDragBase());
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getSurfaceTensionBase());
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getViscosityBase());
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getVorticityBase());
	TwDefine((std::string("TweakBar/FluidModel group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/'Drag force' group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/'Surface tension' group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/Viscosity group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/Vorticity group='") + model->getId() + "'").c_str());
}

void reset()
{
	Utilities::Timing::printAverageTimes();
	Utilities::Timing::reset();

	Simulation::getCurrent()->reset();
	base->getSelectedParticles().clear();
}

void timeStep ()
{
	const Real pauseAt = base->getValue<Real>(DemoBase::PAUSE_AT);
	if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime()))
		base->setValue(DemoBase::PAUSE, true);

	if (base->getValue<bool>(DemoBase::PAUSE))
		return;

	// Simulation code
	const unsigned int numSteps = base->getValue<unsigned int>(DemoBase::NUM_STEPS_PER_RENDER);
	for (unsigned int i = 0; i < numSteps; i++)
	{
		START_TIMING("SimStep");
		Simulation::getCurrent()->getTimeStep()->step();
		STOP_TIMING_AVG;

		base->step();
	}
}

void render()
{
	float gridColor[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
	MiniGL::drawGrid(gridColor);

	MiniGL::coordinateSystem();
	MiniGL::drawTime(TimeManager::getCurrent()->getTime());

	Simulation *sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);

		float fluidColor[4] = { 0.3f, 0.5f, 0.9f, 1.0f };
		MiniGL::hsvToRgb(0.61f-0.1f*i, 0.66f, 0.9f, fluidColor);
		base->renderFluid(model, fluidColor);
	}
	renderBoundary();
}

void renderBoundary()
{
	Simulation *sim = Simulation::getCurrent();
	Shader &shader = base->getShader();
	Shader &meshShader = base->getMeshShader();
	SceneLoader::Scene &scene = base->getScene();
	const int renderWalls = base->getValue<int>(DemoBase::RENDER_WALLS);
	GLint context_major_version = base->getContextMajorVersion();

	float wallColor[4] = { 0.1f, 0.6f, 0.6f, 1.0f };
	if ((renderWalls == 1) || (renderWalls == 2))
	{
		if (context_major_version > 3)
		{
			shader.begin();
			glUniform3fv(shader.getUniform("color"), 1, &wallColor[0]);
			glEnableVertexAttribArray(0);
			for (int body = sim->numberOfBoundaryModels() - 1; body >= 0; body--)
			{
				if ((renderWalls == 1) || (!scene.boundaryModels[body]->isWall))
				{
					BoundaryModel *bm = sim->getBoundaryModel(body);
					glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, &bm->getPosition(0));
					glDrawArrays(GL_POINTS, 0, bm->numberOfParticles());
				}
			}
			glDisableVertexAttribArray(0);
			shader.end();
		}
		else
		{
			glDisable(GL_LIGHTING);
			glPointSize(4.0f);
 
			glBegin(GL_POINTS);
			for (int body = sim->numberOfBoundaryModels() - 1; body >= 0; body--)
			{
				if ((renderWalls == 1) || (!scene.boundaryModels[body]->isWall))
				{
					BoundaryModel *bm = sim->getBoundaryModel(body);
					for (unsigned int i = 0; i < bm->numberOfParticles(); i++)
					{
						glColor3fv(wallColor);
						glVertex3v(&bm->getPosition(i)[0]);
					}
				}
			}
			glEnd();
			glEnable(GL_LIGHTING);
		}
	}
	else if ((renderWalls == 3) || (renderWalls == 4))
	{
		for (int body = sim->numberOfBoundaryModels() - 1; body >= 0; body--)
		{
			if ((renderWalls == 3) || (!scene.boundaryModels[body]->isWall))
			{
				meshShader.begin();
				glUniform1f(meshShader.getUniform("shininess"), 5.0f);
				glUniform1f(meshShader.getUniform("specular_factor"), 0.2f);

				GLfloat matrix[16];
				glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
				glUniformMatrix4fv(meshShader.getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
				GLfloat pmatrix[16];
				glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
				glUniformMatrix4fv(meshShader.getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);

				glUniform3fv(meshShader.getUniform("surface_color"), 1, wallColor);

				BoundaryModel *bm = sim->getBoundaryModel(body);
				MiniGL::drawMesh(((StaticRigidBody*)bm->getRigidBodyObject())->getGeometry(), wallColor);

				meshShader.end();
			}
		}		
	}
}

void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale)
{
	std::vector<OBJLoader::Vec3f> x;
	std::vector<OBJLoader::Vec3f> normals;
	std::vector<MeshFaceIndices> faces;
	OBJLoader::Vec3f s = { (float) scale[0], (float)scale[1], (float)scale[2] };
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

void initBoundaryData()
{
	std::string scene_path = FileSystem::getFilePath(base->getSceneFile());
	std::string scene_file_name = FileSystem::getFileName(base->getSceneFile());
	SceneLoader::Scene &scene = base->getScene();
	const bool useCache = base->getUseParticleCaching();

	string cachePath = scene_path + "/Cache";

	for (unsigned int i = 0; i < scene.boundaryModels.size(); i++)
	{
		string meshFileName = FileSystem::normalizePath(scene_path + "/" + scene.boundaryModels[i]->meshFile);

		std::string md5FileName = FileSystem::normalizePath(cachePath + "/" + FileSystem::getFileNameWithExt(meshFileName) + ".md5");
		bool md5 = false;
		if (useCache)
		{
			string md5Str = FileSystem::getFileMD5(meshFileName);
			if (FileSystem::fileExists(md5FileName))
				md5 = FileSystem::checkMD5(md5Str, md5FileName);
		}

		std::vector<Vector3r> boundaryParticles;
		if (scene.boundaryModels[i]->samplesFile != "")
		{
			string particleFileName = scene_path + "/" + scene.boundaryModels[i]->samplesFile;
			PartioReaderWriter::readParticles(particleFileName, scene.boundaryModels[i]->translation, scene.boundaryModels[i]->rotation, scene.boundaryModels[i]->scale[0], boundaryParticles);
		}

		StaticRigidBody *rb = new StaticRigidBody();
		TriangleMesh &geo = rb->getGeometry();
		loadObj(meshFileName, geo, scene.boundaryModels[i]->scale);

		if (scene.boundaryModels[i]->samplesFile == "")
		{
			// Cache sampling
			std::string mesh_base_path = FileSystem::getFilePath(scene.boundaryModels[i]->meshFile);
			std::string mesh_file_name = FileSystem::getFileName(scene.boundaryModels[i]->meshFile);

			const string resStr = to_string(scene.boundaryModels[i]->scale[0]) + "_" + to_string(scene.boundaryModels[i]->scale[1]) + "_" + to_string(scene.boundaryModels[i]->scale[2]);
			const string particleFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_sbd_" + std::to_string(scene.particleRadius) + "_" + resStr + ".bgeo");

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
				LOG_INFO << "Surface sampling of " << meshFileName;
				START_TIMING("Poisson disk sampling");
				PoissonDiskSampling sampling;
				sampling.sampleMesh(geo.numVertices(), geo.getVertices().data(), geo.numFaces(), geo.getFaces().data(), scene.particleRadius, 10, 1, boundaryParticles);
				STOP_TIMING_AVG;

				// Cache sampling
				if (useCache && (FileSystem::makeDir(cachePath) == 0))
				{
					LOG_INFO << "Save particle sampling: " << particleFileName;
					PartioReaderWriter::writeParticles(particleFileName, (unsigned int)boundaryParticles.size(), boundaryParticles.data(), nullptr, scene.particleRadius);
					FileSystem::writeMD5File(meshFileName, md5FileName);
				}

				// transform particles
				for (unsigned int j = 0; j < (unsigned int)boundaryParticles.size(); j++)
					boundaryParticles[j] = scene.boundaryModels[i]->rotation * boundaryParticles[j] + scene.boundaryModels[i]->translation;
			}
		}
		for (unsigned int j = 0; j < geo.numVertices(); j++)
			geo.getVertices()[j] = scene.boundaryModels[i]->rotation * geo.getVertices()[j] + scene.boundaryModels[i]->translation;

		geo.updateNormals();
		geo.updateVertexNormals();

		Simulation::getCurrent()->addBoundaryModel(rb, static_cast<unsigned int>(boundaryParticles.size()), &boundaryParticles[0]);
	}
	Simulation::getCurrent()->updateBoundaryPsi();
}

void TW_CALL setCurrentFluidModel(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	*((unsigned int*)clientData) = val;

	initParameters();
}

void TW_CALL getCurrentFluidModel(void *value, void *clientData)
{
	*(unsigned int *)(value) = *((unsigned int*)clientData);
}
