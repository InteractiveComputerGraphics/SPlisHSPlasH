#include "SPlisHSPlasH/Common.h"
#include "GL/glew.h"
#include "Visualization/MiniGL.h"
#include "GL/glut.h"
#include "SPlisHSPlasH/TimeManager.h"
#include <Eigen/Dense>
#include <iostream>
#include "Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "PositionBasedDynamicsWrapper/PBDRigidBody.h"
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/Utilities/PoissonDiskSampling.h"
#include "PositionBasedDynamicsWrapper/PBDWrapper.h"
#include "Demos/Common/DemoBase.h"
#include "Utilities/FileSystem.h"
#include "Utilities/Version.h"
#include "Utilities/SystemInfo.h"
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
void updateBoundaryParticles(const bool forceUpdate);
void updateBoundaryForces();
void TW_CALL setCurrentFluidModel(const void *value, void *clientData);
void TW_CALL getCurrentFluidModel(void *value, void *clientData);

DemoBase *base;
PBDWrapper pbdWrapper;
unsigned int currentFluidModel = 0;

// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;

	base = new DemoBase();
	base->init(argc, argv, "DynamicBoundaryDemo");

	Simulation *sim = Simulation::getCurrent();
	sim->init(base->getScene().particleRadius);

	//////////////////////////////////////////////////////////////////////////
	// PBD
	//////////////////////////////////////////////////////////////////////////
	pbdWrapper.initShader();
	pbdWrapper.readScene(base->getSceneFile());

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

	pbdWrapper.initModel(TimeManager::getCurrent()->getTimeStepSize());

	MiniGL::setClientIdleFunc(50, timeStep);
	MiniGL::setKeyFunc(0, 'r', reset);
	MiniGL::setClientSceneFunc(render);

	glutMainLoop ();	

	base->cleanup();

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

	pbdWrapper.initGUI();
}

void reset()
{
	Utilities::Timing::printAverageTimes();
	Utilities::Timing::reset();

	Simulation::getCurrent()->reset();

	//////////////////////////////////////////////////////////////////////////
	// PBD
	//////////////////////////////////////////////////////////////////////////
	pbdWrapper.reset();

	updateBoundaryParticles(true);
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

		updateBoundaryForces();

		//////////////////////////////////////////////////////////////////////////
		// PBD
		//////////////////////////////////////////////////////////////////////////
		START_TIMING("SimStep - PBD");
		pbdWrapper.timeStep();
		STOP_TIMING_AVG;

		updateBoundaryParticles(false);

		base->step();
	}
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
		MiniGL::hsvToRgb(0.61f - 0.1f*i, 0.66f, 0.9f, fluidColor);
		base->renderFluid(model, fluidColor);
	}
	renderBoundary();

	//////////////////////////////////////////////////////////////////////////
	// PBD
	//////////////////////////////////////////////////////////////////////////

	PBD::SimulationModel &model = pbdWrapper.getSimulationModel();
	PBD::SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	const int renderWalls = base->getValue<int>(DemoBase::RENDER_WALLS);
	SceneLoader::Scene &scene = base->getScene();
	if ((renderWalls == 3) || (renderWalls == 4))
	{
		for (size_t i = 0; i < rb.size(); i++)
		{
			const PBD::VertexData &vd = rb[i]->getGeometry().getVertexData();
			const Utilities::IndexedFaceMesh &mesh = rb[i]->getGeometry().getMesh();
			if ((renderWalls == 3) || (!scene.boundaryModels[i]->isWall))
			{
				float *col = &scene.boundaryModels[i]->color[0];
				if (!scene.boundaryModels[i]->isWall)
				{
					base->meshShaderBegin(col);
					pbdWrapper.drawMesh(vd, mesh, 0, col);
					base->meshShaderEnd();
				}
				else
				{
					base->meshShaderBegin(col);
					pbdWrapper.drawMesh(vd, mesh, 0, col);
					base->meshShaderEnd();
				}
			}
		}
	}

	pbdWrapper.renderTriangleModels();
	pbdWrapper.renderTetModels();
	pbdWrapper.renderConstraints();
	pbdWrapper.renderBVH();
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
			PartioReaderWriter::readParticles(particleFileName, Vector3r::Zero(), Matrix3r::Identity(), scene.boundaryModels[i]->scale[0], boundaryParticles);
		}

		PBD::SimulationModel &model = pbdWrapper.getSimulationModel();
		PBD::SimulationModel::RigidBodyVector &rigidBodies = model.getRigidBodies();
		PBDRigidBody *rb = new PBDRigidBody(rigidBodies[i]);
		PBD::RigidBodyGeometry &geo = rigidBodies[i]->getGeometry();
		Utilities::IndexedFaceMesh &mesh = geo.getMesh();
		PBD::VertexData &vd = geo.getVertexData();

		if (scene.boundaryModels[i]->samplesFile == "")
		{
			// Cache sampling
			std::string mesh_base_path = FileSystem::getFilePath(scene.boundaryModels[i]->meshFile);
			std::string mesh_file_name = FileSystem::getFileName(scene.boundaryModels[i]->meshFile);
			
			const string resStr = to_string(scene.boundaryModels[i]->scale[0]) + "_" + to_string(scene.boundaryModels[i]->scale[1]) + "_" + to_string(scene.boundaryModels[i]->scale[2]);
			const string particleFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_dbd_" + std::to_string(scene.particleRadius) + "_" + resStr + ".bgeo");

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
				LOG_INFO << "Surface sampling of " << scene.boundaryModels[i]->meshFile;
				START_TIMING("Poisson disk sampling");
				PoissonDiskSampling sampling;
				sampling.sampleMesh(mesh.numVertices(), &vd.getPosition(0), mesh.numFaces(), mesh.getFaces().data(), scene.particleRadius, 10, 1, boundaryParticles);
				STOP_TIMING_AVG;

				// transform back to local coordinates
				for (unsigned int j = 0; j < boundaryParticles.size(); j++)
					boundaryParticles[j] = rb->getRotation().transpose() * (boundaryParticles[j] - rb->getPosition());

				// Cache sampling
				if (useCache && (FileSystem::makeDir(cachePath) == 0))
				{
					LOG_INFO << "Save particle sampling: " << particleFileName;
					PartioReaderWriter::writeParticles(particleFileName, (unsigned int)boundaryParticles.size(), boundaryParticles.data(), nullptr, scene.particleRadius);
					FileSystem::writeMD5File(meshFileName, md5FileName);
				}
			}
		}
		Simulation::getCurrent()->addBoundaryModel(rb, static_cast<unsigned int>(boundaryParticles.size()), &boundaryParticles[0]);
	}
	Simulation::getCurrent()->updateBoundaryPsi();
	updateBoundaryParticles(true);
}


void updateBoundaryParticles(const bool forceUpdate = false)
{
	Simulation *sim = Simulation::getCurrent();
	SceneLoader::Scene &scene = base->getScene();
	const unsigned int nObjects = sim->numberOfBoundaryModels();
	for (unsigned int i = 0; i < nObjects; i++)
	{
		BoundaryModel *bm = sim->getBoundaryModel(i);
		RigidBodyObject *rbo = bm->getRigidBodyObject();
		if (rbo->isDynamic() || forceUpdate)
		{
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static)  
				for (int j = 0; j < (int)bm->numberOfParticles(); j++)
				{
					bm->getPosition(j) = rbo->getRotation() * bm->getPosition0(j) + rbo->getPosition();
					bm->getVelocity(j) = rbo->getAngularVelocity().cross(bm->getPosition(j) - rbo->getPosition()) + rbo->getVelocity();
				}
			}
		}
	}
}

void updateBoundaryForces()
{
	Real h = TimeManager::getCurrent()->getTimeStepSize();
	SceneLoader::Scene &scene = base->getScene();
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nObjects = sim->numberOfBoundaryModels();	
	for (unsigned int i = 0; i < nObjects; i++)
	{
		BoundaryModel *bm = sim->getBoundaryModel(i);
		RigidBodyObject *rbo = bm->getRigidBodyObject();
		if (rbo->isDynamic())
		{
			((PBDRigidBody*)rbo)->updateTimeStepSize();
			Vector3r force, torque;
			force.setZero();
			torque.setZero();

			for (int j = 0; j < (int)bm->numberOfParticles(); j++)
			{
				force += bm->getForce(j);
				torque += (bm->getPosition(j) - rbo->getPosition()).cross(bm->getForce(j));
				bm->getForce(j).setZero();
			}
			rbo->addForce(force);
			rbo->addTorque(torque);
		}
	}
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
