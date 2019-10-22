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
#include "SPlisHSPlasH/Utilities/SurfaceSampling.h"
#include "PositionBasedDynamicsWrapper/PBDWrapper.h"
#include "Simulators/Common/SimulatorBase.h"
#include "Utilities/FileSystem.h"
#include "Utilities/Version.h"
#include "Utilities/SystemInfo.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Simulators/Common/TweakBarParameters.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "Discregrid/All"
#include "SPlisHSPlasH/Utilities/GaussQuadrature.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "GL/freeglut_ext.h"
#include "SPlisHSPlasH/Emitter.h"

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
bool timeStepNoGUI();
void initBoundaryData();
void render ();
void renderBoundary();
void reset();
void initParameters();
void TW_CALL setCurrentFluidModel(const void *value, void *clientData);
void TW_CALL getCurrentFluidModel(void *value, void *clientData);
void TW_CALL setColorField(const void *value, void *clientData);
void TW_CALL getColorField(void *value, void *clientData);
void TW_CALL setRenderMaxValue(const void *value, void *clientData);
void TW_CALL getRenderMaxValue(void *value, void *clientData);
void TW_CALL setRenderMinValue(const void *value, void *clientData);
void TW_CALL getRenderMinValue(void *value, void *clientData);
void TW_CALL setColorMapType(const void *value, void *clientData);
void TW_CALL getColorMapType(void *value, void *clientData);

SimulatorBase *base;
PBDWrapper pbdWrapper;
unsigned int currentFluidModel = 0;
std::vector<std::string> colorFieldNames;

// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;

	base = new SimulatorBase();
	base->init(argc, argv, "DynamicBoundarySimulator");

	if (base->getScene().sim2D)
	{
		LOG_ERR << "DynamicBoundarySimulator does not support 2D simulations, use the StaticBoundarySimulator.";
		exit(1);
	}

	Simulation *sim = Simulation::getCurrent();
	sim->init(base->getScene().particleRadius, base->getScene().sim2D);

	// create additional rigid body information for emitters
	const SceneLoader::Scene &scene = base->getScene();
	std::vector<PBDWrapper::RBData> additionalRBs;
	const std::string scene_path = FileSystem::getFilePath(base->getSceneFile());
	// the last boundary models are the ones that were added for the emitters
	for (auto i = 0; i < scene.emitters.size(); i++)
	{
		SceneLoader::EmitterData *ed = scene.emitters[i];
		PBDWrapper::RBData rb;
		rb.x = ed->x;
		rb.R = ed->rotation;
		rb.scale = Emitter::getSize(ed->width, ed->height, ed->type);
		rb.restitution = 0.6;
		rb.friction = 0.1;
		if (ed->type == 0)
		{
			rb.objFile = FileSystem::normalizePath(base->getDataPath() + "/models/EmitterBox.obj");
			rb.collisionType = 2;
		}
		else if (ed->type == 1)
		{
			rb.objFile = FileSystem::normalizePath(base->getDataPath() + "/models/EmitterCylinder.obj");
			rb.collisionType = 5;
		}
		additionalRBs.push_back(rb);
	}

	//////////////////////////////////////////////////////////////////////////
	// PBD
	//////////////////////////////////////////////////////////////////////////
	pbdWrapper.initShader();
	pbdWrapper.readScene(base->getSceneFile(), additionalRBs);

	base->buildModel();

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		unsigned int nBoundaryParticles = 0;
		for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
			nBoundaryParticles += static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(i))->numberOfParticles();

		LOG_INFO << "Number of boundary particles: " << nBoundaryParticles;
	}

	const bool useGUI = base->getUseGUI();

	if (useGUI)
	{
		for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
		{
			FluidModel *model = sim->getFluidModel(i);
			const std::string &key = model->getId();
			model->setDragMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getDragBase()); });
			model->setSurfaceMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getSurfaceTensionBase()); });
			model->setViscosityMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getViscosityBase()); });
			model->setVorticityMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getVorticityBase()); });
			model->setElasticityMethodChangedCallback([&]() { reset(); initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getElasticityBase()); });
		}
	
		initParameters();
		Simulation::getCurrent()->setSimulationMethodChangedCallback([&]() { reset(); initParameters(); base->getSceneLoader()->readParameterObject("Configuration", Simulation::getCurrent()->getTimeStep()); });
	}

	base->readParameters();

	initBoundaryData();

	pbdWrapper.initModel(TimeManager::getCurrent()->getTimeStepSize());

	if (base->getStateFile() != "")
		base->loadState(base->getStateFile());

	if (!useGUI)
	{
		const Real stopAt = base->getValue<Real>(SimulatorBase::STOP_AT);
		if (stopAt < 0.0)
		{
			LOG_ERR << "StopAt parameter must be set when starting without GUI.";
			exit(1);
		}

		while (true)
		{
			if (!timeStepNoGUI())
				break;
		}
	}
	else
	{
		MiniGL::setClientIdleFunc(50, timeStep);
		MiniGL::addKeyFunc('r', reset);
		MiniGL::setClientSceneFunc(render);
	}

	if (useGUI)
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

	colorFieldNames.clear();
	colorFieldNames.resize(model->numberOfFields());
	TwType enumType = TwDefineEnum("ColorFieldEnum", NULL, 0);
	std::ostringstream oss;
	int idx = 0;
	for (unsigned int j = 0; j < model->numberOfFields(); j++)
	{
		const FieldDescription &field = model->getField(j);
		if ((field.type == FieldType::Scalar) || (field.type == FieldType::Vector3))
		{
			if (idx != 0)
				oss << ", ";
			oss << idx << " {" << field.name.c_str() << "}";
			colorFieldNames[idx] = field.name;
			idx++;
		}
	}
	std::string enumStr = " label='Color field' enum='" + oss.str() + "' group='" + model->getId() + "'";
	enumStr = enumStr + " help='Choose vector or scalar field for particle coloring.'";
	TwAddVarCB(MiniGL::getTweakBar(), "ColorField", enumType, setColorField, getColorField, nullptr, enumStr.c_str());

	TwType enumType2 = TwDefineEnum("ColorMapTypeEnum", NULL, 0);
	std::string str = " label='Color map' enum='0 {None}, 1 {Jet}, 2 {Plasma}' group='" + model->getId() + "' help='Choose a color map.'";
	TwAddVarCB(MiniGL::getTweakBar(), "ColorMapType", enumType2, setColorMapType, getColorMapType, nullptr, str.c_str());
	str = " label='Min. value (shader)' step=0.001 precision=3 group='" + model->getId() + "' help='Minimal value used for color-coding the color field in the rendering process.'";
	TwAddVarCB(MiniGL::getTweakBar(), "RenderMinValue", TW_TYPE_REAL, setRenderMinValue, getRenderMinValue, nullptr, str.c_str());
	str = " label='Max. value (shader)' step=0.001 precision=3 group='" + model->getId() + "' help='Maximal value used for color-coding the color field in the rendering process.'";
	TwAddVarCB(MiniGL::getTweakBar(), "RenderMaxValue", TW_TYPE_REAL, setRenderMaxValue, getRenderMaxValue, nullptr, str.c_str());

	TweakBarParameters::createParameterObjectGUI(model);
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getDragBase());
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getSurfaceTensionBase());
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getViscosityBase());
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getVorticityBase());
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getElasticityBase());
	TwDefine((std::string("TweakBar/FluidModel group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/'Drag force' group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/'Surface tension' group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/Viscosity group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/Vorticity group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/'Elasticity' group='") + model->getId() + "'").c_str());

	pbdWrapper.initGUI();
}

void reset()
{
	Utilities::Timing::printAverageTimes();
	Utilities::Timing::reset();

	Utilities::Counting::printAverageCounts();
	Utilities::Counting::reset();

	Simulation::getCurrent()->reset();
	base->reset();

	//////////////////////////////////////////////////////////////////////////
	// PBD
	//////////////////////////////////////////////////////////////////////////
	pbdWrapper.reset();

	Simulation *sim = Simulation::getCurrent();
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		base->updateBoundaryParticles(true);
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		base->updateDMVelocity();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		base->updateVMVelocity();

	base->getSelectedParticles().clear();
}

void timeStep ()
{
	const Real stopAt = base->getValue<Real>(SimulatorBase::STOP_AT);
	if ((stopAt > 0.0) && (stopAt < TimeManager::getCurrent()->getTime()))
		glutLeaveMainLoop();

	const Real pauseAt = base->getValue<Real>(SimulatorBase::PAUSE_AT);
	if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime()))
		base->setValue(SimulatorBase::PAUSE, true);

	if (base->getValue<bool>(SimulatorBase::PAUSE))
		return;

	// Simulation code
	Simulation *sim = Simulation::getCurrent();
	const bool sim2D = sim->is2DSimulation();
	const unsigned int numSteps = base->getValue<unsigned int>(SimulatorBase::NUM_STEPS_PER_RENDER);
	for (unsigned int i = 0; i < numSteps; i++)
	{
		START_TIMING("SimStep");
		SPH::Simulation::getCurrent()->getTimeStep()->step();
		STOP_TIMING_AVG;

		base->updateBoundaryForces();

		//////////////////////////////////////////////////////////////////////////
		// PBD
		//////////////////////////////////////////////////////////////////////////
		START_TIMING("SimStep - PBD");
		pbdWrapper.timeStep();
		STOP_TIMING_AVG;

		Simulation *sim = Simulation::getCurrent();
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			base->updateBoundaryParticles(false);
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			base->updateDMVelocity();
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			base->updateVMVelocity();

		base->step();

		INCREASE_COUNTER("Time step size", TimeManager::getCurrent()->getTimeStepSize());

		// Make sure that particles stay in xy-plane in a 2D simulation
		if (sim2D)
		{
			for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
			{
				FluidModel *model = sim->getFluidModel(i);
				for (unsigned int i = 0; i < model->numActiveParticles(); i++)
				{
					model->getPosition(i)[2] = 0.0;
					model->getVelocity(i)[2] = 0.0;
				}
			}
		}
	}
}

bool timeStepNoGUI()
{
	const Real stopAt = base->getValue<Real>(SimulatorBase::STOP_AT);
	if ((stopAt > 0.0) && (stopAt < TimeManager::getCurrent()->getTime()))
		return false;

	// Simulation code
	Simulation *sim = Simulation::getCurrent();
	const bool sim2D = sim->is2DSimulation();

	START_TIMING("SimStep");
	SPH::Simulation::getCurrent()->getTimeStep()->step();
	STOP_TIMING_AVG;

	base->updateBoundaryForces();

	//////////////////////////////////////////////////////////////////////////
	// PBD
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("SimStep - PBD");
	pbdWrapper.timeStep();
	STOP_TIMING_AVG;

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		base->updateBoundaryParticles(false);
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		base->updateDMVelocity();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		base->updateVMVelocity();

	base->step();

	INCREASE_COUNTER("Time step size", TimeManager::getCurrent()->getTimeStepSize());

	// Make sure that particles stay in xy-plane in a 2D simulation
	if (sim2D)
	{
		for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
		{
			FluidModel *model = sim->getFluidModel(i);
			for (unsigned int i = 0; i < model->numActiveParticles(); i++)
			{
				model->getPosition(i)[2] = 0.0;
				model->getVelocity(i)[2] = 0.0;
			}
		}
	}
	return true;
}

void renderBoundary()
{
	Simulation *sim = Simulation::getCurrent();
	Shader &shader = base->getShaderScalar();
	Shader &meshShader = base->getMeshShader();
	SceneLoader::Scene &scene = base->getScene();
	const int renderWalls = base->getValue<int>(SimulatorBase::RENDER_WALLS);
	GLint context_major_version = base->getContextMajorVersion();

	if (((renderWalls == 1) || (renderWalls == 2)) && 
		(sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012))
	{
		if (context_major_version > 3)
		{
			for (int body = sim->numberOfBoundaryModels() - 1; body >= 0; body--)
			{
				if ((renderWalls == 1) || (!scene.boundaryModels[body]->isWall))
				{
					base->pointShaderBegin(&shader, scene.boundaryModels[body]->color.data(), 0.0, 100000.0);
					glEnableVertexAttribArray(0);

					BoundaryModel_Akinci2012 *bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(body));
					glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &bm->getPosition(0));
					glDrawArrays(GL_POINTS, 0, bm->numberOfParticles());
					glDisableVertexAttribArray(0);

					base->pointShaderEnd(&shader);
				}
			}
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
					BoundaryModel_Akinci2012 *bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(body));
					for (unsigned int i = 0; i < bm->numberOfParticles(); i++)
					{
						glColor3fv(scene.boundaryModels[body]->color.data());
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
	const bool sim2D = Simulation::getCurrent()->is2DSimulation();
	if (sim2D)
		MiniGL::drawGrid_xy(gridColor);
	else
		MiniGL::drawGrid_xz(gridColor);

	MiniGL::coordinateSystem();
	MiniGL::drawTime(TimeManager::getCurrent()->getTime());

	Simulation *sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		float fluidColor[4] = { 0.3f, 0.5f, 0.9f, 1.0f };
		MiniGL::hsvToRgb(0.61f - 0.1f*i, 0.66f, 0.9f, fluidColor);
		base->renderFluid(i, fluidColor);
	}
	renderBoundary();

	//////////////////////////////////////////////////////////////////////////
	// PBD
	//////////////////////////////////////////////////////////////////////////

	PBD::SimulationModel &model = pbdWrapper.getSimulationModel();
	PBD::SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	const int renderWalls = base->getValue<int>(SimulatorBase::RENDER_WALLS);
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
	pbdWrapper.renderSDF();
}

void initBoundaryData()
{
	std::string scene_path = FileSystem::getFilePath(base->getSceneFile());
	std::string scene_file_name = FileSystem::getFileName(base->getSceneFile());
	SceneLoader::Scene &scene = base->getScene();
	const bool useCache = base->getUseParticleCaching();
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

		// if a samples file is given, use this one
		std::vector<Vector3r> boundaryParticles;

		PBD::SimulationModel &model = pbdWrapper.getSimulationModel();
		PBD::SimulationModel::RigidBodyVector &rigidBodies = model.getRigidBodies();
		PBDRigidBody *rb = new PBDRigidBody(rigidBodies[i]);
		PBD::RigidBodyGeometry &geo = rigidBodies[i]->getGeometry();
		Utilities::IndexedFaceMesh &mesh = geo.getMesh();
		PBD::VertexData &vd = geo.getVertexData();

		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		{
			// if a samples file is given, use this one
			if (scene.boundaryModels[i]->samplesFile != "")
			{
				string particleFileName = scene_path + "/" + scene.boundaryModels[i]->samplesFile;
				PartioReaderWriter::readParticles(particleFileName, Vector3r::Zero(), Matrix3r::Identity(), scene.boundaryModels[i]->scale[0], boundaryParticles);

				// transform back to local coordinates
				for (unsigned int j = 0; j < boundaryParticles.size(); j++)
					boundaryParticles[j] = rb->getRotation().transpose() * (rb->getWorldSpaceRotation() * boundaryParticles[j] + rb->getWorldSpacePosition() - rb->getPosition());
			}
			else		// if no samples file is given, sample the surface model
			{
				// Cache sampling
				std::string mesh_base_path = FileSystem::getFilePath(scene.boundaryModels[i]->meshFile);
				std::string mesh_file_name = FileSystem::getFileName(scene.boundaryModels[i]->meshFile);

				const string resStr = base->real2String(scene.boundaryModels[i]->scale[0]) + "_" + base->real2String(scene.boundaryModels[i]->scale[1]) + "_" + base->real2String(scene.boundaryModels[i]->scale[2]);

				const string modeStr = "_m" + std::to_string(scene.boundaryModels[i]->samplingMode);
				const string particleFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_db_" + base->real2String(scene.particleRadius) + "_" + resStr + modeStr + ".bgeo");

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
					const auto samplePoissonDisk = [&]()
					{
						LOG_INFO << "Poisson disk surface sampling of " << meshFileName;
						START_TIMING("Poisson disk sampling");
						PoissonDiskSampling sampling;
						sampling.sampleMesh(mesh.numVertices(), &vd.getPosition(0), mesh.numFaces(), mesh.getFaces().data(), scene.particleRadius, 10, 1, boundaryParticles);
						STOP_TIMING_AVG;
					};
					const auto sampleRegularTriangle = [&]()
					{
						LOG_INFO << "Regular triangle surface sampling of " << meshFileName;
						START_TIMING("Regular triangle sampling");
						RegularTriangleSampling sampling;
						sampling.sampleMesh(mesh.numVertices(), &vd.getPosition(0), mesh.numFaces(), mesh.getFaces().data(), 1.5f * scene.particleRadius, boundaryParticles);
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

					// transform back to local coordinates
					for (unsigned int j = 0; j < boundaryParticles.size(); j++)
						boundaryParticles[j] = rb->getRotation().transpose() * (boundaryParticles[j] - rb->getPosition());

					// Cache sampling
					if (useCache && (FileSystem::makeDir(cachePath) == 0))
					{
						LOG_INFO << "Save particle sampling: " << particleFileName;
						PartioReaderWriter::writeParticles(particleFileName, (unsigned int)boundaryParticles.size(), boundaryParticles.data(), nullptr, 0.0);
					}
				}
			}

			BoundaryModel_Akinci2012 *bm = new BoundaryModel_Akinci2012();
			bm->initModel(rb, static_cast<unsigned int>(boundaryParticles.size()), &boundaryParticles[0]);
			sim->addBoundaryModel(bm);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		{
			BoundaryModel_Koschier2017 *bm = new BoundaryModel_Koschier2017();
			bm->initModel(rb);
			sim->addBoundaryModel(bm);

			PBD::RigidBodyGeometry &geo = rigidBodies[i]->getGeometry();
			Utilities::IndexedFaceMesh &mesh = geo.getMesh();
			PBD::VertexData &vd = geo.getVertexData();
			// transform back 
			std::vector<Vector3r> xLocal;
			xLocal.resize(vd.size());
			for (unsigned int j = 0; j < vd.size(); j++)
				xLocal[j] = rigidBodies[i]->getRotationMatrix().transpose() * (vd.getPosition(j) - rigidBodies[i]->getPosition());

			base->initDensityMap(xLocal, mesh.getFaces(), scene.boundaryModels[i], md5, true, bm);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		{
			BoundaryModel_Bender2019 *bm = new BoundaryModel_Bender2019();
			bm->initModel(rb);
			sim->addBoundaryModel(bm);

			// transform back 
			std::vector<Vector3r> xLocal;
			xLocal.resize(vd.size());
			for (unsigned int j = 0; j < vd.size(); j++)
				xLocal[j] = rb->getRotation().transpose() * (vd.getPosition(j) - rb->getPosition());
			base->initVolumeMap(xLocal, mesh.getFaces(), scene.boundaryModels[i], md5, true, bm);
		}
		if (useCache && !md5)
			FileSystem::writeMD5File(meshFileName, md5FileName);
	}

	Simulation::getCurrent()->performNeighborhoodSearchSort();
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		base->updateBoundaryParticles(true);
		Simulation::getCurrent()->updateBoundaryVolume();
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		base->updateDMVelocity();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		base->updateVMVelocity();
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

void TW_CALL setColorField(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	base->setColorField(currentFluidModel, colorFieldNames[val]);
}

void TW_CALL getColorField(void *value, void *clientData)
{
	const std::string &fieldName = base->getColorField(currentFluidModel);
	unsigned int index = 0;
	for (auto i = 0; i < colorFieldNames.size(); i++)
	{
		if (colorFieldNames[i] == fieldName)
		{
			index = i;
			break;
		}
	}
	*(unsigned int *)(value) = index;
}

void TW_CALL setRenderMaxValue(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	base->setRenderMaxValue(currentFluidModel, val);
}

void TW_CALL getRenderMaxValue(void *value, void *clientData)
{
	*(Real *)(value) = base->getRenderMaxValue(currentFluidModel);
}

void TW_CALL setRenderMinValue(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	base->setRenderMinValue(currentFluidModel, val);
}

void TW_CALL getRenderMinValue(void *value, void *clientData)
{
	*(Real *)(value) = base->getRenderMinValue(currentFluidModel);
}

void TW_CALL setColorMapType(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	base->setColorMapType(currentFluidModel, val);
}

void TW_CALL getColorMapType(void *value, void *clientData)
{
	*(unsigned int *)(value) = base->getColorMapType(currentFluidModel);
}
