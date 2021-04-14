#include "SimulatorBase.h"
#include "SPlisHSPlasH/Utilities/SceneLoader.h"
#include "SceneConfiguration.h"
#include "Utilities/FileSystem.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/PartioReaderWriter.h"
#include "SPlisHSPlasH/Emitter.h"
#include "SPlisHSPlasH/EmitterSystem.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Vorticity/MicropolarModel_Bender2017.h"
#include "NumericParameter.h"
#include "Utilities/Logger.h"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "Utilities/Version.h"
#include "Utilities/SystemInfo.h"
#include "extern/partio/src/lib/Partio.h"
#include "SPlisHSPlasH/Utilities/GaussQuadrature.h"
#include "SPlisHSPlasH/Utilities/SimpleQuadrature.h"
#include "extern/cxxopts/cxxopts.hpp"
#include "Simulator/GUI/Simulator_GUI_Base.h"
#include "SPlisHSPlasH/Utilities/VolumeSampling.h"
#include "Utilities/OBJLoader.h"
#include "Utilities/BinaryFileReaderWriter.h"
#include "PositionBasedDynamicsWrapper/PBDBoundarySimulator.h"
#include "StaticBoundarySimulator.h"
#include "Exporter/ExporterBase.h"

INIT_LOGGING
INIT_TIMING
INIT_COUNTING

using namespace SPH;
using namespace std;
using namespace GenParam;
using namespace Utilities;

int SimulatorBase::PAUSE = -1;
int SimulatorBase::PAUSE_AT = -1;
int SimulatorBase::STOP_AT = -1;
int SimulatorBase::NUM_STEPS_PER_RENDER = -1;
int SimulatorBase::DATA_EXPORT_FPS = -1;
int SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES = -1;
int SimulatorBase::STATE_EXPORT = -1;
int SimulatorBase::STATE_EXPORT_FPS = -1;
int SimulatorBase::ASYNC_EXPORT = -1;
int SimulatorBase::RENDER_WALLS = -1;
int SimulatorBase::ENUM_WALLS_NONE = -1;
int SimulatorBase::ENUM_WALLS_PARTICLES_ALL = -1;
int SimulatorBase::ENUM_WALLS_PARTICLES_NO_WALLS = -1;
int SimulatorBase::ENUM_WALLS_GEOMETRY_ALL = -1;
int SimulatorBase::ENUM_WALLS_GEOMETRY_NO_WALLS = -1;

 
SimulatorBase::SimulatorBase()
{
	Utilities::logger.addSink(unique_ptr<Utilities::ConsoleSink>(new Utilities::ConsoleSink(Utilities::LogLevel::INFO)));

	m_boundarySimulator = nullptr;
	m_gui = nullptr;
	m_isStaticScene = true;
	m_numberOfStepsPerRenderUpdate = 4;
	m_renderWalls = 4;
	m_doPause = true;
	m_pauseAt = -1.0;
	m_stopAt = -1.0;
	m_useParticleCaching = true;
	m_useGUI = true;
	m_enableRigidBodyVTKExport = false;
	m_enableRigidBodyExport = false;
	m_enableStateExport = false;
	m_enableAsyncExport = false;
	m_framesPerSecond = 25;
	m_framesPerSecondState = 1;
	m_nextFrameTime = 0.0;
	m_nextFrameTimeState = 0.0;
	m_frameCounter = 1;
	m_isFirstFrame = true;
	m_isFirstFrameVTK = true;
	m_firstState = true;
	m_colorField.resize(1, "velocity");
	m_colorMapType.resize(1, 0);
	m_renderMinValue.resize(1, 0.0);
	m_renderMaxValue.resize(1, 5.0);
	m_particleAttributes = "velocity";
	m_timeStepCB = nullptr;
	m_resetCB = nullptr;
#ifdef DL_OUTPUT
	m_nextTiming = 1.0;
#endif
}

SimulatorBase::~SimulatorBase()
{
	Utilities::Timing::printAverageTimes();
	Utilities::Timing::printTimeSums();

	Utilities::Counting::printAverageCounts();
	Utilities::Counting::printCounterSums();

	delete m_boundarySimulator;
	cleanupExporters();
}

void SimulatorBase::initParameters()
{
	ParameterObject::initParameters();

	PAUSE = createBoolParameter("pause", "Pause", &m_doPause);
	setGroup(PAUSE, "General");
	setDescription(PAUSE, "Pause simulation.");
	setHotKey(PAUSE, "space");

	PAUSE_AT = createNumericParameter("pauseAt", "Pause simulation at", &m_pauseAt);
	setGroup(PAUSE_AT, "General");
	setDescription(PAUSE_AT, "Pause simulation at the given time. When the value is negative, the simulation is not paused.");

	STOP_AT = createNumericParameter("stopAt", "Stop simulation at", &m_stopAt);
	setGroup(STOP_AT, "General");
	setDescription(STOP_AT, "Stop simulation at the given time. When the value is negative, the simulation is not stopped.");

	NUM_STEPS_PER_RENDER = createNumericParameter("numberOfStepsPerRenderUpdate", "# time steps / update", &m_numberOfStepsPerRenderUpdate);
	setGroup(NUM_STEPS_PER_RENDER, "Visualization");
	setDescription(NUM_STEPS_PER_RENDER, "Number of simulation steps per rendered frame.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(NUM_STEPS_PER_RENDER))->setMinValue(1);

	RENDER_WALLS = createEnumParameter("renderWalls", "Render walls", &m_renderWalls);
	setGroup(RENDER_WALLS, "Visualization");
	setDescription(RENDER_WALLS, "Make walls visible/invisible.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(RENDER_WALLS));
	enumParam->addEnumValue("None", ENUM_WALLS_NONE);
	enumParam->addEnumValue("Particles (all)", ENUM_WALLS_PARTICLES_ALL);
	enumParam->addEnumValue("Particles (no walls)", ENUM_WALLS_PARTICLES_NO_WALLS);
	enumParam->addEnumValue("Geometry (all)", ENUM_WALLS_GEOMETRY_ALL);
	enumParam->addEnumValue("Geometry (no walls)", ENUM_WALLS_GEOMETRY_NO_WALLS);

	DATA_EXPORT_FPS = createNumericParameter("dataExportFPS", "Export FPS", &m_framesPerSecond);
	setGroup(DATA_EXPORT_FPS, "Export");
	setDescription(DATA_EXPORT_FPS, "Frame rate of partio, vtk and rigid body export.");

	STATE_EXPORT = createBoolParameter("enableStateExport", "Simulation state export", &m_enableStateExport);
	setGroup(STATE_EXPORT, "Export");
	setDescription(STATE_EXPORT, "Enable/disable export of complete simulation state.");

	STATE_EXPORT_FPS = createNumericParameter("stateExportFPS", "State export FPS", &m_framesPerSecondState);
	setGroup(STATE_EXPORT_FPS, "Export");
	setDescription(STATE_EXPORT_FPS, "Frame rate of simulation state export.");

	PARTICLE_EXPORT_ATTRIBUTES = createStringParameter("particleAttributes", "Export attributes", &m_particleAttributes);
	getParameter(PARTICLE_EXPORT_ATTRIBUTES)->setReadOnly(true);
	setGroup(PARTICLE_EXPORT_ATTRIBUTES, "Export");
	setDescription(PARTICLE_EXPORT_ATTRIBUTES, "Attributes that are exported in the partio files (except id and position).");

	ASYNC_EXPORT = createBoolParameter("enableAsyncExport", "Asynchronous export", &m_enableAsyncExport);
	setGroup(ASYNC_EXPORT, "Export");
	setDescription(ASYNC_EXPORT, "Enable/disable asynchronous export of data. The total performance is faster but disable it when measuring the timings of simulation components. \n\n Currently only supported by partio exporter.");

	for (size_t i = 0; i < m_particleExporters.size(); i++)
	{
		m_particleExporters[i].m_id = createBoolParameter(m_particleExporters[i].m_key, m_particleExporters[i].m_name, 
			[i,this]() -> bool { return m_particleExporters[i].m_exporter->getActive(); },
			[i,this](bool active) { m_particleExporters[i].m_exporter->setActive(active); });
		setGroup(m_particleExporters[i].m_id, "Particle exporters");
		setDescription(m_particleExporters[i].m_id, m_particleExporters[i].m_description);
	}

	for (size_t i = 0; i < m_rbExporters.size(); i++)
	{
		m_rbExporters[i].m_id = createBoolParameter(m_rbExporters[i].m_key, m_rbExporters[i].m_name,
			[i, this]() -> bool { return m_rbExporters[i].m_exporter->getActive(); },
			[i, this](bool active) { m_rbExporters[i].m_exporter->setActive(active); });
		setGroup(m_rbExporters[i].m_id, "Rigid body exporters");
		setDescription(m_rbExporters[i].m_id, m_rbExporters[i].m_description);
	}
}

void SimulatorBase::run()
{
	initSimulation();
	runSimulation();
	cleanup();
}

void SimulatorBase::init(std::vector<std::string> argv, const std::string &windowName)
{
	m_argc = static_cast<int>(argv.size());
	m_argv_vec.clear();
	m_argv_vec.reserve(argv.size());
	for (auto & a : argv)
		m_argv_vec.push_back(&a[0]);

	m_argv = m_argv_vec.data();
	init(m_argc, m_argv, windowName);
}

void SimulatorBase::init(int argc, char **argv, const std::string &windowName)
{
	m_argc = argc;
	m_argv = argv;
	m_windowName = windowName;

	createExporters();

	initParameters();
	m_exePath = FileSystem::getProgramPath();
	setUseParticleCaching(true);

	try
	{
		cxxopts::Options options(argv[0], "SPlisHSPlasH - An open-source library for the physically-based simulation of fluids.");
		options
			.positional_help("[scene file]")
			.show_positional_help();

		options.add_options()
			("h,help", "Print help")
			("v,version", "Print version")
			("no-cache", "Disable caching of boundary samples/maps.")
			("state-file", "State file (state_<time>.bin) that should be loaded.", cxxopts::value<std::string>())
			("output-dir", "Output directory for log file and partio files.", cxxopts::value<std::string>())
			("no-initial-pause", "Disable initial pause when starting the simulation.")
			("no-gui", "Disable GUI.")
			("stopAt", "Sets or overwrites the stopAt parameter of the scene.", cxxopts::value<Real>())
			("param", "Sets or overwrites a parameter of the scene.\n\n" 
					  "- Setting a fluid parameter:\n\t<fluid-id>:<parameter-name>:<value>\n"
					  "- Example: --param Fluid:viscosity:0.01\n\n"
					  "- Setting a configuration parameter:\n\t<parameter-name>:<value>\n"
					  "- Example: --param cflMethod:1\n"
					  , cxxopts::value<std::string>())
			;

		options.add_options("invisible")
			("scene-file", "Scene file", cxxopts::value<std::string>());

		options.parse_positional("scene-file");
		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help({ "" }) << std::endl;
			exit(0);
		}

		if (result.count("version"))
		{
			std::cout << SPLISHSPLASH_VERSION << std::endl;
			exit(0);
		}

		if (result.count("no-cache"))
		{
			setUseParticleCaching(false);
		}

		if (result.count("no-gui"))
		{
			setUseGUI(false);
		}

		if (result.count("stopAt"))
		{
			m_stopAt = result["stopAt"].as<Real>();
		}

		if (result.count("param"))
		{
			const string paramStr = result["param"].as<std::string>();
			Utilities::StringTools::tokenize(paramStr, m_paramTokens, ":");

			// add third element to get a unified method for all parameters
			if (m_paramTokens.size() == 2)
				m_paramTokens.insert(m_paramTokens.begin(), "config");
			if (m_paramTokens.size() != 3)
			{
				LOG_ERR << "--param has wrong syntax!";
				LOG_ERR << "Example 1: --param Fluid:viscosity:0.01";
				LOG_ERR << "Example 2: --param cflMethod:1";
				exit(1);
			}
		}

		std::string sceneFile = "";

		if (result.count("scene-file"))
		{
			sceneFile = result["scene-file"].as<std::string>();
			if (FileSystem::isRelativePath(sceneFile))
				sceneFile = FileSystem::normalizePath(m_exePath + "/" + sceneFile);
		}
#ifdef WIN32
		else
		{
			std::string scenePath = FileSystem::normalizePath(m_exePath + "/../data/Scenes/");
			if (!FileSystem::isDirectory(scenePath))
				scenePath = m_exePath;
			std::replace(scenePath.begin(), scenePath.end(), '/', '\\');
			sceneFile = FileSystem::fileDialog(0, scenePath.c_str(), "*.json");

			if (sceneFile == "")
				exit(0);
		}
#endif
		SceneConfiguration::getCurrent()->setSceneFile(sceneFile);

		m_outputPath = FileSystem::normalizePath(getExePath() + "/output/" + FileSystem::getFileName(sceneFile));
		if (result.count("output-dir"))
		{
			m_outputPath = result["output-dir"].as<std::string>();
		}

		if (result.count("no-initial-pause"))
		{
			m_doPause = false;
		}

		setStateFile("");
		if (result.count("state-file"))
		{
			setStateFile(result["state-file"].as<std::string>());
			if (FileSystem::isRelativePath(getStateFile()))
				setStateFile(FileSystem::normalizePath(m_exePath + "/" + getStateFile()));
		}
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}


	std::string logPath = FileSystem::normalizePath(m_outputPath + "/log");
	FileSystem::makeDirs(logPath);
	Utilities::logger.addSink(unique_ptr<Utilities::FileSink>(new Utilities::FileSink(Utilities::LogLevel::DEBUG, logPath + "/SPH_log.txt")));

	LOG_INFO  << "SPlisHSPlasH version: " << SPLISHSPLASH_VERSION;
	LOG_DEBUG << "Git refspec:          " << GIT_REFSPEC;
	LOG_DEBUG << "Git SHA1:             " << GIT_SHA1;
	LOG_DEBUG << "Git status:           " << GIT_LOCAL_STATUS;
	LOG_DEBUG << "Host name:            " << SystemInfo::getHostName();

	if (!getUseParticleCaching())
		LOG_INFO << "Boundary cache disabled.";
	LOG_INFO << "Output directory: " << m_outputPath;

	//////////////////////////////////////////////////////////////////////////
	// read scene
	//////////////////////////////////////////////////////////////////////////
	m_sceneLoader = std::unique_ptr<SceneLoader>(new SceneLoader());
	const std::string& sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
	Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
	if (sceneFile != "")
		m_sceneLoader->readScene(sceneFile.c_str(), scene);
	else
		return;

	//////////////////////////////////////////////////////////////////////////
	// init boundary simulation
	//////////////////////////////////////////////////////////////////////////
	m_isStaticScene = true;
	for (unsigned int i = 0; i < scene.boundaryModels.size(); i++)
	{
		if (scene.boundaryModels[i]->dynamic)
		{
			m_isStaticScene = false;
			break;
		}
	}
	if (m_isStaticScene)
	{
		LOG_INFO << "Initialize static boundary simulation";
		m_boundarySimulator = new StaticBoundarySimulator(this);
	}
	else
	{
		LOG_INFO << "Initialize dynamic boundary simulation";
		m_boundarySimulator = new PBDBoundarySimulator(this);
	}

	initExporters();
}

void SimulatorBase::initSimulation()
{
	const std::string& sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();

#ifdef DL_OUTPUT
		// copy scene files in output so that the simulation can be reproduced
	std::string sceneFilePath = FileSystem::normalizePath(m_outputPath + "/scene");
	FileSystem::makeDirs(sceneFilePath);
	FileSystem::copyFile(sceneFile, sceneFilePath + "/" + FileSystem::getFileNameWithExt(sceneFile));

	std::string modelsFilePath = FileSystem::normalizePath(m_outputPath + "/models");
	FileSystem::makeDirs(modelsFilePath);
	for (unsigned int i = 0; i < scene.boundaryModels.size(); i++)
	{
		std::string meshFileName = scene.boundaryModels[i]->meshFile;
		if (meshFileName != "")
		{
			if (FileSystem::isRelativePath(meshFileName))
				meshFileName = FileSystem::normalizePath(FileSystem::getFilePath(sceneFile) + "/" + meshFileName);
			FileSystem::copyFile(meshFileName, modelsFilePath + "/" + FileSystem::getFileNameWithExt(meshFileName));
		}
		std::string mapFileName = scene.boundaryModels[i]->mapFile;
		if (mapFileName != "")
		{
			if (FileSystem::isRelativePath(mapFileName))
				mapFileName = FileSystem::normalizePath(FileSystem::getFilePath(sceneFile) + "/" + mapFileName);
			FileSystem::copyFile(scene.boundaryModels[i]->mapFile, modelsFilePath + "/" + FileSystem::getFileNameWithExt(scene.boundaryModels[i]->mapFile));
		}
	}

	for (unsigned int i = 0; i < scene.fluidModels.size(); i++)
	{
		std::string samplesFileName = scene.fluidModels[i]->samplesFile;
		if (samplesFileName != "")
		{
			if (FileSystem::isRelativePath(samplesFileName))
				samplesFileName = FileSystem::normalizePath(FileSystem::getFilePath(sceneFile) + "/" + samplesFileName);
			FileSystem::copyFile(samplesFileName, modelsFilePath + "/" + FileSystem::getFileNameWithExt(samplesFileName));
		}
	}

	std::string progFilePath = FileSystem::normalizePath(m_outputPath + "/program");
	FileSystem::makeDirs(progFilePath);
	FileSystem::copyFile(m_argv[0], progFilePath + "/" + FileSystem::getFileNameWithExt(m_argv[0]));
	#endif

	Simulation *sim = Simulation::getCurrent();
	sim->init(scene.particleRadius, scene.sim2D);

	buildModel();

#ifdef USE_DEBUG_TOOLS
	sim->createDebugTools();
#endif

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		unsigned int nBoundaryParticles = 0;
		for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
			nBoundaryParticles += static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(i))->numberOfParticles();

		LOG_INFO << "Number of boundary particles: " << nBoundaryParticles;
	}

	if (m_useGUI)
		m_gui->init(m_argc, m_argv, m_windowName.c_str());

	m_boundarySimulator->init();

	readParameters();
	if (m_useGUI)
	{
		Simulation *sim = Simulation::getCurrent();
		for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
		{
			FluidModel *model = sim->getFluidModel(i);
			const std::string &key = model->getId();
			model->setDragMethodChangedCallback([this, model]() { m_gui->initSimulationParameterGUI(); getSceneLoader()->readMaterialParameterObject(model->getId(), (ParameterObject*)model->getDragBase()); });
			model->setSurfaceMethodChangedCallback([this, model]() { m_gui->initSimulationParameterGUI(); getSceneLoader()->readMaterialParameterObject(model->getId(), (ParameterObject*)model->getSurfaceTensionBase()); });
			model->setViscosityMethodChangedCallback([this, model]() { m_gui->initSimulationParameterGUI(); getSceneLoader()->readMaterialParameterObject(model->getId(), (ParameterObject*)model->getViscosityBase()); });
			model->setVorticityMethodChangedCallback([this, model]() { m_gui->initSimulationParameterGUI(); getSceneLoader()->readMaterialParameterObject(model->getId(), (ParameterObject*)model->getVorticityBase()); });
			model->setElasticityMethodChangedCallback([this, model]() { m_gui->initSimulationParameterGUI(); getSceneLoader()->readMaterialParameterObject(model->getId(), (ParameterObject*)model->getElasticityBase()); });
		}

		m_gui->initSimulationParameterGUI();
		Simulation::getCurrent()->setSimulationMethodChangedCallback([this]() { 
			reset(); 
			m_gui->initSimulationParameterGUI(); 
			getSceneLoader()->readParameterObject("Configuration", Simulation::getCurrent()->getTimeStep()); 
#ifdef USE_DEBUG_TOOLS
			getSceneLoader()->readParameterObject("Configuration", Simulation::getCurrent()->getDebugTools());
#endif
			});
	}
	setCommandLineParameter();
	updateScalarField();

	m_boundarySimulator->initBoundaryData();
}

void SimulatorBase::deferredInit()
{
	Simulation* sim = Simulation::getCurrent();
	if (m_useGUI)
	{
		for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
		{
			FluidModel* model = sim->getFluidModel(i);
			if (model->getDragBase())
				getSceneLoader()->readMaterialParameterObject(model->getId(), (ParameterObject*)model->getDragBase());
			if (model->getSurfaceTensionBase())
				getSceneLoader()->readMaterialParameterObject(model->getId(), (ParameterObject*)model->getSurfaceTensionBase());
			if (model->getViscosityBase())
				getSceneLoader()->readMaterialParameterObject(model->getId(), (ParameterObject*)model->getViscosityBase());
			if (model->getVorticityBase())
				getSceneLoader()->readMaterialParameterObject(model->getId(), (ParameterObject*)model->getVorticityBase());
			if (model->getElasticityBase())
				getSceneLoader()->readMaterialParameterObject(model->getId(), (ParameterObject*)model->getElasticityBase());
			getSceneLoader()->readParameterObject("Configuration", Simulation::getCurrent()->getTimeStep());
		}
		m_gui->initSimulationParameterGUI();
	}
	sim->setSimulationInitialized(true);
	m_boundarySimulator->deferredInit();
}

void SimulatorBase::runSimulation()
{
	deferredInit();

	if (getStateFile() != "")
		loadState(getStateFile());

	if (!m_useGUI)
	{
		const Real stopAt = getValue<Real>(SimulatorBase::STOP_AT);
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

	if (m_useGUI)
		m_gui->run();
}

void SimulatorBase::cleanup()
{
	if (m_useGUI)
		m_gui->cleanup();

	delete SceneConfiguration::getCurrent();
	delete Simulation::getCurrent();
}

void SimulatorBase::readParameters()
{
	m_sceneLoader->readParameterObject("Configuration", this);
	m_sceneLoader->readParameterObject("Configuration", Simulation::getCurrent());
	m_sceneLoader->readParameterObject("Configuration", Simulation::getCurrent()->getTimeStep());
#ifdef USE_DEBUG_TOOLS
	m_sceneLoader->readParameterObject("Configuration", Simulation::getCurrent()->getDebugTools());
#endif

	Simulation *sim = Simulation::getCurrent();
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		const std::string &key = model->getId();
		m_sceneLoader->readMaterialParameterObject(key, model);
		m_sceneLoader->readMaterialParameterObject(key, (ParameterObject*) model->getDragBase());
		m_sceneLoader->readMaterialParameterObject(key, (ParameterObject*) model->getSurfaceTensionBase());
		m_sceneLoader->readMaterialParameterObject(key, (ParameterObject*) model->getViscosityBase());
		m_sceneLoader->readMaterialParameterObject(key, (ParameterObject*) model->getVorticityBase());
		m_sceneLoader->readMaterialParameterObject(key, (ParameterObject*) model->getElasticityBase());

		for (auto material : scene.materials)
		{
			if (material->id == key)
			{
				setColorField(i, material->colorField);
				setColorMapType(i, material->colorMapType);
				setRenderMinValue(i, material->minVal);
				setRenderMaxValue(i, material->maxVal);
			}
		}
	}
}

void SimulatorBase::setCommandLineParameter()
{
	Simulation *sim = Simulation::getCurrent();
	if (m_paramTokens.size() != 3)
		return;

	setCommandLineParameter((ParameterObject*)sim);
	setCommandLineParameter((ParameterObject*)sim->getTimeStep());
#ifdef USE_DEBUG_TOOLS
	setCommandLineParameter((ParameterObject*)sim->getDebugTools());
#endif
	
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		const std::string &key = model->getId();
		
		if (m_paramTokens[0] == key)
		{			
			setCommandLineParameter((ParameterObject*)model);
 			setCommandLineParameter((ParameterObject*)model->getDragBase());
 			setCommandLineParameter((ParameterObject*)model->getSurfaceTensionBase());
			setCommandLineParameter((ParameterObject*)model->getViscosityBase());
 			setCommandLineParameter((ParameterObject*)model->getVorticityBase());
 			setCommandLineParameter((ParameterObject*)model->getElasticityBase());
		}
	}
}

void SimulatorBase::setCommandLineParameter(GenParam::ParameterObject *paramObj)
{
	if (paramObj == nullptr)
		return;

	const unsigned int numParams = paramObj->numParameters();
	for (unsigned int j = 0; j < numParams; j++)
	{
		ParameterBase *paramBase = paramObj->getParameter(j);
		if (m_paramTokens[1] == paramBase->getName())
		{
			if (paramBase->getType() == RealParameterType)
			{
				const Real val = stof(m_paramTokens[2]);
				static_cast<NumericParameter<Real>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::UINT32)
			{
				const unsigned int val = stoi(m_paramTokens[2]);
				static_cast<NumericParameter<unsigned int>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::UINT16)
			{
				const unsigned short val = stoi(m_paramTokens[2]);
				static_cast<NumericParameter<unsigned short>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::UINT8)
			{
				const unsigned char val = stoi(m_paramTokens[2]);
				static_cast<NumericParameter<unsigned char>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::INT32)
			{
				const int val = stoi(m_paramTokens[2]);
				static_cast<NumericParameter<int>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::INT16)
			{
				const short val = stoi(m_paramTokens[2]);
				static_cast<NumericParameter<short>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::INT8)
			{
				const char val = stoi(m_paramTokens[2]);
				static_cast<NumericParameter<char>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::ENUM)
			{
				const int val = stoi(m_paramTokens[2]);
				static_cast<EnumParameter*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::BOOL)
			{
				const bool val = stoi(m_paramTokens[2]);
				static_cast<BoolParameter*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == RealVectorParameterType)
			{
 				if (static_cast<VectorParameter<Real>*>(paramBase)->getDim() == 3)
 				{
					vector<string> tokens;
					Utilities::StringTools::tokenize(m_paramTokens[2], tokens, ",");
					if (tokens.size() == 3)
					{
						Vector3r val(stof(tokens[0]), stof(tokens[1]), stof(tokens[2]));
						static_cast<VectorParameter<Real>*>(paramBase)->setValue(val.data());
					}
 				}
			}
			else if (paramBase->getType() == ParameterBase::STRING)
			{
				const std::string val = m_paramTokens[2];
				static_cast<StringParameter*>(paramBase)->setValue(val);
			}
		}
	}
}

void SimulatorBase::cleanupExporters()
{
	for (size_t i = 0; i < m_particleExporters.size(); i++)
		delete m_particleExporters[i].m_exporter;
	m_particleExporters.clear();
	for (size_t i = 0; i < m_rbExporters.size(); i++)
		delete m_rbExporters[i].m_exporter;
	m_rbExporters.clear();
}

void SimulatorBase::initExporters()
{
	for (size_t i = 0; i < m_particleExporters.size(); i++)
		m_particleExporters[i].m_exporter->init(m_outputPath);
	for (size_t i = 0; i < m_rbExporters.size(); i++)
		m_rbExporters[i].m_exporter->init(m_outputPath);
}


void SimulatorBase::buildModel()
{
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
	TimeManager::getCurrent()->setTimeStepSize(scene.timeStepSize);

	initFluidData();

	createEmitters();
	createAnimationFields();

	Simulation *sim = Simulation::getCurrent();

	if (sim->getTimeStep())
		sim->getTimeStep()->resize();

	if (!sim->is2DSimulation())
	{
		sim->setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC);
		sim->setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}
	else
	{
		sim->setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_CUBIC_2D);
		sim->setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_CUBIC_2D);
	}
}

void SimulatorBase::reset()
{
	Utilities::Timing::printAverageTimes();
	Utilities::Timing::reset();

	Utilities::Counting::printAverageCounts();
	Utilities::Counting::reset();

	Simulation::getCurrent()->reset();
#ifdef USE_DEBUG_TOOLS
	Simulation::getCurrent()->getDebugTools()->reset();
#endif

	m_boundarySimulator->reset();
	if (m_gui)
		m_gui->reset();

	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
	if (Simulation::getCurrent()->getValue<int>(Simulation::CFL_METHOD) != Simulation::ENUM_CFL_NONE)
		TimeManager::getCurrent()->setTimeStepSize(scene.timeStepSize);
	m_nextFrameTime = 0.0;
	m_nextFrameTimeState = 0.0;
	m_frameCounter = 1;
	m_isFirstFrame = true;
	m_isFirstFrameVTK = true;
#ifdef DL_OUTPUT
	m_nextTiming = 1.0;
#endif
	updateScalarField();

	for (size_t i = 0; i < m_particleExporters.size(); i++)
		m_particleExporters[i].m_exporter->reset();
	for (size_t i = 0; i < m_rbExporters.size(); i++)
		m_rbExporters[i].m_exporter->reset();

	if (m_resetCB)
		m_resetCB();
}

void SimulatorBase::singleTimeStep()
{
	const unsigned int numSteps = getValue<unsigned int>(SimulatorBase::NUM_STEPS_PER_RENDER);
	setValue<unsigned int>(SimulatorBase::NUM_STEPS_PER_RENDER, 1);
	setValue<bool>(SimulatorBase::PAUSE, false);

	timeStep();

	setValue<unsigned int>(SimulatorBase::NUM_STEPS_PER_RENDER, numSteps);
	setValue<bool>(SimulatorBase::PAUSE, true);

	const Real stopAt = getValue<Real>(SimulatorBase::STOP_AT);
	if (m_gui && (stopAt > 0.0) && (stopAt < TimeManager::getCurrent()->getTime()))
		m_gui->stop();
}

void SimulatorBase::timeStep()
{
	const Real stopAt = getValue<Real>(SimulatorBase::STOP_AT);
	if (m_gui && (stopAt > 0.0) && (stopAt < TimeManager::getCurrent()->getTime()))
		m_gui->stop();

	const Real pauseAt = getValue<Real>(SimulatorBase::PAUSE_AT);
	if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime()))
		setValue(SimulatorBase::PAUSE, true);

	if (getValue<bool>(SimulatorBase::PAUSE))
		return;

	// Simulation code
	Simulation *sim = Simulation::getCurrent();
	const bool sim2D = sim->is2DSimulation();
	const unsigned int numSteps = getValue<unsigned int>(SimulatorBase::NUM_STEPS_PER_RENDER);
	for (unsigned int i = 0; i < numSteps; i++)
	{
		START_TIMING("SimStep");
		Simulation::getCurrent()->getTimeStep()->step();
		STOP_TIMING_AVG;

		m_boundarySimulator->timeStep();

		step();

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

#ifdef USE_DEBUG_TOOLS
		Simulation::getCurrent()->getDebugTools()->step();
#endif

		if (m_timeStepCB)
			m_timeStepCB();
	}

	updateScalarField();
}

bool SimulatorBase::timeStepNoGUI()
{
	const Real stopAt = getValue<Real>(SimulatorBase::STOP_AT);
	if ((stopAt > 0.0) && (stopAt < TimeManager::getCurrent()->getTime()))
		return false;

	// Simulation code
	Simulation *sim = Simulation::getCurrent();
	const bool sim2D = sim->is2DSimulation();

	START_TIMING("SimStep");
	Simulation::getCurrent()->getTimeStep()->step();
	STOP_TIMING_AVG;

	m_boundarySimulator->timeStep();

	step();

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

#ifdef USE_DEBUG_TOOLS
	Simulation::getCurrent()->getDebugTools()->step();
#endif

	if (m_timeStepCB)
		m_timeStepCB();
	return true;
}

void SimulatorBase::updateScalarField()
{
	Simulation* sim = Simulation::getCurrent();
	m_scalarField.resize(sim->numberOfFluidModels());
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < sim->numberOfFluidModels(); fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		// Draw simulation model
		const unsigned int nParticles = model->numActiveParticles();

		const FieldDescription* field = nullptr;
		field = &model->getField(getColorField(fluidModelIndex));

		if (field != nullptr)
		{
			m_scalarField[fluidModelIndex].resize(model->numActiveParticles());
			for (unsigned int i = 0u; i < model->numActiveParticles(); i++)
			{
				if (field->type == FieldType::Vector3)
				{
					Eigen::Map<Vector3r> vec((Real*)field->getFct(i));
					m_scalarField[fluidModelIndex][i] = static_cast<float>(vec.norm());
				}
				else if (field->type == FieldType::Scalar)
				{
					m_scalarField[fluidModelIndex][i] = static_cast<float>(*(Real*)field->getFct(i));
				}
				else if (field->type == FieldType::UInt)
				{
					m_scalarField[fluidModelIndex][i] = static_cast<float>(*(unsigned int*)field->getFct(i));
				}
				else if (field->type == FieldType::Matrix3)
				{
					Eigen::Map<Matrix3r> m((Real*)field->getFct(i));
					m_scalarField[fluidModelIndex][i] = static_cast<float>(m.norm());
				}
				else if (field->type == FieldType::Vector6)
				{
					Eigen::Map<Vector6r> m((Real*)field->getFct(i));
					m_scalarField[fluidModelIndex][i] = static_cast<float>(m.norm());
				}
				else if (field->type == FieldType::Matrix6)
				{
					Eigen::Map<Matrix6r> m((Real*)field->getFct(i));
					m_scalarField[fluidModelIndex][i] = static_cast<float>(m.norm());
				}
			}
		}
	}
}

void SimulatorBase::loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale)
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

void SimulatorBase::activateExporter(const std::string& exporterName, const bool active)
{
	for (auto i = 0; i < m_particleExporters.size(); i++)
	{
		if (exporterName == m_particleExporters[i].m_name)
		{
			m_particleExporters[i].m_exporter->setActive(active);
			return;
		}
	}
	for (auto i = 0; i < m_rbExporters.size(); i++)
	{
		if (exporterName == m_rbExporters[i].m_name)
		{
			m_rbExporters[i].m_exporter->setActive(active);
			return;
		}
	}
}

void SimulatorBase::setInitialVelocity(const Vector3r& vel, const Vector3r& angVel, const unsigned int numParticles, Vector3r* fluidParticles, Vector3r* fluidVelocities)
{
	Vector3r com;
	com.setZero();
	for (auto i = 0; i < numParticles; i++)
	{
		// all have the same mass
		com += fluidParticles[i];
	}
	com /= (Real)numParticles;
	for (auto i = 0; i < numParticles; i++)
	{
		fluidVelocities[i] = vel + angVel.cross(fluidParticles[i] - com);
	}
}

void SimulatorBase::initFluidData()
{
	LOG_INFO << "Initialize fluid particles";

	Simulation *sim = Simulation::getCurrent();
	const std::string& sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();

	//////////////////////////////////////////////////////////////////////////
	// Determine number of different fluid IDs
	//////////////////////////////////////////////////////////////////////////
	std::map<std::string, unsigned int> fluidIDs;
	unsigned int index = 0;
	for (unsigned int i = 0; i < scene.fluidBlocks.size(); i++)
	{
		if (fluidIDs.find(scene.fluidBlocks[i]->id) == fluidIDs.end())
			fluidIDs[scene.fluidBlocks[i]->id] = index++;
	}
	for (unsigned int i = 0; i < scene.fluidModels.size(); i++)
	{
		if (fluidIDs.find(scene.fluidModels[i]->id) == fluidIDs.end())
			fluidIDs[scene.fluidModels[i]->id] = index++;
	}
	for (unsigned int i = 0; i < scene.emitters.size(); i++)
	{
		if (fluidIDs.find(scene.emitters[i]->id) == fluidIDs.end())
			fluidIDs[scene.emitters[i]->id] = index++;
	}
	const unsigned int numberOfFluidModels = static_cast<unsigned int>(fluidIDs.size());

	std::vector<std::vector<Vector3r>> fluidParticles;
	std::vector<std::vector<Vector3r>> fluidVelocities;
	fluidParticles.resize(numberOfFluidModels);
	fluidVelocities.resize(numberOfFluidModels);

	createFluidBlocks(fluidIDs, fluidParticles, fluidVelocities);

	std::string base_path = FileSystem::getFilePath(sceneFile);

	const bool useCache = getUseParticleCaching();
	std::string scene_path = FileSystem::getFilePath(sceneFile);
	string cachePath = scene_path + "/Cache";

	unsigned int startIndex = 0;
	unsigned int endIndex = 0;
	for (unsigned int i = 0; i < scene.fluidModels.size(); i++)
	{
		const unsigned int fluidIndex = fluidIDs[scene.fluidModels[i]->id];
		const unsigned int startIndex = (unsigned int)fluidParticles[fluidIndex].size();

		std::string fileName;
		if (FileSystem::isRelativePath(scene.fluidModels[i]->samplesFile))
			fileName = base_path + "/" + scene.fluidModels[i]->samplesFile;
		else
			fileName = scene.fluidModels[i]->samplesFile;

		string ext = FileSystem::getFileExt(fileName);
		transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
		if (ext == "OBJ")
		{
			// check if mesh file has changed
			std::string md5FileName = FileSystem::normalizePath(cachePath + "/" + FileSystem::getFileNameWithExt(fileName) + "_fluid.md5");
			bool md5 = false;
			if (useCache)
			{
				string md5Str = FileSystem::getFileMD5(fileName);
				if (FileSystem::fileExists(md5FileName))
					md5 = FileSystem::checkMD5(md5Str, md5FileName);
			}

			// Cache sampling
			std::string mesh_base_path = FileSystem::getFilePath(fileName);
			std::string mesh_file_name = FileSystem::getFileName(fileName);

			std::string invert = to_string(scene.fluidModels[i]->invert);
			std::string mode = to_string(scene.fluidModels[i]->mode);
			const string scaleStr = StringTools::real2String(scene.fluidModels[i]->scale[0]) + "_" + StringTools::real2String(scene.fluidModels[i]->scale[1]) + "_" + StringTools::real2String(scene.fluidModels[i]->scale[2]);
			const string resStr = to_string(scene.fluidModels[i]->resolutionSDF[0]) + "_" + to_string(scene.fluidModels[i]->resolutionSDF[1]) + "_" + to_string(scene.fluidModels[i]->resolutionSDF[2]);
			const string particleFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_fluid_" + StringTools::real2String(scene.particleRadius) + "_i" + invert + "_m" + mode + "_s" + scaleStr + "_r" + resStr + ".bgeo");

			// check MD5 if cache file is available
			bool foundCacheFile = false;

			if (useCache)
				foundCacheFile = FileSystem::fileExists(particleFileName);

			if (useCache && foundCacheFile && md5)
			{
				PartioReaderWriter::readParticles(particleFileName, scene.fluidModels[i]->translation, scene.fluidModels[i]->rotation, 1.0, fluidParticles[fluidIndex], fluidVelocities[fluidIndex]);
				LOG_INFO << "Loaded cached fluid sampling: " << particleFileName;
			}

			if (!useCache || !foundCacheFile || !md5)
			{
				LOG_INFO << "Volume sampling of " << fileName;

				TriangleMesh mesh;
				loadObj(fileName, mesh, scene.fluidModels[i]->scale);

				bool invert = scene.fluidModels[i]->invert;
				int mode = scene.fluidModels[i]->mode;
				std::array<unsigned int, 3> resolutionSDF = scene.fluidModels[i]->resolutionSDF;

				LOG_INFO << "SDF resolution: " << resolutionSDF[0] << ", " << resolutionSDF[1] << ", " << resolutionSDF[2];

				const unsigned int size_before_sampling = fluidParticles[fluidIndex].size();
				START_TIMING("Volume sampling");
				Utilities::VolumeSampling::sampleMesh(mesh.numVertices(), mesh.getVertices().data(), mesh.numFaces(), mesh.getFaces().data(),
					scene.particleRadius, nullptr, resolutionSDF, invert, mode, fluidParticles[fluidIndex]);
				STOP_TIMING_AVG;
                const unsigned int size_after_sampling = fluidParticles[fluidIndex].size();

				fluidVelocities[fluidIndex].resize(fluidParticles[fluidIndex].size());

				// Cache sampling
				if (useCache && (FileSystem::makeDir(cachePath) == 0))
				{
					LOG_INFO << "Save particle sampling: " << particleFileName;
					PartioReaderWriter::writeParticles(particleFileName, size_after_sampling-size_before_sampling, &fluidParticles[fluidIndex][size_before_sampling], &fluidVelocities[fluidIndex][size_before_sampling], 0.0);
                    // PartioReaderWriter::writeParticles(particleFileName, (unsigned int)fluidParticles[fluidIndex].size(), fluidParticles[fluidIndex].data(), fluidVelocities[fluidIndex].data(), 0.0);
					FileSystem::writeMD5File(fileName, md5FileName);
				}

				// transform particles
				for (unsigned int j = size_before_sampling; j < size_after_sampling; j++)
					fluidParticles[fluidIndex][j] = scene.fluidModels[i]->rotation * fluidParticles[fluidIndex][j] + scene.fluidModels[i]->translation;
			}
		}
		else
		{
			if (!PartioReaderWriter::readParticles(fileName, scene.fluidModels[i]->translation, scene.fluidModels[i]->rotation, scene.fluidModels[i]->scale[0], fluidParticles[fluidIndex], fluidVelocities[fluidIndex]))
				LOG_ERR << "File not found: " << fileName;
		}

		const unsigned int numAddedParticles = (unsigned int)fluidParticles[fluidIndex].size() - startIndex;
		setInitialVelocity(scene.fluidModels[i]->initialVelocity, scene.fluidModels[i]->initialAngularVelocity, numAddedParticles, &fluidParticles[fluidIndex][startIndex], &fluidVelocities[fluidIndex][startIndex]);
		Simulation::getCurrent()->setValue(Simulation::PARTICLE_RADIUS, scene.particleRadius);
	}

	unsigned int nParticles = 0;
	for (auto it = fluidIDs.begin(); it != fluidIDs.end(); it++)
	{
		const unsigned int index = it->second;

		unsigned int maxEmitterParticles = 10000;
		for (auto material : scene.materials)
		{
			if (material->id == it->first)
				maxEmitterParticles = material->maxEmitterParticles;
		}
		sim->addFluidModel(it->first, (unsigned int)fluidParticles[index].size(), fluidParticles[index].data(), fluidVelocities[index].data(), maxEmitterParticles);
		nParticles += (unsigned int)fluidParticles[index].size();
	}

	m_colorField.resize(sim->numberOfFluidModels(), "velocity");
	m_colorMapType.resize(sim->numberOfFluidModels(), 0);
	m_renderMinValue.resize(sim->numberOfFluidModels(), 0.0);
	m_renderMaxValue.resize(sim->numberOfFluidModels(), 5.0);

	LOG_INFO << "Number of fluid particles: " << nParticles;
}

void SimulatorBase::createEmitters()
{
	Simulation *sim = Simulation::getCurrent();
	Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();

	//////////////////////////////////////////////////////////////////////////
	// emitters
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int i = 0; i < scene.emitters.size(); i++)
	{
		SceneLoader::EmitterData *ed = scene.emitters[i];

		FluidModel *model = nullptr;
		unsigned int j;
		for (j = 0; j < sim->numberOfFluidModels(); j++)
		{
			model = sim->getFluidModel(j);
			if (model->getId() == ed->id)
				break;
		}

		if (j < sim->numberOfFluidModels())
		{
			if (sim->is2DSimulation())
			{
				// in 2D convert all emitters to box emitters
				// and set width to 1
				if (ed->type == 1)
					ed->height = ed->width;
				ed->width = 1;
				ed->type = 0;
			}
			model->getEmitterSystem()->addEmitter(
				ed->width, ed->height,
				ed->x, ed->rotation,
				ed->velocity,
				ed->type);

			// Generate boundary geometry around emitters
			Emitter *emitter = model->getEmitterSystem()->getEmitters().back();
			SceneLoader::BoundaryData *emitterBoundary = new SceneLoader::BoundaryData();
			emitterBoundary->dynamic = false;
			emitterBoundary->isWall = false;
			emitterBoundary->color = { 0.2f, 0.2f, 0.2f, 1.0f };
			emitterBoundary->rotation = ed->rotation;
			const Real supportRadius = sim->getSupportRadius();
			const Vector3r & emitDir = ed->rotation.col(0);
			emitterBoundary->scale = Emitter::getSize(static_cast<Real>(ed->width), static_cast<Real>(ed->height), ed->type);
			const Vector3r pos = ed->x;
			emitterBoundary->translation = pos;
			emitterBoundary->samplesFile = "";
			emitterBoundary->mapInvert = false;
			emitterBoundary->mapResolution = Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign>(20, 20, 20);
			emitterBoundary->mapThickness = 0.0;

			if (sim->is2DSimulation())
				emitterBoundary->scale[2] = 2 * supportRadius;

			if (ed->type == 0)
				emitterBoundary->meshFile = FileSystem::normalizePath(getExePath() + "/resources/emitter_boundary/EmitterBox.obj");
			else if (ed->type == 1)
				emitterBoundary->meshFile = FileSystem::normalizePath(getExePath() + "/resources/emitter_boundary/EmitterCylinder.obj");
			scene.boundaryModels.push_back(emitterBoundary);
			
			// reuse particles if they are outside of a bounding box
			bool emitterReuseParticles = false;
			Vector3r emitterBoxMin(-1.0, -1.0, -1.0);
			Vector3r emitterBoxMax(1.0, 1.0, 1.0);
			for (auto material : scene.materials)
			{
				if (material->id == model->getId())
				{
					emitterReuseParticles = material->emitterReuseParticles;
					emitterBoxMin = material->emitterBoxMin;
					emitterBoxMax = material->emitterBoxMax;
				}
			}

			if (emitterReuseParticles)
				model->getEmitterSystem()->enableReuseParticles(emitterBoxMin, emitterBoxMax);

			emitter->setEmitStartTime(ed->emitStartTime);
			emitter->setEmitEndTime(ed->emitEndTime);
		}
	}
}

void SimulatorBase::createAnimationFields()
{
	Simulation *sim = Simulation::getCurrent();
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();

	//////////////////////////////////////////////////////////////////////////
	// animation fields
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int i = 0; i < scene.animatedFields.size(); i++)
	{
		SceneLoader::AnimationFieldData *data = scene.animatedFields[i];

		sim->getAnimationFieldSystem()->addAnimationField(
				data->particleFieldName, 
				data->x, data->rotation, data->scale,
				data->expression, 
				data->shapeType);

		AnimationField *field = sim->getAnimationFieldSystem()->getAnimationFields().back();
		field->setStartTime(data->startTime);
		field->setEndTime(data->endTime);
	}
}


void SimulatorBase::createFluidBlocks(std::map<std::string, unsigned int> &fluidIDs, std::vector<std::vector<Vector3r>> &fluidParticles, std::vector<std::vector<Vector3r>> &fluidVelocities)
{
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
	for (unsigned int i = 0; i < scene.fluidBlocks.size(); i++)
	{
		const unsigned int fluidIndex = fluidIDs[scene.fluidBlocks[i]->id];
		const Real diam = static_cast<Real>(2.0)*scene.particleRadius;

		Real xshift = diam;
		Real yshift = diam;
		const Real eps = static_cast<Real>(1.0e-9);
		if (scene.fluidBlocks[i]->mode == 1)
			yshift = sqrt(static_cast<Real>(3.0)) * scene.particleRadius + eps;
		else if (scene.fluidBlocks[i]->mode == 2)
		{
			xshift = sqrt(static_cast<Real>(6.0)) * diam / static_cast<Real>(3.0) + eps;
			yshift = sqrt(static_cast<Real>(3.0)) * scene.particleRadius + eps;
		}

		Vector3r diff = scene.fluidBlocks[i]->box.m_maxX - scene.fluidBlocks[i]->box.m_minX;
		if (scene.fluidBlocks[i]->mode == 1)
		{
			diff[0] -= diam;
			diff[2] -= diam;
		}
		else if (scene.fluidBlocks[i]->mode == 2)
		{
			diff[0] -= xshift;
			diff[2] -= diam;
		}

		const int stepsX = (int)round(diff[0] / xshift) - 1;
		const int stepsY = (int)round(diff[1] / yshift) - 1;
		int stepsZ = (int)round(diff[2] / diam) - 1;

		Vector3r start = scene.fluidBlocks[i]->box.m_minX + static_cast<Real>(2.0)*scene.particleRadius*Vector3r::Ones();
		const unsigned int startIndex = (unsigned int)fluidParticles[fluidIndex].size();
		const unsigned int numAddedParticles = stepsX * stepsY * stepsZ;
		fluidParticles[fluidIndex].reserve(fluidParticles[fluidIndex].size() + numAddedParticles);
		fluidVelocities[fluidIndex].resize(fluidVelocities[fluidIndex].size() + numAddedParticles);

		if (Simulation::getCurrent()->is2DSimulation())
		{
			stepsZ = 1;
			start[2] = 0.0;
		}

		for (int j = 0; j < stepsX; j++)
		{
			for (int k = 0; k < stepsY; k++)
			{
				for (int l = 0; l < stepsZ; l++)
				{
					Vector3r currPos = Vector3r(j*xshift, k*yshift, l*diam) + start;
					if (scene.fluidBlocks[i]->mode == 1)
					{
						if (k % 2 == 0)
							currPos += Vector3r(0, 0, scene.particleRadius);
						else
							currPos += Vector3r(scene.particleRadius, 0, 0);
					}
					else if (scene.fluidBlocks[i]->mode == 2)
					{
						currPos += Vector3r(0, 0, scene.particleRadius);

						Vector3r shift_vec(0, 0, 0);
						if ((j % 2) && !Simulation::getCurrent()->is2DSimulation())
						{
							shift_vec[2] += diam / (static_cast<Real>(2.0) * (k % 2 ? -1 : 1));
						}
						if (k % 2 == 0)
						{
							shift_vec[0] += xshift / static_cast<Real>(2.0);
						}
						currPos += shift_vec;
					}
					fluidParticles[fluidIndex].push_back(currPos);
				}
			}
		}
		setInitialVelocity(scene.fluidBlocks[i]->initialVelocity, scene.fluidBlocks[i]->initialAngularVelocity, numAddedParticles, &fluidParticles[fluidIndex][startIndex], &fluidVelocities[fluidIndex][startIndex]);
	}
}

void SimulatorBase::particleInfo(std::vector<std::vector<unsigned int>> &particles)
{
	Simulation *sim = Simulation::getCurrent();
	const int maxWidth = 25;
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		if (particles[i].size() > 0)
		{
			LOG_INFO << "---------------------------------------------------------------------------";
			LOG_INFO << model->getId();
			LOG_INFO << "---------------------------------------------------------------------------";
		}
		for (unsigned int j = 0; j < particles[i].size(); j++)
		{
			unsigned int index = particles[i][j];
 			LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << "Index:" << index;
 			for (unsigned int k = 0; k < model->numberOfFields(); k++)
 			{
				const FieldDescription &field = model->getField(k);
				if (field.type == Scalar)
 					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << field.name + ":" << *((Real*) field.getFct(index));
				else if (field.type == UInt)
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << field.name + ":" << *((unsigned int*)field.getFct(index));
				else if (field.type == Vector3)
				{
					Eigen::Map<Vector3r> vec((Real*) field.getFct(index));
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << field.name + ":" << vec.transpose();
				}
				else if (field.type == Vector6)
				{
					Eigen::Map<Vector6r> vec((Real*)field.getFct(index));
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << field.name + ":" << vec.transpose();
				}
				else if (field.type == Matrix3)
				{
					Eigen::Map<Matrix3r> mat((Real*)field.getFct(index));
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << field.name + ":" << mat.row(0);
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << " " << mat.row(1);
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << " " << mat.row(2);
				}
				else if (field.type == Matrix6)
				{
					Eigen::Map<Matrix6r> mat((Real*) field.getFct(index));
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << field.name + ":" << mat.row(0);
					for (unsigned int k = 1; k < 6; k++)
						LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << " " << mat.row(k);
				}
 			}
			LOG_INFO << "---------------------------------------------------------------------------\n";
		}
	}
}

void SimulatorBase::step()
{
	if (TimeManager::getCurrent()->getTime() >= m_nextFrameTime)
	{
		m_nextFrameTime += static_cast<Real>(1.0) / m_framesPerSecond;

		for (size_t i = 0; i < m_particleExporters.size(); i++)
		{
			m_particleExporters[i].m_exporter->step(m_frameCounter);
		}
		for (size_t i = 0; i < m_rbExporters.size(); i++)
			m_rbExporters[i].m_exporter->step(m_frameCounter);

		m_frameCounter++;
	}
	if (TimeManager::getCurrent()->getTime() >= m_nextFrameTimeState)
	{
		m_nextFrameTimeState += static_cast<Real>(1.0) / m_framesPerSecondState;
		if (m_enableStateExport)
			saveState();
	}
#ifdef DL_OUTPUT
	if (TimeManager::getCurrent()->getTime() >= m_nextTiming)
	{
		LOG_INFO << "---------------------------------------------------------------------------";
		LOG_INFO << "Time: " << TimeManager::getCurrent()->getTime();
		Timing::printAverageTimes();
		Timing::printTimeSums();
		Counting::printAverageCounts();
		Counting::printCounterSums();
		m_nextTiming += 1.0;
	}
#endif
}

void SimulatorBase::updateBoundaryParticles(const bool forceUpdate = false)
{
	Simulation *sim = Simulation::getCurrent();
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();

	const unsigned int nObjects = sim->numberOfBoundaryModels();
	for (unsigned int i = 0; i < nObjects; i++)
	{
		BoundaryModel_Akinci2012 *bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(i));
		RigidBodyObject *rbo = bm->getRigidBodyObject();
		if (rbo->isDynamic() || forceUpdate)
		{
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static)  
				for (int j = 0; j < (int)bm->numberOfParticles(); j++)
				{
					bm->getPosition(j) = rbo->getRotation() * bm->getPosition0(j) + rbo->getPosition();
					if (rbo->isDynamic())
						bm->getVelocity(j) = rbo->getAngularVelocity().cross(bm->getPosition(j) - rbo->getPosition()) + rbo->getVelocity();
					else
						bm->getVelocity(j).setZero();
				}
			}
			#ifdef GPU_NEIGHBORHOOD_SEARCH
			// copy the particle data to the GPU
			if (forceUpdate)
				sim->getNeighborhoodSearch()->update_point_sets();
			#endif 
		}
	}
}

void SPH::SimulatorBase::updateDMVelocity()
{
	Simulation *sim = Simulation::getCurrent();
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();

	const unsigned int nObjects = sim->numberOfBoundaryModels();
	for (unsigned int i = 0; i < nObjects; i++)
	{
		BoundaryModel_Koschier2017 *bm = static_cast<BoundaryModel_Koschier2017*>(sim->getBoundaryModel(i));
		RigidBodyObject *rbo = bm->getRigidBodyObject();
		if (rbo->isDynamic())
		{
			const Real maxDist = bm->getMaxDist();
			const Vector3r x(maxDist, 0.0, 0.0);
			const Vector3r vel = rbo->getAngularVelocity().cross(x) + rbo->getVelocity();
			bm->setMaxVel(vel.norm());
		}
	}
}

void SPH::SimulatorBase::updateVMVelocity()
{
	Simulation *sim = Simulation::getCurrent();
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
	const unsigned int nObjects = sim->numberOfBoundaryModels();
	for (unsigned int i = 0; i < nObjects; i++)
	{
		BoundaryModel_Bender2019 *bm = static_cast<BoundaryModel_Bender2019*>(sim->getBoundaryModel(i));
		RigidBodyObject *rbo = bm->getRigidBodyObject();
		if (rbo->isDynamic())
		{
			const Real maxDist = bm->getMaxDist();
			const Vector3r x(maxDist, 0.0, 0.0);
			const Vector3r vel = rbo->getAngularVelocity().cross(x) + rbo->getVelocity();
			bm->setMaxVel(vel.norm());
		}
	}
}

void SPH::SimulatorBase::saveState(const std::string& stateFile)
{
	std::string stateFilePath; 
	std::string exportFileName;
	const Real time = TimeManager::getCurrent()->getTime();
	const std::string timeStr = StringTools::real2String(time);
	if (stateFile == "")
	{
		stateFilePath = FileSystem::normalizePath(m_outputPath + "/state");
		exportFileName = FileSystem::normalizePath(stateFilePath + "/state_" + timeStr);
	}
	else
	{
		stateFilePath = FileSystem::getFilePath(stateFile);
		exportFileName = FileSystem::normalizePath(stateFilePath + "/" + FileSystem::getFileName(stateFile));
	}
	FileSystem::makeDirs(stateFilePath);

	const std::string& sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
	string md5Str = FileSystem::getFileMD5(sceneFile);

	Simulation *sim = Simulation::getCurrent();

		// Save additional data
	BinaryFileWriter binWriter;
	binWriter.openFile(exportFileName + ".bin");
	binWriter.write(md5Str);

	binWriter.write(m_nextFrameTime);
	binWriter.write(m_nextFrameTimeState);
	binWriter.write(m_frameCounter);
	binWriter.write(m_isFirstFrame);
	binWriter.write(m_isFirstFrameVTK);
	
	writeParameterState(binWriter);
	TimeManager::getCurrent()->saveState(binWriter);
	Simulation::getCurrent()->saveState(binWriter);

	// fluid models
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		std::string fileName = "particle";
		fileName = fileName + "_" + model->getId(); // +"_" + std::to_string(m_frameCounter);

		// Save particle data
		std::string expFileName = FileSystem::normalizePath(exportFileName + "_" + fileName);
		writeFluidParticlesState(expFileName + ".bgeo", model);
	}

	// boundary models
	if (m_firstState)
	{
		for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
		{
			BoundaryModel *model = sim->getBoundaryModel(i);
			std::string fileName = "boundary";
			fileName = fileName + "_" + to_string(i); // +"_" + std::to_string(m_frameCounter);

			// Save particle data
			std::string expFileName = FileSystem::normalizePath(exportFileName + "_" + fileName);
			writeBoundaryState(expFileName + ".bgeo", model);
		}
		m_firstState = false;
	}

	// dynamic bodies
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel *bm = sim->getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic())
		{
			binWriter.writeMatrix(bm->getRigidBodyObject()->getPosition());
			binWriter.writeMatrix(bm->getRigidBodyObject()->getRotation());
			binWriter.writeMatrix(bm->getRigidBodyObject()->getVelocity());
			binWriter.writeMatrix(bm->getRigidBodyObject()->getAngularVelocity());
		}
	}
	binWriter.closeFile();

	LOG_INFO << "Saved state: " << exportFileName + ".bin";
}

#ifdef WIN32
void SPH::SimulatorBase::loadStateDialog()
{
	const std::string stateFilePath = FileSystem::normalizePath(m_outputPath + "/state");
	const std::string stateFileName = FileSystem::fileDialog(0, stateFilePath, "*.bin");
	if (stateFileName == "")
		return;
	loadState(stateFileName);
}
#endif 

void SPH::SimulatorBase::loadState(const std::string &stateFile)
{
	Simulation *sim = Simulation::getCurrent();
	const std::string& sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
	string md5Str = FileSystem::getFileMD5(sceneFile);

	// Load additional data
	BinaryFileReader binReader;
	if (!binReader.openFile(stateFile))
		return;

	// check if scene file has changed
	std::string md5StrState;
	binReader.read(md5StrState);

	binReader.read(m_nextFrameTime);
	binReader.read(m_nextFrameTimeState);
	binReader.read(m_frameCounter);
	binReader.read(m_isFirstFrame);
	binReader.read(m_isFirstFrameVTK);

	readParameterState(binReader);
	if (md5Str != md5StrState)
		LOG_WARN << "State was stored for another scene file.";
	TimeManager::getCurrent()->loadState(binReader);
	Simulation::getCurrent()->loadState(binReader);

	const std::string importFilePath = FileSystem::getFilePath(stateFile);
	const std::string importFileName = FileSystem::getFileName(stateFile);

	// fluid models
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		std::string fileName = "particle";
		fileName = fileName + "_" + model->getId(); // +"_" + std::to_string(m_frameCounter);

		std::string impFileName = FileSystem::normalizePath(importFilePath + "/" + importFileName + "_" + fileName);
		readFluidParticlesState(impFileName + ".bgeo", model);
	}

	// boundary models
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel *model = sim->getBoundaryModel(i);
		std::string fileName = "boundary";
		fileName = fileName + "_" + to_string(i); // +"_" + std::to_string(m_frameCounter);

		// Save particle data
		std::string impFileName = FileSystem::normalizePath(importFilePath + "/" + importFileName + "_" + fileName);
		readBoundaryState(impFileName + ".bgeo", model);
	}

	// dynamic bodies
	bool dynamic = false;
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel *bm = sim->getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic())
		{
			Vector3r x;
			binReader.readMatrix(x);
			bm->getRigidBodyObject()->setPosition(x);

			Matrix3r R;
			binReader.readMatrix(R);
			bm->getRigidBodyObject()->setRotation(R);

			Vector3r v;
			binReader.readMatrix(v);
			bm->getRigidBodyObject()->setVelocity(v);

			binReader.readMatrix(v);
			bm->getRigidBodyObject()->setAngularVelocity(v);

			dynamic = true;
		}		
	}
	if (dynamic)
	{
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			updateBoundaryParticles(true);
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			updateDMVelocity();
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			updateVMVelocity();
	}
	binReader.closeFile();

	//sim->performNeighborhoodSearchSort();
}

void SimulatorBase::writeParameterState(BinaryFileWriter &binWriter)
{
	writeParameterObjectState(binWriter, this);
	writeParameterObjectState(binWriter, Simulation::getCurrent());
	writeParameterObjectState(binWriter, Simulation::getCurrent()->getTimeStep());
#ifdef USE_DEBUG_TOOLS
	writeParameterObjectState(binWriter, Simulation::getCurrent()->getDebugTools());
#endif

	Simulation *sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		writeParameterObjectState(binWriter, model);
		writeParameterObjectState(binWriter, (ParameterObject*)model->getDragBase());
		writeParameterObjectState(binWriter, (ParameterObject*)model->getSurfaceTensionBase());
		writeParameterObjectState(binWriter, (ParameterObject*)model->getViscosityBase());
		writeParameterObjectState(binWriter, (ParameterObject*)model->getVorticityBase());
		writeParameterObjectState(binWriter, (ParameterObject*)model->getElasticityBase());

		binWriter.write(getColorField(model->getPointSetIndex()));
		binWriter.write(getColorMapType(model->getPointSetIndex()));
		binWriter.write(getRenderMinValue(model->getPointSetIndex()));
		binWriter.write(getRenderMaxValue(model->getPointSetIndex()));
	}
}

void SimulatorBase::writeParameterObjectState(BinaryFileWriter &binWriter, GenParam::ParameterObject *paramObj)
{
	if (paramObj == nullptr)
		return;

	const unsigned int numParams = paramObj->numParameters();

	for (unsigned int i = 0; i < numParams; i++)
	{
		ParameterBase *paramBase = paramObj->getParameter(i);

		if (paramBase->getType() == RealParameterType)
			binWriter.write(static_cast<NumericParameter<Real>*>(paramBase)->getValue());
 		else if (paramBase->getType() == ParameterBase::UINT32)
 			binWriter.write(static_cast<NumericParameter<unsigned int>*>(paramBase)->getValue());
 		else if (paramBase->getType() == ParameterBase::UINT16)
 			binWriter.write(static_cast<NumericParameter<unsigned short>*>(paramBase)->getValue());
 		else if (paramBase->getType() == ParameterBase::UINT8)
 			binWriter.write(static_cast<NumericParameter<unsigned char>*>(paramBase)->getValue());
 		else if (paramBase->getType() == ParameterBase::INT32)
 			binWriter.write(static_cast<NumericParameter<int>*>(paramBase)->getValue());
 		else if (paramBase->getType() == ParameterBase::INT16)
 			binWriter.write(static_cast<NumericParameter<short>*>(paramBase)->getValue());
 		else if (paramBase->getType() == ParameterBase::INT8)
 			binWriter.write(static_cast<NumericParameter<char>*>(paramBase)->getValue());
 		else if (paramBase->getType() == ParameterBase::ENUM)
 			binWriter.write(static_cast<EnumParameter*>(paramBase)->getValue());
 		else if (paramBase->getType() == ParameterBase::BOOL)
 			binWriter.write(static_cast<BoolParameter*>(paramBase)->getValue());
 		else if (paramBase->getType() == RealVectorParameterType)
 		{
 			VectorParameter<Real> *vec = static_cast<VectorParameter<Real>*>(paramBase);
 			binWriter.writeBuffer((char*)vec->getValue(), vec->getDim() * sizeof(Real));;
 		}
 		else if (paramBase->getType() == ParameterBase::STRING)
 			binWriter.write(static_cast<StringParameter*>(paramBase)->getValue());
	}
}

void SimulatorBase::readParameterState(BinaryFileReader &binReader)
{
	readParameterObjectState(binReader, this);
 	readParameterObjectState(binReader, Simulation::getCurrent());
 	readParameterObjectState(binReader, Simulation::getCurrent()->getTimeStep());
#ifdef USE_DEBUG_TOOLS
	readParameterObjectState(binReader, Simulation::getCurrent()->getDebugTools());
#endif

 
 	Simulation *sim = Simulation::getCurrent();
 	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
 	{
 		FluidModel *model = sim->getFluidModel(i);
 		readParameterObjectState(binReader, model);
 		readParameterObjectState(binReader, (ParameterObject*)model->getDragBase());
 		readParameterObjectState(binReader, (ParameterObject*)model->getSurfaceTensionBase());
 		readParameterObjectState(binReader, (ParameterObject*)model->getViscosityBase());
 		readParameterObjectState(binReader, (ParameterObject*)model->getVorticityBase());
 		readParameterObjectState(binReader, (ParameterObject*)model->getElasticityBase());
 
		std::string field;
		binReader.read(field);
		setColorField(model->getPointSetIndex(), field);
		int type;
		binReader.read(type);
		setColorMapType(model->getPointSetIndex(), type);
		Real v;
		binReader.read(v);
		setRenderMinValue(model->getPointSetIndex(), v);
		binReader.read(v);
		setRenderMaxValue(model->getPointSetIndex(), v);
 	}
}

void SimulatorBase::readParameterObjectState(BinaryFileReader &binReader, GenParam::ParameterObject *paramObj)
{
	if (paramObj == nullptr)
		return;

	const unsigned int numParams = paramObj->numParameters();

	for (unsigned int i = 0; i < numParams; i++)
	{
		ParameterBase *paramBase = paramObj->getParameter(i);

		if (paramBase->getType() == RealParameterType)
		{
			Real val;
			binReader.read(val);
			static_cast<NumericParameter<Real>*>(paramBase)->setValue(val);
		}
		else if (paramBase->getType() == ParameterBase::UINT32)
		{
			unsigned int val;
			binReader.read(val);
			static_cast<NumericParameter<unsigned int>*>(paramBase)->setValue(val);
		}
		else if (paramBase->getType() == ParameterBase::UINT16)
		{
			unsigned short val;
			binReader.read(val);
			static_cast<NumericParameter<unsigned short>*>(paramBase)->setValue(val);
		}
		else if (paramBase->getType() == ParameterBase::UINT8)
		{
			unsigned char val;
			binReader.read(val);
			static_cast<NumericParameter<unsigned char>*>(paramBase)->setValue(val);
		}
		else if (paramBase->getType() == ParameterBase::INT32)
		{
			int val;
			binReader.read(val);
			static_cast<NumericParameter<int>*>(paramBase)->setValue(val);
		}
		else if (paramBase->getType() == ParameterBase::INT16)
		{
			short val;
			binReader.read(val);
			static_cast<NumericParameter<short>*>(paramBase)->setValue(val);
		}
		else if (paramBase->getType() == ParameterBase::INT8)
		{
			char val;
			binReader.read(val);
			static_cast<NumericParameter<char>*>(paramBase)->setValue(val);
		}
		else if (paramBase->getType() == ParameterBase::ENUM)
		{
			int val;
			binReader.read(val);
			static_cast<EnumParameter*>(paramBase)->setValue(val);
		}
 		else if (paramBase->getType() == ParameterBase::BOOL)
		{
			bool val;
			binReader.read(val);
			static_cast<BoolParameter*>(paramBase)->setValue(val);
		}
 		else if (paramBase->getType() == RealVectorParameterType)
 		{
 			VectorParameter<Real> *vec = static_cast<VectorParameter<Real>*>(paramBase);
 			binReader.readBuffer((char*)vec->getValue(), vec->getDim() * sizeof(Real));;
 		}
 		else if (paramBase->getType() == ParameterBase::STRING)
		{
			std::string val;
			binReader.read(val);
			static_cast<StringParameter*>(paramBase)->setValue(val);
		}
	}
}


void SimulatorBase::writeFluidParticlesState(const std::string &fileName, FluidModel *model)
{
	Partio::ParticlesDataMutable& particleData = *Partio::create();

	std::vector<std::pair<unsigned int, Partio::ParticleAttribute>> partioAttr;
	for (unsigned int j = 0; j < model->numberOfFields(); j++)
	{
		const FieldDescription &field = model->getField(j);
		if (field.storeData)
		{
//			LOG_INFO << "Store field: " << field.name;
			if (field.type == Scalar)
			{
				partioAttr.push_back({ j, particleData.addAttribute(field.name.c_str(), Partio::FLOAT, 1) });
			}
			else if (field.type == UInt)
			{
				partioAttr.push_back({ j, particleData.addAttribute(field.name.c_str(), Partio::INT, 1) });
			}
			else if (field.type == Vector3)
			{
				partioAttr.push_back({ j, particleData.addAttribute(field.name.c_str(), Partio::VECTOR, 3) });
			}
			else
			{
				LOG_WARN << "Only scalar and vector fields are currently supported by the partio exporter.";
			}
		}
	}

	const unsigned int numParticles = model->numActiveParticles();

	for (unsigned int i = 0; i < numParticles; i++)
	{
		Partio::ParticleIndex index = particleData.addParticle();

		for (unsigned int j = 0; j < partioAttr.size(); j++)
		{
			const unsigned int fieldIndex = partioAttr[j].first;
			const Partio::ParticleAttribute &attr = partioAttr[j].second;

			const FieldDescription &field = model->getField(fieldIndex);
			if (field.type == FieldType::Scalar)
			{
				float* val = particleData.dataWrite<float>(attr, index);
				*val = (float)* ((Real*) field.getFct(i));
			}
			else if (field.type == FieldType::UInt)
			{
				int* val = particleData.dataWrite<int>(attr, index);
				*val = (int)* ((unsigned int*)field.getFct(i));
			}
			else if (field.type == FieldType::Vector3)
			{
				float* val = particleData.dataWrite<float>(attr, index);
				Eigen::Map<Vector3r> vec((Real*) field.getFct(i));
				val[0] = (float)vec[0];
				val[1] = (float)vec[1];
				val[2] = (float)vec[2];
			}
		}
	}

	Partio::write(fileName.c_str(), particleData, true);
	particleData.release();
}


void SimulatorBase::readFluidParticlesState(const std::string &fileName, FluidModel *model)
{
	if (!FileSystem::fileExists(fileName))
	{
		LOG_WARN << "File " << fileName << " does not exist.";
		return;
	}

	Partio::ParticlesDataMutable* data = Partio::read(fileName.c_str());
	if (!data)
	{
		LOG_WARN << "Partio file " << fileName << " not readable.";
		return;
	}

	std::vector<std::pair<unsigned int, Partio::ParticleAttribute>> partioAttr;
	for (int i = 0; i < data->numAttributes(); i++)
	{
		Partio::ParticleAttribute attr;
		data->attributeInfo(i, attr);
		for (unsigned int j = 0; j < model->numberOfFields(); j++)
		{
			const FieldDescription &field = model->getField(j);
			if (field.name == attr.name)
			{
				//LOG_INFO << "Read field: " << field.name;
				partioAttr.push_back({ j, attr});
			}
		}
	}

	for (unsigned int j = 0; j < partioAttr.size(); j++)
	{
		const unsigned int fieldIndex = partioAttr[j].first;
		const Partio::ParticleAttribute &attr = partioAttr[j].second;
 
		const FieldDescription &field = model->getField(fieldIndex);

		for (int i = 0; i < data->numParticles(); i++)
		{			
			if (field.type == FieldType::Scalar)
			{
				const float *value = data->data<float>(attr, i);
				*((Real*) field.getFct(i)) = value[0];
			}
			else if (field.type == FieldType::UInt)
			{
				const int *value = data->data<int>(attr, i);
				*((unsigned int*)field.getFct(i)) = value[0];
			}
			else if (field.type == FieldType::Vector3)
			{
				const float *value = data->data<float>(attr, i);
				Eigen::Map<Vector3r> vec((Real*) field.getFct(i));
				vec[0] = value[0];
				vec[1] = value[1];
				vec[2] = value[2];
			}
		}
	}
	data->release();
}


void SimulatorBase::writeBoundaryState(const std::string &fileName, BoundaryModel *bm)
{
	Simulation *sim = Simulation::getCurrent();
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		BoundaryModel_Akinci2012 *model = static_cast<BoundaryModel_Akinci2012*>(bm);
		Partio::ParticlesDataMutable& particleData = *Partio::create();
		const Partio::ParticleAttribute &attrX0 = particleData.addAttribute("position0", Partio::VECTOR, 3);
		const Partio::ParticleAttribute &attrX = particleData.addAttribute("position", Partio::VECTOR, 3);
		const Partio::ParticleAttribute &attrVel = particleData.addAttribute("velocity", Partio::VECTOR, 3);
		const Partio::ParticleAttribute &attrVol = particleData.addAttribute("volume", Partio::FLOAT, 1);

		const unsigned int numParticles = model->numberOfParticles();

		for (unsigned int i = 0; i < numParticles; i++)
		{
			Partio::ParticleIndex index = particleData.addParticle();

			float* val = particleData.dataWrite<float>(attrX0, index);
			const Vector3r &x0 = model->getPosition0(i);
			val[0] = (float)x0[0];
			val[1] = (float)x0[1];
			val[2] = (float)x0[2];

			val = particleData.dataWrite<float>(attrX, index);
			const Vector3r &x = model->getPosition(i);
			val[0] = (float)x[0];
			val[1] = (float)x[1];
			val[2] = (float)x[2];

			val = particleData.dataWrite<float>(attrVel, index);
			const Vector3r &v = model->getVelocity(i);
			val[0] = (float)v[0];
			val[1] = (float)v[1];
			val[2] = (float)v[2];

			val = particleData.dataWrite<float>(attrVol, index);
			val[0] = (float)model->getVolume(i);
		}

		Partio::write(fileName.c_str(), particleData, true);
		particleData.release();
	}
}


void SimulatorBase::readBoundaryState(const std::string &fileName, BoundaryModel *bm)
{
	Simulation *sim = Simulation::getCurrent();
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		if (!FileSystem::fileExists(fileName))
		{
			LOG_WARN << "File " << fileName << " does not exist.";
			return;
		}

		BoundaryModel_Akinci2012 *model = static_cast<BoundaryModel_Akinci2012*>(bm);
		Partio::ParticlesDataMutable* data = Partio::read(fileName.c_str());
		if (!data)
		{
			LOG_WARN << "Partio file " << fileName << " not readable.";
			return;
		}

		unsigned int pos0Index = 0xffffffff;
		unsigned int posIndex = 0xffffffff;
		unsigned int velIndex = 0xffffffff;
		unsigned int volIndex = 0xffffffff;

		for (int i = 0; i < data->numAttributes(); i++)
		{
			Partio::ParticleAttribute attr;
			data->attributeInfo(i, attr);
			if (attr.name == "position0")
				pos0Index = i;
			else if (attr.name == "position")
				posIndex = i;
			else if (attr.name == "velocity")
				velIndex = i;
			else if (attr.name == "volume")
				volIndex = i;
		}

		if ((pos0Index == 0xffffffff) ||
			(posIndex == 0xffffffff) ||
			(velIndex == 0xffffffff) ||
			(volIndex == 0xffffffff))
		{
			LOG_WARN << "File " << fileName << " does not has the correct attributes.";
			return;
		}

		Partio::ParticleAttribute attrX0;
		Partio::ParticleAttribute attrX;
		Partio::ParticleAttribute attrVel;
		Partio::ParticleAttribute attrVol;

		data->attributeInfo(pos0Index, attrX0);
		data->attributeInfo(posIndex, attrX);
		data->attributeInfo(velIndex, attrVel);
		data->attributeInfo(volIndex, attrVol);

		model->resize(data->numParticles());
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *pos0 = data->data<float>(attrX0, i);
			model->getPosition0(i) = Vector3r(pos0[0], pos0[1], pos0[2]);

			const float *pos = data->data<float>(attrX, i);
			model->getPosition(i) = Vector3r(pos[0], pos[1], pos[2]);

			const float *vel = data->data<float>(attrVel, i);
			model->getVelocity(i) = Vector3r(vel[0], vel[1], vel[2]);

			const float *vol = data->data<float>(attrVol, i);
			model->setVolume(i, vol[0]);
		}

		data->release();

		NeighborhoodSearch *neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();
		neighborhoodSearch->update_point_sets();
		neighborhoodSearch->resize_point_set(model->getPointSetIndex(), &model->getPosition(0)[0], model->numberOfParticles());
	}
}


void SimulatorBase::initDensityMap(std::vector<Vector3r> &x, std::vector<unsigned int> &faces, const Utilities::SceneLoader::BoundaryData *boundaryData, const bool md5, const bool isDynamic, BoundaryModel_Koschier2017 *boundaryModel)
{
	Simulation *sim = Simulation::getCurrent();
	const std::string& sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
	const Real supportRadius = sim->getSupportRadius();
	std::string scene_path = FileSystem::getFilePath(sceneFile);
	std::string scene_file_name = FileSystem::getFileName(sceneFile);
	const bool useCache = getUseParticleCaching();
	Discregrid::CubicLagrangeDiscreteGrid *densityMap;

	// if a map file is given, use this one
	if (boundaryData->mapFile != "")
	{
		std::string mapFileName = boundaryData->mapFile;
		if (FileSystem::isRelativePath(mapFileName))
			mapFileName = FileSystem::normalizePath(scene_path + "/" + mapFileName);
		densityMap = new Discregrid::CubicLagrangeDiscreteGrid(mapFileName);
		boundaryModel->setMap(densityMap);
		LOG_INFO << "Loaded density map: " << mapFileName;
		return;
	}

	string cachePath = scene_path + "/Cache";

	// Cache map
	std::string mesh_base_path = FileSystem::getFilePath(boundaryData->meshFile);
	std::string mesh_file_name = FileSystem::getFileName(boundaryData->meshFile);


	Eigen::Matrix<unsigned int, 3, 1> resolutionSDF = boundaryData->mapResolution;
	const string scaleStr = "s" + StringTools::real2String(boundaryData->scale[0]) + "_" + StringTools::StringTools::real2String(boundaryData->scale[1]) + "_" + StringTools::real2String(boundaryData->scale[2]);
	const string resStr = "r" + to_string(resolutionSDF[0]) + "_" + to_string(resolutionSDF[1]) + "_" + to_string(resolutionSDF[2]);
	const string invertStr = "i" + to_string((int)boundaryData->mapInvert);
	const string thicknessStr = "t" + StringTools::real2String(boundaryData->mapThickness);
	const string kernelStr = "k" + to_string(sim->getKernel());
	string densityMapFileName = "";
	if (isDynamic)
		densityMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_db_dm_" + StringTools::real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + "_" + kernelStr + ".cdm");
	else 
		densityMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_sb_dm_" + StringTools::real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + "_" + kernelStr + ".cdm");

	// check MD5 if cache file is available
	bool foundCacheFile = false;

	if (useCache)
		foundCacheFile = FileSystem::fileExists(densityMapFileName);

	if (useCache && foundCacheFile && md5)
	{
		densityMap = new Discregrid::CubicLagrangeDiscreteGrid(densityMapFileName);
		boundaryModel->setMap(densityMap);
		LOG_INFO << "Loaded cached density map: " << densityMapFileName;
		return;
	}

	if (!useCache || !foundCacheFile || !md5)
	{
		//////////////////////////////////////////////////////////////////////////
		// Generate distance field of object using Discregrid
		//////////////////////////////////////////////////////////////////////////
#ifdef USE_DOUBLE
		Discregrid::TriangleMesh sdfMesh(&x[0][0], faces.data(), x.size(), faces.size() / 3);
#else
		// if type is float, copy vector to double vector
		std::vector<double> doubleVec;
		doubleVec.resize(3 * x.size());
		for (unsigned int i = 0; i < x.size(); i++)
			for (unsigned int j = 0; j < 3; j++)
				doubleVec[3 * i + j] = x[i][j];
		Discregrid::TriangleMesh sdfMesh(&doubleVec[0], faces.data(), x.size(), faces.size() / 3);
#endif

		Discregrid::MeshDistance md(sdfMesh);
		Eigen::AlignedBox3d domain;
		for (auto const& x_ : x)
		{
			domain.extend(x_.cast<double>());
		}
		const Real tolerance = boundaryData->mapThickness;
		domain.max() += (4.0*supportRadius + tolerance) * Eigen::Vector3d::Ones();
		domain.min() -= (4.0*supportRadius + tolerance) * Eigen::Vector3d::Ones();

		LOG_INFO << "Domain - min: " << domain.min()[0] << ", " << domain.min()[1] << ", " << domain.min()[2];
		LOG_INFO << "Domain - max: " << domain.max()[0] << ", " << domain.max()[1] << ", " << domain.max()[2];

		LOG_INFO << "Set SDF resolution: " << resolutionSDF[0] << ", " << resolutionSDF[1] << ", " << resolutionSDF[2];
		densityMap = new Discregrid::CubicLagrangeDiscreteGrid(domain, std::array<unsigned int, 3>({ resolutionSDF[0], resolutionSDF[1], resolutionSDF[2] }));
		auto func = Discregrid::DiscreteGrid::ContinuousFunction{};

		Real sign = 1.0;
		if (boundaryData->mapInvert)
			sign = -1.0;
		func = [&md, &sign, &tolerance](Eigen::Vector3d const& xi) {return sign * (md.signedDistanceCached(xi) - tolerance); };

		LOG_INFO << "Generate SDF";
		START_TIMING("SDF Construction");
		densityMap->addFunction(func, false);
		STOP_TIMING_AVG

		const bool sim2D = sim->is2DSimulation();

		//////////////////////////////////////////////////////////////////////////
		// Generate density map of object using Discregrid
		//////////////////////////////////////////////////////////////////////////
		if (sim2D)
			SimpleQuadrature::determineSamplePointsInCircle(supportRadius, 30);

		auto int_domain = Eigen::AlignedBox3d(Eigen::Vector3d::Constant(-supportRadius), Eigen::Vector3d::Constant(supportRadius));
		Real factor = 5.0;
		if (sim2D)
			factor = 1.75;
		auto density_func = [&](Eigen::Vector3d const& x)
		{
			auto d = densityMap->interpolate(0u, x);
			if (d > (1.0 + 1.0 / factor) * supportRadius)
			{
				return 0.0;
			}

			auto integrand = [&](Eigen::Vector3d const& xi)
			{
				if (xi.squaredNorm() > supportRadius*supportRadius)
					return 0.0;

				auto dist = densityMap->interpolate(0u, x + xi);

				// Linear function gamma
				if (dist > 1.0 / factor * supportRadius)
					return 0.0;
				return static_cast<double>((1.0 - factor * dist / supportRadius) * sim->W(xi.cast<Real>()));
			};

			double res = 0.0;
			if (sim2D)
				res = 0.8 * SimpleQuadrature::integrate(integrand);
			else
				res = 0.8 * GaussQuadrature::integrate(integrand, int_domain, 50);
			
			return res;
		};

		auto cell_diag = densityMap->cellSize().norm();
		std::cout << "Generate density map..." << std::endl;
		const bool no_reduction = true;
		START_TIMING("Density Map Construction");
		densityMap->addFunction(density_func, false, [&](Eigen::Vector3d const& x_)
		{
			if (no_reduction)
			{
				return true;
			}
			auto x = x_.cwiseMax(densityMap->domain().min()).cwiseMin(densityMap->domain().max());
			auto dist = densityMap->interpolate(0u, x);
			if (dist == std::numeric_limits<double>::max())
			{
				return false;
			}

			return fabs(dist) < 2.5 * supportRadius;
		});
		STOP_TIMING_PRINT;

		// reduction
		if (!no_reduction)
		{
			std::cout << "Reduce discrete fields...";
			densityMap->reduceField(0u, [&](const Eigen::Vector3d &, double v)
			{
				return fabs(v) < 2.5 * supportRadius;
			});
			densityMap->reduceField(1u, [&](const Eigen::Vector3d &, double v)->double
			{
				if (v == std::numeric_limits<double>::max())
					return false;
				return true;
			});
			std::cout << "DONE" << std::endl;
		}

		boundaryModel->setMap(densityMap);

		// Store cache file
		if (useCache && (FileSystem::makeDir(cachePath) == 0))
		{
			LOG_INFO << "Save density map: " << densityMapFileName;
			densityMap->save(densityMapFileName);
		}
	}
}

void SimulatorBase::initVolumeMap(std::vector<Vector3r> &x, std::vector<unsigned int> &faces, const Utilities::SceneLoader::BoundaryData *boundaryData, const bool md5, const bool isDynamic, BoundaryModel_Bender2019 *boundaryModel)
{
	Simulation *sim = Simulation::getCurrent();
	const std::string& sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
	const Real supportRadius = sim->getSupportRadius();
	std::string scene_path = FileSystem::getFilePath(sceneFile);
	std::string scene_file_name = FileSystem::getFileName(sceneFile);
	const bool useCache = getUseParticleCaching();
	Discregrid::CubicLagrangeDiscreteGrid *volumeMap;

	// if a map file is given, use this one
	if (boundaryData->mapFile != "")
	{
		std::string mapFileName = boundaryData->mapFile;
		if (FileSystem::isRelativePath(mapFileName))
			mapFileName = FileSystem::normalizePath(scene_path + "/" + mapFileName);
		volumeMap = new Discregrid::CubicLagrangeDiscreteGrid(mapFileName);
		boundaryModel->setMap(volumeMap);
		LOG_INFO << "Loaded volume map: " << mapFileName;
		return;
	}

	string cachePath = scene_path + "/Cache";

	// Cache map
	std::string mesh_base_path = FileSystem::getFilePath(boundaryData->meshFile);
	std::string mesh_file_name = FileSystem::getFileName(boundaryData->meshFile);

	Eigen::Matrix<unsigned int, 3, 1> resolutionSDF = boundaryData->mapResolution;
	const string scaleStr = "s" + StringTools::real2String(boundaryData->scale[0]) + "_" + StringTools::real2String(boundaryData->scale[1]) + "_" + StringTools::real2String(boundaryData->scale[2]);
	const string resStr = "r" + to_string(resolutionSDF[0]) + "_" + to_string(resolutionSDF[1]) + "_" + to_string(resolutionSDF[2]);
	const string invertStr = "i" + to_string((int)boundaryData->mapInvert);
	const string thicknessStr = "t" + StringTools::real2String(boundaryData->mapThickness);
	string volumeMapFileName = "";
	if (isDynamic)
		volumeMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_db_vm_" + StringTools::real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + ".cdm");
	else 
		volumeMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_sb_vm_" + StringTools::real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + ".cdm");

	// check MD5 if cache file is available
	bool foundCacheFile = false;

	if (useCache)
		foundCacheFile = FileSystem::fileExists(volumeMapFileName);

	if (useCache && foundCacheFile && md5)
	{
		volumeMap = new Discregrid::CubicLagrangeDiscreteGrid(volumeMapFileName);
		boundaryModel->setMap(volumeMap);
		LOG_INFO << "Loaded cached volume map: " << volumeMapFileName;
		return;
	}

	if (!useCache || !foundCacheFile || !md5)
	{
		//////////////////////////////////////////////////////////////////////////
		// Generate distance field of object using Discregrid
		//////////////////////////////////////////////////////////////////////////
#ifdef USE_DOUBLE
		Discregrid::TriangleMesh sdfMesh(&x[0][0], faces.data(), x.size(), faces.size() / 3);
#else
		// if type is float, copy vector to double vector
		std::vector<double> doubleVec;
		doubleVec.resize(3 * x.size());
		for (unsigned int i = 0; i < x.size(); i++)
			for (unsigned int j = 0; j < 3; j++)
				doubleVec[3 * i + j] = x[i][j];
		Discregrid::TriangleMesh sdfMesh(&doubleVec[0], faces.data(), x.size(), faces.size() / 3);
#endif

		Discregrid::MeshDistance md(sdfMesh);
		Eigen::AlignedBox3d domain;
		for (auto const& x_ : x)
		{
			domain.extend(x_.cast<double>());
		}
		const Real tolerance = boundaryData->mapThickness;
		domain.max() += (4.0*supportRadius + tolerance) * Eigen::Vector3d::Ones();
		domain.min() -= (4.0*supportRadius + tolerance) * Eigen::Vector3d::Ones();

		LOG_INFO << "Domain - min: " << domain.min()[0] << ", " << domain.min()[1] << ", " << domain.min()[2];
		LOG_INFO << "Domain - max: " << domain.max()[0] << ", " << domain.max()[1] << ", " << domain.max()[2];

		LOG_INFO << "Set SDF resolution: " << resolutionSDF[0] << ", " << resolutionSDF[1] << ", " << resolutionSDF[2];
		volumeMap = new Discregrid::CubicLagrangeDiscreteGrid(domain, std::array<unsigned int, 3>({ resolutionSDF[0], resolutionSDF[1], resolutionSDF[2] }));
		auto func = Discregrid::DiscreteGrid::ContinuousFunction{};

		//volumeMap->setErrorTolerance(0.001);

		Real sign = 1.0;
		if (boundaryData->mapInvert)
			sign = -1.0;
		const Real particleRadius = sim->getParticleRadius();
		// subtract 0.5 * particle radius to prevent penetration of particles and the boundary
		func = [&md, &sign, &tolerance, &particleRadius](Eigen::Vector3d const& xi) {return sign * (md.signedDistanceCached(xi) - tolerance - 0.5*particleRadius); };

		LOG_INFO << "Generate SDF";
		START_TIMING("SDF Construction");
		volumeMap->addFunction(func, false);
		STOP_TIMING_PRINT

		//////////////////////////////////////////////////////////////////////////
		// Generate volume map of object using Discregrid
		//////////////////////////////////////////////////////////////////////////

		Simulation *sim = Simulation::getCurrent();
		const bool sim2D = sim->is2DSimulation();

		if (sim2D)
			SimpleQuadrature::determineSamplePointsInCircle(supportRadius, 30);
		auto int_domain = Eigen::AlignedBox3d(Eigen::Vector3d::Constant(-supportRadius), Eigen::Vector3d::Constant(supportRadius));
		Real factor = 1.0;
		if (sim2D)
			factor = 1.0;
		auto volume_func = [&](Eigen::Vector3d const& x)
		{
			auto dist = volumeMap->interpolate(0u, x);
			if (dist > (1.0 + 1.0 /*/ factor*/) * supportRadius)
			{
				return 0.0;
			}

			auto integrand = [&volumeMap, &x, &supportRadius, &factor, &sim, &sim2D](Eigen::Vector3d const& xi) -> double
			{
				if (xi.squaredNorm() > supportRadius*supportRadius)
					return 0.0;

				auto dist = volumeMap->interpolate(0u, x + xi);

				if (dist <= 0.0)
					return 1.0 - 0.001 * dist / supportRadius;
  				if (dist < 1.0 / factor * supportRadius)
  					return static_cast<double>(CubicKernel::W(factor * static_cast<Real>(dist)) / CubicKernel::W_zero());
 				return 0.0;
			};

			double res = 0.0;
			if (sim2D)
				res = 1.2 * SimpleQuadrature::integrate(integrand);
			else
				res = 1.2 * GaussQuadrature::integrate(integrand, int_domain, 30);

			return res;
		};

		auto cell_diag = volumeMap->cellSize().norm();
		std::cout << "Generate volume map..." << std::endl;
		const bool no_reduction = true;
		START_TIMING("Volume Map Construction");
		volumeMap->addFunction(volume_func, false, [&](Eigen::Vector3d const& x_)
		{
			if (no_reduction)
			{
				return true;
			}
			auto x = x_.cwiseMax(volumeMap->domain().min()).cwiseMin(volumeMap->domain().max());
			auto dist = volumeMap->interpolate(0u, x);
			if (dist == std::numeric_limits<double>::max())
			{
				return false;
			}

			return fabs(dist) < 4.0 * supportRadius;
		});
		STOP_TIMING_PRINT;

		// reduction
		if (!no_reduction)
		{
			std::cout << "Reduce discrete fields...";
			volumeMap->reduceField(0u, [&](const Eigen::Vector3d &, double v)
			{
				return fabs(v) < 4.0 * supportRadius;
			});
			volumeMap->reduceField(1u, [&](const Eigen::Vector3d &, double v)->double
			{
				if (v == std::numeric_limits<double>::max())
					return false;				
				return true;
			});
			std::cout << "DONE" << std::endl;
		}

		boundaryModel->setMap(volumeMap);

		// Store cache file
		if (useCache && (FileSystem::makeDir(cachePath) == 0))
		{
			LOG_INFO << "Save volume map: " << volumeMapFileName;
			volumeMap->save(volumeMapFileName);
		}
	}

	// store maximal distance of a point to center of mass for CFL
	if (boundaryData->dynamic)
	{
		// determine center of mass
		Vector3r com;
		com.setZero();
		for (unsigned int i = 0; i < x.size(); i++)
		{
			com += x[i];
		}
		com /= static_cast<Real>(x.size());

		// determine point with maximal distance to center of mass
		Real maxDist = 0.0;
		for (unsigned int i = 0; i < x.size(); i++)
		{
			const Vector3r diff = x[i] - com;
			const Real dist = diff.norm();
			if (dist > maxDist)
			{
				maxDist = dist;
			}
		}
		boundaryModel->setMaxDist(maxDist);
	}
}

void SimulatorBase::determineMinMaxOfScalarField()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const std::string& colorFieldName = getColorField(fluidModelIndex);

		const FieldDescription* field = nullptr;
		field = &model->getField(colorFieldName);

		Real minValue = REAL_MAX;
		Real maxValue = REAL_MIN;
		for (unsigned int i = 0u; i < model->numActiveParticles(); i++)
		{
			if (field->type == FieldType::Vector3)
			{
				Eigen::Map<Vector3r> vec((Real*)field->getFct(i));
				const Real val = vec.norm();
				minValue = std::min(minValue, val);
				maxValue = std::max(maxValue, val);
			}
			else if (field->type == FieldType::Scalar)
			{
				minValue = std::min(minValue, static_cast<Real>(*(Real*)field->getFct(i)));
				maxValue = std::max(maxValue, static_cast<Real>(*(Real*)field->getFct(i)));
			}
			else if (field->type == FieldType::UInt)
			{
				minValue = std::min(minValue, static_cast<Real>(*(unsigned int*)field->getFct(i)));
				maxValue = std::max(maxValue, static_cast<Real>(*(unsigned int*)field->getFct(i)));
			}
			else if (field->type == FieldType::Matrix3)
			{
				Eigen::Map <Matrix3r> vec((Real*)field->getFct(i));
				const Real val = vec.norm();
				minValue = std::min(minValue, static_cast<Real>(val));
				maxValue = std::max(maxValue, static_cast<Real>(val));
			}
			else if (field->type == FieldType::Vector6)
			{
				Eigen::Map<Vector6r> vec((Real*)field->getFct(i));
				const Real val = vec.norm();
				minValue = std::min(minValue, static_cast<Real>(val));
				maxValue = std::max(maxValue, static_cast<Real>(val));
			}
			else if (field->type == FieldType::Matrix6)
			{
				Eigen::Map<Matrix6r> vec((Real*)field->getFct(i));
				const Real val = vec.norm();
				minValue = std::min(minValue, static_cast<Real>(val));
				maxValue = std::max(maxValue, static_cast<Real>(val));
			}
		}
		setRenderMinValue(fluidModelIndex, minValue);
		setRenderMaxValue(fluidModelIndex, maxValue);
	}
	//updateScalarField();
}

