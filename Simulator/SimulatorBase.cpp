#include "SimulatorBase.h"
#include "SPlisHSPlasH/Utilities/SceneLoader.h"
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
int SimulatorBase::PARTIO_EXPORT = -1;
int SimulatorBase::VTK_EXPORT = -1;
int SimulatorBase::RB_VTK_EXPORT = -1;
int SimulatorBase::RB_EXPORT = -1;
int SimulatorBase::DATA_EXPORT_FPS = -1;
int SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES = -1;
int SimulatorBase::STATE_EXPORT = -1;
int SimulatorBase::STATE_EXPORT_FPS = -1;
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
	m_sceneFile = "";
	m_renderWalls = 4;
	m_doPause = true;
	m_pauseAt = -1.0;
	m_stopAt = -1.0;
	m_useParticleCaching = true;
	m_useGUI = true;
	m_enablePartioExport = false;
	m_enableVTKExport = false;
	m_enableRigidBodyVTKExport = false;
	m_enableRigidBodyExport = false;
	m_enableStateExport = false;
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

	PARTIO_EXPORT = createBoolParameter("enablePartioExport", "Partio export", &m_enablePartioExport);
	setGroup(PARTIO_EXPORT, "Export");
	setDescription(PARTIO_EXPORT, "Enable/disable partio export.");

	RB_EXPORT = createBoolParameter("enableRigidBodyExport", "Rigid body export", &m_enableRigidBodyExport);
	setGroup(RB_EXPORT, "Export");
	setDescription(RB_EXPORT, "Enable/disable rigid body export.");

	VTK_EXPORT = createBoolParameter("enableVTKExport", "VTK export", &m_enableVTKExport);
	setGroup(VTK_EXPORT, "Export");
	setDescription(VTK_EXPORT, "Enable/disable VTK export.");

	RB_VTK_EXPORT = createBoolParameter("enableRigidBodyVTKExport", "Rigid body VTK export", &m_enableRigidBodyVTKExport);
	setGroup(RB_VTK_EXPORT, "Export");
	setDescription(RB_VTK_EXPORT, "Enable/disable rigid body VTK export.");

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
}

void SimulatorBase::run()
{
	initSimulation();
	runSimulation();
	cleanup();
}

void SimulatorBase::init(std::vector<std::string> argv, const std::string &windowName){
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

		m_sceneFile = "";

		if (result.count("scene-file"))
		{
			m_sceneFile = result["scene-file"].as<std::string>();
			if (FileSystem::isRelativePath(m_sceneFile))
				m_sceneFile = FileSystem::normalizePath(m_exePath + "/" + m_sceneFile);
		}
#ifdef WIN32
		else
		{
			std::string scenePath = FileSystem::normalizePath(m_exePath + "/../data/Scenes/");
			if (!FileSystem::isDirectory(scenePath))
				scenePath = m_exePath;
			std::replace(scenePath.begin(), scenePath.end(), '/', '\\');
			m_sceneFile = FileSystem::fileDialog(0, scenePath.c_str(), "*.json");

			if (m_sceneFile == "")
				exit(0);
		}
#endif

		m_outputPath = FileSystem::normalizePath(getExePath() + "/output/" + FileSystem::getFileName(m_sceneFile));
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
	if (m_sceneFile != "")
		m_sceneLoader->readScene(m_sceneFile.c_str(), m_scene);
	else
		return;

	//////////////////////////////////////////////////////////////////////////
	// init boundary simulation
	//////////////////////////////////////////////////////////////////////////
	m_isStaticScene = true;
	for (unsigned int i = 0; i < m_scene.boundaryModels.size(); i++)
	{
		if (m_scene.boundaryModels[i]->dynamic)
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
}

void SimulatorBase::initSimulation()
{
#ifdef DL_OUTPUT
		// copy scene files in output so that the simulation can be reproduced
	std::string sceneFilePath = FileSystem::normalizePath(m_outputPath + "/scene");
	FileSystem::makeDirs(sceneFilePath);
	FileSystem::copyFile(m_sceneFile, sceneFilePath + "/" + FileSystem::getFileNameWithExt(m_sceneFile));

	std::string modelsFilePath = FileSystem::normalizePath(m_outputPath + "/models");
	FileSystem::makeDirs(modelsFilePath);
	for (unsigned int i = 0; i < m_scene.boundaryModels.size(); i++)
	{
		std::string meshFileName = m_scene.boundaryModels[i]->meshFile;
		if (meshFileName != "")
		{
			if (FileSystem::isRelativePath(meshFileName))
				meshFileName = FileSystem::normalizePath(FileSystem::getFilePath(m_sceneFile) + "/" + meshFileName);
			FileSystem::copyFile(meshFileName, modelsFilePath + "/" + FileSystem::getFileNameWithExt(meshFileName));
		}
		std::string mapFileName = m_scene.boundaryModels[i]->mapFile;
		if (mapFileName != "")
		{
			if (FileSystem::isRelativePath(mapFileName))
				mapFileName = FileSystem::normalizePath(FileSystem::getFilePath(m_sceneFile) + "/" + mapFileName);
			FileSystem::copyFile(m_scene.boundaryModels[i]->mapFile, modelsFilePath + "/" + FileSystem::getFileNameWithExt(m_scene.boundaryModels[i]->mapFile));
		}
	}

	for (unsigned int i = 0; i < m_scene.fluidModels.size(); i++)
	{
		std::string samplesFileName = m_scene.fluidModels[i]->samplesFile;
		if (samplesFileName != "")
		{
			if (FileSystem::isRelativePath(samplesFileName))
				samplesFileName = FileSystem::normalizePath(FileSystem::getFilePath(m_sceneFile) + "/" + samplesFileName);
			FileSystem::copyFile(samplesFileName, modelsFilePath + "/" + FileSystem::getFileNameWithExt(samplesFileName));
		}
	}

	std::string progFilePath = FileSystem::normalizePath(m_outputPath + "/program");
	FileSystem::makeDirs(progFilePath);
	FileSystem::copyFile(m_argv[0], progFilePath + "/" + FileSystem::getFileNameWithExt(m_argv[0]));
	#endif

	Simulation *sim = Simulation::getCurrent();
	sim->init(getScene().particleRadius, getScene().sim2D);

	buildModel();

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
			model->setElasticityMethodChangedCallback([this, model]() { reset(); m_gui->initSimulationParameterGUI(); getSceneLoader()->readMaterialParameterObject(model->getId(), (ParameterObject*)model->getElasticityBase()); });
		}

		m_gui->initSimulationParameterGUI();
		Simulation::getCurrent()->setSimulationMethodChangedCallback([this]() { reset(); m_gui->initSimulationParameterGUI(); getSceneLoader()->readParameterObject("Configuration", Simulation::getCurrent()->getTimeStep()); });
	}
	readParameters();
	setCommandLineParameter();
}

void SimulatorBase::runSimulation()
{
	m_boundarySimulator->initBoundaryData();

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
	for (unsigned int i = 0; i < m_scene.boundaryModels.size(); i++)
		delete m_scene.boundaryModels[i];
	m_scene.boundaryModels.clear();

	for (unsigned int i = 0; i < m_scene.fluidModels.size(); i++)
		delete m_scene.fluidModels[i];
	m_scene.fluidModels.clear();

	for (unsigned int i = 0; i < m_scene.fluidBlocks.size(); i++)
		delete m_scene.fluidBlocks[i];
	m_scene.fluidBlocks.clear();

	for (unsigned int i = 0; i < m_scene.emitters.size(); i++)
		delete m_scene.emitters[i];
	m_scene.emitters.clear();

	for (unsigned int i = 0; i < m_scene.animatedFields.size(); i++)
		delete m_scene.animatedFields[i];
	m_scene.animatedFields.clear();

	if (m_useGUI)
		m_gui->cleanup();

	delete Simulation::getCurrent();
}

void SimulatorBase::readParameters()
{
	m_sceneLoader->readParameterObject("Configuration", this);
	m_sceneLoader->readParameterObject("Configuration", Simulation::getCurrent());
	m_sceneLoader->readParameterObject("Configuration", Simulation::getCurrent()->getTimeStep());

	Simulation *sim = Simulation::getCurrent();
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

		for (auto material : m_scene.materials)
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


void SimulatorBase::buildModel()
{
	TimeManager::getCurrent()->setTimeStepSize(m_scene.timeStepSize);

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
	m_boundarySimulator->reset();
	if (m_gui)
		m_gui->reset();

	if (Simulation::getCurrent()->getValue<int>(Simulation::CFL_METHOD) != Simulation::ENUM_CFL_NONE)
		TimeManager::getCurrent()->setTimeStepSize(m_scene.timeStepSize);
	m_nextFrameTime = 0.0;
	m_nextFrameTimeState = 0.0;
	m_frameCounter = 1;
	m_isFirstFrame = true;
	m_isFirstFrameVTK = true;
#ifdef DL_OUTPUT
	m_nextTiming = 1.0;
#endif
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

		if (m_timeStepCB)
			m_timeStepCB();
	}
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
	if (m_timeStepCB)
		m_timeStepCB();
	return true;
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


void SimulatorBase::initFluidData()
{
	LOG_INFO << "Initialize fluid particles";

	Simulation *sim = Simulation::getCurrent();

	//////////////////////////////////////////////////////////////////////////
	// Determine number of different fluid IDs
	//////////////////////////////////////////////////////////////////////////
	std::map<std::string, unsigned int> fluidIDs;
	unsigned int index = 0;
	for (unsigned int i = 0; i < m_scene.fluidBlocks.size(); i++)
	{
		if (fluidIDs.find(m_scene.fluidBlocks[i]->id) == fluidIDs.end())
			fluidIDs[m_scene.fluidBlocks[i]->id] = index++;
	}
	for (unsigned int i = 0; i < m_scene.fluidModels.size(); i++)
	{
		if (fluidIDs.find(m_scene.fluidModels[i]->id) == fluidIDs.end())
			fluidIDs[m_scene.fluidModels[i]->id] = index++;
	}
	for (unsigned int i = 0; i < m_scene.emitters.size(); i++)
	{
		if (fluidIDs.find(m_scene.emitters[i]->id) == fluidIDs.end())
			fluidIDs[m_scene.emitters[i]->id] = index++;
	}
	const unsigned int numberOfFluidModels = static_cast<unsigned int>(fluidIDs.size());

	std::vector<std::vector<Vector3r>> fluidParticles;
	std::vector<std::vector<Vector3r>> fluidVelocities;
	fluidParticles.resize(numberOfFluidModels);
	fluidVelocities.resize(numberOfFluidModels);

	createFluidBlocks(fluidIDs, fluidParticles, fluidVelocities);

	std::string base_path = FileSystem::getFilePath(m_sceneFile);

	const bool useCache = getUseParticleCaching();
	std::string scene_path = FileSystem::getFilePath(getSceneFile());
	string cachePath = scene_path + "/Cache";

	unsigned int startIndex = 0;
	unsigned int endIndex = 0;
	for (unsigned int i = 0; i < m_scene.fluidModels.size(); i++)
	{
		const unsigned int fluidIndex = fluidIDs[m_scene.fluidModels[i]->id];

		std::string fileName;
		if (FileSystem::isRelativePath(m_scene.fluidModels[i]->samplesFile))
			fileName = base_path + "/" + m_scene.fluidModels[i]->samplesFile;
		else
			fileName = m_scene.fluidModels[i]->samplesFile;

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

			std::string mode = to_string(m_scene.fluidModels[i]->mode);
			const string scaleStr = real2String(m_scene.fluidModels[i]->scale[0]) + "_" + real2String(m_scene.fluidModels[i]->scale[1]) + "_" + real2String(m_scene.fluidModels[i]->scale[2]);
			const string resStr = to_string(m_scene.fluidModels[i]->resolutionSDF[0]) + "_" + to_string(m_scene.fluidModels[i]->resolutionSDF[1]) + "_" + to_string(m_scene.fluidModels[i]->resolutionSDF[2]);
			const string particleFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_fluid_" + real2String(m_scene.particleRadius) + "_m" + mode + "_s" + scaleStr + "_r" + resStr + ".bgeo");

			// check MD5 if cache file is available
			bool foundCacheFile = false;

			if (useCache)
				foundCacheFile = FileSystem::fileExists(particleFileName);

			if (useCache && foundCacheFile && md5)
			{
				PartioReaderWriter::readParticles(particleFileName, m_scene.fluidModels[i]->translation, m_scene.fluidModels[i]->rotation, 1.0, fluidParticles[fluidIndex], fluidVelocities[fluidIndex]);
				LOG_INFO << "Loaded cached fluid sampling: " << particleFileName;
			}

			if (!useCache || !foundCacheFile || !md5)
			{
				LOG_INFO << "Volume sampling of " << fileName;

				TriangleMesh mesh;
				loadObj(fileName, mesh, m_scene.fluidModels[i]->scale);

				bool invert = m_scene.fluidModels[i]->invert;
				int mode = m_scene.fluidModels[i]->mode;
				std::array<unsigned int, 3> resolutionSDF = m_scene.fluidModels[i]->resolutionSDF;

				LOG_INFO << "SDF resolution: " << resolutionSDF[0] << ", " << resolutionSDF[1] << ", " << resolutionSDF[2];
					
				START_TIMING("Volume sampling");
				Utilities::VolumeSampling::sampleMesh(mesh.numVertices(), mesh.getVertices().data(), mesh.numFaces(), mesh.getFaces().data(),
					m_scene.particleRadius, nullptr, resolutionSDF, invert, mode, fluidParticles[fluidIndex]);
				STOP_TIMING_AVG;

				fluidVelocities[fluidIndex].resize(fluidParticles[fluidIndex].size(), m_scene.fluidModels[i]->initialVelocity);

				// Cache sampling
				if (useCache && (FileSystem::makeDir(cachePath) == 0))
				{
					LOG_INFO << "Save particle sampling: " << particleFileName;
					PartioReaderWriter::writeParticles(particleFileName, (unsigned int)fluidParticles[fluidIndex].size(), fluidParticles[fluidIndex].data(), fluidVelocities[fluidIndex].data(), 0.0);
					FileSystem::writeMD5File(fileName, md5FileName);
				}

				// transform particles
				for (unsigned int j = 0; j < (unsigned int)fluidParticles[fluidIndex].size(); j++)
					fluidParticles[fluidIndex][j] = m_scene.fluidModels[i]->rotation * fluidParticles[fluidIndex][j] + m_scene.fluidModels[i]->translation;
			}
		}
		else
		{
			if (!PartioReaderWriter::readParticles(fileName, m_scene.fluidModels[i]->translation, m_scene.fluidModels[i]->rotation, m_scene.fluidModels[i]->scale[0], fluidParticles[fluidIndex], fluidVelocities[fluidIndex]))
				LOG_ERR << "File not found: " << fileName;
		}
		Simulation::getCurrent()->setValue(Simulation::PARTICLE_RADIUS, m_scene.particleRadius);
	}

	unsigned int nParticles = 0;
	for (auto it = fluidIDs.begin(); it != fluidIDs.end(); it++)
	{
		const unsigned int index = it->second;

		unsigned int maxEmitterParticles = 10000;
		for (auto material : m_scene.materials)
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

	//////////////////////////////////////////////////////////////////////////
	// emitters
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int i = 0; i < m_scene.emitters.size(); i++)
	{
		SceneLoader::EmitterData *ed = m_scene.emitters[i];

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
			m_scene.boundaryModels.push_back(emitterBoundary);
			
			// reuse particles if they are outside of a bounding box
			bool emitterReuseParticles = false;
			Vector3r emitterBoxMin(-1.0, -1.0, -1.0);
			Vector3r emitterBoxMax(1.0, 1.0, 1.0);
			for (auto material : m_scene.materials)
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

	//////////////////////////////////////////////////////////////////////////
	// animation fields
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int i = 0; i < m_scene.animatedFields.size(); i++)
	{
		SceneLoader::AnimationFieldData *data = m_scene.animatedFields[i];

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
	for (unsigned int i = 0; i < m_scene.fluidBlocks.size(); i++)
	{
		const unsigned int fluidIndex = fluidIDs[m_scene.fluidBlocks[i]->id];
		const Real diam = static_cast<Real>(2.0)*m_scene.particleRadius;

		Real xshift = diam;
		Real yshift = diam;
		const Real eps = static_cast<Real>(1.0e-9);
		if (m_scene.fluidBlocks[i]->mode == 1)
			yshift = sqrt(static_cast<Real>(3.0)) * m_scene.particleRadius + eps;
		else if (m_scene.fluidBlocks[i]->mode == 2)
		{
			xshift = sqrt(static_cast<Real>(6.0)) * diam / static_cast<Real>(3.0) + eps;
			yshift = sqrt(static_cast<Real>(3.0)) * m_scene.particleRadius + eps;
		}

		Vector3r diff = m_scene.fluidBlocks[i]->box.m_maxX - m_scene.fluidBlocks[i]->box.m_minX;
		if (m_scene.fluidBlocks[i]->mode == 1)
		{
			diff[0] -= diam;
			diff[2] -= diam;
		}
		else if (m_scene.fluidBlocks[i]->mode == 2)
		{
			diff[0] -= xshift;
			diff[2] -= diam;
		}

		const int stepsX = (int)round(diff[0] / xshift) - 1;
		const int stepsY = (int)round(diff[1] / yshift) - 1;
		int stepsZ = (int)round(diff[2] / diam) - 1;

		Vector3r start = m_scene.fluidBlocks[i]->box.m_minX + static_cast<Real>(2.0)*m_scene.particleRadius*Vector3r::Ones();
		fluidParticles[fluidIndex].reserve(fluidParticles[fluidIndex].size() + stepsX*stepsY*stepsZ);
		fluidVelocities[fluidIndex].resize(fluidVelocities[fluidIndex].size() + stepsX*stepsY*stepsZ, m_scene.fluidBlocks[i]->initialVelocity);

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
					if (m_scene.fluidBlocks[i]->mode == 1)
					{
						if (k % 2 == 0)
							currPos += Vector3r(0, 0, m_scene.particleRadius);
						else
							currPos += Vector3r(m_scene.particleRadius, 0, 0);
					}
					else if (m_scene.fluidBlocks[i]->mode == 2)
					{
						currPos += Vector3r(0, 0, m_scene.particleRadius);

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

void SimulatorBase::rigidBodyExport()
{
	std::string exportPath = FileSystem::normalizePath(m_outputPath + "/rigid_bodies");
	std::string exportPathVTK = FileSystem::normalizePath(m_outputPath + "/vtk");
	if (m_enableRigidBodyExport)
	{
		FileSystem::makeDirs(exportPath);
		writeRigidBodiesBIN(exportPath);
	}
	if (m_enableRigidBodyVTKExport)
	{
		FileSystem::makeDirs(exportPathVTK);
		writeRigidBodiesVTK(exportPathVTK);
	}
}

void SimulatorBase::particleExport()
{	
	std::string partioExportPath = FileSystem::normalizePath(m_outputPath + "/partio");
	std::string vtkExportPath = FileSystem::normalizePath(m_outputPath + "/vtk");
	if (m_enablePartioExport)
		FileSystem::makeDirs(partioExportPath);
	if (m_enableVTKExport)
		FileSystem::makeDirs(vtkExportPath);

	Simulation *sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		std::string fileName = "ParticleData";
		fileName = fileName + "_" + model->getId() + "_" + std::to_string(m_frameCounter);

		if (m_enablePartioExport)
		{
			std::string exportFileName = FileSystem::normalizePath(partioExportPath + "/" + fileName);
			writeParticlesPartio(exportFileName + ".bgeo", model);
		}
		if (m_enableVTKExport)
		{
			std::string exportFileName = FileSystem::normalizePath(vtkExportPath + "/" + fileName);
			writeParticlesVTK(exportFileName + ".vtk", model);
		}
	}
}


void SimulatorBase::writeParticlesPartio(const std::string &fileName, FluidModel *model)
{
	Partio::ParticlesDataMutable& particleData = *Partio::create();
	Partio::ParticleAttribute posAttr = particleData.addAttribute("position", Partio::VECTOR, 3);
	Partio::ParticleAttribute idAttr = particleData.addAttribute("id", Partio::INT, 1);

	// add attributes
	std::vector<std::string> attributes;
	StringTools::tokenize(m_particleAttributes, attributes, ";");

 	std::map<unsigned int, int> attrMap;
 	std::map<unsigned int, Partio::ParticleAttribute> partioAttrMap;
 	for (unsigned int i = 0; i < attributes.size(); i++)
 	{
 		// position is exported anyway
		if (attributes[i] == "position")
		{
			attrMap[i] = -1;
			continue;
		}
 
 		bool found = false;
 		for (unsigned int j = 0; j < model->numberOfFields(); j++)
 		{
 			const FieldDescription &field = model->getField(j);
 			if (field.name == attributes[i])
 			{
 				found = true;
 				if (field.type == Scalar)
 				{
 					attrMap[i] = j;
 					partioAttrMap[i] = particleData.addAttribute(attributes[i].c_str(), Partio::FLOAT, 1);
 				}
				else if (field.type == UInt)
				{
					attrMap[i] = j;
					partioAttrMap[i] = particleData.addAttribute(attributes[i].c_str(), Partio::INT, 1);
				}
 				else if (field.type == Vector3)
 				{
 					attrMap[i] = j;
 					partioAttrMap[i] = particleData.addAttribute(attributes[i].c_str(), Partio::VECTOR, 3);
 				}
 				else
 				{
 					attrMap[i] = -1;
 					LOG_WARN << "Only scalar and vector fields are currently supported by the partio exporter.";
 				}
 				break;
 			}
 		}
		if (!found)
		{
			attrMap[i] = -1;
			LOG_WARN << "Unknown field cannot be exported in partio file: " << attributes[i];
		}
 	}

	const unsigned int numParticles = model->numActiveParticles();

	for (unsigned int i = 0; i < numParticles; i++)
	{
		Partio::ParticleIndex index = particleData.addParticle();
		float* pos = particleData.dataWrite<float>(posAttr, index);
		int* id = particleData.dataWrite<int>(idAttr, index);

		const Vector3r &x = model->getPosition(i);
		pos[0] = (float)x[0];
		pos[1] = (float)x[1];
		pos[2] = (float)x[2];
	
		id[0] = model->getParticleId(i);

 		for (unsigned int j = 0; j < attributes.size(); j++)
 		{
 			const int fieldIndex = attrMap[j];
 			if (fieldIndex != -1)
 			{
 				const FieldDescription &field = model->getField(fieldIndex);
 				if (field.type == FieldType::Scalar)
 				{
 					float* val = particleData.dataWrite<float>(partioAttrMap[j], index);
 					*val = (float) *((Real*) field.getFct(i));
 				}
				else if (field.type == FieldType::UInt)
				{
					int* val = particleData.dataWrite<int>(partioAttrMap[j], index);
					*val = (int) *((unsigned int*)field.getFct(i));
				}
 				else if (field.type == FieldType::Vector3)
 				{
 					float* val = particleData.dataWrite<float>(partioAttrMap[j], index);
					Eigen::Map<Vector3r> vec((Real*) field.getFct(i));
 					val[0] = (float)vec[0];
 					val[1] = (float)vec[1];
 					val[2] = (float)vec[2];
 				}
 			}
 		}
	}

	Partio::write(fileName.c_str(), particleData, true);
	particleData.release();
}

void SPH::SimulatorBase::writeParticlesVTK(const std::string &fileName, FluidModel *model)
{
	const unsigned int numParticles = model->numActiveParticles();
	if (0 == numParticles)
		return;

#ifdef USE_DOUBLE
	const char * real_str = " double\n";
#else 
	const char * real_str = " float\n";
#endif

	// Open the file
	std::ofstream outfile{ fileName, std::ios::binary };
	if (!outfile.is_open()) 
	{
		LOG_WARN << "Cannot open a file to save VTK particles.";
		return;
	}

	outfile << "# vtk DataFile Version 4.1\n";
	outfile << "SPlisHSPlasH particle data\n"; // title of the data set, (any string up to 256 characters+\n)
	outfile << "BINARY\n";
	outfile << "DATASET UNSTRUCTURED_GRID\n";

	// add attributes
	std::vector<std::string> attributes;
	StringTools::tokenize(m_particleAttributes, attributes, ";");

	//////////////////////////////////////////////////////////////////////////
	// positions and ids exported anyways
	attributes.erase(
		std::remove_if(attributes.begin(), attributes.end(), [](const std::string&s) { return (s == "position" || s == "id"); }),
		attributes.end());

	//////////////////////////////////////////////////////////////////////////
	// export position attribute as POINTS
	{
		std::vector<Vector3r> positions;
		positions.reserve(numParticles);
		for (unsigned int i = 0u; i < numParticles; i++)
			positions.emplace_back(model->getPosition(i));
		// swap endianess
		for (unsigned int i = 0; i < numParticles; i++)
			for (unsigned int c = 0; c < 3; c++)
				swapByteOrder(&positions[i][c]);
		// export to vtk
		outfile << "POINTS " << numParticles << real_str;
		outfile.write(reinterpret_cast<char*>(positions[0].data()), 3 * numParticles * sizeof(Real));
		outfile << "\n";
	}

	//////////////////////////////////////////////////////////////////////////
	// export particle IDs as CELLS
	{
		std::vector<Eigen::Vector2i> cells;
		cells.reserve(numParticles);
		unsigned int nodes_per_cell_swapped = 1;
		swapByteOrder(&nodes_per_cell_swapped);
		for (unsigned int i = 0u; i < numParticles; i++)
		{
			unsigned int idSwapped = model->getParticleId(i);
			swapByteOrder(&idSwapped);
			cells.emplace_back(nodes_per_cell_swapped, idSwapped);
		}

		// particles are cells with one element and the index of the particle
		outfile << "CELLS " << numParticles << " " << 2 * numParticles << "\n";
		outfile.write(reinterpret_cast<char*>(cells[0].data()), 2 * numParticles * sizeof(unsigned int));
		outfile << "\n";
	}
	//////////////////////////////////////////////////////////////////////////
	// export cell types
	{
		// the type of a particle cell is always 1
		std::vector<int> cellTypes;
		int cellTypeSwapped = 1;
		swapByteOrder(&cellTypeSwapped);
		cellTypes.resize(numParticles, cellTypeSwapped);
		outfile << "CELL_TYPES " << numParticles << "\n";
		outfile.write(reinterpret_cast<char*>(cellTypes.data()), numParticles * sizeof(int));
		outfile << "\n";
	}

	//////////////////////////////////////////////////////////////////////////
	// write additional attributes as per-particle data
	{
		outfile << "POINT_DATA " << numParticles << "\n";
		// write IDs
		outfile << "SCALARS id unsigned_int 1\n";
		outfile << "LOOKUP_TABLE id_table\n";
		// copy data
		std::vector<unsigned int> attrData;
		attrData.reserve(numParticles);
		for (unsigned int i = 0u; i < numParticles; i++)
			attrData.emplace_back(model->getParticleId(i));
		// swap endianess
		for (unsigned int i = 0; i < numParticles; i++)
			swapByteOrder(&attrData[i]);
		// export to vtk
		outfile.write(reinterpret_cast<char*>(attrData.data()), numParticles * sizeof(unsigned int));
		outfile << "\n";
	}

	//////////////////////////////////////////////////////////////////////////
	// per point fields (all attributes except for positions)
	const auto numFields = attributes.size();
	outfile << "FIELD FieldData " << std::to_string(numFields) << "\n";

	// iterate over attributes
	for (const std::string & a : attributes)
	{
		const FieldDescription & field = model->getField(a);

		std::string attrNameVTK;
		std::regex_replace(std::back_inserter(attrNameVTK), a.begin(), a.end(), std::regex("\\s+"), "_");

		if (field.type == FieldType::Scalar)
		{
			// write header information
			outfile << attrNameVTK << " 1 " << numParticles << real_str;

			// copy data
			std::vector<Real> attrData;
			attrData.reserve(numParticles);
			for (unsigned int i = 0u; i < numParticles; i++)
				attrData.emplace_back(*((Real*) field.getFct(i)));
			// swap endianess
			for (unsigned int i = 0; i < numParticles; i++)
				swapByteOrder(&attrData[i]);
			// export to vtk
			outfile.write(reinterpret_cast<char*>(attrData.data()), numParticles * sizeof(Real));
		}
		else if (field.type == FieldType::Vector3)
		{
			// write header information
			outfile << attrNameVTK << " 3 " << numParticles << real_str;

			// copy from partio data
			std::vector<Vector3r> attrData;
			attrData.reserve(numParticles);
			for (unsigned int i = 0u; i < numParticles; i++)
				attrData.emplace_back((Real*) field.getFct(i));
			// swap endianess
			for (unsigned int i = 0; i < numParticles; i++)
				for (unsigned int c = 0; c < 3; c++)
					swapByteOrder(&attrData[i][c]);
			// export to vtk
			outfile.write(reinterpret_cast<char*>(attrData[0].data()), 3 * numParticles * sizeof(Real));
		}
		// TODO support other field types
		else
		{
			LOG_WARN << "Skipping attribute " << a << ", because it is of unsupported type\n";
			continue;
		}
		// end of block
		outfile << "\n";
	}
	outfile.close();
}

void SimulatorBase::writeRigidBodiesBIN(const std::string &exportPath)
{
	std::string fileName = "rb_data_";
	fileName = fileName + std::to_string(m_frameCounter) + ".bin";
	std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nBoundaryModels = sim->numberOfBoundaryModels();

	std::string scene_path = FileSystem::getFilePath(getSceneFile());

	// check if we have a static model
	bool isStatic = true;
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel *bm = sim->getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic())
		{
			isStatic = false;
			break;
		}
	}

	BinaryFileWriter binWriter;
	if (m_isFirstFrame || !isStatic)
		binWriter.openFile(exportFileName.c_str());

	if (m_isFirstFrame)
	{
		binWriter.write(nBoundaryModels);

		for (unsigned int i = 0; i < m_scene.boundaryModels.size(); i++)
		{
			std::string meshFileName = m_scene.boundaryModels[i]->meshFile;
			if (FileSystem::isRelativePath(meshFileName))
				meshFileName = FileSystem::normalizePath(scene_path + "/" + meshFileName);

			const string fileName = Utilities::FileSystem::getFileNameWithExt(meshFileName);
			binWriter.write(fileName);
			Eigen::Vector3f s = m_scene.boundaryModels[i]->scale.template cast<float>();
			binWriter.writeMatrix(s);
			std::string targetFilePath = exportPath + "/" + fileName;
			if (!Utilities::FileSystem::fileExists(targetFilePath))
			{
				Utilities::FileSystem::copyFile(meshFileName, targetFilePath);
			}
			binWriter.write((char)m_scene.boundaryModels[i]->isWall);
			binWriter.writeMatrix(m_scene.boundaryModels[i]->color);
		}
	}

	if (m_isFirstFrame || !isStatic)
	{
		for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
		{
			BoundaryModel *bm = sim->getBoundaryModel(i);
			const Vector3r &x = bm->getRigidBodyObject()->getWorldSpacePosition();
			const Eigen::Vector3f x_f = x.template cast<float>();
			binWriter.writeMatrix(x_f);

			const Matrix3r &R = bm->getRigidBodyObject()->getWorldSpaceRotation();
			//const Eigen::Matrix3f RT = R.transpose().template cast<float>();
			binWriter.writeMatrix(R);
		}
		binWriter.closeFile();
	}

	m_isFirstFrame = false;
}

void SimulatorBase::writeRigidBodiesVTK(const std::string &exportPath)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nBoundaryModels = sim->numberOfBoundaryModels();

	// check if we have a static model
	bool isStatic = true;
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel *bm = sim->getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic())
		{
			isStatic = false;
			break;
		}
	}

#ifdef USE_DOUBLE
	const char * real_str = " double\n";
#else 
	const char * real_str = " float\n";
#endif

	if (m_isFirstFrameVTK || !isStatic)
	{
		for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
		{
			std::string fileName = "rb_data_";
			fileName = fileName + std::to_string(i) + "_" + std::to_string(m_frameCounter) + ".vtk";
			std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

			// Open the file
			std::ofstream outfile(exportFileName, std::ios::binary);
			if (!outfile)
			{
				LOG_WARN << "Cannot open a file to save VTK mesh.";
				return;
			}

			// Header
			outfile << "# vtk DataFile Version 4.2\n";
			outfile << "SPlisHSPlasH mesh data\n";
			outfile << "BINARY\n";
			outfile << "DATASET UNSTRUCTURED_GRID\n";

			BoundaryModel *bm = sim->getBoundaryModel(i);
			const std::vector<Vector3r> &vertices = bm->getRigidBodyObject()->getVertices();
			const std::vector<unsigned int> &faces = bm->getRigidBodyObject()->getFaces();
			int n_vertices = (int)vertices.size();
			int n_triangles = (int)faces.size() / 3;

			// Vertices
			{
				std::vector<Vector3r> positions;
				positions.reserve(n_vertices);
				for (int j = 0u; j < n_vertices; j++)
				{
					Vector3r x = vertices[j];
					swapByteOrder(&x[0]);
					swapByteOrder(&x[1]);
					swapByteOrder(&x[2]);
					positions.emplace_back(x);
				}
				// export to vtk
				outfile << "POINTS " << n_vertices << real_str;
				outfile.write(reinterpret_cast<char*>(positions[0].data()), 3 * n_vertices * sizeof(Real));
				outfile << "\n";
			}

			// Connectivity
			{
				std::vector<int> connectivity_to_write;
				connectivity_to_write.reserve(4 * n_triangles);
				for (int tri_i = 0; tri_i < n_triangles; tri_i++)
				{
					int val = 3;
					swapByteOrder(&val);
					connectivity_to_write.push_back(val);
					val = faces[3 * tri_i + 0];
					swapByteOrder(&val);
					connectivity_to_write.push_back(val);
					val = faces[3 * tri_i + 1];
					swapByteOrder(&val);
					connectivity_to_write.push_back(val);
					val = faces[3 * tri_i + 2];
					swapByteOrder(&val);
					connectivity_to_write.push_back(val);
				}
				// export to vtk
				outfile << "CELLS " << n_triangles << " " << 4 * n_triangles << "\n";
				outfile.write(reinterpret_cast<char*>(&connectivity_to_write[0]), connectivity_to_write.size() * sizeof(int));
				outfile << "\n";
			}

			// Cell types
			{
				outfile << "CELL_TYPES " << n_triangles << "\n";
				int cell_type_swapped = 5;
				swapByteOrder(&cell_type_swapped);
				std::vector<int> cell_type_arr(n_triangles, cell_type_swapped);
				outfile.write(reinterpret_cast<char*>(&cell_type_arr[0]), cell_type_arr.size() * sizeof(int));
				outfile << "\n";
			}
			outfile.close();
		}
	}

	m_isFirstFrameVTK = false;
}

void SimulatorBase::step()
{
	if (TimeManager::getCurrent()->getTime() >= m_nextFrameTime)
	{
		m_nextFrameTime += static_cast<Real>(1.0) / m_framesPerSecond;
		if (m_enablePartioExport || m_enableVTKExport)
			particleExport();
		if (m_enableRigidBodyExport || m_enableRigidBodyVTKExport)
			rigidBodyExport();
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
	SceneLoader::Scene &scene = getScene();
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
	SceneLoader::Scene &scene = getScene();
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
	SceneLoader::Scene &scene = getScene();
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


std::string SimulatorBase::real2String(const Real r)
{
	string str = to_string(r);
	str.erase(str.find_last_not_of('0') + 1, std::string::npos);
	str.erase(str.find_last_not_of('.') + 1, std::string::npos);
	return str;
}

void SPH::SimulatorBase::saveState()
{
	std::string stateFilePath = FileSystem::normalizePath(m_outputPath + "/state");
	FileSystem::makeDirs(stateFilePath);

	string md5Str = FileSystem::getFileMD5(m_sceneFile);

	Simulation *sim = Simulation::getCurrent();

	const Real time = TimeManager::getCurrent()->getTime();
	const std::string timeStr = real2String(time);

	// Save additional data
	BinaryFileWriter binWriter;
	std::string exportFileName = FileSystem::normalizePath(stateFilePath + "/state_" + timeStr);
	binWriter.openFile(exportFileName + ".bin");
	binWriter.write(md5Str);

	binWriter.write(m_nextFrameTime);
	binWriter.write(m_nextFrameTimeState);
	binWriter.write(m_frameCounter);
	binWriter.write(m_isFirstFrame);
	
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
			std::string fileName = "state_boundary";
			fileName = fileName + "_" + to_string(i); // +"_" + std::to_string(m_frameCounter);

			// Save particle data
			std::string expFileName = FileSystem::normalizePath(stateFilePath + "/" + fileName);
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
	string md5Str = FileSystem::getFileMD5(m_sceneFile);

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
		std::string fileName = "state_boundary";
		fileName = fileName + "_" + to_string(i); // +"_" + std::to_string(m_frameCounter);

		// Save particle data
		std::string impFileName = FileSystem::normalizePath(importFilePath + "/" + fileName);
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

	sim->performNeighborhoodSearchSort();
}

void SimulatorBase::writeParameterState(BinaryFileWriter &binWriter)
{
	writeParameterObjectState(binWriter, this);
	writeParameterObjectState(binWriter, Simulation::getCurrent());
	writeParameterObjectState(binWriter, Simulation::getCurrent()->getTimeStep());

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
	const Real supportRadius = sim->getSupportRadius();
	std::string scene_path = FileSystem::getFilePath(getSceneFile());
	std::string scene_file_name = FileSystem::getFileName(getSceneFile());
	SceneLoader::Scene &scene = getScene();
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
	const string scaleStr = "s" + real2String(boundaryData->scale[0]) + "_" + real2String(boundaryData->scale[1]) + "_" + real2String(boundaryData->scale[2]);
	const string resStr = "r" + to_string(resolutionSDF[0]) + "_" + to_string(resolutionSDF[1]) + "_" + to_string(resolutionSDF[2]);
	const string invertStr = "i" + to_string((int)boundaryData->mapInvert);
	const string thicknessStr = "t" + real2String(boundaryData->mapThickness);
	const string kernelStr = "k" + to_string(sim->getKernel());
	string densityMapFileName = "";
	if (isDynamic)
		densityMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_db_dm_" + real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + "_" + kernelStr + ".cdm");
	else 
		densityMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_sb_dm_" + real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + "_" + kernelStr + ".cdm");

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
	const Real supportRadius = sim->getSupportRadius();
	std::string scene_path = FileSystem::getFilePath(getSceneFile());
	std::string scene_file_name = FileSystem::getFileName(getSceneFile());
	SceneLoader::Scene &scene = getScene();
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
	const string scaleStr = "s" + real2String(boundaryData->scale[0]) + "_" + real2String(boundaryData->scale[1]) + "_" + real2String(boundaryData->scale[2]);
	const string resStr = "r" + to_string(resolutionSDF[0]) + "_" + to_string(resolutionSDF[1]) + "_" + to_string(resolutionSDF[2]);
	const string invertStr = "i" + to_string((int)boundaryData->mapInvert);
	const string thicknessStr = "t" + real2String(boundaryData->mapThickness);
	string volumeMapFileName = "";
	if (isDynamic)
		volumeMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_db_vm_" + real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + ".cdm");
	else 
		volumeMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_sb_vm_" + real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + ".cdm");

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
			factor = 1.75;
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
					return 1.0 - 0.1 * dist / supportRadius;
  				if (dist < 1.0 / factor * supportRadius)
  					return static_cast<double>(CubicKernel::W(factor * static_cast<Real>(dist)) / CubicKernel::W_zero());
 				return 0.0;
			};

			double res = 0.0;
			if (sim2D)
				res = 0.8 * SimpleQuadrature::integrate(integrand);
			else
				res = 0.8 * GaussQuadrature::integrate(integrand, int_domain, 30);

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
