#include "SimulatorBase.h"
#include "Visualization/MiniGL.h"
#include "SPlisHSPlasH/Utilities/SceneLoader.h"
#include "Utilities/FileSystem.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/PartioReaderWriter.h"
#include "Visualization/Selection.h"
#include "GL/glut.h"
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
#include "Visualization/colormaps/colormap_jet.h"
#include "Visualization/colormaps/colormap_plasma.h"
#include "extern/partio/src/lib/Partio.h"
#include "SPlisHSPlasH/Utilities/GaussQuadrature.h"
#include "SPlisHSPlasH/Utilities/SimpleQuadrature.h"
#include "extern/cxxopts/cxxopts.hpp"

#include "SPlisHSPlasH/Utilities/VolumeSampling.h"
#include "Utilities/OBJLoader.h"
#include "Utilities/BinaryFileReaderWriter.h"

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
	m_enableRigidBodyExport = false;
	m_enableStateExport = false;
	m_framesPerSecond = 25;
	m_framesPerSecondState = 1;
	m_nextFrameTime = 0.0;
	m_nextFrameTimeState = 0.0;
	m_frameCounter = 1;
	m_isFirstFrame = true;
	m_firstState = true;
	m_colorField.resize(1, "velocity");
	m_colorMapType.resize(1, 0);
	m_renderMinValue.resize(1, 0.0);
	m_renderMaxValue.resize(1, 5.0);
	m_particleAttributes = "velocity";
#ifdef DL_OUTPUT
	m_nextTiming = 1.0;
#endif
}

SimulatorBase::~SimulatorBase()
{

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

	VTK_EXPORT = createBoolParameter("enableVTKExport", "VTK export", &m_enableVTKExport);
	setGroup(VTK_EXPORT, "Export");
	setDescription(VTK_EXPORT, "Enable/disable VTK export.");

	RB_EXPORT = createBoolParameter("enableRigidBodyExport", "Rigid body export", &m_enableRigidBodyExport);
	setGroup(RB_EXPORT, "Export");
	setDescription(RB_EXPORT, "Enable/disable rigid body export.");

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

void SimulatorBase::init(int argc, char **argv, const char *simName)
{
	initParameters();
	m_exePath = FileSystem::getProgramPath();
	m_dataPath = FileSystem::normalizePath(getExePath() + "/" + std::string(SPH_DATA_PATH));
	setUseParticleCaching(true);

	try
	{
		cxxopts::Options options(argv[0], "SPlisHSPlasH - An open-source library for the physically-based simulation of fluids.");
		options
			.positional_help("[scene file]")
			.show_positional_help();

		options.add_options()
			("h,help", "Print help")
			("no-cache", "Disable caching of boundary samples/maps.")	
			("state-file", "State file (state_<time>.bin) that should be loaded.", cxxopts::value<std::string>())
			("data-path", "Path of the data directory.", cxxopts::value<std::string>())
			("output-dir", "Output directory for log file and partio files.", cxxopts::value<std::string>())
			("no-initial-pause", "Disable caching of boundary samples/maps.")
			("no-gui", "Disable GUI.")
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

		if (result.count("no-cache"))
		{
			setUseParticleCaching(false);
		}

		if (result.count("no-gui"))
		{
			setUseGUI(false);
		}

		m_dataPath = FileSystem::normalizePath(getExePath() + "/" + std::string(SPH_DATA_PATH));
		if (result.count("data-path"))
		{
			m_dataPath = result["data-path"].as<std::string>();
		}

		m_sceneFile = getDataPath() + "/Scenes/DoubleDamBreak.json";

		if (result.count("scene-file"))
		{
			m_sceneFile = result["scene-file"].as<std::string>();
			if (FileSystem::isRelativePath(m_sceneFile))
				m_sceneFile = FileSystem::normalizePath(m_exePath + "/" + m_sceneFile);
		}
#ifdef WIN32
		else
		{
			std::string scenePath = FileSystem::normalizePath(m_dataPath + "/Scenes/");
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

	LOG_DEBUG << "Git refspec: " << GIT_REFSPEC;
	LOG_DEBUG << "Git SHA1:    " << GIT_SHA1;
	LOG_DEBUG << "Git status:  " << GIT_LOCAL_STATUS;
	LOG_DEBUG << "Host name:   " << SystemInfo::getHostName();

	if (!getUseParticleCaching())
		LOG_INFO << "Boundary cache disabled.";
	LOG_INFO << "Output directory: " << m_outputPath;

	m_sceneLoader = std::unique_ptr<SceneLoader>(new SceneLoader());
	if (m_sceneFile != "")
		m_sceneLoader->readScene(m_sceneFile.c_str(), m_scene);
	else
		return;

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
	FileSystem::copyFile(argv[0], progFilePath + "/" + FileSystem::getFileNameWithExt(argv[0]));
#endif

	// OpenGL
	if (getUseGUI())
	{
		MiniGL::init(argc, argv, 1280, 960, 0, 0, simName);
		MiniGL::initLights();
		MiniGL::getOpenGLVersion(m_context_major_version, m_context_minor_version);
		const bool sim2D = getScene().sim2D;
		if (sim2D)
			MiniGL::setViewport(40.0, 0.1f, 500.0, m_scene.camPosition, m_scene.camLookat);
		else
			MiniGL::setViewport(40.0, 0.1f, 500.0, m_scene.camPosition, m_scene.camLookat);
		MiniGL::setSelectionFunc(selection, this);
		MiniGL::addKeyFunc('i', std::bind(&SimulatorBase::particleInfo, this));
		MiniGL::addKeyFunc('s', std::bind(&SimulatorBase::saveState, this));
#ifdef WIN32
		MiniGL::addKeyFunc('l', std::bind(&SimulatorBase::loadStateDialog, this));
#endif 

		if (MiniGL::checkOpenGLVersion(3, 3))
			initShaders();
	}
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
}

void SimulatorBase::initShaders()
{
	string vertFile = getExePath() + "/resources/shaders/vs_points_vector.glsl";
	string fragFile = getExePath() + "/resources/shaders/fs_points.glsl";
	m_shader_vector.compileShaderFile(GL_VERTEX_SHADER, vertFile);
	m_shader_vector.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
	m_shader_vector.createAndLinkProgram();
	m_shader_vector.begin();
	m_shader_vector.addUniform("modelview_matrix");
	m_shader_vector.addUniform("projection_matrix");
	m_shader_vector.addUniform("radius");
	m_shader_vector.addUniform("viewport_width");
	m_shader_vector.addUniform("color");
	m_shader_vector.addUniform("min_scalar");
	m_shader_vector.addUniform("max_scalar");
	m_shader_vector.end();

	string vertFileScalar = getExePath() + "/resources/shaders/vs_points_scalar.glsl";
	m_shader_scalar.compileShaderFile(GL_VERTEX_SHADER, vertFileScalar);
	m_shader_scalar.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
	m_shader_scalar.createAndLinkProgram();
	m_shader_scalar.begin();
	m_shader_scalar.addUniform("modelview_matrix");
	m_shader_scalar.addUniform("projection_matrix");
	m_shader_scalar.addUniform("radius");
	m_shader_scalar.addUniform("viewport_width");
	m_shader_scalar.addUniform("color");
	m_shader_scalar.addUniform("min_scalar");
	m_shader_scalar.addUniform("max_scalar");
	m_shader_scalar.end();

	string fragFileMap = getExePath() + "/resources/shaders/fs_points_colormap.glsl";
	m_shader_vector_map.compileShaderFile(GL_VERTEX_SHADER, vertFile);
	m_shader_vector_map.compileShaderFile(GL_FRAGMENT_SHADER, fragFileMap);
	m_shader_vector_map.createAndLinkProgram();
	m_shader_vector_map.begin();
	m_shader_vector_map.addUniform("modelview_matrix");
	m_shader_vector_map.addUniform("projection_matrix");
	m_shader_vector_map.addUniform("radius");
	m_shader_vector_map.addUniform("viewport_width");
	m_shader_vector_map.addUniform("color");
	m_shader_vector_map.addUniform("min_scalar");
	m_shader_vector_map.addUniform("max_scalar");
	m_shader_vector_map.end();

	m_shader_scalar_map.compileShaderFile(GL_VERTEX_SHADER, vertFileScalar);
	m_shader_scalar_map.compileShaderFile(GL_FRAGMENT_SHADER, fragFileMap);
	m_shader_scalar_map.createAndLinkProgram();
	m_shader_scalar_map.begin();
	m_shader_scalar_map.addUniform("modelview_matrix");
	m_shader_scalar_map.addUniform("projection_matrix");
	m_shader_scalar_map.addUniform("radius");
	m_shader_scalar_map.addUniform("viewport_width");
	m_shader_scalar_map.addUniform("color");
	m_shader_scalar_map.addUniform("min_scalar");
	m_shader_scalar_map.addUniform("max_scalar");
	m_shader_scalar_map.end();

	vertFile = getExePath() + "/resources/shaders/vs_smooth.glsl";
	fragFile = getExePath() + "/resources/shaders/fs_smooth.glsl";
	m_meshShader.compileShaderFile(GL_VERTEX_SHADER, vertFile);
	m_meshShader.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
	m_meshShader.createAndLinkProgram();
	m_meshShader.begin();
	m_meshShader.addUniform("modelview_matrix");
	m_meshShader.addUniform("projection_matrix");
	m_meshShader.addUniform("surface_color");
	m_meshShader.addUniform("shininess");
	m_meshShader.addUniform("specular_factor");
	m_meshShader.end();

	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &m_textureMap);
}


void SimulatorBase::meshShaderBegin(const float *col)
{
	m_meshShader.begin();
	glUniform1f(m_meshShader.getUniform("shininess"), 5.0f);
	glUniform1f(m_meshShader.getUniform("specular_factor"), 0.2f);

	GLfloat matrix[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
	glUniformMatrix4fv(m_meshShader.getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
	GLfloat pmatrix[16];
	glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
	glUniformMatrix4fv(m_meshShader.getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);
	glUniform3fv(m_meshShader.getUniform("surface_color"), 1, col);
}

void SimulatorBase::meshShaderEnd()
{
	m_meshShader.end();
}

void SimulatorBase::pointShaderBegin(Shader *shader, const float *col, const Real minVal, const Real maxVal, const bool useTexture, float const* color_map)
{
	shader->begin();

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	const Real radius = Simulation::getCurrent()->getValue<Real>(Simulation::PARTICLE_RADIUS);
	glUniform1f(shader->getUniform("viewport_width"), (float)viewport[2]);
	glUniform1f(shader->getUniform("radius"), (float)radius);
	glUniform1f(shader->getUniform("min_scalar"), (GLfloat)minVal);
	glUniform1f(shader->getUniform("max_scalar"), (GLfloat)maxVal);
	glUniform3fv(shader->getUniform("color"), 1, col);

	if (useTexture)
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_1D, m_textureMap);
		glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 256u, 0, GL_RGB, GL_FLOAT, color_map);
		
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glGenerateMipmap(GL_TEXTURE_1D);
	}


	GLfloat matrix[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
	glUniformMatrix4fv(shader->getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
	GLfloat pmatrix[16];
	glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
	glUniformMatrix4fv(shader->getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);

	glEnable(GL_DEPTH_TEST);
	// Point sprites do not have to be explicitly enabled since OpenGL 3.2 where
	// they are enabled by default. Moreover GL_POINT_SPRITE is deprecate and only
	// supported before OpenGL 3.2 or with compatibility profile enabled.
	glEnable(GL_POINT_SPRITE);
	glEnable(GL_PROGRAM_POINT_SIZE);
	glPointParameterf(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);
}

void SimulatorBase::pointShaderEnd(Shader *shader, const bool useTexture)
{
	glBindTexture(GL_TEXTURE_1D, 0);
	shader->end();
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
		m_sceneLoader->readParameterObject(model->getId(), model);
		m_sceneLoader->readParameterObject(key, (ParameterObject*) model->getDragBase());
		m_sceneLoader->readParameterObject(key, (ParameterObject*) model->getSurfaceTensionBase());
		m_sceneLoader->readParameterObject(key, (ParameterObject*) model->getViscosityBase());
		m_sceneLoader->readParameterObject(key, (ParameterObject*) model->getVorticityBase());
		m_sceneLoader->readParameterObject(key, (ParameterObject*) model->getElasticityBase());

		SceneLoader::ColoringData colorData = m_sceneLoader->readColoringInfo(key);
		setColorField(i, colorData.colorField);
		setColorMapType(i, colorData.colorMapType);
		setRenderMinValue(i, colorData.minVal);
		setRenderMaxValue(i, colorData.maxVal);
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
			const string resStr = real2String(m_scene.fluidModels[i]->resolutionSDF[0]) + "_" + real2String(m_scene.fluidModels[i]->resolutionSDF[1]) + "_" + real2String(m_scene.fluidModels[i]->resolutionSDF[2]);
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

		unsigned int maxEmitterParticles = 1000;
		m_sceneLoader->readValue(it->first, "maxEmitterParticles", maxEmitterParticles);
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
			emitterBoundary->scale = Emitter::getSize(ed->width, ed->height, ed->type);
			const Vector3r pos = ed->x;
			emitterBoundary->translation = pos;
			emitterBoundary->samplesFile = "";
			emitterBoundary->mapInvert = false;
			emitterBoundary->mapResolution = Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign>(20, 20, 20);
			emitterBoundary->mapThickness = 0.0;

			if (sim->is2DSimulation())
				emitterBoundary->scale[2] = 2 * supportRadius;

			if (ed->type == 0)
				emitterBoundary->meshFile = FileSystem::normalizePath(getDataPath() + "/models/EmitterBox.obj");
			else if (ed->type == 1)
				emitterBoundary->meshFile = FileSystem::normalizePath(getDataPath() + "/models/EmitterCylinder.obj");
			m_scene.boundaryModels.push_back(emitterBoundary);
			
			// reuse particles if they are outside of a bounding box
			bool emitterReuseParticles = false;
			m_sceneLoader->readValue(model->getId(), "emitterReuseParticles", emitterReuseParticles);

			if (emitterReuseParticles)
			{
				// boxMin
				Vector3r emitterBoxMin(-1.0, -1.0, -1.0);
				m_sceneLoader->readVector(model->getId(), "emitterBoxMin", emitterBoxMin);

				// boxMax
				Vector3r emitterBoxMax(1.0, 1.0, 1.0);
				m_sceneLoader->readVector(model->getId(), "emitterBoxMax", emitterBoxMax);

				model->getEmitterSystem()->enableReuseParticles(emitterBoxMin, emitterBoxMax);
			}
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

void SimulatorBase::renderFluid(const unsigned int fluidModelIndex, float *fluidColor)
{
	// Draw simulation model
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nParticles = model->numActiveParticles();
	if (nParticles == 0)
		return;

	float surfaceColor[4] = { 0.2f, 0.6f, 0.8f, 1 };
	float speccolor[4] = { 1.0, 1.0, 1.0, 1.0 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, speccolor);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);
	glColor3fv(surfaceColor);
	
	const Real supportRadius = sim->getSupportRadius();
	Real vmax = static_cast<Real>(0.4*2.0)*supportRadius / TimeManager::getCurrent()->getTimeStepSize();
	Real vmin = 0.0;

	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		Shader *shader_vector = &m_shader_vector_map;
		Shader *shader_scalar = &m_shader_scalar_map;
		float const *color_map = nullptr;
		if (m_colorMapType[fluidModelIndex] == 1)
			color_map = reinterpret_cast<float const*>(colormap_jet);
		else if (m_colorMapType[fluidModelIndex] == 2)
			color_map = reinterpret_cast<float const*>(colormap_plasma);

		if (m_colorMapType[fluidModelIndex] == 0)
		{
			shader_vector = &m_shader_vector;
			shader_scalar = &m_shader_scalar;
		}

		const FieldDescription *field = nullptr;
		const std::string &colorFieldName = m_colorField[fluidModelIndex];
		field = &model->getField(colorFieldName);


		if (field == nullptr) 
			pointShaderBegin(shader_scalar, &fluidColor[0], m_renderMinValue[fluidModelIndex], m_renderMaxValue[fluidModelIndex],  false);
		else if (field->type == FieldType::Vector3)
			pointShaderBegin(shader_vector, &fluidColor[0], m_renderMinValue[fluidModelIndex], m_renderMaxValue[fluidModelIndex], true, color_map);
		else if (field->type == FieldType::Scalar)
			pointShaderBegin(shader_scalar, &fluidColor[0], m_renderMinValue[fluidModelIndex], m_renderMaxValue[fluidModelIndex], true, color_map);

		if (model->numActiveParticles() > 0)
		{
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &model->getPosition(0));

			if (field != nullptr)
			{
				if (field->type == FieldType::Vector3)
				{
					glEnableVertexAttribArray(1);
					glVertexAttribPointer(1, 3, GL_REAL, GL_FALSE, 0, field->getFct(0));
				}
				else if (field->type == FieldType::Scalar)
				{
					glEnableVertexAttribArray(1);
					glVertexAttribPointer(1, 1, GL_REAL, GL_FALSE, 0, field->getFct(0));
				}
			}	

			glDrawArrays(GL_POINTS, 0, model->numActiveParticles());
			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
		}

		if (field == nullptr)
			pointShaderEnd(shader_scalar, false);
		else if (field->type == FieldType::Vector3)
			pointShaderEnd(shader_vector, true);
		else if (field->type == FieldType::Scalar)
			pointShaderEnd(shader_scalar, true);
	}
	else
	{
		glPointSize(4.0);
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < nParticles; i++)
		{
			Real v = model->getVelocity(i).norm();
			v = static_cast<Real>(0.5)*((v - vmin) / (vmax - vmin));
			v = min(static_cast<Real>(128.0)*v*v, static_cast<Real>(0.5));
			float fluidColor[4] = { 0.2f, 0.2f, 0.2f, 1.0 };
			MiniGL::hsvToRgb(0.55f, 1.0f, 0.5f + (float)v, fluidColor);

			glColor3fv(fluidColor);
			glVertex3v(&model->getPosition(i)[0]);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}


	float red[4] = { 0.8f, 0.0f, 0.0f, 1 };
	const unsigned int fluidIndex = model->getPointSetIndex();
	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		if ((getSelectedParticles().size() > 0) && ((getSelectedParticles()[fluidIndex].size() > 0)))
		{
			pointShaderBegin(&m_shader_vector, &red[0], m_renderMinValue[fluidModelIndex], m_renderMaxValue[fluidModelIndex]);
			const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
			glUniform1f(m_shader_vector.getUniform("radius"), (float)radius*1.05f);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &model->getPosition(0));
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_REAL, GL_FALSE, 0, &model->getVelocity(0));
			glDrawElements(GL_POINTS, (GLsizei) getSelectedParticles()[fluidIndex].size(), GL_UNSIGNED_INT, getSelectedParticles()[fluidIndex].data());
			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
			pointShaderEnd(&m_shader_vector);
		}
	}
	else
	{
		if (getSelectedParticles().size() > 0)
		{
			glPointSize(4.0);
			glDisable(GL_LIGHTING);
			glBegin(GL_POINTS);
			for (unsigned int i = 0; i < getSelectedParticles()[fluidIndex].size(); i++)
			{
				glColor3fv(red);
				glVertex3v(&model->getPosition(getSelectedParticles()[fluidIndex][i])[0]);
			}
			glEnd();
			glEnable(GL_LIGHTING);
		}
	}

}

void SimulatorBase::mouseMove(int x, int y, void *clientData)
{
	SimulatorBase *base = (SimulatorBase*)clientData;
	Simulation *sim = Simulation::getCurrent();

	Vector3r mousePos;
	MiniGL::unproject(x, y, mousePos);
	const Vector3r diff = mousePos - base->m_oldMousePos;

	TimeManager *tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		for (unsigned int j = 0; j < base->m_selectedParticles[i].size(); j++)
		{
			model->getVelocity(base->m_selectedParticles[i][j]) += 5.0*diff / h;
		}
	}
	base->m_oldMousePos = mousePos;
}

void SimulatorBase::selection(const Vector2i &start, const Vector2i &end, void *clientData)
{
	SimulatorBase *base = (SimulatorBase*)clientData;
	Simulation *sim = Simulation::getCurrent();
	base->m_selectedParticles.resize(sim->numberOfFluidModels());
	bool selected = false;
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);

		const unsigned int nParticles = model->numActiveParticles();
 		if (nParticles != 0)
 		{
 			std::vector<unsigned int> hits;
 			base->m_selectedParticles[i].clear();
 			Selection::selectRect(start, end, &model->getPosition(0),
 				&model->getPosition(model->numActiveParticles() - 1),
 				base->m_selectedParticles[i]);
			if (base->m_selectedParticles[i].size() > 0)
				selected = true;
 		}
	}
	if (selected)
		MiniGL::setMouseMoveFunc(GLUT_MIDDLE_BUTTON, mouseMove);
	else
		MiniGL::setMouseMoveFunc(-1, NULL);

	MiniGL::unproject(end[0], end[1], base->m_oldMousePos);
}

void SimulatorBase::particleInfo()
{
	Simulation *sim = Simulation::getCurrent();
	const int maxWidth = 25;
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		if (m_selectedParticles[i].size() > 0)
		{
			LOG_INFO << "---------------------------------------------------------------------------";
			LOG_INFO << model->getId();
			LOG_INFO << "---------------------------------------------------------------------------";
		}
		for (unsigned int j = 0; j < m_selectedParticles[i].size(); j++)
		{
			unsigned int index = m_selectedParticles[i][j];
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
	if (m_enableRigidBodyExport)
		FileSystem::makeDirs(exportPath);

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
			binWriter.write((char) m_scene.boundaryModels[i]->isWall);
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
			const Eigen::Matrix3f RT = R.transpose().template cast<float>();
			binWriter.writeMatrix(RT);
		}
		binWriter.closeFile();
	}

	m_isFirstFrame = false;
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
	if (!outfile.is_open()) {
		LOG_WARN << "Cannot open a file to save a VTK mesh.";
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

void SimulatorBase::step()
{
	if (TimeManager::getCurrent()->getTime() >= m_nextFrameTime)
	{
		m_nextFrameTime += static_cast<Real>(1.0) / m_framesPerSecond;
		if (m_enablePartioExport || m_enableVTKExport)
			particleExport();
		if (m_enableRigidBodyExport)
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

void SimulatorBase::reset()
{
	if (Simulation::getCurrent()->getValue<int>(Simulation::CFL_METHOD) != Simulation::ENUM_CFL_NONE)
		TimeManager::getCurrent()->setTimeStepSize(m_scene.timeStepSize);
	m_nextFrameTime = 0.0;
	m_nextFrameTimeState = 0.0;
	m_frameCounter = 1;
	m_isFirstFrame = true;
#ifdef DL_OUTPUT
	m_nextTiming = 1.0;
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

void SimulatorBase::updateBoundaryForces()
{
	Real h = TimeManager::getCurrent()->getTimeStepSize();
	SceneLoader::Scene &scene = getScene();
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nObjects = sim->numberOfBoundaryModels();
	for (unsigned int i = 0; i < nObjects; i++)
	{
		BoundaryModel *bm = sim->getBoundaryModel(i);
		RigidBodyObject *rbo = bm->getRigidBodyObject();
		if (rbo->isDynamic())
		{
			Vector3r force, torque;
			bm->getForceAndTorque(force, torque);
			rbo->addForce(force);
			rbo->addTorque(torque);
			bm->clearForceAndTorque();
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
	const string resStr = "r" + real2String(resolutionSDF[0]) + "_" + real2String(resolutionSDF[1]) + "_" + real2String(resolutionSDF[2]);
	const string invertStr = "i" + to_string((int)boundaryData->mapInvert);
	const string thicknessStr = "t" + real2String(boundaryData->mapThickness);
	string densityMapFileName = "";
	if (isDynamic)
		densityMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_db_dm_" + real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + ".cdm");
	else 
		densityMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_sb_dm_" + real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + ".cdm");

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

			return -6.0 * supportRadius < dist + cell_diag && dist - cell_diag < 2.0 * supportRadius;
		});
		STOP_TIMING_PRINT;

		// reduction
		if (!no_reduction)
		{
			std::cout << "Reduce discrete fields...";
			densityMap->reduceField(0u, [&](Eigen::Vector3d const &, double v)->double
			{
				return 0.0 <= v && v <= 3.0;
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
	const string resStr = "r" + real2String(resolutionSDF[0]) + "_" + real2String(resolutionSDF[1]) + "_" + real2String(resolutionSDF[2]);
	const string invertStr = "i" + to_string((int)boundaryData->mapInvert);
	const string thicknessStr = "t" + real2String(boundaryData->mapThickness);
	const string kernelStr = "k" + to_string(sim->getKernel());
	string volumeMapFileName = "";
	if (isDynamic)
		volumeMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_db_vm_" + real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + "_" + kernelStr + ".cdm");
	else 
		volumeMapFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_sb_vm_" + real2String(scene.particleRadius) + "_" + scaleStr + "_" + resStr + "_" + invertStr + "_" + thicknessStr + "_" + kernelStr + ".cdm");

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

			return -6.0 * supportRadius < dist + cell_diag && dist - cell_diag < 2.0 * supportRadius;
		});
		STOP_TIMING_PRINT;

		// reduction
		if (!no_reduction)
		{
			std::cout << "Reduce discrete fields...";
			volumeMap->reduceField(0u, [&](Eigen::Vector3d const &, double v)->double
			{
				return 0.0 <= v && v <= 3.0;
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
