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
int SimulatorBase::PARTIO_EXPORT_FPS = -1;
int SimulatorBase::PARTIO_EXPORT_ATTRIBUTES = -1;
int SimulatorBase::RENDER_WALLS = -1;
int SimulatorBase::ENUM_WALLS_NONE = -1;
int SimulatorBase::ENUM_WALLS_PARTICLES_ALL = -1;
int SimulatorBase::ENUM_WALLS_PARTICLES_NO_WALLS = -1;
int SimulatorBase::ENUM_WALLS_GEOMETRY_ALL = -1;
int SimulatorBase::ENUM_WALLS_GEOMETRY_NO_WALLS = -1;
int SimulatorBase::ENUM_RENDER_NONE = -1;
int SimulatorBase::ENUM_RENDER_VELOCITY = -1;
int SimulatorBase::ENUM_RENDER_ANGULAR_VELOCITY = -1;
int SimulatorBase::ENUM_RENDER_DENSITY = -1;

 
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
	m_enablePartioExport = false;
	m_framesPerSecond = 25;
	m_nextFrameTime = 0.0;
	m_frameCounter = 1;
	m_colorField.resize(1, "velocity");
	m_colorMapType.resize(1, 0);
	m_renderMinValue.resize(1, 0.0);
	m_renderMaxValue.resize(1, 5.0);
	m_partioAttributes = "velocity";
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

	PARTIO_EXPORT_FPS = createNumericParameter("partioFPS", "Export FPS", &m_framesPerSecond);
	setGroup(PARTIO_EXPORT_FPS, "Export");
	setDescription(PARTIO_EXPORT_FPS, "Frame rate of partio export.");

	PARTIO_EXPORT_ATTRIBUTES = createStringParameter("partioAttributes", "Export attributes", &m_partioAttributes);
	getParameter(PARTIO_EXPORT_ATTRIBUTES)->setReadOnly(true);
	setGroup(PARTIO_EXPORT_ATTRIBUTES, "Export");
	setDescription(PARTIO_EXPORT_ATTRIBUTES, "Attributes that are exported in the partio files (except id and position).");
}

void SimulatorBase::init(int argc, char **argv, const char *simName)
{
	initParameters();
	m_exePath = FileSystem::getProgramPath();
	m_dataPath = FileSystem::normalizePath(getExePath() + "/" + std::string(SPH_DATA_PATH));

	m_sceneFile = getDataPath() + "/Scenes/DoubleDamBreak.json";
	setUseParticleCaching(true);
	for (int i = 1; i < argc; i++)
	{
		string argStr = argv[i];
		if (argStr == "--no-cache")
			setUseParticleCaching(false);
		else
		{
			m_sceneFile = string(argv[i]);
			if (FileSystem::isRelativePath(m_sceneFile))
				m_sceneFile = FileSystem::normalizePath(m_exePath + "/" + m_sceneFile);
		}
	}

	m_outputPath = FileSystem::normalizePath(getExePath() + "/output/" + FileSystem::getFileName(m_sceneFile));

#ifdef DL_OUTPUT
	std::string sceneFilePath = FileSystem::normalizePath(m_outputPath + "/scene");
	FileSystem::makeDirs(sceneFilePath);
	FileSystem::copyFile(m_sceneFile, sceneFilePath + "/" + FileSystem::getFileNameWithExt(m_sceneFile));

	std::string progFilePath = FileSystem::normalizePath(m_outputPath + "/program");
	FileSystem::makeDirs(progFilePath);
	FileSystem::copyFile(argv[0], progFilePath + "/" + FileSystem::getFileNameWithExt(argv[0]));
#endif

	std::string logPath = FileSystem::normalizePath(m_outputPath + "/log");
	FileSystem::makeDirs(logPath);
	Utilities::logger.addSink(unique_ptr<Utilities::FileSink>(new Utilities::FileSink(Utilities::LogLevel::DEBUG, logPath + "/SPH.log")));

	LOG_DEBUG << "Git refspec: " << GIT_REFSPEC;
	LOG_DEBUG << "Git SHA1:    " << GIT_SHA1;
	LOG_DEBUG << "Git status:  " << GIT_LOCAL_STATUS;
	LOG_DEBUG << "Host name:   " << SystemInfo::getHostName();

	m_sceneLoader = std::unique_ptr<SceneLoader>(new SceneLoader());
	if (m_sceneFile != "")
		m_sceneLoader->readScene(m_sceneFile.c_str(), m_scene);
	else
		return;

	// OpenGL
	MiniGL::init(argc, argv, 1280, 960, 0, 0, simName);
	MiniGL::initLights();
	MiniGL::getOpenGLVersion(m_context_major_version, m_context_minor_version);
	const bool sim2D = getScene().sim2D;
	if (sim2D)
		MiniGL::setViewport(40.0, 0.1f, 500.0, Vector3r(0.0, 0.0, 8.0), Vector3r(0.0, 0.0, 0.0));
	else
		MiniGL::setViewport(40.0, 0.1f, 500.0, Vector3r(0.0, 3.0, 8.0), Vector3r(0.0, 0.0, 0.0));
	MiniGL::setSelectionFunc(selection, this);
	MiniGL::addKeyFunc('i', std::bind(&SimulatorBase::particleInfo, this));

	if (MiniGL::checkOpenGLVersion(3, 3))
		initShaders();
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

		PartioReaderWriter::readParticles(fileName, m_scene.fluidModels[i]->translation, m_scene.fluidModels[i]->rotation, m_scene.fluidModels[i]->scale, fluidParticles[fluidIndex], fluidVelocities[fluidIndex]);
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
			model->getEmitterSystem()->addEmitter(ed->width, ed->height,
				ed->x, ed->dir, ed->v, ed->emitsPerSecond, ed->type);

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
		}
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
						if (j % 2)
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
		pointShaderBegin(&m_shader_scalar, &red[0], m_renderMinValue[fluidModelIndex], m_renderMaxValue[fluidModelIndex]);
		if ((getSelectedParticles().size() > 0) && ((getSelectedParticles()[fluidIndex].size() > 0)))
		{
			const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
			glUniform1f(m_shader_scalar.getUniform("radius"), (float)radius*1.05f);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &model->getPosition(0));
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_REAL, GL_FALSE, 0, &model->getVelocity(0));
			glDrawElements(GL_POINTS, (GLsizei) getSelectedParticles()[fluidIndex].size(), GL_UNSIGNED_INT, getSelectedParticles()[fluidIndex].data());
			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
		}
		pointShaderEnd(&m_shader_scalar);
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

void SimulatorBase::selection(const Eigen::Vector2i &start, const Eigen::Vector2i &end, void *clientData)
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
	
	if (selected && (MiniGL::getModifierKey() == 3))
	{

	}
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
 					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << field.name + ":" << *field.getFct(index);
				else if (field.type == Vector3)
				{
					Eigen::Map<Vector3r> vec(field.getFct(index));
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << field.name + ":" << vec.transpose();
				}
				else if (field.type == Vector6)
				{
					Eigen::Map<Vector6r> vec(field.getFct(index));
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << field.name + ":" << vec.transpose();
				}
				else if (field.type == Matrix3)
				{
					Eigen::Map<Matrix3r> mat(field.getFct(index));
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << field.name + ":" << mat.row(0);
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << " " << mat.row(1);
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << " " << mat.row(2);
				}
				else if (field.type == Matrix6)
				{
					Eigen::Map<Matrix6r> mat(field.getFct(index));
					LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << field.name + ":" << mat.row(0);
					for (unsigned int k = 1; k < 6; k++)
						LOG_INFO << std::left << std::setw(maxWidth) << std::setfill(' ') << " " << mat.row(k);
				}
 			}
			LOG_INFO << "---------------------------------------------------------------------------\n";
		}
	}
}

void SimulatorBase::partioExport()
{	
	std::string exportPath = FileSystem::normalizePath(m_outputPath + "/partio");
	FileSystem::makeDirs(exportPath);

	Simulation *sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		std::string fileName = "ParticleData";
		fileName = fileName + "_" + model->getId() + "_" + std::to_string(m_frameCounter) + ".bgeo";
		std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

		writeParticles(exportFileName, model);
	}
}


void SimulatorBase::writeParticles(const std::string &fileName, FluidModel *model)
{
	Partio::ParticlesDataMutable& particleData = *Partio::create();
	Partio::ParticleAttribute posAttr = particleData.addAttribute("position", Partio::VECTOR, 3);
	Partio::ParticleAttribute idAttr = particleData.addAttribute("id", Partio::INT, 1);

	// add attributes
	std::vector<std::string> attributes;
	StringTools::tokenize(m_partioAttributes, attributes, ";");

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
	
		id[0] = i;

 		for (unsigned int j = 0; j < attributes.size(); j++)
 		{
 			const int fieldIndex = attrMap[j];
 			if (fieldIndex != -1)
 			{
 				const FieldDescription &field = model->getField(fieldIndex);
 				if (field.type == FieldType::Scalar)
 				{
 					float* val = particleData.dataWrite<float>(partioAttrMap[j], index);
 					*val = (float) *field.getFct(i);
 				}
 				else if (field.type == FieldType::Vector3)
 				{
 					float* val = particleData.dataWrite<float>(partioAttrMap[j], index);
					Eigen::Map<Vector3r> vec(field.getFct(i));
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

void SimulatorBase::step()
{
	if (m_enablePartioExport)
	{
		if (TimeManager::getCurrent()->getTime() >= m_nextFrameTime)
		{
			m_nextFrameTime += static_cast<Real>(1.0) / m_framesPerSecond;
			partioExport();
			m_frameCounter++;
		}
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
	//TimeManager::getCurrent()->setTimeStepSize(m_scene.timeStepSize);
	m_nextFrameTime = 0.0;
	m_frameCounter = 1;
#ifdef DL_OUTPUT
	m_nextTiming = 1.0;
#endif
}

