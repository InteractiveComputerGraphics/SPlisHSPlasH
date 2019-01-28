#include "DemoBase.h"
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

INIT_LOGGING
INIT_TIMING
INIT_COUNTING

using namespace SPH;
using namespace std;
using namespace GenParam;
using namespace Utilities;

int DemoBase::PAUSE = -1;
int DemoBase::PAUSE_AT = -1;
int DemoBase::STOP_AT = -1;
int DemoBase::NUM_STEPS_PER_RENDER = -1;
int DemoBase::PARTIO_EXPORT = -1;
int DemoBase::PARTIO_EXPORT_FPS = -1;
int DemoBase::RENDER_MIN_VALUE = -1;
int DemoBase::RENDER_MAX_VALUE = -1;
int DemoBase::RENDER_WALLS = -1;
int DemoBase::RENDER_COLOR_FIELD = -1;
int DemoBase::RENDER_COLOR_MAP_TYPE = -1;
int DemoBase::ENUM_WALLS_NONE = -1;
int DemoBase::ENUM_WALLS_PARTICLES_ALL = -1;
int DemoBase::ENUM_WALLS_PARTICLES_NO_WALLS = -1;
int DemoBase::ENUM_WALLS_GEOMETRY_ALL = -1;
int DemoBase::ENUM_WALLS_GEOMETRY_NO_WALLS = -1;
int DemoBase::ENUM_RENDER_NONE = -1;
int DemoBase::ENUM_RENDER_VELOCITY = -1;
int DemoBase::ENUM_RENDER_ANGULAR_VELOCITY = -1;
int DemoBase::ENUM_RENDER_DENSITY = -1;
int DemoBase::ENUM_COLORMAP_NONE = -1;
int DemoBase::ENUM_COLORMAP_JET = -1;
int DemoBase::ENUM_COLORMAP_PLASMA = -1;

 
DemoBase::DemoBase()
{
	Utilities::logger.addSink(unique_ptr<Utilities::ConsoleSink>(new Utilities::ConsoleSink(Utilities::LogLevel::INFO)));

	m_numberOfStepsPerRenderUpdate = 4;
	m_colorField = 1;
	setColorMapType(0);
	m_renderMinValue = 0.0;
	m_renderMaxValue = 5.0;
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
#ifdef DL_OUTPUT
	m_nextTiming = 1.0;
#endif
}

DemoBase::~DemoBase()
{

}

void DemoBase::initParameters()
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

	RENDER_MIN_VALUE = createNumericParameter("renderMinValue", "Min. value (shader)", &m_renderMinValue);
	setGroup(RENDER_MIN_VALUE, "Visualization");
	setDescription(RENDER_MIN_VALUE, "Minimal value used for color-coding the color field in the rendering process.");

	RENDER_MAX_VALUE = createNumericParameter("renderMaxValue", "Max. value (shader)", &m_renderMaxValue);
	setGroup(RENDER_MAX_VALUE, "Visualization");
	setDescription(RENDER_MAX_VALUE, "Maximal value used for color-coding the color field in the rendering process.");

	RENDER_COLOR_FIELD = createEnumParameter("colorField", "Color field", &m_colorField);
	setGroup(RENDER_COLOR_FIELD, "Visualization");
	setDescription(RENDER_COLOR_FIELD, "Choose vector or scalar field for particle coloring.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(RENDER_COLOR_FIELD));
	enumParam->addEnumValue("None", ENUM_RENDER_NONE);
	enumParam->addEnumValue("Velocity", ENUM_RENDER_VELOCITY);
	enumParam->addEnumValue("Angular velocity (micropolar model)", ENUM_RENDER_ANGULAR_VELOCITY);
	enumParam->addEnumValue("Density", ENUM_RENDER_DENSITY);

	ParameterBase::GetFunc<int> getColorMapTypeFct = std::bind(&DemoBase::getColorMapType, this);
	ParameterBase::SetFunc<int> setColorMapTypeFct = std::bind(&DemoBase::setColorMapType, this, std::placeholders::_1);
	RENDER_COLOR_MAP_TYPE = createEnumParameter("colorMapType", "Color map", getColorMapTypeFct, setColorMapTypeFct);
	setGroup(RENDER_COLOR_MAP_TYPE, "Visualization");
	setDescription(RENDER_COLOR_MAP_TYPE, "Choose a color map.");
	enumParam = static_cast<EnumParameter*>(getParameter(RENDER_COLOR_MAP_TYPE));
	enumParam->addEnumValue("None", ENUM_COLORMAP_NONE);
	enumParam->addEnumValue("Jet", ENUM_COLORMAP_JET);
	enumParam->addEnumValue("Plasma", ENUM_COLORMAP_PLASMA);

	RENDER_WALLS = createEnumParameter("renderWalls", "Render walls", &m_renderWalls);
	setGroup(RENDER_WALLS, "Visualization");
	setDescription(RENDER_WALLS, "Make walls visible/invisible.");
	enumParam = static_cast<EnumParameter*>(getParameter(RENDER_WALLS));
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
}

void DemoBase::init(int argc, char **argv, const char *demoName)
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
	MiniGL::init(argc, argv, 1280, 960, 0, 0, demoName);
	MiniGL::initLights();
	MiniGL::getOpenGLVersion(m_context_major_version, m_context_minor_version);
	MiniGL::setViewport(40.0, 0.1f, 500.0, Vector3r(0.0, 3.0, 8.0), Vector3r(0.0, 0.0, 0.0));
	MiniGL::setSelectionFunc(selection, this);
	MiniGL::addKeyFunc('i', std::bind(&DemoBase::particleInfo, this));

	if (MiniGL::checkOpenGLVersion(3, 3))
		initShaders();
}

void DemoBase::cleanup()
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

void DemoBase::initShaders()
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


void DemoBase::meshShaderBegin(const float *col)
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

void DemoBase::meshShaderEnd()
{
	m_meshShader.end();
}

void DemoBase::pointShaderBegin(Shader *shader, const float *col, const bool useTexture)
{
	shader->begin();

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	const Real radius = Simulation::getCurrent()->getValue<Real>(Simulation::PARTICLE_RADIUS);
	glUniform1f(shader->getUniform("viewport_width"), (float)viewport[2]);
	glUniform1f(shader->getUniform("radius"), (float)radius);
	glUniform1f(shader->getUniform("min_scalar"), (GLfloat)m_renderMinValue);
	glUniform1f(shader->getUniform("max_scalar"), (GLfloat)m_renderMaxValue);
	glUniform3fv(shader->getUniform("color"), 1, col);

	if (useTexture)
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_1D, m_textureMap);
		glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, m_colorMapLength, 0, GL_RGB, GL_FLOAT,
			reinterpret_cast<float const*>(m_colorMapBuffer));
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

void DemoBase::pointShaderEnd(Shader *shader, const bool useTexture)
{
	glBindTexture(GL_TEXTURE_1D, 0);
	shader->end();
}

void DemoBase::readParameters()
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
	}
}


void DemoBase::buildModel()
{
	TimeManager::getCurrent()->setTimeStepSize(m_scene.timeStepSize);

	initFluidData();

	createEmitters();

	Simulation *sim = Simulation::getCurrent();

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


void DemoBase::initFluidData()
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

	LOG_INFO << "Number of fluid particles: " << nParticles;
}

void DemoBase::createEmitters()
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


void DemoBase::createFluidBlocks(std::map<std::string, unsigned int> &fluidIDs, std::vector<std::vector<Vector3r>> &fluidParticles, std::vector<std::vector<Vector3r>> &fluidVelocities)
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

void DemoBase::renderFluid(FluidModel *model, float *fluidColor)
{
	// Draw simulation model
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

	Simulation *sim = Simulation::getCurrent();
	const Real supportRadius = sim->getSupportRadius();
	Real vmax = static_cast<Real>(0.4*2.0)*supportRadius / TimeManager::getCurrent()->getTimeStepSize();
	Real vmin = 0.0;

	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		Shader *shader_vector = &m_shader_vector_map;
		Shader *shader_scalar = &m_shader_scalar_map;
		if (m_colorMapType == ENUM_COLORMAP_NONE)
		{
			shader_vector = &m_shader_vector;
			shader_scalar = &m_shader_scalar;
		}

		if ((m_colorField == ENUM_RENDER_VELOCITY) || (m_colorField == ENUM_RENDER_ANGULAR_VELOCITY))
			pointShaderBegin(shader_vector, &fluidColor[0], true);
		else if (m_colorField == ENUM_RENDER_DENSITY)
			pointShaderBegin(shader_scalar, &fluidColor[0], true);
		else 
			pointShaderBegin(shader_scalar, &fluidColor[0], false);

		if (model->numActiveParticles() > 0)
		{
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &model->getPosition(0));
			
			if (m_colorField == ENUM_RENDER_VELOCITY)
			{
				glEnableVertexAttribArray(1);
				glVertexAttribPointer(1, 3, GL_REAL, GL_FALSE, 0, &model->getVelocity(0));
			}
			else if ((m_colorField == ENUM_RENDER_ANGULAR_VELOCITY) && ((VorticityMethods)model->getVorticityMethod() == VorticityMethods::Micropolar))
			{
				glEnableVertexAttribArray(1);
				float fluidColor2[4] = { 0.3f, 0.9f, 0.5f, 1.0f };
				glUniform3fv(shader_vector->getUniform("color"), 1, fluidColor2);
				glVertexAttribPointer(1, 3, GL_REAL, GL_FALSE, 0, &((MicropolarModel_Bender2017*)model->getVorticityBase())->getAngularVelocity(0)[0]);
			}
			else if (m_colorField == ENUM_RENDER_DENSITY)
			{
				glEnableVertexAttribArray(1);
				float fluidColor2[4] = { 0.9f, 0.3f, 0.3f, 1.0f };
				glUniform3fv(shader_scalar->getUniform("color"), 1, fluidColor2);
				glVertexAttribPointer(1, 1, GL_REAL, GL_FALSE, 0, &model->getDensity(0));
			}				

			glDrawArrays(GL_POINTS, 0, model->numActiveParticles());
			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
		}

		if ((m_colorField == ENUM_RENDER_VELOCITY) || (m_colorField == ENUM_RENDER_ANGULAR_VELOCITY))
			pointShaderEnd(shader_vector, true);
		else if (m_colorField == ENUM_RENDER_DENSITY)
			pointShaderEnd(shader_scalar, true);
		else
			pointShaderEnd(shader_scalar, false);
		
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
		pointShaderBegin(&m_shader_scalar, &red[0]);
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

void DemoBase::mouseMove(int x, int y, void *clientData)
{
	DemoBase *base = (DemoBase*)clientData;
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

void DemoBase::selection(const Eigen::Vector2i &start, const Eigen::Vector2i &end, void *clientData)
{
	DemoBase *base = (DemoBase*)clientData;
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

void DemoBase::particleInfo()
{
	Simulation *sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		for (unsigned int j = 0; j < m_selectedParticles[i].size(); j++)
		{
			unsigned int index = m_selectedParticles[i][j];
			LOG_INFO << index;
			LOG_INFO << "x:       " << model->getPosition(index).transpose();
			LOG_INFO << "v:       " << model->getVelocity(index).transpose();
			LOG_INFO << "density: " << model->getDensity(index);
		}
	}
}

void DemoBase::partioExport()
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

		PartioReaderWriter::writeParticles(exportFileName, model->numActiveParticles(), &model->getPosition(0), &model->getVelocity(0), 0.0);
	}
}

void DemoBase::step()
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

void DemoBase::reset()
{
	//TimeManager::getCurrent()->setTimeStepSize(m_scene.timeStepSize);
	m_nextFrameTime = 0.0;
	m_frameCounter = 1;
#ifdef DL_OUTPUT
	m_nextTiming = 1.0;
#endif
}

void DemoBase::setColorMapType(const int v)
{
	m_colorMapType = v;
	if (m_colorMapType == ENUM_COLORMAP_JET)
	{
		m_colorMapBuffer = colormap_jet[0];
		m_colorMapLength = 256u;
	}
	else if (m_colorMapType == ENUM_COLORMAP_PLASMA)
	{
		m_colorMapBuffer = colormap_plasma[0];
		m_colorMapLength = 256u;
	}
}
