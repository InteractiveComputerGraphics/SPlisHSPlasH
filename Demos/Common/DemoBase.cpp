#include "DemoBase.h"
#include "Visualization/MiniGL.h"
#include "SPlisHSPlasH/Utilities/SceneLoader.h"
#include "Utilities/FileSystem.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/PartioReaderWriter.h"
#include "Visualization/Selection.h"
#include "GL/glut.h"
#include "SPlisHSPlasH/Emitter.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Vorticity/MicropolarModel_Bender2017.h"
#include "NumericParameter.h"
#include "Utilities/Logger.h"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "Utilities/Version.h"
#include "Utilities/SystemInfo.h"

INIT_LOGGING
INIT_TIMING
INIT_COUNTING

using namespace SPH;
using namespace std;
using namespace GenParam;
using namespace Utilities;

int DemoBase::PAUSE = -1;
int DemoBase::PAUSE_AT = -1;
int DemoBase::NUM_STEPS_PER_RENDER = -1;
int DemoBase::PARTIO_EXPORT = -1;
int DemoBase::PARTIO_EXPORT_FPS = -1;
int DemoBase::RENDER_OMEGA = -1;
int DemoBase::RENDER_MAX_VEL = -1;
int DemoBase::RENDER_WALLS = -1;
int DemoBase::ENUM_WALLS_NONE = -1;
int DemoBase::ENUM_WALLS_PARTICLES_ALL = -1;
int DemoBase::ENUM_WALLS_PARTICLES_NO_WALLS = -1;
int DemoBase::ENUM_WALLS_GEOMETRY_ALL = -1;
int DemoBase::ENUM_WALLS_GEOMETRY_NO_WALLS = -1;

 
DemoBase::DemoBase()
{
	Utilities::logger.addSink(unique_ptr<Utilities::ConsoleSink>(new Utilities::ConsoleSink(Utilities::LogLevel::INFO)));

	m_numberOfStepsPerRenderUpdate = 4;
	m_renderAngularVelocities = false;
	m_renderMaxVelocity = 25.0;
	m_sceneFile = "";
	m_renderWalls = 4;
	m_doPause = true;
	m_pauseAt = -1.0;
	m_useParticleCaching = true;
	m_enablePartioExport = false;
	m_framesPerSecond = 25;
	m_nextFrameTime = 0.0;
	m_frameCounter = 1;
}

DemoBase::~DemoBase()
{

}

void DemoBase::initParameters()
{
	ParameterObject::initParameters();

	PAUSE = createBoolParameter("pause", "Pause", &m_doPause);
	setGroup(PAUSE, "Simulation");
	setDescription(PAUSE, "Pause simulation.");
	setHotKey(PAUSE, "space");

	PAUSE_AT = createNumericParameter("pauseAt", "Pause simulation at", &m_pauseAt);
	setGroup(PAUSE_AT, "Simulation");
	setDescription(PAUSE_AT, "Pause simulation at the given time. When the value is negative, the simulation is not paused.");

	PARTIO_EXPORT = createBoolParameter("enablePartioExport", "Partio export", &m_enablePartioExport);
	setGroup(PARTIO_EXPORT, "Export");
	setDescription(PARTIO_EXPORT, "Enable/disable partio export.");

	PARTIO_EXPORT_FPS = createNumericParameter("partioFPS", "Export FPS", &m_framesPerSecond);
	setGroup(PARTIO_EXPORT_FPS, "Export");
	setDescription(PARTIO_EXPORT_FPS, "Frame rate of partio export.");

	NUM_STEPS_PER_RENDER = createNumericParameter("numberOfStepsPerRenderUpdate", "# time steps / update", &m_numberOfStepsPerRenderUpdate);
	setGroup(NUM_STEPS_PER_RENDER, "Visualization");
	setDescription(NUM_STEPS_PER_RENDER, "Pause simulation at the given time. When the value is negative, the simulation is not paused.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(NUM_STEPS_PER_RENDER))->setMinValue(1);

	RENDER_MAX_VEL = createNumericParameter("renderMaxVelocity", "Max. velocity (shader)", &m_renderMaxVelocity);
	setGroup(RENDER_MAX_VEL, "Visualization");
	setDescription(RENDER_MAX_VEL, "Maximal velocity used for color-coding the velocity in the rendering process.");
	static_cast<RealParameter*>(getParameter(RENDER_MAX_VEL))->setMinValue(0.00001);

	RENDER_OMEGA = createBoolParameter("renderAngularVelocities", "Render angular velocities", &m_renderAngularVelocities);
	setGroup(RENDER_OMEGA, "Visualization");
	setDescription(RENDER_OMEGA, "Render angular velocities (micropolar model).");

	RENDER_WALLS = createEnumParameter("renderWalls", "Render walls", &m_renderWalls);
	setGroup(RENDER_WALLS, "Visualization");
	setDescription(RENDER_WALLS, "Make walls visible/invisible.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(RENDER_WALLS));
	enumParam->addEnumValue("None", ENUM_WALLS_NONE);
	enumParam->addEnumValue("Particles (all)", ENUM_WALLS_PARTICLES_ALL);
	enumParam->addEnumValue("Particles (no walls)", ENUM_WALLS_PARTICLES_NO_WALLS);
	enumParam->addEnumValue("Geometry (all)", ENUM_WALLS_GEOMETRY_ALL);
	enumParam->addEnumValue("Geometry (no walls)", ENUM_WALLS_GEOMETRY_NO_WALLS);
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
}

void DemoBase::initShaders()
{
	string vertFile = getExePath() + "/resources/shaders/vs_points.glsl";
	string fragFile = getExePath() + "/resources/shaders/fs_points.glsl";
	m_shader.compileShaderFile(GL_VERTEX_SHADER, vertFile);
	m_shader.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
	m_shader.createAndLinkProgram();
	m_shader.begin();
	m_shader.addUniform("modelview_matrix");
	m_shader.addUniform("projection_matrix");
	m_shader.addUniform("radius");
	m_shader.addUniform("viewport_width");
	m_shader.addUniform("color");
	m_shader.addUniform("projection_radius");
	m_shader.addUniform("max_velocity");
	m_shader.end();

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

void DemoBase::pointShaderBegin(const float *col)
{
	m_shader.begin();

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	const Real radius = Simulation::getCurrent()->getModel()->getValue<Real>(FluidModel::PARTICLE_RADIUS);
	glUniform1f(m_shader.getUniform("viewport_width"), (float)viewport[2]);
	glUniform1f(m_shader.getUniform("radius"), (float)radius);
	glUniform1f(m_shader.getUniform("max_velocity"), (GLfloat) m_renderMaxVelocity);
	glUniform3fv(m_shader.getUniform("color"), 1, col);

	GLfloat matrix[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
	glUniformMatrix4fv(m_shader.getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
	GLfloat pmatrix[16];
	glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
	glUniformMatrix4fv(m_shader.getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);

	glEnable(GL_DEPTH_TEST);
	// Point sprites do not have to be explicitly enabled since OpenGL 3.2 where
	// they are enabled by default. Moreover GL_POINT_SPRITE is deprecate and only
	// supported before OpenGL 3.2 or with compatibility profile enabled.
	glEnable(GL_POINT_SPRITE);
	glEnable(GL_PROGRAM_POINT_SIZE);
	glPointParameterf(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);
}

void DemoBase::pointShaderEnd()
{
	m_shader.end();
}

void DemoBase::readParameters()
{
	m_sceneLoader->readParameterObject(this);
	m_sceneLoader->readParameterObject(Simulation::getCurrent());
	m_sceneLoader->readParameterObject(Simulation::getCurrent()->getModel());
	m_sceneLoader->readParameterObject(Simulation::getCurrent()->getTimeStep());
	m_sceneLoader->readParameterObject(Simulation::getCurrent()->getDragBase());
	m_sceneLoader->readParameterObject(Simulation::getCurrent()->getSurfaceTensionBase());
	m_sceneLoader->readParameterObject(Simulation::getCurrent()->getViscosityBase());
	m_sceneLoader->readParameterObject(Simulation::getCurrent()->getVorticityBase());
}


void DemoBase::buildModel()
{
	TimeManager::getCurrent()->setTimeStepSize(m_scene.timeStepSize);

	std::vector<Vector3r> fluidParticles;
	std::vector<Vector3r> fluidVelocities;
	initFluidData(fluidParticles, fluidVelocities);

	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getModel();
	model->setValue(FluidModel::PARTICLE_RADIUS, m_scene.particleRadius);

	model->initModel((unsigned int)fluidParticles.size(), fluidParticles.data(), fluidVelocities.data(), m_scene.maxEmitterParticles);
	sim->getTimeStep()->resize();

	LOG_INFO << "Number of fluid particles: " << fluidParticles.size();

	unsigned int nBoundaryParticles = 0;
	for (unsigned int i = 0; i < model->numberOfRigidBodyParticleObjects(); i++)
		nBoundaryParticles += model->getRigidBodyParticleObject(i)->numberOfParticles();

	LOG_INFO << "Number of boundary particles: " << nBoundaryParticles;

	//////////////////////////////////////////////////////////////////////////
	// emitters
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int i = 0; i < m_scene.emitters.size(); i++)
	{
		SceneLoader::EmitterData *ed = m_scene.emitters[i];
		model->getEmitterSystem()->addEmitter(ed->width, ed->height,
			ed->x, ed->dir, ed->v, ed->emitsPerSecond, ed->type);
	}
	if (m_scene.emitterReuseParticles)
		model->getEmitterSystem()->enableReuseParticles(m_scene.emitterBoxMin, m_scene.emitterBoxMax);

	model->setValue(FluidModel::KERNEL_METHOD, FluidModel::ENUM_KERNEL_PRECOMPUTED_CUBIC);
	model->setValue(FluidModel::GRAD_KERNEL_METHOD, FluidModel::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
}


void DemoBase::initFluidData(std::vector<Vector3r> &fluidParticles, std::vector<Vector3r> &fluidVelocities)
{
	LOG_INFO << "Initialize fluid particles";
	createFluidBlocks(fluidParticles, fluidVelocities);

	std::string base_path = FileSystem::getFilePath(m_sceneFile);

	unsigned int startIndex = 0;
	unsigned int endIndex = 0;
	for (unsigned int i = 0; i < m_scene.fluidModels.size(); i++)
	{
		std::string fileName;
		if (FileSystem::isRelativePath(m_scene.fluidModels[i]->samplesFile))
			fileName = base_path + "/" + m_scene.fluidModels[i]->samplesFile;
		else
			fileName = m_scene.fluidModels[i]->samplesFile;

		PartioReaderWriter::readParticles(fileName, m_scene.fluidModels[i]->translation, m_scene.fluidModels[i]->rotation, m_scene.fluidModels[i]->scale, fluidParticles, fluidVelocities);
		Simulation::getCurrent()->getModel()->setValue(FluidModel::PARTICLE_RADIUS, m_scene.particleRadius);
	}
}


void DemoBase::createFluidBlocks(std::vector<Vector3r> &fluidParticles, std::vector<Vector3r> &fluidVelocities)
{
	for (unsigned int i = 0; i < m_scene.fluidBlocks.size(); i++)
	{
		const Real diam = 2.0*m_scene.particleRadius;

		Real xshift = diam;
		Real yshift = diam;
		const Real eps = 1.0e-9;
		if (m_scene.fluidBlocks[i]->mode == 1)
			yshift = sqrt(3.0) * m_scene.particleRadius + eps;
		else if (m_scene.fluidBlocks[i]->mode == 2)
		{
			xshift = sqrt(6.0) * diam / 3.0 + eps;
			yshift = sqrt(3.0) * m_scene.particleRadius + eps;
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
		const int stepsZ = (int)round(diff[2] / diam) - 1;

		Vector3r start = m_scene.fluidBlocks[i]->box.m_minX + 2.0*m_scene.particleRadius*Vector3r::Ones();
		fluidParticles.reserve(fluidParticles.size() + stepsX*stepsY*stepsZ);
		fluidVelocities.resize(fluidParticles.size() + stepsX*stepsY*stepsZ, m_scene.fluidBlocks[i]->initialVelocity);
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
							shift_vec[2] += diam / (2.0 * (k % 2 ? -1 : 1));
						}
						if (k % 2 == 0)
						{
							shift_vec[0] += xshift / 2.0;
						}
						currPos += shift_vec;
					}
					fluidParticles.push_back(currPos);
				}
			}
		}
	}
}

void DemoBase::renderFluid()
{
	// Draw simulation model
	MiniGL::drawTime(TimeManager::getCurrent()->getTime());
	FluidModel *model = Simulation::getCurrent()->getModel();
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


	const Real supportRadius = model->getSupportRadius();
	Real vmax = 0.4*2.0*supportRadius / TimeManager::getCurrent()->getTimeStepSize();
	Real vmin = 0.0;

	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		float fluidColor[4] = { 0.3f, 0.5f, 0.9f, 1.0f };
		float fluidColor2[4] = { 0.3f, 0.9f, 0.5f, 1.0f };
		pointShaderBegin(&fluidColor[0]);

		if (model->numActiveParticles() > 0)
		{
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, &model->getPosition(0, 0));
			glEnableVertexAttribArray(1);
			if (m_renderAngularVelocities && ((VorticityMethods)Simulation::getCurrent()->getVorticityMethod() == VorticityMethods::Micropolar))
			{
				glUniform3fv(m_shader.getUniform("color"), 1, fluidColor2);
				glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 0, &((MicropolarModel_Bender2017*)Simulation::getCurrent()->getVorticityBase())->getAngularVelocity(0)[0]);
			}
			else
				glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 0, &model->getVelocity(0, 0));

			glDrawArrays(GL_POINTS, 0, model->numActiveParticles());
			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
		}

		pointShaderEnd();
	}
	else
	{
		glPointSize(4.0);
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < nParticles; i++)
		{
			Real v = model->getVelocity(0, i).norm();
			v = 0.5*((v - vmin) / (vmax - vmin));
			v = min(128.0*v*v, 0.5);
			float fluidColor[4] = { 0.2f, 0.2f, 0.2f, 1.0 };
			MiniGL::hsvToRgb(0.55f, 1.0f, 0.5f + (float)v, fluidColor);

			glColor3fv(fluidColor);
			glVertex3v(&model->getPosition(0, i)[0]);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}


	float red[4] = { 0.8f, 0.0f, 0.0f, 1 };
	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		pointShaderBegin(&red[0]);
		if (getSelectedParticles().size() > 0)
		{
			const Real radius = model->getValue<Real>(FluidModel::PARTICLE_RADIUS);
			glUniform1f(m_shader.getUniform("radius"), (float)radius*1.05f);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, &model->getPosition(0, 0));
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 0, &model->getVelocity(0, 0));
			glDrawElements(GL_POINTS, (GLsizei) getSelectedParticles().size(), GL_UNSIGNED_INT, getSelectedParticles().data());
			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
		}
		pointShaderEnd();
	}
	else
	{
		glPointSize(4.0);
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < getSelectedParticles().size(); i++)
		{			
			glColor3fv(red);
			glVertex3v(&model->getPosition(0, getSelectedParticles()[i])[0]);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}

}

void DemoBase::mouseMove(int x, int y, void *clientData)
{
	DemoBase *base = (DemoBase*)clientData;
	FluidModel *model = Simulation::getCurrent()->getModel();

	Vector3r mousePos;
	MiniGL::unproject(x, y, mousePos);
	const Vector3r diff = mousePos - base->m_oldMousePos;

	TimeManager *tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	for (unsigned int j = 0; j < base->m_selectedParticles.size(); j++)
	{
		model->getVelocity(0, base->m_selectedParticles[j]) += 5.0*diff / h;
	}
	base->m_oldMousePos = mousePos;
}

void DemoBase::selection(const Eigen::Vector2i &start, const Eigen::Vector2i &end, void *clientData)
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	DemoBase *base = (DemoBase*)clientData;
	const unsigned int nParticles = model->numActiveParticles();
	if (nParticles == 0)
		return;

	std::vector<unsigned int> hits;
	base->m_selectedParticles.clear();
	Selection::selectRect(start, end, &model->getPosition(0, 0),
		&model->getPosition(0, model->numActiveParticles() - 1),
		base->m_selectedParticles);
	if (base->m_selectedParticles.size() > 0)
		MiniGL::setMouseMoveFunc(GLUT_MIDDLE_BUTTON, mouseMove);
	else
		MiniGL::setMouseMoveFunc(-1, NULL);

	MiniGL::unproject(end[0], end[1], base->m_oldMousePos);
}

void DemoBase::partioExport()
{
	FluidModel *model = Simulation::getCurrent()->getModel();
	std::string exportPath = FileSystem::normalizePath(m_outputPath + "/partio");
	FileSystem::makeDirs(exportPath);

	std::string fileName = "ParticleData";
	fileName = fileName + std::to_string(m_frameCounter) + ".bgeo";
	std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

	PartioReaderWriter::writeParticles(exportFileName, model->numActiveParticles(), &model->getPosition(0, 0), &model->getVelocity(0, 0), 0.0);
}

void DemoBase::step()
{
	if (m_enablePartioExport)
	{
		if (TimeManager::getCurrent()->getTime() >= m_nextFrameTime)
		{
			m_nextFrameTime += 1.0 / m_framesPerSecond;
			partioExport();
			m_frameCounter++;
		}
	}
}

void DemoBase::reset()
{
	m_nextFrameTime = 0.0;
	m_frameCounter = 1;
}