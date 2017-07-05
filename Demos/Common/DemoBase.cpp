#include "DemoBase.h"
#include "Visualization/MiniGL.h"
#include "Utilities/SceneLoader.h"
#include "Utilities/FileSystem.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/WCSPH/TimeStepWCSPH.h"
#include "SPlisHSPlasH/PCISPH/TimeStepPCISPH.h"
#include "SPlisHSPlasH/PBF/TimeStepPBF.h"
#include "SPlisHSPlasH/IISPH/TimeStepIISPH.h"
#include "SPlisHSPlasH/DFSPH/TimeStepDFSPH.h"
#include "SPlisHSPlasH/PF/TimeStepPF.h"
#include "Utilities/PartioReaderWriter.h"
#include "Visualization/Selection.h"
#include "GL/glut.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_Bender2017.h"


using namespace SPH;
using namespace std;
 
DemoBase::DemoBase()
{
	m_numberOfStepsPerRenderUpdate = 4;
	m_sceneFile = "";
	m_renderWalls = 4;
	m_doPause = true;
	m_pauseAt = -1.0;
	m_useParticleCaching = true;
	m_enablePartioExport = false;
	m_framesPerSecond = 25;
	m_simulationMethodChangedFct = NULL;
}

DemoBase::~DemoBase()
{

}

void DemoBase::init(int argc, char **argv, const char *demoName)
{
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

	if (m_sceneFile != "")
		SceneLoader::readScene(m_sceneFile.c_str(), m_scene);
	else
		return;

	// OpenGL
	MiniGL::init(argc, argv, 1024, 768, 0, 0, demoName);
	MiniGL::initLights();
	MiniGL::getOpenGLVersion(m_context_major_version, m_context_minor_version);
	MiniGL::setViewport(40.0, 0.1f, 500.0, Vector3r(0.0, 3.0, 8.0), Vector3r(0.0, 0.0, 0.0));
	MiniGL::setSelectionFunc(selection, this);

	if (MiniGL::checkOpenGLVersion(3, 3))
		initShaders();
}

void DemoBase::cleanup()
{
	delete m_simulationMethod.simulation;
	delete TimeManager::getCurrent();

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
	string vertFile = getDataPath() + "/shaders/vs_points.glsl";
	string fragFile = getDataPath() + "/shaders/fs_points.glsl";
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

	vertFile = getDataPath() + "/shaders/vs_smooth.glsl";
	fragFile = getDataPath() + "/shaders/fs_smooth.glsl";
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
	glUniform1f(m_shader.getUniform("viewport_width"), (float)viewport[2]);
	glUniform1f(m_shader.getUniform("radius"), (float)m_scene.particleRadius);
	glUniform1f(m_shader.getUniform("max_velocity"), 25.0);
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

void DemoBase::initParameters()
{
	TwRemoveAllVars(MiniGL::getTweakBar());
	m_parameters.clear();

	MiniGL::initTweakBarParameters();
	TwAddVarRW(MiniGL::getTweakBar(), "Pause", TW_TYPE_BOOLCPP, &m_doPause, " label='Pause' group=Simulation key=SPACE ");
	TwAddVarRW(MiniGL::getTweakBar(), "PauseAt", TW_TYPE_REAL, &m_pauseAt, " label='Pause simulation at' step=0.001 precision=4 group=Simulation ");
	TwAddVarRW(MiniGL::getTweakBar(), "numberOfStepsPerRenderUpdate", TW_TYPE_UINT32, &m_numberOfStepsPerRenderUpdate, " label='# time steps / update' min=1 group=Simulation ");
	TwAddVarRW(MiniGL::getTweakBar(), "EnablePartioExport", TW_TYPE_BOOLCPP, &m_enablePartioExport, " label='Partio export' group=Simulation ");
	TwAddVarRW(MiniGL::getTweakBar(), "FramesPerSecond", TW_TYPE_UINT32, &m_framesPerSecond, " label='Export FPS' group=Simulation ");
	TwAddSeparator(MiniGL::getTweakBar(), NULL, " group=Simulation");

	m_parameters.push_back(Parameter(ParameterIDs::TimeStepSize, "TimeStepSize", TW_TYPE_REAL, " label='Time step size'  min=0.0 max = 0.1 step=0.001 precision=4 group=Simulation ", this));

	m_parameters.push_back(Parameter(ParameterIDs::IterationCount, "IterationCount", TW_TYPE_UINT32, " label='Iterations' readonly=true group=Simulation ", this));

	m_parameters.push_back(Parameter(ParameterIDs::Gravitation, "Gravitation", TW_TYPE_DIR3R, " label='Gravitation' group=Simulation", this));

	TwType enumType = TwDefineEnum("SimulationMethodType", NULL, 0);
	m_parameters.push_back(Parameter(ParameterIDs::SimMethod, "SimulationMethod", enumType, " label='Simulation method' enum='0 {WCSPH}, 1 {PCISPH}, 2 {PBF}, 3 {IISPH}, 4 {DFSPH}, 5 {PF}' group=Simulation", this));

	m_parameters.push_back(Parameter(ParameterIDs::MaxIterations, "MaxIterations", TW_TYPE_UINT32, " label='Max. iterations' group=Simulation ", this));
	m_parameters.push_back(Parameter(ParameterIDs::MaxError, "MaxError", TW_TYPE_REAL, " label='Max.density error(%)'  min=0.00001 precision=4 group=Simulation ", this));

	TwType enumTypeVisco = TwDefineEnum("ViscosityMethod", NULL, 0);
	m_parameters.push_back(Parameter(ParameterIDs::ViscosityMethod, "ViscosityMethod", enumTypeVisco, " label='Viscosity' enum='0 {None}, 1 {Standard}, 2 {XSPH}, 3 {Bender2017}' group=Simulation", this));
	m_parameters.push_back(Parameter(ParameterIDs::Viscosity, "Viscosity", TW_TYPE_REAL, " label='Viscosity coefficient'  min=0.0 step=0.001 precision=4 group=Simulation ", this));

	if (m_simulationMethod.simulation->getViscosityMethod() == ViscosityMethods::Bender2017)
	{
		m_parameters.push_back(Parameter(ParameterIDs::ViscoMaxIter, "ViscoMaxIterations", TW_TYPE_UINT32, " label='Max. iterations (visco)' group=Simulation ", this));
		m_parameters.push_back(Parameter(ParameterIDs::ViscoMaxError, "ViscoMaxError", TW_TYPE_REAL, " label='Max. visco error'  min=0.001 precision=3 group=Simulation ", this));
	}
	
	TwType enumTypeST = TwDefineEnum("SurfaceTensionMethod", NULL, 0);
	m_parameters.push_back(Parameter(ParameterIDs::SurfaceTensionMethod, "SurfaceTensionMethod", enumTypeST, " label='Surface tension' enum='0 {None}, 1 {Becker & Teschner 2007}, 2 {Akinci et al. 2013}, 3 {He et al. 2014}' group=Simulation", this));
	m_parameters.push_back(Parameter(ParameterIDs::SurfaceTension, "SurfaceTension", TW_TYPE_REAL, " label='Surface tension coefficient'  min=0.0 step=0.001 precision=4 group=Simulation ", this));

	TwType enumType3 = TwDefineEnum("CFL_Method", NULL, 0);
	m_parameters.push_back(Parameter(ParameterIDs::CFL_Method, "CFL_Method", enumType3, " label='CFL - method' enum='0 {None}, 1 {CFL}, 2 {CFL - iterations}' group=CFL ", this));

	m_parameters.push_back(Parameter(ParameterIDs::CFL_Factor, "CFL_Factor", TW_TYPE_REAL, " label='CFL - factor'  min=0.001 precision=3 group=CFL ", this));
	m_parameters.push_back(Parameter(ParameterIDs::CFL_MaxTimeStepSize, "CFL_MaxTimeStepSize", TW_TYPE_REAL, " label='CFL - max. time step size' min=0.0 precision=4 group=CFL ", this));

	TwType enumType4 = TwDefineEnum("Kernel_Method", NULL, 0);
	m_parameters.push_back(Parameter(ParameterIDs::Kernel_Method, "Kernel_Method", enumType4, " label='Kernel' enum='0 {Cubic spline}, 1 {Poly6}, 2 {Spiky}, 3 {Precomputed cubic spline}' group=Kernel ", this));

	TwType enumType5 = TwDefineEnum("GradKernel_Method", NULL, 0);
	m_parameters.push_back(Parameter(ParameterIDs::GradKernel_Method, "GradKernel_Method", enumType5, " label='Gradient of kernel' enum='0 {Cubic spline}, 1 {Poly6}, 2 {Spiky}, 3 {Precomputed cubic spline}' group=Kernel ", this));

	if (m_simulationMethod.simulationMethod == SimulationMethods::PBF)
	{
		TwType enumType2 = TwDefineEnum("VelocityUpdateMethodType", NULL, 0);
		m_parameters.push_back(Parameter(ParameterIDs::VelocityUpdateMethod, "VelocityUpdateMethod", enumType2, " label='Velocity update method' enum='0 {First Order Update}, 1 {Second Order Update}' group=PBF", this));
	}

	if (m_simulationMethod.simulationMethod == SimulationMethods::WCSPH)
	{
		m_parameters.push_back(Parameter(ParameterIDs::WCSPH_Stiffness, "WCSPH_Stiffness", TW_TYPE_REAL, " label='Stiffness (B)' min=0.0 group=WCSPH", this));
		m_parameters.push_back(Parameter(ParameterIDs::WCSPH_Exponent, "WCSPH_Exponent", TW_TYPE_REAL, " label='Exponent (gamma)' min=0.0 group=WCSPH", this));
	}

	if (m_simulationMethod.simulationMethod == SimulationMethods::DFSPH)
	{
		m_parameters.push_back(Parameter(ParameterIDs::IterationCountV, "IterationCountV", TW_TYPE_UINT32, " label='Iterations (divergence)' readonly=true group=DFSPH ", this));
		m_parameters.push_back(Parameter(ParameterIDs::DFSPH_EnableDivergenceSolver, "DFSPH_EnableDivergenceSolver", TW_TYPE_BOOL32, " label='Enable divergence solver' group=DFSPH", this));
		m_parameters.push_back(Parameter(ParameterIDs::MaxIterationsV, "MaxIterationsV", TW_TYPE_UINT32, " label='Max. iterations (divergence)' group=DFSPH ", this));
		m_parameters.push_back(Parameter(ParameterIDs::MaxErrorV, "MaxErrorV", TW_TYPE_REAL, " label='Max. divergence error(%)'  min=0.00001 precision=4 group=DFSPH ", this));
	}

	if (m_simulationMethod.simulationMethod == SimulationMethods::PF)
	{
		// no special parameters yet
	}


	for (auto &p : m_parameters)
	{
		TwAddVarCB(MiniGL::getTweakBar(), p.name.c_str(), p.type, setParameter, getParameter, &p, p.tweakBarDefinition.c_str());
	}

	TwType enumTypeWalls = TwDefineEnum("RenderWalls", NULL, 0);
	TwAddVarRW(MiniGL::getTweakBar(), "RenderWalls", enumTypeWalls, &m_renderWalls, " label='Render walls' enum='0 {None}, 1 {Particles (all)}, 2 {Particles (no walls)}, 3 {Geometry (all)}, 4 {Geometry (no walls)}' group=Visualization ");
}


void DemoBase::buildModel()
{
	TimeManager::getCurrent()->setTimeStepSize(m_scene.timeStepSize);

	std::vector<Vector3r> fluidParticles;
	std::vector<Vector3r> fluidVelocities;
	initFluidData(fluidParticles, fluidVelocities);

	m_simulationMethod.model.setParticleRadius(m_scene.particleRadius);
	setPauseAt(m_scene.pauseAt);
	setNumberOfStepsPerRenderUpdate(m_scene.numberOfStepsPerRenderUpdate);

	m_simulationMethod.model.initModel((unsigned int)fluidParticles.size(), fluidParticles.data());

	std::cout << "Number of fluid particles: " << fluidParticles.size() << "\n";

	unsigned int nBoundaryParticles = 0;
	for (unsigned int i = 0; i < m_simulationMethod.model.numberOfRigidBodyParticleObjects(); i++)
		nBoundaryParticles += m_simulationMethod.model.getRigidBodyParticleObject(i)->numberOfParticles();

	std::cout << "Number of boundary particles: " << nBoundaryParticles << "\n";

	m_simulationMethod.simulation = new TimeStepDFSPH(&m_simulationMethod.model);
	m_simulationMethod.simulationMethod = SimulationMethods::DFSPH;
	m_simulationMethod.model.setKernel(3);
	m_simulationMethod.model.setGradKernel(3);

	setSimulationMethod((SimulationMethods) m_scene.simulationMethod);

	m_simulationMethod.simulation->setCflMethod(m_scene.cflMethod);
	m_simulationMethod.simulation->setCflFactor(m_scene.cflFactor);
	m_simulationMethod.simulation->setCflMaxTimeStepSize(m_scene.cflMaxTimeStepSize);
	m_simulationMethod.simulation->setMaxIterations(m_scene.maxIterations);
	m_simulationMethod.simulation->setMaxError(m_scene.maxError);
	m_simulationMethod.simulation->setMaxIterationsV(m_scene.maxIterationsV);
	m_simulationMethod.simulation->setMaxErrorV(m_scene.maxErrorV);
	m_simulationMethod.simulation->setViscosityMethod((ViscosityMethods) m_scene.viscosityMethod);
	m_simulationMethod.simulation->setSurfaceTensionMethod((SurfaceTensionMethods)m_scene.surfaceTensionMethod);
	
	m_simulationMethod.model.setEnableDivergenceSolver(m_scene.enableDivergenceSolver);
	m_simulationMethod.model.setViscosity(m_scene.viscosity);
	m_simulationMethod.model.setSurfaceTension(m_scene.surfaceTension);
	m_simulationMethod.model.setDensity0(m_scene.density0);
	m_simulationMethod.model.setGravitation(m_scene.gravitation);
	m_simulationMethod.model.setVelocityUpdateMethod(m_scene.velocityUpdateMethod);
	m_simulationMethod.model.setStiffness(m_scene.stiffness);
	m_simulationMethod.model.setExponent(m_scene.exponent);

	setEnablePartioExport(m_scene.enablePartioExport);
	setFramesPerSecond(m_scene.partioFPS);

	if (m_simulationMethod.simulation->getViscosityMethod() == ViscosityMethods::Bender2017)
	{
		((Viscosity_Bender2017*)m_simulationMethod.simulation->getViscosityBase())->setMaxError(m_scene.viscoMaxError);
		((Viscosity_Bender2017*)m_simulationMethod.simulation->getViscosityBase())->setMaxIter(m_scene.viscoMaxIter);
	}

	initParameters();
}


void DemoBase::initFluidData(std::vector<Vector3r> &fluidParticles, std::vector<Vector3r> &fluidVelocities)
{
	std::cout << "Initialize fluid particles\n";
	createFluidBlocks(fluidParticles);

	std::string base_path = FileSystem::getFilePath(m_sceneFile);

	unsigned int startIndex = 0;
	unsigned int endIndex = 0;
	for (unsigned int i = 0; i < m_scene.fluidModels.size(); i++)
	{
		string fileName = base_path + "/" + m_scene.fluidModels[i]->samplesFile;
		PartioReaderWriter::readParticles(fileName, m_scene.fluidModels[i]->translation, m_scene.fluidModels[i]->rotation, m_scene.fluidModels[i]->scale, fluidParticles, fluidVelocities);
		m_simulationMethod.model.setParticleRadius(m_scene.particleRadius);
	}
}


void DemoBase::createFluidBlocks(std::vector<Vector3r> &fluidParticles)
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

void TW_CALL DemoBase::setParameter(const void *value, void *clientData)
{
	Parameter *p = ((Parameter*)clientData);
	DemoBase *base = p->base;
	SimulationMethod &sm = base->getSimulationMethod();

	if (p->id == ParameterIDs::TimeStepSize)
	{
		const Real val = *(const Real *)(value);
		TimeManager::getCurrent()->setTimeStepSize(val);
	}
	else if (p->id == ParameterIDs::Gravitation)
	{
		const Real *val = (const Real*)(value);
		const Vector3r v(val[0], val[1], val[2]);
		sm.model.setGravitation(v);
	}
	else if (p->id == ParameterIDs::SimMethod)
	{
		const short val = *(const short *)(value);
		base->setSimulationMethod((SimulationMethods)val);
	}
	else if (p->id == ParameterIDs::VelocityUpdateMethod)
	{
		const short val = *(const short *)(value);
		sm.model.setVelocityUpdateMethod((unsigned int)val);
	}
	else if (p->id == ParameterIDs::Viscosity)
	{
		const Real val = *(const Real *)(value);
		sm.model.setViscosity(val);
	}
	else if (p->id == ParameterIDs::ViscoMaxIter)
	{
		const unsigned int val = *(const unsigned int *)(value);
		if (sm.simulation->getViscosityMethod() == ViscosityMethods::Bender2017)
			((Viscosity_Bender2017*)sm.simulation->getViscosityBase())->setMaxIter(val);
	}
	else if (p->id == ParameterIDs::ViscoMaxError)
	{
		const Real val = *(const Real *)(value);
		if (sm.simulation->getViscosityMethod() == ViscosityMethods::Bender2017)
			((Viscosity_Bender2017*)sm.simulation->getViscosityBase())->setMaxError(val);
	}
	else if (p->id == ParameterIDs::SurfaceTension)
	{
		const Real val = *(const Real *)(value);
		sm.model.setSurfaceTension(val);
	}
	else if (p->id == ParameterIDs::SurfaceTensionMethod)
	{
		const short val = *(const short *)(value);
		sm.simulation->setSurfaceTensionMethod((SurfaceTensionMethods)val);
	}
	else if (p->id == ParameterIDs::ViscosityMethod)
	{
		const short val = *(const short *)(value);
		sm.simulation->setViscosityMethod((ViscosityMethods)val);
		if (sm.simulation->getViscosityMethod() == ViscosityMethods::Bender2017)
		{
			((Viscosity_Bender2017*)sm.simulation->getViscosityBase())->setMaxIter(base->m_scene.viscoMaxIter);
			((Viscosity_Bender2017*)sm.simulation->getViscosityBase())->setMaxError(base->m_scene.viscoMaxError);
			base->initParameters();
		}		
	}
	else if (p->id == ParameterIDs::WCSPH_Stiffness)
	{
		const Real val = *(const Real *)(value);
		sm.model.setStiffness(val);
	}
	else if (p->id == ParameterIDs::WCSPH_Exponent)
	{
		const Real val = *(const Real *)(value);
		sm.model.setExponent(val);
	}
	else if (p->id == ParameterIDs::DFSPH_EnableDivergenceSolver)
	{
		const bool val = *(const bool *)(value);
		sm.model.setEnableDivergenceSolver(val);
	}
	else if (p->id == ParameterIDs::CFL_Method)
	{
		const short val = *(const short *)(value);
		sm.simulation->setCflMethod((unsigned int)val);
	}
	else if (p->id == ParameterIDs::CFL_Factor)
	{
		const Real val = *(const Real *)(value);
		sm.simulation->setCflFactor(val);
	}
	else if (p->id == ParameterIDs::CFL_MaxTimeStepSize)
	{
		const Real val = *(const Real *)(value);
		sm.simulation->setCflMaxTimeStepSize(val);
	}
	else if (p->id == ParameterIDs::Kernel_Method)
	{
		const short val = *(const short *)(value);
		if (val != sm.model.getKernel())
		{
			sm.model.setKernel((unsigned int)val);
			sm.model.updateBoundaryPsi();
		}
	}
	else if (p->id == ParameterIDs::GradKernel_Method)
	{
		const short val = *(const short *)(value);
		sm.model.setGradKernel((unsigned int)val);
	}
	else if (p->id == ParameterIDs::MaxIterations)
	{
		sm.simulation->setMaxIterations(*(unsigned int *)(value));
	}
	else if (p->id == ParameterIDs::MaxError)
	{
		sm.simulation->setMaxError(*(Real *)(value));
	}
	else if (p->id == ParameterIDs::MaxIterationsV)
	{
		sm.simulation->setMaxIterationsV(*(unsigned int *)(value));
	}
	else if (p->id == ParameterIDs::MaxErrorV)
	{
		sm.simulation->setMaxErrorV(*(Real *)(value));
	}
}

void TW_CALL DemoBase::getParameter(void *value, void *clientData)
{
	Parameter *p = ((Parameter*)clientData);
	DemoBase *base = p->base;
	SimulationMethod &sm = base->getSimulationMethod();

	if (p->id == ParameterIDs::TimeStepSize)
	{
		*(Real *)(value) = TimeManager::getCurrent()->getTimeStepSize();
	}
	else if (p->id == ParameterIDs::IterationCount)
	{
		if (sm.simulation != NULL)
			*((unsigned int*)value) = sm.simulation->getIterationCount();
		else
			*((unsigned int*)value) = 0;
	}
	else if (p->id == ParameterIDs::IterationCountV)
	{
		if (sm.simulation != NULL)
			*((unsigned int*)value) = sm.simulation->getIterationCountV();
		else
			*((unsigned int*)value) = 0;
	}
	else if (p->id == ParameterIDs::Gravitation)
	{
		const Vector3r &val = sm.model.getGravitation();
		((Real*)value)[0] = val[0];
		((Real*)value)[1] = val[1];
		((Real*)value)[2] = val[2];
	}
	else if (p->id == ParameterIDs::SimMethod)
	{
		*(short *)(value) = sm.simulationMethod;
	}
	else if (p->id == ParameterIDs::VelocityUpdateMethod)
	{
		*(short *)(value) = (short)sm.model.getVelocityUpdateMethod();
	}
	else if (p->id == ParameterIDs::Viscosity)
	{
		*(Real *)(value) = sm.model.getViscosity();
	}
	else if (p->id == ParameterIDs::ViscoMaxIter)
	{
		if (sm.simulation->getViscosityMethod() == ViscosityMethods::Bender2017)
			*(unsigned int *)(value) = ((Viscosity_Bender2017*)sm.simulation->getViscosityBase())->getMaxIter();
	}
	else if (p->id == ParameterIDs::ViscoMaxError)
	{
		if (sm.simulation->getViscosityMethod() == ViscosityMethods::Bender2017)
			*(Real *)(value) = ((Viscosity_Bender2017*)sm.simulation->getViscosityBase())->getMaxError();
	}
	else if (p->id == ParameterIDs::SurfaceTension)
	{
		*(Real *)(value) = sm.model.getSurfaceTension();
	}
	else if (p->id == ParameterIDs::SurfaceTensionMethod)
	{
		*(short *)(value) = (short)sm.simulation->getSurfaceTensionMethod();
	}
	else if (p->id == ParameterIDs::ViscosityMethod)
	{
		*(short *)(value) = (short)sm.simulation->getViscosityMethod();
	}
	else if (p->id == ParameterIDs::WCSPH_Stiffness)
	{
		*(Real *)(value) = sm.model.getStiffness();
	}
	else if (p->id == ParameterIDs::WCSPH_Exponent)
	{
		*(Real *)(value) = sm.model.getExponent();
	}
	else if (p->id == ParameterIDs::DFSPH_EnableDivergenceSolver)
	{
		*(bool *)(value) = sm.model.getEnableDivergenceSolver();
	}
	else if (p->id == ParameterIDs::CFL_Method)
	{
		*(short *)(value) = (short)sm.simulation->getCflMethod();
	}
	else if (p->id == ParameterIDs::CFL_Factor)
	{
		*(Real *)(value) = sm.simulation->getCflFactor();
	}
	else if (p->id == ParameterIDs::CFL_MaxTimeStepSize)
	{
		*(Real *)(value) = sm.simulation->getCflMaxTimeStepSize();
	}
	else if (p->id == ParameterIDs::Kernel_Method)
	{
		*(short *)(value) = (short)sm.model.getKernel();
	}
	else if (p->id == ParameterIDs::GradKernel_Method)
	{
		*(short *)(value) = (short)sm.model.getGradKernel();
	}
	else if (p->id == ParameterIDs::MaxIterations)
	{
		*(unsigned int *)(value) = sm.simulation->getMaxIterations();
	}
	else if (p->id == ParameterIDs::MaxError)
	{
		*(Real *)(value) = sm.simulation->getMaxError();
	}
	else if (p->id == ParameterIDs::MaxIterationsV)
	{
		*(unsigned int *)(value) = sm.simulation->getMaxIterationsV();
	}
	else if (p->id == ParameterIDs::MaxErrorV)
	{
		*(Real *)(value) = sm.simulation->getMaxErrorV();
	}
}

void DemoBase::renderFluid()
{
	// Draw simulation model

	const unsigned int nParticles = m_simulationMethod.model.numParticles();

	float surfaceColor[4] = { 0.2f, 0.6f, 0.8f, 1 };
	float speccolor[4] = { 1.0, 1.0, 1.0, 1.0 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, speccolor);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);
	glColor3fv(surfaceColor);


	const Real supportRadius = m_simulationMethod.model.getSupportRadius();
	Real vmax = 0.4*2.0*supportRadius / TimeManager::getCurrent()->getTimeStepSize();
	Real vmin = 0.0;

	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		float fluidColor[4] = { 0.3f, 0.5f, 0.9f, 1.0f };
		pointShaderBegin(&fluidColor[0]);

		if (m_simulationMethod.model.numParticles() > 0)
		{
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, &m_simulationMethod.model.getPosition(0, 0));
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 0, &m_simulationMethod.model.getVelocity(0, 0));
			glDrawArrays(GL_POINTS, 0, m_simulationMethod.model.numParticles());
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
			Real v = m_simulationMethod.model.getVelocity(0, i).norm();
			v = 0.5*((v - vmin) / (vmax - vmin));
			v = min(128.0*v*v, 0.5);
			float fluidColor[4] = { 0.2f, 0.2f, 0.2f, 1.0 };
			MiniGL::hsvToRgb(0.55f, 1.0f, 0.5f + (float)v, fluidColor);

			glColor3fv(fluidColor);
			glVertex3v(&m_simulationMethod.model.getPosition(0, i)[0]);
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
			glUniform1f(m_shader.getUniform("radius"), (float)m_scene.particleRadius*1.05f);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, &m_simulationMethod.model.getPosition(0, 0));
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 0, &m_simulationMethod.model.getVelocity(0, 0));
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
			glVertex3v(&getSimulationMethod().model.getPosition(0, getSelectedParticles()[i])[0]);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}

	MiniGL::drawTime(TimeManager::getCurrent()->getTime());
}


void DemoBase::mouseMove(int x, int y, void *clientData)
{
	DemoBase *base = (DemoBase*)clientData;

	Vector3r mousePos;
	MiniGL::unproject(x, y, mousePos);
	const Vector3r diff = mousePos - base->m_oldMousePos;

	TimeManager *tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	for (unsigned int j = 0; j < base->m_selectedParticles.size(); j++)
	{
		base->m_simulationMethod.model.getVelocity(0, base->m_selectedParticles[j]) += 5.0*diff / h;
	}
	base->m_oldMousePos = mousePos;
}

void DemoBase::selection(const Eigen::Vector2i &start, const Eigen::Vector2i &end, void *clientData)
{
	DemoBase *base = (DemoBase*)clientData;

	std::vector<unsigned int> hits;
	base->m_selectedParticles.clear();
	Selection::selectRect(start, end, &base->m_simulationMethod.model.getPosition(0, 0), 
		&base->m_simulationMethod.model.getPosition(0, base->m_simulationMethod.model.numParticles() - 1), 
		base->m_selectedParticles);
	if (base->m_selectedParticles.size() > 0)
		MiniGL::setMouseMoveFunc(GLUT_MIDDLE_BUTTON, mouseMove);
	else
		MiniGL::setMouseMoveFunc(-1, NULL);

	MiniGL::unproject(end[0], end[1], base->m_oldMousePos);
}

void DemoBase::setSimulationMethod(SimulationMethods method)
{
	if ((method < 0) || (method >= SimulationMethods::NUM_METHODS))
		method = SimulationMethods::DFSPH;

	if (method != m_simulationMethod.simulationMethod)
	{
		delete m_simulationMethod.simulation;

		if (method == SimulationMethods::WCSPH)
		{
			m_simulationMethod.simulation = new TimeStepWCSPH(&m_simulationMethod.model);
			m_simulationMethod.simulation->setCflMethod(0);
			m_simulationMethod.model.setGradKernel(0);
			if (m_simulationMethod.model.getKernel() != 0)
			{
				m_simulationMethod.model.setKernel(0);
				m_simulationMethod.model.updateBoundaryPsi();
			}
			TimeManager::getCurrent()->setTimeStepSize(0.001);
		}
		else if (method == SimulationMethods::PCISPH)
		{
			m_simulationMethod.simulation = new TimeStepPCISPH(&m_simulationMethod.model);
			m_simulationMethod.model.setGradKernel(0);
			if (m_simulationMethod.model.getKernel() != 0)
			{
				m_simulationMethod.model.setKernel(0);
				m_simulationMethod.model.updateBoundaryPsi();
			}
		}
		else if (method == SimulationMethods::PBF)
		{
			m_simulationMethod.simulation = new TimeStepPBF(&m_simulationMethod.model);
			m_simulationMethod.model.setGradKernel(2);
			if (m_simulationMethod.model.getKernel() != 1)
			{
				m_simulationMethod.model.setKernel(1);
				m_simulationMethod.model.updateBoundaryPsi();
			}
		}
		else if (method == SimulationMethods::IISPH)
		{
			m_simulationMethod.simulation = new TimeStepIISPH(&m_simulationMethod.model);
			m_simulationMethod.model.setGradKernel(0);
			if (m_simulationMethod.model.getKernel() != 0)
			{
				m_simulationMethod.model.setKernel(0);
				m_simulationMethod.model.updateBoundaryPsi();
			}
		}
		else if (method == SimulationMethods::DFSPH)
		{
			m_simulationMethod.simulation = new TimeStepDFSPH(&m_simulationMethod.model);
			m_simulationMethod.model.setGradKernel(3);
			if (m_simulationMethod.model.getKernel() != 3)
			{
				m_simulationMethod.model.setKernel(3);
				m_simulationMethod.model.updateBoundaryPsi();
			}
		}
		else if (method == SimulationMethods::PF)
		{
			m_simulationMethod.simulation = new TimeStepPF(&m_simulationMethod.model);
			m_simulationMethod.model.setGradKernel(3);
			if (m_simulationMethod.model.getKernel() != 3)
			{
				m_simulationMethod.model.setKernel(3);
				m_simulationMethod.model.updateBoundaryPsi();
			}
		}
		m_simulationMethod.simulationMethod = method;

		initParameters();
		if (m_simulationMethodChangedFct)
			m_simulationMethodChangedFct();	
	}
}
