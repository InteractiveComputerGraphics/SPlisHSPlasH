#include "Simulator_GUI_TweakBar.h"
#include "GUI/OpenGL/MiniGL.h"
#include "GUI/TweakBar/TweakBarParameters.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "../OpenGL/Simulator_OpenGL.h"
#include "SPlisHSPlasH/Utilities/SceneLoader.h"
#include "GUI/OpenGL/Selection.h"

using namespace SPH;
using namespace Utilities;

Simulator_GUI_TweakBar::Simulator_GUI_TweakBar(SimulatorBase *base) :
	Simulator_GUI_Base(base)
{	
	m_tweakBar = nullptr;
	m_currentFluidModel = 0;
}

Simulator_GUI_TweakBar::~Simulator_GUI_TweakBar(void)
{	
}

void Simulator_GUI_TweakBar::init(int argc, char **argv, const char *name)
{
	MiniGL::init(argc, argv, 1280, 960, name);
	MiniGL::initLights();

	auto scene = m_simulatorBase->getScene();
	const bool sim2D = scene.sim2D;
	if (sim2D)
		MiniGL::setViewport(40.0, 0.1f, 500.0, scene.camPosition, scene.camLookat);
	else
		MiniGL::setViewport(40.0, 0.1f, 500.0, scene.camPosition, scene.camLookat);
	MiniGL::setSelectionFunc(selection, this);
	MiniGL::addKeyFunc('i', std::bind(&Simulator_GUI_TweakBar::particleInfo, this));
	MiniGL::addKeyFunc('s', std::bind(&SimulatorBase::saveState, m_simulatorBase));
#ifdef WIN32
	MiniGL::addKeyFunc('l', std::bind(&SimulatorBase::loadStateDialog, m_simulatorBase));
#endif 

	if (MiniGL::checkOpenGLVersion(3, 3))
		Simulator_OpenGL::initShaders(m_simulatorBase->getExePath() + "/resources/shaders");

	const int width = MiniGL::getWidth();
	const int height = MiniGL::getHeight();

	// Initialize AntTweakBar
	if (!TwInit(TW_OPENGL, NULL))
	{
		// A fatal error occured    
		fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
		exit(1);
	}
	TwWindowSize(width, height);
	initTweakBar();

	MiniGL::addReshapeFunc([](int width, int height) { TwWindowSize(width, height); });
	MiniGL::addKeyboardFunc([](int key, int scancode, int action, int mods) -> bool { return TwEventKeyGLFW(key, action); });
	MiniGL::addCharFunc([](int key, int action) -> bool { return TwEventCharGLFW(key, action); });
	MiniGL::addMousePressFunc([](int button, int action, int mods) -> bool { return TwEventMouseButtonGLFW(button, action); });
	MiniGL::addMouseMoveFunc([](int x, int y) -> bool { return TwEventMousePosGLFW(x, y); });
	MiniGL::addMouseWheelFunc([](int pos, double xoffset, double yoffset) -> bool { return TwEventMouseWheelGLFW(pos); });

	MiniGL::setClientIdleFunc(std::bind(&SimulatorBase::timeStep, m_simulatorBase));
	MiniGL::addKeyFunc('r', std::bind(&SimulatorBase::reset, m_simulatorBase));
	MiniGL::setClientSceneFunc(std::bind(&Simulator_GUI_TweakBar::render, this));
}


void Simulator_GUI_TweakBar::initTweakBar()
{
	// Create a tweak bar
	m_tweakBar = TwNewBar("TweakBar");
	TwDefine(" GLOBAL help='MiniGL TweakBar.' "); // Message added to the help bar.
	//TwDefine(" TweakBar size='300 900' valueswidth=120 position='5 5' color='96 200 224' text=dark "); // change default tweak bar size and color
	TwDefine(" TweakBar size='300 900' valueswidth=120 position='5 5'"); // change default tweak bar size and color

	initTweakBarParameters();
}


void Simulator_GUI_TweakBar::initTweakBarParameters()
{
	TwAddVarCB(getTweakBar(), "Time", TW_TYPE_REAL, setTimeCB, getTimeCB, nullptr, " label='Time' precision=5 group='General'");

	// Add callback to toggle auto-rotate mode (callback functions are defined above).
	TwAddVarCB(getTweakBar(), "Wireframe", TW_TYPE_BOOL32, setWireframeCB, getWireframeCB, nullptr,
		" label='Wireframe' key=w help='Toggle wireframe mode.' group='Visualization' ");

	TwAddVarCB(getTweakBar(), "Rotation", TW_TYPE_QUAT4F, setRotationCB, getRotationCB, nullptr,
		" label='Rotation' open help='Change the rotation.' group='Visualization' ");
}

void Simulator_GUI_TweakBar::initSimulationParameterGUI()
{
	TwRemoveAllVars(getTweakBar());
	TweakBarParameters::cleanup();

	initTweakBarParameters();

	Simulation *sim = Simulation::getCurrent();
	TweakBarParameters::createParameterGUI(getTweakBar());
	if (m_simulatorBase)
		TweakBarParameters::createParameterObjectGUI(getTweakBar(), m_simulatorBase);
	TweakBarParameters::createParameterObjectGUI(getTweakBar(), sim);
	TweakBarParameters::createParameterObjectGUI(getTweakBar(), (GenParam::ParameterObject*) sim->getTimeStep());
#ifdef USE_DEBUG_TOOLS
	TweakBarParameters::createParameterObjectGUI(getTweakBar(), (GenParam::ParameterObject*) sim->getDebugTools());
#endif


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
		TwAddVarCB(getTweakBar(), "CurrentFluidModel", enumType, setCurrentFluidModel, getCurrentFluidModel, this, enumStr.c_str());
	}

	// show GUI only for currently selected fluid model
	unsigned int i = m_currentFluidModel;
	FluidModel *model = sim->getFluidModel(m_currentFluidModel);

	m_colorFieldNames.clear();
	m_colorFieldNames.resize(model->numberOfFields());
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
			m_colorFieldNames[idx] = field.name;
			idx++;
		}
	}
	std::string enumStr = " label='Color field' enum='" + oss.str() + "' group='" + model->getId() + "'";
	enumStr = enumStr + " help='Choose vector or scalar field for particle coloring.'";
	TwAddVarCB(getTweakBar(), "ColorField", enumType, setColorField, getColorField, this, enumStr.c_str());

	TwType enumType2 = TwDefineEnum("ColorMapTypeEnum", NULL, 0);
	std::string str = " label='Color map' enum='0 {None}, 1 {Jet}, 2 {Plasma}' group='" + model->getId() + "' help='Choose a color map.'";
	TwAddVarCB(getTweakBar(), "ColorMapType", enumType2, setColorMapType, getColorMapType, this, str.c_str());
	str = " label='Min. value (shader)' step=0.001 precision=3 group='" + model->getId() + "' help='Minimal value used for color-coding the color field in the rendering process.'";
	TwAddVarCB(getTweakBar(), "RenderMinValue", TW_TYPE_REAL, setRenderMinValue, getRenderMinValue, this, str.c_str());
	str = " label='Max. value (shader)' step=0.001 precision=3 group='" + model->getId() + "' help='Maximal value used for color-coding the color field in the rendering process.'";
	TwAddVarCB(getTweakBar(), "RenderMaxValue", TW_TYPE_REAL, setRenderMaxValue, getRenderMaxValue, this, str.c_str());

	TweakBarParameters::createParameterObjectGUI(getTweakBar(), model);
	TweakBarParameters::createParameterObjectGUI(getTweakBar(), (GenParam::ParameterObject*) model->getDragBase());
	TweakBarParameters::createParameterObjectGUI(getTweakBar(), (GenParam::ParameterObject*) model->getSurfaceTensionBase());
	TweakBarParameters::createParameterObjectGUI(getTweakBar(), (GenParam::ParameterObject*) model->getViscosityBase());
	TweakBarParameters::createParameterObjectGUI(getTweakBar(), (GenParam::ParameterObject*) model->getVorticityBase());
	TweakBarParameters::createParameterObjectGUI(getTweakBar(), (GenParam::ParameterObject*) model->getElasticityBase());
	TwDefine((std::string("TweakBar/FluidModel group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/'Drag force' group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/'Surface tension' group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/Viscosity group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/Vorticity group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/'Elasticity' group='") + model->getId() + "'").c_str());
}

void Simulator_GUI_TweakBar::initParameterGUI()
{
	TwRemoveAllVars(getTweakBar());
}

void Simulator_GUI_TweakBar::update()
{
	TwRefreshBar(getTweakBar());
	TwDraw();
}

void Simulator_GUI_TweakBar::cleanup()
{
	TwDeleteBar(getTweakBar());
	TwTerminate();
	MiniGL::getKeyFunc().clear();
}

void Simulator_GUI_TweakBar::render()
{
	float gridColor[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
	const bool sim2D = Simulation::getCurrent()->is2DSimulation();
	if (sim2D)
		MiniGL::drawGrid_xy(gridColor);
	else
		MiniGL::drawGrid_xz(gridColor);

	MiniGL::coordinateSystem();

	Simulation *sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		float fluidColor[4] = { 0.3f, 0.5f, 0.9f, 1.0f };
		MiniGL::hsvToRgb(0.61f - 0.1f*i, 0.66f, 0.9f, fluidColor);
		FluidModel *model = sim->getFluidModel(i);
		SimulatorBase *base = getSimulatorBase();
		Simulator_OpenGL::renderFluid(model, fluidColor, base->getColorMapType(i),
			base->getColorField(i), base->getRenderMinValue(i), base->getRenderMaxValue(i));
		Simulator_OpenGL::renderSelectedParticles(model, getSelectedParticles(), base->getColorMapType(i),
			base->getColorField(i), base->getRenderMinValue(i), base->getRenderMaxValue(i));
	}
	renderBoundary();
	update();
}

void Simulator_GUI_TweakBar::renderBoundary()
{
	Simulation *sim = Simulation::getCurrent();
	SimulatorBase *base = getSimulatorBase();
	SceneLoader::Scene &scene = base->getScene();
	const int renderWalls = base->getValue<int>(SimulatorBase::RENDER_WALLS);

	if (((renderWalls == 1) || (renderWalls == 2)) &&
		(sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012))
	{
		for (int body = sim->numberOfBoundaryModels() - 1; body >= 0; body--)
		{
			if ((renderWalls == 1) || (!scene.boundaryModels[body]->isWall))
			{
				BoundaryModel_Akinci2012 *bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(body));
				Simulator_OpenGL::renderBoundaryParticles(bm, scene.boundaryModels[body]->color.data());
			}
		}
	}
	else if ((renderWalls == 3) || (renderWalls == 4))
	{
		for (int body = sim->numberOfBoundaryModels() - 1; body >= 0; body--)
		{
			if ((renderWalls == 3) || (!scene.boundaryModels[body]->isWall))
			{
				BoundaryModel *bm = sim->getBoundaryModel(body);
				Simulator_OpenGL::renderBoundary(bm, scene.boundaryModels[body]->color.data());
			}
		}
	}
}

void Simulator_GUI_TweakBar::reset()
{
	m_selectedParticles.clear();
}

void Simulator_GUI_TweakBar::selection(const Vector2i &start, const Vector2i &end, void *clientData)
{
	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	Simulation *sim = Simulation::getCurrent();
	std::vector<std::vector<unsigned int>> &selectedParticles = gui->getSelectedParticles();
	selectedParticles.resize(sim->numberOfFluidModels());
	bool selected = false;
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);

		const unsigned int nParticles = model->numActiveParticles();
		if (nParticles != 0)
		{
			std::vector<unsigned int> hits;
			selectedParticles[i].clear();
			Selection::selectRect(start, end, &model->getPosition(0),
				&model->getPosition(model->numActiveParticles() - 1),
				selectedParticles[i]);
			if (selectedParticles[i].size() > 0)
				selected = true;
		}
	}
	if (selected)
		MiniGL::setMouseMoveFunc(2, mouseMove);
	else
		MiniGL::setMouseMoveFunc(-1, NULL);

	MiniGL::unproject(end[0], end[1], gui->m_oldMousePos);
}


void Simulator_GUI_TweakBar::mouseMove(int x, int y, void *clientData)
{
	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	Simulation *sim = Simulation::getCurrent();
	std::vector<std::vector<unsigned int>> &selectedParticles = gui->getSelectedParticles();

	Vector3r mousePos;
	MiniGL::unproject(x, y, mousePos);
	const Vector3r diff = mousePos - gui->m_oldMousePos;

	TimeManager *tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		for (unsigned int j = 0; j < selectedParticles[i].size(); j++)
		{
			model->getVelocity(selectedParticles[i][j]) += 5.0*diff / h;
		}
	}
	gui->m_oldMousePos = mousePos;
}

void Simulator_GUI_TweakBar::particleInfo()
{
	SimulatorBase::particleInfo(m_selectedParticles);
}

void Simulator_GUI_TweakBar::run()
{
	MiniGL::mainLoop();
}

void Simulator_GUI_TweakBar::stop()
{
	MiniGL::leaveMainLoop();
}

void Simulator_GUI_TweakBar::addKeyFunc(char k, std::function<void()> const& func)
{
	MiniGL::addKeyFunc(k, func);
}

void TW_CALL Simulator_GUI_TweakBar::setTimeCB(const void *value, void *clientData)
{
	// read-only: nothing happens here
}

void TW_CALL Simulator_GUI_TweakBar::getTimeCB(void *value, void *clientData)
{
	*(Real *)(value) = TimeManager::getCurrent()->getTime();
}

void TW_CALL Simulator_GUI_TweakBar::setWireframeCB(const void *value, void *clientData)
{
	const int val = *(const int *)(value);
	if (val == 0)
		MiniGL::setDrawMode(GL_FILL);
	else
		MiniGL::setDrawMode(GL_LINE);
}

void TW_CALL Simulator_GUI_TweakBar::getWireframeCB(void *value, void *clientData)
{
	*(int *)(value) = MiniGL::getDrawMode() == GL_LINE;
}

void TW_CALL Simulator_GUI_TweakBar::setRotationCB(const void *value, void *clientData)
{
	const float *val = (const float *)(value);

	Quaternionr q;
	q.x() = (Real)val[0];
	q.y() = (Real)val[1];
	q.z() = (Real)val[2];
	q.w() = -(Real)val[3];
	MiniGL::setRotation(q);
}

void TW_CALL Simulator_GUI_TweakBar::getRotationCB(void *value, void *clientData)
{
	Quaternionr q = MiniGL::getRotation();
	float *val = (float*)(value);
	val[0] = (float)q.x();
	val[1] = (float)q.y();
	val[2] = (float)q.z();
	val[3] = -(float)q.w();
}

void TW_CALL Simulator_GUI_TweakBar::setCurrentFluidModel(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	gui->m_currentFluidModel = val;
	gui->initSimulationParameterGUI();
}

void TW_CALL Simulator_GUI_TweakBar::getCurrentFluidModel(void *value, void *clientData)
{
	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	*(unsigned int *)(value) = gui->m_currentFluidModel;
}

void TW_CALL Simulator_GUI_TweakBar::setColorField(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);

	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	gui->getSimulatorBase()->setColorField(gui->m_currentFluidModel, gui->m_colorFieldNames[val]);
}

void TW_CALL Simulator_GUI_TweakBar::getColorField(void *value, void *clientData)
{
	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	const std::string &fieldName = gui->getSimulatorBase()->getColorField(gui->m_currentFluidModel);
	unsigned int index = 0;
	for (auto i = 0; i < gui->m_colorFieldNames.size(); i++)
	{
		if (gui->m_colorFieldNames[i] == fieldName)
		{
			index = i;
			break;
		}
	}
	*(unsigned int *)(value) = index;
}

void TW_CALL Simulator_GUI_TweakBar::setRenderMaxValue(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);

	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	gui->getSimulatorBase()->setRenderMaxValue(gui->m_currentFluidModel, val);
}

void TW_CALL Simulator_GUI_TweakBar::getRenderMaxValue(void *value, void *clientData)
{
	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	*(Real *)(value) = gui->getSimulatorBase()->getRenderMaxValue(gui->m_currentFluidModel);
}

void TW_CALL Simulator_GUI_TweakBar::setRenderMinValue(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);

	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	gui->getSimulatorBase()->setRenderMinValue(gui->m_currentFluidModel, val);
}

void TW_CALL Simulator_GUI_TweakBar::getRenderMinValue(void *value, void *clientData)
{
	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	*(Real *)(value) = gui->getSimulatorBase()->getRenderMinValue(gui->m_currentFluidModel);
}

void TW_CALL Simulator_GUI_TweakBar::setColorMapType(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);

	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	gui->getSimulatorBase()->setColorMapType(gui->m_currentFluidModel, val);
}

void TW_CALL Simulator_GUI_TweakBar::getColorMapType(void *value, void *clientData)
{
	Simulator_GUI_TweakBar *gui = (Simulator_GUI_TweakBar*)clientData;
	*(unsigned int *)(value) = gui->getSimulatorBase()->getColorMapType(gui->m_currentFluidModel);
}
