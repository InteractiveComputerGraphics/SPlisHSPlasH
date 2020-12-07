#include "Simulator_GUI_imgui.h"
#include "GUI/OpenGL/MiniGL.h"
#include "GUI/imgui/imguiParameters.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "../OpenGL/Simulator_OpenGL.h"
#include "SPlisHSPlasH/Utilities/SceneLoader.h"
#include "GUI/OpenGL/Selection.h"
#include "Utilities/FileSystem.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

using namespace SPH;
using namespace Utilities;

Simulator_GUI_imgui::Simulator_GUI_imgui(SimulatorBase *base) :
	Simulator_GUI_Base(base)
{	
	m_currentFluidModel = 0;
}

Simulator_GUI_imgui::~Simulator_GUI_imgui(void)
{	
	imguiParameters::cleanup();
}

void Simulator_GUI_imgui::init(int argc, char **argv, const char *name)
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
	MiniGL::addKeyFunc('i', std::bind(&Simulator_GUI_imgui::particleInfo, this));
	MiniGL::addKeyFunc('s', std::bind(&SimulatorBase::saveState, m_simulatorBase, ""));
#ifdef WIN32
	MiniGL::addKeyFunc('l', std::bind(&SimulatorBase::loadStateDialog, m_simulatorBase));
#endif 

	if (MiniGL::checkOpenGLVersion(3, 3))
		Simulator_OpenGL::initShaders(m_simulatorBase->getExePath() + "/resources/shaders");

	const int width = MiniGL::getWidth();
	const int height = MiniGL::getHeight();

	initImgui();
	initImguiParameters();

	MiniGL::addKeyboardFunc([](int key, int scancode, int action, int mods) -> bool { ImGui_ImplGlfw_KeyCallback(MiniGL::getWindow(), key, scancode, action, mods); return ImGui::GetIO().WantCaptureKeyboard; });
	MiniGL::addCharFunc([](int key, int action) -> bool { ImGui_ImplGlfw_CharCallback(MiniGL::getWindow(), key); return ImGui::GetIO().WantCaptureKeyboard; });
	MiniGL::addMousePressFunc([](int button, int action, int mods) -> bool { ImGui_ImplGlfw_MouseButtonCallback(MiniGL::getWindow(), button, action, mods); return ImGui::GetIO().WantCaptureMouse; });
	MiniGL::addMouseWheelFunc([](int pos, double xoffset, double yoffset) -> bool { ImGui_ImplGlfw_ScrollCallback(MiniGL::getWindow(), xoffset, yoffset); return ImGui::GetIO().WantCaptureMouse; });

	MiniGL::setClientIdleFunc(std::bind(&SimulatorBase::timeStep, m_simulatorBase));
	MiniGL::setClientDestroyFunc(std::bind(&Simulator_GUI_imgui::destroy, this));
	MiniGL::addKeyFunc('r', std::bind(&SimulatorBase::reset, m_simulatorBase));
	MiniGL::addKeyFunc('w', Simulator_GUI_imgui::switchDrawMode);
	MiniGL::addKeyFunc(' ', std::bind(&Simulator_GUI_imgui::switchPause, this));
	MiniGL::addKeyFunc('m', std::bind(&SimulatorBase::determineMinMaxOfScalarField, m_simulatorBase));
	MiniGL::setClientSceneFunc(std::bind(&Simulator_GUI_imgui::render, this));
}

void Simulator_GUI_imgui::initImgui()
{
	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();

	ImGuiStyle* style = &ImGui::GetStyle();
	ImVec4* colors = style->Colors;
	colors[ImGuiCol_Text] = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
	colors[ImGuiCol_WindowBg] = ImVec4(0.1f, 0.1f, 0.1f, 0.8f);
	style->FrameBorderSize = 0.5f;
	style->FrameRounding = 3.0f;
	style->TabBorderSize = 1.0f;
	
	std::string font = Utilities::FileSystem::normalizePath(m_simulatorBase->getExePath() + "/resources/fonts/Roboto-Medium.ttf");
	//std::string font = Utilities::FileSystem::normalizePath(m_simulatorBase->getExePath() + "/resources/fonts/DroidSans.ttf");
	io.Fonts->AddFontFromFileTTF(font.c_str(), 15.0f);

	// Setup Platform/Renderer bindings
	ImGui_ImplGlfw_InitForOpenGL(MiniGL::getWindow(), false);
	const char* glsl_version = "#version 330";
	ImGui_ImplOpenGL3_Init(glsl_version);
}

void Simulator_GUI_imgui::initImguiParameters()
{
	imguiParameters::imguiNumericParameter<Real>* timeParam = new imguiParameters::imguiNumericParameter<Real>();
	timeParam->description = "Current simulation time";
	timeParam->label = "Time";
	timeParam->readOnly = true;
	timeParam->getFct = []() -> Real { return TimeManager::getCurrent()->getTime(); };
	imguiParameters::addParam("General", "General", timeParam);

	imguiParameters::imguiNumericParameter<Real>* timeStepSizeParam = new imguiParameters::imguiNumericParameter<Real>();
	timeStepSizeParam->description = "Set time step size";
	timeStepSizeParam->label = "Time step size";
	timeStepSizeParam->minValue = static_cast<Real>(0.00001);
	timeStepSizeParam->maxValue = static_cast<Real>(0.1);
	timeStepSizeParam->getFct = []() -> Real { return TimeManager::getCurrent()->getTimeStepSize(); };
	timeStepSizeParam->setFct = [](Real v) { TimeManager::getCurrent()->setTimeStepSize(v); };
	imguiParameters::addParam("General", "General", timeStepSizeParam);

	imguiParameters::imguiBoolParameter* wireframeParam = new imguiParameters::imguiBoolParameter();
	wireframeParam->description = "Switch wireframe mode";
	wireframeParam->label = "Wireframe";
	wireframeParam->readOnly = false;
	wireframeParam->getFct = []() -> bool { return MiniGL::getDrawMode() == GL_LINE; };
	wireframeParam->setFct = [](bool v) {
		if (!v)
			MiniGL::setDrawMode(GL_FILL);
		else
			MiniGL::setDrawMode(GL_LINE);
	};
	imguiParameters::addParam("Visualization", "", wireframeParam);
}

void Simulator_GUI_imgui::createSimulationParameterGUI()
{
	ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(390, 900), ImGuiCond_FirstUseEver);

	ImGui::Begin("Settings");  
	ImGui::PushItemWidth(175);

	imguiParameters::createParameterGUI();

	ImGui::PopItemWidth();
	ImGui::End();
}

void Simulator_GUI_imgui::initSimulationParameterGUI()
{
	imguiParameters::cleanup();

	initImguiParameters();

	Simulation *sim = Simulation::getCurrent();
	if (m_simulatorBase)
		imguiParameters::createParameterObjectGUI(m_simulatorBase);
	imguiParameters::createParameterObjectGUI(sim);
	imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) sim->getTimeStep());
#ifdef USE_DEBUG_TOOLS
	imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) sim->getDebugTools());
#endif

	// Enum for all fluid models
	if (sim->numberOfFluidModels() > 0)
	{
		FluidModel* model = sim->getFluidModel(m_currentFluidModel);

		// Select fluid model
		{
			imguiParameters::imguiEnumParameter* param = new imguiParameters::imguiEnumParameter();
			param->description = "Select a fluid model to set its parameters below.";
			param->label = "Current fluid model";
			param->readOnly = false;
			for (unsigned int j = 0; j < sim->numberOfFluidModels(); j++)
			{
				param->items.push_back(sim->getFluidModel(j)->getId());
			}
			param->getFct = [this]() -> int { return m_currentFluidModel; };
			param->setFct = [this](int v) { m_currentFluidModel = v; initSimulationParameterGUI(); };
			imguiParameters::addParam("Fluid Model", "", param);
		}

		// show GUI only for currently selected fluid model
		unsigned int i = m_currentFluidModel;

		m_colorFieldNames.clear();
		m_colorFieldNames.resize(model->numberOfFields());

		// Select color field
		{
			imguiParameters::imguiEnumParameter* param = new imguiParameters::imguiEnumParameter();
			param->description = "Choose vector or scalar field for particle coloring.";
			param->label = "Color field";
			param->readOnly = false;
			int idx = 0;
			for (unsigned int j = 0; j < model->numberOfFields(); j++)
			{
				const FieldDescription& field = model->getField(j);
				if ((field.type == FieldType::Scalar) || (field.type == FieldType::Vector3) || 
					(field.type == FieldType::UInt) || (field.type == FieldType::Matrix3) ||
					(field.type == FieldType::Vector6) || (field.type == FieldType::Matrix6))
				{
					param->items.push_back(field.name);
					m_colorFieldNames[idx] = field.name;
					idx++;
				}
			}
			param->getFct = [this]() -> int { 
				const std::string& fieldName = getSimulatorBase()->getColorField(m_currentFluidModel);
				for (auto i = 0; i < m_colorFieldNames.size(); i++)
				{
					if (m_colorFieldNames[i] == fieldName)
						return i;
				}
				return 0;
			};
			param->setFct = [this](int v) { 
				getSimulatorBase()->setColorField(m_currentFluidModel, m_colorFieldNames[v]); 
				getSimulatorBase()->determineMinMaxOfScalarField();
				getSimulatorBase()->updateScalarField();
			};
			imguiParameters::addParam("Fluid Model", model->getId(), param);
		}

		// Select color map type
		{
			imguiParameters::imguiEnumParameter* param = new imguiParameters::imguiEnumParameter();
			param->description = "Choose a color map.";
			param->label = "Color map";
			param->readOnly = false;
			param->items.push_back("None");
			param->items.push_back("Jet");
			param->items.push_back("Plasma");
			param->items.push_back("CoolWarm");
			param->items.push_back("BlueWhiteRed");
			param->items.push_back("Seismic");
			param->getFct = [this]() -> int { return getSimulatorBase()->getColorMapType(m_currentFluidModel); };
			param->setFct = [this](int v) { getSimulatorBase()->setColorMapType(m_currentFluidModel, v); };
			imguiParameters::addParam("Fluid Model", model->getId(), param);
		}

		// Select color min/max value
		{
			imguiParameters::imguiNumericParameter<Real>* param1 = new imguiParameters::imguiNumericParameter<Real>();
			param1->description = "Minimal value used for color-coding the color field in the rendering process.";
			param1->label = "Min. value (shader)";
			param1->getFct = [this]() -> Real { return getSimulatorBase()->getRenderMinValue(m_currentFluidModel); };
			param1->setFct = [this](Real v) { getSimulatorBase()->setRenderMinValue(m_currentFluidModel, v); };
			imguiParameters::addParam("Fluid Model", model->getId(), param1);

			imguiParameters::imguiNumericParameter<Real>* param2 = new imguiParameters::imguiNumericParameter<Real>();
			param2->description = "Maximal value used for color-coding the color field in the rendering process.";
			param2->label = "Max. value (shader)";
			param2->getFct = [this]() -> Real { return getSimulatorBase()->getRenderMaxValue(m_currentFluidModel); };
			param2->setFct = [this](Real v) { getSimulatorBase()->setRenderMaxValue(m_currentFluidModel, v); };
			imguiParameters::addParam("Fluid Model", model->getId(), param2);

			imguiParameters::imguiFunctionParameter* param3 = new imguiParameters::imguiFunctionParameter();
			param3->description = "Recompute min and max values for color-coding the color field in the rendering process.";
			param3->label = "Rescale";
			param3->readOnly = false;
			param3->function = [this]() { getSimulatorBase()->determineMinMaxOfScalarField(); };
			imguiParameters::addParam("Fluid Model", model->getId(), param3);
		}

		imguiParameters::createParameterObjectGUI(model);
		imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getDragBase());
		imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getSurfaceTensionBase());
		imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getViscosityBase());
		imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getVorticityBase());
		imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getElasticityBase());
	}
}

void Simulator_GUI_imgui::initParameterGUI()
{
}

void Simulator_GUI_imgui::update()
{
	// Start the Dear ImGui frame
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();
	
	createSimulationParameterGUI();

	// Rendering
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void Simulator_GUI_imgui::destroy()
{
	// Cleanup
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
}

void Simulator_GUI_imgui::cleanup()
{
	MiniGL::getKeyFunc().clear();
}

void Simulator_GUI_imgui::render()
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

		const FieldDescription* field = nullptr;
		field = &model->getField(base->getColorField(i));

		bool useScalarField = true;
		if ((field == nullptr) || (base->getScalarField().size() == 0))
			useScalarField = false;
		Simulator_OpenGL::renderFluid(model, fluidColor, base->getColorMapType(i),
			useScalarField, base->getScalarField(), base->getRenderMinValue(i), base->getRenderMaxValue(i));
		Simulator_OpenGL::renderSelectedParticles(model, getSelectedParticles(), base->getColorMapType(i),
			base->getRenderMinValue(i), base->getRenderMaxValue(i));
	}
	renderBoundary();
	update();
}

void Simulator_GUI_imgui::renderBoundary()
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

void Simulator_GUI_imgui::reset()
{
	m_selectedParticles.clear();
}

void Simulator_GUI_imgui::selection(const Vector2i &start, const Vector2i &end, void *clientData)
{
	Simulator_GUI_imgui *gui = (Simulator_GUI_imgui*)clientData;
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


void Simulator_GUI_imgui::mouseMove(int x, int y, void *clientData)
{
	Simulator_GUI_imgui *gui = (Simulator_GUI_imgui*)clientData;
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

void Simulator_GUI_imgui::particleInfo()
{
	SimulatorBase::particleInfo(m_selectedParticles);
}

void Simulator_GUI_imgui::run()
{
	MiniGL::mainLoop();
}

void Simulator_GUI_imgui::stop()
{
	MiniGL::leaveMainLoop();
}

void Simulator_GUI_imgui::addKeyFunc(char k, std::function<void()> const& func)
{
	MiniGL::addKeyFunc(k, func);
}

void Simulator_GUI_imgui::switchPause()
{
	m_simulatorBase->setValue(SimulatorBase::PAUSE, !m_simulatorBase->getValue<bool>(SimulatorBase::PAUSE));
}

void Simulator_GUI_imgui::switchDrawMode()
{
	if (MiniGL::getDrawMode() == GL_LINE)
		MiniGL::setDrawMode(GL_FILL);
	else
		MiniGL::setDrawMode(GL_LINE);
}
