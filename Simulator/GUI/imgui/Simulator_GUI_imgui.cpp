#include "Simulator_GUI_imgui.h"
#include "GUI/OpenGL/MiniGL.h"
#include "GUI/imgui/imguiParameters.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "../OpenGL/Simulator_OpenGL.h"
#include "SPlisHSPlasH/Utilities/SceneLoader.h"
#include "GUI/OpenGL/Selection.h"
#include "Utilities/FileSystem.h"
#include "Simulator/SceneConfiguration.h"
#include "LogWindow.h"

#include "imgui.h"
#include "imgui_internal.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>


using namespace SPH;
using namespace Utilities;

Simulator_GUI_imgui::Simulator_GUI_imgui(SimulatorBase *base) :
	Simulator_GUI_Base(base)
{	
	m_currentFluidModel = 0;
	m_currentScaleIndex = 0;
	m_vsync = false;
	m_iniFound = false;
	m_showLogWindow = true;
}

Simulator_GUI_imgui::~Simulator_GUI_imgui(void)
{	
	imguiParameters::cleanup();
	delete m_logWindow;
}

void Simulator_GUI_imgui::init(const char *name)
{
	m_logWindow = new LogWindow();

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	m_context = ImGui::CreateContext();

	// Add .ini handle for UserData type
	ImGuiSettingsHandler ini_handler;
	ini_handler.TypeName = "SPHSimulator";
	ini_handler.TypeHash = ImHashStr("SPHSimulator");
	ini_handler.ReadOpenFn = readOpenIni;
	ini_handler.ReadLineFn = readIni;
	ini_handler.WriteAllFn = writeIni;
	ini_handler.ApplyAllFn = applySettings;
	ini_handler.UserData = this;
	m_context->SettingsHandlers.push_back(ini_handler);

	// load ini file before window is created
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;           // Enable Docking
	ImGui::LoadIniSettingsFromDisk(io.IniFilename);

	MiniGL::init(1280, 960, name, m_userSettings.vsync, m_userSettings.maximized);
	MiniGL::initLights();

	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
	const bool sim2D = scene.sim2D;
	if (sim2D)
		MiniGL::setViewport(40.0, 0.1f, 500.0, m_simulatorBase->getCameraPosition(), m_simulatorBase->getCameraLookAt());
	else
		MiniGL::setViewport(40.0, 0.1f, 500.0, m_simulatorBase->getCameraPosition(), m_simulatorBase->getCameraLookAt());
	MiniGL::setSelectionFunc(selection, this);
	MiniGL::addKeyFunc(GLFW_KEY_I, 0, std::bind(&Simulator_GUI_imgui::particleInfo, this));
	MiniGL::addKeyFunc(GLFW_KEY_S, GLFW_MOD_CONTROL, std::bind(&SimulatorBase::saveState, m_simulatorBase, ""));
#ifdef USE_NFD_FILE_DIALOG
	MiniGL::addKeyFunc(GLFW_KEY_L, GLFW_MOD_CONTROL, std::bind(&SimulatorBase::loadStateDialog, m_simulatorBase));
#endif
	MiniGL::addKeyFunc(GLFW_KEY_W, GLFW_MOD_CONTROL, std::bind(&SimulatorBase::writeScene, m_simulatorBase));
	MiniGL::addKeyFunc(GLFW_KEY_KP_ADD, 0, std::bind(&SimulatorBase::singleTimeStep, m_simulatorBase));

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
	MiniGL::addKeyFunc(GLFW_KEY_R, 0, std::bind(&SimulatorBase::reset, m_simulatorBase));
	MiniGL::addKeyFunc(GLFW_KEY_W, 0, Simulator_GUI_imgui::switchDrawMode);
	MiniGL::addKeyFunc(GLFW_KEY_SPACE, 0, std::bind(&Simulator_GUI_imgui::switchPause, this));
	MiniGL::addKeyFunc(GLFW_KEY_M, 0, std::bind(&SimulatorBase::determineMinMaxOfScalarField, m_simulatorBase));
	MiniGL::setClientSceneFunc(std::bind(&Simulator_GUI_imgui::render, this));

	// apply user settings from ini file 
	if (m_iniFound)
	{
		MiniGL::setWindowPos(m_userSettings.win_x, m_userSettings.win_y);
		MiniGL::setWindowSize(m_userSettings.win_width, m_userSettings.win_height);
	}
}

void Simulator_GUI_imgui::initStyle()
{
	m_context->Style = ImGuiStyle();

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();

	ImGuiStyle* style = &ImGui::GetStyle();
	ImVec4* colors = style->Colors;
	colors[ImGuiCol_Text] = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
	colors[ImGuiCol_WindowBg] = ImVec4(0.1f, 0.1f, 0.1f, 1.0f);
	style->FrameBorderSize = 0.5f;
	style->FrameRounding = 3.0f;
	style->TabBorderSize = 1.0f;
	style->WindowRounding = 6.0f;
	style->ScaleAllSizes(1);
}

void* Simulator_GUI_imgui::readOpenIni(ImGuiContext* ctx, ImGuiSettingsHandler* handler, const char* name)
{
	Simulator_GUI_imgui* gui = (Simulator_GUI_imgui*)handler->UserData;
	return &gui->m_userSettings;
}

void Simulator_GUI_imgui::readIni(ImGuiContext* ctx, ImGuiSettingsHandler* handler, void* entry, const char* line)
{ 
	Simulator_GUI_imgui* gui = (Simulator_GUI_imgui*)handler->UserData;
	UserSettings* settings = (UserSettings*)entry;
	int x, y, w, h;
	int i;
	if (sscanf(line, "pos=%d,%d", &x, &y) == 2) { settings->win_x = x; settings->win_y = y; }
	else if (sscanf(line, "size=%d,%d", &w, &h) == 2) { settings->win_width = w; settings->win_height = h; }
	else if (sscanf(line, "scale=%d", &i) == 1) { settings->scaleIndex = i; }
	else if (sscanf(line, "maximized=%d", &i) == 1) { settings->maximized = (i != 0); }
	else if (sscanf(line, "vsync=%d", &i) == 1) { settings->vsync = (i != 0); }
	else if (sscanf(line, "show_log_window=%d", &i) == 1) { settings->show_log_window = (i != 0); }
	else if (sscanf(line, "log_filter=%d", &i) == 1) { settings->log_filter = i; }
}

void Simulator_GUI_imgui::writeIni(ImGuiContext* ctx, ImGuiSettingsHandler* handler, ImGuiTextBuffer* out_buf)
{ 
	Simulator_GUI_imgui* gui = (Simulator_GUI_imgui*)handler->UserData;
	out_buf->reserve(out_buf->size() + 200);
	out_buf->appendf("[%s][%s]\n", handler->TypeName, "Settings");
	out_buf->appendf("scale=%d\n", gui->m_currentScaleIndex);

	int x, y;
	MiniGL::getWindowPos(x, y);
	out_buf->appendf("pos=%d,%d\n", x, y);

	int w, h;
	MiniGL::getWindowSize(w, h);
	out_buf->appendf("size=%d,%d\n", w, h);
	out_buf->appendf("maximized=%d\n", MiniGL::getWindowMaximized());

	out_buf->appendf("vsync=%d\n", gui->m_vsync);
	out_buf->appendf("show_log_window=%d\n", gui->m_showLogWindow);
	out_buf->appendf("log_filter=%d\n", gui->m_logWindow->getSelectedFilter());
}

void Simulator_GUI_imgui::applySettings(ImGuiContext* ctx, ImGuiSettingsHandler* handler)
{
	Simulator_GUI_imgui* gui = (Simulator_GUI_imgui*)handler->UserData;
	UserSettings* settings = (UserSettings*) &gui->m_userSettings;
	gui->m_currentScaleIndex = settings->scaleIndex;
	gui->m_vsync = settings->vsync;
	gui->m_showLogWindow = settings->show_log_window;
	gui->m_iniFound = true;
	gui->m_logWindow->setSelectedFilter(settings->log_filter);
}

void Simulator_GUI_imgui::initImgui()
{
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	
	std::string font = Utilities::FileSystem::normalizePath(m_simulatorBase->getExePath() + "/resources/fonts/Roboto-Medium.ttf");
	std::string font2 = Utilities::FileSystem::normalizePath(m_simulatorBase->getExePath() + "/resources/fonts/Cousine-Regular.ttf");

	m_scales.push_back(1.0f);
	m_scales.push_back(1.25f);
	m_scales.push_back(1.5f);
	m_scales.push_back(1.75f);
	m_scales.push_back(2.0f);

	for(int i=0; i < 5; i++)
		m_fonts.push_back(io.Fonts->AddFontFromFileTTF(font.c_str(), m_baseSize * m_scales[i]));
	for (int i = 0; i < 5; i++)
		m_fonts2.push_back(io.Fonts->AddFontFromFileTTF(font2.c_str(), m_baseSize * m_scales[i]));

	initStyle();

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

bool Simulator_GUI_imgui::alignedButton(const char* label, float alignment)
{
	ImGuiStyle& style = ImGui::GetStyle();

	const float size = ImGui::CalcTextSize(label).x + style.FramePadding.x * 2.0f;
	const float avail = ImGui::GetContentRegionAvail().x;

	const float offset = (avail - size) * alignment;
	if (offset > 0.0f)
		ImGui::SetCursorPosX(ImGui::GetCursorPosX() + offset);

	return ImGui::Button(label);
}

void Simulator_GUI_imgui::createMenuBar()
{
	bool openpopup = false;

	ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0.0f, 10.0f));

	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("Write scene file", "Ctrl+w"))
				getSimulatorBase()->writeScene();
			ImGui::Separator();
			if (ImGui::MenuItem("Save state", "Ctrl+s"))
				getSimulatorBase()->saveState();
			if (ImGui::MenuItem("Load state", "Ctrl+l"))
				getSimulatorBase()->loadStateDialog();
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("Simulation"))
		{
			if (ImGui::MenuItem("Pause/run simulation", "Space"))
				switchPause();
			if (ImGui::MenuItem("Reset simulation", "r"))
				getSimulatorBase()->reset();
			if (ImGui::MenuItem("Single step", "+"))
				getSimulatorBase()->singleTimeStep();
			ImGui::Separator();
			if (ImGui::MenuItem("Print particle info", "i", false, (m_selectedParticles.size() > 0)))
				particleInfo();
			if (ImGui::MenuItem("Rescale scalar field", "m"))
				getSimulatorBase()->determineMinMaxOfScalarField();
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("GUI"))
		{
			if (ImGui::MenuItem("V-Sync", "", m_vsync))
			{
				m_vsync = !m_vsync;
				openpopup = true;
			}
			if (ImGui::MenuItem("Scale - 100%", "", m_currentScaleIndex == 0))
			{
				m_currentScaleIndex = 0;
				initStyle();
			}
			if (ImGui::MenuItem("Scale - 125%", "", m_currentScaleIndex == 1))
			{
				m_currentScaleIndex = 1;
				initStyle();
			}
			if (ImGui::MenuItem("Scale - 150%", "", m_currentScaleIndex == 2))
			{
				m_currentScaleIndex = 2;
				initStyle();
			}
			if (ImGui::MenuItem("Scale - 175%", "", m_currentScaleIndex == 3))
			{
				m_currentScaleIndex = 3;
				initStyle();
			}
			if (ImGui::MenuItem("Scale - 200%", "", m_currentScaleIndex == 4))
			{
				m_currentScaleIndex = 4;
				initStyle();
			}
			ImGui::Separator();
			if (ImGui::MenuItem("Show log window", "", m_showLogWindow))
			{
				m_showLogWindow = !m_showLogWindow;
			}
			ImGui::EndMenu();
		}
		ImGui::EndMainMenuBar();
	}

	ImGui::PopStyleVar(1);

	if (openpopup)
	{
		ImGui::OpenPopup("Info");
		openpopup = false;
	}
	bool open = true;
	ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(5.0f, 5.0f));
	ImGui::PushStyleVar(ImGuiStyleVar_WindowTitleAlign, ImVec2(0.5f, 0.5f));
	ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(5.0f, 15.0f));
	ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(10.0f, 10.0f));
	if (ImGui::BeginPopupModal("Info", &open, ImGuiWindowFlags_AlwaysAutoResize))
	{
		ImGui::Text("To turn on/off the vertical sync, \nyou have to restart the simulator.");
		if (alignedButton("Close"))
			ImGui::CloseCurrentPopup();
		ImGui::EndPopup();
	}
	ImGui::PopStyleVar(4);
}

void Simulator_GUI_imgui::createSimulationParameterGUI()
{
	ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(390, 900), ImGuiCond_FirstUseEver);

	float alpha = 0.8f;
	if (ImGui::IsWindowDocked())
		alpha = 1.0f;

	ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.15f, 0.15f, 0.15f, alpha));

	ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(5.0f, 10.0f));
	ImGui::Begin("Settings");

	ImGui::PushItemWidth(175 * m_scales[m_currentScaleIndex]);

	imguiParameters::createParameterGUI();

	ImGui::PopItemWidth();
	ImGui::End();
	ImGui::PopStyleVar(1);
	ImGui::PopStyleColor(1);

	if (m_showLogWindow)
		m_logWindow->drawWindow(m_fonts2[m_currentScaleIndex]);
}

void Simulator_GUI_imgui::initSimulationParameterGUI()
{
	imguiParameters::cleanup();

	Simulation *sim = Simulation::getCurrent();

	imguiParameters::createParameterObjectGUI(sim);
	imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) sim->getTimeStep());

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

			imguiParameters::imguiBoolParameter* param4 = new imguiParameters::imguiBoolParameter();
			param4->description = "Define if fluid phase is rendered.";
			param4->label = "Visible";
			param4->getFct = [this]() -> bool { return getSimulatorBase()->getVisible(m_currentFluidModel); };
			param4->setFct = [this](bool v) { getSimulatorBase()->setVisible(m_currentFluidModel, v); };
			imguiParameters::addParam("Fluid Model", model->getId(), param4);
		}

		imguiParameters::createParameterObjectGUI(model);
		imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getXSPH());
		imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getDragBase());
		imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getSurfaceTensionBase());
		imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getViscosityBase());
		imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getVorticityBase());
		imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getElasticityBase());
	}
	
	initImguiParameters();

	if (m_simulatorBase)
	{
		imguiParameters::createParameterObjectGUI(m_simulatorBase);
#ifdef USE_EMBEDDED_PYTHON
		if (m_simulatorBase->getScriptObject())
			imguiParameters::createParameterObjectGUI(m_simulatorBase->getScriptObject());
#endif
	}

#ifdef USE_DEBUG_TOOLS
	imguiParameters::createParameterObjectGUI((GenParam::ParameterObject*)sim->getDebugTools());
#endif
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

	// init dock space
	static ImGuiDockNodeFlags dockspaceFlags = ImGuiDockNodeFlags_PassthruCentralNode | ImGuiDockNodeFlags_NoDockingInCentralNode | ImGuiDockNodeFlags_AutoHideTabBar;
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;
	ImGuiViewport* viewport = ImGui::GetMainViewport();

	// add enough space for the menubar
	ImVec2 pos = viewport->Pos;
	ImVec2 size = viewport->Size;
	size.y -= m_baseSize * m_scales[m_currentScaleIndex];
	pos.y += m_baseSize * m_scales[m_currentScaleIndex];

	ImGui::SetNextWindowPos(pos);
	ImGui::SetNextWindowSize(size);
	ImGui::SetNextWindowViewport(viewport->ID);
	windowFlags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
	windowFlags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;
	if (dockspaceFlags & ImGuiDockNodeFlags_PassthruCentralNode) 
		windowFlags |= ImGuiWindowFlags_NoBackground;

	ImGui::PushFont(m_fonts[m_currentScaleIndex]);
	ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
	ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
	ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
	ImGui::Begin("DockSpace", nullptr, windowFlags);

	ImGuiID dockspaceID = ImGui::GetID("DockSpace");
	ImGui::DockSpace(dockspaceID, ImVec2(0.0f, 0.0f), dockspaceFlags);
	static int first = true;
	if (!m_iniFound && first)
	{
		first = false;
		ImGui::DockBuilderRemoveNode(dockspaceID);
		ImGui::DockBuilderAddNode(dockspaceID, dockspaceFlags | ImGuiDockNodeFlags_DockSpace);
		ImGui::DockBuilderSetNodeSize(dockspaceID, viewport->Size);

		ImGuiID dock_id_top = -1;
		auto dock_id_bottom = ImGui::DockBuilderSplitNode(dockspaceID, ImGuiDir_Down, 0.3f, nullptr, &dock_id_top);
		ImGuiID dock_id_right = -1;
		auto dock_id_left = ImGui::DockBuilderSplitNode(dock_id_top, ImGuiDir_Left, 0.3f, nullptr, &dock_id_right);
		ImGui::DockBuilderDockWindow("Settings", dock_id_left);
		ImGui::DockBuilderDockWindow("Log", dock_id_bottom);

		ImGui::DockBuilderFinish(dockspaceID);
	}

	ImGui::End();
	ImGui::PopStyleVar(3);

	createMenuBar();
	
	createSimulationParameterGUI();
	ImGui::PopFont();

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

		if (!base->getVisible(i))
			continue;

		const FieldDescription* field = nullptr;
		field = &model->getField(base->getColorField(i));

		bool useScalarField = true;
		if ((field == nullptr) || (base->getScalarField(i).size() == 0))
			useScalarField = false;
		Simulator_OpenGL::renderFluid(model, fluidColor, base->getColorMapType(i),
			useScalarField, base->getScalarField(i), base->getRenderMinValue(i), base->getRenderMaxValue(i));
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
	const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
	const int renderWalls = base->getValue<int>(SimulatorBase::RENDER_WALLS);

	if (((renderWalls == 1) || (renderWalls == 2)) &&
		(sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012))
	{
		for (int body = sim->numberOfBoundaryModels() - 1; body >= 0; body--)
		{
			if ((renderWalls == 1) || (!scene.boundaryModels[body]->isWall))
			{
				BoundaryModel_Akinci2012 *bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(body));
				Eigen::Matrix<float, 4, 1, Eigen::DontAlign> col = scene.boundaryModels[body]->color.cast<float>();
				Simulator_OpenGL::renderBoundaryParticles(bm, col.data(), base->getRenderMinValue(0), base->getRenderMaxValue(0));
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
				Eigen::Matrix<float, 4, 1, Eigen::DontAlign> col = scene.boundaryModels[body]->color.cast<float>();
				Simulator_OpenGL::renderBoundary(bm, col.data());
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
	if (m_selectedParticles.size() > 0)
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

void Simulator_GUI_imgui::addKeyFunc(int key, int modifiers, std::function<void()> const& func)
{
	MiniGL::addKeyFunc(key, modifiers, func);
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
