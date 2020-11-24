#include "PartioViewer_GUI_imgui.h"
#include "GUI/OpenGL/MiniGL.h"
#include "GUI/imgui/imguiParameters.h"
#include "SPlisHSPlasH/Simulation.h"
#include "../../PartioViewer.h"
#include "../OpenGL/PartioViewer_OpenGL.h"
#include "GUI/OpenGL/Selection.h"
#include "Utilities/FileSystem.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"


using namespace SPH;


std::istream& operator >> (std::istream& istream, Vector3r& v)
{
	return istream >> std::skipws >> v[0] >> v[1] >> v[2];
}


PartioViewer_GUI_imgui::PartioViewer_GUI_imgui(PartioViewer *viewer) :
	PartioViewer_GUI_Base()
{	
	m_viewer = viewer;
	m_simulatorBase = nullptr;
	m_currentFluidModel = 0;
	m_renderWalls = false;
	m_showBBox = false;
	m_camPos = Vector3r(0.0, 3.0, 10.0);
	m_camLookat = Vector3r(0, 0, 0);
}

PartioViewer_GUI_imgui::~PartioViewer_GUI_imgui(void)
{	
	imguiParameters::cleanup();
}

void PartioViewer_GUI_imgui::init()
{
	// OpenGL
	MiniGL::init(m_viewer->getArgc(), m_viewer->getArgv(), m_viewer->getWidth(), m_viewer->getHeight(), "Partio Viewer");
	MiniGL::initLights();
	MiniGL::setViewport(40.0, 0.1f, 500.0, m_camPos, m_camLookat);

	MiniGL::setSelectionFunc(selection, m_viewer);

	MiniGL::addKeyFunc('i', std::bind(&PartioViewer::particleInfo, m_viewer));
	MiniGL::addKeyFunc('+', std::bind(&PartioViewer::nextFrame, m_viewer));
	MiniGL::addKeyFunc('-', std::bind(&PartioViewer::prevFrame, m_viewer));
	MiniGL::addKeyFunc('s', std::bind(&PartioViewer::saveFrame, m_viewer));
	MiniGL::addKeyFunc('v', std::bind(&PartioViewer::generateVideo, m_viewer));
	MiniGL::addKeyFunc('j', std::bind(&PartioViewer::generateSequence, m_viewer));
	MiniGL::addKeyFunc(' ', [&] { m_viewer->setPause(!m_viewer->getPause()); });
	MiniGL::addKeyFunc('r', [&] { m_viewer->reset(); update();  });
	MiniGL::addKeyFunc('m', [&] { determineMinMaxValues(); });

	MiniGL::setClientIdleFunc(std::bind(&PartioViewer::timeStep, m_viewer));	
	MiniGL::setClientDestroyFunc(std::bind(&PartioViewer_GUI_imgui::destroy, this));

	const int width = MiniGL::getWidth();
	const int height = MiniGL::getHeight();

	initImgui();
	initImguiParameters();

	MiniGL::addKeyboardFunc([](int key, int scancode, int action, int mods) -> bool { ImGui_ImplGlfw_KeyCallback(MiniGL::getWindow(), key, scancode, action, mods); return ImGui::GetIO().WantCaptureKeyboard; });
	MiniGL::addCharFunc([](int key, int action) -> bool { ImGui_ImplGlfw_CharCallback(MiniGL::getWindow(), key); return ImGui::GetIO().WantCaptureKeyboard; });
	MiniGL::addMousePressFunc([](int button, int action, int mods) -> bool { ImGui_ImplGlfw_MouseButtonCallback(MiniGL::getWindow(), button, action, mods); return ImGui::GetIO().WantCaptureMouse; });
	MiniGL::addMouseWheelFunc([](int pos, double xoffset, double yoffset) -> bool { ImGui_ImplGlfw_ScrollCallback(MiniGL::getWindow(), xoffset, yoffset); return ImGui::GetIO().WantCaptureMouse; });
}

void PartioViewer_GUI_imgui::run()
{
	MiniGL::mainLoop();
}

void PartioViewer_GUI_imgui::stop()
{
	MiniGL::leaveMainLoop();
}

void PartioViewer_GUI_imgui::addOptions(cxxopts::Options &options)
{
 	options.add_options()
		("camPos", "Camera position (e.g. --camPos \"0 1 5\")", cxxopts::value<Vector3r>()->default_value("0 3 10"))
		("camLookat", "Camera lookat (e.g. --camLookat \"0 0 0\")", cxxopts::value<Vector3r>()->default_value("0 0 0"))
		("showBBox", "Show bounding box.")
 		;
}

void PartioViewer_GUI_imgui::parseOptions(cxxopts::ParseResult &result)
{
	if (result.count("camPos"))
		m_camPos = result["camPos"].as<Vector3r>();

	if (result.count("camLookat"))
		m_camLookat = result["camLookat"].as<Vector3r>();

	m_showBBox = result.count("showBBox") != 0;
}

void PartioViewer_GUI_imgui::initImgui()
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

	std::string font = Utilities::FileSystem::normalizePath(m_viewer->getExePath() + "/resources/fonts/Roboto-Medium.ttf");
	//std::string font = Utilities::FileSystem::normalizePath(m_simulatorBase->getExePath() + "/resources/fonts/DroidSans.ttf");
	io.Fonts->AddFontFromFileTTF(font.c_str(), 15.0f);

	// Setup Platform/Renderer bindings
	ImGui_ImplGlfw_InitForOpenGL(MiniGL::getWindow(), false);
	const char* glsl_version = "#version 330";
	ImGui_ImplOpenGL3_Init(glsl_version);
}


void PartioViewer_GUI_imgui::initImguiParameters()
{
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

void PartioViewer_GUI_imgui::initParameterGUI()
{
	imguiParameters::cleanup();

	initImguiParameters();

	imguiParameters::imguiNumericParameter<unsigned int>* uparam = new imguiParameters::imguiNumericParameter<unsigned int>();
	uparam->description = "Current frame index";
	uparam->label = "Frame index";
	uparam->minValue = 0;
	uparam->getFct = [this]() -> int { return m_viewer->getFrameIndex(); };
	uparam->setFct = [this](int v) { m_viewer->setFrame(v); };
	imguiParameters::addParam("General", "General", uparam);


	imguiParameters::imguiNumericParameter<Real>* rparam = new imguiParameters::imguiNumericParameter<Real>();
	rparam->description = "Particle radius";
	rparam->label = "Particle radius";
	rparam->minValue = 0.001;
	rparam->getFct = [this]() -> Real { return m_viewer->getParticleRadius(); };
	rparam->setFct = [this](Real v) { m_viewer->setParticleRadius(v); };
	imguiParameters::addParam("General", "General", rparam);

	imguiParameters::imguiBoolParameter *bparam = new imguiParameters::imguiBoolParameter();
	bparam->description = "Render walls";
	bparam->label = "Render walls";
	bparam->getFct = [this]() -> bool { return m_renderWalls; };
	bparam->setFct = [this](bool v) { m_renderWalls = v; };
	imguiParameters::addParam("Visualization", "", bparam);

	bparam = new imguiParameters::imguiBoolParameter();
	bparam->description = "Show bounding box";
	bparam->label = "Show bounding box";
	bparam->getFct = [this]() -> bool { return m_showBBox; };
	bparam->setFct = [this](bool v) { m_showBBox = v; };
	imguiParameters::addParam("Visualization", "", bparam);

	bparam = new imguiParameters::imguiBoolParameter();
	bparam->description = "Use plane";
	bparam->label = "Use plane";
	bparam->getFct = [this]() -> bool { return m_viewer->getUsePlane(); };
	bparam->setFct = [this](bool v) { m_viewer->setUsePlane(v); if (v) m_viewer->updateData();  };
	imguiParameters::addParam("Visualization", "", bparam);

	imguiParameters::imguiVec3fParameter* vparam = new imguiParameters::imguiVec3fParameter();
	vparam->description = "Plane point";
	vparam->label = "Plane point";
	vparam->getFct = [this]() -> Eigen::Vector3f { return m_viewer->getPlanePoint(); };
	vparam->setFct = [this](Eigen::Vector3f &v) { m_viewer->setPlanePoint(v); if (m_viewer->getUsePlane()) m_viewer->updateData(); };
	imguiParameters::addParam("Visualization", "", vparam);

	vparam = new imguiParameters::imguiVec3fParameter();
	vparam->description = "Plane normal";
	vparam->label = "Plane normal";
	vparam->getFct = [this]() -> Eigen::Vector3f { return m_viewer->getPlaneNormal(); };
	vparam->setFct = [this](Eigen::Vector3f& v) { m_viewer->setPlaneNormal(v); if (m_viewer->getUsePlane()) m_viewer->updateData(); };
	imguiParameters::addParam("Visualization", "", vparam);
 

	std::vector<SPH::Fluid>& fluids = m_viewer->getFluids();

	if (fluids.size() > 0)
	{
		// Select fluid model
		{
			imguiParameters::imguiEnumParameter* param = new imguiParameters::imguiEnumParameter();
			param->description = "Select a fluid model to set its parameters below.";
			param->label = "Current fluid model";
			for (unsigned int j = 0; j < fluids.size(); j++)
			{
				param->items.push_back(Utilities::FileSystem::getFileName(fluids[j].inputFile));
			}
			param->getFct = [this]() -> int { return m_currentFluidModel; };
			param->setFct = [this](int v) { m_currentFluidModel = v; };
			imguiParameters::addParam("Visualization", "", param);
		}
	
		// show GUI only for currently selected fluid model
		if (fluids[m_currentFluidModel].partioData)
		{
			m_mapColorField2Attr.clear();
			imguiParameters::imguiEnumParameter* eparam = new imguiParameters::imguiEnumParameter();
			eparam->description = "Choose vector or scalar field for particle coloring.";
			eparam->label = "Color field";

			bool found = false;
			int idx = 0;
			for (int i = 0; i < fluids[m_currentFluidModel].partioData->numAttributes(); i++)
			{
				Partio::ParticleAttribute attr;
				fluids[m_currentFluidModel].partioData->attributeInfo(i, attr);
				if ((attr.type == Partio::FLOAT) || (attr.type == Partio::VECTOR) || (attr.type == Partio::INT))
				{
					m_mapColorField2Attr[idx] = i;
					eparam->items.push_back(attr.name);
					idx++;
				}
			}

			eparam->getFct = [this]() -> int { return m_viewer->getFluids()[m_currentFluidModel].m_colorField; };
			eparam->setFct = [this](int v) { 
				m_viewer->getFluids()[m_currentFluidModel].m_colorField = v; 
				determineMinMaxValues();
				PartioViewer_OpenGL::updateScalarField();
			};
			imguiParameters::addParam("Visualization", "", eparam);
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
			param->getFct = [this]() -> int { return m_viewer->getFluids()[m_currentFluidModel].m_colorMapType; };
			param->setFct = [this](int v) { m_viewer->getFluids()[m_currentFluidModel].m_colorMapType = v; };
			imguiParameters::addParam("Visualization", "", param);
		}

		// Select color min/max value
		{
			imguiParameters::imguiNumericParameter<Real>* param1 = new imguiParameters::imguiNumericParameter<Real>();
			param1->description = "Minimal value used for color-coding the color field in the rendering process.";
			param1->label = "Min. value (shader)";
			param1->getFct = [this]() -> Real { return m_viewer->getFluids()[m_currentFluidModel].m_renderMinValue; };
			param1->setFct = [this](Real v) { m_viewer->getFluids()[m_currentFluidModel].m_renderMinValue = v; };
			imguiParameters::addParam("Visualization", "", param1);

			imguiParameters::imguiNumericParameter<Real>* param2 = new imguiParameters::imguiNumericParameter<Real>();
			param2->description = "Maximal value used for color-coding the color field in the rendering process.";
			param2->label = "Max. value (shader)";
			param2->getFct = [this]() -> Real { return  m_viewer->getFluids()[m_currentFluidModel].m_renderMaxValue; };
			param2->setFct = [this](Real v) { m_viewer->getFluids()[m_currentFluidModel].m_renderMaxValue = v; };
			imguiParameters::addParam("Visualization", "", param2);

			imguiParameters::imguiFunctionParameter* param3 = new imguiParameters::imguiFunctionParameter();
			param3->description = "Recompute min and max values for color-coding the color field in the rendering process.";
			param3->label = "Rescale";
			param3->readOnly = false;
			param3->function = [this]() { determineMinMaxValues(); };
			imguiParameters::addParam("Visualization", "", param3);
		}
	}

	imguiParameters::imguiNumericParameter<int>* iparam = new imguiParameters::imguiNumericParameter<int>();
	iparam->description = "Start frame";
	iparam->label = "Start frame";
	iparam->minValue = 0;
	iparam->getFct = [this]() -> int { return m_viewer->getStartFrame(); };
	iparam->setFct = [this](int v) { m_viewer->setStartFrame(v); };
	imguiParameters::addParam("Export", "", iparam);

	iparam = new imguiParameters::imguiNumericParameter<int>();
	iparam->description = "End frame";
	iparam->label = "End frame";
	iparam->minValue = 0;
	iparam->getFct = [this]() -> int { return m_viewer->getEndFrame(); };
	iparam->setFct = [this](int v) { m_viewer->setEndFrame(v); };
	imguiParameters::addParam("Export", "", iparam);

	iparam = new imguiParameters::imguiNumericParameter<int>();
	iparam->description = "FPS";
	iparam->label = "FPS";
	iparam->minValue = 1;
	iparam->getFct = [this]() -> int { return m_viewer->getFPS(); };
	iparam->setFct = [this](int v) { m_viewer->setFPS(v); };
	imguiParameters::addParam("Export", "", iparam);

	if (MiniGL::checkOpenGLVersion(3, 3))
		PartioViewer_OpenGL::initShaders(m_viewer->getExePath() + "/resources/shaders");

	MiniGL::setClientSceneFunc(std::bind(&PartioViewer_GUI_imgui::renderScene, this));
}

void PartioViewer_GUI_imgui::renderScene()
{
	PartioViewer_OpenGL::renderGrid();

	std::vector<SPH::Fluid>& fluids = m_viewer->getFluids();
	for (size_t i = 0; i < fluids.size(); i++)
	{
		float fluidColor[4] = { 0.3f, 0.5f, 0.9f, 1.0f };
		PartioViewer_OpenGL::hsvToRgb(0.61f - 0.1f*i, 0.66f, 0.9f, fluidColor);
		PartioViewer_OpenGL::render(fluids[i], m_viewer->getParticleRadius(), fluidColor, 
			fluids[i].m_colorMapType, m_mapColorField2Attr[fluids[i].m_colorField],
			fluids[i].m_renderMinValue, fluids[i].m_renderMaxValue, m_viewer->getUsePlane());
	}

	float boundaryColor[4] = { 0.4f, 0.4f, 0.4f, 1.0f };
	std::vector<SPH::Boundary>& boundaries = m_viewer->getBoundaries();
	for (size_t i = 0; i < boundaries.size(); i++)
	{
		PartioViewer_OpenGL::render(boundaries[i], m_renderWalls, boundaryColor);
	}
	// Render bounding box - fluid
	float col[4] = { 0.8,0.8,0.8,1 };
	if (m_showBBox)
		PartioViewer_OpenGL::renderAABB(m_viewer->getFluidBoundingBox(), col);

	update();
}

void PartioViewer_GUI_imgui::determineMinMaxValues()
{
	if (m_viewer->getFluids().size() > 0)
	{
		auto &fluid = m_viewer->getFluids()[m_currentFluidModel];
		const unsigned int nParticles = (unsigned int)fluid.partioData->numParticles();

		Partio::ParticleAttribute attr;
		fluid.partioData->attributeInfo(fluid.m_colorField, attr);

		float minValue = FLT_MAX;
		float maxValue = FLT_MIN;
		for (unsigned int i = 0u; i < nParticles; i++)
		{
			if (attr.type == Partio::VECTOR)
			{
				const Eigen::Map<const Eigen::Vector3f> vec(fluid.partioData->data<float>(attr, i));
				minValue = std::min(minValue, vec.norm());
				maxValue = std::max(maxValue, vec.norm());
			}
			else if (attr.type == Partio::FLOAT)
			{
				minValue = std::min(minValue, *fluid.partioData->data<float>(attr, i));
				maxValue = std::max(maxValue, *fluid.partioData->data<float>(attr, i));
			}
			else if (attr.type == Partio::INT)
			{
				minValue = std::min(minValue, static_cast<float>(*fluid.partioData->data<int>(attr, i)));
				maxValue = std::max(maxValue, static_cast<float>(*fluid.partioData->data<int>(attr, i)));
			}
		}
		fluid.m_renderMinValue = minValue;
		fluid.m_renderMaxValue = maxValue;
	}
}

void PartioViewer_GUI_imgui::render()
{
	MiniGL::viewport();
	renderScene();
	MiniGL::swapBuffers();
}

void PartioViewer_GUI_imgui::update()
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

void PartioViewer_GUI_imgui::cleanup()
{
	MiniGL::getKeyFunc().clear();
}

void PartioViewer_GUI_imgui::destroy()
{
	// Cleanup
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
}


void PartioViewer_GUI_imgui::selection(const Vector2i &start, const Vector2i &end, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	std::vector<SPH::Fluid>& fluids = viewer->getFluids();
	for (size_t i = 0; i < fluids.size(); i++)
	{
		size_t nParticles;
		if (viewer->getUsePlane())
			nParticles = fluids[i].visibleParticles.size();
		else
			nParticles = fluids[i].partioData->numParticles();

		if (nParticles != 0)
		{
			Partio::ParticleAttribute posAttr;
			fluids[i].partioData->attributeInfo(fluids[i].posIndex, posAttr);
			const float* partioX = fluids[i].partioData->data<float>(posAttr, 0);
			std::vector<Vector3r> x;
			x.resize(nParticles);
			for (int j = 0; j < nParticles; j++)
			{
				if (viewer->getUsePlane())
					x[j] = Eigen::Map<const Eigen::Vector3f>(&partioX[3 * fluids[i].visibleParticles[j]]).cast<Real>();
				else
					x[j] = Eigen::Map<const Eigen::Vector3f>(&partioX[3 * j]).cast<Real>();
			}

			std::vector<unsigned int> hits;
			fluids[i].selectedParticles.clear();
			Selection::selectRect(start, end, x.begin(),
				x.end(),
				fluids[i].selectedParticles);

			for (int j = 0; j < fluids[i].selectedParticles.size(); j++)
			{
				if (viewer->getUsePlane())
					fluids[i].selectedParticles[j] = fluids[i].visibleParticles[fluids[i].selectedParticles[j]];
				else
					fluids[i].selectedParticles[j] = fluids[i].selectedParticles[j];
			}
		}
	}
}

void PartioViewer_GUI_imgui::createSimulationParameterGUI()
{
	ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(390, 900), ImGuiCond_FirstUseEver);

	ImGui::Begin("Settings");
	ImGui::PushItemWidth(175);

	imguiParameters::createParameterGUI();

	ImGui::PopItemWidth();
	ImGui::End();
}

unsigned int PartioViewer_GUI_imgui::getWidth() const
{
	return MiniGL::getWidth();
}

unsigned int PartioViewer_GUI_imgui::getHeight() const
{
	return MiniGL::getHeight();
}
