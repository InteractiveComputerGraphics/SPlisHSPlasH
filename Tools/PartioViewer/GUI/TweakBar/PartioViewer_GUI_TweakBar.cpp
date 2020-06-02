#include "PartioViewer_GUI_TweakBar.h"
#include "GUI/OpenGL/MiniGL.h"
#include "GUI/TweakBar/TweakBarParameters.h"
#include "SPlisHSPlasH/Simulation.h"
#include "../../PartioViewer.h"
#include "../OpenGL/PartioViewer_OpenGL.h"
#include "GUI/OpenGL/Selection.h"
#include "Utilities/FileSystem.h"

using namespace SPH;


std::istream& operator >> (std::istream& istream, Vector3r& v)
{
	return istream >> std::skipws >> v[0] >> v[1] >> v[2];
}


PartioViewer_GUI_TweakBar::PartioViewer_GUI_TweakBar(PartioViewer *viewer) :
	PartioViewer_GUI_Base()
{	
	m_viewer = viewer;
	m_tweakBar = nullptr;
	m_simulatorBase = nullptr;
	m_currentFluidModel = 0;
	m_renderWalls = false;
	m_showBBox = false;
	m_camPos = Vector3r(0.0, 3.0, 10.0);
	m_camLookat = Vector3r(0, 0, 0);
}

PartioViewer_GUI_TweakBar::~PartioViewer_GUI_TweakBar(void)
{	
}

void PartioViewer_GUI_TweakBar::init()
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

	MiniGL::setClientIdleFunc(std::bind(&PartioViewer::timeStep, m_viewer));	

	const int width = MiniGL::getWidth();
	const int height = MiniGL::getHeight();

	// Initialize AntTweakBar
	// (note that AntTweakBar could also be initialized after GLUT, no matter)
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
}

void PartioViewer_GUI_TweakBar::run()
{
	MiniGL::mainLoop();
}

void PartioViewer_GUI_TweakBar::stop()
{
	MiniGL::leaveMainLoop();
}

void PartioViewer_GUI_TweakBar::addOptions(cxxopts::Options &options)
{
 	options.add_options()
		("camPos", "Camera position (e.g. --camPos \"0 1 5\")", cxxopts::value<Vector3r>()->default_value("0 3 10"))
		("camLookat", "Camera lookat (e.g. --camLookat \"0 0 0\")", cxxopts::value<Vector3r>()->default_value("0 0 0"))
 		;
}

void PartioViewer_GUI_TweakBar::parseOptions(cxxopts::ParseResult &result)
{
	if (result.count("camPos"))
		m_camPos = result["camPos"].as<Vector3r>();

	if (result.count("camLookat"))
		m_camLookat = result["camLookat"].as<Vector3r>();
}

void PartioViewer_GUI_TweakBar::initTweakBar()
{
	// Create a tweak bar
	m_tweakBar = TwNewBar("TweakBar");
	TwDefine(" GLOBAL help='MiniGL TweakBar.' "); // Message added to the help bar.
	//TwDefine(" TweakBar size='300 900' valueswidth=120 position='5 5' color='96 200 224' text=dark "); // change default tweak bar size and color
	TwDefine(" TweakBar size='300 900' valueswidth=120 position='5 5'"); // change default tweak bar size and color

	initTweakBarParameters();
}


void PartioViewer_GUI_TweakBar::initTweakBarParameters()
{
	// Add callback to toggle auto-rotate mode (callback functions are defined above).
	TwAddVarCB(getTweakBar(), "Wireframe", TW_TYPE_BOOLCPP, setWireframeCB, getWireframeCB, nullptr,
		" label='Wireframe' key=w help='Toggle wireframe mode.' group='Visualization' ");

	TwAddVarCB(getTweakBar(), "Rotation", TW_TYPE_QUAT4F, setRotationCB, getRotationCB, nullptr,
		" label='Rotation' open help='Change the rotation.' group='Visualization' ");
}

void PartioViewer_GUI_TweakBar::initParameterGUI()
{
	TwRemoveAllVars(getTweakBar());

	initTweakBarParameters();

	TwAddVarCB(getTweakBar(), "frameIndex", TW_TYPE_UINT32, setFrameIndex, getFrameIndex, m_viewer, " label='Frame index' min=0 group=General");
	TwAddVarCB(getTweakBar(), "particleRadius", TW_TYPE_REAL, setParticleRadius, getParticleRadius, m_viewer, " label='Particle radius' min=0.001 group=General");
	
	TwAddVarRW(getTweakBar(), "renderWalls", TW_TYPE_BOOLCPP, &m_renderWalls, " label='Render walls' group=Visualization");	
	TwAddVarRW(getTweakBar(), "showBBox", TW_TYPE_BOOLCPP, &m_showBBox, " label='Show bounding box' group=Visualization");

	TwAddVarCB(getTweakBar(), "usePlane", TW_TYPE_BOOLCPP, setUsePlane, getUsePlane, m_viewer, " label='Use plane' group=Visualization");
	TwAddVarCB(getTweakBar(), "planePoint", TW_TYPE_DIR3F, setPlanePoint, getPlanePoint, m_viewer, " label='Plane point' group=Visualization");
	TwAddVarCB(getTweakBar(), "planeNormal", TW_TYPE_DIR3F, setPlaneNormal, getPlaneNormal, m_viewer, " label='Plane normal' group=Visualization");
	  	
	    
	std::vector<SPH::Fluid>& fluids = m_viewer->getFluids();
	if (fluids.size() > 0)
	{
		// Select fluid model
		{
			TwType enumType = TwDefineEnum("CurrentFluidModel", NULL, 0);
			std::ostringstream oss;
			int idx = 0;
			for (unsigned int j = 0; j < fluids.size(); j++)
			{
				if (idx != 0)
					oss << ", ";
				oss << idx << " {" << Utilities::FileSystem::getFileName(fluids[j].inputFile) << "}";
				idx++;
			}
			std::string enumStr = " label='Current fluid model' enum='" + oss.str() + "' group='Visualization'";
			enumStr = enumStr + " help='Select a fluid model to set its parameters below.'";
			TwAddVarCB(getTweakBar(), "CurrentFluidModel", enumType, setCurrentFluidModel, getCurrentFluidModel, this, enumStr.c_str());
		}

		// show GUI only for currently selected fluid model
		if (fluids[m_currentFluidModel].partioData)
		{
			TwType enumType = TwDefineEnum("ColorFieldEnum", NULL, 0);
			std::ostringstream oss;
			m_mapColorField2Attr.clear();
			int idx = 0;
			for (int i = 0; i < fluids[m_currentFluidModel].partioData->numAttributes(); i++)
			{
				Partio::ParticleAttribute attr;
				fluids[m_currentFluidModel].partioData->attributeInfo(i, attr);
				if ((attr.type == Partio::FLOAT) || (attr.type == Partio::VECTOR))
				{
					m_mapColorField2Attr[idx] = i;
					if (idx != 0)
						oss << ", ";
					oss << idx << " {" << attr.name.c_str() << "}";
					idx++;
				}
			}

			std::string enumStr = " label='Color field' enum='" + oss.str() + "' group='Visualization'";
			enumStr = enumStr + " help='Choose vector or scalar field for particle coloring.'";
			TwAddVarRW(getTweakBar(), "ColorField", enumType, &fluids[m_currentFluidModel].m_colorField, enumStr.c_str());
		}
	}
		  
	TwType enumType2 = TwDefineEnum("ColorMapTypeEnum", NULL, 0);
	std::string str = " label='Color map' enum='0 {None}, 1 {Jet}, 2 {Plasma}' group='Visualization' help='Choose a color map.'";
	TwAddVarRW(getTweakBar(), "ColorMapType", enumType2, &fluids[m_currentFluidModel].m_colorMapType, str.c_str());
	str = " label='Min. value (shader)' step=0.001 precision=3 group='Visualization' help='Minimal value used for color-coding the color field in the rendering process.'";
	TwAddVarRW(getTweakBar(), "RenderMinValue", TW_TYPE_FLOAT, &fluids[m_currentFluidModel].m_renderMinValue, str.c_str());
	str = " label='Max. value (shader)' step=0.001 precision=3 group='Visualization' help='Maximal value used for color-coding the color field in the rendering process.'";
	TwAddVarRW(getTweakBar(), "RenderMaxValue", TW_TYPE_FLOAT, &fluids[m_currentFluidModel].m_renderMaxValue, str.c_str());
	  
	TwAddVarCB(getTweakBar(), "StartFrame", TW_TYPE_INT32, setStartFrame, getStartFrame, m_viewer, " label='Start frame' min=0 group=Export");
	TwAddVarCB(getTweakBar(), "EndFrame", TW_TYPE_INT32, setEndFrame, getEndFrame, m_viewer, " label='End frame' min=0 group=Export");
	TwAddVarCB(getTweakBar(), "FPS", TW_TYPE_INT32, setFPS, getFPS, m_viewer, " label='FPS' min=1 group=Export");

	if (MiniGL::checkOpenGLVersion(3, 3))
		PartioViewer_OpenGL::initShaders(m_viewer->getExePath() + "/resources/shaders");

	MiniGL::setClientSceneFunc(std::bind(&PartioViewer_GUI_TweakBar::renderScene, this));
}

void PartioViewer_GUI_TweakBar::renderScene()
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

void PartioViewer_GUI_TweakBar::render()
{
	MiniGL::viewport();
	renderScene();
	MiniGL::swapBuffers();
}

void PartioViewer_GUI_TweakBar::update()
{
	TwRefreshBar(getTweakBar());
	TwDraw();
}

void PartioViewer_GUI_TweakBar::cleanup()
{
	TwDeleteBar(getTweakBar());
	TwTerminate();
}


void PartioViewer_GUI_TweakBar::selection(const Vector2i &start, const Vector2i &end, void *clientData)
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


void TW_CALL PartioViewer_GUI_TweakBar::setWireframeCB(const void *value, void *clientData)
{
	const int val = *(const int *)(value);
	if (val == 0)
		MiniGL::setDrawMode(GL_FILL);
	else
		MiniGL::setDrawMode(GL_LINE);
}

void TW_CALL PartioViewer_GUI_TweakBar::getWireframeCB(void *value, void *clientData)
{
	*(int *)(value) = MiniGL::getDrawMode() == GL_LINE;
}

void TW_CALL PartioViewer_GUI_TweakBar::setRotationCB(const void *value, void *clientData)
{
	const float *val = (const float *)(value);

	Quaternionr q;
	q.x() = (Real)val[0];
	q.y() = (Real)val[1];
	q.z() = (Real)val[2];
	q.w() = -(Real)val[3];
	MiniGL::setRotation(q);
}

void TW_CALL PartioViewer_GUI_TweakBar::getRotationCB(void *value, void *clientData)
{
	Quaternionr q = MiniGL::getRotation();
	float *val = (float*)(value);
	val[0] = (float)q.x();
	val[1] = (float)q.y();
	val[2] = (float)q.z();
	val[3] = -(float)q.w();
}

void TW_CALL PartioViewer_GUI_TweakBar::setFrameIndex(const void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	const unsigned int val = *(const unsigned int *)(value);
	viewer->setFrame(val);
}

void TW_CALL PartioViewer_GUI_TweakBar::getFrameIndex(void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	*(unsigned int *)(value) = viewer->getFrameIndex();
}

void TW_CALL PartioViewer_GUI_TweakBar::setUsePlane(const void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	const bool val = *(const bool *)(value);
	viewer->setUsePlane(val);
	if (val)
		viewer->updateData();
}

void TW_CALL PartioViewer_GUI_TweakBar::getUsePlane(void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	*(bool *)(value) = viewer->getUsePlane();
}

void TW_CALL PartioViewer_GUI_TweakBar::setPlaneNormal(const void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	float * const val = (float* const)(value);
	viewer->setPlaneNormal(Vector3f(val[0], val[1], val[2]));
	if (viewer->getUsePlane())
		viewer->updateData();
}

void TW_CALL PartioViewer_GUI_TweakBar::getPlaneNormal(void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	auto planeNormal = viewer->getPlaneNormal();
	((float*)value)[0] = planeNormal[0];
	((float*)value)[1] = planeNormal[1];
	((float*)value)[2] = planeNormal[2];
}

void TW_CALL PartioViewer_GUI_TweakBar::setPlanePoint(const void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	float * const val = (float* const)(value);
	viewer->setPlanePoint(Vector3f(val[0], val[1], val[2]));
	if (viewer->getUsePlane())
		viewer->updateData();
}

void TW_CALL PartioViewer_GUI_TweakBar::getPlanePoint(void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	auto planePoint = viewer->getPlanePoint();
	((float*)value)[0] = planePoint[0];
	((float*)value)[1] = planePoint[1];
	((float*)value)[2] = planePoint[2];
}

void TW_CALL PartioViewer_GUI_TweakBar::setParticleRadius(const void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	const Real val = *(const Real *)(value);
	viewer->setParticleRadius(val);
	viewer->updateBoundingBox();
}

void TW_CALL PartioViewer_GUI_TweakBar::getParticleRadius(void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	*(Real *)(value) = viewer->getParticleRadius();
}

void TW_CALL PartioViewer_GUI_TweakBar::setStartFrame(const void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	const int val = *(const int *)(value);
	viewer->setStartFrame(val);
}

void TW_CALL PartioViewer_GUI_TweakBar::getStartFrame(void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	*(int *)(value) = viewer->getStartFrame();
}

void TW_CALL PartioViewer_GUI_TweakBar::setEndFrame(const void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	const int val = *(const int *)(value);
	viewer->setEndFrame(val);
}

void TW_CALL PartioViewer_GUI_TweakBar::getEndFrame(void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	*(int *)(value) = viewer->getEndFrame();
}

void TW_CALL PartioViewer_GUI_TweakBar::setFPS(const void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	const int val = *(const int *)(value);
	viewer->setFPS(val);
}

void TW_CALL PartioViewer_GUI_TweakBar::getFPS(void *value, void *clientData)
{
	PartioViewer *viewer = (PartioViewer*)clientData;
	*(int *)(value) = viewer->getFPS();
}

void TW_CALL PartioViewer_GUI_TweakBar::setCurrentFluidModel(const void* value, void* clientData)
{
	PartioViewer_GUI_TweakBar* pv = (PartioViewer_GUI_TweakBar*)clientData;
	const int val = *(const int*)(value);
	pv->m_currentFluidModel = val;
	pv->initParameterGUI();
}

void TW_CALL PartioViewer_GUI_TweakBar::getCurrentFluidModel(void* value, void* clientData)
{
	PartioViewer_GUI_TweakBar* pv = (PartioViewer_GUI_TweakBar*)clientData;
	*(int*)(value) = pv->m_currentFluidModel;
}

unsigned int PartioViewer_GUI_TweakBar::getWidth() const
{
	return MiniGL::getWidth();
}

unsigned int PartioViewer_GUI_TweakBar::getHeight() const
{
	return MiniGL::getHeight();
}
