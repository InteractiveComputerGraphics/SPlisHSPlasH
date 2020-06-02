#include "SPlisHSPlasH/Common.h"
#include "Simulator/SimulatorBase.h"
#include "Simulator/GUI/OpenGL/Simulator_OpenGL.h"
#include "PositionBasedDynamicsWrapper/PBDBoundarySimulator.h"
#ifdef USE_IMGUI
	#include "Simulator/GUI/imgui/Simulator_GUI_imgui.h"
	#include "PositionBasedDynamicsWrapper/PBD_Simulator_GUI_imgui.h"
#else
	#include "Simulator/GUI/TweakBar/Simulator_GUI_TweakBar.h"
	#include "PositionBasedDynamicsWrapper/PBD_Simulator_GUI_TweakBar.h"
#endif

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
	#define new DEBUG_NEW 
#endif
#endif

using namespace SPH;
using namespace std;

SimulatorBase *base = nullptr;
Simulator_GUI_Base *gui = nullptr;

// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;
	base = new SimulatorBase();
	base->init(argc, argv, "SPlisHSPlasH");

	if (base->getUseGUI())
	{
#ifdef USE_IMGUI
		if (base->isStaticScene())
			gui = new Simulator_GUI_imgui(base);
		else
			gui = new PBD_Simulator_GUI_imgui(base, ((PBDBoundarySimulator*)base->getBoundarySimulator())->getPBDWrapper());
#else
		if (base->isStaticScene())
			gui = new Simulator_GUI_TweakBar(base);
		else
			gui = new PBD_Simulator_GUI_TweakBar(base, ((PBDBoundarySimulator*)base->getBoundarySimulator())->getPBDWrapper());
#endif
		base->setGui(gui);
	}
	base->run();

	delete base;
	delete gui;
	
	return 0;
}

