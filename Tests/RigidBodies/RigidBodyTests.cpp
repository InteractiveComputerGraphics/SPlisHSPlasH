#include "SPlisHSPlasH/Common.h"
#include "Simulator/SimulatorBase.h"
#include "Simulator/SceneConfiguration.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Utilities/FileSystem.h"
#include "Utilities/Logger.h"
#include <fstream>

// Let Catch provide main():
#define CATCH_CONFIG_MAIN

#include "catch.hpp"


// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
#define new DEBUG_NEW
#endif
#endif

using namespace SPH;
using namespace std;

/** Write a scene file with dynamic boxes, an animated box and, in the end, a static tank.
*/
void writeScene(const std::string &sceneFile, const int boundaryHandlingMethod, const unsigned int numDynamicBodies)
{
	std::ofstream file(sceneFile);
	auto writeBox = [&](const double x, const double size, const bool dynamic, const bool animated, const bool wall)
	{
		file << "\t\t{ \"geometryFile\": \"../data/models/UnitBox.obj\", \"translation\": [" << x << ", 0, 0], "
			<< "\"scale\": [" << size << ", " << size << ", " << size << "], \"isDynamic\": " << (dynamic ? "true" : "false")
			<< ", \"isAnimated\": " << (animated ? "true" : "false") << ", \"isWall\": " << (wall ? "true" : "false")
			<< ", \"mapInvert\": " << (wall ? "true" : "false") << ", \"mapResolution\": [5, 5, 5] }" << (wall ? "" : ",") << "\n";
	};

	file << "{\n\t\"Configuration\": { \"particleRadius\": 0.05, \"simulationMethod\": 4, \"boundaryHandlingMethod\": " << boundaryHandlingMethod << " },\n";
	file << "\t\"RigidBodies\": [\n";
	writeBox(1.0, 0.25, false, true, false);
	for (unsigned int i = 0; i < numDynamicBodies; i++)
		writeBox(2.0 + i, 0.25, true, false, false);
	writeBox(-2.0, 1.0, false, false, true);
	file << "\t],\n\t\"FluidBlocks\": [ { \"start\": [-2.2, -0.4, -0.2], \"end\": [-1.8, 0, 0.2] } ]\n}\n";
}

TEST_CASE("Boundary models of a dynamic scene", "[rigid_bodies]")
{
	REPORT_MEMORY_LEAKS;

	Utilities::logger.activate(false);

	const std::string exePath = Utilities::FileSystem::getProgramPath();
	const std::string sceneFile = Utilities::FileSystem::normalizePath(exePath + "/RigidBodyTests.json");

	for (int method = 0; method < static_cast<int>(BoundaryHandlingMethods::NumSimulationMethods); method++)
	{
		for (const unsigned int numDynamicBodies : { 2u, 20u })
		{
			writeScene(sceneFile, method, numDynamicBodies);

			SimulatorBase* base = new SimulatorBase();
			base->init({ exePath, sceneFile }, "SPlisHSPlasH");
			base->setUseGUI(false);
			base->setUseParticleCaching(false);
			base->initSimulation();
			base->deferredInit();

			// each boundary model is initialized with the parameters of its own rigid body
			const Utilities::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();
			Simulation* sim = Simulation::getCurrent();
			REQUIRE(sim->numberOfBoundaryModels() == numDynamicBodies + 2);
			REQUIRE(scene.boundaryModels.size() == numDynamicBodies + 2);
			for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
			{
				RigidBodyObject* rb = sim->getBoundaryModel(i)->getRigidBodyObject();
				REQUIRE(rb->isDynamic() == scene.boundaryModels[i]->dynamic);
				REQUIRE((rb->getPosition() - scene.boundaryModels[i]->translation).norm() < 1.0e-5);
			}

			base->cleanup();
			delete base;
		}
	}
}
