#include "SPlisHSPlasH/Common.h"
#include "Utilities/Logger.h"
#include "JsonSchemaGenerator.h"
#include "extern/cxxopts/cxxopts.hpp"
#include "SceneExampleGenerator.h"
#include "Utilities/FileSystem.h"

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
	#define new DEBUG_NEW 
#endif
#endif

using namespace SPH;
using namespace Eigen;
using namespace Utilities;
using namespace std;
using namespace GenParam;

string outputFile = "";
int mode = 0;

// main 
int main(int argc, char** argv)
{
	REPORT_MEMORY_LEAKS;

	try
	{
		cxxopts::Options options(argv[0], "ParameterParser - Parses all generic parameters of the simulation framework and generates an example scene with all parameters or a json schema.");

		options.add_options()
			("h,help", "Print help")
			("o,output", "Output file", cxxopts::value<std::string>())
			("m,mode", "Mode:\n  - 0: generate example scene with all parameters\n  - 1: generate json schema\n", cxxopts::value<int>()->default_value("0"))
			;

		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
		}

		if (result.count("mode"))
			mode = result["mode"].as<int>();

		// set default file name if none is defined
		std::string exePath = FileSystem::getProgramPath();
		outputFile = FileSystem::normalizePath(exePath + "/all-parameters.json");
		if (mode == 1)
			outputFile = FileSystem::normalizePath(exePath + "/scene-file-schema.json");

		if (result.count("output"))
			outputFile = result["output"].as<std::string>();
	}
	catch (const cxxopts::exceptions::exception& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}

	if (mode == 0)
	{
		SceneExampleGenerator sceneGen;
		sceneGen.generateExampleSceneFile(outputFile);
	}
	else
	{
		JsonSchemaGenerator schemaGen;
		schemaGen.generateSchemaFile(outputFile);
	}
	return 0;
}