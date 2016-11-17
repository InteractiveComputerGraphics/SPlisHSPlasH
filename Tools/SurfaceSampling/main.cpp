#include "SPlisHSPlasH/Common.h"
#include <Eigen/Dense>
#include <iostream>
#include "SPlisHSPlasH/Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/Utilities/PoissonDiskSampling.h"
#include "Utilities/FileSystem.h"
#include "Utilities/StringTools.h"

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
	#define new DEBUG_NEW 
#endif
#endif

using namespace SPH;
using namespace Eigen;
using namespace std;

string inputFile = "";
string outputFile = "";
Real particleRadius = 0.025;
Vector3r scale = Vector3r::Ones();

// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;

	for (int i = 1; i < argc; i++)
	{
		string argStr = argv[i];
		string type_str = argStr.substr(0, 2);
		if ((type_str == "-r") && (i+1 < argc))
			particleRadius = stof(argv[++i]);
		else if ((type_str == "-s") && (i + 1 < argc))
		{
			vector<string> tokens;
			StringTools::tokenize(argv[++i], tokens, ",");
			scale[0] = stof(tokens[0]);
			scale[1] = stof(tokens[1]);
			scale[2] = stof(tokens[2]);
		}
		else if (i + 1 < argc)
		{
			inputFile = argv[i];
			outputFile = argv[++i];
		}
		else
		{
			std::cerr << "Not enough parameters!\n";
			std::cerr << "Usage: SurfaceSampling.exe [-r particle_radius] [-s scaleX,scaleY,scaleZ] in.obj out.bgeo\n";
			return -1;
		}
	}

	TriangleMesh mesh;
	OBJLoader::loadObj(inputFile, mesh, scale);

	std::cout << "Surface sampling of " << inputFile << "\n";
	START_TIMING("Poisson disk sampling");
	PoissonDiskSampling sampling;
	std::vector<Vector3r> samplePoints;
	sampling.sampleMesh(mesh.numVertices(), mesh.getVertices().data(), mesh.numFaces(), mesh.getFaces().data(), particleRadius, 10, 1, samplePoints);
	STOP_TIMING_AVG;
	std::cout << "Number of sample points: " << samplePoints.size() << "\n";


	PartioReaderWriter::writeParticles(outputFile, (unsigned int) samplePoints.size(), samplePoints.data(), NULL, particleRadius);

	Timing::printAverageTimes();
	Timing::printTimeSums();
	
	return 0;
}

