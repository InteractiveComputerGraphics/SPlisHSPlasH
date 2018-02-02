#include "SPlisHSPlasH/Common.h"
#include <Eigen/Dense>
#include <iostream>
#include "Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/Utilities/PoissonDiskSampling.h"
#include "Utilities/FileSystem.h"
#include "Utilities/StringTools.h"
#include "Utilities/Version.h"
#include "SPlisHSPlasH/TriangleMesh.h"
#include "extern/cxxopts/cxxopts.hpp"

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
	#define new DEBUG_NEW 
#endif
#endif

INIT_TIMING
INIT_LOGGING

using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;

void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale);

string inputFile = "";
string outputFile = "";
Real particleRadius = 0.025;
Vector3r scale = Vector3r::Ones();

std::istream& operator >> (std::istream& istream, Vector3r& v)
{
	return istream >> std::skipws >> v[0] >> v[1] >> v[2];
}


// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;

	std::cout << "Git refspec: " << GIT_REFSPEC << std::endl;
	std::cout << "Git SHA1: " << GIT_SHA1 << std::endl;
	std::cout << "Git status: " << GIT_LOCAL_STATUS << std::endl;

	try
	{
		cxxopts::Options options(argv[0], "SurfaceSampling - Sample a surface geometry given by an OBJ file.");

		options.add_options()
			("h,help", "Print help")
			("i,input", "Input file (obj)", cxxopts::value<std::string>())
			("o,output", "Output file (bgeo)", cxxopts::value<std::string>())
			("r,radius", "Particle radius", cxxopts::value<Real>()->default_value("0.025"))
			("s,scale", "Scaling of input geometry (e.g. --scale \"1 2 3\")", cxxopts::value<Vector3r>())
			;

		options.parse(argc, argv);

		if (options.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
		}

		if (options.count("input") && options.count("output"))
		{
			inputFile = options["input"].as<std::string>();
			std::cout << "Input = " << inputFile << std::endl;
			outputFile = options["output"].as<std::string>();
			std::cout << "Output = " << outputFile << std::endl;
		}
		else
		{
			std::cout << "Input or output missing!" << std::endl;
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(1);
		}

		if (options.count("radius"))
			particleRadius = options["radius"].as<Real>();
		cout << "Radius: " << particleRadius << endl;

		if (options.count("scale"))
			scale = options["scale"].as<Vector3r>();
		cout << "Scale: " << scale << endl;

		options.parse(argc, argv);
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}


	TriangleMesh mesh;
	loadObj(inputFile, mesh, scale);

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

void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale)
{
	std::vector<OBJLoader::Vec3f> x;
	std::vector<OBJLoader::Vec3f> normals;
	std::vector<MeshFaceIndices> faces;
	OBJLoader::Vec3f s = { (float)scale[0], (float)scale[1], (float)scale[2] };
	OBJLoader::loadObj(filename, &x, &faces, &normals, nullptr, s);

	mesh.release();
	const unsigned int nPoints = (unsigned int)x.size();
	const unsigned int nFaces = (unsigned int)faces.size();
	mesh.initMesh(nPoints, nFaces);
	for (unsigned int i = 0; i < nPoints; i++)
	{
		mesh.addVertex(Vector3r(x[i][0], x[i][1], x[i][2]));
	}
	for (unsigned int i = 0; i < nFaces; i++)
	{
		// Reduce the indices by one
		int posIndices[3];
		for (int j = 0; j < 3; j++)
		{
			posIndices[j] = faces[i].posIndices[j] - 1;
		}

		mesh.addFace(&posIndices[0]);
	}

	LOG_INFO << "Number of triangles: " << nFaces;
	LOG_INFO << "Number of vertices: " << nPoints;
}
