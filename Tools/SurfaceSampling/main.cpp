#include "SPlisHSPlasH/Common.h"
#include <Eigen/Dense>
#include <iostream>
#include "Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/Utilities/SurfaceSampling.h"
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
using namespace Utilities;

void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale);

std::string inputFile = "";
std::string outputFile = "";
Real particleRadius = 0.025;
Vector3r scale = Vector3r::Ones();
unsigned int samplingMode = 1;
Vector3r translation = Vector3r::Zero();
Vector3r rotationAxis = Vector3r::UnitY();
Real angle = 0;

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
			("m,mode", "Sampling mode 0 Poisson disk, 1 Regular, 2 2D sampling", cxxopts::value<unsigned int>()->default_value("1"))
			("t,translation", "Translation for 2D sampling (default: \"0 0 0\")", cxxopts::value<Vector3r>())
			("rotationAxis", "Rotation axis for 2D sampling (default: \"0 1 0\")", cxxopts::value<Vector3r>())
			("a,angle", "Rotation angle for 2D simulation", cxxopts::value<Real>()->default_value("0"))
			;

		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
		}

		if (result.count("input") && result.count("output"))
		{
			inputFile = result["input"].as<std::string>();
			std::cout << "Input = " << inputFile << std::endl;
			outputFile = result["output"].as<std::string>();
			std::cout << "Output = " << outputFile << std::endl;
		}
		else
		{
			std::cout << "Input or output missing!" << std::endl;
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(1);
		}

		if (result.count("radius"))
			particleRadius = result["radius"].as<Real>();
		std::cout << "Radius: " << particleRadius << std::endl;

		if (result.count("scale"))
			scale = result["scale"].as<Vector3r>();
		std::cout << "Scale: [" << scale.transpose() << "]^T" << std::endl;

		if (result.count("mode"))
			samplingMode = result["mode"].as<unsigned int>();
		std::cout << "Sampling mode: " << samplingMode << std::endl;

		if (result.count("translation"))
			translation = result["translation"].as<Vector3r>();
		if (result.count("rotationAxis"))
			rotationAxis = result["rotationAxis"].as<Vector3r>();
		if (result.count("angle"))
			angle = result["angle"].as<Real>();

		if (samplingMode == 2)
		{
			std::cout << "Translation:    [" << translation.transpose() << "]^T" << std::endl;
			std::cout << "Rotation axis:  [" << rotationAxis.transpose() << "]^T" << std::endl;
			std::cout << "Rotation angle: " << angle << std::endl;
		}
		else if (result.count("translation") || result.count("rotationAxis") || result.count("angle"))
		{
			std::cout << "--translation, --rotationAxis and --angle only affect 2D sampling and are ignored." << std::endl;
		}
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}


	TriangleMesh mesh;
	loadObj(inputFile, mesh, scale);

	std::vector<Vector3r> samplePoints;

	const auto samplePoissonDisk = [&]()
	{
		std::cout << "Poisson disk surface sampling of " << inputFile << std::endl;
		START_TIMING("Poisson disk sampling");
		PoissonDiskSampling sampling;
		sampling.sampleMesh(mesh.numVertices(), mesh.getVertices().data(), mesh.numFaces(), mesh.getFaces().data(), particleRadius, 10, 1, samplePoints);
		STOP_TIMING_PRINT;
	};
	const auto sampleRegularTriangle = [&]()
	{
		std::cout << "Regular triangle surface sampling of " << inputFile << std::endl;
		START_TIMING("Regular triangle sampling");
		RegularTriangleSampling sampling;
		sampling.sampleMesh(mesh.numVertices(), mesh.getVertices().data(), mesh.numFaces(), mesh.getFaces().data(), 1.5f * particleRadius, samplePoints);
		STOP_TIMING_PRINT;
	};
	const auto sampleRegular2D = [&]()
	{
		std::cout << "2D regular sampling of " << inputFile << std::endl;
		START_TIMING("2D regular sampling");
		RegularSampling2D sampling;
		sampling.sampleMesh(AngleAxisr(angle, rotationAxis).toRotationMatrix(), translation,
			mesh.numVertices(), mesh.getVertices().data(), mesh.numFaces(),
			mesh.getFaces().data(), 1.75f * particleRadius, samplePoints);
		STOP_TIMING_AVG;
	};
	if (SurfaceSamplingMode::PoissonDisk == samplingMode)
		samplePoissonDisk();
	else if (SurfaceSamplingMode::RegularTriangle == samplingMode)
		sampleRegularTriangle();
	else if (SurfaceSamplingMode::Regular2D == samplingMode)
		sampleRegular2D();
	else
	{
		std::cout << "Unknown surface sampling method: " << samplingMode;
		std::cout << "Falling back to:";
		sampleRegularTriangle();
	}
	
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
