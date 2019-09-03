#include "SPlisHSPlasH/Common.h"
#include <Eigen/Dense>
#include "extern/cxxopts/cxxopts.hpp"
#include "Utilities/Timing.h"
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/Utilities/VolumeSampling.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/Version.h"
#include "SPlisHSPlasH/TriangleMesh.h"
#include "Discregrid/All"

using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;

INIT_TIMING
INIT_LOGGING

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
#define new DEBUG_NEW 
#endif
#endif

void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale);

string exePath;
string inputFile = "";
string outputFile = "";
Real radius = 0.025f;
Real scale = 1.0;
Real diameter = radius*2.0;
int mode = 0;
AlignedBox3r region;
bool useRegion = false;
std::vector<Vector3r> particles;
std::array<unsigned int, 3> resolutionSDF = { 30, 30, 30 };
Vector3r bbmin, bbmax;
std::shared_ptr<Discregrid::CubicLagrangeDiscreteGrid> distanceField;
bool invert = false;

std::istream& operator >> (std::istream& istream, AlignedBox3r &r)
{
	return istream >> std::skipws >> r.min()[0] >> r.min()[1] >> r.min()[2] >> r.max()[0] >> r.max()[1] >> r.max()[2];
}

std::ostream& operator << (std::ostream& out, const AlignedBox3r& r)
{
	out << r.min()[0] << ", " << r.min()[1] << ", " << r.min()[2] << ", " << r.max()[0] << ", " << r.max()[1] << ", " << r.max()[2];
	return out;
}

std::istream& operator >> (std::istream& istream, std::array<unsigned int, 3>& r)
{
	return istream >> std::skipws >> r[0] >> r[1] >> r[2];
}

std::ostream& operator << (std::ostream& out, const std::array<unsigned int, 3>& r)
{
	out << r[0] << ", " << r[1] << ", " << r[2];
	return out;
}

// main 
int main(int argc, char **argv)
{
	REPORT_MEMORY_LEAKS;

	Utilities::logger.addSink(unique_ptr<Utilities::ConsoleSink>(new Utilities::ConsoleSink(Utilities::LogLevel::INFO)));

	LOG_INFO << "Git refspec: " << GIT_REFSPEC;
	LOG_INFO << "Git SHA1: " << GIT_SHA1;
	LOG_INFO << "Git status: " << GIT_LOCAL_STATUS;

	Timing::m_dontPrintTimes = true;

	exePath = string(argv[0]);
	const size_t npos = exePath.rfind("\\");
	if (npos != string::npos)
		exePath = exePath.substr(0, exePath.rfind("\\") + 1);
	else
		exePath = ".";

	try
	{
		cxxopts::Options options(argv[0], "VolumeSampling - Sample a volumetric geometry given by an OBJ file.");

		options.add_options()
			("h,help", "Print help")
			("i,input", "Input file (obj)", cxxopts::value<std::string>())
			("o,output", "Output file (bgeo)", cxxopts::value<std::string>())
			("r,radius", "Particle radius", cxxopts::value<Real>()->default_value("0.025"))
			("s,scale", "Scaling of input geometry", cxxopts::value<Real>()->default_value("1.0"))
			("m,mode", "Mode (regular=0, almost dense=1, dense=2)", cxxopts::value<int>()->default_value("0"))
			("region", "Region to fill with particles (e.g. --region \"0 0 0 1 1 1\"", cxxopts::value<AlignedBox3r>())
			("res", "Resolution of the Signed Distance Field (e.g. --res \"30 30 30\"", cxxopts::value<std::array<unsigned int, 3>>())
			("invert", "Invert the SDF to sample the outside of the object in the bounding box/region")
			;

		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			LOG_INFO << options.help({ "", "Group" });
			exit(0);
		}

		if (result.count("input") && result.count("output"))
		{
			inputFile = result["input"].as<std::string>();
			LOG_INFO << "Input = " << inputFile;
			outputFile = result["output"].as<std::string>();
			LOG_INFO << "Output = " << outputFile;
		}
		else
		{
			LOG_INFO << "Input or output missing!";
			LOG_INFO << options.help({ "", "Group" });
			exit(1);
		}

		if (result.count("radius"))
			radius = result["radius"].as<Real>();
		LOG_INFO << "Radius: " << radius;

		if (result.count("scale"))
			scale = result["scale"].as<Real>();
		LOG_INFO << "Scale: " << scale;

		if (result.count("mode"))
			mode = result["mode"].as<int>();
		LOG_INFO << "Mode: " << mode;

		if (result.count("region"))
		{
			region = result["region"].as<AlignedBox3r>();
			useRegion = true;
			LOG_INFO << "Region: " << region;
		}

		if (result.count("res"))
		{
			resolutionSDF = result["res"].as<std::array<unsigned int, 3>>();
			LOG_INFO << "SDF resolution: " << resolutionSDF;
		}

		if (result.count("invert"))
		{
			invert = true;
		}
	}
	catch (const cxxopts::OptionException& e)
	{
		LOG_INFO << "error parsing options: " << e.what();
		exit(1);
	}

	TriangleMesh mesh;
	loadObj(inputFile, mesh, scale*Vector3r::Ones());

	START_TIMING("Volume sampling");
	if (useRegion)
	{
		region.min() = scale * region.min();
		region.max() = scale * region.max();
		Utilities::VolumeSampling::sampleMesh(mesh.numVertices(), mesh.getVertices().data(), mesh.numFaces(), mesh.getFaces().data(),
			radius, &region, resolutionSDF, invert, mode, particles);
	}
	else
	{
		Utilities::VolumeSampling::sampleMesh(mesh.numVertices(), mesh.getVertices().data(), mesh.numFaces(), mesh.getFaces().data(),
			radius, nullptr, resolutionSDF, invert, mode, particles);
	}
	STOP_TIMING_PRINT;
	PartioReaderWriter::writeParticles(outputFile, (unsigned int)particles.size(), particles.data(), NULL, radius);

	LOG_INFO << "Generated particles: " << particles.size();

	Timing::printAverageTimes();

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
