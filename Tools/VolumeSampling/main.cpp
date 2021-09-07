#include "SPlisHSPlasH/Common.h"
#include <Eigen/Dense>
#include "extern/cxxopts/cxxopts.hpp"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/TriangleMesh.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/SystemInfo.h"
#include "Utilities/Version.h"
#include "Discregrid/All"
#include "Utilities/FileSystem.h"
#include "SPHVolumeSampling.h"
#include "SPHVolumeSampling_Jiang2015.h"

using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;

INIT_TIMING
INIT_COUNTING
INIT_LOGGING

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
#define new DEBUG_NEW 
#endif
#endif


void sampleObject(TriangleMesh &mesh);
void generateSDF(SPH::TriangleMesh &mesh);
void computeBoundingBox(TriangleMesh &mesh);

double distance(const Vector3r &x, const Real tolerance);

string exePath;
string inputFile = "";
string outputFile = "";
int output_format = 0; // 0: partio, 1: vtk
Real radius = static_cast<Real>(0.025);
Real viscosity = 0.25;
Real cohesion = 0.0;
Real adhesion = 0.0;
Real cflFactor = 0.25;
Real stiffness = 10000.0;
Real dt = 0.0005;
Vector3r scale(1.0, 1.0, 1.0);
Real diameter = radius* static_cast<Real>(2.0);
int mode = 4;
unsigned int steps = 100;
Eigen::Matrix<unsigned int, 3, 1> resolutionSDF(30, 30, 30);
SamplingBase::Region region;
bool invert = false;
bool useRegion = false;
std::vector<Vector3r> particles;
Vector3r bbmin, bbmax;
std::string outputPath = "";
bool useCache = true;
std::shared_ptr<Discregrid::CubicLagrangeDiscreteGrid> distanceField;


std::istream& operator >> (std::istream& istream, SamplingBase::Region& r)
{
	return istream >> std::skipws >> r.m_min[0] >> r.m_min[1] >> r.m_min[2] >> r.m_max[0] >> r.m_max[1] >> r.m_max[2];
}

std::ostream& operator << (std::ostream& out, const SamplingBase::Region& r)
{
	out << r.m_min[0] << ", " << r.m_min[1] << ", " << r.m_min[2] << ", " << r.m_max[0] << ", " << r.m_max[1] << ", " << r.m_max[2];
	return out;
}

std::istream& operator >> (std::istream& istream, Eigen::Matrix<unsigned int, 3, 1>& r)
{
	return istream >> std::skipws >> r[0] >> r[1] >> r[2];
}

std::ostream& operator << (std::ostream& out, const Eigen::Matrix<unsigned int, 3, 1>& r)
{
	out << r[0] << ", " << r[1] << ", " << r[2];
	return out;
}

std::istream& operator >> (std::istream& istream, Vector3r& r)
{
	return istream >> std::skipws >> r[0] >> r[1] >> r[2];
}

std::ostream& operator << (std::ostream& out, const Vector3r& r)
{
	out << r[0] << ", " << r[1] << ", " << r[2];
	return out;
}

// main 
int main(int argc, char **argv)
{
	REPORT_MEMORY_LEAKS;

	exePath = FileSystem::getProgramPath();

	Timing::m_dontPrintTimes = true;

	try
	{
		cxxopts::Options options(argv[0], "VolumeSampling - Sample a volumetric geometry given by an OBJ file.");

		options.add_options()
			("h,help", "Print help")
			("i,input", "Input file (obj)", cxxopts::value<std::string>())
			("o,output", "Output file (bgeo or vtk)", cxxopts::value<std::string>())
			("r,radius", "Particle radius", cxxopts::value<Real>()->default_value("0.025"))
			("s,scale", "Scaling of input geometry (e.g. --scale \"2 1 2\")", cxxopts::value<Vector3r>()->default_value("1 1 1"))
			("m,mode", "Mode (regular=0, almost dense=1, dense=2, Jiang et al. 2015=3, Kugelstadt et al. 2021=4)", cxxopts::value<int>()->default_value("4"))
			("region", "Region to fill with particles (e.g. --region \"0 0 0 1 1 1\")", cxxopts::value<SamplingBase::Region>())
			("steps", "SPH time steps", cxxopts::value<unsigned int>()->default_value("100"))
			("cflFactor", "CFL factor", cxxopts::value<Real>()->default_value("0.25"))
			("viscosity", "Viscosity coefficient (XSPH)", cxxopts::value<Real>()->default_value("0.25"))
			("cohesion", "Cohesion coefficient", cxxopts::value<Real>())
			("adhesion", "Adhesion coefficient", cxxopts::value<Real>())
			("stiffness", "Stiffness coefficient (only mode 3)", cxxopts::value<Real>()->default_value("10000.0"))
			("dt", "Time step size (only mode 3)", cxxopts::value<Real>()->default_value("0.0005"))
			("res", "Resolution of the Signed Distance Field (e.g. --res \"30 30 30\")", cxxopts::value<Eigen::Matrix<unsigned int, 3, 1>>())
			("invert", "Invert the SDF to sample the outside of the object in the bounding box/region")
			("no-cache", "Disable caching of SDF.")
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
			outputFile = result["output"].as<std::string>();
			output_format = 0;
			if (Utilities::StringTools::to_upper(FileSystem::getFileExt(outputFile)) == "VTK")
				output_format = 1;
		}
		else
		{
			std::cout << "Input or output missing!" << std::endl;
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(1);
		}

		if (result.count("mode"))
			mode = result["mode"].as<int>();

		outputPath = FileSystem::normalizePath(exePath + "/output/" + FileSystem::getFileName(inputFile));

		Utilities::logger.addSink(unique_ptr<Utilities::ConsoleSink>(new Utilities::ConsoleSink(Utilities::LogLevel::INFO)));
		std::string logPath = FileSystem::normalizePath(outputPath + "/log");
		FileSystem::makeDirs(logPath);
		Utilities::logger.addSink(unique_ptr<Utilities::FileSink>(new Utilities::FileSink(Utilities::LogLevel::DEBUG, logPath + "/SPH_log.txt")));

		LOG_INFO << "SPlisHSPlasH version: " << SPLISHSPLASH_VERSION;
		LOG_DEBUG << "Git refspec:          " << GIT_REFSPEC;
		LOG_DEBUG << "Git SHA1:             " << GIT_SHA1;
		LOG_DEBUG << "Git status:           " << GIT_LOCAL_STATUS;
		LOG_DEBUG << "Host name:            " << Utilities::SystemInfo::getHostName();

		LOG_INFO << "Input = " << inputFile;
		LOG_INFO << "Output = " << outputFile;
		LOG_INFO << "Mode: " << mode;

#ifdef DL_OUTPUT
		std::string modelsFilePath = FileSystem::normalizePath(outputPath + "/models");
		FileSystem::makeDirs(modelsFilePath);
		FileSystem::copyFile(inputFile, modelsFilePath + "/" + FileSystem::getFileNameWithExt(inputFile));
#endif 

		if (result.count("radius"))
			radius = result["radius"].as<Real>();
		LOG_INFO << "Radius: " << radius;

		if (result.count("scale"))
			scale = result["scale"].as<Vector3r>();
		LOG_INFO << "Scale: " << scale;

		if (result.count("steps"))
			steps = result["steps"].as<unsigned int>();
		LOG_INFO << "SPH steps: " << steps;

		if (result.count("cflFactor"))
			cflFactor = result["cflFactor"].as<Real>();
		LOG_INFO << "CFL Factor: " << cflFactor;

		if (result.count("viscosity"))
			viscosity = result["viscosity"].as<Real>();
		LOG_INFO << "Viscosity coefficient: " << viscosity;

		if (result.count("cohesion"))
			cohesion = result["cohesion"].as<Real>();

		if (result.count("adhesion"))
			adhesion = result["adhesion"].as<Real>();

		if (result.count("stiffness"))
			stiffness = result["stiffness"].as<Real>();
		LOG_INFO << "Stiffness coefficient: " << stiffness;

		if (result.count("dt"))
			dt = result["dt"].as<Real>();
		LOG_INFO << "Time step size: " << dt;

		if (result.count("res"))
		{
			resolutionSDF = result["res"].as<Eigen::Matrix<unsigned int, 3, 1>>();
			LOG_INFO << "SDF resolution: " << resolutionSDF;
		}

		if (result.count("invert"))
		{
			invert = true;
		}

		if (result.count("no-cache"))
			useCache = false;

		if (result.count("region"))
		{
			region = result["region"].as<SamplingBase::Region>();
			useRegion = true;
			LOG_INFO << "Region: " << region;
		}
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}

	SPHSamplingBase* sampling = nullptr;
	if (mode == 3)
	{
		sampling = new SPHVolumeSampling_Jiang2015();
		((SPHVolumeSampling_Jiang2015*)sampling)->setStiffness(stiffness);
		((SPHVolumeSampling_Jiang2015*)sampling)->setTimeStepSize(dt);
	}
	else if (mode == 4)
		sampling = new SPHVolumeSampling();
	else if (mode > 4)
		LOG_ERR << "Mode not supported!";

	if (sampling != nullptr)
	{
		sampling->setOutputPath(outputPath);
		sampling->setScale(scale);
		sampling->setUseRegion(useRegion);
		sampling->setRegion(region);
		sampling->setResolutionSdf(resolutionSDF);
		sampling->setViscosity(viscosity);
		sampling->setCFLFactor(cflFactor);
		sampling->setSteps(steps);
		sampling->setRadius(radius);
		sampling->setInvert(invert);
		sampling->setUseCache(useCache);
		if (cohesion != 0.0)
			sampling->setCohesion(cohesion);
		if (adhesion != 0.0)
			sampling->setAdhesion(adhesion);

		LOG_INFO << "Cohesion coefficient: " << sampling->getCohesion();
		LOG_INFO << "Adhesion coefficient: " << sampling->getAdhesion();
			
		sampling->generateSampling(inputFile, outputFile);
		Timing::printTimeSums();
		Timing::printAverageTimes();

		return 0;
	}

	diameter = static_cast<Real>(2.0) * radius;
	TriangleMesh mesh;
	SamplingBase::loadObj(inputFile, mesh, scale);
	computeBoundingBox(mesh);
	START_TIMING("generateSDF");
	generateSDF(mesh);
	STOP_TIMING_AVG_PRINT;

	sampleObject(mesh);

	if (output_format == 0)
		PartioReaderWriter::writeParticles(outputFile, (unsigned int)particles.size(), particles.data(), NULL, 0.0);
	else
		SamplingBase::writeParticlesVTK(outputFile, particles);

	LOG_INFO << "Generated particles: " << particles.size();

	Timing::printTimeSums();
	Timing::printAverageTimes();

	return 0;
}

 void sampleObject(TriangleMesh &mesh)
 {
	 int surface = 0;

 	// sample object
 	const unsigned int numberOfSamplePoints = (((unsigned int)((1.0f / diameter) * (bbmax[2] - bbmin[2]))) + 1) *
		(((unsigned int)((1.0f / diameter) * (bbmax[1] - bbmin[1]))) + 1) *
		(((unsigned int)((1.0f / diameter) * (bbmax[0] - bbmin[0]))) + 1);
 	unsigned int currentSample = 0;
 	Real currentPercent = 0.01;
 	int counter_x = 0;
 	int counter_y = 0;
 	Real xshift = diameter;
 	Real yshift = diameter;
 
 	if ((mode == 1) || (mode == 4))
 		yshift = sqrt(static_cast<Real>(3.0)) * radius;
 	else if (mode == 2)
 	{
 		xshift = sqrt(static_cast<Real>(3.0)) * radius;
 		yshift = sqrt(static_cast<Real>(6.0)) * diameter / static_cast<Real>(3.0);
 	}
 	for (Real z = bbmin[2]; z <= bbmax[2]; z += diameter)
 	{
 		for (Real y = bbmin[1]; y <= bbmax[1]; y += yshift)
 		{
 			for (Real x = bbmin[0]; x <= bbmax[0]; x += xshift)
 			{
 				Vector3r particlePosition;
				if ((mode == 1) || (mode == 4))
 				{					
 					if (counter_y % 2 == 0)
 						particlePosition = Vector3r(x, y + radius, z + radius);
 					else
 						particlePosition = Vector3r(x + radius, y + radius, z);
 				}
 				else if (mode == 2)
 				{
 					particlePosition = Vector3r(x, y + radius, z + radius);
 
 					Vector3r shift_vec(0, 0, 0);
 					if (counter_x % 2)
 					{
 						shift_vec[2] += diameter / (static_cast<Real>(2.0) * (counter_y % 2 ? -1 : 1));
 					}
 					if (counter_y % 2)
 					{
 						shift_vec[0] += xshift / static_cast<Real>(2.0);
 						shift_vec[2] += diameter / static_cast<Real>(2.0);
 					}
 					particlePosition += shift_vec;
 				}
 				else
 				{
 					// Use center of voxel
 					particlePosition = Vector3r(x + radius, y + radius, z + radius);
 				}
 
				const double dist = distance(particlePosition, 0.0);
				if ((dist != std::numeric_limits<double>::max()) && (dist > 0.0))
				{
					particles.push_back(particlePosition);
				}
 				currentSample++;
 
 				if ((Real)currentSample / (Real)numberOfSamplePoints > currentPercent)
 				{
					std::cout << "\r" << "Generating samples " << std::setw(14)	<< currentPercent * 100.0 << "%";
 					currentPercent += 0.01;
 				}
 				counter_x++;
 			}
 			counter_x = 0;
 			counter_y++;
 		}
 		counter_y = 0;
 	}
	std::cout << "\r" << "Generating samples " << std::setw(14) << 100.0 << "%\n";
	LOG_INFO << "Surface particles: " << surface << "\n";
}

void computeBoundingBox(TriangleMesh &mesh)
{
	Vector3r *v = mesh.getVertices().data();

	// compute bounding box	 
	bbmin = v[0];
	bbmax = bbmin;
	for (unsigned int i = 1; i < mesh.numVertices(); ++i)
	{
		const Vector3r& p = v[i];
		for (unsigned int j = 0; j < 3; ++j) {
			bbmin[j] = std::min(bbmin[j], p[j]);
			bbmax[j] = std::max(bbmax[j], p[j]);
		}
	}

	if (useRegion)
	{
		for (unsigned int i = 0; i < 3; i++)
		{
			bbmin[i] = std::max(scale[i] * region.m_min[i], bbmin[i]);
			bbmax[i] = std::min(scale[i] * region.m_max[i], bbmax[i]);
		}
	}
}
 
// Determine distance of a point x to the surface of the mesh 
double distance(const Vector3r &x, const Real tolerance)
{
	const double dist = distanceField->interpolate(0, x.template cast<double>());
	if (dist == std::numeric_limits<double>::max())
		return dist;
	return dist - tolerance;
}

void generateSDF(SPH::TriangleMesh &mesh)
{
	//////////////////////////////////////////////////////////////////////////
	// Generate distance field of object using Discregrid
	//////////////////////////////////////////////////////////////////////////
#ifdef USE_DOUBLE
	Discregrid::TriangleMesh sdfMesh(&mesh.getVertices()[0][0], mesh.getFaces().data(), mesh.numVertices(), mesh.numFaces());
#else
	// if type is float, copy vector to double vector
	std::vector<double> doubleVec;
	doubleVec.resize(3 * mesh.numVertices());
	for (unsigned int i = 0; i < mesh.numVertices(); i++)
		for (unsigned int j = 0; j < 3; j++)
			doubleVec[3 * i + j] = mesh.getVertices()[i][j];
	Discregrid::TriangleMesh sdfMesh(&doubleVec[0], mesh.getFaces().data(), mesh.numVertices(), mesh.numFaces());
#endif

	Discregrid::MeshDistance md(sdfMesh);
	Eigen::AlignedBox3d domain;
	domain.extend(bbmin.cast<double>());
	domain.extend(bbmax.cast<double>());
	//domain.max() += 1.0e-3 * domain.diagonal().norm() * Eigen::Vector3d::Ones();
	//domain.min() -= 1.0e-3 * domain.diagonal().norm() * Eigen::Vector3d::Ones();
	domain.max() += Eigen::Vector3d::Ones() * radius * 4.0;
	domain.min() -= Eigen::Vector3d::Ones() * radius * 4.0;

	LOG_INFO << "Set SDF resolution: " << resolutionSDF[0] << ", " << resolutionSDF[1] << ", " << resolutionSDF[2];
	distanceField = std::make_shared<Discregrid::CubicLagrangeDiscreteGrid>(domain, std::array<unsigned int, 3>({ resolutionSDF[0], resolutionSDF[1], resolutionSDF[2] }));
	auto func = Discregrid::DiscreteGrid::ContinuousFunction{};

	// invert the distance field since the particles should stay inside
	func = [&md](Eigen::Vector3d const& xi) {return -md.signedDistanceCached(xi); };

	LOG_INFO << "Generate SDF";
	distanceField->addFunction(func, false);
}


