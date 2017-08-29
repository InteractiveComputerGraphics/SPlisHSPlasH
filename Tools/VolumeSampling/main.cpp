#include "SPlisHSPlasH/Common.h"
#include <Eigen/Dense>
#include "extern/cxxopts/cxxopts.hpp"
#include "SPlisHSPlasH/Utilities/Timing.h"
#include "Utilities/OBJLoader.h"
#include "WindingNumbers.h"
#include "Utilities/PartioReaderWriter.h"

using namespace SPH;
using namespace Eigen;
using namespace std;


// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
#define new DEBUG_NEW 
#endif
#endif


void sampleObject(TriangleMesh &mesh);
void partioExport();

struct Region
{
public:
	typedef Real value_type;


	Region() {}
	Region(Real minx, Real miny, Real minz,
		Real maxx, Real maxy, Real maxz)
	{
		m_min = Vector3r(minx, miny, minz);
		m_max = Vector3r(maxx, maxy, maxz);
	}

	Vector3r m_min;
	Vector3r m_max;
};

string exePath;
string inputFile = "";
string outputFile = "";
Real radius = 0.05f;
Real scale = 1.0;
Real diameter = radius*2.0;
int mode = 0;
Region region;
bool useRegion = false;
std::vector<Vector3r> particles;


std::istream& operator >> (std::istream& istream, Region& r)
{
	return istream >> std::skipws >> r.m_min[0] >> r.m_min[1] >> r.m_min[2] >> r.m_max[0] >> r.m_max[1] >> r.m_max[2];
}

std::ostream& operator << (std::ostream& out, const Region& r)
{
	out << r.m_min[0] << ", " << r.m_min[1] << ", " << r.m_min[2] << ", " << r.m_max[0] << ", " << r.m_max[1] << ", " << r.m_max[2];
	return out;
}



// main 
int main(int argc, char **argv)
{
	REPORT_MEMORY_LEAKS;

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
			("region", "Region to fill with particles (e.g. --region \"0 0 0 1 1 1\"", cxxopts::value<Region>())
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
			std::cout << "Input or output missing!" << std::endl;;
			exit(1);
		}

		if (options.count("radius"))
			radius = options["radius"].as<Real>();
		cout << "Radius: " << radius << endl;

		if (options.count("scale"))
			scale = options["scale"].as<Real>();
		cout << "Scale: " << scale << endl;

		if (options.count("mode"))
			mode = options["mode"].as<int>();
		cout << "Mode: " << mode << endl;

		if (options.count("region"))
		{
			region = options["region"].as<Region>();
			useRegion = true;
			cout << "Region: " << region << endl;
		}
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}


	diameter = 2.0 * radius;
	TriangleMesh mesh;
	OBJLoader::loadObj(inputFile, mesh, scale*Vector3r::Ones());
 	sampleObject(mesh);
	PartioReaderWriter::writeParticles(outputFile, (unsigned int)particles.size(), particles.data(), NULL, radius);

	std::cout << "Generated particles: " << particles.size() << std::endl;

	Timing::printAverageTimes();

	return 0;
}

 void sampleObject(TriangleMesh &mesh)
 {
	 Vector3r *v = mesh.getVertices().data();

 	// compute bounding box
 	Vector3r min, max;
 	min = v[0];
 	max = min;
 	for (unsigned int i = 1; i < mesh.numVertices(); ++i)
 	{
 		const Vector3r& p = v[i];
 		for (unsigned int j = 0; j < 3; ++j) {
 			min[j] = std::min(min[j], p[j]);
 			max[j] = std::max(max[j], p[j]);
 		}
 	}

	if (useRegion)
 	{
		for (unsigned int i=0; i < 3; i++)
		{
			min[i] = std::max(scale * region.m_min[i], min[i]);
			max[i] = std::min(scale * region.m_max[i], max[i]);
		}
 	}
 
 	// sample object
 	const unsigned int numberOfSamplePoints = (((unsigned int)((1.0f / diameter) * (max[2] - min[2]))) + 1) *
 		(((unsigned int)((1.0f / diameter) * (max[1] - min[1]))) + 1) *
 		(((unsigned int)((1.0f / diameter) * (max[0] - min[0]))) + 1);
 	unsigned int currentSample = 0;
 	Real currentPercent = 0.01;
 	int counter_x = 0;
 	int counter_y = 0;
 	Real xshift = diameter;
 	Real yshift = diameter;
 
 	if (mode == 1)
 		yshift = sqrt(3.0) * radius;
 	else if (mode == 2)
 	{
 		xshift = sqrt(3.0) * radius;
 		yshift = sqrt(6.0) * diameter / 3.0;
 	}
 	for (Real z = min[2]; z <= max[2]; z += diameter)
 	{
 		for (Real y = min[1]; y <= max[1]; y += yshift)
 		{
 			for (Real x = min[0]; x <= max[0]; x += xshift)
 			{
 				Vector3r particlePosition;
 				if (mode == 1)
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
 						shift_vec[2] += diameter / (2.0 * (counter_y % 2 ? -1 : 1));
 					}
 					if (counter_y % 2)
 					{
 						shift_vec[0] += xshift / 2.0;
 						shift_vec[2] += diameter / 2.0;
 					}
 					particlePosition += shift_vec;
 				}
 				else
 				{
 					// Use center of voxel
 					particlePosition = Vector3r(x + radius, y + radius, z + radius);
 				}
 
 				const Real w_p = WindingNumbers::computeGeneralizedWindingNumber(particlePosition, mesh);
 				if (w_p > 0.5f)
 					particles.push_back(particlePosition);
 				currentSample++;
 
 				if ((Real)currentSample / (Real)numberOfSamplePoints > currentPercent)
 				{
 					std::cout << currentPercent * 100.0 << "%\n";
 					currentPercent += 0.01;
 				}
 				counter_x++;
 			}
 			counter_x = 0;
 			counter_y++;
 		}
 		counter_y = 0;
 	}
 	std::cout << "100%\n";
}
 
