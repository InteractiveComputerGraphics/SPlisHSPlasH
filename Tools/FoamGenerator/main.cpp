#include "SPlisHSPlasH/Common.h"
#include <Eigen/Dense>
#include <iostream>
#include "Utilities/Timing.h"
#include "Utilities/FileSystem.h"
#include "extern/cxxopts/cxxopts.hpp"
#include "CompactNSearch.h"
#include "Utilities/StringTools.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "FoamKernel.h"
#include "Utilities/PartioReaderWriter.h"
#include "extern/partio/src/lib/Partio.h"
#include <array>
#include <future>
#include <random>

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

void generateFoamFiles();
void queryValues(Real &vDiffAvgMax, Real & omegaDiffAvgMax, Real &curvatureAvgMax, Real &energyAvgMax);
void determineValues();
void generateFoam();
bool readCurrentFrame();
void writeCurrentFrame();
void writeCurrentFrame_bgeo();
void writeCurrentFrame_vtk();

bool readParticles(const std::string &fileName, std::vector<Vector3r> &positions, std::vector<Vector3r> &velocities, std::vector<Vector3r> &angularVelocities);
void computeNormals();
void computeDensities();
void generateFoamParticles(const unsigned int index, const unsigned int numParticles, const unsigned int numTrappedAir, const unsigned int numWaveCrest, const unsigned int numVorticity, const Real I_ke);
void getOrthogonalVectors(const Vector3r &vec, Vector3r &x, Vector3r &y);
Real clampAndNormalize(const Real val, const Real minVal, const Real maxVal);
void advectFoamParticles();
void removeParticles();

std::vector<Vector3r> x_copy;
std::vector<unsigned char> particleType;
enum ParticleType {Foam, Spray, Bubbles, NUM_PARTICLE_TYPES};
enum GeneratedType {TrappedAir, WaveCrest, Vorticity, NUM_GENERATOR_TYPES};
constexpr unsigned int maxNumTypes = (NUM_PARTICLE_TYPES > NUM_GENERATOR_TYPES) ? NUM_PARTICLE_TYPES : NUM_GENERATOR_TYPES;
unsigned int numTypes = 0;
std::array<unsigned int, maxNumTypes> numParticlesOfType;
std::array<std::string, maxNumTypes> nameExtensions;
std::array<std::string, maxNumTypes> extendedFileNames;
std::array<std::vector<Vector3r>, maxNumTypes> particleCopies;
std::vector<Real> v_diff, curvature, omega_diff, energy;
std::string fileName_copy;
std::future<void> handle;
std::array<std::future<void>, maxNumTypes> handlesType;
const Real density0 = 1000.0;
string input = "";
string output = "";
int output_format = 0; // 0: bgeo, 1: vtk
unsigned int startFrame = 1;
unsigned int endFrame = 0xffffffff;
unsigned int currentFrame = 1;
unsigned int skipframes = 0;
string exePath;
Real particleRadius = 0.025;
Real timeStepSize = 0.02;
Real invDt;
Real supportRadius;
Real mass = 1.0;
Real foam_scale = 1000.0;
Real taMin, taMax, wcMin, wcMax, keMin, keMax;
Real k_ta = 4000;
Real k_wc = 100000;
Real k_buoyancy = 2.0;
Real k_drag = 0.8;
Real lifetimeMin = 2.0;
Real lifetimeMax = 5.0;
Vector3r bbMin(Vector3r::Constant(-std::numeric_limits<Real>::max()));
Vector3r bbMax(Vector3r::Constant(std::numeric_limits<Real>::max()));
enum class BbType { Kill, Lifesteal, Clamp };
BbType bbType = BbType::Lifesteal;
CompactNSearch::NeighborhoodSearch *neighborhoodSearch;
std::vector<Vector3r> x0, v0, normals;
std::vector<Real> densities;
std::vector<Vector3r> fx, fv;
std::vector<Real> flifetime;
bool queryMode = false;
bool automaticMode = true;
bool splitTypes = false;
bool splitGenerators = false;
enum KernelType { CubicSpline = 0, Ihmsen2012, NumKernelTypes };
KernelType kernelType = KernelType::CubicSpline;
Real(*kernelFct)(const Vector3r &) = CubicKernel::W;

std::random_device rand_device;
std::mt19937_64 random_generator(rand_device());
std::uniform_real_distribution<Real> uniform_distr_0_1(Real(0.0), Real(1.0));
Real max_v_diff = 0.0;
Real sum_v_diff = 0.0;
Real sum_max_v_diff = 0.0;
Real max_curvature = 0.0;
Real sum_curvature = 0.0;
Real sum_max_curvature = 0.0;
Real max_energy = 0.0;
Real sum_energy = 0.0;
Real sum_max_energy = 0.0;
unsigned int numberOfFrames;
Real inertia = 2.0;
Real k_vo = 4000;
Real voMin, voMax;
std::vector<Vector3r> omega0;
Real max_omega_diff = 0.0;
Real sum_omega_diff = 0.0;
Real sum_max_omega_diff = 0.0;

std::vector<std::vector<Vector3r>> fxPerThread, fvPerThread;
std::vector<std::vector<Real>> flifetimePerThread;
std::vector<std::vector<unsigned char>> particleTypePerThread;
std::vector<std::array<unsigned int, maxNumTypes>> numParticlesOfTypePerThread;

inline unsigned int numberOfPointSets() 
{
	return static_cast<unsigned int>(neighborhoodSearch->n_point_sets());
}

inline unsigned int numberOfNeighbors(const unsigned int pointSetIndex, const unsigned int index)
{
	return static_cast<unsigned int>(neighborhoodSearch->point_set(0).n_neighbors(pointSetIndex, index));
}

inline unsigned int getNeighbor(const unsigned int pointSetIndex, const unsigned int index, const unsigned int k)
{
	return neighborhoodSearch->point_set(0).neighbor(pointSetIndex, index, k);
}

inline unsigned int numberOfNeighbors(const unsigned int pointSetIndex, const unsigned int index, const unsigned int neighborPointSetIndex)
{
	return static_cast<unsigned int>(neighborhoodSearch->point_set(pointSetIndex).n_neighbors(neighborPointSetIndex, index)); 
} 			

inline unsigned int getNeighbor(const unsigned int pointSetIndex, const unsigned int index, const unsigned int neighborPointSetIndex, const unsigned int k)
{ 
	return neighborhoodSearch->point_set(pointSetIndex).neighbor(neighborPointSetIndex, index, k); 
}


// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;

	Utilities::logger.addSink(shared_ptr<Utilities::ConsoleSink>(new Utilities::ConsoleSink(Utilities::LogLevel::INFO)));
	exePath = FileSystem::getProgramPath();

	try
	{
		cxxopts::Options options(argv[0], "FoamGen - An implementation of Bender et al. \"Turbulent Micropolar SPH Fluids with Foam\", 2018\n\n"
			"By default the limits and the factors (ta, wc, vo) are determined automatically. The amount of generated foam particles "
			"is solely determined by the parameter foamscale. If the -no_auto flag is set, all parameters can be set manually.\n"
		);

		options.add_options()
			("i,input", "Input file (partio)", cxxopts::value<std::string>())
			("o,output", "Output file (partio or vtk)", cxxopts::value<std::string>())
			("q,query", "Query mode: determines max/avg values ")
			("no-auto", "Disable automatic mode. Limits and factors ta, wc, vo must be set manually.")
			("splittypes", "Output each foam type to a different file")
			("splitgenerators", "Output different foam files depending on which potential generated the foam. Overrides --splittypes.")
			("s,startframe", "Start frame", cxxopts::value<unsigned int>()->default_value("1"))
			("e,endframe", "End frame", cxxopts::value<unsigned int>())
			("r,radius", "Particle radius", cxxopts::value<Real>()->default_value("0.025"))
			("t,timestepsize", "Time step size", cxxopts::value<Real>()->default_value("0.02"))
			("k,kernel", "0: Cubic spline, 1: Ihmsen et al. 2012", cxxopts::value<int>()->default_value("0"))
			("l,limits", "Limits (min/max) for potentials (trapped air, wave crest, vorticity, kinetic energy)", cxxopts::value<std::string>()->default_value("5,20,2,8,5,20,5,50"))
			("lifetime", "Lifetime (min/max)", cxxopts::value<std::string>()->default_value("2.0,5.0"))
			("b,buoyancy", "Buoyancy", cxxopts::value<Real>()->default_value("2.0"))
			("d,drag", "Drag", cxxopts::value<Real>()->default_value("0.8"))
			("ta", "Trapped air factor", cxxopts::value<Real>()->default_value("4000"))
			("wc", "Wave crest factor", cxxopts::value<Real>()->default_value("50000"))
			("vo", "Vorticity factor", cxxopts::value<Real>()->default_value("4000"))
			("bbsize", "minimum and maximum coordinates of and axis aligned bounding-box (minX, minY, minZ, maxX, maxY, maxZ)", cxxopts::value<std::string>())
			("bbtype", "chose how the bounding-box is used [kill | lifesteal | clamp]. Use in combination with --bbsize.", cxxopts::value<std::string>())
			("skipframes", "number of frames to skip when writing foam", cxxopts::value<unsigned int>()->default_value("0"))
			("f,foamscale", "Global multiplier for number of generated foam particles", cxxopts::value<Real>()->default_value("1000"))
			("h,help", "Print help")
			;

		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
		}

		queryMode = false;
		if (result.count("query"))
		{
			queryMode = true;
		}

		automaticMode = true;
		if (result.count("no-auto"))
		{
			automaticMode = false;
		}

		if (result.count("splitgenerators"))
		{
			splitGenerators = true;
			nameExtensions[GeneratedType::TrappedAir] = std::string("_TrappedAir");
			nameExtensions[GeneratedType::WaveCrest]  = std::string("_WaveCrest") ;
			nameExtensions[GeneratedType::Vorticity]  = std::string("_Vorticity") ;
			numTypes = GeneratedType::NUM_GENERATOR_TYPES;
		}
		if (result.count("splittypes") && !splitGenerators)
		{
			splitTypes = true;
			nameExtensions[ParticleType::Foam] = std::string("_foam");
			nameExtensions[ParticleType::Spray] = std::string("_spray");
			nameExtensions[ParticleType::Bubbles] = std::string("_bubbles");
			numTypes = ParticleType::NUM_PARTICLE_TYPES;
		}

		if (result.count("input"))
		{
			input = result["input"].as<std::string>();
			LOG_INFO << "Input = " << input;
		}
		if (result.count("output"))
		{
			output = result["output"].as<std::string>();
			LOG_INFO << "Output = " << output;

			output_format = 0;
			if (Utilities::StringTools::to_upper(FileSystem::getFileExt(output)) == "VTK")
				output_format = 1;
		}
		if (input == "")
		{
			LOG_ERR << "Input is missing!";
			exit(1);
		}
		if ((output == "") && !queryMode)
		{
			LOG_ERR << "Output is missing!";
			exit(1);
		}

		if (result.count("startframe"))
			startFrame = result["startframe"].as<unsigned int>();

		if (result.count("endframe"))
			endFrame = result["endframe"].as<unsigned int>();

		if (result.count("radius"))
			particleRadius = result["radius"].as<Real>();

		if (result.count("timestepsize"))
			timeStepSize = result["timestepsize"].as<Real>();

		if (result.count("kernel"))
		{
			int kt = result["kernel"].as<int>();
			if ((kt < 0) || (kt >= KernelType::NumKernelTypes))
			{
				kt = 0;
				LOG_WARN << "Wrong kernel type. Using cubic spline kernel.";
			}
			kernelType = static_cast<KernelType>(kt);

			if (kernelType == KernelType::CubicSpline)
				kernelFct = CubicKernel::W;
			else if (kernelType == KernelType::Ihmsen2012)
				kernelFct = FoamKernel::W;
		}

		if (result.count("buoyancy"))
			k_buoyancy = result["buoyancy"].as<Real>();

		if (result.count("drag"))
			k_drag = result["drag"].as<Real>();

		if (result.count("skipframes"))
			skipframes = result["skipframes"].as<unsigned int>();

		foam_scale = result["foamscale"].as<Real>();

		k_ta = result["ta"].as<Real>();
		k_wc = result["wc"].as<Real>();
		k_vo = result["vo"].as<Real>();
		{
			std::vector<string> tokens;
			std::string limitsStr = result["limits"].as<std::string>();
			StringTools::tokenize(limitsStr, tokens, ",");
			if (tokens.size() != 8)
			{
				LOG_ERR << "Wrong number of limit values!";
				exit(1);
			}
			taMin = stof(tokens[0]);
			taMax = stof(tokens[1]);
			wcMin = stof(tokens[2]);
			wcMax = stof(tokens[3]);
			voMin = stof(tokens[4]);
			voMax = stof(tokens[5]);
			keMin = stof(tokens[6]);
			keMax = stof(tokens[7]);
			LOG_INFO << "Limits: " << taMin << ", " << taMax << ", "
				<< wcMin << ", " << wcMax << ", "
				<< voMin << ", " << voMax << ", "
				<< keMin << ", " << keMax;
		}

		{
			std::vector<string> tokens;
			std::string ltStr = result["lifetime"].as<std::string>();
			StringTools::tokenize(ltStr, tokens, ",");
			if (tokens.size() != 2)
			{
				LOG_ERR << "Wrong number of lifetime values!";
				exit(1);
			}
			lifetimeMin = stof(tokens[0]);
			lifetimeMax = stof(tokens[1]);
			LOG_INFO << "Lifetime (min/max): " << lifetimeMin << ", " << lifetimeMax;
		}

		if (result.count("bbsize"))
		{
			std::vector<std::string> tokens;
			std::string limitsStr = result["bbsize"].as<std::string>();
			StringTools::tokenize(limitsStr, tokens, ",");
			if (tokens.size() != 6)
			{
				std::cerr << "Wrong number of bounding box values!\n";
				exit(1);
			}
			bbMin = Vector3r(stof(tokens[0]), stof(tokens[1]), stof(tokens[2]));
			bbMax = Vector3r(stof(tokens[3]), stof(tokens[4]), stof(tokens[5]));
			if ((bbMin.array() > bbMax.array()).any())
			{
				std::cerr << "Bounding-box: A min value is larger than its corresponding max value!\n";
				exit(1);
			}
		}
		if (result.count("bbtype"))
		{
			if (!result.count("bbsize"))
			{
				std::cerr << "No bounding-box given! Use --bbsize to set the size of the bounding box.\n";
				exit(1);
			}
			auto s = result["bbtype"].as<std::string>();
				if (s.compare("kill") == 0)
			{
				bbType = BbType::Kill;
			}
			else if (s.compare("lifesteal") == 0)
			{
				bbType = BbType::Lifesteal;
			}
			else if (s.compare("clamp") == 0)
			{
				bbType = BbType::Clamp;
			}
			else
			{
				std::cerr << "Unknown bounding-box type: " << s << "\n";
				exit(1);
			}
		}
		LOG_INFO << "Particle radius: " << particleRadius;
		LOG_INFO << "Drag: " << k_drag;
		LOG_INFO << "Buoyancy: " << k_buoyancy;
		LOG_INFO << "Trapped air factor: " << k_ta;
		LOG_INFO << "Wave crest factor: " << k_wc;
		LOG_INFO << "Vorticity factor: " << k_vo;
		LOG_INFO << "Bounding box: " << bbMin.transpose() << ", " << bbMax.transpose();
		LOG_INFO << "Bounding box type: " << (int) bbType;
		LOG_INFO << "Foam scale: " << foam_scale;
		LOG_INFO << "Time step size: " << timeStepSize;
		if (skipframes > 0)
			std::cout << "Skipping " << skipframes << " frame(s) between writes." << std::endl;
	}
	catch (const cxxopts::exceptions::exception& e)
	{
		LOG_ERR << "error parsing options: " << e.what();
		exit(1);
	}

#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
#else
	const int maxThreads = 1;
#endif

	fxPerThread.resize(maxThreads); 
	fvPerThread.resize(maxThreads);
	flifetimePerThread.resize(maxThreads); 
	particleTypePerThread.resize(maxThreads);
	numParticlesOfTypePerThread.resize(maxThreads);

	// init particle mass
	const Real diam = static_cast<Real>(2.0) * particleRadius;
	mass = static_cast<Real>(0.8) * diam*diam*diam * density0;

	// if in query mode, determine the max values per frame 
	Real vDiffAvgMax, omegaDiffAvgMax, curvatureAvgMax, energyAvgMax;
	if (queryMode)
		queryValues(vDiffAvgMax, omegaDiffAvgMax, curvatureAvgMax, energyAvgMax);
	else if (automaticMode)
	{
		queryValues(vDiffAvgMax, omegaDiffAvgMax, curvatureAvgMax, energyAvgMax);

		taMin = static_cast<Real>(0.1)*vDiffAvgMax;
		taMax = vDiffAvgMax;
		wcMin = static_cast<Real>(0.1)*curvatureAvgMax;
		wcMax = curvatureAvgMax;
		voMin = static_cast<Real>(0.1)*omegaDiffAvgMax;
		voMax = omegaDiffAvgMax;
		keMin = static_cast<Real>(0.1)*energyAvgMax;
		keMax = energyAvgMax;
		LOG_INFO << "Limits: " << taMin << ", " << taMax << ", "
			<< wcMin << ", " << wcMax << ", "
			<< voMin << ", " << voMax << ", "
			<< keMin << ", " << keMax;

		k_ta = 1.0;
		k_wc = 1.0;
		k_vo = 1.0;
		
		generateFoamFiles();
	}
	else
		generateFoamFiles();

	Timing::printAverageTimes();
	Timing::printTimeSums();

	return 0;
}

/** Read the sequence of particle data files and determine for each frame the values 
* that are required to decide how many foam particles should be generated. 
* Finally, sum up the max values and compute the average wrt the number offrames. */
void queryValues(Real &vDiffAvgMax, Real & omegaDiffAvgMax, Real &curvatureAvgMax, Real &energyAvgMax)
{
	currentFrame = startFrame;
	numberOfFrames = 0;

	// Initialize neighborhood search
	supportRadius = static_cast<Real>(4.0) * particleRadius;

	// Initialize kernel
	FoamKernel::setRadius(supportRadius);
	CubicKernel::setRadius(supportRadius);

	invDt = static_cast<Real>(1.0) / timeStepSize;

	bool chk = true;
	bool pointSetAdded = false;
	START_TIMING("total");
	while (chk)
	{
		LOG_INFO << "Reading frame: " << currentFrame;
		chk = readCurrentFrame();
		if (!chk)
			LOG_ERR << "Failed to read frame";

		START_TIMING("iteration");
		if (neighborhoodSearch == nullptr)
		{
			neighborhoodSearch = new CompactNSearch::NeighborhoodSearch(supportRadius, false);
			neighborhoodSearch->set_radius(supportRadius);
			// Fluids 
			neighborhoodSearch->add_point_set(&x0[0][0], x0.size(), true, true);
		}
		else
			neighborhoodSearch->resize_point_set(0, &x0[0][0], x0.size());

		START_TIMING("neighborhoodSearch");
		neighborhoodSearch->find_neighbors();
		STOP_TIMING_AVG

		START_TIMING("determineValues");
		determineValues();
		STOP_TIMING_AVG

		STOP_TIMING_AVG

		numberOfFrames++;

		if (numberOfFrames % 10 == 0)
		{
			LOG_INFO << "v_diff     - max.: " << max_v_diff << ", avg.: " << sum_v_diff / (Real)numberOfFrames;
			LOG_INFO << "Curvature  - max.: " << max_curvature << ", avg.: " << sum_curvature / (Real)numberOfFrames;
			LOG_INFO << "Omega_diff - max.: " << max_omega_diff << ", avg.: " << sum_omega_diff / (Real)numberOfFrames;
			LOG_INFO << "Energy     - max.: " << max_energy << ", avg.: " << sum_energy / (Real)numberOfFrames;
		}

		currentFrame++;
		if (currentFrame > endFrame)
			break;
	}
	delete neighborhoodSearch;
	neighborhoodSearch = nullptr;

	STOP_TIMING_AVG;
	// average maxima
	vDiffAvgMax = sum_max_v_diff / (Real)numberOfFrames;
	omegaDiffAvgMax = sum_max_omega_diff / (Real)numberOfFrames;
	curvatureAvgMax = sum_max_curvature / (Real)numberOfFrames;
	energyAvgMax = sum_max_energy / (Real)numberOfFrames;
	LOG_INFO << "v_diff     - avg. max.: " << vDiffAvgMax;
	LOG_INFO << "Curvature  - avg. max.: " << curvatureAvgMax;
	LOG_INFO << "Omega_diff - avg. max.: " << omegaDiffAvgMax;
	LOG_INFO << "Energy     - avg. max.: " << energyAvgMax;
}


/** Compute values that are required to decide how many foam particles should be generated.\n\n
* 
* See:\n
* - Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012 \n
* - Bender et al., "Turbulent Micropolar SPH Fluids with Foam", 2018 
*/
void determineValues()
{
	// compute the density for each particle
	computeDensities();
	// compute the particle normals which are required for the curvature term
	computeNormals();

	Real sumVDiff = 0.0;
	Real sumCurvature = 0.0;
	Real sumEnergy = 0.0;
	Real maxVDiff = 0.0;
	Real maxCurvature = 0.0;
	Real maxEnergy = 0.0;
	Real sumOmegaDiff = 0.0;
	Real maxOmegaDiff = 0.0;

	v_diff.resize(x0.size());
	omega_diff.resize(x0.size());
	curvature.resize(x0.size());
	energy.resize(x0.size());

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)x0.size(); i++)
		{
			const Vector3r &xi = x0[i];
			const Vector3r &vi = v0[i];

			v_diff[i] = 0.0;
			curvature[i] = 0.0;
			omega_diff[i] = 0.0;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(0, i, j);
				const Vector3r &xj = x0[neighborIndex];
				const Vector3r &vj = v0[neighborIndex];

				//////////////////////////////////////////////////////////////////////////
				// Trapped air potential 
				// Compute Eq. 2 in Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
				//////////////////////////////////////////////////////////////////////////
				Vector3r vivj = vi - vj;
				const Real vivjmag = vivj.norm();
				if (vivjmag > 1.0e-6)
					vivj = (1.0 / vivjmag) * vivj;

				Vector3r xixj = xi - xj;
				xixj.normalize();

				if (kernelType != KernelType::Ihmsen2012)
					v_diff[i] += mass / densities[j] * vivjmag * (static_cast<Real>(1.0) - vivj.dot(xixj)) * kernelFct(xi - xj);
				else
					v_diff[i] += vivjmag * (static_cast<Real>(1.0) - vivj.dot(xixj)) * FoamKernel::W(xi - xj);

				//////////////////////////////////////////////////////////////////////////
				// Wave crest curvature
				// Compute Eq. 4 in Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
				//////////////////////////////////////////////////////////////////////////
				if (-xixj.dot(normals[i]) < 0)
				{
					if (kernelType != KernelType::Ihmsen2012)
						curvature[i] += mass / densities[j] * (static_cast<Real>(1.0) - normals[i].dot(normals[neighborIndex])) * kernelFct(xi - xj);
					else
						curvature[i] += (static_cast<Real>(1.0) - normals[i].dot(normals[neighborIndex])) * FoamKernel::W(xi - xj);
				}

				//////////////////////////////////////////////////////////////////////////
				// Vorticity
				// Compute omaga^diff in Bender et al., "Turbulent Micropolar SPH Fluids with Foam", 2018
				//////////////////////////////////////////////////////////////////////////
				if (omega0.size() > 0)
				{
					const Vector3r &omegai = omega0[i];
					const Vector3r &omegaj = omega0[neighborIndex];
					if (kernelType != KernelType::Ihmsen2012)
						omega_diff[i] += mass / densities[j] * (omegai - omegaj).norm() * kernelFct(xi - xj);
					else
						omega_diff[i] += (omegai - omegaj).norm() * FoamKernel::W(xi - xj);
				}
			}

			//////////////////////////////////////////////////////////////////////////
			// Wave crest
			// Compute Eq. 7 in Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
			//////////////////////////////////////////////////////////////////////////
			Real delta = 0.0;
			Vector3r vi_normalized = vi;
			vi_normalized.normalize();
			if (vi_normalized.dot(normals[i]) >= 0.6)
				delta = 1.0;

			curvature[i] = delta * curvature[i];

			//////////////////////////////////////////////////////////////////////////
			// Kinetic energic
			//////////////////////////////////////////////////////////////////////////
			if (omega0.size() > 0)
				energy[i] = static_cast<Real>(0.5)*mass*vi.squaredNorm() + static_cast<Real>(0.5) * inertia * omega0[i].squaredNorm();
			else
				energy[i] = static_cast<Real>(0.5)*mass*vi.squaredNorm();
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Compute sum and max of all values
	//////////////////////////////////////////////////////////////////////////
	for (auto i = 0; i < x0.size(); i++)
	{
		sumVDiff += v_diff[i];
		maxVDiff = std::max(maxVDiff, v_diff[i]);

		if (omega0.size() > 0)
		{
			sumOmegaDiff += omega_diff[i];
			maxOmegaDiff = std::max(maxOmegaDiff, omega_diff[i]);
		}

		sumCurvature += curvature[i];
		maxCurvature = std::max(maxCurvature, curvature[i]);

		sumEnergy += energy[i];
		maxEnergy = std::max(maxEnergy, energy[i]);
	}

	sum_v_diff += sumVDiff / (Real)x0.size();
	max_v_diff = std::max(max_v_diff, maxVDiff);
	sum_max_v_diff += maxVDiff;

	if (omega0.size() > 0)
	{
		sum_omega_diff += sumOmegaDiff / (Real)x0.size();
		max_omega_diff = std::max(max_omega_diff, maxOmegaDiff);
		sum_max_omega_diff += maxOmegaDiff;
	}

	sum_curvature += sumCurvature / (Real)x0.size();
	max_curvature = std::max(max_curvature, maxCurvature);
	sum_max_curvature += maxCurvature;

	sum_energy += sumEnergy / (Real)x0.size();
	max_energy = std::max(max_energy, maxEnergy);
	sum_max_energy += maxEnergy;
}

/** Generate a foam file for each frame.
 */
void generateFoamFiles()
{
	currentFrame = startFrame;

	// Initialize kernels
	supportRadius = static_cast<Real>(4.0) * particleRadius;
	FoamKernel::setRadius(supportRadius);
	CubicKernel::setRadius(supportRadius);

	invDt = static_cast<Real>(1.0) / timeStepSize;

#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
#else
	const int maxThreads = 1;
#endif


	bool chk = true;
	bool first = true;
	bool pointSetAdded = false;
	START_TIMING("total");
	while (chk)
	{
		// Read current frame of fluid particles
		LOG_INFO << "Reading frame: " << currentFrame;
		chk = readCurrentFrame();

		if (first)
		{
			// init foam particle vectors
			fx.reserve(10 * x0.size());
			fv.reserve(10 * x0.size());
			flifetime.reserve(10 * x0.size());
			first = false;

			for (int i = 0; i < maxThreads; i++)
			{
				fxPerThread.reserve(10 * x0.size());
				fvPerThread.reserve(10 * x0.size());
				flifetimePerThread.reserve(10 * x0.size());
			}
		}

		// init neighborhood search
		if (neighborhoodSearch == NULL)
		{
			neighborhoodSearch = new CompactNSearch::NeighborhoodSearch(supportRadius, false);
			neighborhoodSearch->set_radius(supportRadius);
			// Fluid 
			neighborhoodSearch->add_point_set(&x0[0][0], x0.size(), true, true);
			neighborhoodSearch->set_active(0u, 0u, true);
			pointSetAdded = true;
		}
		else
			neighborhoodSearch->resize_point_set(0, &x0[0][0], x0.size());

		// find the fluid neighbors of each fluid particle
		START_TIMING("neighborhoodSearch - advection");
		neighborhoodSearch->find_neighbors();
		STOP_TIMING_AVG

		// remove foam particles which exceeded their lifetime
		START_TIMING("removeParticles");
		removeParticles();
		STOP_TIMING_AVG;

		// advect each foam particle by performing the time integration
		START_TIMING("advectFoamParticles");
		advectFoamParticles();
		STOP_TIMING_AVG

		// write the foam particles
		START_TIMING("writeCurrentFrame");
		writeCurrentFrame();
		STOP_TIMING_AVG;

		// generate new foam particles 
		START_TIMING("generateFoam");
		generateFoam();
		STOP_TIMING_AVG;

		currentFrame++;
		if (currentFrame > endFrame)
			break;
	}

	delete neighborhoodSearch;
	neighborhoodSearch = NULL;

	STOP_TIMING_AVG;
}

/** Generate new foam particles using the algorithm described in :\n
* - Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012 \n
* - Bender et al., "Turbulent Micropolar SPH Fluids with Foam", 2018 
*/
void generateFoam()
{
	// Compute the density for each fluid particle
	computeDensities();
	// Compute the normals for the fluid particles which is needed to compute the curvature
	computeNormals();

#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
#else
	const int maxThreads = 1;
#endif

	if (splitGenerators)
	{
		for (int g = 0; g < NUM_GENERATOR_TYPES; ++g)
		{
			numParticlesOfType[g] = 0;
			for (int i=0; i < maxThreads; i++)
			{
				numParticlesOfTypePerThread[i][g] = 0;
			}
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)x0.size(); i++)
		{
			const Vector3r &xi = x0[i];
			const Vector3r &vi = v0[i]; 

			Real v_diff = 0.0;
			Real curvature = 0.0;
			Real omega_diff = 0.0;

			// Only generate foam if we have enough neighbors
			if (numberOfNeighbors(0, i) < 15)
				continue;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(0, i, j);
				const Vector3r &xj = x0[neighborIndex];
				const Vector3r &vj = v0[neighborIndex];

				//////////////////////////////////////////////////////////////////////////
				// Trapped air potential 
				// Compute Eq. 2 in Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
				//////////////////////////////////////////////////////////////////////////
				Vector3r vivj = vi - vj;
				const Real vivjmag = vivj.norm();
				if (vivjmag > 1.0e-6)
					vivj = (1.0 / vivjmag) * vivj;

				Vector3r xixj = xi - xj;
				xixj.normalize();

				if (kernelType != KernelType::Ihmsen2012)
					v_diff += mass / densities[j] * vivjmag * (static_cast<Real>(1.0) - vivj.dot(xixj)) * kernelFct(xi - xj);
				else
					v_diff += vivjmag * (static_cast<Real>(1.0) - vivj.dot(xixj)) * FoamKernel::W(xi - xj);

				//////////////////////////////////////////////////////////////////////////
				// Wave crest curvature
				// Compute Eq. 4 in Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
				//////////////////////////////////////////////////////////////////////////
				if (-xixj.dot(normals[i]) < 0)
				{
					if (kernelType != KernelType::Ihmsen2012)
						curvature += mass / densities[j] * (static_cast<Real>(1.0) - normals[i].dot(normals[neighborIndex])) * kernelFct(xi - xj);
					else
						curvature += (static_cast<Real>(1.0) - normals[i].dot(normals[neighborIndex])) * FoamKernel::W(xi - xj);
				}

				//////////////////////////////////////////////////////////////////////////
				// Vorticity
				// Compute omaga^diff in Bender et al., "Turbulent Micropolar SPH Fluids with Foam", 2018
				//////////////////////////////////////////////////////////////////////////
				if (omega0.size() > 0)
				{
					const Vector3r &omegai = omega0[i];
					const Vector3r &omegaj = omega0[neighborIndex];
					if (kernelType != KernelType::Ihmsen2012)
						omega_diff += mass / densities[j] * (omegai - omegaj).norm() * kernelFct(xi - xj);
					else
						omega_diff += (omegai - omegaj).norm() * FoamKernel::W(xi - xj);
				}
			}

			//////////////////////////////////////////////////////////////////////////
			// Trapped air potential 
			//////////////////////////////////////////////////////////////////////////
			const Real I_ta = clampAndNormalize(v_diff, taMin, taMax);

			//////////////////////////////////////////////////////////////////////////
			// Vorticity potential 
			//////////////////////////////////////////////////////////////////////////
			Real I_vo = 0.0;
			if (omega0.size() > 0)
				I_vo = clampAndNormalize(omega_diff, voMin, voMax);

			//////////////////////////////////////////////////////////////////////////
			// Wave crest
			// Compute Eq. 7 in Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
			//////////////////////////////////////////////////////////////////////////
			Real delta = 0.0;
			Vector3r vi_normalized = vi;
			vi_normalized.normalize();
			if (vi_normalized.dot(normals[i]) >= 0.6)
				delta = 1.0;
			const Real I_wc = clampAndNormalize(delta * curvature, wcMin, wcMax);

			//////////////////////////////////////////////////////////////////////////
			// Kinetic energy
			//////////////////////////////////////////////////////////////////////////
			Real I_ke = 0.0;
			unsigned int nd = 0;
			unsigned int nv = 0;
			if (omega0.size() > 0)
			{
				I_ke = clampAndNormalize(static_cast<Real>(0.5)*mass*vi.squaredNorm() + static_cast<Real>(0.5) * inertia * omega0[i].squaredNorm(), keMin, keMax);
				nd = (unsigned int)(foam_scale * I_ke * (k_ta * I_ta + k_wc * I_wc + k_vo * I_vo) * timeStepSize + 0.5);
				nv = (unsigned int)(foam_scale * I_ke * k_vo * I_vo * timeStepSize + 0.5);
			}
			else
			{
				I_ke = clampAndNormalize(static_cast<Real>(0.5)*mass*vi.squaredNorm(), keMin, keMax);
				nd = (unsigned int) std::max((foam_scale * I_ke * (k_ta * I_ta + k_wc * I_wc) * timeStepSize + static_cast<Real>(0.5)), static_cast<Real>(0.0));
				nv = 0;
			}

			const unsigned int nt = (unsigned int)(foam_scale * I_ke * k_ta * I_ta * timeStepSize + 0.5);
			const unsigned int nw = (unsigned int)(foam_scale * I_ke * k_wc * I_wc * timeStepSize + 0.5);

			//////////////////////////////////////////////////////////////////////////
			// Generate foam particles
			////////////////////////////////////////////////////////////////////////// 
			if (splitGenerators)
				generateFoamParticles(i, nt + nw + nv, nt, nw, nv, I_ke);
			else
				generateFoamParticles(i, nd, -1, -1, -1, I_ke);
		}
	}

	// Merge particles that have been generated in different threads
	for (int j = 0; j < maxThreads; j++)
	{
		for (int i = 0; i < fxPerThread[j].size(); i++)
		{
			fx.push_back(fxPerThread[j][i]);
			fv.push_back(fvPerThread[j][i]);
			flifetime.push_back(flifetimePerThread[j][i]);
			particleType.push_back(particleTypePerThread[j][i]);
		}
		fxPerThread[j].clear();
		fvPerThread[j].clear();
		flifetimePerThread[j].clear();
		particleTypePerThread[j].clear();

		for (int i = 0; i < 3; i++)
		{
			numParticlesOfType[i] += numParticlesOfTypePerThread[j][i];
		}
	}

	LOG_INFO << "# foam particles: " << fx.size();
}

std::string zeroPadding(const unsigned int number, const unsigned int length) 
{
	ostringstream out;
	out << std::internal << std::setfill('0') << std::setw(length) << number;
	return out.str();
}

/** Substitute the placeholder in the file name with the current frame number.
*/
std::string convertFileName(const std::string &inputFileName, const unsigned int currentFrame)
{
	std::string fileName = inputFileName;
	std::string::size_type pos1 = fileName.find_first_of("#", 0);
	if (pos1 == std::string::npos)
	{
		LOG_ERR << "# missing in file name.";
		exit(1);
	}
	std::string::size_type pos2 = fileName.find_first_not_of("#", pos1);
	std::string::size_type length = pos2 - pos1;

	std::string numberStr = zeroPadding(currentFrame, (unsigned int)length);
	fileName.replace(pos1, length, numberStr);
	return fileName;
}

/** Read fluid particles of current frame. 
*/
bool readCurrentFrame()
{
	std::string fileName = convertFileName(input, currentFrame);

	x0.clear();
	v0.clear();
	omega0.clear();
	return readParticles(fileName, x0, v0, omega0);
}

void writeCurrentFrame()
{
	if (output_format == 1)
		writeCurrentFrame_vtk();
	else
		writeCurrentFrame_bgeo();
}

/** Write foam particles of current frame.
*/
void writeCurrentFrame_bgeo()
{
	if (currentFrame % (1 + skipframes))
	{
		std::cout << "Skipping write of frame: " << currentFrame << std::endl;
		return;
	}
	else
	{
		std::cout << "Writing frame: " << currentFrame << std::endl;
	}

	std::string fileName = output;
	FileSystem::makeDir(FileSystem::getFilePath(FileSystem::normalizePath(output)));

	// local references for lambda capture
	Utilities::Logger & logger = Utilities::logger;
	const auto & numOfType = numParticlesOfType;
	auto & handles = handlesType;
	const auto & particleData = fx;
	const auto & particleTypes = particleType;
	const auto & nameExt = nameExtensions;
	// helper lambda to write particles depending on type
	auto writeByType = [&](const char type) {		
		extendedFileNames[type] = convertFileName(fileName, currentFrame / (1 + skipframes));
		const std::string::size_type pos = extendedFileNames[type].rfind(".");
		if (pos != std::string::npos)
			extendedFileNames[type].insert(pos, nameExt[type]);

		LOG_INFO << "Writing: " << extendedFileNames[type];
		if (handles[type].valid())
			handles[type].wait();
		particleCopies[type].clear();
		particleCopies[type].reserve(numOfType[type]);
		for (int i = 0; i < (int)particleData.size(); i++)
		{
			if (particleTypes[i] == type)
				particleCopies[type].push_back(particleData[i]);
		}
		handles[type] = std::async(std::launch::async, [type] { PartioReaderWriter::writeParticles(extendedFileNames[type], (unsigned int)particleCopies[type].size(), particleCopies[type].data(), nullptr, 0.0); });
	};

	if (splitGenerators)
	{
		writeByType(GeneratedType::TrappedAir);
		writeByType(GeneratedType::WaveCrest );
		writeByType(GeneratedType::Vorticity );
	}
	else if (splitTypes)
	{
		writeByType(ParticleType::Foam   );
		writeByType(ParticleType::Spray  );
		writeByType(ParticleType::Bubbles);
	}
	else
	{
		std::string fileName = convertFileName(output, currentFrame / (1 + skipframes));

		LOG_INFO << "Writing: " << fileName;

		if (handle.valid())
			handle.wait();

		x_copy.resize(fx.size());
		//v_copy.resize(fv.size());
		fileName_copy = fileName;
		for (int i = 0; i < (int)fx.size(); i++)
		{
			x_copy[i] = fx[i];
		}

		handle = std::async(std::launch::async, [] { PartioReaderWriter::writeParticles(fileName_copy, (unsigned int)x_copy.size(), x_copy.data(), nullptr, 0.0); });
	}
}

// VTK expects big endian
template<typename T>
inline void swapByteOrder(T* v)
{
	constexpr size_t n = sizeof(T);
	uint8_t* bytes = reinterpret_cast<uint8_t*>(v);
	for (unsigned int c = 0u; c < n / 2; c++)
		std::swap(bytes[c], bytes[n - c - 1]);
}

void writeParticlesVTK(const std::string& fileName, const unsigned int numParticles, const Vector3r* particlePositions)
{
	if (0 == numParticles)
		return;

#ifdef USE_DOUBLE
	const char* real_str = " double\n";
#else 
	const char* real_str = " float\n";
#endif

	// Open the file
	std::ofstream outfile{ fileName, std::ios::binary };
	if (!outfile.is_open())
	{
		LOG_WARN << "Cannot open a file to save VTK particles.";
		return;
	}

	outfile << "# vtk DataFile Version 4.1\n";
	outfile << "SPlisHSPlasH particle data\n"; // title of the data set, (any string up to 256 characters+\n)
	outfile << "BINARY\n";
	outfile << "DATASET UNSTRUCTURED_GRID\n";

	//////////////////////////////////////////////////////////////////////////
	// export position attribute as POINTS
	{
		std::vector<Vector3r> positions;
		positions.reserve(numParticles);
		for (unsigned int i = 0u; i < numParticles; i++)
			positions.emplace_back(particlePositions[i]);
		// swap endianess
		for (unsigned int i = 0; i < numParticles; i++)
			for (unsigned int c = 0; c < 3; c++)
				swapByteOrder(&positions[i][c]);
		// export to vtk
		outfile << "POINTS " << numParticles << real_str;
		outfile.write(reinterpret_cast<char*>(positions[0].data()), 3 * numParticles * sizeof(Real));
		outfile << "\n";
	}

	//////////////////////////////////////////////////////////////////////////
	// export particle IDs as CELLS
	{
		std::vector<Eigen::Vector2i> cells;
		cells.reserve(numParticles);
		unsigned int nodes_per_cell_swapped = 1;
		swapByteOrder(&nodes_per_cell_swapped);
		for (unsigned int i = 0u; i < numParticles; i++)
		{
			unsigned int idSwapped = i;
			swapByteOrder(&idSwapped);
			cells.emplace_back(nodes_per_cell_swapped, idSwapped);
		}

		// particles are cells with one element and the index of the particle
		outfile << "CELLS " << numParticles << " " << 2 * numParticles << "\n";
		outfile.write(reinterpret_cast<char*>(cells[0].data()), 2 * numParticles * sizeof(unsigned int));
		outfile << "\n";
	}
	//////////////////////////////////////////////////////////////////////////
	// export cell types
	{
		// the type of a particle cell is always 1
		std::vector<int> cellTypes;
		int cellTypeSwapped = 1;
		swapByteOrder(&cellTypeSwapped);
		cellTypes.resize(numParticles, cellTypeSwapped);
		outfile << "CELL_TYPES " << numParticles << "\n";
		outfile.write(reinterpret_cast<char*>(cellTypes.data()), numParticles * sizeof(int));
		outfile << "\n";
	}

	//////////////////////////////////////////////////////////////////////////
	// write additional attributes as per-particle data
	{
		outfile << "POINT_DATA " << numParticles << "\n";
		// write IDs
		outfile << "SCALARS id unsigned_int 1\n";
		outfile << "LOOKUP_TABLE id_table\n";
		// copy data
		std::vector<unsigned int> attrData;
		attrData.reserve(numParticles);
		for (unsigned int i = 0u; i < numParticles; i++)
			attrData.emplace_back(i);
		// swap endianess
		for (unsigned int i = 0; i < numParticles; i++)
			swapByteOrder(&attrData[i]);
		// export to vtk
		outfile.write(reinterpret_cast<char*>(attrData.data()), numParticles * sizeof(unsigned int));
		outfile << "\n";
	}

	outfile.close();
}

/** Write foam particles of current frame.
*/
void writeCurrentFrame_vtk()
{
	if (currentFrame % (1 + skipframes))
	{
		std::cout << "Skipping write of frame: " << currentFrame << std::endl;
		return;
	}
	else
	{
		std::cout << "Writing frame: " << currentFrame << std::endl;
	}

	std::string fileName = output;
	FileSystem::makeDir(FileSystem::getFilePath(FileSystem::normalizePath(output)));

	// local references for lambda capture
	Utilities::Logger& logger = Utilities::logger;
	const auto& numOfType = numParticlesOfType;
	auto& handles = handlesType;
	const auto& particleData = fx;
	const auto& particleTypes = particleType;
	const auto& nameExt = nameExtensions;
	// helper lambda to write particles depending on type
	auto writeByType = [&](const char type) {
		extendedFileNames[type] = convertFileName(fileName, currentFrame / (1 + skipframes));
		const std::string::size_type pos = extendedFileNames[type].rfind(".");
		if (pos != std::string::npos)
			extendedFileNames[type].insert(pos, nameExt[type]);

		LOG_INFO << "Writing: " << extendedFileNames[type];
		if (handles[type].valid())
			handles[type].wait();
		particleCopies[type].clear();
		particleCopies[type].reserve(numOfType[type]);
		for (int i = 0; i < (int)particleData.size(); i++)
		{
			if (particleTypes[i] == type)
				particleCopies[type].push_back(particleData[i]);
		}
		handles[type] = std::async(std::launch::async, [type] { writeParticlesVTK(extendedFileNames[type], (unsigned int)particleCopies[type].size(), particleCopies[type].data()); });
	};

	if (splitGenerators)
	{
		writeByType(GeneratedType::TrappedAir);
		writeByType(GeneratedType::WaveCrest);
		writeByType(GeneratedType::Vorticity);
	}
	else if (splitTypes)
	{
		writeByType(ParticleType::Foam);
		writeByType(ParticleType::Spray);
		writeByType(ParticleType::Bubbles);
	}
	else
	{
		std::string fileName = convertFileName(output, currentFrame / (1 + skipframes));

		LOG_INFO << "Writing: " << fileName;

		if (handle.valid())
			handle.wait();

		x_copy.resize(fx.size());
		//v_copy.resize(fv.size());
		fileName_copy = fileName;
		for (int i = 0; i < (int)fx.size(); i++)
		{
			x_copy[i] = fx[i];
		}

		handle = std::async(std::launch::async, [] { writeParticlesVTK(fileName_copy, (unsigned int)x_copy.size(), x_copy.data()); });
	}
}

/** Clamp and normalize value. Compute Eq. 1 in 
* Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
*/
Real clampAndNormalize(const Real val, const Real minVal, const Real maxVal)
{
	return (std::min(val, maxVal) - std::min(val, minVal)) / (maxVal - minVal);
}

/** Compute normals of fluid particles using a color field.
*/
void computeNormals()
{
	normals.resize(x0.size());

	// Compute normals
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)x0.size(); i++)
		{
			const Vector3r &xi = x0[i];
			Vector3r &ni = normals[i];
			ni.setZero();
			
			// We are only interested in surface particles
			if (numberOfNeighbors(0, i) > 20)
				continue;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(0, i, j);
				const Vector3r &xj = x0[neighborIndex];
				const Real density_j = densities[neighborIndex];
				ni -= mass / density_j * CubicKernel::gradW(xi - xj);
			}
			ni.normalize();
		}
	}

}

/** Compute densities of the fluid particles.
*/
void computeDensities()
{
	densities.resize(x0.size());

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)x0.size(); i++)
		{
			Real &density = densities[i];

			// Compute current density for particle i
			density = mass * CubicKernel::W_zero();
			const Vector3r &xi = x0[i];

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(0, i, j);
				const Vector3r &xj = x0[neighborIndex];
				density += mass * CubicKernel::W(xi - xj);
			}
		}
	}
}

/** Returns two orthogonal vectors to vec which are also orthogonal to each other.
*/
void getOrthogonalVectors(const Vector3r &vec, Vector3r &x, Vector3r &y)
{
	// Get plane vectors x, y
	Vector3r v(1, 0, 0);

	// Check, if v has same direction as vec
	if (fabs(v.dot(vec)) > 0.999)
		v = Vector3r(0, 1, 0);

	x = vec.cross(v);
	y = vec.cross(x);
	x.normalize();
	y.normalize();
}

/** Generate new foam particles for the fluid particle with the given index. 
* The particles are generated in a cylinder as described in 
* Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
*/
void generateFoamParticles(const unsigned int index, const unsigned int numParticles, const unsigned int numTrappedAir, const unsigned int numWaveCrest, const unsigned int numVorticity, const Real I_ke)
{
	Vector3r e1, e2;
	const Vector3r v = v0[index]; 
	Vector3r vn = v; 
	vn.normalize();
	getOrthogonalVectors(vn, e1, e2);

	e1 = particleRadius * e1;
	e2 = particleRadius * e2;

	for (unsigned int i = 0; i < numParticles; i++)
	{
		// Generate a random distribution of the foam particles in a cylinder.
		const Real Xr = uniform_distr_0_1(random_generator);
		const Real Xtheta = uniform_distr_0_1(random_generator);
		const Real Xh = uniform_distr_0_1(random_generator);

		const Real r = particleRadius * sqrt(Xr);
		const Real theta = Xtheta * static_cast<Real>(2.0 * M_PI);
		const Real h = (Xh- static_cast<Real>(0.5)) * timeStepSize * v.norm();

		const Vector3r xd = x0[index] + r * cos(theta) * e1 + r * sin(theta) * e2 + h * vn;
		const Vector3r vd = r * cos(theta) * e1 + r * sin(theta) * e2 + v;
		
		unsigned char generatorType = GeneratedType::TrappedAir;
		if (i >= numTrappedAir)
			generatorType = GeneratedType::WaveCrest;
		if (i >= (numTrappedAir + numWaveCrest))
			generatorType = GeneratedType::Vorticity;


#ifdef _OPENMP
		int tid = omp_get_thread_num();
#else
		int tid = 0;
#endif
		fxPerThread[tid].push_back(xd);
		fvPerThread[tid].push_back(vd);
		if (splitGenerators)
		{
			particleTypePerThread[tid].push_back(generatorType);
			numParticlesOfTypePerThread[tid][generatorType]++;
		}
		else
		{
			particleTypePerThread[tid].push_back(ParticleType::Foam); // default, actual value will be set during advection
		}
		Real lt = lifetimeMin + I_ke / keMax * uniform_distr_0_1(random_generator) *(lifetimeMax - lifetimeMin);
		flifetimePerThread[tid].push_back(lt);
	}
}

/** Avect all foam particles according to 
* Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
*/
void advectFoamParticles()
{
	Vector3r g(0.0, -9.81, 0.0);
	if (splitTypes)
	{
		numParticlesOfType[ParticleType::Foam] = 0;
		numParticlesOfType[ParticleType::Spray] = 0;
		numParticlesOfType[ParticleType::Bubbles] = 0;
	}

	std::vector<std::vector<unsigned int>> neighbors;
	neighbors.reserve(100);

	#pragma omp parallel default(shared), private(neighbors)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) fx.size(); i++)
		{
			const Vector3r &xi = fx[i];

			// determine particle type
			neighborhoodSearch->find_neighbors(xi.data(), neighbors);
			const unsigned int numFluidNeighbors = (unsigned int) neighbors[0].size();

			unsigned char ftype = ParticleType::Foam;
			if (numFluidNeighbors < 6)
				ftype = ParticleType::Spray;
			else if (numFluidNeighbors > 20)
				ftype = ParticleType::Bubbles;

			// Fix particle type if classified as foam because of boundary.
			// If a particle is close to the boundary the number of fluid neighbors will be lower than for a particle inside the fluid.
			// A corrected number of fluid neighbors is estimated by adding the current number of fluid neighbors scaled by the volume
			// fraction of the cap of the sphere making up the neighborhood that lies outside of the boundary.
			if (ftype == ParticleType::Foam)
			{
				auto correctedFluidNeighbors = numFluidNeighbors;
				// Used to compute the volume of the cap of the neighborhood sphere that lies outside the boundary.
				auto unitSphereCapVolume = [&](const Real h)
				{
					return h * h * (M_PI / 3) * (3 - h);
				};
				for (auto j = 0u; j < 3; ++j)
				{
					Real d;
					d = std::abs(xi[j] - bbMin[j]);
					if (d < supportRadius)
						correctedFluidNeighbors += (decltype(correctedFluidNeighbors))(std::ceil(correctedFluidNeighbors * unitSphereCapVolume(1 - d / supportRadius)));
					d = std::abs(xi[j] - bbMax[j]);
					if (d < supportRadius)
						correctedFluidNeighbors += (decltype(correctedFluidNeighbors))(std::ceil(correctedFluidNeighbors * unitSphereCapVolume(1 - d / supportRadius)));
					if (correctedFluidNeighbors > 20)
					{
						ftype = ParticleType::Bubbles;
						break;
					}
				}
			}

			if (splitTypes)
			{
				particleType[i] = ftype;
				numParticlesOfType[ftype]++;
			}

			// handle bounding box
			if (((xi.array() < bbMin.array() || xi.array() > bbMax.array()).any()))
			{
				if (bbType == BbType::Kill)
				{
					flifetime[i] = 0;
				}
				else if (bbType == BbType::Lifesteal)
				{
					flifetime[i] -= static_cast<Real>(1000.0) * timeStepSize;
				}
				// clamp particle position to bounding box and reflect velocity on box wall
				else if (bbType == BbType::Clamp)
				{
					const Vector3r bbNormal = ((fx[i].array() < bbMin.array()).cast<Real>() - (fx[i].array() > bbMax.array()).cast<Real>()).matrix().normalized();
					const Vector3r vReflect = fv[i] - 2 * bbNormal.dot(fv[i]) * bbNormal;
					fv[i] = Real(0.25) * vReflect;
					fx[i] = (fx[i].array() < bbMin.array()).select(bbMin, fx[i]);
					fx[i] = (fx[i].array() > bbMax.array()).select(bbMax, fx[i]);
					// also lifesteal
					//flifetime[i] -= 1000.0 * timeStepSize;
				}
			}

			// advect particle
			if (ftype == ParticleType::Spray)
			{
				// spray
				fv[i] += timeStepSize * g;
				//fv[i] *= 0.99;
				fx[i] += timeStepSize * fv[i];

// 				if (numFluidNeighbors < 2)
// 					flifetime[i] -= 15.0*timeStepSize;
// 				else
// 					flifetime[i] -= 2.0*timeStepSize;
			}
			else if ((ftype == ParticleType::Foam) || (ftype == ParticleType::Bubbles))
			{
				Vector3r vf = Vector3r::Zero();
				Real sumK = 0.0;
				for (unsigned int j = 0; j < numFluidNeighbors; j++)
				{
					const unsigned int neighborIndex = neighbors[0][j];
					const Vector3r &xj = x0[neighborIndex];
					const Real K = CubicKernel::W(xi - xj);
					const Vector3r v = v0[neighborIndex]; 
					vf += v * K;
					sumK += K;
				}
				vf = (1.0 / sumK) * vf;

				if (ftype == ParticleType::Foam)
				{
					// foam
					//fv[i] = vf;
					fx[i] += timeStepSize * vf;
					flifetime[i] -= timeStepSize;
				}
				else if (ftype == ParticleType::Bubbles)
				{
					// bubbles
					fv[i] += timeStepSize * (-k_buoyancy * g + k_drag * invDt * (vf - fv[i]));
					fx[i] += timeStepSize * fv[i];

					//flifetime[i] -= 0.333*timeStepSize;
				}
			}
		}
	}
	if (splitTypes || splitGenerators)
		for (unsigned int i = 0; i < numTypes; i++)
			LOG_INFO << "#" << nameExtensions[i] << ": " << numParticlesOfType[i];
}

/** Remove foam particles which are at the end of their lifetime.
*/
void removeParticles()
{
	unsigned int removedParticles = 0;
	for (int i = 0; i < (int)fx.size(); i++)
	{
		if (flifetime[i] <= 0.0)
		{
			removedParticles++;
		}
		else
		{
			fx[i - removedParticles] = fx[i];
			fv[i - removedParticles] = fv[i];
			flifetime[i - removedParticles] = flifetime[i];
			particleType[i - removedParticles] = particleType[i];
		}
	}
	if (removedParticles > 0)
	{
		fx.resize(fx.size() - removedParticles);
		fv.resize(fv.size() - removedParticles);
		flifetime.resize(flifetime.size() - removedParticles);
		particleType.resize(particleType.size() - removedParticles);
	}
}

/** Read particle data from a partio file.
*/
bool readParticles(const std::string &fileName, std::vector<Vector3r> &positions, std::vector<Vector3r> &velocities, std::vector<Vector3r> &angularVelocities)
{
	if (!FileSystem::fileExists(fileName))
		return false;

	Partio::ParticlesDataMutable* data = Partio::read(fileName.c_str());
	if (!data)
		return false;

	unsigned int posIndex = 0xffffffff;
	unsigned int velIndex = 0xffffffff;
	unsigned int omegaIndex = 0xffffffff;

	for (int i = 0; i < data->numAttributes(); i++)
	{
		Partio::ParticleAttribute attr;
		data->attributeInfo(i, attr);
		if (attr.name == "position")
			posIndex = i;
		else if (attr.name == "velocity")
			velIndex = i;
		else if (attr.name == "angularVelocity")
			omegaIndex = i;
	}

	Partio::ParticleAttribute attr;

	if (posIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int)positions.size();
		positions.resize(fSize + data->numParticles());
		data->attributeInfo(posIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *pos = data->data<float>(attr, i);
			Vector3r x(pos[0], pos[1], pos[2]);
			positions[i + fSize] = x;
		}
	}

	if (velIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int)velocities.size();
		velocities.resize(fSize + data->numParticles());
		data->attributeInfo(velIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *vel = data->data<float>(attr, i);
			Vector3r v(vel[0], vel[1], vel[2]);
			velocities[i + fSize] = v;
		}
	}
	else
	{
		unsigned int fSize = (unsigned int)velocities.size();
		velocities.resize(fSize + data->numParticles());
		for (int i = 0; i < data->numParticles(); i++)
			velocities[i + fSize].setZero();
	}

	if (omegaIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int)angularVelocities.size();
		angularVelocities.resize(fSize + data->numParticles());
		data->attributeInfo(omegaIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *vel = data->data<float>(attr, i);
			Vector3r v(vel[0], vel[1], vel[2]);
			angularVelocities[i + fSize] = v;
		}
	}

	data->release();
	return true;
}
