#include "SPlisHSPlasH/Common.h"
#include <Eigen/Dense>
#include "extern/cxxopts/cxxopts.hpp"
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/TriangleMesh.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/Version.h"
#include "Utilities/FileSystem.h"
#include "SPlisHSPlasH/NeighborhoodSearch.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Utilities/SceneLoader.h"
#include "extern/happly/happly.h"
#include "extern/partio/src/lib/Partio.h"
#include "SPlisHSPlasH/Utilities/MeshImport.h"

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

void meshSkinning();
void computeMatrixL();
void computeF();
void deformMesh();
void performMeshSkinning();
bool readFrame(std::vector<Vector3r>& x, const unsigned int frame);
void writeCurrentFrame();
void init();
void precomputeValues();
void exportOBJ();
void exportVTK();
void exportPLY();

struct Comparator {
	Comparator(const Vector3r &xi, std::vector<Vector3r>* x) : m_xi(xi), m_x(x) {};
	bool operator()(unsigned int a, unsigned int b)
	{
		return ((*m_x)[a] - m_xi).squaredNorm() < ((*m_x)[b] - m_xi).squaredNorm();
	}

	Vector3r m_xi;
	std::vector<Vector3r> *m_x;
};


const Real density0 = 1000.0;
string input = "";
string output = "";
string meshFile = "";
string sceneFile = "";
bool overwriteExistingFiles = false;
Vector3r scale(1.0, 1.0, 1.0);
Vector3r translation(1.0, 1.0, 1.0);
Matrix3r rotation = Matrix3r::Identity();
unsigned int startFrame = 1;
unsigned int endFrame = 0xffffffff;
unsigned int currentFrame = 1;
string exePath;
bool useSplitting = false;
int output_format = 0; // 0: OBJ, 1: vtk, 2: ply
Real radius = static_cast<Real>(0.025);
Real diameter = radius * static_cast<Real>(2.0);
std::string partioPath = "";
Real supportRadius;
Real mass = 1.0;
TriangleMesh mesh;
std::vector<Vector3r> x0;
std::vector<Vector3r> x;
std::vector<int> id;
std::vector<Vector3r> meshX;
std::vector<Matrix3r> L;
std::vector<Matrix3r> F;
std::vector<Real> densities;
std::vector<Real> restVolumes;
std::vector<Real> shepard;
Real supportRadiusFactor = 6.0;
unsigned int maxNeighbors = 30;
Real W_zero;
Real(*kernelFct)(const Vector3r&);
Vector3r(*gradKernelFct)(const Vector3r& r);
NeighborhoodSearch* neighborhoodSearch;
std::vector<unsigned int> activeParticles;
std::vector<std::vector<unsigned int>> initialNeighbors;
std::vector<std::vector<unsigned int>> initialMeshNeighbors;
#ifdef USE_AVX
std::vector<Vector3f8, Eigen::aligned_allocator<Vector3f8>> precomp_V_gradW8;
std::vector<unsigned int> precomputed_indices8;
#else
std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> precomp_V_gradW;
std::vector<unsigned int> precomputed_indices;
#endif


unsigned int numberOfNeighbors(const unsigned int pointSetIndex, const unsigned int index)
{
	return static_cast<unsigned int>(neighborhoodSearch->point_set(0).n_neighbors(pointSetIndex, index));
}

unsigned int getNeighbor(const unsigned int pointSetIndex, const unsigned int index, const unsigned int k)
{
	return neighborhoodSearch->point_set(0).neighbor(pointSetIndex, index, k);
}

std::ostream& operator << (std::ostream& out, const Vector3r& r)
{
	out << r[0] << ", " << r[1] << ", " << r[2];
	return out;
}


// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;

	Utilities::logger.addSink(shared_ptr<Utilities::ConsoleSink>(new Utilities::ConsoleSink(Utilities::LogLevel::INFO)));
	exePath = FileSystem::getProgramPath();

	try
	{
		cxxopts::Options options(argv[0], "MeshSkinning");

		options.add_options()
			("i,input", "Input file ", cxxopts::value<std::string>())
			("o,output", "Output file", cxxopts::value<std::string>())
			("m,mesh", "Mesh file ", cxxopts::value<std::string>())
			("scene", "Scene file (all settings are imported from the scene file)", cxxopts::value<std::string>())
			("partioPath", "Path of the partio files (when using a scene file). If not set, it is assumed that the files are in the standard output path.", cxxopts::value<std::string>())
			("scale", "Scaling of input geometry (e.g. --scale 2,1,2)", cxxopts::value<std::vector<Real>>()->default_value("1,1,1"))
			("t,translation", "Translation of input geometry (e.g. --translation 2,1,2)", cxxopts::value<std::vector<Real>>()->default_value("0,0,0"))
			("axis", "Rotation axis of input geometry (e.g. --axis 1,0,0)", cxxopts::value<std::vector<Real>>()->default_value("1,0,0"))
			("angle", "Angle of input geometry (e.g. --angle 1)", cxxopts::value<Real>()->default_value("0.0"))
			("s,startframe", "Start frame", cxxopts::value<unsigned int>()->default_value("1"))
			("e,endframe", "End frame", cxxopts::value<unsigned int>())
			("r,radius", "Particle radius", cxxopts::value<Real>()->default_value("0.025"))			
			("supportRadiusFactor", "The support radius is defined as factor*particleRadius", cxxopts::value<Real>()->default_value("6.0"))
			("maxNeighbors", "The maximum number of neighbors that are used for the interpolation.", cxxopts::value<unsigned int>()->default_value("60"))
			("splitting", "Read a scene which used the object splitting export option.")
			("overwrite", "Overwrite existing files.")
			("h,help", "Print help")
			;

		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
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
			else if (Utilities::StringTools::to_upper(FileSystem::getFileExt(output)) == "PLY")
				output_format = 2;
		}
		if (result.count("mesh"))
		{
			meshFile = result["mesh"].as<std::string>();
			LOG_INFO << "mesh = " << meshFile;
		}
		if (result.count("scene"))
		{
			sceneFile = result["scene"].as<std::string>();
		}

		if (result.count("partioPath"))
		{
			partioPath = result["partioPath"].as<std::string>();
			if (FileSystem::isRelativePath(partioPath))
				partioPath = FileSystem::normalizePath(exePath + "/" + partioPath);
			LOG_INFO << "Partio path = " << partioPath;
		}

		if (result.count("splitting"))
			useSplitting = true;

		if (result.count("overwrite"))
			overwriteExistingFiles = true;

		if (sceneFile == "")
		{
			if (input == "")
			{
				LOG_ERR << "Input is missing!";
				exit(1);
			}
			if (output == "")
			{
				LOG_ERR << "Output is missing!";
				exit(1);
			}
			if (meshFile == "")
			{
				LOG_ERR << "Mesh file is missing!";
				exit(1);
			}

			if (result.count("scale"))
				scale = Vector3r(result["scale"].as<std::vector<Real>>().data()); 
			LOG_INFO << "Scale: " << scale;

			if (result.count("translation"))
				translation = Vector3r(result["translation"].as<std::vector<Real>>().data()); 
			LOG_INFO << "Translation: " << translation;

			Vector3r axis = Vector3r::Zero();
			Real angle = 0.0;
			rotation = Matrix3r::Identity();
			if (result.count("axis"))
				axis = Vector3r(result["axis"].as<std::vector<Real>>().data()); 
			if (result.count("angle"))
			{
				angle = result["angle"].as<Real>();
				axis.normalize();
				rotation = AngleAxisr(angle, axis);

				LOG_INFO << "Axis: " << axis;
				LOG_INFO << "Angle: " << angle;
			}

			if (result.count("radius"))
				radius = result["radius"].as<Real>();

			LOG_INFO << "Particle radius: " << radius;
		}

		if (result.count("startframe"))
			startFrame = result["startframe"].as<unsigned int>();

		if (result.count("endframe"))
			endFrame = result["endframe"].as<unsigned int>();

		if (result.count("supportRadiusFactor"))
			supportRadiusFactor = result["supportRadiusFactor"].as<Real>();

		if (result.count("maxNeighbors"))
			maxNeighbors = result["maxNeighbors"].as<unsigned int>();
	}
	catch (const cxxopts::exceptions::exception& e)
	{
		LOG_ERR << "error parsing options: " << e.what();
		exit(1);
	}

	performMeshSkinning();

	Timing::printAverageTimes();
	Timing::printTimeSums();

	return 0;
}


void init()
{
	diameter = static_cast<Real>(2.0) * radius;
	mass = static_cast<Real>(1.0) * diameter * diameter * diameter * density0;

	restVolumes.resize(x0.size());
	shepard.resize(mesh.numVertices());
	densities.resize(x0.size());
	L.resize(x0.size());
	F.resize(x0.size());

	supportRadius = supportRadiusFactor * radius;

	CubicKernel::setRadius(supportRadius);
#ifdef USE_AVX
	CubicKernel_AVX::setRadius(supportRadius);
#endif
	W_zero = CubicKernel::W_zero();
	kernelFct = CubicKernel::W;
	gradKernelFct = CubicKernel::gradW;

	// Init neighborhood search
	neighborhoodSearch = new NeighborhoodSearch(supportRadius, false);
	neighborhoodSearch->set_radius(supportRadius);
	neighborhoodSearch->add_point_set(&x0[0][0], x0.size(), true, true);

	// find initial neighbors
	START_TIMING("neighborhoodSearch");
	neighborhoodSearch->find_neighbors();
	STOP_TIMING_AVG;

	// store initial neighbors and init rest volumes
	int numParticles = (int) x0.size();
	initialNeighbors.resize(numParticles);
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < numParticles; i++)
		{
			// only neighbors in same phase will influence elasticity
			const unsigned int numNeighbors = numberOfNeighbors(0, i);
			initialNeighbors[i].resize(numNeighbors);
			for (unsigned int j = 0; j < numNeighbors; j++)
				initialNeighbors[i][j] = getNeighbor(0, i, j);

			std::sort(initialNeighbors[i].begin(), initialNeighbors[i].end(), Comparator(x0[i], &x0));

			if (initialNeighbors[i].size() > maxNeighbors)
				initialNeighbors[i].resize(maxNeighbors);

			// compute volume
			Real d = W_zero;
			const Vector3r& xi = x0[i];
			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(0, i, j);
				const Vector3r& xj = x0[neighborIndex];
				d += kernelFct(xi - xj);
			}
			restVolumes[i] = static_cast<Real>(1.0) / d;
		}
	}

	// compute the neighbors for each vertex in the mesh
	unsigned int numVertices = mesh.numVertices();
	initialMeshNeighbors.resize(numVertices);
	Vector3r* v = mesh.getVertices().data();

	std::vector<std::vector<unsigned int>> neighbors;
	neighbors.reserve(200);
	meshX.resize(numVertices);
	#pragma omp parallel default(shared), private(neighbors)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int) numVertices; i++)
		{
			const Vector3r& xi = v[i];
			// copy mesh x
			meshX[i] = xi;

			// determine particle type
			neighborhoodSearch->find_neighbors(xi.data(), neighbors);

			// only neighbors in same phase will influence elasticity
			const unsigned int numNeighbors = (unsigned int) neighbors[0].size();
			initialMeshNeighbors[i].resize(numNeighbors);
			for (unsigned int j = 0; j < numNeighbors; j++)
				initialMeshNeighbors[i][j] = neighbors[0][j];

			std::sort(initialMeshNeighbors[i].begin(), initialMeshNeighbors[i].end(), Comparator(meshX[i], &x0));

			if (initialMeshNeighbors[i].size() > maxNeighbors)
				initialMeshNeighbors[i].resize(maxNeighbors);
		}

		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)mesh.numVertices(); i++)
		{
			shepard[i] = 0.0;
			const Vector3r& xi = meshX[i];
			for (unsigned int j = 0; j < (unsigned int) initialMeshNeighbors[i].size(); j++)
			{
				const unsigned int neighborIndex = initialMeshNeighbors[i][j];
				const Vector3r& xj = x0[neighborIndex];
				shepard[i] += restVolumes[neighborIndex] * kernelFct(xi - xj);
			}

			shepard[i] = static_cast<Real>(1.0) / shepard[i];
		}
	}

	int counter = 0;
	std::vector<bool> active;
	active.resize(numParticles, false);
	unsigned int minNeighbors = 0xffffffff;
	unsigned int maxNeighbors = 0;
	for (unsigned int i = 0; i < numVertices; i++)
	{
		const unsigned int numNeighbors = (unsigned int)initialMeshNeighbors[i].size();
		minNeighbors = std::min(minNeighbors, numNeighbors);
		maxNeighbors = std::max(maxNeighbors, numNeighbors);

		if (numNeighbors == 0)
			LOG_INFO << i;

		for (unsigned int j = 0; j < numNeighbors; j++)
		{
			unsigned int neighbor = initialMeshNeighbors[i][j];
			if (!active[neighbor])
			{
				active[neighbor] = true;
				counter++;
			}			
		}
	}
	activeParticles.clear();
	activeParticles.reserve(counter);
	for (int i = 0; i < numParticles; i++)
	{
		if (active[i])
			activeParticles.push_back(i);
	}
	LOG_INFO << "Active particles: " << activeParticles.size();
	LOG_INFO << "Min. # neighbors: " << minNeighbors;
	LOG_INFO << "Max. # neighbors: " << maxNeighbors;

	START_TIMING("computeMatrixL");
	computeMatrixL();
	STOP_TIMING_AVG;

	START_TIMING("precomputeValues");
	precomputeValues();
	STOP_TIMING_AVG;

	delete neighborhoodSearch;
	neighborhoodSearch = nullptr;
}

void performMeshSkinning(const std::string &id, const int objId)
{
	const std::string outputPath = FileSystem::normalizePath(exePath + "/output/" + FileSystem::getFileName(sceneFile));
	if (partioPath == "")
		partioPath = FileSystem::normalizePath(outputPath + "/partio");

	std::string fileName;
	if (useSplitting)
	{
		fileName = FileSystem::normalizePath(partioPath + "/ParticleData_" + id + "_" + std::to_string(objId) + "_" + std::to_string(startFrame) + ".bgeo");
		input = FileSystem::normalizePath(partioPath + "/ParticleData_" + id + "_" + std::to_string(objId) + "_#.bgeo");
		output = FileSystem::normalizePath(partioPath + "/../meshes/ParticleData_" + id + "_" + std::to_string(objId) + "_#.ply");
	}
	else
	{
		fileName = FileSystem::normalizePath(partioPath + "/ParticleData_" + id + "_" + std::to_string(startFrame) + ".bgeo");
		input = FileSystem::normalizePath(partioPath + "/ParticleData_" + id + "_#.bgeo");
		output = FileSystem::normalizePath(partioPath + "/../meshes/ParticleData_" + id + "_#.ply");
	}

	if (FileSystem::fileExists(fileName))
	{
		output_format = 2;

		if (!FileSystem::fileExists(meshFile))
		{
			LOG_ERR << "Visualization mesh file defined in scene file is missing: " << meshFile;
			return;
		}

		LOG_INFO << "Input: " << input;
		LOG_INFO << "Output: " << output;
		LOG_INFO << "Mesh file: " << meshFile;
		LOG_INFO << "Particle radius: " << radius;
		LOG_INFO << "Scale: " << scale;
		LOG_INFO << "Translation: " << translation;

		MeshImport::importMesh(meshFile, mesh, translation, rotation, scale);

		// read reference configuration
		readFrame(x0, 1);
		LOG_INFO << "# particles: " << x0.size();

		init();

		meshSkinning();
	}
	else
	{
		LOG_ERR << "Input file not found: " << fileName;
	}
}

void performMeshSkinning()
{
	if (sceneFile != "")
	{
		SceneLoader* sceneLoader = new SceneLoader();
		Utilities::SceneLoader::Scene scene;
		sceneLoader->readScene(sceneFile.c_str(), scene);

		const std::string outputPath = FileSystem::normalizePath(exePath + "/output/" + FileSystem::getFileName(sceneFile));
		if (partioPath == "")
			partioPath = FileSystem::normalizePath(outputPath + "/partio");
		const std::string scene_path = FileSystem::getFilePath(sceneFile);

		int objId = 0;
		for (size_t i = 0; i < scene.fluidBlocks.size(); i++)
		{
			translation = 0.5 * (scene.fluidBlocks[i]->boxMax + scene.fluidBlocks[i]->boxMin);
			rotation = Matrix3r::Identity();
			scale = scene.fluidBlocks[i]->boxMax - scene.fluidBlocks[i]->boxMin;
			radius = scene.particleRadius;
			std::string id = scene.fluidBlocks[i]->id;
			meshFile = FileSystem::normalizePath(FileSystem::getFilePath(sceneFile) + "/" + scene.fluidBlocks[i]->visMeshFile);

			performMeshSkinning(id, objId);
			objId++;
		}


		for (size_t i = 0; i < scene.fluidModels.size(); i++)
		{
			translation = scene.fluidModels[i]->translation;
			rotation = AngleAxisr(scene.fluidModels[i]->angle, scene.fluidModels[i]->axis).toRotationMatrix();
			scale = scene.fluidModels[i]->scale;
			radius = scene.particleRadius;
			std::string id = scene.fluidModels[i]->id;
			meshFile = FileSystem::normalizePath(FileSystem::getFilePath(sceneFile) + "/" + scene.fluidModels[i]->visMeshFile);
			LOG_INFO << FileSystem::getFilePath(sceneFile) << ", " << sceneFile;
			LOG_INFO << meshFile;

			performMeshSkinning(id, objId);
			objId++;
		}	

		delete sceneLoader;
	}
	else
	{
		MeshImport::importMesh(meshFile, mesh, translation, rotation, scale);

		// read reference configuration
		readFrame(x0, 1);
		LOG_INFO << "# particles: " << x0.size();

		init();

		meshSkinning();
	}
}

std::string zeroPadding(const unsigned int number, const unsigned int length) 
{
	ostringstream out;
	out << std::internal << std::setfill('0') << std::setw(length) << number;
	return out.str();
}

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

bool readFrame(std::vector<Vector3r>& x, const unsigned int frame)
{
	std::string fileName = convertFileName(input, frame);

	x.clear();
	id.clear();

	// check if file exists
	if (!FileSystem::fileExists(fileName))
		return false;

	// read partio file
	Partio::ParticlesDataMutable* data = Partio::read(fileName.c_str());
	if (!data)
		return false;

	// find position and id attribute
	unsigned int posIndex = 0xffffffff;
	unsigned int idIndex = 0xffffffff;
	for (int i = 0; i < data->numAttributes(); i++)
	{
		Partio::ParticleAttribute attr;
		data->attributeInfo(i, attr);
		if (attr.name == "position")
			posIndex = i;
		if (attr.name == "id")
			idIndex = i;
	}

	// read IDs
	Partio::ParticleAttribute attr;
	if (idIndex != 0xffffffff)
	{
		id.resize(data->numParticles());
		data->attributeInfo(idIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const int* id_ = data->data<int>(attr, i);
			id[i] = id_[0];
		}
	}

	// compute the first id since for multiple objects, the ids might not start at 0
	int firstId = INT_MAX;
	for (int i = 0; i < id.size(); i++)
		firstId = std::min(id[i], firstId);


	// read positions
	if (posIndex != 0xffffffff)
	{
		x.resize(data->numParticles());
		data->attributeInfo(posIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float* pos = data->data<float>(attr, i);

			// use id as index if possible to consider z-sorting in the particle data
			if (idIndex != 0xffffffff)
				x[id[i]- firstId] = Vector3r(pos[0], pos[1], pos[2]);
			else 
				x[i] = Vector3r(pos[0], pos[1], pos[2]);
		}
	}

	data->release();
	return true;
}


void computeMatrixL()
{
	const int numParticles = (int) x0.size();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			const Vector3r &xi0 = x0[i];
			Matrix3r L_;
			L_.setZero();

			const size_t numNeighbors = initialNeighbors[i].size();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numNeighbors; j++)
			{
				const unsigned int neighborIndex = initialNeighbors[i][j];

				const Vector3r &xj0 = x0[neighborIndex];
				const Vector3r xj_xi_0 = xj0 - xi0;
				const Vector3r gradW = gradKernelFct(xj_xi_0);

				// minus because gradW(xij0) == -gradW(xji0)
				L_ -= restVolumes[neighborIndex] * gradW * xj_xi_0.transpose();
			}

			//// add 1 to z-component. otherwise we get a singular matrix in 2D
			//if (sim->is2DSimulation())
			//	L(2, 2) = 1.0;

			bool invertible = false;
			L_.computeInverseWithCheck(L[i], invertible, 1e-9);
			if (!invertible)
			{
				//MathFunctions::pseudoInverse(L, m_L[i]);
				L[i] = Matrix3r::Identity();
			}
		}
	}
}

#ifdef USE_AVX

void computeF()
{
	const int numParticles = (int)activeParticles.size();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int index = activeParticles[i];
			const Vector3r &xi = x[index];
			const Vector3r& xi0 = x0[index];

			const Vector3f8 xi_avx(xi);
			const Vector3f8 xi0_avx(xi0);

			const unsigned int numNeighbors = (unsigned int) initialNeighbors[index].size();

			Matrix3f8 F_avx;
			F_avx.setZero();
			for (unsigned int j = 0; j < numNeighbors; j += 8)
			{
				const unsigned int count = std::min(numNeighbors - j, 8u);

				const Vector3f8 xj_avx = convertVec_zero(&initialNeighbors[index][j], &x[0], count);
				const Vector3f8 xj_xi = xj_avx - xi_avx;
				const Vector3f8& V_gradW = precomp_V_gradW8[precomputed_indices8[i] + j / 8];

				const Scalarf8 restVolume_j_avx = convert_zero(&initialNeighbors[index][j], &restVolumes[0], count);
				const Vector3f8 xj0_avx = convertVec_zero(&initialNeighbors[index][j], &x0[0], count);
				//const Vector3f8& V_gradW = Matrix3f8(L[index]) * CubicKernel_AVX::gradW(xi0_avx - xj0_avx) * restVolume_j_avx;
				Matrix3f8 dyad;
				dyadicProduct(xj_xi, V_gradW, dyad);
				F_avx += dyad;
			}
			F[index] = F_avx.reduce();
		}
	}
}


void deformMesh()
{
	const int numVertices = (int)mesh.numVertices();
	auto& meshX0 = mesh.getVertices();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numVertices; i++)
		{
			Vector3r &mesh_xi = meshX[i];
			const Vector3r& mesh_xi0 = meshX0[i];
			const Vector3f8 xi0_avx(mesh_xi0);

			const unsigned int numNeighbors = (unsigned int)initialMeshNeighbors[i].size();

			Vector3f8 xNew;
			xNew.setZero();
			for (unsigned int j = 0; j < numNeighbors; j += 8)
			{
				const unsigned int count = std::min(numNeighbors - j, 8u);

				const Vector3f8 xj_avx = convertVec_zero(&initialMeshNeighbors[i][j], &x[0], count);
				const Vector3f8 xj0_avx = convertVec_zero(&initialMeshNeighbors[i][j], &x0[0], count);
				const Matrix3f8 Fj_avx = convertMat_zero(&initialMeshNeighbors[i][j], &F[0], count);

				const Vector3f8 xi_xj_0 = xi0_avx - xj0_avx;
				const Scalarf8 restVolume_j_avx = convert_zero(&initialMeshNeighbors[i][j], &restVolumes[0], count);
				const Scalarf8 W = CubicKernel_AVX::W(xi_xj_0);

				xNew += (Fj_avx * xi_xj_0 + xj_avx) * restVolume_j_avx * W;
			}
			mesh_xi = xNew.reduce();
			mesh_xi *= shepard[i];
		}
	}
}

#else

void computeF()
{
	const int numParticles = (int)activeParticles.size();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int index = activeParticles[i];
			const Vector3r &xi = x[index];
			const Vector3r& xi0 = x0[index];
			F[index].setZero();

			const size_t numNeighbors = initialNeighbors[index].size();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numNeighbors; j++)
			{
				const unsigned int neighborIndex = initialNeighbors[index][j];

				const Vector3r &xj = x[neighborIndex];
				const Vector3r &xj0 = x0[neighborIndex];
				const Vector3r xj_xi = xj - xi;
				const Vector3r xi_xj_0 = xi0 - xj0;
				//const Vector3r correctedKernel = L[index] * gradKernelFct(xi_xj_0);
				const Vector3r V_gradW = precomp_V_gradW[precomputed_indices[i] + j];
				F[index] += xj_xi * V_gradW.transpose();
			}
		}
	}
}

void deformMesh()
{
	const int numVertices = (int)mesh.numVertices();
	auto& meshX0 = mesh.getVertices();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numVertices; i++)
		{
			Vector3r &mesh_xi = meshX[i];
			const Vector3r& mesh_xi0 = meshX0[i];

			const size_t numNeighbors = initialMeshNeighbors[i].size();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			Vector3r deltaX;
			deltaX.setZero();
			for (unsigned int j = 0; j < numNeighbors; j++)
			{
				const unsigned int neighborIndex = initialMeshNeighbors[i][j];

				const Vector3r &xj = x[neighborIndex];
				const Vector3r &xj0 = x0[neighborIndex];
				const Vector3r xi_xj_0 = mesh_xi0 - xj0;
				deltaX += shepard[i] * restVolumes[neighborIndex] * (F[neighborIndex] * xi_xj_0 + xj) * kernelFct(xi_xj_0);
			}
			mesh_xi = deltaX;
		}
	}
}
#endif




void precomputeValues()
{
	const int numParticles = (int)activeParticles.size();
	unsigned int numVertices = mesh.numVertices();
	auto& meshX0 = mesh.getVertices();

#ifdef USE_AVX
	precomputed_indices8.clear();
	precomp_V_gradW8.clear();
	precomputed_indices8.resize(numParticles);
#else
	precomputed_indices.clear();
	precomp_V_gradW.clear();
	precomputed_indices.resize(numParticles);
#endif


	unsigned int sumNeighborParticles = 0;
	unsigned int sumNeighborParticles8 = 0;
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int index = activeParticles[i];
		const size_t numNeighbors = initialNeighbors[index].size();

#ifdef USE_AVX
		precomputed_indices8[i] = sumNeighborParticles8;

		// steps of 8 values due to avx
		sumNeighborParticles8 += (unsigned int) numNeighbors / 8u;
		if (numNeighbors % 8 != 0)
			sumNeighborParticles8++;
#else
		precomputed_indices[i] = sumNeighborParticles;
		sumNeighborParticles += (int) numNeighbors;
#endif
	}

#ifdef USE_AVX
	precomp_V_gradW8.resize(sumNeighborParticles8);
#else
	precomp_V_gradW.resize(sumNeighborParticles);
#endif

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int index = activeParticles[i];
			const Vector3r& xi0 = x0[index];
			const unsigned int numNeighbors = (unsigned int)initialNeighbors[index].size();

#ifdef USE_AVX
			const Vector3f8 xi0_avx(xi0);
			unsigned int base8 = precomputed_indices8[i];
			unsigned int idx = 0;
			Matrix3f8 L_avx(L[index]);
			for (unsigned int j = 0; j < numNeighbors; j += 8)
			{
				const unsigned int count = std::min(numNeighbors - j, 8u);
				const Scalarf8 restVolume_j_avx = convert_zero(&initialNeighbors[index][j], &restVolumes[0], count);
				const Vector3f8 xj0_avx = convertVec_zero(&initialNeighbors[index][j], &x0[0], count);

				const Vector3f8 gradW = CubicKernel_AVX::gradW(xi0_avx - xj0_avx);
				precomp_V_gradW8[base8 + idx] = (L_avx * gradW) * restVolume_j_avx;
				idx++;
			}
#else
			unsigned int base = precomputed_indices[i];

			for (unsigned int j = 0; j < numNeighbors; j++)
			{
				const unsigned int neighborIndex = initialNeighbors[index][j];
				const Vector3r& xj0 = x0[neighborIndex];
				Vector3r gradW = gradKernelFct(xi0 - xj0);
				precomp_V_gradW[base + j] = restVolumes[neighborIndex] * L[index] * gradW;
			}
#endif
		}
	}
}

void meshSkinning()
{
	currentFrame = startFrame;
	std::string fileName = convertFileName(output, currentFrame);
	if (overwriteExistingFiles || !FileSystem::fileExists(fileName))
		writeCurrentFrame();
	else 
		LOG_INFO << "Skipping file: " << fileName;

	bool chk = true;
	while (chk)
	{
		std::string fileName = convertFileName(output, currentFrame);
		if (!overwriteExistingFiles && FileSystem::fileExists(fileName))
		{
			LOG_INFO << "Skipping file: " << fileName;
			currentFrame++;
			if (currentFrame > endFrame)
				break;
			continue;
		}

		LOG_INFO << "Reading frame: " << currentFrame;
		chk = readFrame(x, currentFrame);
		if (!chk)
			break;

		START_TIMING("computeF");
		computeF();
		STOP_TIMING_AVG;

		START_TIMING("deformMesh");
		deformMesh();
		STOP_TIMING_AVG;

		writeCurrentFrame();

		currentFrame++;
		if (currentFrame > endFrame)
			break;
	}
}


void exportOBJ()
{
	std::string fileName = convertFileName(output, currentFrame);
	LOG_INFO << "Writing: " << fileName;

	FileSystem::makeDir(FileSystem::getFilePath(FileSystem::normalizePath(output)));

	// Open the file
	std::ofstream outfile(fileName);
	if (!outfile)
	{
		LOG_WARN << "Cannot open a file to save OBJ mesh.";
		return;
	}

	// Header
	outfile << "# Created by SPlisHSPlasH version " << SPLISHSPLASH_VERSION << "\n";
	outfile << "g default\n";

	const std::vector<Vector3r>& vertices = meshX;
	const std::vector<unsigned int>& faces = mesh.getFaces();
	int n_vertices = (int)vertices.size();
	int n_triangles = (int)faces.size() / 3;

	// Vertices
	{
		for (int j = 0u; j < n_vertices; j++)
		{
			Vector3r x = vertices[j];
			outfile << "v " << x[0] << " " << x[1] << " " << x[2] << "\n";
		}
	}

	// faces
	{
		for (int j = 0; j < n_triangles; j++)
		{
			outfile << "f " << faces[3 * j + 0] + 1 << " " << faces[3 * j + 1] + 1 << " " << faces[3 * j + 2] + 1 << "\n";
		}
	}
	outfile.close();
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

void exportVTK()
{
	std::string fileName = convertFileName(output, currentFrame);
	LOG_INFO << "Writing: " << fileName;

	FileSystem::makeDir(FileSystem::getFilePath(FileSystem::normalizePath(output)));

	// Open the file
	std::ofstream outfile(fileName, std::ios::binary);
	if (!outfile)
	{
		LOG_WARN << "Cannot open a file to save VTK mesh.";
		return;
	}

#ifdef USE_DOUBLE
	const char* real_str = " double\n";
#else 
	const char* real_str = " float\n";
#endif

	// Header
	outfile << "# vtk DataFile Version 4.2\n";
	outfile << "Created by SPlisHSPlasH version " << SPLISHSPLASH_VERSION << "\n";
	outfile << "BINARY\n";
	outfile << "DATASET UNSTRUCTURED_GRID\n";

	const std::vector<Vector3r>& vertices = meshX;
	const std::vector<unsigned int>& faces = mesh.getFaces();
	int n_vertices = (int)vertices.size();
	int n_triangles = (int)faces.size() / 3;

	// Vertices
	{
		std::vector<Vector3r> positions;
		positions.reserve(n_vertices);
		for (int j = 0u; j < n_vertices; j++)
		{
			Vector3r x = vertices[j];
			swapByteOrder(&x[0]);
			swapByteOrder(&x[1]);
			swapByteOrder(&x[2]);
			positions.emplace_back(x);
		}
		// export to vtk
		outfile << "POINTS " << n_vertices << real_str;
		outfile.write(reinterpret_cast<char*>(positions[0].data()), 3 * n_vertices * sizeof(Real));
		outfile << "\n";
	}

	// Connectivity
	{
		std::vector<int> connectivity_to_write;
		connectivity_to_write.reserve(4 * n_triangles);
		for (int tri_i = 0; tri_i < n_triangles; tri_i++)
		{
			int val = 3;
			swapByteOrder(&val);
			connectivity_to_write.push_back(val);
			val = faces[3 * tri_i + 0];
			swapByteOrder(&val);
			connectivity_to_write.push_back(val);
			val = faces[3 * tri_i + 1];
			swapByteOrder(&val);
			connectivity_to_write.push_back(val);
			val = faces[3 * tri_i + 2];
			swapByteOrder(&val);
			connectivity_to_write.push_back(val);
		}
		// export to vtk
		outfile << "CELLS " << n_triangles << " " << 4 * n_triangles << "\n";
		outfile.write(reinterpret_cast<char*>(&connectivity_to_write[0]), connectivity_to_write.size() * sizeof(int));
		outfile << "\n";
	}

	// Cell types
	{
		outfile << "CELL_TYPES " << n_triangles << "\n";
		int cell_type_swapped = 5;
		swapByteOrder(&cell_type_swapped);
		std::vector<int> cell_type_arr(n_triangles, cell_type_swapped);
		outfile.write(reinterpret_cast<char*>(&cell_type_arr[0]), cell_type_arr.size() * sizeof(int));
		outfile << "\n";
	}

	outfile.close();
}

void exportPLY()
{
	std::string fileName = convertFileName(output, currentFrame);
	LOG_INFO << "Writing: " << fileName;

	FileSystem::makeDir(FileSystem::getFilePath(FileSystem::normalizePath(output)));


	// Suppose these hold your data
	std::vector<std::array<double, 3>> meshVertexPositions;
	std::vector<std::vector<size_t>> meshFaceIndices;

	const std::vector<Vector3r>& vertices = meshX;
	const std::vector<unsigned int>& faces = mesh.getFaces();
	int n_vertices = (int)vertices.size();
	int n_triangles = (int)faces.size() / 3;

	meshVertexPositions.resize(n_vertices);
	for (int i = 0; i < n_vertices; i++)
	{
		meshVertexPositions[i] = { vertices[i][0], vertices[i][1], vertices[i][2] };
	}

	meshFaceIndices.resize(n_triangles);
	for (int i = 0; i < n_triangles; i++)
	{
		meshFaceIndices[i].resize(3);
		meshFaceIndices[i] = { faces[3*i], faces[3 * i+1], faces[3 * i+2] };
	}

	// Create an empty object
	happly::PLYData plyOut;

	// Add mesh data (elements are created automatically)
	plyOut.addVertexPositions(meshVertexPositions);
	plyOut.addFaceIndices(meshFaceIndices);


	// Write the object to file
	plyOut.write(fileName, happly::DataFormat::Binary);
}

void writeCurrentFrame()
{
	if (output_format == 1)
		exportVTK();
	else if (output_format == 2)
		exportPLY();
	else 
		exportOBJ();
}