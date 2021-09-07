#include "SPlisHSPlasH/Common.h"
#include "SamplingBase.h"
#include "Utilities/Timing.h"
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/TriangleMesh.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/Version.h"
#include "Utilities/FileSystem.h"


using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;

SamplingBase::SamplingBase()
{
	m_radius = static_cast<Real>(0.025);
	m_scale = Vector3r::Ones();
	m_diameter = m_radius * static_cast<Real>(2.0);
	m_resolutionSDF = Eigen::Matrix<unsigned int, 3, 1>(30, 30, 30);
	m_useRegion = false;
	m_output_format = 0;
}

SamplingBase::~SamplingBase()
{

}

void SamplingBase::generateSampling(const std::string& inputFile, const std::string& outputFile)
{
	loadObj(inputFile, m_mesh, m_scale);
	computeBoundingBox(m_mesh);

	bool md5 = false;
	string exePath = FileSystem::getProgramPath();
	string cachePath = exePath + "/Cache";
	if (m_useCache)
	{
		string md5Str = FileSystem::getFileMD5(inputFile);
		std::string mesh_file_name = FileSystem::getFileName(inputFile);
		const string scaleStr = "s" + StringTools::real2String(m_scale[0]) + "_" + StringTools::real2String(m_scale[1]) + "_" + StringTools::real2String(m_scale[2]);
		const string resStr = "r" + to_string(m_resolutionSDF[0]) + "_" + to_string(m_resolutionSDF[1]) + "_" + to_string(m_resolutionSDF[2]);
		const string invertStr = "i" + to_string((int)m_invert);

		string cacheFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_" + md5Str + "_" + scaleStr + "_" + resStr + "_" + invertStr  + ".csdf");

		// check MD5 if cache file is available
		bool foundCacheFile = FileSystem::fileExists(cacheFileName);

		if (foundCacheFile)
		{
			m_distanceField = std::make_shared<Discregrid::CubicLagrangeDiscreteGrid>(cacheFileName);
			LOG_INFO << "Loaded SDF cache file: " << cacheFileName;
		}
		else
		{
			if (FileSystem::makeDir(cachePath) == 0)
			{
				START_TIMING("generateSDF");
				generateSDF(m_mesh);
				STOP_TIMING_AVG_PRINT;

				m_distanceField->save(cacheFileName);
				LOG_INFO << "Saved SDF cache file: " << cacheFileName;
			}
			else
				LOG_ERR << "Failed to make directory: " << cachePath;
		}
	}
	else
	{
		START_TIMING("generateSDF");
		generateSDF(m_mesh);
		STOP_TIMING_AVG_PRINT;
	}

	m_output_format = 0;
	if (Utilities::StringTools::to_upper(FileSystem::getFileExt(outputFile)) == "VTK")
		m_output_format = 1;

	generateSamples();

	if (m_output_format == 0)
		PartioReaderWriter::writeParticles(outputFile, (unsigned int)m_x.size(), m_x.data(), NULL, 0.0);
	else
		writeParticlesVTK(outputFile, m_x);

	std::cout << "Generated particles: " << m_x.size() << std::endl;
}

void SamplingBase::writeParticlesVTK(const std::string& fileName, std::vector<Vector3r>& x)
{
	const unsigned int numParticles = (int)x.size();
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
			positions.emplace_back(x[i]);
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

void SamplingBase::loadObj(const std::string& filename, TriangleMesh& mesh, const Vector3r& scale)
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

void SamplingBase::sampleObject(const int mode)
{
	// sample object
	const unsigned int numberOfSamplePoints = (((unsigned int)((1.0f / m_diameter) * (m_bbmax[2] - m_bbmin[2]))) + 1) *
		(((unsigned int)((1.0f / m_diameter) * (m_bbmax[1] - m_bbmin[1]))) + 1) *
		(((unsigned int)((1.0f / m_diameter) * (m_bbmax[0] - m_bbmin[0]))) + 1);
	unsigned int currentSample = 0;
	Real currentPercent = 0.01;
	int counter_x = 0;
	int counter_y = 0;
	Real xshift = m_diameter;
	Real yshift = m_diameter;

	if (mode == 1)
		yshift = sqrt(static_cast<Real>(3.0)) * m_radius;
	else if (mode == 2)
	{
		xshift = sqrt(static_cast<Real>(3.0)) * m_radius;
		yshift = sqrt(static_cast<Real>(6.0)) * m_diameter / static_cast<Real>(3.0);
	}
	for (Real z = m_bbmin[2]; z <= m_bbmax[2]; z += m_diameter)
	{
		for (Real y = m_bbmin[1]; y <= m_bbmax[1]; y += yshift)
		{
			for (Real x = m_bbmin[0]; x <= m_bbmax[0]; x += xshift)
			{
				Vector3r particlePosition;
				if (mode == 1)
				{
					if (counter_y % 2 == 0)
						particlePosition = Vector3r(x, y + m_radius, z + m_radius);
					else
						particlePosition = Vector3r(x + m_radius, y + m_radius, z);
				}
				else if (mode == 2)
				{
					particlePosition = Vector3r(x, y + m_radius, z + m_radius);

					Vector3r shift_vec(0, 0, 0);
					if (counter_x % 2)
					{
						shift_vec[2] += m_diameter / (static_cast<Real>(2.0) * (counter_y % 2 ? -1 : 1));
					}
					if (counter_y % 2)
					{
						shift_vec[0] += xshift / static_cast<Real>(2.0);
						shift_vec[2] += m_diameter / static_cast<Real>(2.0);
					}
					particlePosition += shift_vec;
				}
				else
				{
					// Use center of voxel
					particlePosition = Vector3r(x + m_radius, y + m_radius, z + m_radius);
				}

				const double dist = distance(particlePosition, 0.0);
				if ((dist != std::numeric_limits<double>::max()) && (dist > 0.0))
				{
					m_x.push_back(particlePosition);
				}
				currentSample++;

				if ((Real)currentSample / (Real)numberOfSamplePoints > currentPercent)
				{
					std::cout << "\r" << "Generating samples " << std::setw(14) << currentPercent * 100.0 << "%";
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
}


void SamplingBase::computeBoundingBox(TriangleMesh& mesh)
{
	Vector3r* v = mesh.getVertices().data();

	// compute bounding box	 
	m_bbmin = v[0];
	m_bbmax = m_bbmin;
	for (unsigned int i = 1; i < mesh.numVertices(); ++i)
	{
		const Vector3r& p = v[i];
		for (unsigned int j = 0; j < 3; ++j) {
			m_bbmin[j] = std::min(m_bbmin[j], p[j]);
			m_bbmax[j] = std::max(m_bbmax[j], p[j]);
		}
	}

	if (m_useRegion)
	{
		for (unsigned int i = 0; i < 3; i++)
		{
			m_bbmin[i] = std::max(m_scale[i] * m_region.m_min[i], m_bbmin[i]);
			m_bbmax[i] = std::min(m_scale[i] * m_region.m_max[i], m_bbmax[i]);
		}
	}
}


// Determine distance of a point x to the surface of the mesh and corresponding surface normal and 
// next point on the surface.
double SamplingBase::distance(const Vector3r& x, const Real tolerance, Vector3r& normal, Vector3r& nextSurfacePoint)
{
	Eigen::Vector3d n;
	const double dist = m_distanceField->interpolate(0, x.template cast<double>(), &n);
	if (dist == std::numeric_limits<double>::max())
		return dist;
	n.normalize();
	normal = n.template cast<Real>();

	nextSurfacePoint = (x - dist * normal);

	return dist - tolerance;
}

// Determine distance of a point x to the surface of the mesh 
double SamplingBase::distance(const Vector3r& x, const Real tolerance)
{
	const double dist = m_distanceField->interpolate(0, x.template cast<double>());
	if (dist == std::numeric_limits<double>::max())
		return dist;
	return dist - tolerance;
}

void SamplingBase::generateSDF(SPH::TriangleMesh& mesh)
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
	domain.extend(m_bbmin.cast<double>());
	domain.extend(m_bbmax.cast<double>());
	domain.max() += Eigen::Vector3d::Ones() * m_radius * 4.0;
	domain.min() -= Eigen::Vector3d::Ones() * m_radius * 4.0;

	LOG_INFO << "Set SDF resolution: " << m_resolutionSDF[0] << ", " << m_resolutionSDF[1] << ", " << m_resolutionSDF[2];
	m_distanceField = std::make_shared<Discregrid::CubicLagrangeDiscreteGrid>(domain, std::array<unsigned int, 3>({ m_resolutionSDF[0], m_resolutionSDF[1], m_resolutionSDF[2] }));
	auto func = Discregrid::DiscreteGrid::ContinuousFunction{};

	// invert the distance field since the particles should stay inside
	if (!m_invert)
		func = [&md](Eigen::Vector3d const& xi) {return -md.signedDistanceCached(xi); };
	else
		func = [&md](Eigen::Vector3d const& xi) {return md.signedDistanceCached(xi); };

	LOG_INFO << "Generate SDF";
	m_distanceField->addFunction(func, false);
}

