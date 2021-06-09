#include "RigidBodyExporter_VTK.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;
using namespace Utilities;

RigidBodyExporter_VTK::RigidBodyExporter_VTK(SimulatorBase* base) :
	ExporterBase(base)
{
	m_isFirstFrame = true;
}

RigidBodyExporter_VTK::~RigidBodyExporter_VTK(void)
{
}

void RigidBodyExporter_VTK::init(const std::string& outputPath)
{
	m_exportPath = FileSystem::normalizePath(outputPath + "/vtk");
}

void RigidBodyExporter_VTK::step(const unsigned int frame)
{
	if (!m_active)
		return;

	writeRigidBodies(frame);
}

void RigidBodyExporter_VTK::reset()
{
	m_isFirstFrame = true;
}

void RigidBodyExporter_VTK::setActive(const bool active)
{
	ExporterBase::setActive(active);
	if (m_active)
		FileSystem::makeDirs(m_exportPath);
}


void RigidBodyExporter_VTK::writeRigidBodies(const unsigned int frame)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nBoundaryModels = sim->numberOfBoundaryModels();

	// check if we have a static model
	bool isStatic = true;
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel* bm = sim->getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic() || bm->getRigidBodyObject()->isAnimated())
		{
			isStatic = false;
			break;
		}
	}

#ifdef USE_DOUBLE
	const char* real_str = " double\n";
#else 
	const char* real_str = " float\n";
#endif

	if (m_isFirstFrame || !isStatic)
	{
		for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
		{
			std::string fileName = "rb_data_";
			fileName = fileName + std::to_string(i) + "_" + std::to_string(frame) + ".vtk";
			std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);

			// Open the file
			std::ofstream outfile(exportFileName, std::ios::binary);
			if (!outfile)
			{
				LOG_WARN << "Cannot open a file to save VTK mesh.";
				return;
			}

			// Header
			outfile << "# vtk DataFile Version 4.2\n";
			outfile << "SPlisHSPlasH mesh data\n";
			outfile << "BINARY\n";
			outfile << "DATASET UNSTRUCTURED_GRID\n";

			BoundaryModel* bm = sim->getBoundaryModel(i);
			const std::vector<Vector3r>& vertices = bm->getRigidBodyObject()->getVertices();
			const std::vector<unsigned int>& faces = bm->getRigidBodyObject()->getFaces();
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
	}

	m_isFirstFrame = false;
}