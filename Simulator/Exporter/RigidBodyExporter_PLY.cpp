#include "RigidBodyExporter_PLY.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"
#include "extern/happly/happly.h"
#include <Utilities/Version.h>

using namespace SPH;
using namespace Utilities;

RigidBodyExporter_PLY::RigidBodyExporter_PLY(SimulatorBase* base) :
	ExporterBase(base)
{
	m_isFirstFrame = true;
}

RigidBodyExporter_PLY::~RigidBodyExporter_PLY(void)
{
}

void RigidBodyExporter_PLY::init(const std::string& outputPath)
{
	// define output path for the data
	m_exportPath = FileSystem::normalizePath(outputPath + "/ply");
}

void RigidBodyExporter_PLY::step(const unsigned int frame)
{
	// check if the exporter is active
	if (!m_active)
		return;

	writeRigidBodies(frame);
}

void RigidBodyExporter_PLY::reset()
{
	m_isFirstFrame = true;
}

void RigidBodyExporter_PLY::setActive(const bool active)
{
	ExporterBase::setActive(active);
	// create output folder
	if (m_active)
		FileSystem::makeDirs(m_exportPath);
}


void RigidBodyExporter_PLY::writeRigidBodies(const unsigned int frame)
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

	// If we have a static model, write the data only for the first frame.
	// Otherwise each frame is exported.
	if (m_isFirstFrame || !isStatic)
	{
		for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
		{
			std::string fileName = "rb_data_";
			fileName = fileName + std::to_string(i) + "_" + std::to_string(frame) + ".ply";
			std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);

			BoundaryModel* bm = sim->getBoundaryModel(i);
			const std::vector<Vector3r>& vertices = bm->getRigidBodyObject()->getVertices();
			const std::vector<unsigned int>& faces = bm->getRigidBodyObject()->getFaces();
			int n_vertices = (int)vertices.size();
			int n_triangles = (int)faces.size() / 3;

			// Suppose these hold your data
			std::vector<std::array<double, 3>> meshVertexPositions;
			std::vector<std::vector<size_t>> meshFaceIndices;
			
			// vertices
			meshVertexPositions.resize(n_vertices);
			for (int j = 0u; j < n_vertices; j++)
				meshVertexPositions[j] = { vertices[j][0], vertices[j][1], vertices[j][2] };

			// faces
			meshFaceIndices.resize(n_triangles);
			for (int j = 0; j < n_triangles; j++)
			{
				meshFaceIndices[j].resize(3);
				meshFaceIndices[j] = { faces[3 * j], faces[3 * j + 1], faces[3 * j + 2] };
			}

			// Create an empty object
			happly::PLYData plyOut;

			// Add mesh data (elements are created automatically)
			plyOut.addVertexPositions(meshVertexPositions);
			plyOut.addFaceIndices(meshFaceIndices);


			// Write the object to file
			plyOut.write(exportFileName, happly::DataFormat::Binary);
		}
	}

	m_isFirstFrame = false;
}