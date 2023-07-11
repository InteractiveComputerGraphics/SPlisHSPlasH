#include "RigidBodyExporter_OBJ.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"
#include <Utilities/Version.h>

using namespace SPH;
using namespace Utilities;

RigidBodyExporter_OBJ::RigidBodyExporter_OBJ(SimulatorBase* base) :
	ExporterBase(base)
{
	m_isFirstFrame = true;
}

RigidBodyExporter_OBJ::~RigidBodyExporter_OBJ(void)
{
}

void RigidBodyExporter_OBJ::init(const std::string& outputPath)
{
	// define output path for the data
	m_exportPath = FileSystem::normalizePath(outputPath + "/obj");
}

void RigidBodyExporter_OBJ::step(const unsigned int frame)
{
	// check if the exporter is active
	if (!m_active)
		return;

	writeRigidBodies(frame);
}

void RigidBodyExporter_OBJ::reset()
{
	m_isFirstFrame = true;
}

void RigidBodyExporter_OBJ::setActive(const bool active)
{
	ExporterBase::setActive(active);
	// create output folder
	if (m_active)
		FileSystem::makeDirs(m_exportPath);
}


void RigidBodyExporter_OBJ::writeRigidBodies(const unsigned int frame)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nBoundaryModels = sim->numberOfBoundaryModels();

	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel* bm = sim->getBoundaryModel(i);

		// If we have a static body, write the data only for the first frame.
		// Otherwise each frame is exported.
		if (m_isFirstFrame || bm->getRigidBodyObject()->isDynamic() || bm->getRigidBodyObject()->isAnimated())
		{
			std::string fileName = "rb_data_";
			fileName = fileName + std::to_string(i) + "_" + std::to_string(frame) + ".obj";
			std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);

			// Open the file
			std::ofstream outfile(exportFileName);
			if (!outfile)
			{
				LOG_WARN << "Cannot open a file to save OBJ mesh.";
				return;
			}

			// Header
			outfile << "# Created by SPlisHSPlasH version " << SPLISHSPLASH_VERSION << "\n";
			outfile << "g default\n";

			const std::vector<Vector3r>& vertices = bm->getRigidBodyObject()->getVertices();
			const std::vector<unsigned int>& faces = bm->getRigidBodyObject()->getFaces();
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
	}

	m_isFirstFrame = false;
}