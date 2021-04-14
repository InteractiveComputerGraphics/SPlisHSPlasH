#include "ParticleExporter_VTK.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"
#include <regex>

using namespace SPH;
using namespace Utilities;

ParticleExporter_VTK::ParticleExporter_VTK(SimulatorBase *base) :
	ExporterBase(base)
{
}

ParticleExporter_VTK::~ParticleExporter_VTK(void)
{
}

void ParticleExporter_VTK::init(const std::string& outputPath)
{
	m_exportPath = FileSystem::normalizePath(outputPath + "/vtk");
}

void ParticleExporter_VTK::step(const unsigned int frame)
{
	if (!m_active)
		return;

	Simulation* sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel* model = sim->getFluidModel(i);
		std::string fileName = "ParticleData";
		fileName = fileName + "_" + model->getId() + "_" + std::to_string(frame);

		std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);
		writeParticles(exportFileName + ".vtk", model);
	}
}

void ParticleExporter_VTK::reset()
{
}

void ParticleExporter_VTK::setActive(const bool active)
{
	ExporterBase::setActive(active);
	if (m_active)
		FileSystem::makeDirs(m_exportPath);
}


void ParticleExporter_VTK::writeParticles(const std::string& fileName, FluidModel* model)
{
	const unsigned int numParticles = model->numActiveParticles();
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

	// add attributes
	std::vector<std::string> attributes;
	StringTools::tokenize(m_base->getValue<std::string>(SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES), attributes, ";");

	//////////////////////////////////////////////////////////////////////////
	// positions and ids exported anyways
	attributes.erase(
		std::remove_if(attributes.begin(), attributes.end(), [](const std::string& s) { return (s == "position" || s == "id"); }),
		attributes.end());

	//////////////////////////////////////////////////////////////////////////
	// export position attribute as POINTS
	{
		std::vector<Vector3r> positions;
		positions.reserve(numParticles);
		for (unsigned int i = 0u; i < numParticles; i++)
			positions.emplace_back(model->getPosition(i));
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
			unsigned int idSwapped = model->getParticleId(i);
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
			attrData.emplace_back(model->getParticleId(i));
		// swap endianess
		for (unsigned int i = 0; i < numParticles; i++)
			swapByteOrder(&attrData[i]);
		// export to vtk
		outfile.write(reinterpret_cast<char*>(attrData.data()), numParticles * sizeof(unsigned int));
		outfile << "\n";
	}

	//////////////////////////////////////////////////////////////////////////
	// per point fields (all attributes except for positions)
	const auto numFields = attributes.size();
	outfile << "FIELD FieldData " << std::to_string(numFields) << "\n";

	// iterate over attributes
	for (const std::string& a : attributes)
	{
		const FieldDescription& field = model->getField(a);

		std::string attrNameVTK;
		std::regex_replace(std::back_inserter(attrNameVTK), a.begin(), a.end(), std::regex("\\s+"), "_");

		if (field.type == FieldType::Scalar)
		{
			// write header information
			outfile << attrNameVTK << " 1 " << numParticles << real_str;

			// copy data
			std::vector<Real> attrData;
			attrData.reserve(numParticles);
			for (unsigned int i = 0u; i < numParticles; i++)
				attrData.emplace_back(*((Real*)field.getFct(i)));
			// swap endianess
			for (unsigned int i = 0; i < numParticles; i++)
				swapByteOrder(&attrData[i]);
			// export to vtk
			outfile.write(reinterpret_cast<char*>(attrData.data()), numParticles * sizeof(Real));
		}
		else if (field.type == FieldType::Vector3)
		{
			// write header information
			outfile << attrNameVTK << " 3 " << numParticles << real_str;

			// copy from partio data
			std::vector<Vector3r> attrData;
			attrData.reserve(numParticles);
			for (unsigned int i = 0u; i < numParticles; i++)
				attrData.emplace_back((Real*)field.getFct(i));
			// swap endianess
			for (unsigned int i = 0; i < numParticles; i++)
				for (unsigned int c = 0; c < 3; c++)
					swapByteOrder(&attrData[i][c]);
			// export to vtk
			outfile.write(reinterpret_cast<char*>(attrData[0].data()), 3 * numParticles * sizeof(Real));
		}
		else if (field.type == FieldType::Matrix3)
		{
			// write header information
			outfile << attrNameVTK << " 9 " << numParticles << real_str;

			// copy from partio data
			std::vector<Matrix3r> attrData;
			attrData.reserve(numParticles);
			for (unsigned int i = 0u; i < numParticles; i++)
			{
				Eigen::Map<Matrix3r> m((Real*)field.getFct(i));
				attrData.emplace_back(m);
			}
			// swap endianess
			for (unsigned int i = 0; i < numParticles; i++)
				for (unsigned int j = 0; j < 3; j++)
					for (unsigned int k = 0; k < 3; k++)
						swapByteOrder(&attrData[i](j, k));
			// export to vtk
			outfile.write(reinterpret_cast<char*>(attrData[0].data()), 9 * numParticles * sizeof(Real));
		}
		else if (field.type == FieldType::UInt)
		{
			// write header information
			outfile << attrNameVTK << " 1 " << numParticles << " unsigned_int\n";

			// copy data
			std::vector<unsigned int> attrData;
			attrData.reserve(numParticles);
			for (unsigned int i = 0u; i < numParticles; i++)
				attrData.emplace_back(*((unsigned int*)field.getFct(i)));
			// swap endianess
			for (unsigned int i = 0; i < numParticles; i++)
				swapByteOrder(&attrData[i]);
			// export to vtk
			outfile.write(reinterpret_cast<char*>(attrData.data()), numParticles * sizeof(unsigned int));
		}
		// TODO support other field types
		else
		{
			LOG_WARN << "Skipping attribute " << a << ", because it is of unsupported type\n";
			continue;
		}
		// end of block
		outfile << "\n";
	}
	outfile.close();
}