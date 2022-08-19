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
	m_outfile = nullptr;
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
		if (!m_base->getValue<bool>(SimulatorBase::EXPORT_OBJECT_SPLITTING))
		{
			fileName = fileName + "_" + model->getId() + "_" + std::to_string(frame);
			std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);
			writeParticles(exportFileName + ".vtk", model);
		}
		else
		{
			// object splitting
			for (auto j = 0u; j < m_base->getLastObjectId(); j++)
			{
				std::string fileName2 = fileName + "_" + model->getId() + "_" + std::to_string(j) + "_" + std::to_string(frame);
				std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName2);
				writeParticles(exportFileName + ".vtk", model, j);
			}
		}
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

void ParticleExporter_VTK::createParticleFile(const std::string& fileName, FluidModel* model)
{
	// Open the file
	m_outfile = new std::ofstream(fileName, std::ios::binary);
	if (!m_outfile->is_open())
	{
		LOG_WARN << "Cannot open a file to save VTK particles.";
		return;
	}

	*m_outfile << "# vtk DataFile Version 4.1\n";
	*m_outfile << "SPlisHSPlasH particle data\n"; // title of the data set, (any string up to 256 characters+\n)
	*m_outfile << "BINARY\n";
	*m_outfile << "DATASET UNSTRUCTURED_GRID\n";

	// add attributes
	m_attributes.clear();
	StringTools::tokenize(m_base->getValue<std::string>(SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES), m_attributes, ";");


	//////////////////////////////////////////////////////////////////////////
	// positions and ids exported anyways
	m_attributes.erase(
		std::remove_if(m_attributes.begin(), m_attributes.end(), [](const std::string& s) { return (s == "position" || s == "id"); }),
		m_attributes.end());
}


void ParticleExporter_VTK::writeParticles(const std::string& fileName, FluidModel* model, const unsigned int objId)
{
#ifdef USE_DOUBLE
	const char* real_str = " double\n";
#else 
	const char* real_str = " float\n";
#endif

	const unsigned int numParticles = model->numActiveParticles();

	std::vector<Vector3r> positions;
	positions.reserve(numParticles);
	std::vector<Eigen::Vector2i> cells;
	cells.reserve(numParticles);
	std::vector<int> cellTypes;
	cellTypes.reserve(numParticles);
	std::vector<unsigned int> attrData;
	attrData.reserve(numParticles);

	unsigned int counter = 0;
	for (unsigned int i = 0; i < numParticles; i++)
	{
		if ((objId != 0xffffffff) && (model->getObjectId(i) != objId))
			continue;

		//////////////////////////////////////////////////////////////////////////
		// export position attribute as POINTS
		{
			positions.push_back(model->getPosition(i));
			// swap endianess
			for (unsigned int c = 0; c < 3; c++)
				swapByteOrder(&positions[positions.size()-1][c]);
		}

		//////////////////////////////////////////////////////////////////////////
		// export particle IDs as CELLS
		{
			unsigned int nodes_per_cell_swapped = 1;
			swapByteOrder(&nodes_per_cell_swapped);
			unsigned int idSwapped = counter++;
			swapByteOrder(&idSwapped);
			cells.push_back({ nodes_per_cell_swapped, idSwapped });
		}
		//////////////////////////////////////////////////////////////////////////
		// export cell types
		{
			// the type of a particle cell is always 1
			int cellTypeSwapped = 1;
			swapByteOrder(&cellTypeSwapped);
			cellTypes.push_back(cellTypeSwapped);
		}

		//////////////////////////////////////////////////////////////////////////
		// write additional attributes as per-particle data
		{
			unsigned int id = model->getParticleId(i);
			swapByteOrder(&id);
			attrData.push_back(id);
		}
	}

	createParticleFile(fileName, model);

	if (m_outfile != nullptr)
	{
		// export to vtk
		const unsigned int nPoints = (unsigned int) positions.size();
		*m_outfile << "POINTS " << nPoints << real_str;
		m_outfile->write(reinterpret_cast<char*>(positions[0].data()), 3 * nPoints * sizeof(Real));
		*m_outfile << "\n";

		// particles are cells with one element and the index of the particle
		*m_outfile << "CELLS " << nPoints << " " << 2 * nPoints << "\n";
		m_outfile->write(reinterpret_cast<char*>(cells[0].data()), 2 * nPoints * sizeof(unsigned int));
		*m_outfile << "\n";

		*m_outfile << "CELL_TYPES " << nPoints << "\n";
		m_outfile->write(reinterpret_cast<char*>(cellTypes.data()), nPoints * sizeof(int));
		*m_outfile << "\n";

		*m_outfile << "POINT_DATA " << nPoints << "\n";
		// write IDs
		*m_outfile << "SCALARS id unsigned_int 1\n";
		*m_outfile << "LOOKUP_TABLE id_table\n";
		m_outfile->write(reinterpret_cast<char*>(attrData.data()), nPoints * sizeof(unsigned int));
		*m_outfile << "\n";



		//////////////////////////////////////////////////////////////////////////
		// per point fields (all attributes except for positions)
		const auto numFields = m_attributes.size();
		*m_outfile << "FIELD FieldData " << std::to_string(numFields) << "\n";

		// iterate over attributes
		for (const std::string& a : m_attributes)
		{
			const FieldDescription& field = model->getField(a);

			std::string attrNameVTK;
			std::regex_replace(std::back_inserter(attrNameVTK), a.begin(), a.end(), std::regex("\\s+"), "_");

			if (field.type == FieldType::Scalar)
			{
				// write header information
				*m_outfile << attrNameVTK << " 1 " << nPoints << real_str;

				// copy data
				std::vector<Real> attrData;
				attrData.reserve(nPoints);
				for (unsigned int i = 0u; i < numParticles; i++)
				{
					if ((objId != 0xffffffff) && (model->getObjectId(i) != objId))
						continue;
					Real val = *((Real*)field.getFct(i));
					swapByteOrder(&val);		// swap endianess
					attrData.emplace_back(val);
				}

				// export to vtk
				m_outfile->write(reinterpret_cast<char*>(attrData.data()), nPoints * sizeof(Real));
			}
			else if (field.type == FieldType::Vector3)
			{
				// write header information
				*m_outfile << attrNameVTK << " 3 " << nPoints << real_str;

				// copy from partio data
				std::vector<Vector3r> attrData;
				attrData.reserve(nPoints);
				for (unsigned int i = 0u; i < numParticles; i++)
				{
					if ((objId != 0xffffffff) && (model->getObjectId(i) != objId))
						continue;
					Vector3r val((Real*)field.getFct(i));
					for (unsigned int c = 0; c < 3; c++)
						swapByteOrder(&val[c]);		// swap endianess
					attrData.emplace_back(val);
				}
				// export to vtk
				m_outfile->write(reinterpret_cast<char*>(attrData[0].data()), 3 * nPoints * sizeof(Real));
			}
			else if (field.type == FieldType::UInt)
			{
				// write header information
				*m_outfile << attrNameVTK << " 1 " << nPoints << " unsigned_int\n";

				// copy data
				std::vector<unsigned int> attrData;
				attrData.reserve(nPoints);
				for (unsigned int i = 0u; i < numParticles; i++)
				{
					if ((objId != 0xffffffff) && (model->getObjectId(i) != objId))
						continue;
					unsigned int val = *((unsigned int*)field.getFct(i));
					swapByteOrder(&val);		// swap endianess
					attrData.emplace_back(val);
				}
				// export to vtk
				m_outfile->write(reinterpret_cast<char*>(attrData.data()), nPoints * sizeof(unsigned int));
			}
			// TODO support other field types
			else
			{
				LOG_WARN << "Skipping attribute " << a << ", because it is of unsupported type\n";
				continue;
			}
			// end of block
			*m_outfile << "\n";

		}
        m_outfile->close();
        delete m_outfile;
        m_outfile = nullptr;
	}
}