#include "RigidBodyParticleExporter_VTK.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"
#include <regex>
using namespace SPH;
using namespace Utilities;

RigidBodyParticleExporter_VTK::RigidBodyParticleExporter_VTK(SimulatorBase* base) : ExporterBase(base) {
	m_isFirstFrame = true;
}

RigidBodyParticleExporter_VTK::~RigidBodyParticleExporter_VTK(void) {}

void RigidBodyParticleExporter_VTK::init(const std::string& outputPath) {
	m_exportPath = FileSystem::normalizePath(outputPath + "/vtk");
}

void RigidBodyParticleExporter_VTK::step(const unsigned int frame) {
	if (!m_active) {
		return;
	}
	writeRigidBodies(frame);
}

void RigidBodyParticleExporter_VTK::reset() {
	m_isFirstFrame = true;
}

void RigidBodyParticleExporter_VTK::setActive(const bool active) {
	ExporterBase::setActive(active);
	if (m_active) {
		FileSystem::makeDirs(m_exportPath);
	}
}

void RigidBodyParticleExporter_VTK::writeRigidBodies(const unsigned int frame) {
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nBoundaryModels = sim->numberOfBoundaryModels();

	// check if we have a static model
	bool isStatic = true;
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++) {
		BoundaryModel* bm = sim->getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic() || bm->getRigidBodyObject()->isAnimated()) {
			isStatic = false;
			break;
		}
	}

#ifdef USE_DOUBLE
	const char* real_str = " double\n";
#else 
	const char* real_str = " float\n";
#endif

	if (m_isFirstFrame || !isStatic) {
		for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++) {
			std::string fileName = "rb_particle_data_";
			fileName = fileName + std::to_string(i) + "_" + std::to_string(frame) + ".vtk";
			std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);

			// Open the file
			std::ofstream outfile(exportFileName, std::ios::binary);
			if (!outfile) {
				LOG_WARN << "Cannot open a file to save VTK mesh.";
				return;
			}

			// Header
			outfile << "# vtk DataFile Version 4.2\n";
			outfile << "SPlisHSPlasH rigid body particle data\n";
			outfile << "BINARY\n";
			outfile << "DATASET UNSTRUCTURED_GRID\n";

			// add attributes
			m_attributes.clear();
			StringTools::tokenize(m_base->getValue<std::string>(SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES), m_attributes, ";");

			//////////////////////////////////////////////////////////////////////////
			// positions and ids exported anyways
			m_attributes.erase(
				std::remove_if(m_attributes.begin(), m_attributes.end(), [](const std::string& s) { return (s == "position" || s == "id"); }),
				m_attributes.end());

			BoundaryModel_Akinci2012* model = dynamic_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(i));
			if (model != nullptr) {
				Simulation* sim = Simulation::getCurrent();
				const unsigned int nFluids = sim->numberOfFluidModels();

				const unsigned int numBoundaryParticles = model->numberOfParticles();

				std::vector<Vector3r> positions;
				positions.reserve(numBoundaryParticles);
				std::vector<Eigen::Vector2i> cells;
				cells.reserve(numBoundaryParticles);
				std::vector<int> cellTypes;
				cellTypes.reserve(numBoundaryParticles);
				std::vector<unsigned int> attrData;
				attrData.reserve(numBoundaryParticles);

				unsigned int counter = 0;
				for (unsigned int i = 0; i < numBoundaryParticles; i++) {


					//////////////////////////////////////////////////////////////////////////
					// export position attribute as POINTS
					{
						positions.push_back(model->getPosition(i));
						// swap endianess
						for (unsigned int c = 0; c < 3; c++)
							swapByteOrder(&positions[positions.size() - 1][c]);
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
					//{
					//	unsigned int id = model->getParticleId(i);
					//	swapByteOrder(&id);
					//	attrData.push_back(id);
					//}
				}

				//createParticleFile(fileName, model);

				if (outfile.good()) {
					// export to vtk
					const unsigned int nPoints = (unsigned int)positions.size();
					outfile << "POINTS " << nPoints << real_str;
					outfile.write(reinterpret_cast<char*>(positions[0].data()), 3 * nPoints * sizeof(Real));
					outfile << "\n";

					// particles are cells with one element and the index of the particle
					outfile << "CELLS " << nPoints << " " << 2 * nPoints << "\n";
					outfile.write(reinterpret_cast<char*>(cells[0].data()), 2 * nPoints * sizeof(unsigned int));
					outfile << "\n";

					outfile << "CELL_TYPES " << nPoints << "\n";
					outfile.write(reinterpret_cast<char*>(cellTypes.data()), nPoints * sizeof(int));
					outfile << "\n";

					outfile << "POINT_DATA " << nPoints << "\n";
					// write IDs
					outfile << "SCALARS id unsigned_int 1\n";
					outfile << "LOOKUP_TABLE id_table\n";
					outfile.write(reinterpret_cast<char*>(attrData.data()), nPoints * sizeof(unsigned int));
					outfile << "\n";



					//////////////////////////////////////////////////////////////////////////
					// per point fields (all attributes except for positions)
					const auto numFields = m_attributes.size();
					outfile << "FIELD FieldData " << std::to_string(numFields) << "\n";

					// iterate over attributes.
					for (const std::string& a : m_attributes) {
						const FieldDescription& field = model->getField(a);

						std::string attrNameVTK;
						std::regex_replace(std::back_inserter(attrNameVTK), a.begin(), a.end(), std::regex("\\s+"), "_");

						if (field.type == FieldType::Scalar) {
							// write header information
							outfile << attrNameVTK << " 1 " << nPoints << real_str;

							// copy data
							std::vector<Real> attrData;
							attrData.reserve(nPoints);
							for (unsigned int i = 0u; i < numBoundaryParticles; i++) {
								Real val = *((Real*)field.getFct(i));
								swapByteOrder(&val);		// swap endianess
								attrData.emplace_back(val);
							}

							// export to vtk
							outfile.write(reinterpret_cast<char*>(attrData.data()), nPoints * sizeof(Real));
						} else if (field.type == FieldType::Vector3) {
							// write header information
							outfile << attrNameVTK << " 3 " << nPoints << real_str;

							// copy from partio data
							std::vector<Vector3r> attrData;
							attrData.reserve(nPoints);
							for (unsigned int i = 0u; i < numBoundaryParticles; i++) {
								Vector3r val((Real*)field.getFct(i));
								for (unsigned int c = 0; c < 3; c++)
									swapByteOrder(&val[c]);		// swap endianess
								attrData.emplace_back(val);
							}
							// export to vtk
							outfile.write(reinterpret_cast<char*>(attrData[0].data()), 3 * nPoints * sizeof(Real));
						} else if (field.type == FieldType::UInt) {
							// write header information
							outfile << attrNameVTK << " 1 " << nPoints << " unsigned_int\n";

							// copy data
							std::vector<unsigned int> attrData;
							attrData.reserve(nPoints);
							for (unsigned int i = 0u; i < numBoundaryParticles; i++) {
								unsigned int val = *((unsigned int*)field.getFct(i));
								swapByteOrder(&val);		// swap endianess
								attrData.emplace_back(val);
							}
							// export to vtk
							outfile.write(reinterpret_cast<char*>(attrData.data()), nPoints * sizeof(unsigned int));
						}
						// TODO support other field types
						else {
							LOG_WARN << "Skipping attribute " << a << ", because it is of unsupported type\n";
							continue;
						}
						// end of block
						outfile << "\n";

					}
					outfile.close();
				}

			}

		}

		m_isFirstFrame = false;

	}
}