#include "Exporter.h"

#include <regex>

#include "extern/partio/src/lib/Partio.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/FileSystem.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Logger.h"

using namespace SPH;
using namespace Utilities;

int Exporter::PARTIO_EXPORT = -1;
int Exporter::VTK_EXPORT = -1;
int Exporter::RB_EXPORT = -1;
int Exporter::RB_VTK_EXPORT = -1;

bool Exporter::enablePartioExport = false;
bool Exporter::enableVTKExport = false;
bool Exporter::enableRigidBodyExport = false;
bool Exporter::enableRigidBodyVTKExport = false;

int Exporter::DATA_EXPORT_FPS = -1;
int Exporter::PARTICLE_EXPORT_ATTRIBUTES = -1;

Real Exporter::framesPerSecond = 25;
std::string Exporter::particleAttributes = "velocity";

Exporter::Exporter() {
	m_isFirstFrameVTK = true;
}

Exporter::~Exporter() {
	
}

void Exporter::rigidBodyExport(std::string scenePath, std::string temporalIdentifier, bool isFirstFrame, SceneLoader::Scene scene)
{
	std::string exportPath = FileSystem::normalizePath(m_outputPath + "/rigid_bodies");
	std::string exportPathVTK = FileSystem::normalizePath(m_outputPath + "/vtk");
	if (enableRigidBodyExport)
	{
		FileSystem::makeDirs(exportPath);
		writeRigidBodiesBIN(exportPath, scenePath, temporalIdentifier, isFirstFrame, scene);
	}
	if (enableRigidBodyVTKExport)
	{
		FileSystem::makeDirs(exportPathVTK);
		writeRigidBodiesVTK(exportPathVTK, temporalIdentifier);
	}
}

void Exporter::particleExport(std::string exportName, std::string temporalIdentifier, std::string folder, bool partio, bool vtk)
{
	std::string partioExportPath = FileSystem::normalizePath(m_outputPath + "/" + folder + "/partio");
	std::string vtkExportPath = FileSystem::normalizePath(m_outputPath + "/" + folder + "/vtk");
	if (partio)
		FileSystem::makeDirs(partioExportPath);
	if (vtk)
		FileSystem::makeDirs(vtkExportPath);

	Simulation *sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		std::string fileName;
		if (exportName.empty())
			fileName = model->getId() + "_" + temporalIdentifier;
		else
			fileName = exportName + "_" + model->getId() + "_" + temporalIdentifier;

		if (partio)
		{
			std::string exportFileName = FileSystem::normalizePath(partioExportPath + "/" + fileName);
			writeParticlesPartio(exportFileName + ".bgeo", model);
		}
		if (vtk)
		{
			std::string exportFileName = FileSystem::normalizePath(vtkExportPath + "/" + fileName);
			writeParticlesVTK(exportFileName + ".vtk", model);
		}
	}
}


void Exporter::writeParticlesPartio(const std::string& fileName, FluidModel* model)
{
	Partio::ParticlesDataMutable& particleData = *Partio::create();
	Partio::ParticleAttribute posAttr = particleData.addAttribute("position", Partio::VECTOR, 3);
	Partio::ParticleAttribute idAttr = particleData.addAttribute("id", Partio::INT, 1);

	// add attributes
	std::vector<std::string> attributes;
	StringTools::tokenize(particleAttributes, attributes, ";");

	std::map<unsigned int, int> attrMap;
	std::map<unsigned int, Partio::ParticleAttribute> partioAttrMap;
	for (unsigned int i = 0; i < attributes.size(); i++)
	{
		// position is exported anyway
		if (attributes[i] == "position")
		{
			attrMap[i] = -1;
			continue;
		}

		bool found = false;
		for (unsigned int j = 0; j < model->numberOfFields(); j++)
		{
			const FieldDescription& field = model->getField(j);
			if (field.name == attributes[i])
			{
				found = true;
				if (field.type == Scalar)
				{
					attrMap[i] = j;
					partioAttrMap[i] = particleData.addAttribute(attributes[i].c_str(), Partio::FLOAT, 1);
				}
				else if (field.type == UInt)
				{
					attrMap[i] = j;
					partioAttrMap[i] = particleData.addAttribute(attributes[i].c_str(), Partio::INT, 1);
				}
				else if (field.type == Vector3)
				{
					attrMap[i] = j;
					partioAttrMap[i] = particleData.addAttribute(attributes[i].c_str(), Partio::VECTOR, 3);
				}
				else
				{
					attrMap[i] = -1;
					LOG_WARN << "Only scalar and vector fields are currently supported by the partio exporter.";
				}
				break;
			}
		}
		if (!found)
		{
			attrMap[i] = -1;
			LOG_WARN << "Unknown field cannot be exported in partio file: " << attributes[i];
		}
	}

	const unsigned int numParticles = model->numActiveParticles();

	for (unsigned int i = 0; i < numParticles; i++)
	{
		Partio::ParticleIndex index = particleData.addParticle();
		float* pos = particleData.dataWrite<float>(posAttr, index);
		int* id = particleData.dataWrite<int>(idAttr, index);

		const Vector3r& x = model->getPosition(i);
		pos[0] = (float)x[0];
		pos[1] = (float)x[1];
		pos[2] = (float)x[2];

		id[0] = model->getParticleId(i);

		for (unsigned int j = 0; j < attributes.size(); j++)
		{
			const int fieldIndex = attrMap[j];
			if (fieldIndex != -1)
			{
				const FieldDescription& field = model->getField(fieldIndex);
				if (field.type == FieldType::Scalar)
				{
					float* val = particleData.dataWrite<float>(partioAttrMap[j], index);
					*val = (float)*((Real*)field.getFct(i));
				}
				else if (field.type == FieldType::UInt)
				{
					int* val = particleData.dataWrite<int>(partioAttrMap[j], index);
					*val = (int)*((unsigned int*)field.getFct(i));
				}
				else if (field.type == FieldType::Vector3)
				{
					float* val = particleData.dataWrite<float>(partioAttrMap[j], index);
					Eigen::Map<Vector3r> vec((Real*)field.getFct(i));
					val[0] = (float)vec[0];
					val[1] = (float)vec[1];
					val[2] = (float)vec[2];
				}
			}
		}
	}

	Partio::write(fileName.c_str(), particleData, true);
	particleData.release();
}

void Exporter::writeParticlesVTK(const std::string& fileName, FluidModel* model)
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
	StringTools::tokenize(particleAttributes, attributes, ";");

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

void Exporter::writeRigidBodiesBIN(const std::string& exportPath, std::string scene_path, std::string temporalIdentifier, bool isFirstFrame, SceneLoader::Scene scene)
{
	std::string fileName = "rb_data_";
	fileName = fileName + temporalIdentifier + ".bin";
	std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

	Simulation* sim = Simulation::getCurrent();
	const unsigned int nBoundaryModels = sim->numberOfBoundaryModels();

	// check if we have a static model
	bool isStatic = true;
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel* bm = sim->getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic())
		{
			isStatic = false;
			break;
		}
	}

	BinaryFileWriter binWriter;
	if (isFirstFrame || !isStatic)
		binWriter.openFile(exportFileName.c_str());

	if (isFirstFrame)
	{
		binWriter.write(nBoundaryModels);

		for (unsigned int i = 0; i < scene.boundaryModels.size(); i++)
		{
			std::string meshFileName = scene.boundaryModels[i]->meshFile;
			if (FileSystem::isRelativePath(meshFileName))
				meshFileName = FileSystem::normalizePath(scene_path + "/" + meshFileName);

			const std::string fileName = FileSystem::getFileNameWithExt(meshFileName);
			binWriter.write(fileName);
			Eigen::Vector3f s = scene.boundaryModels[i]->scale.template cast<float>();
			binWriter.writeMatrix(s);
			std::string targetFilePath = exportPath + "/" + fileName;
			if (!FileSystem::fileExists(targetFilePath))
			{
				FileSystem::copyFile(meshFileName, targetFilePath);
			}
			binWriter.write((char)scene.boundaryModels[i]->isWall);
			binWriter.writeMatrix(scene.boundaryModels[i]->color);
		}
	}

	if (isFirstFrame || !isStatic)
	{
		for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
		{
			BoundaryModel* bm = sim->getBoundaryModel(i);
			const Vector3r& x = bm->getRigidBodyObject()->getWorldSpacePosition();
			const Eigen::Vector3f x_f = x.template cast<float>();
			binWriter.writeMatrix(x_f);

			const Matrix3r& R = bm->getRigidBodyObject()->getWorldSpaceRotation();
			//const Eigen::Matrix3f RT = R.transpose().template cast<float>();
			binWriter.writeMatrix(R);
		}
		binWriter.closeFile();
	}
}

void Exporter::writeRigidBodiesVTK(const std::string& exportPath, std::string temporalIdentifier)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nBoundaryModels = sim->numberOfBoundaryModels();

	// check if we have a static model
	bool isStatic = true;
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel* bm = sim->getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic())
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

	if (m_isFirstFrameVTK || !isStatic)
	{
		for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
		{
			std::string fileName = "rb_data_";
			fileName = fileName + std::to_string(i) + "_" + temporalIdentifier + ".vtk";
			std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

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

	m_isFirstFrameVTK = false;
}