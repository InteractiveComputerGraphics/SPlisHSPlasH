#include "ParticleExporter_Partio.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"
#include "extern/partio/src/lib/Partio.h"

using namespace SPH;
using namespace Utilities;

ParticleExporter_Partio::ParticleExporter_Partio(SimulatorBase *base) :
	ExporterBase(base)
{
}

ParticleExporter_Partio::~ParticleExporter_Partio(void)
{
}

void ParticleExporter_Partio::init(const std::string& outputPath)
{
	m_exportPath = FileSystem::normalizePath(outputPath + "/partio");
}

void ParticleExporter_Partio::step(const unsigned int frame)
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
		writeParticlesPartio(exportFileName + ".bgeo", model);
	}
}

void ParticleExporter_Partio::reset()
{
}

void ParticleExporter_Partio::setActive(const bool active)
{
	ExporterBase::setActive(active);
	if (m_active)
		FileSystem::makeDirs(m_exportPath);
}


void ParticleExporter_Partio::writeParticlesPartio(const std::string& fileName, FluidModel* model)
{
	const bool async = m_base->getValue<bool>(SimulatorBase::ASYNC_EXPORT);
	if (async)
	{
		if (m_handle.valid())
			m_handle.wait();
	}

	m_particleData = Partio::create();
	Partio::ParticlesDataMutable& particleData = *m_particleData;

	Partio::ParticleAttribute posAttr = particleData.addAttribute("position", Partio::VECTOR, 3);
	Partio::ParticleAttribute idAttr = particleData.addAttribute("id", Partio::INT, 1);

	// add attributes
	std::vector<std::string> attributes;
	StringTools::tokenize(m_base->getValue<std::string>(SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES), attributes, ";");

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

	m_particleFile = fileName;
	if (async)
		m_handle = std::async(std::launch::async, [&] { Partio::write(m_particleFile.c_str(), particleData, true); particleData.release(); });
	else
	{
		Partio::write(m_particleFile.c_str(), particleData, true);
		particleData.release();
	}
}