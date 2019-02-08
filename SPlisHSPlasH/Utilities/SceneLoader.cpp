#include "SceneLoader.h"
#include <iostream>
#include "extern/json/json.hpp"
#include <fstream>
#include "Utilities/FileSystem.h"
#include "Utilities/Logger.h"

using namespace Utilities;
using namespace GenParam;
using namespace std;


void SceneLoader::readScene(const char *fileName, Scene &scene)
{
	LOG_INFO << "Load scene file: " << fileName;

	std::ifstream input_file(fileName);
	if (!input_file.is_open())
	{
		LOG_ERR << "Cannot open file!";
		return;
	}
	try
	{
		m_jsonData << input_file;
	}
	catch (const std::exception& e)
	{
		LOG_ERR << e.what();
		exit(1);
	}	

	std::string base_path = FileSystem::getFilePath(fileName);

	//////////////////////////////////////////////////////////////////////////
	// read configuration 
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find("Configuration") != m_jsonData.end())
	{
		nlohmann::json config = m_jsonData["Configuration"];

		scene.timeStepSize = 0.001;
		readValue(config["timeStepSize"], scene.timeStepSize);

		scene.particleRadius = 0.025;
		readValue(config["particleRadius"], scene.particleRadius);

		scene.sim2D = false;
		readValue(config["sim2D"], scene.sim2D);
	}

	//////////////////////////////////////////////////////////////////////////
	// read boundary models
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find("RigidBodies") != m_jsonData.end())
	{
		nlohmann::json boundaryModels = m_jsonData["RigidBodies"];
		for (auto& boundaryModel : boundaryModels)
		{
			std::string particleFile = "";
			std::string meshFile = "";
			const bool bMesh = readValue<std::string>(boundaryModel["geometryFile"], meshFile);
			const bool bSamples = readValue<std::string>(boundaryModel["particleFile"], particleFile);

			if (bMesh || bSamples)
			{
				BoundaryData *data = new BoundaryData();
				data->meshFile = meshFile;
				data->samplesFile = particleFile;

				// translation
				data->translation = Vector3r::Zero();
				readVector(boundaryModel["translation"], data->translation);

				// rotation axis
				Vector3r axis = Vector3r::Zero();
				Real angle = 0.0;
				data->rotation = Matrix3r::Identity();
				if (readVector(boundaryModel["rotationAxis"], axis) &&
					readValue<Real>(boundaryModel["rotationAngle"], angle))
					data->rotation = AngleAxisr(angle, axis);

				// scale
				data->scale = Vector3r::Ones();
				readVector(boundaryModel["scale"], data->scale);

				data->dynamic = false;
				readValue<bool>(boundaryModel["isDynamic"], data->dynamic);

				data->isWall = false;
				readValue<bool>(boundaryModel["isWall"], data->isWall);

				data->color = Eigen::Vector4f(1.0f, 0.0f, 0.0f, 0.0f);
				readVector(boundaryModel["color"], data->color);

				scene.boundaryModels.push_back(data);
			}
		}
	}


	//////////////////////////////////////////////////////////////////////////
	// read fluid models
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find("FluidModels") != m_jsonData.end())
	{
		nlohmann::json fluidModels = m_jsonData["FluidModels"];
		for (auto& fluidModel : fluidModels)
		{
			std::string particleFile;
			if (readValue<std::string>(fluidModel["particleFile"], particleFile))
			{
				FluidData *data = new FluidData();
				data->samplesFile = particleFile;

				// id
				data->id = "Fluid";
				readValue(fluidModel["id"], data->id);

				// translation
				data->translation = Vector3r::Zero();
				readVector(fluidModel["translation"], data->translation);

				// rotation axis
				Vector3r axis = Vector3r::Zero();
				Real angle = 0.0;
				data->rotation = Matrix3r::Identity();
				if (readVector(fluidModel["rotationAxis"], axis) &&
					readValue<Real>(fluidModel["rotationAngle"], angle))
					data->rotation = AngleAxisr(angle, axis);

				// scale
				data->scale = 1.0;
				readValue(fluidModel["scale"], data->scale);

				scene.fluidModels.push_back(data);
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// read fluid blocks
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find("FluidBlocks") != m_jsonData.end())
	{
		nlohmann::json fluidBlocks = m_jsonData["FluidBlocks"];
		for (auto& fluidBlock : fluidBlocks)
		{
			// translation
			Vector3r translation = Vector3r::Zero();
			readVector(fluidBlock["translation"], translation);

			// scale
			Vector3r scale = Vector3r::Ones();
			readVector(fluidBlock["scale"], scale);

			Vector3r minX, maxX;
			if (readVector(fluidBlock["start"], minX) &&
				readVector(fluidBlock["end"], maxX))
			{
				FluidBlock *block = new FluidBlock();
				block->box.m_minX[0] = scale[0] * minX[0] + translation[0];
				block->box.m_minX[1] = scale[1] * minX[1] + translation[1];
				block->box.m_minX[2] = scale[2] * minX[2] + translation[2];
				block->box.m_maxX[0] = scale[0] * maxX[0] + translation[0];
				block->box.m_maxX[1] = scale[1] * maxX[1] + translation[1];
				block->box.m_maxX[2] = scale[2] * maxX[2] + translation[2];

				// id
				block->id = "Fluid";
				readValue(fluidBlock["id"], block->id);

				readValue(fluidBlock["denseMode"], block->mode);

				// velocity
				block->initialVelocity = Vector3r::Zero();
				readVector(fluidBlock["initialVelocity"], block->initialVelocity);

				scene.fluidBlocks.push_back(block);
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// read emitters
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find("Emitters") != m_jsonData.end())
	{
		nlohmann::json emitters = m_jsonData["Emitters"];
		for (auto& emitter : emitters)
		{
			EmitterData *data = new EmitterData();

			// id
			data->id = "Fluid";
			readValue(emitter["id"], data->id);

			// width
			data->width = 5;
			readValue(emitter["width"], data->width);

			// height
			data->height = 5;
			readValue(emitter["height"], data->height);

			// translation
			data->x = Vector3r::Zero();
			readVector(emitter["translation"], data->x);

			// direction
			data->dir = Vector3r(1.0, 0.0, 0.0);
			readVector(emitter["direction"], data->dir);

			// velocity
			data->v = Vector3r(1.0, 0.0, 0.0);
			readVector(emitter["velocity"], data->v);

			// emits per second
			data->emitsPerSecond = 10;
			readValue(emitter["emitsPerSecond"], data->emitsPerSecond);

			// type: 0 = rectangular, 1 = circle
			data->type = 0;
			readValue(emitter["type"], data->type);

			scene.emitters.push_back(data);
		}
	}
}


template <>
bool SceneLoader::readValue(const nlohmann::json &j, bool &v)
{
	if (j.is_null())
		return false;

	if (j.is_number_integer())
	{
		int val = j.get<int>();
		v = val != 0;
	}
	else
		v = j.get<bool>();
	return true;
}

void SceneLoader::readParameterObject(const std::string &key, ParameterObject *paramObj)
{
	if (paramObj == nullptr)
		return;

	const unsigned int numParams = paramObj->numParameters();


	//////////////////////////////////////////////////////////////////////////
	// read configuration 
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find(key) != m_jsonData.end())
	{
		nlohmann::json config = m_jsonData[key];
		std::vector<std::string> newParamList;

		for (unsigned int i = 0; i < numParams; i++)
		{
			ParameterBase *paramBase = paramObj->getParameter(i);

			if (paramBase->getType() == RealParameterType)
			{
				Real val;
				if (readValue(config[paramBase->getName()], val))
					static_cast<NumericParameter<Real>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::UINT32)
			{
				unsigned int val;
				if (readValue(config[paramBase->getName()], val))
					static_cast<NumericParameter<unsigned int>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::UINT16)
			{
				unsigned short val;
				if (readValue(config[paramBase->getName()], val))
					static_cast<NumericParameter<unsigned short>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::UINT8)
			{
				unsigned char val;
				if (readValue(config[paramBase->getName()], val))
					static_cast<NumericParameter<unsigned char>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::INT32)
			{
				int val;
				if (readValue(config[paramBase->getName()], val))
					static_cast<NumericParameter<int>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::INT16)
			{
				short val;
				if (readValue(config[paramBase->getName()], val))
					static_cast<NumericParameter<short>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::INT8)
			{
				char val;
				if (readValue(config[paramBase->getName()], val))
					static_cast<NumericParameter<char>*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::ENUM)
			{
				int val;
				if (readValue(config[paramBase->getName()], val))
					static_cast<EnumParameter*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == ParameterBase::BOOL)
			{
				bool val;
				if (readValue(config[paramBase->getName()], val))
					static_cast<BoolParameter*>(paramBase)->setValue(val);
			}
			else if (paramBase->getType() == RealVectorParameterType)
			{
				if (static_cast<VectorParameter<Real>*>(paramBase)->getDim() == 3)
				{
					Vector3r val;
					if (readVector(config[paramBase->getName()], val))
						static_cast<VectorParameter<Real>*>(paramBase)->setValue(val.data());
				}
			}
			else if (paramBase->getType() == ParameterBase::STRING)
			{
				std::string val;
				if (readValue(config[paramBase->getName()], val))
					static_cast<StringParameter*>(paramBase)->setValue(val);
			}
		}
	}
}

Utilities::SceneLoader::ColoringData Utilities::SceneLoader::readColoringInfo(const std::string &key)
{
	ColoringData data;
	data.colorField = "velocity";
	data.colorMapType = 0;
	data.minVal = 0.0;
	data.maxVal = 5.0;

	//////////////////////////////////////////////////////////////////////////
	// read configuration 
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find(key) != m_jsonData.end())
	{
		nlohmann::json config = m_jsonData[key];
#
		readValue(config["renderMinValue"], data.minVal);
		readValue(config["renderMaxValue"], data.maxVal);
		readValue(config["colorField"], data.colorField);
		readValue(config["colorMapType"], data.colorMapType);
	}
	return data;
}
