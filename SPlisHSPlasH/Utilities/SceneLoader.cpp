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

		if (scene.sim2D)
			scene.camPosition = Vector3r(0.0, 0.0, 8.0);
		else
			scene.camPosition = Vector3r(0.0, 3.0, 8.0);
		readVector(config["cameraPosition"], scene.camPosition);
		scene.camLookat = Vector3r(0.0, 0.0, 0.0);
		readVector(config["cameraLookat"], scene.camLookat);
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
			std::string mapFile = "";
			const bool bMesh = readValue<std::string>(boundaryModel["geometryFile"], meshFile);
			const bool bSamples = readValue<std::string>(boundaryModel["particleFile"], particleFile);
			const bool bMap = readValue<std::string>(boundaryModel["mapFile"], mapFile);

			if (bMesh || bSamples)
			{
				BoundaryData *data = new BoundaryData();
				data->meshFile = meshFile;
				data->samplesFile = particleFile;
				data->mapFile = mapFile;

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

				data->samplingMode = 0;
				readValue<unsigned int>(boundaryModel["samplingMode"], data->samplingMode);

				// Maps
				data->mapInvert = false;
				readValue(boundaryModel["mapInvert"], data->mapInvert);

				data->mapThickness = 0.0;
				readValue(boundaryModel["mapThickness"], data->mapThickness);

				data->mapResolution = Eigen::Matrix<unsigned int, 3, 1>(20, 20, 20);
				readVector(boundaryModel["mapResolution"], data->mapResolution);


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
				data->scale = Vector3r::Ones();
				readVector(fluidModel["scale"], data->scale);

				// velocity
				data->initialVelocity = Vector3r::Zero();
				readVector(fluidModel["initialVelocity"], data->initialVelocity);

				data->invert = false;
				readValue(fluidModel["invert"], data->invert);

				data->resolutionSDF = { 20, 20, 20 };
				Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign> res(20,20,20);
				readVector(fluidModel["resolutionSDF"], res);
				data->resolutionSDF[0] = res[0];
				data->resolutionSDF[1] = res[1];
				data->resolutionSDF[2] = res[2];

				data->mode = 0;
				readValue(fluidModel["denseMode"], data->mode);

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
			
			// rotation
			// default direction without rotation is +x
			Vector3r axis(0,0,1);
			Real angle = 0.0;
			data->rotation = Matrix3r::Identity();
			if (readVector(emitter["rotationAxis"], axis) &&
				readValue(emitter["rotationAngle"], angle))
			{
				// in 2D simulations always rotate around z-axis
				if (scene.sim2D)
					axis = { 0.0, 0.0, 1.0 };
				axis.normalize();
				data->rotation = AngleAxisr(angle, axis).toRotationMatrix();
 			}

			// emission velocity
			data->velocity = 1;
			readValue(emitter["velocity"], data->velocity);

			// time when emission starts and stops
			data->emitStartTime = 0;
			readValue(emitter["emitStartTime"], data->emitStartTime);
			data->emitEndTime = std::numeric_limits<Real>::max();
			readValue(emitter["emitEndTime"], data->emitEndTime);

			// type: 0 = rectangular, 1 = circle
			data->type = 0;
			readValue(emitter["type"], data->type);

			scene.emitters.push_back(data);
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// read animation fields
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find("AnimationFields") != m_jsonData.end())
	{
		nlohmann::json fields = m_jsonData["AnimationFields"];
		for (auto& field : fields)
		{
			AnimationFieldData * data = new AnimationFieldData();

			data->particleFieldName = "";
			readValue(field["particleField"], data->particleFieldName);

			data->expression[0] = "";
			readValue(field["expression_x"], data->expression[0]);

			data->expression[1] = "";
			readValue(field["expression_y"], data->expression[1]);

			data->expression[2] = "";
			readValue(field["expression_z"], data->expression[2]);

			// 0=Box, 1=Cylinder
			data->shapeType = 0;
			readValue(field["shapeType"], data->shapeType);

			// time when emission starts and stops
			data->startTime = 0;
			readValue(field["startTime"], data->startTime);
			data->endTime = std::numeric_limits<Real>::max();
			readValue(field["endTime"], data->endTime);

			data->scale= Vector3r::Ones();
			readVector(field["scale"], data->scale);
			
			// shape position
			data->x = Vector3r::Zero();
			readVector(field["translation"], data->x);

			// rotation
			// default direction without rotation is +x
			Vector3r axis = Vector3r::Zero();
			Real angle = 0.0;
			data->rotation = Matrix3r::Identity();
			if (readVector(field["rotationAxis"], axis) &&
				readValue(field["rotationAngle"], angle))
			{
				axis.normalize();
				data->rotation = AngleAxisr(angle, axis).toRotationMatrix();
			}
			scene.animatedFields.push_back(data);
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
