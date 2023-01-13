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
		m_jsonData = nlohmann::json::parse(input_file);
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
			BoundaryParameterObject* data = new BoundaryParameterObject();
			data->initParameters();
			readParameterObject(boundaryModel, data);
			data->axis.normalize();
			if ((data->meshFile != "") || (data->samplesFile != ""))
				scene.boundaryModels.push_back(data);
			else
				delete data;
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
			FluidModelParameterObject* data = new FluidModelParameterObject();
			data->initParameters();
			readParameterObject(fluidModel, data);
			data->axis.normalize();
			scene.fluidModels.push_back(data);
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
			FluidBlockParameterObject *block = new FluidBlockParameterObject();
			block->initParameters();
			readParameterObject(fluidBlock, block);
			scene.fluidBlocks.push_back(block);
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
			EmitterParameterObject *data = new EmitterParameterObject();
			data->initParameters();
			readParameterObject(emitter, data);

			// in 2D simulations always rotate around z-axis
			if (scene.sim2D)
				data->axis = { 0.0, 0.0, 1.0 };
			data->axis.normalize();
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
			AnimationFieldParameterObject* data = new AnimationFieldParameterObject();
			data->initParameters();
			readParameterObject(field, data);
			data->axis.normalize();
			scene.animatedFields.push_back(data);
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// read materials
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find("Materials") != m_jsonData.end())
	{
		nlohmann::json materials = m_jsonData["Materials"];
		for (auto& material : materials)
		{
			MaterialParameterObject* data = new MaterialParameterObject();
			data->initParameters();
			readParameterObject(material, data);
			scene.materials.push_back(data);
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

void SceneLoader::readMaterialParameterObject(const std::string &key, ParameterObject *paramObj)
{
	if (paramObj == nullptr)
		return;

	if (m_jsonData.find("Materials") != m_jsonData.end())
	{
		nlohmann::json materials = m_jsonData["Materials"];
		for (auto& material : materials)
		{
			string id = "";
			readValue(material["id"], id);

			if (key == id)
			{
				readParameterObject(material, paramObj);
			}
		}
	}
}

void SceneLoader::readParameterObject(const std::string& key, ParameterObject* paramObj)
{
	if (paramObj == nullptr)
		return;

	//////////////////////////////////////////////////////////////////////////
	// read configuration 
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find(key) != m_jsonData.end())
	{
		nlohmann::json config = m_jsonData[key];
		readParameterObject(config, paramObj);
	}
}

void SceneLoader::readParameterObject(nlohmann::json& config, ParameterObject* paramObj)
{
	if (paramObj == nullptr)
		return;

	const unsigned int numParams = paramObj->numParameters();
	for (unsigned int i = 0; i < numParams; i++)
	{
		ParameterBase* paramBase = paramObj->getParameter(i);

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
			else if (static_cast<VectorParameter<Real>*>(paramBase)->getDim() == 4)
			{
				Vector4r val;
				if (readVector(config[paramBase->getName()], val))
					static_cast<VectorParameter<Real>*>(paramBase)->setValue(val.data());
			}
		}
		else if (paramBase->getType() == ParameterBase::VEC_UINT32)
		{
			if (static_cast<VectorParameter<unsigned int>*>(paramBase)->getDim() == 3)
			{
				Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign> val;
				if (readVector(config[paramBase->getName()], val))
					static_cast<VectorParameter<unsigned int>*>(paramBase)->setValue(val.data());
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




