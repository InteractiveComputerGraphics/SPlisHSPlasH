#include "SceneWriter.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include <iostream>
#include "extern/json/json.hpp"
#include <fstream>
#include "Utilities/FileSystem.h"
#include "Utilities/Logger.h"

using namespace Utilities;
using namespace GenParam;
using namespace std;
using namespace SPH;


void SceneWriter::writeScene(const char* fileName)
{
	LOG_INFO << "Write scene file: " << fileName;


	//////////////////////////////////////////////////////////////////////////
	// update configuration 
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find("Configuration") != m_jsonData.end())
	{
		nlohmann::json &config = m_jsonData["Configuration"];
		config["timeStepSize"] = TimeManager::getCurrent()->getTimeStepSize();
	}


	std::ofstream output_file(fileName);
	if (!output_file.is_open())
	{
		LOG_ERR << "Cannot open file!";
		return;
	}
	try
	{
		output_file << std::setw(4) << m_jsonData << std::endl;
	}
	catch (const std::exception& e)
	{
		LOG_ERR << e.what();
		exit(1);
	}
	return;
}

void SceneWriter::updateMaterialParameterConfig(const std::string& key, ParameterObject* paramObj)
{
	if (paramObj == nullptr)
		return;

	if (m_jsonData.find("Materials") != m_jsonData.end())
	{
		nlohmann::json &materials = m_jsonData["Materials"];
		for (auto& material : materials)
		{
			string id = material["id"];
			if (key == id)
			{
				writeParameterObject(material, paramObj);
			}
		}
	}
}

void SceneWriter::writeParameterObject(const std::string& key, ParameterObject* paramObj)
{
	if (paramObj == nullptr)
		return;

	//////////////////////////////////////////////////////////////////////////
	// write configuration 
	//////////////////////////////////////////////////////////////////////////
	if (m_jsonData.find(key) != m_jsonData.end())
	{
		nlohmann::json &config = m_jsonData[key];
		writeParameterObject(config, paramObj);
	}
}


void SceneWriter::writeParameterObject(nlohmann::json& config, ParameterObject* paramObj)
{
	if (paramObj == nullptr)
		return;

	const unsigned int numParams = paramObj->numParameters();
	for (unsigned int i = 0; i < numParams; i++)
	{
		ParameterBase* paramBase = paramObj->getParameter(i);

		// only update writable values
		if (paramBase->getReadOnly())
			continue;

		if (paramBase->getType() == RealParameterType)
		{
			const Real val = static_cast<NumericParameter<Real>*>(paramBase)->getValue();
			writeValue(config, paramBase->getName(), val);
		}
		else if (paramBase->getType() == ParameterBase::UINT32)
		{
			const unsigned int val = static_cast<NumericParameter<Real>*>(paramBase)->getValue();
			writeValue(config, paramBase->getName(), val);
		}
		else if (paramBase->getType() == ParameterBase::UINT16)
		{
			const unsigned short val = static_cast<NumericParameter<Real>*>(paramBase)->getValue();
			writeValue(config, paramBase->getName(), val);
		}
		else if (paramBase->getType() == ParameterBase::UINT8)
		{
			const unsigned char val = static_cast<NumericParameter<Real>*>(paramBase)->getValue();
			writeValue(config, paramBase->getName(), val);
		}
		else if (paramBase->getType() == ParameterBase::INT32)
		{
			const int val = static_cast<NumericParameter<Real>*>(paramBase)->getValue();
			writeValue(config, paramBase->getName(), val);
		}
		else if (paramBase->getType() == ParameterBase::INT16)
		{
			const short val = static_cast<NumericParameter<Real>*>(paramBase)->getValue();
			writeValue(config, paramBase->getName(), val);
		}
		else if (paramBase->getType() == ParameterBase::INT8)
		{
			const char val = static_cast<NumericParameter<Real>*>(paramBase)->getValue();
			writeValue(config, paramBase->getName(), val);
		}
		else if (paramBase->getType() == ParameterBase::ENUM)
		{
			const int val = static_cast<EnumParameter*>(paramBase)->getValue();
			writeValue(config, paramBase->getName(), val);
		}
		else if (paramBase->getType() == ParameterBase::BOOL)
		{
			const bool val = static_cast<BoolParameter*>(paramBase)->getValue();
			writeValue(config, paramBase->getName(), val);
		}
		else if (paramBase->getType() == RealVectorParameterType)
		{
			if (static_cast<VectorParameter<Real>*>(paramBase)->getDim() == 3)
			{
				const Vector3r val(static_cast<VectorParameter<Real>*>(paramBase)->getValue());
				writeVector(config, paramBase->getName(), val);
			}
		}
	}
}
