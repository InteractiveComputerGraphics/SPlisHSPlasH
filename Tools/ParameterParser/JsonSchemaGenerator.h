#pragma once

#include "SPlisHSPlasH/Common.h"
#include "ParameterObject.h"
#include "extern/json/json.hpp"

namespace SPH
{
	class JsonSchemaGenerator
	{
	protected:
		nlohmann::json m_jsonData;
		nlohmann::json* m_currentData;

		void jsonBoolParam(GenParam::BoolParameter* param);
		void jsonStringParam(GenParam::StringParameter* param);
		void jsonEnumParam(GenParam::EnumParameter* param);
		void jsonVecUintParam(GenParam::VectorParameter<unsigned int>* param);
		void jsonVecParam(GenParam::RealVectorParameter* param);

		template<typename T>
		void jsonNumericParam(T* param)
		{
			std::string str = param->getDescription() + " (Default: " + std::to_string(param->getValue()) + ", Type: " + typeid(param->getValue()).name() + ")";

			(*m_currentData)[param->getName()] = nullptr;
			(*m_currentData)[param->getName()]["description"] = str;
			if ((typeid(param->getValue()) == typeid(unsigned int)) ||
				(typeid(param->getValue()) == typeid(unsigned short)) ||
				(typeid(param->getValue()) == typeid(unsigned char)) ||
				(typeid(param->getValue()) == typeid(int)) ||
				(typeid(param->getValue()) == typeid(short)) ||
				(typeid(param->getValue()) == typeid(char)))
				(*m_currentData)[param->getName()]["type"] = "integer";
			else if (typeid(param->getValue()) == typeid(Real))
				(*m_currentData)[param->getName()]["type"] = "number";
			(*m_currentData)[param->getName()]["minimum"] = param->getMinValue();
			(*m_currentData)[param->getName()]["maximum"] = param->getMaxValue();
		}

		void jsonParameterObject(const GenParam::ParameterObject* obj);

	public:
		void generateSchemaFile(const std::string &fileName);
	};
}
