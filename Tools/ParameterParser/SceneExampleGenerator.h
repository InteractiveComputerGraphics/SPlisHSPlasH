#pragma once

#include "SPlisHSPlasH/Common.h"
#include "ParameterObject.h"
#include "extern/json/json.hpp"

namespace SPH
{
	class SceneExampleGenerator
	{
	protected:
		nlohmann::json m_jsonData;
		nlohmann::json* m_currentData;

		void jsonEnumParam(GenParam::EnumParameter* param);
		void jsonStringParam(GenParam::StringParameter* param);
		void jsonVecParam(GenParam::RealVectorParameter* param);
		void jsonVecUintParam(GenParam::VectorParameter<unsigned int>* param);

		template<typename T>
		void jsonParam(T* param)
		{
			std::string str1 = param->getName() + " - comment";
			std::string str2 = param->getDescription() + "(Default: " + std::to_string(param->getValue()) + ", Type: " + typeid(param->getValue()).name() + ")";
			(*m_currentData)[str1] = str2;
			(*m_currentData)[param->getName()] = param->getValue();
		}

		void jsonParameterObject(const GenParam::ParameterObject* obj);

	public:
		void generateExampleSceneFile(const std::string& fileName);
	};
}
