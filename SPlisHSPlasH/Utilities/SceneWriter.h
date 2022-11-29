#ifndef __SceneWriter_h__
#define __SceneWriter_h__

#include "SPlisHSPlasH/Common.h"
#include "extern/json/json.hpp"
#include <vector>
#include "ParameterObject.h"

namespace Utilities
{
	/** \brief Importer of SPlisHSPlasH scene files. 
	*/
	class SceneWriter
	{
	protected:
		nlohmann::json m_jsonData;

		void writeParameterObject(nlohmann::json& config, GenParam::ParameterObject* paramObj);

	public:
		SceneWriter(const nlohmann::json& config) : m_jsonData(config) {}

		void writeScene(const char *fileName);

		template <typename T>
		bool writeValue(nlohmann::json &j, const std::string& key, const T &v)
		{
			if (j.is_null())
				return false;

			j[key] = v;
			return true;
		}

		template <typename T>
		bool writeVector(nlohmann::json& j, const std::string& key, const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& vec)
		{
			if (j.is_null())
				return false;

			j[key] = { vec[0], vec[1], vec[2] };
			return true;
		}

		template <typename T, int size>
		bool readVector(const std::string &section, const std::string &key, Eigen::Matrix<T, size, 1, Eigen::DontAlign> &vec)
		{
			if (m_jsonData.find(section) != m_jsonData.end())
			{
				nlohmann::json j = m_jsonData[section];
				if (j.is_null())
					return false;

				nlohmann::json j2 = j[key];
				if (j2.is_null())
					return false;

				std::vector<T> values = j2.get<std::vector<T>>();
				for (unsigned int i = 0; i < values.size(); i++)
					vec[i] = values[i];
				return true;
			}
			return false;
		}

		void updateMaterialParameterConfig(const std::string& key, GenParam::ParameterObject* paramObj);
		template <typename T>
		void updateMaterialParameterConfig(const std::string& id, const std::string& key, const T &v)
		{
			if (m_jsonData.find("Materials") != m_jsonData.end())
			{
				nlohmann::json& materials = m_jsonData["Materials"];
				for (auto& material : materials)
				{
					std::string mid = material["id"];
					if (mid == id)
					{
						writeValue(material, key, v);
					}
				}
			}
		}
		void writeParameterObject(const std::string &key, GenParam::ParameterObject *paramObj);
	};

}

#endif
