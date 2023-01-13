#ifndef __SceneLoader_h__
#define __SceneLoader_h__

#include "SPlisHSPlasH/Common.h"
#include "extern/json/json.hpp"
#include <vector>
#include "SceneParameterObjects.h"
#include "Utilities/Logger.h"

namespace Utilities
{
	/** \brief Importer of SPlisHSPlasH scene files. 
	*/
	class SceneLoader
	{
	protected:
		nlohmann::json m_jsonData;

		void readParameterObject(nlohmann::json& config, GenParam::ParameterObject* paramObj);

	public:
		/** \brief Struct to store scene information */
		struct Scene
		{
			std::vector<BoundaryParameterObject*> boundaryModels;
			std::vector<FluidModelParameterObject*> fluidModels;
			std::vector<FluidBlockParameterObject*> fluidBlocks;
			std::vector<EmitterParameterObject*> emitters;
			std::vector<AnimationFieldParameterObject*> animatedFields;
			std::vector<MaterialParameterObject*> materials;
			Real particleRadius;
			bool sim2D;
		};

		nlohmann::json& getJSONData() { return m_jsonData; };
		void readScene(const char *fileName, Scene &scene);

		template <typename T>
		bool readValue(const nlohmann::json &j, T &v)
		{
			if (j.is_null())
				return false;

			try
			{
				v = j.get<T>();
			}
			catch (const std::exception& e)
			{
				LOG_ERR << e.what();
				exit(1);
			}
			return true;
		}

		template <typename T, int size>
		bool readVector(const nlohmann::json &j, Eigen::Matrix<T, size, 1, Eigen::DontAlign> &vec)
		{
			if (j.is_null())
				return false;

			try
			{
				std::vector<T> values = j.get<std::vector<T>>();
				for (unsigned int i = 0; i < values.size(); i++)
					vec[i] = values[i];
			}
			catch (const std::exception& e)
			{
				LOG_ERR << e.what();
				exit(1);
			}
			return true;
		}

		template <typename T>
		bool readValue(const std::string &section, const std::string &key, T &v)
		{
			if (m_jsonData.find(section) != m_jsonData.end())
			{
				nlohmann::json j = m_jsonData[section];
				if (j.is_null())
					return false;

				nlohmann::json j2 = j[key];
				if (j2.is_null())
					return false;

				try
				{
					v = j2.get<T>();
				}
				catch (const std::exception& e)
				{
					LOG_ERR << e.what();
					exit(1);
				}
				return true;
			}
			return false;
		}

		bool hasValue(const std::string& section, const std::string& key)
		{
			if (m_jsonData.find(section) != m_jsonData.end())
			{
				nlohmann::json j = m_jsonData[section];
				if (j.is_null())
					return false;

				nlohmann::json j2 = j[key];
				if (j2.is_null())
					return false;

				return true;
			}
			return false;
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

		void readMaterialParameterObject(const std::string& key, GenParam::ParameterObject* paramObj);
		void readParameterObject(const std::string &key, GenParam::ParameterObject *paramObj);
	};

	template <>
	bool SceneLoader::readValue<bool>(const nlohmann::json &j, bool &v);

}

#endif
