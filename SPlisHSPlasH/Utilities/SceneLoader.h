#ifndef __SceneLoader_h__
#define __SceneLoader_h__

#include "SPlisHSPlasH/Common.h"
#include "extern/json/json.hpp"
#include <vector>
#include "ParameterObject.h"

namespace Utilities
{
	/** \brief Importer of SPlisHSPlasH scene files. 
	*/
	class SceneLoader
	{
	protected:
		nlohmann::json m_jsonData;

	public:
		/** \brief Struct for an AABB */
		struct Box
		{
			Vector3r m_minX;
			Vector3r m_maxX;
		};

		/** \brief Struct to store a boundary object */
		struct BoundaryData
		{
			std::string samplesFile;
			std::string meshFile;
			Vector3r translation;
			Matrix3r rotation;
			Vector3r scale;
			Real density;
			bool dynamic;
			bool isWall;
			Eigen::Vector4f color;
			void *rigidBody;
		};

		/** \brief Struct to store a fluid object */
		struct FluidData
		{
			std::string id;
			std::string samplesFile;
			Vector3r translation;
			Matrix3r rotation;
			Real scale;
		};

		/** \brief Struct to store a fluid block */
		struct FluidBlock
		{
			std::string id;
			Box box;
			unsigned char mode;
			Vector3r initialVelocity;
		};

		/** \brief Struct to store an emitter object */
		struct EmitterData
		{
			std::string id;
			unsigned int width;
			unsigned int height;
			Vector3r x;
			Vector3r dir;
			Vector3r v;
			Real emitsPerSecond;
			unsigned int type;
		};

		/** \brief Struct to store scene information */
		struct Scene
		{
			std::vector<BoundaryData*> boundaryModels;
			std::vector<FluidData*> fluidModels;
			std::vector<FluidBlock*> fluidBlocks;
			std::vector<EmitterData*> emitters;
			Real particleRadius;
			Real timeStepSize;
		};


		void readScene(const char *fileName, Scene &scene);

		template <typename T>
		bool readValue(const nlohmann::json &j, T &v)
		{
			if (j.is_null())
				return false;

			v = j.get<T>();
			return true;
		}

		template <typename T, int size>
		bool readVector(const nlohmann::json &j, Eigen::Matrix<T, size, 1> &vec)
		{
			if (j.is_null())
				return false;

			std::vector<T> values = j.get<std::vector<T>>();
			for (unsigned int i = 0; i < values.size(); i++)
				vec[i] = values[i];
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

				v = j2.get<T>();
				return true;
			}
			return false;
		}

		template <typename T, int size>
		bool readVector(const std::string &section, const std::string &key, Eigen::Matrix<T, size, 1> &vec)
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

		void readParameterObject(const std::string &key, GenParam::ParameterObject *paramObj);
	};

	template <>
	bool SceneLoader::readValue<bool>(const nlohmann::json &j, bool &v);

}

#endif
