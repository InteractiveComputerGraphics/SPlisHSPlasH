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
			Eigen::Matrix<float, 4, 1, Eigen::DontAlign> color;
			void *rigidBody;

			std::string mapFile;
			bool mapInvert;
			Real mapThickness;	
			Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign> mapResolution;
			unsigned int samplingMode;
		};

		/** \brief Struct to store a fluid object */
		struct FluidData
		{
			std::string id;
			std::string samplesFile;
			Vector3r translation;
			Matrix3r rotation;
			Vector3r scale;
			Vector3r initialVelocity;
			unsigned char mode;
			bool invert;
			std::array<unsigned int, 3> resolutionSDF;
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
			Real velocity; // emission velocity
			Matrix3r rotation;
			Real emitStartTime;
			Real emitEndTime;
			unsigned int type;
		};

		/** \brief Struct to store an animation field object
		 */
		struct AnimationFieldData
		{
			std::string particleFieldName;
			std::string expression[3];
			unsigned int shapeType;
			Vector3r x;
			Matrix3r rotation;
			Vector3r scale; 
			Real startTime;
			Real endTime;
		};

		/** \brief Struct to store scene information */
		struct Scene
		{
			std::vector<BoundaryData*> boundaryModels;
			std::vector<FluidData*> fluidModels;
			std::vector<FluidBlock*> fluidBlocks;
			std::vector<EmitterData*> emitters;
			std::vector<AnimationFieldData*> animatedFields;
			Real particleRadius;
			bool sim2D;
			Real timeStepSize;
			Vector3r camPosition;
			Vector3r camLookat;
		};

		/** \brief Struct to store particle coloring information */
		struct ColoringData
		{
			std::string colorField;
			unsigned int colorMapType;
			Real minVal;
			Real maxVal;
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
		bool readVector(const nlohmann::json &j, Eigen::Matrix<T, size, 1, Eigen::DontAlign> &vec)
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

		void readParameterObject(const std::string &key, GenParam::ParameterObject *paramObj);
		ColoringData readColoringInfo(const std::string &key);
	};

	template <>
	bool SceneLoader::readValue<bool>(const nlohmann::json &j, bool &v);

}

#endif
