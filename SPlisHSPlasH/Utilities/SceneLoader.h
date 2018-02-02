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
			std::string samplesFile;
			Vector3r translation;
			Matrix3r rotation;
			Real scale;
		};

		/** \brief Struct to store a fluid block */
		struct FluidBlock
		{
			Box box;
			unsigned char mode;
			Vector3r initialVelocity;
		};

		/** \brief Struct to store an emitter object */
		struct EmitterData
		{
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
			unsigned int maxEmitterParticles;
			bool emitterReuseParticles;
			Vector3r emitterBoxMin;
			Vector3r emitterBoxMax;
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
			unsigned int index = 0;
			if (j.is_null())
				return false;

			std::vector<T> values = j.get<std::vector<T>>();
			for (unsigned int i = 0; i < values.size(); i++)
				vec[i] = values[i];
			return true;
		}

		void readParameterObject(GenParam::ParameterObject *paramObj);
	};

	template <>
	bool SceneLoader::readValue<bool>(const nlohmann::json &j, bool &v);

}

#endif
