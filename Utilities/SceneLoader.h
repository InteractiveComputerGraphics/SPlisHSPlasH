#ifndef __SceneLoader_h__
#define __SceneLoader_h__

#include "SPlisHSPlasH/Common.h"
#include "extern/json/json.hpp"
#include <vector>

namespace SPH
{
	/** \brief Importer of SPlisHSPlasH scene files. 
	*/
	class SceneLoader
	{
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
		};

		/** \brief Struct to store scene information */
		struct Scene
		{
			std::vector<BoundaryData*> boundaryModels;
			std::vector<FluidData*> fluidModels;
			std::vector<FluidBlock*> fluidBlocks;
			Real particleRadius;
			Real pauseAt;
			unsigned int numberOfStepsPerRenderUpdate;
			unsigned int cflMethod;
			Real cflFactor;
			Real cflMaxTimeStepSize;
			Real maxError;
			unsigned int maxIterations;
			Real maxErrorV;
			unsigned int maxIterationsV;
			Real viscosity;
			Real surfaceTension;
			Real density0;
			unsigned int velocityUpdateMethod;
			Real stiffness;
			Real exponent;
			bool enableDivergenceSolver;
			Vector3r gravitation;
			Real timeStepSize;
			unsigned int viscosityMethod;
			unsigned int surfaceTensionMethod;
			unsigned int simulationMethod;
		};


		static void readScene(const char *fileName, Scene &scene);

		template <typename T>
		static bool readValue(const nlohmann::json &j, T &v)
		{
			if (j.is_null())
				return false;

			v = j.get<T>();
			return true;
		}

		template <typename T, int size>
		static bool readVector(const nlohmann::json &j, Eigen::Matrix<T, size, 1> &vec)
		{
			unsigned int index = 0;
			if (j.is_null())
				return false;

			std::vector<T> values = j.get<std::vector<T>>();
			for (unsigned int i = 0; i < values.size(); i++)
				vec[i] = values[i];
			return true;
		}
	};

	template <>
	bool SceneLoader::readValue<bool>(const nlohmann::json &j, bool &v);

}

#endif
