#ifndef __SceneParameterObjects_h__
#define __SceneParameterObjects_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include <array>
#include "ParameterObject.h"

namespace Utilities
{
	class FluidBlockParameterObject : public GenParam::ParameterObject
	{
	public:
		std::string id;
		Vector3r boxMin;
		Vector3r boxMax;
		Vector3r translation;
		Vector3r scale;
		std::string visMeshFile;
		unsigned char mode;
		Vector3r initialVelocity;
		Vector3r initialAngularVelocity;

		FluidBlockParameterObject()
		{
			// Default values
			id = "Fluid";
			visMeshFile = "";
			boxMin = Vector3r::Zero();
			boxMax = Vector3r::Zero();
			mode = 0;
			translation = Vector3r::Zero();
			scale = Vector3r::Ones();
			initialVelocity = Vector3r::Zero();
			initialAngularVelocity = Vector3r::Zero();
		}

		FluidBlockParameterObject(std::string id_, std::string visMeshFile_, Vector3r boxMin_, Vector3r boxMax_, unsigned char mode_, Vector3r translation_, Vector3r scale_,
								Vector3r initialVelocity_, Vector3r initialAngularVelocity_)
		{
			// Default values
			id = id_;
			boxMin = boxMin_;
			boxMax = boxMax_;
			translation = translation_;
			scale = scale_;
			visMeshFile = visMeshFile_;
			mode = mode_;
			initialVelocity = initialVelocity_;
			initialAngularVelocity = initialAngularVelocity_;
		}

		static int FLUID_BLOCK_ID;
		static int FLUID_BLOCK_BOX_MINX;
		static int FLUID_BLOCK_BOX_MAXX;
		static int FLUID_BLOCK_TRANSLATION;
		static int FLUID_BLOCK_SCALE;
		static int FLUID_BLOCK_VISMESH;
		static int FLUID_BLOCK_MODE;
		static int FLUID_BLOCK_INITIAL_VEL;
		static int FLUID_BLOCK_INITIAL_ANGVEL;

		virtual void initParameters();
	};

	class FluidModelParameterObject : public GenParam::ParameterObject
	{
	public:
		std::string id;
		std::string visMeshFile;
		std::string samplesFile;
		Vector3r translation;
		Vector3r scale;		
		unsigned char mode;
		Vector3r initialVelocity;
		Vector3r initialAngularVelocity;
		Vector3r axis;
		Real angle;
		bool invert;
		std::array<unsigned int, 3> resolutionSDF;


		FluidModelParameterObject()
		{
			// Default values
			id = "Fluid";
			translation = Vector3r::Zero();
			axis = Vector3r(1, 0, 0);
			angle = 0.0;
			samplesFile = "";
			invert = false;
			resolutionSDF = { 20, 20, 20 };
			scale = Vector3r::Ones();
			visMeshFile = "";
			mode = 0;
			initialVelocity = Vector3r::Zero();
			initialAngularVelocity = Vector3r::Zero();
		}

		FluidModelParameterObject(std::string id_, std::string visMeshFile_, std::string samplesFile_, Vector3r translation_, Vector3r axis_, Real angle_, Vector3r scale_, 
								Vector3r initialVelocity_, Vector3r initialAngularVelocity_, unsigned char mode_, bool invert_, std::array<unsigned int, 3> resolutionSDF_)
		{
			// Default values
			id = id_;
			translation = translation_;
			axis = axis_;
			angle = angle_;
			samplesFile = samplesFile_;
			invert = invert_;
			resolutionSDF = resolutionSDF_;
			scale = scale_;
			visMeshFile = visMeshFile_;
			mode = mode_;
			initialVelocity = initialVelocity_;
			initialAngularVelocity = initialAngularVelocity_;
		}

		static int FLUID_MODEL_ID;
		static int FLUID_MODEL_TRANSLATION;
		static int FLUID_MODEL_SCALE;
		static int FLUID_MODEL_VISMESH;
		static int FLUID_MODEL_MODE;
		static int FLUID_MODEL_INITIAL_VEL;
		static int FLUID_MODEL_INITIAL_ANGVEL;
		static int FLUID_MODEL_SAMPLES_FILE;
		static int FLUID_MODEL_ROTAXIS;
		static int FLUID_MODEL_ROTANGLE;
		static int FLUID_MODEL_INVERT;
		static int FLUID_MODEL_RESSDF;

		virtual void initParameters();
	};

	class EmitterParameterObject : public GenParam::ParameterObject
	{
	public:
		std::string id;
		unsigned int width;
		unsigned int height;
		Vector3r x;
		Real velocity; // emission velocity
		Vector3r axis;
		Real angle;
		Real emitStartTime;
		Real emitEndTime;
		unsigned int type;	// type: 0 = rectangular, 1 = circle

		EmitterParameterObject()
		{
			// Default values
			id = "Fluid";
			width = 5;
			height = 5;
			x = Vector3r::Zero();
			velocity = 1.0;
			axis = Vector3r(0, 0, 1);
			angle = 0.0;
			emitStartTime = 0.0;
			emitEndTime = std::numeric_limits<Real>::max();
			type = 0;
		}

		EmitterParameterObject(std::string id_, unsigned int width_, unsigned int height_, Vector3r x_, Real velocity_, Vector3r axis_, 
							Real angle_, Real emitStartTime_, Real emitEndTime_, unsigned int type_)
		{
			id = id_;
			width = width_;
			height = height_;
			x = x_;
			velocity = velocity_;
			axis = axis_;
			angle = angle_;
			emitStartTime = emitStartTime_;
			emitEndTime = emitEndTime_;
			type = type_;
		}

		static int EMITTER_ID;
		static int EMITTER_WIDTH;
		static int EMITTER_HEIGHT;
		static int EMITTER_POSITION;
		static int EMITTER_VELOCITY;
		static int EMITTER_ROTAXIS;
		static int EMITTER_ROTANGLE;
		static int EMITTER_STARTTIME;
		static int EMITTER_ENDTIME;
		static int EMITTER_TYPE;

		virtual void initParameters();
	};


	/** \brief class to store an animation field object
	 */
	class AnimationFieldParameterObject : public GenParam::ParameterObject
	{
	public:
		std::string particleFieldName;
		std::string expression[3];
		unsigned int shapeType;
		Vector3r translation;
		Vector3r axis;
		Real angle;
		Vector3r scale;
		Real startTime;
		Real endTime;

		AnimationFieldParameterObject()
		{
			// Default values
			particleFieldName = "";
			expression[0] = "";
			expression[1] = "";
			expression[2] = "";
			// 0=Box, 1=Cylinder
			shapeType = 0;

			// time when emission starts and stops
			startTime = 0;
			endTime = std::numeric_limits<Real>::max();
			scale = Vector3r::Ones();

			// shape position
			translation = Vector3r::Zero();

			// rotation
			// default direction without rotation is +x
			axis = Vector3r(1.0, 0.0, 0.0);
			angle = 0.0;
		}

		AnimationFieldParameterObject(std::string particleFieldName_, std::string expressionX_, std::string expressionY_, std::string expressionZ_,
								unsigned int shapeType_, Vector3r translation_, Vector3r axis_, Real angle_, Vector3r scale_, Real startTime_, Real endTime_)
		{
			particleFieldName = particleFieldName_;
			expression[0] = expressionX_;
			expression[1] = expressionY_;
			expression[2] = expressionZ_;
			shapeType = shapeType_;
			startTime = startTime_;
			endTime = endTime_;
			scale = scale_;
			translation = translation_;
			axis = axis_;
			angle = angle_;
		}

		static int ANIMATIONFIELD_PARTICLE_FIELD;
		static int ANIMATIONFIELD_EXPR_X;
		static int ANIMATIONFIELD_EXPR_Y;
		static int ANIMATIONFIELD_EXPR_Z;
		static int ANIMATIONFIELD_SHAPETYPE;
		static int ANIMATIONFIELD_TRANSLATION;
		static int ANIMATIONFIELD_ROTAXIS;
		static int ANIMATIONFIELD_ROTANGLE;
		static int ANIMATIONFIELD_SCALE;
		static int ANIMATIONFIELD_STARTTIME;
		static int ANIMATIONFIELD_ENDTIME;

		virtual void initParameters();
	};

	/** \brief Class to store particle coloring information */
	class MaterialParameterObject : public GenParam::ParameterObject
	{
	public:
		std::string id;
		std::string colorField;
		unsigned int colorMapType;
		Real minVal;
		Real maxVal;
		unsigned int maxEmitterParticles;
		bool emitterReuseParticles;
		Vector3r emitterBoxMin;
		Vector3r emitterBoxMax;

		MaterialParameterObject()
		{
			// Default values
			id = "Fluid";
			minVal = 0.0;
			maxVal = 10.0;
			colorField = "velocity";
			colorMapType = 1;
			maxEmitterParticles = 10000;
			emitterReuseParticles = false;
			emitterBoxMin = Vector3r(-1.0, -1.0, -1.0);
			emitterBoxMax = Vector3r(1.0, 1.0, 1.0);
		}

		MaterialParameterObject(std::string id_, std::string colorField_, unsigned int colorMapType_, Real minVal_, Real maxVal_,
							unsigned int maxEmitterParticles_, bool emitterReuseParticles_, Vector3r emitterBoxMin_, Vector3r emitterBoxMax_)
		{
			id = id_;
			minVal = minVal_;
			maxVal = maxVal_;
			colorField = colorField_;
			colorMapType = colorMapType_;
			maxEmitterParticles = maxEmitterParticles_;
			emitterReuseParticles = emitterReuseParticles_;
			emitterBoxMin = emitterBoxMin_;
			emitterBoxMax = emitterBoxMax_;
		}

		static int MATERIAL_ID;
		static int MATERIAL_MIN_VAL;
		static int MATERIAL_MAX_VAL;
		static int MATERIAL_COLOR_FIELD;
		static int MATERIAL_COLOR_MAP;
		static int MATERIAL_MAX_EMITTER_PARTICLES;
		static int MATERIAL_EMITTER_REUSE;
		static int MATERIAL_EMITTER_BOX_MIN;
		static int MATERIAL_EMITTER_BOX_MAX;

		virtual void initParameters();
	};

	/** \brief Struct to store a boundary object */
	class BoundaryParameterObject : public GenParam::ParameterObject
	{
	public:
		std::string samplesFile;
		std::string meshFile;
		Vector3r translation;
		Vector3r axis;
		Real angle;
		Vector3r scale;
		bool dynamic;
		bool isWall;
		Vector4r color;
		std::string mapFile;
		bool mapInvert;
		Real mapThickness;
		Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign> mapResolution;
		unsigned int samplingMode;
		bool isAnimated;

		BoundaryParameterObject()
		{
			// Default values
			samplesFile = "";
			meshFile = "";
			translation = Vector3r::Zero();
			axis = Vector3r(1, 0, 0);
			angle = 0.0;
			scale = Vector3r::Ones();
			dynamic = false;
			isWall = false;
			color = Vector4r(1.0, 0.0, 0.0, 0.0);
			samplingMode = 0;
			isAnimated = false;
			// Maps
			mapFile = "";
			mapInvert = false;
			mapThickness = 0.0;
			mapResolution = Eigen::Matrix<unsigned int, 3, 1>(20, 20, 20);
		}

		BoundaryParameterObject(std::string samplesFile_, std::string meshFile_, Vector3r translation_, Vector3r axis_, Real angle_, Vector3r scale_,
								bool dynamic_, bool isWall_, Vector4r color_, std::string mapFile_, bool mapInvert_,
								Real mapThickness_, Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign> mapResolution_, unsigned int samplingMode_, bool isAnimated_)
		{
			samplesFile = samplesFile_;
			meshFile = meshFile_;
			translation = translation_;
			axis = axis_;
			angle = angle_;
			scale = scale_;
			dynamic = dynamic_;
			isWall = isWall_;
			color = color_;
			samplingMode = samplingMode_;
			isAnimated = isAnimated_;
			// Maps
			mapFile = mapFile_;
			mapInvert = mapInvert_;
			mapThickness = mapThickness_;
			mapResolution = mapResolution_;
		}

		static int BOUNDARY_SAMPLES_FILE;
		static int BOUNDARY_MESH_FILE;
		static int BOUNDARY_TRANSLATION;
		static int BOUNDARY_AXIS;
		static int BOUNDARY_ANGLE;
		static int BOUNDARY_SCALE;
		static int BOUNDARY_DYNAMIC;
		static int BOUNDARY_IS_WALL;
		static int BOUNDARY_COLOR;
		static int BOUNDARY_MAP_FILE;
		static int BOUNDARY_MAP_INVERT;
		static int BOUNDARY_MAP_THICKNESS;
		static int BOUNDARY_MAP_RESOLUTION;
		static int BOUNDARY_SAMPLING_MODE;
		static int BOUNDARY_IS_ANIMATED;

		virtual void initParameters();
	};
}

#endif
