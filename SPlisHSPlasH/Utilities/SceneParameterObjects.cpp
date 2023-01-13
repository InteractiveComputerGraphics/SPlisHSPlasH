#include "SceneParameterObjects.h"

using namespace Utilities;
using namespace GenParam;
using namespace std;


//////////////////////////////////////////////////////////////////////////
// FluidBlockParameterObject
//////////////////////////////////////////////////////////////////////////
int FluidBlockParameterObject::FLUID_BLOCK_ID = -1;
int FluidBlockParameterObject::FLUID_BLOCK_BOX_MINX = -1;
int FluidBlockParameterObject::FLUID_BLOCK_BOX_MAXX = -1;
int FluidBlockParameterObject::FLUID_BLOCK_TRANSLATION = -1;
int FluidBlockParameterObject::FLUID_BLOCK_SCALE = -1;
int FluidBlockParameterObject::FLUID_BLOCK_VISMESH = -1;
int FluidBlockParameterObject::FLUID_BLOCK_MODE = -1;
int FluidBlockParameterObject::FLUID_BLOCK_INITIAL_VEL = -1;
int FluidBlockParameterObject::FLUID_BLOCK_INITIAL_ANGVEL = -1;

void FluidBlockParameterObject::initParameters()
{
	FLUID_BLOCK_ID = createStringParameter("id", "ID", &id);
	setGroup(FLUID_BLOCK_ID, "FluidBlock");
	setDescription(FLUID_BLOCK_ID, "This id is used to define the properties of the fluid block. If no id is defined, then the standard id \"Fluid\" is used.");

	FLUID_BLOCK_BOX_MINX = createVectorParameter("start", "Start", 3u, boxMin.data());
	setGroup(FLUID_BLOCK_BOX_MINX, "FluidBlock");
	setDescription(FLUID_BLOCK_BOX_MINX, "Minimum coordinate of the box which defines the fluid block.");

	FLUID_BLOCK_BOX_MAXX = createVectorParameter("end", "End", 3u, boxMax.data());
	setGroup(FLUID_BLOCK_BOX_MAXX, "FluidBlock");
	setDescription(FLUID_BLOCK_BOX_MAXX, "Maximum coordinate of the box which defines the fluid block.");

	FLUID_BLOCK_TRANSLATION = createVectorParameter("translation", "Translation", 3u, translation.data());
	setGroup(FLUID_BLOCK_TRANSLATION, "FluidBlock");
	setDescription(FLUID_BLOCK_TRANSLATION, "Translation vector of the block.");

	FLUID_BLOCK_SCALE = createVectorParameter("scale", "Scale", 3u, scale.data());
	setGroup(FLUID_BLOCK_SCALE, "FluidBlock");
	setDescription(FLUID_BLOCK_SCALE, "Scaling vector of the fluid block.");

	FLUID_BLOCK_VISMESH = createStringParameter("visMesh", "Visualization mesh", &visMeshFile);
	setGroup(FLUID_BLOCK_VISMESH, "FluidBlock");
	setDescription(FLUID_BLOCK_VISMESH, "Path of an OBJ/PLY file containing a high resolution mesh which is used by the tool MeshSkinning to generate a sequence of deformed meshes (more info about this can be found in the documentation of the tool).");

	FLUID_BLOCK_MODE = createNumericParameter<unsigned char>("denseMode", "Dense mode", &mode);
	setGroup(FLUID_BLOCK_MODE, "FluidBlock");
	setDescription(FLUID_BLOCK_MODE, "Sampling mode: 0: regular sampling, 1: more dense sampling, 2 : dense sampling");

	FLUID_BLOCK_INITIAL_VEL = createVectorParameter("initialVelocity", "Initial velocity", 3u, initialVelocity.data());
	setGroup(FLUID_BLOCK_INITIAL_VEL, "FluidBlock");
	setDescription(FLUID_BLOCK_INITIAL_VEL, "The initial velocity is set for all particles in the block.");

	FLUID_BLOCK_INITIAL_ANGVEL = createVectorParameter("initialAngularVelocity", "Initial angular velocity", 3u, initialAngularVelocity.data());
	setGroup(FLUID_BLOCK_INITIAL_ANGVEL, "FluidBlock");
	setDescription(FLUID_BLOCK_INITIAL_ANGVEL, "The initial angular velocity of the block.");
}

//////////////////////////////////////////////////////////////////////////
// FluidModelParameterObject
//////////////////////////////////////////////////////////////////////////
int FluidModelParameterObject::FLUID_MODEL_ID = -1;
int FluidModelParameterObject::FLUID_MODEL_TRANSLATION = -1;
int FluidModelParameterObject::FLUID_MODEL_SCALE = -1;
int FluidModelParameterObject::FLUID_MODEL_VISMESH = -1;
int FluidModelParameterObject::FLUID_MODEL_MODE = -1;
int FluidModelParameterObject::FLUID_MODEL_INITIAL_VEL = -1;
int FluidModelParameterObject::FLUID_MODEL_INITIAL_ANGVEL = -1;
int FluidModelParameterObject::FLUID_MODEL_SAMPLES_FILE = -1;
int FluidModelParameterObject::FLUID_MODEL_ROTAXIS = -1;
int FluidModelParameterObject::FLUID_MODEL_ROTANGLE = -1;
int FluidModelParameterObject::FLUID_MODEL_INVERT = -1;
int FluidModelParameterObject::FLUID_MODEL_RESSDF = -1;

void FluidModelParameterObject::initParameters()
{
	FLUID_MODEL_ID = createStringParameter("id", "ID", &id);
	setGroup(FLUID_MODEL_ID, "FluidModel");
	setDescription(FLUID_MODEL_ID, "This id is used to define the properties of the fluid model. If no id is defined, then the standard id \"Fluid\" is used.");

	FLUID_MODEL_TRANSLATION = createVectorParameter("translation", "Translation", 3u, translation.data());
	setGroup(FLUID_MODEL_TRANSLATION, "FluidModel");
	setDescription(FLUID_MODEL_TRANSLATION, "Translation vector of the fluid model.");

	FLUID_MODEL_SCALE = createVectorParameter("scale", "Scale", 3u, scale.data());
	setGroup(FLUID_MODEL_SCALE, "FluidModel");
	setDescription(FLUID_MODEL_SCALE, "Scaling vector of the fluid model.");

	FLUID_MODEL_VISMESH = createStringParameter("visMesh", "Visualization mesh", &visMeshFile);
	setGroup(FLUID_MODEL_VISMESH, "FluidModel");
	setDescription(FLUID_MODEL_VISMESH, "Path of an OBJ/PLY file containing a high resolution mesh which is used by the tool MeshSkinning to generate a sequence of deformed meshes (more info about this can be found in the documentation of the tool).");

	FLUID_MODEL_MODE = createNumericParameter<unsigned char>("denseMode", "Dense mode", &mode);
	setGroup(FLUID_MODEL_MODE, "FluidModel");
	setDescription(FLUID_MODEL_MODE, "Sampling mode: 0: regular sampling, 1: more dense sampling, 2 : dense sampling");

	FLUID_MODEL_INITIAL_VEL = createVectorParameter("initialVelocity", "Initial velocity", 3u, initialVelocity.data());
	setGroup(FLUID_MODEL_INITIAL_VEL, "FluidModel");
	setDescription(FLUID_MODEL_INITIAL_VEL, "The initial velocity is set for all particles in the fluid model.");

	FLUID_MODEL_INITIAL_ANGVEL = createVectorParameter("initialAngularVelocity", "Initial angular velocity", 3u, initialAngularVelocity.data());
	setGroup(FLUID_MODEL_INITIAL_ANGVEL, "FluidModel");
	setDescription(FLUID_MODEL_INITIAL_ANGVEL, "The initial angular velocity of the fluid model.");

	FLUID_MODEL_SAMPLES_FILE = createStringParameter("particleFile", "Samples file", &samplesFile);
	setGroup(FLUID_MODEL_SAMPLES_FILE, "FluidModel");
	setDescription(FLUID_MODEL_SAMPLES_FILE, "Path of the partio file which contains the particle sampling or path of an OBJ/PLY file containing a closed mesh which is automatically sampled by SPlisHSPlasH. If you choose the latter option, you can define the sampling mode (denseMode), the resolution of the signed distance field which is used for the sampling (resolutionSDF) and if the signed distance field should be inverted (invert).");

	FLUID_MODEL_ROTAXIS = createVectorParameter("rotationAxis", "Rotation axis", 3u, axis.data());
	setGroup(FLUID_MODEL_ROTAXIS, "FluidModel");
	setDescription(FLUID_MODEL_ROTAXIS, "Axis used to rotate the particle data after loading.");

	FLUID_MODEL_ROTANGLE = createNumericParameter<Real>("rotationAngle", "Rotation angle", &angle);
	setGroup(FLUID_MODEL_ROTANGLE, "FluidModel");
	setDescription(FLUID_MODEL_ROTANGLE, "Rotation angle for the initial rotation of the particle data.");

	FLUID_MODEL_INVERT = createBoolParameter("invert", "Invert geomtry", &invert);
	setGroup(FLUID_MODEL_INVERT, "FluidModel");
	setDescription(FLUID_MODEL_INVERT, "Invert the signed distance field, flips inside/outside.");

	FLUID_MODEL_RESSDF = createVectorParameter("resolutionSDF", "SDF resolution", 3u, resolutionSDF.data());
	setGroup(FLUID_MODEL_RESSDF, "FluidModel");
	setDescription(FLUID_MODEL_RESSDF, "Resolution of the signed distance field.");
}


//////////////////////////////////////////////////////////////////////////
// EmitterParameterObject
//////////////////////////////////////////////////////////////////////////
int EmitterParameterObject::EMITTER_ID = -1;
int EmitterParameterObject::EMITTER_WIDTH = -1;
int EmitterParameterObject::EMITTER_HEIGHT = -1;
int EmitterParameterObject::EMITTER_POSITION = -1;
int EmitterParameterObject::EMITTER_VELOCITY = -1;
int EmitterParameterObject::EMITTER_ROTAXIS = -1;
int EmitterParameterObject::EMITTER_ROTANGLE = -1;
int EmitterParameterObject::EMITTER_STARTTIME = -1;
int EmitterParameterObject::EMITTER_ENDTIME = -1;
int EmitterParameterObject::EMITTER_TYPE = -1;

void EmitterParameterObject::initParameters()
{
	EMITTER_ID = createStringParameter("id", "ID", &id);
	setGroup(EMITTER_ID, "Emitter");
	setDescription(EMITTER_ID, "This id is used to define the properties of the fluid model. If no id is defined, then the standard id \"Fluid\" is used.");

	EMITTER_WIDTH = createNumericParameter<unsigned int>("width", "Width", &width);
	setGroup(EMITTER_WIDTH, "Emitter");
	setDescription(EMITTER_WIDTH, "Width of the box or radius of the circle emitter in number of particles.");

	EMITTER_HEIGHT = createNumericParameter<unsigned int>("height", "Height", &height);
	setGroup(EMITTER_HEIGHT, "Emitter");
	setDescription(EMITTER_HEIGHT, "Height of the box in number of particles (is only used for type 0).");

	EMITTER_POSITION = createVectorParameter("translation", "Translation", 3u, x.data());
	setGroup(EMITTER_POSITION, "Emitter");
	setDescription(EMITTER_POSITION, "Translation vector of the emitter.");

	EMITTER_VELOCITY = createNumericParameter<Real>("velocity", "Velocity", &velocity);
	setGroup(EMITTER_VELOCITY, "Emitter");
	setDescription(EMITTER_VELOCITY, "Initial velocity of the emitted particles in direction of the emitter.");

	EMITTER_ROTAXIS = createVectorParameter("rotationAxis", "Rotation axis", 3u, axis.data());
	setGroup(EMITTER_ROTAXIS, "Emitter");
	setDescription(EMITTER_ROTAXIS, "Axis used to rotate the emitter. Note that in 2D simulations the axis is always set to [0,0,1].");

	EMITTER_ROTANGLE = createNumericParameter<Real>("rotationAngle", "Rotation angle", &angle);
	setGroup(EMITTER_ROTANGLE, "Emitter");
	setDescription(EMITTER_ROTANGLE, "Rotation angle for the initial rotation of the emitter.");

	EMITTER_STARTTIME = createNumericParameter<Real>("emitStartTime", "Start time", &emitStartTime);
	setGroup(EMITTER_STARTTIME, "Emitter");
	setDescription(EMITTER_STARTTIME, "Start time of the emitter.");

	EMITTER_ENDTIME = createNumericParameter<Real>("emitEndTime", "End time", &emitEndTime);
	setGroup(EMITTER_ENDTIME, "Emitter");
	setDescription(EMITTER_ENDTIME, "End time of the emitter.");

	EMITTER_TYPE = createNumericParameter<unsigned int>("type", "Emitter type", &type);
	setGroup(EMITTER_TYPE, "Emitter");
	setDescription(EMITTER_TYPE, "Defines the shape of the emitter: 0: box, 1: circle.");
}


//////////////////////////////////////////////////////////////////////////
// AnimationFieldParameterObject
//////////////////////////////////////////////////////////////////////////
int AnimationFieldParameterObject::ANIMATIONFIELD_PARTICLE_FIELD = -1;
int AnimationFieldParameterObject::ANIMATIONFIELD_EXPR_X = -1;
int AnimationFieldParameterObject::ANIMATIONFIELD_EXPR_Y = -1;
int AnimationFieldParameterObject::ANIMATIONFIELD_EXPR_Z = -1;
int AnimationFieldParameterObject::ANIMATIONFIELD_SHAPETYPE = -1;
int AnimationFieldParameterObject::ANIMATIONFIELD_TRANSLATION = -1;
int AnimationFieldParameterObject::ANIMATIONFIELD_ROTAXIS = -1;
int AnimationFieldParameterObject::ANIMATIONFIELD_ROTANGLE = -1;
int AnimationFieldParameterObject::ANIMATIONFIELD_SCALE = -1;
int AnimationFieldParameterObject::ANIMATIONFIELD_STARTTIME = -1;
int AnimationFieldParameterObject::ANIMATIONFIELD_ENDTIME = -1;


void AnimationFieldParameterObject::initParameters()
{
	ANIMATIONFIELD_PARTICLE_FIELD = createStringParameter("particleField", "Particle field", &particleFieldName);
	setGroup(ANIMATIONFIELD_PARTICLE_FIELD, "Animation field");
	setDescription(ANIMATIONFIELD_PARTICLE_FIELD, "Defines the field quantity that should be modified by the field (e.g. velocity, angular velocity, position).");
	
	ANIMATIONFIELD_EXPR_X = createStringParameter("expression_x", "Expression X", &expression[0]);
	setGroup(ANIMATIONFIELD_EXPR_X, "Animation field");
	setDescription(ANIMATIONFIELD_EXPR_X, "Math expression for the x-component of the field quantity.");

	ANIMATIONFIELD_EXPR_Y = createStringParameter("expression_y", "Expression Y", &expression[1]);
	setGroup(ANIMATIONFIELD_EXPR_Y, "Animation field");
	setDescription(ANIMATIONFIELD_EXPR_Y, "Math expression for the y-component of the field quantity.");
	
	ANIMATIONFIELD_EXPR_Z = createStringParameter("expression_z", "Expression Z", &expression[2]);
	setGroup(ANIMATIONFIELD_EXPR_Z, "Animation field");
	setDescription(ANIMATIONFIELD_EXPR_Z, "Math expression for the z-component of the field quantity.");

	ANIMATIONFIELD_SHAPETYPE = createNumericParameter<unsigned int>("shapeType", "Shape type", &shapeType);
	setGroup(ANIMATIONFIELD_SHAPETYPE, "Animation field");
	setDescription(ANIMATIONFIELD_SHAPETYPE, "Defines the shape of the animation field: 0: box, 1 : sphere, 2 : cylinder.");

	ANIMATIONFIELD_TRANSLATION = createVectorParameter("translation", "Translation", 3u, translation.data());
	setGroup(ANIMATIONFIELD_TRANSLATION, "Animation field");
	setDescription(ANIMATIONFIELD_TRANSLATION, "Translation vector of the animation field.");

	ANIMATIONFIELD_ROTAXIS = createVectorParameter("rotationAxis", "Rotation axis", 3u, axis.data());
	setGroup(ANIMATIONFIELD_ROTAXIS, "Animation field");
	setDescription(ANIMATIONFIELD_ROTAXIS, "Axis used to rotate the animation field .");

	ANIMATIONFIELD_ROTANGLE = createNumericParameter<Real>("rotationAngle", "Rotation angle", &angle);
	setGroup(ANIMATIONFIELD_ROTANGLE, "Animation field");
	setDescription(ANIMATIONFIELD_ROTANGLE, "Rotation angle for the initial rotation of the animation field.");

	ANIMATIONFIELD_SCALE = createVectorParameter("scale", "Scale", 3u, scale.data());
	setGroup(ANIMATIONFIELD_SCALE, "Animation field");
	setDescription(ANIMATIONFIELD_SCALE, "Scaling vector of the animation field (shapeType=0 (box): This vector defines the width, height, depth of the box. shapeType = 1 (sphere) : The x - component of the vector defines the radius of the sphere.The other components are ignored. shapeType = 2 (cylinder) : The x - and y - component of the vector defines the height and radius of the cylinder, repectively.The z - component is ignored.)");

	ANIMATIONFIELD_STARTTIME = createNumericParameter<Real>("startTime", "Start time", &startTime);
	setGroup(ANIMATIONFIELD_STARTTIME, "Animation field");
	setDescription(ANIMATIONFIELD_STARTTIME, "Start time of the animation field.");

	ANIMATIONFIELD_ENDTIME = createNumericParameter<Real>("endTime", "End time", &endTime);
	setGroup(ANIMATIONFIELD_ENDTIME, "Animation field");
	setDescription(ANIMATIONFIELD_ENDTIME, "End time of the animation field.");	
}

//////////////////////////////////////////////////////////////////////////
// MaterialParameterObject
//////////////////////////////////////////////////////////////////////////

int MaterialParameterObject::MATERIAL_ID = -1;
int MaterialParameterObject::MATERIAL_MIN_VAL = -1;
int MaterialParameterObject::MATERIAL_MAX_VAL = -1;
int MaterialParameterObject::MATERIAL_COLOR_FIELD = -1;
int MaterialParameterObject::MATERIAL_COLOR_MAP = -1;
int MaterialParameterObject::MATERIAL_MAX_EMITTER_PARTICLES = -1;
int MaterialParameterObject::MATERIAL_EMITTER_REUSE = -1;
int MaterialParameterObject::MATERIAL_EMITTER_BOX_MIN = -1;
int MaterialParameterObject::MATERIAL_EMITTER_BOX_MAX = -1;

void MaterialParameterObject::initParameters()
{
	MATERIAL_ID = createStringParameter("id", "ID", &id);
	setGroup(MATERIAL_ID, "Material");
	setDescription(MATERIAL_ID, "Defines the id of the material. You have to give the same id to a FluidBlock, a FluidModel or an Emitter if they should have the defined material behavior.");

	MATERIAL_MIN_VAL = createNumericParameter<Real>("renderMinValue", "Min. value", &minVal);
	setGroup(MATERIAL_MIN_VAL, "Material");
	setDescription(MATERIAL_MIN_VAL, "Minimum value used for color-coding the color field in the rendering process.");

	MATERIAL_MAX_VAL = createNumericParameter<Real>("renderMaxValue", "Max. value", &maxVal);
	setGroup(MATERIAL_MAX_VAL, "Material");
	setDescription(MATERIAL_MAX_VAL, "Maximum value used for color-coding the color field in the rendering process.");

	MATERIAL_COLOR_FIELD = createStringParameter("colorField", "Color field", &colorField);
	setGroup(MATERIAL_COLOR_FIELD, "Material");
	setDescription(MATERIAL_COLOR_FIELD, "Choose vector or scalar field for particle coloring.");

	MATERIAL_COLOR_MAP = createNumericParameter<unsigned int>("colorMapType", "Color map type", &colorMapType);
	setGroup(MATERIAL_COLOR_MAP, "Material");
	setDescription(MATERIAL_COLOR_MAP, "Selection of a color map for coloring the scalar/vector field: 0: None, 1 : Jet, 2 : Plasma, 3 : CoolWarm, 4 : BlueWhiteRed, 5 : Seismic");

	MATERIAL_MAX_EMITTER_PARTICLES = createNumericParameter<unsigned int>("maxEmitterParticles", "Max. emitter particles", &maxEmitterParticles);
	setGroup(MATERIAL_MAX_EMITTER_PARTICLES, "Material");
	setDescription(MATERIAL_MAX_EMITTER_PARTICLES, "Maximum number of particles the emitter generates. Note that reused particles are not counted here.");

	MATERIAL_EMITTER_REUSE = createBoolParameter("emitterReuseParticles", "Reuse particles (emitter)", &emitterReuseParticles);
	setGroup(MATERIAL_EMITTER_REUSE, "Material");
	setDescription(MATERIAL_EMITTER_REUSE, "Reuse particles if they are outside of the bounding box defined by emitterBoxMin, emitterBoxMax.");

	MATERIAL_EMITTER_BOX_MIN = createVectorParameter("emitterBoxMin", "Emitter box min.", 3u, emitterBoxMin.data());
	setGroup(MATERIAL_EMITTER_BOX_MIN, "Material");
	setDescription(MATERIAL_EMITTER_BOX_MIN, "Minimum coordinates of an axis-aligned box (used in combination with emitterReuseParticles).");

	MATERIAL_EMITTER_BOX_MAX = createVectorParameter("emitterBoxMax", "Emitter box max.", 3u, emitterBoxMax.data());
	setGroup(MATERIAL_EMITTER_BOX_MAX, "Material");
	setDescription(MATERIAL_EMITTER_BOX_MAX, "Maximum coordinates of an axis-aligned box (used in combination with emitterReuseParticles).");
}

//////////////////////////////////////////////////////////////////////////
// BoundaryParameterObject
//////////////////////////////////////////////////////////////////////////

int BoundaryParameterObject::BOUNDARY_SAMPLES_FILE = -1;
int BoundaryParameterObject::BOUNDARY_MESH_FILE = -1;
int BoundaryParameterObject::BOUNDARY_TRANSLATION = -1;
int BoundaryParameterObject::BOUNDARY_AXIS = -1;
int BoundaryParameterObject::BOUNDARY_ANGLE = -1;
int BoundaryParameterObject::BOUNDARY_SCALE = -1;
int BoundaryParameterObject::BOUNDARY_DYNAMIC = -1;
int BoundaryParameterObject::BOUNDARY_IS_WALL = -1;
int BoundaryParameterObject::BOUNDARY_COLOR = -1;
int BoundaryParameterObject::BOUNDARY_MAP_FILE = -1;
int BoundaryParameterObject::BOUNDARY_MAP_INVERT = -1;
int BoundaryParameterObject::BOUNDARY_MAP_THICKNESS = -1;
int BoundaryParameterObject::BOUNDARY_MAP_RESOLUTION = -1;
int BoundaryParameterObject::BOUNDARY_SAMPLING_MODE = -1;
int BoundaryParameterObject::BOUNDARY_IS_ANIMATED = -1;

void BoundaryParameterObject::initParameters()
{
	BOUNDARY_SAMPLES_FILE = createStringParameter("particleFile", "Particle file", &samplesFile);
	setGroup(BOUNDARY_SAMPLES_FILE, "Boundary");
	setDescription(BOUNDARY_SAMPLES_FILE, "Path to a partio file which contains a surface sampling of the body. Note that the surface sampling is done automatically if this parameter is missing.");

	BOUNDARY_MESH_FILE = createStringParameter("geometryFile", "Geometry file", &meshFile);
	setGroup(BOUNDARY_MESH_FILE, "Boundary");
	setDescription(BOUNDARY_MESH_FILE, "Path to a OBJ/PLY file which contains the geometry of the body.");

	BOUNDARY_TRANSLATION = createVectorParameter("translation", "Translation", 3u, translation.data());
	setGroup(BOUNDARY_TRANSLATION, "Boundary");
	setDescription(BOUNDARY_TRANSLATION, "Translation vector of the rigid body.");

	BOUNDARY_AXIS = createVectorParameter("rotationAxis", "Rotation axis", 3u, axis.data());
	setGroup(BOUNDARY_AXIS, "Boundary");
	setDescription(BOUNDARY_AXIS, "Axis used to rotate the rigid body after loading.");

	BOUNDARY_ANGLE = createNumericParameter<Real>("rotationAngle", "Rotation angle", &angle);
	setGroup(BOUNDARY_ANGLE, "Boundary");
	setDescription(BOUNDARY_ANGLE, "Rotation angle for the initial rotation of the rigid body.");

	BOUNDARY_SCALE = createVectorParameter("scale", "Scale", 3u, scale.data());
	setGroup(BOUNDARY_SCALE, "Boundary");
	setDescription(BOUNDARY_SCALE, "Scaling vector of the rigid body.");

	BOUNDARY_DYNAMIC = createBoolParameter("isDynamic", "Dynamic", &dynamic);
	setGroup(BOUNDARY_DYNAMIC, "Boundary");
	setDescription(BOUNDARY_DYNAMIC, "Defines if the body is static or dynamic.");

	BOUNDARY_IS_WALL = createBoolParameter("isWall", "Wall", &isWall);
	setGroup(BOUNDARY_IS_WALL, "Boundary");
	setDescription(BOUNDARY_IS_WALL, "Defines if this is a wall. Walls are typically not rendered. This is the only difference.");

	BOUNDARY_COLOR = createVectorParameter("color", "Color", 4u, color.data());
	setGroup(BOUNDARY_COLOR, "Boundary");
	setDescription(BOUNDARY_COLOR, "RGBA color of the body.");

	BOUNDARY_MAP_FILE = createStringParameter("mapFile", "Map file", &mapFile);
	setGroup(BOUNDARY_MAP_FILE, "Boundary");
	setDescription(BOUNDARY_MAP_FILE, "Path to a volume/density map file which contains a volume/density map of the body. Note that the map is generated automatically if this parameter is missing.");

	BOUNDARY_MAP_INVERT = createBoolParameter("mapInvert", "Invert map", &mapInvert);
	setGroup(BOUNDARY_MAP_INVERT, "Boundary");
	setDescription(BOUNDARY_MAP_INVERT, "Invert the map when using density or volume maps, flips inside/outside.");

	BOUNDARY_MAP_THICKNESS = createNumericParameter<Real>("mapThickness", "Map thickness", &mapThickness);
	setGroup(BOUNDARY_MAP_THICKNESS, "Boundary");
	setDescription(BOUNDARY_MAP_THICKNESS, "Additional thickness of a volume or density map.");

	BOUNDARY_MAP_RESOLUTION = createVectorParameter("mapResolution", "Map resolution", 3u, mapResolution.data());
	setGroup(BOUNDARY_MAP_RESOLUTION, "Boundary");
	setDescription(BOUNDARY_MAP_RESOLUTION, " Resolution of a volume or density map.");

	BOUNDARY_SAMPLING_MODE = createNumericParameter<unsigned int>("samplingMode", "Sampling mode", &samplingMode);
	setGroup(BOUNDARY_SAMPLING_MODE, "Boundary");
	setDescription(BOUNDARY_SAMPLING_MODE, "Surface sampling mode. 0 Poisson disk sampling, 1 Regular triangle sampling.");

	BOUNDARY_IS_ANIMATED = createBoolParameter("isAnimated", "Animated", &isAnimated);
	setGroup(BOUNDARY_IS_ANIMATED, "Boundary");
	setDescription(BOUNDARY_IS_ANIMATED, "Defines if the body is animated (e.g. by a script).");	
}