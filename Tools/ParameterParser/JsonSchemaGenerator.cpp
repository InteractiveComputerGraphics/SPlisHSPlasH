#include "JsonSchemaGenerator.h"
#include "Utilities/StringTools.h"
#include "Simulator/SimulatorBase.h"
#include "SPlisHSPlasH/Simulation.h"
#include "ParameterObjectParser.h"

using namespace GenParam;
using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;


void JsonSchemaGenerator::jsonBoolParam(BoolParameter* param)
{
	(*m_currentData)[param->getName()]["description"] = param->getDescription() + " (Default: " + std::to_string(param->getValue()) + ", Type: " + typeid(param->getValue()).name() + ")";
	(*m_currentData)[param->getName()]["type"] = "boolean";
}

void JsonSchemaGenerator::jsonStringParam(StringParameter* param)
{
	(*m_currentData)[param->getName()]["description"] = param->getDescription() + " (Default: " + param->getValue() + ", Type: string)";
	(*m_currentData)[param->getName()]["type"] = "string";
}

void JsonSchemaGenerator::jsonEnumParam(EnumParameter* param)
{
	auto& vals = param->getEnumValues();
	std::string str = param->getDescription() + " (Default: " + std::to_string(param->getValue()) + ", Type: enum)";
	str = str + ", Options:";
	for (size_t i = 0; i < vals.size(); i++)
		str = str + " (" + std::to_string(vals[i].id) + ", " + vals[i].name + ")";

	(*m_currentData)[param->getName()]["description"] = str;
	(*m_currentData)[param->getName()]["type"] = "integer";
}

void JsonSchemaGenerator::jsonVecUintParam(VectorParameter<unsigned int>* param)
{
	if (param->getDim() == 3u)
	{
		Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign> vec(param->getValue());
		(*m_currentData)[param->getName()]["description"] = param->getDescription() + " (Default: (" + std::to_string(vec[0]) + "," + std::to_string(vec[1]) + "," + std::to_string(vec[2]) + "), Type: vec3uint)";
		(*m_currentData)[param->getName()]["$ref"] = "#/definitions/vec3int";
	}
}

void JsonSchemaGenerator::jsonVecParam(RealVectorParameter* param)
{
	if (param->getDim() == 3u)
	{
		Vector3r vec(param->getValue());
		(*m_currentData)[param->getName()]["description"] = param->getDescription() + " (Default: (" + std::to_string(vec[0]) + "," + std::to_string(vec[1]) + "," + std::to_string(vec[2]) + "), Type: vec3real)";
		(*m_currentData)[param->getName()]["$ref"] = "#/definitions/vec3real";
	}
	else if (param->getDim() == 4u)
	{
		Vector4r vec(param->getValue());
		(*m_currentData)[param->getName()]["description"] = param->getDescription() + " (Default: (" + StringTools::real2String(vec[0]) + "," + StringTools::real2String(vec[1]) + "," + StringTools::real2String(vec[2]) + "," + StringTools::real2String(vec[3]) + "), Type: vec4float)";
		(*m_currentData)[param->getName()]["$ref"] = "#/definitions/vec4real";
	}
}

void JsonSchemaGenerator::jsonParameterObject(const ParameterObject* obj)
{
	if ((dynamic_cast<const SimulatorBase*>(obj) != nullptr) ||
		(dynamic_cast<const Simulation*>(obj) != nullptr))
		m_currentData = &m_jsonData["properties"]["Configuration"]["properties"];
	else if ((dynamic_cast<const NonPressureForceBase*>(obj) != nullptr) ||
		(dynamic_cast<const MaterialParameterObject*>(obj) != nullptr) ||
		(dynamic_cast<const FluidModel*>(obj) != nullptr))
	{
		m_jsonData["definitions"]["Material"]["type"] = "object";
		m_currentData = &m_jsonData["definitions"]["Material"]["properties"];
	}
	else if (dynamic_cast<const FluidBlockParameterObject*>(obj) != nullptr)
	{
		m_jsonData["definitions"]["FluidBlock"]["type"] = "object";
		m_currentData = &m_jsonData["definitions"]["FluidBlock"]["properties"];
	}
	else if (dynamic_cast<const FluidModelParameterObject*>(obj) != nullptr)
	{
		m_jsonData["definitions"]["FluidModel"]["type"] = "object";
		m_currentData = &m_jsonData["definitions"]["FluidModel"]["properties"];
	}
	else if (dynamic_cast<const EmitterParameterObject*>(obj) != nullptr)
	{
		m_jsonData["definitions"]["Emitter"]["type"] = "object";
		m_currentData = &m_jsonData["definitions"]["Emitter"]["properties"];
	}
	else if (dynamic_cast<const BoundaryParameterObject*>(obj) != nullptr)
	{
		m_jsonData["definitions"]["RigidBody"]["type"] = "object";
		m_currentData = &m_jsonData["definitions"]["RigidBody"]["properties"];
	}
	else if (dynamic_cast<const AnimationFieldParameterObject*>(obj) != nullptr)
	{
		m_jsonData["definitions"]["AnimationField"]["type"] = "object";
		m_currentData = &m_jsonData["definitions"]["AnimationField"]["properties"];
	}
}

void JsonSchemaGenerator::generateSchemaFile(const std::string& fileName)
{
	std::cout << "Write schema file: " << fileName << "\n";

	//////////////////////////////////////////////////////////////////////////
	// update configuration 
	//////////////////////////////////////////////////////////////////////////
	m_jsonData["$schema"] = "http://json-schema.org/schema#";
	m_jsonData["title"] = "SPlisHSPlasH File Format";
	m_jsonData["description"] = "A json schema for SPlisHSPlasH's scene files";

	m_jsonData["type"] = "object";

	nlohmann::json& config = m_jsonData["properties"]["Configuration"];
	config["type"] = "object";
	config["description"] = "Contains the general settings of the simulation and the pressure solver";

	// Basic type definitions
	nlohmann::json& defs = m_jsonData["definitions"];
	nlohmann::json& def_vec3real = defs["vec3real"];
	def_vec3real["type"] = "array";
	def_vec3real["minItems"] = 3;
	def_vec3real["maxItems"] = 3;
	nlohmann::json& def_vec3real_items = def_vec3real["items"];
	def_vec3real_items["type"] = "number";

	nlohmann::json& def_vec3int = defs["vec3int"];
	def_vec3int["type"] = "array";
	def_vec3int["minItems"] = 3;
	def_vec3int["maxItems"] = 3;
	nlohmann::json& def_vec3uint_items = def_vec3int["items"];
	def_vec3uint_items["type"] = "integer";

	nlohmann::json& def_vec4real = defs["vec4real"];
	def_vec4real["type"] = "array";
	def_vec4real["minItems"] = 4;
	def_vec4real["maxItems"] = 4;
	nlohmann::json& def_vec4real_items = def_vec4real["items"];
	def_vec4real_items["type"] = "number";
	def_vec4real_items["minItems"] = 0;
	def_vec4real_items["maxItems"] = 1;


	ParameterObjectParser parser;
	parser.addRealParamCB(std::bind(&JsonSchemaGenerator::jsonNumericParam<NumericParameter<Real>>, this, std::placeholders::_1));
	parser.addUInt32ParamCB(std::bind(&JsonSchemaGenerator::jsonNumericParam<NumericParameter<unsigned int>>, this, std::placeholders::_1));
	parser.addUInt16ParamCB(std::bind(&JsonSchemaGenerator::jsonNumericParam<NumericParameter<unsigned short>>, this, std::placeholders::_1));
	parser.addUInt8ParamCB(std::bind(&JsonSchemaGenerator::jsonNumericParam<NumericParameter<unsigned char>>, this, std::placeholders::_1));
	parser.addInt32ParamCB(std::bind(&JsonSchemaGenerator::jsonNumericParam<NumericParameter<int>>, this, std::placeholders::_1));
	parser.addInt16ParamCB(std::bind(&JsonSchemaGenerator::jsonNumericParam<NumericParameter<short>>, this, std::placeholders::_1));
	parser.addInt8ParamCB(std::bind(&JsonSchemaGenerator::jsonNumericParam<NumericParameter<char>>, this, std::placeholders::_1));
	parser.addBoolParamCB(std::bind(&JsonSchemaGenerator::jsonBoolParam, this, std::placeholders::_1));
	parser.addEnumParamCB(std::bind(&JsonSchemaGenerator::jsonEnumParam, this, std::placeholders::_1));
	parser.addStringParamCB(std::bind(&JsonSchemaGenerator::jsonStringParam, this, std::placeholders::_1));
	parser.addVecRealParamCB(std::bind(&JsonSchemaGenerator::jsonVecParam, this, std::placeholders::_1));
	parser.addVecUintParamCB(std::bind(&JsonSchemaGenerator::jsonVecUintParam, this, std::placeholders::_1));
	parser.addParameterObjectCB(std::bind(&JsonSchemaGenerator::jsonParameterObject, this, std::placeholders::_1));

	parser.parseParameters();

	// fluid blocks properties
	nlohmann::json& props = m_jsonData["properties"];
	nlohmann::json& blocks = props["FluidBlocks"];
	blocks["type"] = "array";
	blocks["description"] = "Definition of axis-aligned blocks of fluid particles.";
	nlohmann::json& blocks_items = blocks["items"];
	blocks_items["$ref"] = "#/definitions/FluidBlock";

	// fluid models properties
	nlohmann::json& fmodels = props["FluidModels"];
	fmodels["type"] = "array";
	fmodels["description"] = "Definition of fluid models using a sampled geometry.";
	nlohmann::json& fmodels_items = fmodels["items"];
	fmodels_items["$ref"] = "#/definitions/FluidModel";

	// emitters properties
	nlohmann::json& emitters = props["Emitters"];
	emitters["type"] = "array";
	emitters["description"] = "Definition of fluid emitters.";
	nlohmann::json& emitters_items = emitters["items"];
	emitters_items["$ref"] = "#/definitions/Emitter";

	// rigid bodies properties
	nlohmann::json& rbs = props["RigidBodies"];
	rbs["type"] = "array";
	rbs["description"] = "Definition of static and dynamic rigid bodies. The latter are handled by PositionBasedDynamic, thus you may want to take a look at its json format as well.";
	nlohmann::json& rbs_items = rbs["items"];
	rbs_items["$ref"] = "#/definitions/RigidBody";

	// material properties
	nlohmann::json& mat = props["Materials"];
	mat["type"] = "array";
	mat["description"] = "Describes the material properties of the associated fluid models.";
	nlohmann::json& mat_items = mat["items"];
	mat_items["$ref"] = "#/definitions/Material";

	// animation fields properties
	nlohmann::json& afs = props["AnimationFields"];
	afs["type"] = "array";
	afs["description"] = "Definition of animation fields with the following terms: typical math terms: cos, sin; current time: t; Positions: x, y, z; Velocities: vx, vy,vz; Field quantities: valuex, valuey, valuez.";
	nlohmann::json& afs_items = afs["items"];
	afs_items["$ref"] = "#/definitions/AnimationField";

	std::ofstream output_file(fileName);
	if (!output_file.is_open())
	{
		LOG_ERR << "Cannot open file!";
		return;
	}
	try
	{
		output_file << std::setw(4) << m_jsonData << std::endl;
	}
	catch (const std::exception& e)
	{
		LOG_ERR << e.what();
		exit(1);
	}
}