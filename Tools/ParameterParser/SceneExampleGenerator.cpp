#include "SceneExampleGenerator.h"
#include "Utilities/StringTools.h"
#include "Simulator/SimulatorBase.h"
#include "SPlisHSPlasH/Simulation.h"
#include "ParameterObjectParser.h"

using namespace GenParam;
using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;


void SceneExampleGenerator::jsonEnumParam(EnumParameter* param)
{
	auto& vals = param->getEnumValues();
	std::string str1 = param->getName() + " - comment";
	std::string str2 = param->getDescription() + "(Default: " + std::to_string(param->getValue()) + ", Type: enum)";
	str2 = str2 + ", Options:";
	for (size_t i = 0; i < vals.size(); i++)
		str2 = str2 + " (" + std::to_string(vals[i].id) + ", " + vals[i].name + ")";

	(*m_currentData)[str1] = str2;
	(*m_currentData)[param->getName()] = param->getValue();
}

void SceneExampleGenerator::jsonStringParam(StringParameter* param)
{
	std::string str1 = param->getName() + " - comment";
	std::string str2 = param->getDescription() + "(Default: " + param->getValue() + ", Type: string)";
	(*m_currentData)[str1] = str2;
	(*m_currentData)[param->getName()] = param->getValue();
}

void SceneExampleGenerator::jsonVecParam(RealVectorParameter* param)
{
	std::string str1 = param->getName() + " - comment";
	std::string str2 = param->getDescription();
	if (param->getDim() == 3u)
	{
		Vector3r vec(param->getValue());
		str2 = str2 + "(Default: (" + StringTools::real2String(vec[0]) + "," + StringTools::real2String(vec[1]) + "," + StringTools::real2String(vec[2]) + "), Type: vec3float)";
		// convert to std vector 
		std::vector<Real> stdVec{ vec.data(), vec.data() + vec.size() };
		(*m_currentData)[param->getName()] = stdVec;
		(*m_currentData)[str1] = str2;
	}
	else if (param->getDim() == 4u)
	{
		Vector4r vec(param->getValue());
		str2 = str2 + "(Default: (" + StringTools::real2String(vec[0]) + "," + StringTools::real2String(vec[1]) + "," + StringTools::real2String(vec[2]) + "," + StringTools::real2String(vec[3]) + "), Type: vec4float)";
		// convert to std vector 
		std::vector<Real> stdVec{ vec.data(), vec.data() + vec.size() };
		(*m_currentData)[param->getName()] = stdVec;
		(*m_currentData)[str1] = str2;
	}
}

void SceneExampleGenerator::jsonVecUintParam(VectorParameter<unsigned int>* param)
{
	std::string str1 = param->getName() + " - comment";
	if (param->getDim() == 3u)
	{
		Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign> vec(param->getValue());
		std::string str2 = param->getDescription() + "(Default: (" + std::to_string(vec[0]) + "," + std::to_string(vec[1]) + "," + std::to_string(vec[2]) + "), Type: vec3uint)";
		(*m_currentData)[str1] = str2;
		// convert to std vector 
		std::vector<unsigned int> stdVec{ vec.data(), vec.data() + vec.size() };
		(*m_currentData)[param->getName()] = stdVec;
	}
}

void SceneExampleGenerator::jsonParameterObject(const ParameterObject* obj)
{
	if ((dynamic_cast<const SimulatorBase*>(obj) != nullptr) ||
		(dynamic_cast<const Simulation*>(obj) != nullptr))
		m_currentData = &m_jsonData["Configuration"];
	else if ((dynamic_cast<const NonPressureForceBase*>(obj) != nullptr) ||
		(dynamic_cast<const MaterialParameterObject*>(obj) != nullptr) ||
		(dynamic_cast<const FluidModel*>(obj) != nullptr))
	{
		m_currentData = &m_jsonData["Materials"];
		m_currentData = &(*m_currentData)["Fluid"];
	}
	else if (dynamic_cast<const FluidBlockParameterObject*>(obj) != nullptr)
		m_currentData = &m_jsonData["FluidBlocks"];
	else if (dynamic_cast<const FluidModelParameterObject*>(obj) != nullptr)
		m_currentData = &m_jsonData["FluidModels"];
	else if (dynamic_cast<const EmitterParameterObject*>(obj) != nullptr)
		m_currentData = &m_jsonData["Emitters"];
	else if (dynamic_cast<const BoundaryParameterObject*>(obj) != nullptr)
		m_currentData = &m_jsonData["RigidBodies"];
	else if (dynamic_cast<const AnimationFieldParameterObject*>(obj) != nullptr)
		m_currentData = &m_jsonData["AnimationFields"];
}

void SceneExampleGenerator::generateExampleSceneFile(const std::string& fileName)
{
	std::cout << "Write scene file: " << fileName << "\n";

	//////////////////////////////////////////////////////////////////////////
	// update configuration 
	//////////////////////////////////////////////////////////////////////////
	ParameterObjectParser parser;
	parser.addRealParamCB(std::bind(&SceneExampleGenerator::jsonParam<NumericParameter<Real>>, this, std::placeholders::_1));
	parser.addUInt32ParamCB(std::bind(&SceneExampleGenerator::jsonParam<NumericParameter<unsigned int>>, this, std::placeholders::_1));
	parser.addUInt16ParamCB(std::bind(&SceneExampleGenerator::jsonParam<NumericParameter<unsigned short>>, this, std::placeholders::_1));
	parser.addUInt8ParamCB(std::bind(&SceneExampleGenerator::jsonParam<NumericParameter<unsigned char>>, this, std::placeholders::_1));
	parser.addInt32ParamCB(std::bind(&SceneExampleGenerator::jsonParam<NumericParameter<int>>, this, std::placeholders::_1));
	parser.addInt16ParamCB(std::bind(&SceneExampleGenerator::jsonParam<NumericParameter<short>>, this, std::placeholders::_1));
	parser.addInt8ParamCB(std::bind(&SceneExampleGenerator::jsonParam<NumericParameter<char>>, this, std::placeholders::_1));
	parser.addBoolParamCB(std::bind(&SceneExampleGenerator::jsonParam<BoolParameter>, this, std::placeholders::_1));
	parser.addEnumParamCB(std::bind(&SceneExampleGenerator::jsonEnumParam, this, std::placeholders::_1));
	parser.addStringParamCB(std::bind(&SceneExampleGenerator::jsonStringParam, this, std::placeholders::_1));
	parser.addVecRealParamCB(std::bind(&SceneExampleGenerator::jsonVecParam, this, std::placeholders::_1));
	parser.addVecUintParamCB(std::bind(&SceneExampleGenerator::jsonVecUintParam, this, std::placeholders::_1));
	parser.addParameterObjectCB(std::bind(&SceneExampleGenerator::jsonParameterObject, this, std::placeholders::_1));

	parser.parseParameters();

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
