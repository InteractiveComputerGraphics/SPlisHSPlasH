#include "ScriptObject.h"
#include "pySPlisHSPlasH/Embedded.h"
#include "Utilities/FileSystem.h"
#include "Utilities/Logger.h"
#include "SceneConfiguration.h"
#include "SimulatorBase.h"
#include "Simulator/GUI/Simulator_GUI_Base.h"

using namespace SPH;
using namespace std;
using namespace GenParam;
using namespace Utilities;

int ScriptObject::SCRIPT_FILE = -1;

ScriptObject::ScriptObject(SimulatorBase *base) 
{
	m_base = base;
	m_scriptLoaded = false;
	m_scriptModule = "";
	m_scriptFile = "";
}

ScriptObject::~ScriptObject(void)
{
	if (m_scriptLoaded)
	{
		Embedded::getCurrent()->releaseModule(m_scriptModule);
	}
}

void ScriptObject::initParameters()
{
	ParameterObject::initParameters();

	auto tmp = createFunctionParameter("reload", "Reload script", [&]() { 
		if (m_scriptModule != "")
			Embedded::getCurrent()->reloadModule(m_scriptModule);
		updateFunctionParameters();
		m_base->updateGUI();
		});
	setGroup(tmp, "Scripts");
	setDescription(tmp, "Reload the Python script.");

	SCRIPT_FILE = createStringParameter("scriptFile", "Script file", [&]() { return m_scriptFile; },
		[&](std::string str)
		{
			m_scriptFile = str;
			if (FileSystem::isRelativePath(m_scriptFile))
			{
				const std::string& sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
				m_scriptFile = FileSystem::normalizePath(FileSystem::getFilePath(sceneFile) + "/" + str);
			}

			m_scriptModule = loadScriptFile(m_scriptFile);
			if (m_scriptModule != "")
			{
				m_scriptLoaded = true;
				Embedded::getCurrent()->exec_init(m_scriptModule, m_base);
				updateFunctionParameters();

#ifdef DL_OUTPUT
				// copy script files in output so that the simulation can be reproduced
				std::string sceneFilePath = FileSystem::normalizePath(m_base->getOutputPath() + "/scene");
				FileSystem::makeDirs(sceneFilePath);
				FileSystem::copyFile(m_scriptFile, sceneFilePath + "/" + FileSystem::getFileNameWithExt(m_scriptFile));
#endif
			}
			else
			{
				m_scriptLoaded = false;
				removeFunctionParameters();
			}
			m_base->updateGUI();
		});
	setGroup(SCRIPT_FILE, "Scripts");
	setDescription(SCRIPT_FILE, "Embedded Python script that is executed after each time step of the simulator.");

}

void ScriptObject::updateFunctionParameters()
{
	removeFunctionParameters();
	addFunctionParameters();
}

void ScriptObject::addFunctionParameters()
{
	auto functions = Embedded::getCurrent()->getFunctions();
	int counter = 1;
	for (auto iter = functions.begin(); iter != functions.end(); iter++)
	{
		std::string command_name = *iter;
		std::string name = "embedded_function_" + std::to_string(counter);
		auto tmp = createFunctionParameter(name, command_name, [&, command_name]() {
			if (m_scriptLoaded)
				Embedded::getCurrent()->exec_fct(m_scriptModule, command_name);
			});
		//setHotKey(tmp, "x");
		setGroup(tmp, "Scripts");
		setDescription(tmp, "Execute a function in the script.");
		counter++;
	}
}

void ScriptObject::removeFunctionParameters()
{
	while (m_parameters.size() > 2)
		m_parameters.pop_back();
}

void ScriptObject::init()
{
	m_parameters.clear();
	initParameters();
}

std::string ScriptObject::loadScriptFile(const std::string& fileName)
{
	std::ifstream fp;
	fp.open(fileName.c_str(), std::ios_base::in);
	if (fp)
	{
		std::ostringstream output;
		output << fp.rdbuf();
		fp.close();

		return Embedded::getCurrent()->import_script(fileName);
	}
	else
	{
		LOG_ERR << "Error occurred while loading script: " << fileName;
		return "";
	}
}

void ScriptObject::execResetFct()
{
	if (m_scriptLoaded)
		Embedded::getCurrent()->exec_fct(m_scriptModule, "reset");
}

void ScriptObject::execStepFct()
{
	if (m_scriptLoaded)
		Embedded::getCurrent()->exec_fct(m_scriptModule, "step");
}

