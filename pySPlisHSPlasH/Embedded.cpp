#include "Embedded.h"
#include <pybind11/embed.h>
#include <iostream>

#include "common.h"

#include <SPlisHSPlasH/Common.h>
#include <Utilities/Logger.h>
#include <Utilities/Timing.h>
#include <Utilities/Counting.h>
#include <Utilities/FileSystem.h>

namespace py = pybind11;

using namespace SPH;
using namespace Utilities;

Embedded* Embedded::current = nullptr;

Embedded::Embedded () 
{
}

Embedded::~Embedded () 
{
    current->m_mainScope.release();
    py::finalize_interpreter();
    current = 0;
}

Embedded* Embedded::getCurrent ()
{
	if (current == nullptr)
	{
		current = new Embedded ();
        py::initialize_interpreter();
        current->m_mainScope = py::module::import("__main__").attr("__dict__");
	}
	return current;
}

void Embedded::setCurrent (Embedded* tm)
{
	current = tm;
}

bool Embedded::hasCurrent()
{
	return (current != 0);
}


void AnimationFieldModule(py::module);
void ParameterObjectModule(py::module m);
void BoundaryModelModule(py::module);
void EmitterModule(py::module);
void FluidModelModule(py::module);
void GUIModule(py::module);
void ExporterModule(py::module);
void NonPressureForceBaseModule(py::module);
void UtilitiesModule(py::module);
void SimulationModule(py::module);
void SPHKernelsModule(py::module);
void RigidBodyModule(py::module m_sub);
void TimeModule(py::module);
void TriangleMeshModule(py::module);
void PBDModule(py::module);
void DFSPHModule(py::module);
void DragModule(py::module);
void ElasticityModule(py::module);
void IISPHModule(py::module);
void PBFModule(py::module);
void PCISPHModule(py::module);
void PFModule(py::module);
void SurfaceTensionModule(py::module);
void ViscosityModule(py::module);
void VorticityModule(py::module);
void WCSPHModule(py::module);

void ExtrasModule(py::module);


PYBIND11_EMBEDDED_MODULE(splishsplash, m) {
    ParameterObjectModule(m);
    SPHKernelsModule(m);
    AnimationFieldModule(m);
    UtilitiesModule(m);
    BoundaryModelModule(m);
    EmitterModule(m);
    FluidModelModule(m);
    NonPressureForceBaseModule(m);
    SimulationModule(m);
    RigidBodyModule(m);
    TimeModule(m);
    TriangleMeshModule(m);
    PBDModule(m);
    DFSPHModule(m);
    DragModule(m);
    ElasticityModule(m);
    IISPHModule(m);
    PBFModule(m);
    PCISPHModule(m);
    PFModule(m);
    SurfaceTensionModule(m);
    ViscosityModule(m);
    VorticityModule(m);
    WCSPHModule(m);
    GUIModule(m);
    ExporterModule(m);

    ExtrasModule(m);
}

void Embedded::releaseModule(const std::string &moduleName)
{
    try
    {
        m_mainScope[moduleName.c_str()].cast<py::module_>().release();
    }
    catch (py::error_already_set const& pythonErr) { std::cout << pythonErr.what(); }
}

void Embedded::reloadModule(const std::string &moduleName)
{  
    try
    {
        m_mainScope[moduleName.c_str()].cast<py::module_>().reload();
        findFunctions(moduleName);
    }
    catch (py::error_already_set const& pythonErr) { std::cout << pythonErr.what(); }
}

void Embedded::exec_script_file(const std::string& fileName)
{
    std::ifstream fp;
    fp.open(fileName.c_str(), std::ios_base::in);
    if (fp)
    {
        std::ostringstream output;
        output << fp.rdbuf();
        fp.close();

        exec_script_str(output.str());
    }
    else
    {
        LOG_ERR << "Error occurred while loading script: " << fileName;
    }
}

std::string Embedded::import_script(const std::string& fileName)
{
    try
    {
        std::string moduleName = FileSystem::getFileName(fileName);

        // check if module was already imported 
        if (m_mainScope.contains(moduleName.c_str()))
        {
            // reload module
            m_mainScope[moduleName.c_str()].cast<py::module_>().reload();
            findFunctions(moduleName);
        }
        else
        {
            py::module_ sys = py::module_::import("sys");
            sys.attr("path").attr("append")(FileSystem::getFilePath(fileName).c_str());
            auto mod = py::module_::import(moduleName.c_str());
            m_mainScope[moduleName.c_str()] = mod;
            findFunctions(moduleName);
        }
        return moduleName;
    }
    catch (py::error_already_set const& pythonErr) { std::cout << pythonErr.what(); return ""; }
}

void Embedded::findFunctions(const std::string &moduleName)
{
    m_functions.clear();

    py::module_ mod = (*m_mainScope)[moduleName.c_str()].cast<py::module_>();
    if (mod.attr("__dict__").contains("function_list"))
    {
        py::list function_list = mod.attr("function_list");
        for (auto iter = function_list.begin(); iter != function_list.end(); iter++)
        {
            auto value = *iter;
            m_functions.push_back(value.cast<std::string>());
        }
    }

 /*   if (mod.attr("__dict__").contains("function_map"))
    {
        py::dict function_map = mod.attr("function_map");
        for (auto iter = function_map.begin(); iter != function_map.end(); iter++)
        {
            auto value = *iter;
            m_functions.push_back({ value.first.cast<std::string>(),  value.second.cast<std::string>() });
        }
    }*/
}

void Embedded::exec_script_str(const std::string &str)
{
	try 
    {
        py::exec(py::str(str), py::globals(), m_mainScope);
	}
	catch (py::error_already_set const& pythonErr) { std::cout << pythonErr.what(); }
}

void Embedded::exec_fct(const std::string &moduleName, const std::string& fctName)
{
    try
    {
        m_mainScope[moduleName.c_str()].attr(fctName.c_str())();
    }
    catch (py::error_already_set const& pythonErr) { std::cout << pythonErr.what(); }
}

#include "Simulator/SimulatorBase.h"

void Embedded::exec_init(const std::string& moduleName, SimulatorBase *base)
{
    try
    {
        m_mainScope[moduleName.c_str()].attr("init")(base);
    }
    catch (py::error_already_set const& pythonErr) { std::cout << pythonErr.what(); }
}