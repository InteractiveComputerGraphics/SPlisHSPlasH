//
// Created by sjeske on 2/3/20.
//
#include "common.h"

#include <Simulator/GUI/Simulator_GUI_Base.h>
#ifdef USE_IMGUI
#include <Simulator/GUI/imgui/Simulator_GUI_imgui.h>
#else
#include <Simulator/GUI/TweakBar/Simulator_GUI_TweakBar.h>
#endif
#include <Simulator/GUI/OpenGL/Simulator_OpenGL.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

#include <iostream>

namespace py = pybind11;

void GUIModule(py::module m) {
    // ---------------------------------------
    // GUI Submodule
    // ---------------------------------------
    auto m_sub = m.def_submodule("GUI");

    // ---------------------------------------
    // Simulator GUI Base Class
    // ---------------------------------------
    py::class_<SPH::Simulator_GUI_Base>(m_sub, "Simulator_GUI_Base")
            .def(py::init<SPH::SimulatorBase*>())
            .def("init", [](SPH::Simulator_GUI_Base& obj, std::vector<std::string> argv, std::string windowName){
                std::vector<const char *> cargv;
                cargv.reserve(argv.size());
                for (auto & elem : argv) cargv.push_back(elem.c_str());
                obj.init(static_cast<int>(argv.size()), const_cast<char**>(cargv.data()), windowName.c_str());
            })
            .def("initSimulationParameterGUI", &SPH::Simulator_GUI_Base::initSimulationParameterGUI)
            .def("initParameterGUI", &SPH::Simulator_GUI_Base::initParameterGUI)
            .def("render", &SPH::Simulator_GUI_Base::render)
            .def("reset", &SPH::Simulator_GUI_Base::reset)
            .def("update", &SPH::Simulator_GUI_Base::update)
            .def("cleanup", &SPH::Simulator_GUI_Base::cleanup)
            .def("run", &SPH::Simulator_GUI_Base::run)
            .def("stop", &SPH::Simulator_GUI_Base::stop)
            .def("addKeyFunc", &SPH::Simulator_GUI_Base::addKeyFunc)
            .def("getSimulatorBase", &SPH::Simulator_GUI_Base::getSimulatorBase, py::return_value_policy::reference_internal);

#ifdef USE_IMGUI
    // ---------------------------------------
    // imgui GUI Class
    // ---------------------------------------
            py::class_<SPH::Simulator_GUI_imgui, SPH::Simulator_GUI_Base>(m_sub, "Simulator_GUI_imgui")
                .def(py::init<SPH::SimulatorBase*>());
#else
    // ---------------------------------------
    // Tweakbar GUI Class
    // ---------------------------------------
    py::class_<SPH::Simulator_GUI_TweakBar, SPH::Simulator_GUI_Base>(m_sub, "Simulator_GUI_TweakBar")
            .def(py::init<SPH::SimulatorBase*>());
            // .def("getTweakBar", &SPH::Simulator_GUI_TweakBar::getTweakBar); //TODO: Wrap tweak bar for this
#endif 

    // ---------------------------------------
    // OpenGL GUI Class // TODO implement missing functions
    // ---------------------------------------
    py::class_<SPH::Simulator_OpenGL>(m_sub, "Simulator_OpenGL")
            .def(py::init<>())
            .def("initShaders", &SPH::Simulator_OpenGL::initShaders)
            .def("renderFluid", &SPH::Simulator_OpenGL::renderFluid)
            .def("renderSelectedParticles", &SPH::Simulator_OpenGL::renderSelectedParticles)
            .def("renderBoundary", &SPH::Simulator_OpenGL::renderBoundary)
            .def("renderBoundaryParticles", &SPH::Simulator_OpenGL::renderBoundaryParticles);
}
