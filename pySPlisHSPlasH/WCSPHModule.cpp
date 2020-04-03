//
// Created by sjeske on 1/28/20.
//
#include "common.h"

#include <pybind11/pybind11.h>
#include <SPlisHSPlasH/WCSPH/SimulationDataWCSPH.h>
#include <SPlisHSPlasH/WCSPH/TimeStepWCSPH.h>

namespace py = pybind11;

void WCSPHModule(py::module m_sub) {
    // ---------------------------------------
    // Class Simulation Data WCSPH
    // ---------------------------------------
    py::class_<SPH::SimulationDataWCSPH>(m_sub, "SimulationDataWCSPH")
            .def(py::init<>())
            .def("init", &SPH::SimulationDataWCSPH::init)
            .def("cleanup", &SPH::SimulationDataWCSPH::cleanup)
            .def("reset", &SPH::SimulationDataWCSPH::reset)
            .def("performNeighborhoodSearchSort", &SPH::SimulationDataWCSPH::performNeighborhoodSearchSort)
            .def("emittedParticles", &SPH::SimulationDataWCSPH::emittedParticles)
            .def("getPressure", (const Real (SPH::SimulationDataWCSPH::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataWCSPH::getPressure)
            // .def("getPressure", (Real& (SPH::SimulationDataWCSPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataWCSPH::getPressure) // TODO: wont work by reference
            .def("setPressure", &SPH::SimulationDataWCSPH::setPressure)
            .def("getPressureAccel", (const Vector3r& (SPH::SimulationDataWCSPH::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataWCSPH::getPressureAccel)
            // .def("getPressureAccel", (Vector3r& (SPH::SimulationDataWCSPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataWCSPH::getPressureAccel) // TODO: wont work by reference
            .def("setPressureAccel", &SPH::SimulationDataWCSPH::setPressureAccel);

    // ---------------------------------------
    // Time Step WCSPH
    // ---------------------------------------
    py::class_<SPH::TimeStepWCSPH, SPH::TimeStep>(m_sub, "TimeStepWCSPH")
            .def_readwrite_static("STIFFNESS", &SPH::TimeStepWCSPH::STIFFNESS)
            .def_readwrite_static("EXPONENT", &SPH::TimeStepWCSPH::EXPONENT)
            .def(py::init<>());
}