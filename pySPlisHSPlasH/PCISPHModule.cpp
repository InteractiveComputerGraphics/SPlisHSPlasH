//
// Created by sjeske on 1/28/20.
//
#include "common.h"

#include <SPlisHSPlasH/PCISPH/SimulationDataPCISPH.h>
#include <SPlisHSPlasH/PCISPH/TimeStepPCISPH.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

void PCISPHModule(py::module m_sub) {
    // ---------------------------------------
    // Class Simulation Data PCISPH
    // ---------------------------------------
    py::class_<SPH::SimulationDataPCISPH>(m_sub, "SimulationDataPCISPH")
            .def(py::init<>())
            .def("init", &SPH::SimulationDataPCISPH::init)
            .def("cleanup", &SPH::SimulationDataPCISPH::cleanup)
            .def("reset", &SPH::SimulationDataPCISPH::reset)
            .def("performNeighborhoodSearchSort", &SPH::SimulationDataPCISPH::performNeighborhoodSearchSort)
            .def("getPCISPH_ScalingFactor", &SPH::SimulationDataPCISPH::getPCISPH_ScalingFactor)
            .def("emittedParticles", &SPH::SimulationDataPCISPH::emittedParticles)

            .def("getLastPosition", (const Vector3r& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned int)const)(&SPH::SimulationDataPCISPH::getLastPosition))
            // .def("getLastPosition", (Vector3r& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataPCISPH::getLastPosition) // TODO: wont work by reference
            .def("setLastPosition", &SPH::SimulationDataPCISPH::setLastPosition)

            .def("getLastVelocity", (const Vector3r& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataPCISPH::getLastVelocity)
            // .def("getLastVelocity", (Vector3r& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataPCISPH::getLastVelocity) // TODO: wont work by reference
            .def("setLastVelocity", &SPH::SimulationDataPCISPH::setLastVelocity)

            .def("getDensityAdv", (const Real (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataPCISPH::getDensityAdv)
            // .def("getDensityAdv", (Real& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataPCISPH::getDensityAdv) // TODO: wont work by reference
            .def("setDensityAdv", &SPH::SimulationDataPCISPH::setDensityAdv)

            .def("getPressure", (const Real (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataPCISPH::getPressure)
            // .def("getPressure", (Real& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataPCISPH::getPressure) // TODO: wont work by reference
            .def("setPressure", &SPH::SimulationDataPCISPH::setPressure)

            .def("getPressureAccel", (const Vector3r& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataPCISPH::getPressureAccel)
            // .def("getPressureAccel", (Vector3r& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataPCISPH::getPressureAccel) // TODO: wont work by reference
            .def("setPressureAccel", &SPH::SimulationDataPCISPH::setPressureAccel);

    // ---------------------------------------
    // Class Time Step PCISPH
    // ---------------------------------------
    py::class_<SPH::TimeStepPCISPH, SPH::TimeStep>(m_sub, "TimeStepPCISPH")
            .def(py::init<>());
}