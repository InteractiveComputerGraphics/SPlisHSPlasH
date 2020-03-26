//
// Created by sjeske on 1/22/20.
//
#include "common.h"

#include <SPlisHSPlasH/TimeStep.h>
#include <SPlisHSPlasH/IISPH/SimulationDataIISPH.h>
#include <SPlisHSPlasH/IISPH/TimeStepIISPH.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

void IISPHModule(py::module m_sub) {
    // ---------------------------------------
    // Class Simulation Data IISPH
    // ---------------------------------------
    py::class_<SPH::SimulationDataIISPH>(m_sub, "SimulationDataIISPH")
            .def(py::init<>())
            .def("init", &SPH::SimulationDataIISPH::init)
            .def("cleanup", &SPH::SimulationDataIISPH::cleanup)
            .def("reset", &SPH::SimulationDataIISPH::reset)
            .def("performNeighborhoodSearchSort", &SPH::SimulationDataIISPH::performNeighborhoodSearchSort)
            .def("emittedParticles", &SPH::SimulationDataIISPH::emittedParticles)

            .def("getAii", (const Real (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const)(&SPH::SimulationDataIISPH::getAii))
            // .def("getAii", (Real& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int))(&SPH::SimulationDataIISPH::getAii)) // TODO: wont work by reference
            .def("setAii", &SPH::SimulationDataIISPH::setAii)

            .def("getDii", (const Vector3r &(SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const)&SPH::SimulationDataIISPH::getDii)
            // .def("getDii", (Vector3r& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataIISPH::getDii) // TODO: wont work by reference
            .def("setDii", &SPH::SimulationDataIISPH::setDii)

            .def("getDij_pj", (const Vector3r &(SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const)&SPH::SimulationDataIISPH::getDij_pj)
            // .def("getDij_pj", (Vector3r& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataIISPH::getDij_pj) // TODO: wont work by reference
            .def("setDij_pj", &SPH::SimulationDataIISPH::setDij_pj)

            .def("getDensityAdv", (const Real (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const)&SPH::SimulationDataIISPH::getDensityAdv)
            // .def("getDensityAdv", (Real& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataIISPH::getDensityAdv) // TODO: wont work by reference
            .def("setDensityAdv", &SPH::SimulationDataIISPH::setDensityAdv)

            .def("getPressure", (const Real (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const)&SPH::SimulationDataIISPH::getPressure)
            // .def("getPressure", (Real& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataIISPH::getPressure) // TODO: wont work by reference
            .def("setPressure", &SPH::SimulationDataIISPH::setPressure)

            .def("getLastPressure", (const Real (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const)&SPH::SimulationDataIISPH::getLastPressure)
            // .def("getLastPressure", (Real& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataIISPH::getLastPressure) // TODO: wont work by reference
            .def("setLastPressure", &SPH::SimulationDataIISPH::setLastPressure)

            .def("getPressureAccel", (const Vector3r &(SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const)&SPH::SimulationDataIISPH::getPressureAccel)
            // .def("getPressureAccel", (Vector3r& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataIISPH::getPressureAccel) // TODO: wont work by reference
            .def("setPressureAccel", &SPH::SimulationDataIISPH::setPressureAccel);

    // ---------------------------------------
    // Class Simulation Data IISPH
    // ---------------------------------------
    py::class_<SPH::TimeStepIISPH, SPH::TimeStep>(m_sub, "TimeStepIISPH")
            .def("getSimulationData", &SPH::TimeStepIISPH::getSimulationData)
            .def(py::init<>());
}

