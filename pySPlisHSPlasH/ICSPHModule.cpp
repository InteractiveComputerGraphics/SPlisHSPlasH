#include "common.h"

#include <SPlisHSPlasH/TimeStep.h>
#include <SPlisHSPlasH/ICSPH/SimulationDataICSPH.h>
#include <SPlisHSPlasH/ICSPH/TimeStepICSPH.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

void ICSPHModule(py::module m_sub) {
    // ---------------------------------------
    // Class Simulation Data IISPH
    // ---------------------------------------
    py::class_<SPH::SimulationDataICSPH>(m_sub, "SimulationDataICSPH")
            .def(py::init<>())
            .def("init", &SPH::SimulationDataICSPH::init)
            .def("cleanup", &SPH::SimulationDataICSPH::cleanup)
            .def("reset", &SPH::SimulationDataICSPH::reset)
            .def("performNeighborhoodSearchSort", &SPH::SimulationDataICSPH::performNeighborhoodSearchSort)
            .def("emittedParticles", &SPH::SimulationDataICSPH::emittedParticles)

            .def("getAii", (const Real (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned int) const)(&SPH::SimulationDataICSPH::getAii))
            // .def("getAii", (Real& (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned int))(&SPH::SimulationDataICSPH::getAii)) // TODO: wont work by reference
            .def("setAii", &SPH::SimulationDataICSPH::setAii)

            .def("getDensityAdv", (const Real (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned int) const)&SPH::SimulationDataICSPH::getDensityAdv)
            // .def("getDensityAdv", (Real& (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataICSPH::getDensityAdv) // TODO: wont work by reference
            .def("setDensityAdv", &SPH::SimulationDataICSPH::setDensityAdv)

            .def("getPressure", (const Real (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned int) const)&SPH::SimulationDataICSPH::getPressure)
            // .def("getPressure", (Real& (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataICSPH::getPressure) // TODO: wont work by reference
            .def("setPressure", &SPH::SimulationDataICSPH::setPressure)

            .def("getPressureAccel", (const Vector3r &(SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned int) const)&SPH::SimulationDataICSPH::getPressureAccel)
            // .def("getPressureAccel", (Vector3r& (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataICSPH::getPressureAccel) // TODO: wont work by reference
            .def("setPressureAccel", &SPH::SimulationDataICSPH::setPressureAccel)

            .def("getPressureGradient", (const Vector3r & (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned int) const) & SPH::SimulationDataICSPH::getPressureGradient)
            // .def("getPressureGradient", (Vector3r& (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned int))&SPH::SimulationDataICSPH::getPressureAccel) // TODO: wont work by reference
            .def("setPressureGradient", &SPH::SimulationDataICSPH::setPressureGradient);

    // ---------------------------------------
    // Class Simulation Data ICSPH
    // ---------------------------------------
    py::class_<SPH::TimeStepICSPH, SPH::TimeStep>(m_sub, "TimeStepICSPH")
            .def_readwrite_static("LAMBDA", &SPH::TimeStepICSPH::LAMBDA)
            .def_readwrite_static("PRESSURE_CLAMPING", &SPH::TimeStepICSPH::PRESSURE_CLAMPING)
            .def("getSimulationData", &SPH::TimeStepICSPH::getSimulationData)
            .def(py::init<>());
}

