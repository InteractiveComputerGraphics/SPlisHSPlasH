//
// Created by sjeske on 1/28/20.
//
#include "common.h"

#include <SPlisHSPlasH/PBF/SimulationDataPBF.h>
#include <SPlisHSPlasH/PBF/TimeIntegration.h>
#include <SPlisHSPlasH/PBF/TimeStepPBF.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

void PBFModule(py::module m_sub) {
    // ---------------------------------------
    // Class Simulation Data PBF
    // ---------------------------------------
    py::class_<SPH::SimulationDataPBF>(m_sub, "SimulationDataPBF")
            .def(py::init<>())
            .def("init", &SPH::SimulationDataPBF::init)
            .def("cleanup", &SPH::SimulationDataPBF::cleanup)
            .def("reset", &SPH::SimulationDataPBF::reset)
            .def("performNeighborhoodSearchSort", &SPH::SimulationDataPBF::performNeighborhoodSearchSort)
            .def("emittedParticles", &SPH::SimulationDataPBF::emittedParticles)

            .def("getLambda", (const Real & (SPH::SimulationDataPBF::*)(const unsigned int, const unsigned int)const)(&SPH::SimulationDataPBF::getLambda))
            // .def("getLambda", (Real & (SPH::SimulationDataPBF::*)(const unsigned int, const unsigned int))(&SPH::SimulationDataPBF::getLambda)) // TODO: wont work by reference
            .def("setLambda", &SPH::SimulationDataPBF::setLambda)

            .def("getDeltaX", (const Vector3r & (SPH::SimulationDataPBF::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataPBF::getDeltaX)
            // .def("getDeltaX", (Vector3r & (SPH::SimulationDataPBF::*)(const unsigned int, const unsigned int))&SPH::SimulationDataPBF::getDeltaX) // TODO: wont work by reference
            .def("setDeltaX", &SPH::SimulationDataPBF::setDeltaX)

            .def("getLastPosition", (const Vector3r & (SPH::SimulationDataPBF::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataPBF::getLastPosition)
            // .def("getLastPosition", (Vector3r & (SPH::SimulationDataPBF::*)(const unsigned int, const unsigned int))&SPH::SimulationDataPBF::getLastPosition) // TODO: wont work by reference
            .def("setLastPosition", &SPH::SimulationDataPBF::setLastPosition)

            .def("getOldPosition", (const Vector3r & (SPH::SimulationDataPBF::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataPBF::getOldPosition)
            // .def("getOldPosition", (Vector3r & (SPH::SimulationDataPBF::*)(const unsigned int, const unsigned int))&SPH::SimulationDataPBF::getOldPosition) // TODO: wont work by reference
            .def("setOldPosition", &SPH::SimulationDataPBF::setOldPosition);

    // ---------------------------------------
    // Class Time Integration
    // ---------------------------------------
    py::class_<SPH::TimeIntegration>(m_sub, "TimeIntegration")
            .def_static("semiImplicitEuler", &SPH::TimeIntegration::semiImplicitEuler)
            .def_static("velocityUpdateFirstOrder", &SPH::TimeIntegration::velocityUpdateFirstOrder)
            .def_static("velocityUpdateSecondOrder", &SPH::TimeIntegration::velocityUpdateSecondOrder);

    // ---------------------------------------
    // Class Time Step PBF
    // ---------------------------------------
    py::class_<SPH::TimeStepPBF, SPH::TimeStep>(m_sub, "TimeStepPBF")
            .def_readwrite_static("VELOCITY_UPDATE_METHOD", &SPH::TimeStepPBF::VELOCITY_UPDATE_METHOD)
            .def_readwrite_static("ENUM_PBF_FIRST_ORDER", &SPH::TimeStepPBF::ENUM_PBF_FIRST_ORDER)
            .def_readwrite_static("ENUM_PBF_SECOND_ORDER", &SPH::TimeStepPBF::ENUM_PBF_SECOND_ORDER)
            .def(py::init<>());
}