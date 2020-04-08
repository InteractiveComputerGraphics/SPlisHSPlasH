//
// Created by sjeske on 1/28/20.
//
#include "common.h"

#include <SPlisHSPlasH/PF/SimulationDataPF.h>
#include <SPlisHSPlasH/PF/TimeStepPF.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

void PFModule(py::module m_sub) {
    // ---------------------------------------
    // Class Simulation Data PF
    // ---------------------------------------
    py::class_<SPH::SimulationDataPF>(m_sub, "SimulationDataPF")
            .def(py::init<>())
            .def("init", &SPH::SimulationDataPF::init)
            .def("cleanup", &SPH::SimulationDataPF::cleanup)
            .def("reset", &SPH::SimulationDataPF::reset)
            .def("performNeighborhoodSearchSort", &SPH::SimulationDataPF::performNeighborhoodSearchSort)
            .def("emittedParticles", &SPH::SimulationDataPF::emittedParticles)

            .def("getOldPosition", (const Vector3r (SPH::SimulationDataPF::*)(const unsigned int, const unsigned int)const)(&SPH::SimulationDataPF::getOldPosition))
            // .def("getOldPosition", (Vector3r& (SPH::SimulationDataPF::*)(const unsigned int, const unsigned int))(&SPH::SimulationDataPF::getOldPosition)) // TODO: wont work by reference
            .def("setOldPosition", &SPH::SimulationDataPF::setOldPosition)

            .def("getNumFluidNeighbors", (const unsigned int (SPH::SimulationDataPF::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataPF::getNumFluidNeighbors)
            // .def("getNumFluidNeighbors", (unsigned int& (SPH::SimulationDataPF::*)(const unsigned int, const unsigned int))&SPH::SimulationDataPF::getNumFluidNeighbors) // TODO: wont work by reference
            .def("setNumFluidNeighbors", &SPH::SimulationDataPF::setNumFluidNeighbors)

            .def("getS", (const Vector3r& (SPH::SimulationDataPF::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataPF::getS)
            // .def("getS", (Vector3r& (SPH::SimulationDataPF::*)(const unsigned int, const unsigned int))&SPH::SimulationDataPF::getS) // TODO: wont work by reference
            .def("setS", &SPH::SimulationDataPF::setS)

            .def("getDiag", (const Vector3r& (SPH::SimulationDataPF::*)(const unsigned int, const unsigned int)const)&SPH::SimulationDataPF::getDiag)
            // .def("getDiag", (Vector3r& (SPH::SimulationDataPF::*)(const unsigned int, const unsigned int))&SPH::SimulationDataPF::getDiag) // TODO: wont work by reference
            .def("setDiag", &SPH::SimulationDataPF::setDiag)

            .def("getParticleOffset", &SPH::SimulationDataPF::getParticleOffset);

    // ---------------------------------------
    // Class Time Step PF
    // ---------------------------------------
    py::class_<SPH::TimeStepPF, SPH::TimeStep>(m_sub, "TimeStepPF")
            .def_readwrite_static("STIFFNESS", &SPH::TimeStepPF::STIFFNESS)
            .def(py::init<>())
            .def_static("matrixVecProd", &SPH::TimeStepPF::matrixVecProd);
}