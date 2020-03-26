//
// Created by sjeske on 1/22/20.
//
#include "common.h"

#include <SPlisHSPlasH/TimeStep.h>
#include <SPlisHSPlasH/DFSPH/SimulationDataDFSPH.h>
#include <SPlisHSPlasH/DFSPH/TimeStepDFSPH.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

void DFSPHModule(py::module m_sub) {
    // ---------------------------------------
    // Class Simulation Data DFSPH
    // ---------------------------------------
    py::class_<SPH::SimulationDataDFSPH>(m_sub, "SimulationDataDFSPH")
            .def(py::init<>())
            .def("init", &SPH::SimulationDataDFSPH::init)
            .def("cleanup", &SPH::SimulationDataDFSPH::cleanup)
            .def("reset", &SPH::SimulationDataDFSPH::reset)
            .def("performNeighborhoodSearchSort", &SPH::SimulationDataDFSPH::performNeighborhoodSearchSort)

            .def("getFactor", (const Real (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned int)const)(&SPH::SimulationDataDFSPH::getFactor))
            // .def("getFactor", (Real& (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned int))(&SPH::SimulationDataDFSPH::getFactor)) // TODO: wont work by reference
            .def("setFactor", &SPH::SimulationDataDFSPH::setFactor)

            .def("getKappa", (const Real (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned int)const)(&SPH::SimulationDataDFSPH::getKappa))
            // .def("getKappa", (Real& (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned int))(&SPH::SimulationDataDFSPH::getKappa)) // TODO: wont work by reference
            .def("setKappa", &SPH::SimulationDataDFSPH::setKappa)

            .def("getKappaV", (const Real (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned int)const)(&SPH::SimulationDataDFSPH::getKappaV))
            // .def("getKappaV", (Real& (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned int))(&SPH::SimulationDataDFSPH::getKappaV)) // TODO: wont work by reference
            .def("setKappaV", &SPH::SimulationDataDFSPH::setKappaV)

            .def("getDensityAdv", (const Real (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned int)const)(&SPH::SimulationDataDFSPH::getDensityAdv))
            // .def("getDensityAdv", (Real& (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned int))(&SPH::SimulationDataDFSPH::getDensityAdv)) // TODO: wont work by reference
            .def("setDensityAdv", &SPH::SimulationDataDFSPH::setDensityAdv);

    // ---------------------------------------
    // Class Time Step DFSPH
    // ---------------------------------------
    py::class_<SPH::TimeStepDFSPH, SPH::TimeStep>(m_sub, "TimeStepDFSPH")
            .def_readwrite_static("SOLVER_ITERATIONS_V", &SPH::TimeStepDFSPH::SOLVER_ITERATIONS_V)
            .def_readwrite_static("MAX_ITERATIONS_V", &SPH::TimeStepDFSPH::MAX_ITERATIONS_V)
            .def_readwrite_static("MAX_ERROR_V", &SPH::TimeStepDFSPH::MAX_ERROR_V)
            .def_readwrite_static("USE_DIVERGENCE_SOLVER", &SPH::TimeStepDFSPH::USE_DIVERGENCE_SOLVER)

            .def(py::init<>());
}