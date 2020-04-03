//
// Created by sjeske on 1/22/20.
//
#include "common.h"

#include <SPlisHSPlasH/Elasticity/ElasticityBase.h>
#include <SPlisHSPlasH/Elasticity/Elasticity_Becker2009.h>
#include <SPlisHSPlasH/Elasticity/Elasticity_Peer2018.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

void ElasticityModule(py::module m_sub) {
    // ---------------------------------------
    // Class Elasticity Base
    // ---------------------------------------
    py::class_<SPH::ElasticityBase, SPH::NonPressureForceBase>(m_sub, "ElasticityBase")
            .def_readwrite_static("YOUNGS_MODULUS", &SPH::ElasticityBase::YOUNGS_MODULUS)
            .def_readwrite_static("POISSON_RATIO", &SPH::ElasticityBase::POISSON_RATIO);

    // ---------------------------------------
    // Class Elasticity Becker 2009
    // ---------------------------------------
    py::class_<SPH::Elasticity_Becker2009, SPH::ElasticityBase>(m_sub, "Elasticity_Becker2009")
            .def_readwrite_static("ALPHA", &SPH::Elasticity_Becker2009::ALPHA)
            .def(py::init<SPH::FluidModel*>());

    // ---------------------------------------
    // Class Elasticity Peer 2018
    // ---------------------------------------
    py::class_<SPH::Elasticity_Peer2018, SPH::ElasticityBase>(m_sub, "Elasticity_Peer2018")
            .def_readwrite_static("ITERATIONS", &SPH::Elasticity_Peer2018::ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &SPH::Elasticity_Peer2018::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &SPH::Elasticity_Peer2018::MAX_ERROR)
            .def_readwrite_static("ALPHA", &SPH::Elasticity_Peer2018::ALPHA)

            .def_static("matrixVecProd", &SPH::Elasticity_Peer2018::matrixVecProd)
            .def(py::init<SPH::FluidModel*>());
}
