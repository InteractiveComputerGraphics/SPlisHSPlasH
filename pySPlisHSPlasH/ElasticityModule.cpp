//
// Created by sjeske on 1/22/20.
//
#include "common.h"

#include <SPlisHSPlasH/Elasticity/Elasticity_Becker2009.h>
#include <SPlisHSPlasH/Elasticity/Elasticity_Peer2018.h>
#include <SPlisHSPlasH/Elasticity/Elasticity_Kugelstadt2021.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

void ElasticityModule(py::module m_sub) {

    // ---------------------------------------
    // Class Elasticity Becker 2009
    // ---------------------------------------
    py::class_<SPH::Elasticity_Becker2009, SPH::NonPressureForceBase>(m_sub, "Elasticity_Becker2009")
            .def_readwrite_static("YOUNGS_MODULUS", &SPH::Elasticity_Becker2009::YOUNGS_MODULUS)
            .def_readwrite_static("POISSON_RATIO", &SPH::Elasticity_Becker2009::POISSON_RATIO)
            .def_readwrite_static("FIXED_BOX_MIN", &SPH::Elasticity_Becker2009::FIXED_BOX_MIN)
            .def_readwrite_static("FIXED_BOX_MAX", &SPH::Elasticity_Becker2009::FIXED_BOX_MAX)
            .def_readwrite_static("ALPHA", &SPH::Elasticity_Becker2009::ALPHA)
            .def(py::init<SPH::FluidModel*>())
            .def("getMethodName", &SPH::Elasticity_Becker2009::getMethodName);

    // ---------------------------------------
    // Class Elasticity Peer 2018
    // ---------------------------------------
    py::class_<SPH::Elasticity_Peer2018, SPH::NonPressureForceBase>(m_sub, "Elasticity_Peer2018")
            .def_readwrite_static("YOUNGS_MODULUS", &SPH::Elasticity_Peer2018::YOUNGS_MODULUS)
            .def_readwrite_static("POISSON_RATIO", &SPH::Elasticity_Peer2018::POISSON_RATIO)
            .def_readwrite_static("FIXED_BOX_MIN", &SPH::Elasticity_Peer2018::FIXED_BOX_MIN)
            .def_readwrite_static("FIXED_BOX_MAX", &SPH::Elasticity_Peer2018::FIXED_BOX_MAX)
            .def_readwrite_static("ITERATIONS", &SPH::Elasticity_Peer2018::ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &SPH::Elasticity_Peer2018::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &SPH::Elasticity_Peer2018::MAX_ERROR)
            .def_readwrite_static("ALPHA", &SPH::Elasticity_Peer2018::ALPHA)

            .def_static("matrixVecProd", &SPH::Elasticity_Peer2018::matrixVecProd)
            .def(py::init<SPH::FluidModel*>())
            .def("getMethodName", &SPH::Elasticity_Peer2018::getMethodName);


    // ---------------------------------------
    // Class Elasticity Kugelstadt 2021
    // ---------------------------------------

    py::class_<SPH::Elasticity_Kugelstadt2021, SPH::NonPressureForceBase>(m_sub, "Elasticity_Kugelstadt2021")
            .def_readwrite_static("YOUNGS_MODULUS", &SPH::Elasticity_Kugelstadt2021::YOUNGS_MODULUS)
            .def_readwrite_static("POISSON_RATIO", &SPH::Elasticity_Kugelstadt2021::POISSON_RATIO)
            .def_readwrite_static("FIXED_BOX_MIN", &SPH::Elasticity_Kugelstadt2021::FIXED_BOX_MIN)
            .def_readwrite_static("FIXED_BOX_MAX", &SPH::Elasticity_Kugelstadt2021::FIXED_BOX_MAX)
            .def_readwrite_static("ITERATIONS_V", &SPH::Elasticity_Kugelstadt2021::ITERATIONS_V)
            .def_readwrite_static("MAX_ITERATIONS_V", &SPH::Elasticity_Kugelstadt2021::MAX_ITERATIONS_V)
            .def_readwrite_static("MAX_ERROR_V", &SPH::Elasticity_Kugelstadt2021::MAX_ERROR_V)
            .def_readwrite_static("ALPHA", &SPH::Elasticity_Kugelstadt2021::ALPHA)
            .def_readwrite_static("MAX_NEIGHBORS", &SPH::Elasticity_Kugelstadt2021::MAX_NEIGHBORS)

            .def_static("matrixVecProd", &SPH::Elasticity_Kugelstadt2021::matrixVecProd)
            .def("computeRotations", &SPH::Elasticity_Kugelstadt2021::computeRotations)
            .def("getMethodName", &SPH::Elasticity_Kugelstadt2021::getMethodName)
            .def(py::init<SPH::FluidModel*>());
}
