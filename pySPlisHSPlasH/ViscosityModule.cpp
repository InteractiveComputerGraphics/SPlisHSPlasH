//
// Created by sjeske on 1/28/20.
//
#include "common.h"

#include <pybind11/pybind11.h>

#include <SPlisHSPlasH/Viscosity/Viscosity_Bender2017.h>
#include <SPlisHSPlasH/Viscosity/Viscosity_Peer2015.h>
#include <SPlisHSPlasH/Viscosity/Viscosity_Peer2016.h>
#include <SPlisHSPlasH/Viscosity/Viscosity_Standard.h>
#include <SPlisHSPlasH/Viscosity/Viscosity_Takahashi2015.h>
#include <SPlisHSPlasH/Viscosity/Viscosity_Weiler2018.h>

namespace py = pybind11;

void ViscosityModule(py::module m_sub) {

    // ---------------------------------------
    // Viscosity Bender 2017
    // ---------------------------------------
    py::class_<SPH::Viscosity_Bender2017, SPH::NonPressureForceBase>(m_sub, "Viscosity_Bender2017")
            .def_readwrite_static("VISCOSITY_COEFFICIENT", &SPH::Viscosity_Bender2017::VISCOSITY_COEFFICIENT)
            .def_readwrite_static("ITERATIONS", &SPH::Viscosity_Bender2017::ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &SPH::Viscosity_Bender2017::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &SPH::Viscosity_Bender2017::MAX_ERROR)

            .def(py::init<SPH::FluidModel*>())
            .def("computeTargetStrainRate", &SPH::Viscosity_Bender2017::computeTargetStrainRate)
            .def("computeViscosityFactor", &SPH::Viscosity_Bender2017::computeViscosityFactor)
            .def("viscoGradientMultTransposeRightOpt", &SPH::Viscosity_Bender2017::viscoGradientMultTransposeRightOpt)
            .def("getTargetStrainRate", (const Vector6r& (SPH::Viscosity_Bender2017::*)(const unsigned int)const)(&SPH::Viscosity_Bender2017::getTargetStrainRate))
            // .def("getTargetStrainRate", (Vector6r& (SPH::Viscosity_Bender2017::*)(const unsigned int))&SPH::Viscosity_Bender2017::getTargetStrainRate) // TODO: wont work by reference
            .def("setTargetStrainRate", &SPH::Viscosity_Bender2017::setTargetStrainRate)
            .def("getViscosityFactor", (const Matrix6r& (SPH::Viscosity_Bender2017::*)(const unsigned int)const)&SPH::Viscosity_Bender2017::getViscosityFactor)
            // .def("getViscosityFactor", (Matrix6r& (SPH::Viscosity_Bender2017::*)(const unsigned int))&SPH::Viscosity_Bender2017::getViscosityFactor) // TODO: wont work by reference
            .def("setViscosityFactor", &SPH::Viscosity_Bender2017::setViscosityFactor)
            .def("getViscosityLambda", (const Vector6r& (SPH::Viscosity_Bender2017::*)(const unsigned int)const)&SPH::Viscosity_Bender2017::getViscosityLambda)
            // .def("getViscosityLambda", (Vector6r& (SPH::Viscosity_Bender2017::*)(const unsigned int))&SPH::Viscosity_Bender2017::getViscosityLambda) // TODO: wont work by reference
            .def("setViscosityLambda", &SPH::Viscosity_Bender2017::setViscosityLambda)
            .def("getMethodName", &SPH::Viscosity_Bender2017::getMethodName);

    // ---------------------------------------
    // Viscosity Peer 2015
    // ---------------------------------------
    py::class_<SPH::Viscosity_Peer2015, SPH::NonPressureForceBase>(m_sub, "Viscosity_Peer2015")
            .def_readwrite_static("VISCOSITY_COEFFICIENT", &SPH::Viscosity_Peer2015::VISCOSITY_COEFFICIENT)
            .def_readwrite_static("ITERATIONS", &SPH::Viscosity_Peer2015::ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &SPH::Viscosity_Peer2015::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &SPH::Viscosity_Peer2015::MAX_ERROR)
            .def(py::init<SPH::FluidModel*>())
            .def_static("matrixVecProd", &SPH::Viscosity_Peer2015::matrixVecProd)
            .def_static("diagonalMatrixElement", &SPH::Viscosity_Peer2015::diagonalMatrixElement)
            .def("getTargetNablaV", (const Matrix3r& (SPH::Viscosity_Peer2015::*)(const unsigned int)const)(&SPH::Viscosity_Peer2015::getTargetNablaV))
            // .def("getTargetNablaV", (Matrix3r& (SPH::Viscosity_Peer2015::*)(const unsigned int))&SPH::Viscosity_Peer2015::getTargetNablaV) // TODO: wont work by reference
            .def("setTargetNablaV", &SPH::Viscosity_Peer2015::setTargetNablaV)
            .def("getMethodName", &SPH::Viscosity_Peer2015::getMethodName);

    // ---------------------------------------
    // Viscosity Peer2016
    // ---------------------------------------
    py::class_<SPH::Viscosity_Peer2016, SPH::NonPressureForceBase>(m_sub, "Viscosity_Peer2016")
            .def_readwrite_static("VISCOSITY_COEFFICIENT", &SPH::Viscosity_Peer2016::VISCOSITY_COEFFICIENT)
            .def_readwrite_static("ITERATIONS_V", &SPH::Viscosity_Peer2016::ITERATIONS_V)
            .def_readwrite_static("ITERATIONS_OMEGA", &SPH::Viscosity_Peer2016::ITERATIONS_OMEGA)
            .def_readwrite_static("MAX_ITERATIONS_V", &SPH::Viscosity_Peer2016::MAX_ITERATIONS_V)
            .def_readwrite_static("MAX_ERROR_V", &SPH::Viscosity_Peer2016::MAX_ERROR_V)
            .def_readwrite_static("MAX_ITERATIONS_OMEGA", &SPH::Viscosity_Peer2016::MAX_ITERATIONS_OMEGA)
            .def_readwrite_static("MAX_ERROR_OMEGA", &SPH::Viscosity_Peer2016::MAX_ERROR_OMEGA)

            .def(py::init<SPH::FluidModel*>())
            .def_static("matrixVecProdV", &SPH::Viscosity_Peer2016::matrixVecProdV)
            .def_static("diagonalMatrixElementV", &SPH::Viscosity_Peer2016::diagonalMatrixElementV)
            .def_static("matrixVecProdOmega", &SPH::Viscosity_Peer2016::matrixVecProdOmega)
            .def_static("diagonalMatrixElementOmega", &SPH::Viscosity_Peer2016::diagonalMatrixElementOmega)
            .def("getTargetNablaV", (const Matrix3r& (SPH::Viscosity_Peer2016::*)(const unsigned int)const)(&SPH::Viscosity_Peer2016::getTargetNablaV))
            // .def("getTargetNablaV", (Matrix3r& (SPH::Viscosity_Peer2016::*)(const unsigned int))&SPH::Viscosity_Peer2016::getTargetNablaV) // TODO: wont work by reference
            .def("setTargetNablaV", &SPH::Viscosity_Peer2016::setTargetNablaV)
            .def("getOmega", (const Vector3r& (SPH::Viscosity_Peer2016::*)(const unsigned int)const)&SPH::Viscosity_Peer2016::getOmega)
            // .def("getOmega", ( Vector3r& (SPH::Viscosity_Peer2016::*)(const unsigned int))&SPH::Viscosity_Peer2016::getOmega) // TODO: wont work by reference
            .def("setOmega", &SPH::Viscosity_Peer2016::setOmega)
            .def("getMethodName", &SPH::Viscosity_Peer2016::getMethodName);

    // ---------------------------------------
    // Viscosity Standard
    // ---------------------------------------
    py::class_<SPH::Viscosity_Standard, SPH::NonPressureForceBase>(m_sub, "Viscosity_Standard")
            .def_readwrite_static("VISCOSITY_COEFFICIENT", &SPH::Viscosity_Standard::VISCOSITY_COEFFICIENT)
            .def_readwrite_static("VISCOSITY_COEFFICIENT_BOUNDARY", &SPH::Viscosity_Standard::VISCOSITY_COEFFICIENT_BOUNDARY)
            .def(py::init<SPH::FluidModel*>())
            .def("getMethodName", &SPH::Viscosity_Standard::getMethodName);

    // ---------------------------------------
    // Viscosity
    // ---------------------------------------
    py::class_<SPH::Viscosity_Takahashi2015, SPH::NonPressureForceBase>(m_sub, "Viscosity_Takahashi2015")
            .def_readwrite_static("VISCOSITY_COEFFICIENT", &SPH::Viscosity_Takahashi2015::VISCOSITY_COEFFICIENT)
            .def_readwrite_static("ITERATIONS", &SPH::Viscosity_Takahashi2015::ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &SPH::Viscosity_Takahashi2015::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &SPH::Viscosity_Takahashi2015::MAX_ERROR)

            .def(py::init<SPH::FluidModel*>())
            .def_static("matrixVecProd", &SPH::Viscosity_Takahashi2015::matrixVecProd)
            .def("getViscousStress", (const Matrix3r& (SPH::Viscosity_Takahashi2015::*)(const unsigned int)const) &SPH::Viscosity_Takahashi2015::getViscousStress)
            // .def("getViscousStress", (Matrix3r& (SPH::Viscosity_Takahashi2015::*)(const unsigned int))&SPH::Viscosity_Takahashi2015::getViscousStress) // TODO: wont work by reference
            .def("setViscousStress", &SPH::Viscosity_Takahashi2015::setViscousStress)
            .def("getAccel", (const Vector3r& (SPH::Viscosity_Takahashi2015::*)(const unsigned int)const)&SPH::Viscosity_Takahashi2015::getAccel)
            // .def("getAccel", (Vector3r& (SPH::Viscosity_Takahashi2015::*)(const unsigned int))&SPH::Viscosity_Takahashi2015::getAccel) // TODO: wont work by reference
            .def("setAccel", &SPH::Viscosity_Takahashi2015::setAccel)
            .def("getMethodName", &SPH::Viscosity_Takahashi2015::getMethodName);

    // ---------------------------------------
    // Viscosity Weiler 2018
    // ---------------------------------------
    py::class_<SPH::Viscosity_Weiler2018, SPH::NonPressureForceBase>(m_sub, "Viscosity_Weiler2018")
            .def_readwrite_static("VISCOSITY_COEFFICIENT", &SPH::Viscosity_Weiler2018::VISCOSITY_COEFFICIENT)
            .def_readwrite_static("ITERATIONS", &SPH::Viscosity_Weiler2018::ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &SPH::Viscosity_Weiler2018::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &SPH::Viscosity_Weiler2018::MAX_ERROR)
            .def_readwrite_static("VISCOSITY_COEFFICIENT_BOUNDARY", &SPH::Viscosity_Weiler2018::VISCOSITY_COEFFICIENT_BOUNDARY)

            .def(py::init<SPH::FluidModel*>())
            .def_static("matrixVecProd", &SPH::Viscosity_Weiler2018::matrixVecProd)
            .def("getVDiff", (const Vector3r& (SPH::Viscosity_Weiler2018::*)(const unsigned int)const)& SPH::Viscosity_Weiler2018::getVDiff)
            .def("setVDiff", &SPH::Viscosity_Weiler2018::setVDiff)
            .def("getMethodName", &SPH::Viscosity_Weiler2018::getMethodName);
}