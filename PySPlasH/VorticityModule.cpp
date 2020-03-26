//
// Created by sjeske on 1/28/20.
//
#include "common.h"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SPlisHSPlasH/Common.h>
#include <SPlisHSPlasH/NonPressureForceBase.h>
#include <SPlisHSPlasH/Vorticity/VorticityBase.h>
#include <SPlisHSPlasH/Viscosity/Viscosity_Bender2017.h>
#include <SPlisHSPlasH/Vorticity/MicropolarModel_Bender2017.h>
#include <SPlisHSPlasH/Vorticity/VorticityConfinement.h>

namespace py = pybind11;

void VorticityModule(py::module m_sub) {
    // ---------------------------------------
    // Vorticity Base
    // ---------------------------------------
    py::class_<SPH::VorticityBase, SPH::NonPressureForceBase>(m_sub, "VorticityBase")
            .def_readwrite_static("VORTICITY_COEFFICIENT", &SPH::VorticityBase::VORTICITY_COEFFICIENT);

    // ---------------------------------------
    // Vorticity Bender 2017
    // ---------------------------------------
    py::class_<SPH::MicropolarModel_Bender2017, SPH::VorticityBase>(m_sub, "MicropolarModel_Bender2017")
            .def_readwrite_static("VISCOSITY_OMEGA", &SPH::MicropolarModel_Bender2017::VISCOSITY_OMEGA)
            .def_readwrite_static("INERTIA_INVERSE", &SPH::MicropolarModel_Bender2017::INERTIA_INVERSE)

            .def(py::init<SPH::FluidModel*>())
            .def("getAngularAcceleration", (const Vector3r& (SPH::MicropolarModel_Bender2017::*)(const unsigned int)const)&SPH::MicropolarModel_Bender2017::getAngularAcceleration)
            // .def("getAngularAcceleration", (Vector3r& (SPH::MicropolarModel_Bender2017::*)(const unsigned int))&SPH::MicropolarModel_Bender2017::getAngularAcceleration) // TODO: wont work by reference
            .def("setAngularAcceleration", &SPH::MicropolarModel_Bender2017::setAngularAcceleration)
            .def("getAngularVelocity", (const Vector3r& (SPH::MicropolarModel_Bender2017::*)(const unsigned int)const)&SPH::MicropolarModel_Bender2017::getAngularVelocity)
            // .def("getAngularVelocity", (Vector3r& (SPH::MicropolarModel_Bender2017::*)(const unsigned int))&SPH::MicropolarModel_Bender2017::getAngularVelocity) // TODO: wont work by reference
            .def("setAngularVelocity", &SPH::MicropolarModel_Bender2017::setAngularVelocity);

    // ---------------------------------------
    // Vorticity Bender 2017
    // ---------------------------------------
    py::class_<SPH::VorticityConfinement, SPH::VorticityBase>(m_sub, "VorticityConfinement")
            .def(py::init<SPH::FluidModel*>())
            .def("getAngularVelocity", (const Vector3r& (SPH::VorticityConfinement::*)(const unsigned int)const)&SPH::VorticityConfinement::getAngularVelocity)
            // .def("getAngularVelocity", (Vector3r& (SPH::VorticityConfinement::*)(const unsigned int))&SPH::VorticityConfinement::getAngularVelocity) // TODO: wont work by reference
            .def("setAngularVelocity", &SPH::VorticityConfinement::setAngularVelocity);
}