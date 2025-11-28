//
// Created by sjeske on 1/28/20.
//
#include "common.h"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SPlisHSPlasH/Common.h>
#include <SPlisHSPlasH/NonPressureForceBase.h>
#include <SPlisHSPlasH/Vorticity/MicropolarModel_Bender2017.h>
#include <SPlisHSPlasH/Vorticity/VorticityConfinement.h>
#include <SPlisHSPlasH/Vorticity/VorticityRefinement_Liu2021.h>

namespace py = pybind11;

void VorticityModule(py::module m_sub) {

    // ---------------------------------------
    // Vorticity Bender 2017
    // ---------------------------------------
    py::class_<SPH::MicropolarModel_Bender2017, SPH::NonPressureForceBase>(m_sub, "MicropolarModel_Bender2017")
            .def_readwrite_static("VORTICITY_COEFFICIENT", &SPH::MicropolarModel_Bender2017::VORTICITY_COEFFICIENT)
            .def_readwrite_static("VISCOSITY_OMEGA", &SPH::MicropolarModel_Bender2017::VISCOSITY_OMEGA)
            .def_readwrite_static("INERTIA_INVERSE", &SPH::MicropolarModel_Bender2017::INERTIA_INVERSE)

            .def(py::init<SPH::FluidModel*>())
            .def("getAngularAcceleration", (const Vector3r& (SPH::MicropolarModel_Bender2017::*)(const unsigned int)const)&SPH::MicropolarModel_Bender2017::getAngularAcceleration)
            // .def("getAngularAcceleration", (Vector3r& (SPH::MicropolarModel_Bender2017::*)(const unsigned int))&SPH::MicropolarModel_Bender2017::getAngularAcceleration) // TODO: wont work by reference
            .def("setAngularAcceleration", &SPH::MicropolarModel_Bender2017::setAngularAcceleration)
            .def("getAngularVelocity", (const Vector3r& (SPH::MicropolarModel_Bender2017::*)(const unsigned int)const)&SPH::MicropolarModel_Bender2017::getAngularVelocity)
            // .def("getAngularVelocity", (Vector3r& (SPH::MicropolarModel_Bender2017::*)(const unsigned int))&SPH::MicropolarModel_Bender2017::getAngularVelocity) // TODO: wont work by reference
            .def("setAngularVelocity", &SPH::MicropolarModel_Bender2017::setAngularVelocity)
            .def("getMethodName", &SPH::MicropolarModel_Bender2017::getMethodName);

    // ---------------------------------------
    // Vorticity confinement
    // ---------------------------------------
    py::class_<SPH::VorticityConfinement, SPH::NonPressureForceBase>(m_sub, "VorticityConfinement")
            .def_readwrite_static("VORTICITY_COEFFICIENT", &SPH::VorticityConfinement::VORTICITY_COEFFICIENT)
            .def(py::init<SPH::FluidModel*>())
            .def("getAngularVelocity", (const Vector3r& (SPH::VorticityConfinement::*)(const unsigned int)const)&SPH::VorticityConfinement::getAngularVelocity)
            // .def("getAngularVelocity", (Vector3r& (SPH::VorticityConfinement::*)(const unsigned int))&SPH::VorticityConfinement::getAngularVelocity) // TODO: wont work by reference
            .def("setAngularVelocity", &SPH::VorticityConfinement::setAngularVelocity)
            .def("getMethodName", &SPH::VorticityConfinement::getMethodName);

    // ---------------------------------------
    // Liu et al. 2021
    // ---------------------------------------
    py::class_<SPH::VorticityRefinement_Liu2021, SPH::NonPressureForceBase>(m_sub, "VorticityRefinement_Liu2021")
            .def_readwrite_static("VORTICITY_COEFFICIENT", &SPH::VorticityRefinement_Liu2021::VORTICITY_COEFFICIENT)
            .def_readwrite_static("KINEMATIC_VISCOSITY", &SPH::VorticityRefinement_Liu2021::KINEMATIC_VISCOSITY)
            .def(py::init<SPH::FluidModel*>())
            .def("getVorticity", (const Vector3r& (SPH::VorticityRefinement_Liu2021::*)(const unsigned int)const)&SPH::VorticityRefinement_Liu2021::getVorticity)
            // .def("getVorticity", (Vector3r& (SPH::VorticityRefinement_Liu2021::*)(const unsigned int))&SPH::VorticityRefinement_Liu2021::getVorticity) // TODO: wont work by reference
            .def("setVorticity", &SPH::VorticityRefinement_Liu2021::setVorticity)
            .def("getDissipatedVorticity", (const Vector3r& (SPH::VorticityRefinement_Liu2021::*)(const unsigned int)const)&SPH::VorticityRefinement_Liu2021::getDissipatedVorticity)
            // .def("getDissipatedVorticity", (Vector3r& (SPH::VorticityRefinement_Liu2021::*)(const unsigned int))&SPH::VorticityRefinement_Liu2021::getDissipatedVorticity) // TODO: wont work by reference
            .def("setDissipatedVorticity", &SPH::VorticityRefinement_Liu2021::setDissipatedVorticity)
            .def("getMethodName", &SPH::VorticityRefinement_Liu2021::getMethodName);
}