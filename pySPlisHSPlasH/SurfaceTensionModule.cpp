//
// Created by sjeske on 1/28/20.
//
#include "common.h"

#include <SPlisHSPlasH/NonPressureForceBase.h>
#include <SPlisHSPlasH/SurfaceTension/SurfaceTensionBase.h>
#include <SPlisHSPlasH/SurfaceTension/SurfaceTension_Akinci2013.h>
#include <SPlisHSPlasH/SurfaceTension/SurfaceTension_Becker2007.h>
#include <SPlisHSPlasH/SurfaceTension/SurfaceTension_He2014.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

void SurfaceTensionModule(py::module m_sub) {
    // ---------------------------------------
    // Class Surface Tension Base
    // ---------------------------------------
    py::class_<SPH::SurfaceTensionBase, SPH::NonPressureForceBase>(m_sub, "SurfaceTensionBase")
            .def_readwrite_static("SURFACE_TENSION", &SPH::SurfaceTensionBase::SURFACE_TENSION)
            .def_readwrite_static("SURFACE_TENSION_BOUNDARY", &SPH::SurfaceTensionBase::SURFACE_TENSION_BOUNDARY);

    // ---------------------------------------
    // Class Surface Akinci 2013
    // ---------------------------------------
    py::class_<SPH::SurfaceTension_Akinci2013, SPH::SurfaceTensionBase>(m_sub, "SurfaceTension_Akinci2013")
            .def(py::init<SPH::FluidModel*>())
            .def("computeNormals", &SPH::SurfaceTension_Akinci2013::computeNormals)
            // .def("getNormal", (Vector3r& (SPH::SurfaceTension_Akinci2013::*)(const unsigned int))&SPH::SurfaceTension_Akinci2013::getNormal) // TODO: wont work by reference
            .def("getNormal", (const Vector3r& (SPH::SurfaceTension_Akinci2013::*)(const unsigned int)const)&SPH::SurfaceTension_Akinci2013::getNormal)
            .def("setNormal", &SPH::SurfaceTension_Akinci2013::setNormal);

    // ---------------------------------------
    // Class Surface Becker 2007
    // ---------------------------------------
    py::class_<SPH::SurfaceTension_Becker2007, SPH::SurfaceTensionBase>(m_sub, "SurfaceTension_Becker2007")
            .def(py::init<SPH::FluidModel*>());

    // ---------------------------------------
    // Class Surface He 2014
    // ---------------------------------------
    py::class_<SPH::SurfaceTension_He2014, SPH::SurfaceTensionBase>(m_sub, "SurfaceTension_He2014")
            .def(py::init<SPH::FluidModel*>())
            .def("getColor", (const Real (SPH::SurfaceTension_He2014::*)(const unsigned int)const)(&SPH::SurfaceTension_He2014::getColor))
            // .def("getColor", (Real& (SPH::SurfaceTension_He2014::*)(const unsigned int))(&SPH::SurfaceTension_He2014::getColor)) // TODO: wont work by reference
            .def("setColor", &SPH::SurfaceTension_He2014::setColor)
            .def("getGradC2", (const Real (SPH::SurfaceTension_He2014::*)(const unsigned int)const)(&SPH::SurfaceTension_He2014::getGradC2))
            // .def("getGradC2", (Real& (SPH::SurfaceTension_He2014::*)(const unsigned int))(&SPH::SurfaceTension_He2014::getGradC2)) // TODO: wont work by reference
            .def("setGradC2", &SPH::SurfaceTension_He2014::setGradC2);
}