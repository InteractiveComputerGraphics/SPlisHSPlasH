//
// Created by sjeske on 1/22/20.
//
#include "common.h"

#include <SPlisHSPlasH/Drag/DragForce_Gissler2017.h>

#include <pybind11/pybind11.h>
#include <SPlisHSPlasH/Drag/DragForce_Macklin2014.h>

namespace py = pybind11;

void DragModule(py::module m_sub) {

    // ---------------------------------------
    // Class Drag Force Gissler
    // ---------------------------------------
    py::class_<SPH::DragForce_Gissler2017, SPH::NonPressureForceBase>(m_sub, "DragForce_Gissler2017")
            .def(py::init<SPH::FluidModel*>())
            .def_readwrite_static("DRAG_COEFFICIENT", &SPH::DragForce_Gissler2017::DRAG_COEFFICIENT)
            .def("getMethodName", &SPH::DragForce_Gissler2017::getMethodName);

    // ---------------------------------------
    // Class Drag Force Macklin
    // ---------------------------------------
    py::class_<SPH::DragForce_Macklin2014, SPH::NonPressureForceBase>(m_sub, "DragForce_Macklin2014")
            .def(py::init<SPH::FluidModel*>())
            .def_readwrite_static("DRAG_COEFFICIENT", &SPH::DragForce_Macklin2014::DRAG_COEFFICIENT)
            .def("getMethodName", &SPH::DragForce_Macklin2014::getMethodName);
}
