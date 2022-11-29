#include "common.h"

#include <pybind11/pybind11.h>

#include <SPlisHSPlasH/NonPressureForceBase.h>
#include <SPlisHSPlasH/XSPH.h>

namespace py = pybind11;

void XSPHModule(py::module m_sub) {
    // ---------------------------------------
    // XSPH
    // ---------------------------------------
    py::class_<SPH::XSPH, SPH::NonPressureForceBase>(m_sub, "XSPH")
        .def_readwrite_static("FLUID_COEFFICIENT", &SPH::XSPH::FLUID_COEFFICIENT)
        .def_readwrite_static("BOUNDARY_COEFFICIENT", &SPH::XSPH::BOUNDARY_COEFFICIENT)

        .def(py::init<SPH::FluidModel*>());
}