#include "common.h"

#include <Simulator/Exporter/ExporterBase.h>
#include <Simulator/SimulatorBase.h>

#include <pybind11/pybind11.h>


#include <iostream>

namespace py = pybind11;

class pyExporterBase : public SPH::ExporterBase
{
    using SPH::ExporterBase::ExporterBase;

    virtual void init(const std::string& outputPath) { PYBIND11_OVERRIDE(void, SPH::ExporterBase, init, outputPath); };
    virtual void reset() { PYBIND11_OVERRIDE(void, SPH::ExporterBase, reset, ); };
    virtual void step(const unsigned int frame) { PYBIND11_OVERRIDE_PURE(void, SPH::ExporterBase, step, frame); };
    virtual void setActive(const bool active) { PYBIND11_OVERRIDE(void, SPH::ExporterBase, setActive, active); }
};

void ExporterModule(py::module m_sub) {
    // ---------------------------------------
    // Exporter Base
    // ---------------------------------------
    py::class_<SPH::ExporterBase, pyExporterBase>(m_sub, "ExporterBase")
        .def(py::init<SPH::SimulatorBase*>())
        .def("step", &SPH::ExporterBase::step)
        .def("reset", &SPH::ExporterBase::reset)
        .def("setActive", &SPH::ExporterBase::setActive)
        .def("getActive", &SPH::ExporterBase::getActive)
        .def("init", &SPH::ExporterBase::init);

}
