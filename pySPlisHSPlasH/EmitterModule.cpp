//
// Created by sjeske on 1/24/20.
//
#include "common.h"

#include <SPlisHSPlasH/Emitter.h>
#include <SPlisHSPlasH/EmitterSystem.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>

#include "bind_pointer_vector.h"

namespace py = pybind11;

void EmitterModule(py::module m_sub){
    // ---------------------------------------
    // Emitter Class
    // ---------------------------------------
    py::class_<SPH::Emitter>(m_sub, "Emitter")
            .def(py::init<>([](SPH::FluidModel *model,
                               const unsigned int width, const unsigned int height,
                               const Vector3r &pos, const Matrix3r & rotation,
                               const Real velocity,
                               const unsigned int type = 0){
                return SPH::Emitter(model, width, height, pos, rotation, velocity, type);
            }))
            .def("emitParticles", &SPH::Emitter::emitParticles)
            .def("emitParticlesCircle", &SPH::Emitter::emitParticlesCircle)
            .def("getNextEmitTime", &SPH::Emitter::getNextEmitTime)
            .def("setNextEmitTime", &SPH::Emitter::setNextEmitTime)
            .def("setEmitStartTime", &SPH::Emitter::setEmitStartTime)
            .def("setEmitEndTime", &SPH::Emitter::setEmitEndTime)
            .def("getPosition", &SPH::Emitter::getPosition)
            .def("setPosition", &SPH::Emitter::setPosition)
            .def("getRotation", &SPH::Emitter::getRotation)
            .def("setRotation", &SPH::Emitter::setRotation)
            .def("getVelocity", &SPH::Emitter::getVelocity)
            .def("setVelocity", &SPH::Emitter::setVelocity)
            .def_static("getSize", &SPH::Emitter::getSize)
            .def("step", &SPH::Emitter::step)
            .def("reset", &SPH::Emitter::reset)
            .def("saveState", &SPH::Emitter::saveState)
            .def("loadState", &SPH::Emitter::loadState);

    py::bind_pointer_vector<std::vector<SPH::Emitter*>>(m_sub, "EmitterVector");

    // ---------------------------------------
    // Emitter System Class
    // ---------------------------------------
    py::class_<SPH::EmitterSystem>(m_sub, "EmitterSystem")
            .def(py::init<>([](SPH::FluidModel* model){
                return SPH::EmitterSystem(model);
            }))
            .def("enableReuseParticles", &SPH::EmitterSystem::enableReuseParticles)
            .def("disableReuseParticles", &SPH::EmitterSystem::disableReuseParticles)
            .def("addEmitter", &SPH::EmitterSystem::addEmitter)
            .def("numEmitters", &SPH::EmitterSystem::numEmitters)
            .def("getEmitters", &SPH::EmitterSystem::getEmitters, py::return_value_policy::reference_internal)
            .def("numReusedParticles", &SPH::EmitterSystem::numReusedParticles)
            .def("numEmittedParticles", &SPH::EmitterSystem::numEmittedParticles)
            .def("step", &SPH::EmitterSystem::step)
            .def("reset", &SPH::EmitterSystem::reset)
            .def("saveState", &SPH::EmitterSystem::saveState)
            .def("loadState", &SPH::EmitterSystem::loadState);

}
