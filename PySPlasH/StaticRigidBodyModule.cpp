//
// Created by sjeske on 1/22/20.
//
#include "common.h"

#include <SPlisHSPlasH/StaticRigidBody.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

void StaticRigidBodyModule(py::module m_sub){
    // ---------------------------------------
    // Class Static Rigid Body
    // ---------------------------------------
    py::class_<SPH::StaticRigidBody>(m_sub, "StaticRigidBody")
            .def(py::init<>([](){
                return SPH::StaticRigidBody();
            }))
            .def("isDynamic", &SPH::StaticRigidBody::isDynamic)
            .def("getMass", &SPH::StaticRigidBody::getMass)
            .def("getPosition", &SPH::StaticRigidBody::getPosition)
            .def("setPosition", &SPH::StaticRigidBody::setPosition)
            .def("getWorldSpacePosition", &SPH::StaticRigidBody::getWorldSpacePosition)
            .def("getVelocity", &SPH::StaticRigidBody::getVelocity)
            .def("setVelocity", &SPH::StaticRigidBody::setVelocity)
            .def("getRotation", &SPH::StaticRigidBody::getRotation)
            .def("setRotation", &SPH::StaticRigidBody::setRotation)
            .def("getWorldSpaceRotation", &SPH::StaticRigidBody::getWorldSpaceRotation)
            .def("getAngularVelocity", &SPH::StaticRigidBody::getAngularVelocity)
            .def("setAngularVelocity", &SPH::StaticRigidBody::setAngularVelocity)
            .def("addForce", &SPH::StaticRigidBody::addForce)
            .def("addTorque", &SPH::StaticRigidBody::addTorque)

            .def("setWorldSpacePosition", &SPH::StaticRigidBody::setWorldSpacePosition)
            .def("setWorldSpaceRotation", &SPH::StaticRigidBody::setWorldSpaceRotation)
            .def("getGeometry", &SPH::StaticRigidBody::getGeometry);

}
