//
// Created by sjeske on 1/22/20.
//
#include "common.h"

#include <SPlisHSPlasH/StaticRigidBody.h>
#include <SPlisHSPlasH/RigidBodyObject.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

void RigidBodyModule(py::module m_sub){
    // ---------------------------------------
    // Class Rigid Body Object
    // ---------------------------------------
    py::class_<SPH::RigidBodyObject>(m_sub, "RigidBodyObject")
            .def("isDynamic", &SPH::RigidBodyObject::isDynamic)
            .def("isAnimated", &SPH::RigidBodyObject::isAnimated)
            .def("setIsAnimated", &SPH::RigidBodyObject::setIsAnimated)
            .def("getMass", &SPH::RigidBodyObject::getMass)
            .def("getPosition", &SPH::RigidBodyObject::getPosition)
            .def("setPosition", &SPH::RigidBodyObject::setPosition)
            .def("getWorldSpacePosition", &SPH::RigidBodyObject::getWorldSpacePosition)
            .def("getVelocity", &SPH::RigidBodyObject::getVelocity)
            .def("setVelocity", &SPH::RigidBodyObject::setVelocity)
            .def("getRotation", [](SPH::RigidBodyObject& obj) -> Vector4r {
                return obj.getRotation().coeffs();
            })
            .def("setRotation", [](SPH::RigidBodyObject& obj, const Vector4r &qVec) {
                    Quaternionr q;
                    q.coeffs() = qVec;
                    obj.setRotation(q);
                })
            .def("getWorldSpaceRotation", &SPH::RigidBodyObject::getWorldSpaceRotation)
            .def("getAngularVelocity", &SPH::RigidBodyObject::getAngularVelocity)
            .def("setAngularVelocity", &SPH::RigidBodyObject::setAngularVelocity)
            .def("addForce", &SPH::RigidBodyObject::addForce)
            .def("addTorque", &SPH::RigidBodyObject::addTorque)
            .def("getVertices", &SPH::RigidBodyObject::getVertices)
            .def("getVertexNormals", &SPH::RigidBodyObject::getVertexNormals)
            .def("getFaces", &SPH::RigidBodyObject::getFaces)
            .def("updateMeshTransformation", &SPH::RigidBodyObject::updateMeshTransformation);

    // ---------------------------------------
    // Class Static Rigid Body
    // ---------------------------------------
    py::class_<SPH::StaticRigidBody, SPH::RigidBodyObject>(m_sub, "StaticRigidBody")
            .def(py::init<>())
            .def("setWorldSpacePosition", &SPH::StaticRigidBody::setWorldSpacePosition)
            .def("setWorldSpaceRotation", &SPH::StaticRigidBody::setWorldSpaceRotation)
            .def("animate", &SPH::StaticRigidBody::animate)
            .def("getGeometry", &SPH::StaticRigidBody::getGeometry)
            .def("getVertexBuffer", [](SPH::StaticRigidBody &obj) -> py::memoryview {
                auto vertices = obj.getVertices();
                void* base_ptr = &vertices[0][0];
                int num_vert = vertices.size();
                return py::memoryview::from_buffer((Real*)base_ptr, { num_vert, 3 }, { sizeof(Real) * 3, sizeof(Real) });
             })
            .def("getFaceBuffer", [](SPH::StaticRigidBody &obj) -> py::memoryview {
                    auto faces = obj.getFaces();
                    unsigned int* base_ptr = faces.data();
                    return py::memoryview::from_buffer(base_ptr, { faces.size() }, { sizeof(unsigned int) });
             });

}
