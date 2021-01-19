//
// Created by stefan on 13.01.21.
//

#include "common.h"
#include <pybind11/pybind11.h>
#include <Simulator/PositionBasedDynamicsWrapper/PBDBoundarySimulator.h>
#include <Simulator/PositionBasedDynamicsWrapper/PBDRigidBody.h>
#include <Simulator/PositionBasedDynamicsWrapper/PBDWrapper.h>

namespace py = pybind11;

void PBDModule(py::module m) {
    auto m_sub = m.def_submodule("PBD");

    py::class_<PBDWrapper>(m_sub, "PBDWrapper")
            .def(py::init<>())
            .def("reset", &PBDWrapper::reset)
            .def("initModel", &PBDWrapper::initModel)
            .def("readScene", &PBDWrapper::readScene)
            .def("initTriangleModelConstraints", &PBDWrapper::initTriangleModelConstraints)
            .def("initTetModelConstraints", &PBDWrapper::initTetModelConstraints)
            .def("timeStep", &PBDWrapper::timeStep)
            .def("updateVisModels", &PBDWrapper::updateVisModels)
            .def("loadObj", &PBDWrapper::loadObj)
            // .def("getSimulationModel", &PBDWrapper::getSimulationModel)  //TODO: make this work
            // .def("getCollisionDetection", &PBDWrapper::getCollisionDetection)  //TODO: make this work
            // .def("getTimeStepController", &PBDWrapper::getTimeStepController)  //TODO: make this work
            .def("getDampingCoeff", &PBDWrapper::getDampingCoeff)
            .def("setDampingCoeff", &PBDWrapper::setDampingCoeff)
            .def("getClothSimulationMethod", &PBDWrapper::getClothSimulationMethod)
            .def("setClothSimulationMethod", &PBDWrapper::setClothSimulationMethod)
            .def("getSolidSimulationMethod", &PBDWrapper::getSolidSimulationMethod)
            .def("setSolidSimulationMethod", &PBDWrapper::setSolidSimulationMethod)
            .def("getBendingMethod", &PBDWrapper::getBendingMethod)
            .def("setBendingMethod", &PBDWrapper::setBendingMethod);

    py::class_<SPH::PBDBoundarySimulator, SPH::BoundarySimulator>(m_sub, "PBDBoundarySimulator")
            .def(py::init<SPH::SimulatorBase*>());
            //.def("getPBDWrapper", &SPH::PBDBoundarySimulator::getPBDWrapper, py::return_value_policy::reference_internal);  //TODO: make this work

    py::class_<SPH::PBDRigidBody, SPH::RigidBodyObject>(m_sub, "PBDRigidBody")
            .def(py::init<PBD::RigidBody*>())
            .def("getVertices", &SPH::PBDRigidBody::getVertices)
            .def("getFaces", &SPH::PBDRigidBody::getFaces)
            .def("getVertexBuffer", [](SPH::PBDRigidBody &obj) -> py::memoryview {
                const std::vector<Vector3r>& vertices = obj.getVertices();
                void* base_ptr = const_cast<Real*>(&vertices[0][0]);
                int num_vert = vertices.size();
                return py::memoryview::from_buffer((Real*)base_ptr, { num_vert, 3 }, { sizeof(Real) * 3, sizeof(Real) }, true);
             })
            .def("getFaceBuffer", [](SPH::PBDRigidBody &obj) -> py::memoryview {
                const std::vector<unsigned int>& faces = obj.getFaces();
                unsigned int* base_ptr = const_cast<unsigned int*>(&faces[0]);
                return py::memoryview::from_buffer(base_ptr, { (int) faces.size() }, { sizeof(unsigned int) }, true);
             });

}
