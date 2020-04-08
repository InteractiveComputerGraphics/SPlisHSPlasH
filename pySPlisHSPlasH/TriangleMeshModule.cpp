#include "common.h"
#include <SPlisHSPlasH/TimeManager.h>

#include <pybind11/pybind11.h>
#include <SPlisHSPlasH/TriangleMesh.h>

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void TriangleMeshModule(py::module m_sub) {
    // ---------------------------------------
    // Class Time Manager
    // ---------------------------------------
    using Faces = SPH::TriangleMesh::Faces;
    using Normals = SPH::TriangleMesh::Normals;
    using Vertices = SPH::TriangleMesh::Vertices;

    py::class_<SPH::TriangleMesh>(m_sub, "TriangleMesh")
            .def(py::init<>())
            .def("release", &SPH::TriangleMesh::release)
            .def("initMesh", &SPH::TriangleMesh::initMesh)
            .def("addFace", overload_cast_<const unsigned int* const>()(&SPH::TriangleMesh::addFace))
            .def("addFace", overload_cast_<const unsigned int* const>()(&SPH::TriangleMesh::addFace))
            .def("addVertex", &SPH::TriangleMesh::addVertex)

            .def("getFaces", (const Faces & (SPH::TriangleMesh::*)()const)(&SPH::TriangleMesh::getFaces))
            // .def("getFaces", (Faces & (SPH::TriangleMesh::*)())(&SPH::TriangleMesh::getFaces)) // TODO: wont work by reference
            .def("getFaceNormals", (const Normals & (SPH::TriangleMesh::*)()const)(&SPH::TriangleMesh::getFaceNormals))
            // .def("getFaceNormals", (Normals & (SPH::TriangleMesh::*)())(&SPH::TriangleMesh::getFaceNormals)) // TODO: wont work by reference
            .def("getVertexNormals", (const Normals & (SPH::TriangleMesh::*)()const)(&SPH::TriangleMesh::getVertexNormals))
            // .def("getVertexNormals", (Normals & (SPH::TriangleMesh::*)())(&SPH::TriangleMesh::getVertexNormals)) // TODO: wont work by reference
            .def("getVertices", (const Vertices & (SPH::TriangleMesh::*)()const)(&SPH::TriangleMesh::getVertices))
            // .def("getVertices", (Vertices & (SPH::TriangleMesh::*)())(&SPH::TriangleMesh::getVertices)) // TODO: wont work by reference

            .def("numVertices", &SPH::TriangleMesh::numVertices)
            .def("numFaces", &SPH::TriangleMesh::numFaces)

            .def("updateNormals", &SPH::TriangleMesh::updateNormals)
            .def("updateVertexNormals", &SPH::TriangleMesh::updateVertexNormals);
}
