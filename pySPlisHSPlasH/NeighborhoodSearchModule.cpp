//
// Created by lwesthofen on 12/8/22.
//
#include "common.h"

#include <SPlisHSPlasH/NeighborhoodSearch.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void NeighborhoodSearchModule(py::module m_sub) {
    // ---------------------------------------
    // NeighborhoodSearch Module
    // ---------------------------------------
    py::class_<SPH::NeighborhoodSearchWrapper>(m_sub, "NeighborhoodSearch")
            .def(py::init<const Real>())
            .def("getMapPointSet_UserData", &SPH::NeighborhoodSearchWrapper::getMapPointSet_UserData)
            .def("setNeighborhoodSearchActive", overload_cast_<const bool>()(&SPH::NeighborhoodSearchWrapper::setNeighborhoodSearchActive))
            .def("setNeighborhoodSearchActive", overload_cast_<unsigned int, unsigned int, bool>()(&SPH::NeighborhoodSearchWrapper::setNeighborhoodSearchActive))
            .def("setNeighborhoodSearchActive", overload_cast_<unsigned int, bool, bool>()(&SPH::NeighborhoodSearchWrapper::setNeighborhoodSearchActive))
            //.def("zSort", &SPH::NeighborhoodSearchWrapper::zSort) // TODO: May not work indicated by the particle selection
            .def("applyZSort", &SPH::NeighborhoodSearchWrapper::applyZSort<bool>)
            .def("applyZSort", &SPH::NeighborhoodSearchWrapper::applyZSort<unsigned int>)
            .def("applyZSort", &SPH::NeighborhoodSearchWrapper::applyZSort<Real>)
            .def("applyZSort", &SPH::NeighborhoodSearchWrapper::applyZSort<Vector3r>) // TODO: May need further implementations
            .def("updatePointSets", &SPH::NeighborhoodSearchWrapper::updatePointSets)
            .def("findNeighbors", &SPH::NeighborhoodSearchWrapper::findNeighbors)
            .def("addPointSet", overload_cast_<Real*, std::size_t, bool, bool, bool, void*>()(&SPH::NeighborhoodSearchWrapper::addPointSet))
            .def("resizeSet", &SPH::NeighborhoodSearchWrapper::resizeSet)
            .def("setSearchRadius", &SPH::NeighborhoodSearchWrapper::setSearchRadius)
            .def("numberOfPointsInSet", &SPH::NeighborhoodSearchWrapper::numberOfPointsInSet)
            .def("numberOfPointSets", &SPH::NeighborhoodSearchWrapper::numberOfPointSets)
            .def("numberOfNeighbors", &SPH::NeighborhoodSearchWrapper::numberOfNeighbors)
            .def("getNeighbor", &SPH::NeighborhoodSearchWrapper::getNeighbor)
            .def("getNeighborList", &SPH::NeighborhoodSearchWrapper::getNeighborList);          
}