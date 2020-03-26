//
// Created by sjeske on 1/23/20.
//
#include "common.h"

#include <SPlisHSPlasH/BoundaryModel.h>
#include <SPlisHSPlasH/BoundaryModel_Akinci2012.h>
#include <SPlisHSPlasH/BoundaryModel_Bender2019.h>
#include <SPlisHSPlasH/BoundaryModel_Koschier2017.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void BoundaryModelModule(py::module m_sub){
    // ---------------------------------------
    // Boundary Model
    // ---------------------------------------
    py::class_<SPH::BoundaryModel>(m_sub, "BoundaryModel")
            .def(py::init<>())
            .def("reset", &SPH::BoundaryModel::reset)
            .def("performNeighborhoodSearchSort", &SPH::BoundaryModel::performNeighborhoodSearchSort)
            .def("saveState", &SPH::BoundaryModel::saveState)
            .def("loadState", &SPH::BoundaryModel::loadState)
            .def("getRigidBodyObject", &SPH::BoundaryModel::getRigidBodyObject, py::return_value_policy::reference_internal)
            .def("addForce", overload_cast_<const Vector3r&, const Vector3r&>()(&SPH::BoundaryModel::addForce))
#ifdef USE_AVX
            .def("addForce", overload_cast_<const Vector3f8&, const Vector3f8&, const unsigned int>()(&SPH::BoundaryModel::addForce))
#endif
            .def("getPointVelocity", &SPH::BoundaryModel::getPointVelocity) // TODO: remove or fix
            .def("getForceAndTorque", &SPH::BoundaryModel::getForceAndTorque) // TODO: remove or fix
            .def("clearForceAndTorque", &SPH::BoundaryModel::clearForceAndTorque);

    // ---------------------------------------
    // Boundary Model Akinci 2012
    // ---------------------------------------
    py::class_<SPH::BoundaryModel_Akinci2012, SPH::BoundaryModel>(m_sub, "BoundaryModelAkinci2012")
            .def(py::init<>())
            .def("numberOfParticles", &SPH::BoundaryModel_Akinci2012::numberOfParticles)
            .def("getPointSetIndex", &SPH::BoundaryModel_Akinci2012::getPointSetIndex)
            .def("computeBoundaryVolume", &SPH::BoundaryModel_Akinci2012::computeBoundaryVolume)
            .def("resize", &SPH::BoundaryModel_Akinci2012::resize)
            .def("reset", &SPH::BoundaryModel_Akinci2012::reset)
            .def("performNeighborhoodSearchSort", &SPH::BoundaryModel_Akinci2012::performNeighborhoodSearchSort)
            .def("saveState", &SPH::BoundaryModel_Akinci2012::saveState)
            .def("loadState", &SPH::BoundaryModel_Akinci2012::loadState)
            .def("initModel", &SPH::BoundaryModel_Akinci2012::initModel)
            // .def("getPosition0", (Vector3r& (SPH::BoundaryModel_Akinci2012::*)(const unsigned int))(&SPH::BoundaryModel_Akinci2012::getPosition0)) // TODO: wont work by reference
            .def("getPosition0", (const Vector3r& (SPH::BoundaryModel_Akinci2012::*)(const unsigned int) const)(&SPH::BoundaryModel_Akinci2012::getPosition0))
            .def("setPosition0", &SPH::BoundaryModel_Akinci2012::setPosition0)
            // .def("getPosition", (Vector3r& (SPH::BoundaryModel_Akinci2012::*)(const unsigned int))(&SPH::BoundaryModel_Akinci2012::getPosition)) TODO: wont work by reference
            .def("getPosition", (const Vector3r& (SPH::BoundaryModel_Akinci2012::*)(const unsigned int) const)(&SPH::BoundaryModel_Akinci2012::getPosition))
            .def("setPosition", &SPH::BoundaryModel_Akinci2012::setPosition)
            // .def("getVelocity", (Vector3r& (SPH::BoundaryModel_Akinci2012::*)(const unsigned int))(&SPH::BoundaryModel_Akinci2012::getVelocity)) // TODO: wont work by reference
            .def("getVelocity", (const Vector3r& (SPH::BoundaryModel_Akinci2012::*)(const unsigned int) const)(&SPH::BoundaryModel_Akinci2012::getVelocity))
            .def("setVelocity", &SPH::BoundaryModel_Akinci2012::setVelocity)
            // .def("getVolume", (Real& (SPH::BoundaryModel_Akinci2012::*)(const unsigned int))(&SPH::BoundaryModel_Akinci2012::getVolume)) // TODO: might work by reference, but not intended behaviour. Use setter instead
            .def("getVolume", (const Real& (SPH::BoundaryModel_Akinci2012::*)(const unsigned int) const)(&SPH::BoundaryModel_Akinci2012::getVolume))
            .def("setVolume", &SPH::BoundaryModel_Akinci2012::setVolume);

    // ---------------------------------------
    // Boundary Model Bender 2019
    // ---------------------------------------
    py::class_<SPH::BoundaryModel_Bender2019, SPH::BoundaryModel>(m_sub, "BoundaryModelBender2019")
            .def(py::init<>())
            .def("initModel", &SPH::BoundaryModel_Bender2019::initModel)
            .def("reset", &SPH::BoundaryModel_Bender2019::reset)
            .def("getMap", &SPH::BoundaryModel_Bender2019::getMap, py::return_value_policy::reference_internal)
            .def("setMap", &SPH::BoundaryModel_Bender2019::setMap)
            .def("getMaxDist", &SPH::BoundaryModel_Bender2019::getMaxDist)
            .def("setMaxDist", &SPH::BoundaryModel_Bender2019::setMaxDist)
            .def("getMaxVel", &SPH::BoundaryModel_Bender2019::getMaxVel)
            .def("setMaxVel", &SPH::BoundaryModel_Bender2019::setMaxVel)
            .def("getBoundaryVolume", (const Real& (SPH::BoundaryModel_Bender2019::*)(const unsigned int, const unsigned int)const)(&SPH::BoundaryModel_Bender2019::getBoundaryVolume))
            // .def("getBoundaryVolume", (Real& (SPH::BoundaryModel_Bender2019::*)(const unsigned int, const unsigned int))(&SPH::BoundaryModel_Bender2019::getBoundaryVolume)) // TODO: might work by reference, but not intended behaviour. Use setter instead
            .def("setBoundaryVolume", &SPH::BoundaryModel_Bender2019::setBoundaryVolume)
            .def("getBoundaryXj", (const Vector3r& (SPH::BoundaryModel_Bender2019::*)(const unsigned int, const unsigned int)const)(&SPH::BoundaryModel_Bender2019::getBoundaryXj))
            // .def("getBoundaryXj", (Vector3r& (SPH::BoundaryModel_Bender2019::*)(const unsigned int, const unsigned int))(&SPH::BoundaryModel_Bender2019::getBoundaryXj)) // TODO: wont work by reference
            .def("setBoundaryXj", &SPH::BoundaryModel_Bender2019::setBoundaryXj);

    // ---------------------------------------
    // Boundary Model Koschier 2017
    // ---------------------------------------
    py::class_<SPH::BoundaryModel_Koschier2017, SPH::BoundaryModel>(m_sub, "BoundaryModelKoschier2017")
            .def(py::init<>())
            .def("initModel", &SPH::BoundaryModel_Koschier2017::initModel)
            .def("reset", &SPH::BoundaryModel_Koschier2017::reset)
            .def("getMap", &SPH::BoundaryModel_Koschier2017::getMap, py::return_value_policy::reference_internal)
            .def("setMap", &SPH::BoundaryModel_Koschier2017::setMap)
            .def("getMaxDist", &SPH::BoundaryModel_Koschier2017::getMaxDist)
            .def("setMaxDist", &SPH::BoundaryModel_Koschier2017::setMaxDist)
            .def("getMaxVel", &SPH::BoundaryModel_Koschier2017::getMaxVel)
            .def("setMaxVel", &SPH::BoundaryModel_Koschier2017::setMaxVel)
            .def("getBoundaryDensity", (const Real & (SPH::BoundaryModel_Koschier2017::*)(const unsigned int, const unsigned int)const)(&SPH::BoundaryModel_Koschier2017::getBoundaryDensity))
            // .def("getBoundaryDensity", (Real & (SPH::BoundaryModel_Koschier2017::*)(const unsigned int, const unsigned int))(&SPH::BoundaryModel_Koschier2017::getBoundaryDensity)) // TODO: wont work by reference
            .def("setBoundaryDensity", &SPH::BoundaryModel_Koschier2017::setBoundaryDensity)
            .def("getBoundaryDensityGradient", (const Vector3r& (SPH::BoundaryModel_Koschier2017::*)(const unsigned int, const unsigned int)const)(&SPH::BoundaryModel_Koschier2017::getBoundaryDensityGradient))
            // .def("getBoundaryDensityGradient", (Vector3r& (SPH::BoundaryModel_Koschier2017::*)(const unsigned int, const unsigned int))(&SPH::BoundaryModel_Koschier2017::getBoundaryDensityGradient)) // TODO: wont work by reference
            .def("setBoundaryDensityGradient", &SPH::BoundaryModel_Koschier2017::setBoundaryDensityGradient)
            .def("getBoundaryXj", (const Vector3r& (SPH::BoundaryModel_Koschier2017::*)(const unsigned int, const unsigned int)const)(&SPH::BoundaryModel_Koschier2017::getBoundaryXj))
            // .def("getBoundaryXj", (Vector3r& (SPH::BoundaryModel_Koschier2017::*)(const unsigned int, const unsigned int))(&SPH::BoundaryModel_Koschier2017::getBoundaryXj)) // TODO: wont work by reference
            .def("setBoundaryXj", &SPH::BoundaryModel_Koschier2017::setBoundaryXj);
}
