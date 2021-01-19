//
// Created by sjeske on 1/22/20.
//
#include "common.h"

#include <SPlisHSPlasH/FluidModel.h>
#include <SPlisHSPlasH/EmitterSystem.h>
#include <SPlisHSPlasH/SurfaceTension/SurfaceTensionBase.h>
#include <SPlisHSPlasH/Viscosity/ViscosityBase.h>
#include <SPlisHSPlasH/Vorticity/VorticityBase.h>
#include <SPlisHSPlasH/Drag/DragBase.h>
#include <SPlisHSPlasH/Elasticity/ElasticityBase.h>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include <string>

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

// TODO: remove reference getters

std::function<void* (const unsigned int)> makeVoidPointerFct(py::array_t<Real, py::array::c_style | py::array::forcecast> arr)
{
	return [arr](const unsigned int i) mutable -> void* {	
		return arr.mutable_unchecked().mutable_data(i);
	};
}

void FluidModelModule(py::module m_sub){
    // ---------------------------------------
    // Enum Field Type
    // ---------------------------------------
    py::enum_<SPH::FieldType>(m_sub, "FieldType")
            .value("Scalar", SPH::FieldType::Scalar)
            .value("Vector3", SPH::FieldType::Vector3)
            .value("Vector6", SPH::FieldType::Vector6)
            .value("Matrix3", SPH::FieldType::Matrix3)
            .value("Matrix6", SPH::FieldType::Matrix6)
            .value("UInt", SPH::FieldType::UInt);

    // ---------------------------------------
    // Struct Field Description
    // ---------------------------------------
    py::class_<SPH::FieldDescription>(m_sub, "FieldDescription")
            .def(py::init<>([](const std::string &n, const SPH::FieldType &t,
                               const std::function<void*(const unsigned int)> &fct, const bool s = false){
                return SPH::FieldDescription(n, t, fct, s);
            }))
            .def_readwrite("name", &SPH::FieldDescription::name)
            .def_readwrite("type", &SPH::FieldDescription::type)
            .def_readwrite("getFct", &SPH::FieldDescription::getFct)
            .def_readwrite("storeData", &SPH::FieldDescription::storeData);

    py::bind_vector<std::vector<SPH::FieldDescription>>(m_sub, "FieldDescriptionVector");

    // ---------------------------------------
    // Enum class Particle State
    // ---------------------------------------
    py::enum_<SPH::ParticleState>(m_sub, "ParticleState")
            .value("Active", SPH::ParticleState::Active)
            .value("AnimatedByEmitter", SPH::ParticleState::AnimatedByEmitter);

	// ---------------------------------------
	// Class Fluid Model
	// ---------------------------------------
	py::class_<SPH::FluidModel, GenParam::ParameterObject>(m_sub, "FluidModel")
		.def_readwrite_static("NUM_PARTICLES", &SPH::FluidModel::NUM_PARTICLES)
		.def_readwrite_static("NUM_REUSED_PARTICLES", &SPH::FluidModel::NUM_REUSED_PARTICLES)
		.def_readwrite_static("DENSITY0", &SPH::FluidModel::DENSITY0)

		.def_readwrite_static("DRAG_METHOD", &SPH::FluidModel::DRAG_METHOD)
		.def_readwrite_static("SURFACE_TENSION_METHOD", &SPH::FluidModel::SURFACE_TENSION_METHOD)
		.def_readwrite_static("VISCOSITY_METHOD", &SPH::FluidModel::VISCOSITY_METHOD)
		.def_readwrite_static("VORTICITY_METHOD", &SPH::FluidModel::VORTICITY_METHOD)
		.def_readwrite_static("ELASTICITY_METHOD", &SPH::FluidModel::ELASTICITY_METHOD)

		.def(py::init<>())
		.def("init", &SPH::FluidModel::init)
		.def("getId", &SPH::FluidModel::getId)
		.def("getDensity0", &SPH::FluidModel::getDensity0)
		.def("setDensity0", &SPH::FluidModel::setDensity0)
		.def("getPointSetIndex", &SPH::FluidModel::getPointSetIndex)
		.def("addField", &SPH::FluidModel::addField)
		.def("getFields", &SPH::FluidModel::getFields, py::return_value_policy::reference_internal) // TODO: Bind return vector?
		.def("getFieldBuffer", [](SPH::FluidModel& obj, const unsigned int i) -> py::memoryview {
		    auto &field = obj.getField(i);
            void * base_ptr = field.getFct(0);
            int num_particles = obj.numParticles();
            switch (field.type){
                case SPH::FieldType::Scalar:
                    return py::memoryview::from_buffer((Real*)base_ptr, {num_particles}, {sizeof(Real)});
                case SPH::FieldType::Vector3:
                    return py::memoryview::from_buffer((Real*)base_ptr, {num_particles, 3}, {sizeof(Real) * 3, sizeof(Real)});
                case SPH::FieldType::UInt:
                    return py::memoryview::from_buffer((unsigned int*)base_ptr, {num_particles}, {sizeof(unsigned int)});
                default:
                    break;
            }
            return py::memoryview(py::buffer_info());
		})
        .def("getFieldBuffer", [](SPH::FluidModel& obj, const std::string& name) -> py::memoryview {
            auto &field = obj.getField(name);
            void * base_ptr = field.getFct(0);
            int num_particles = obj.numParticles();
            switch (field.type){
                case SPH::FieldType::Scalar:
                    return py::memoryview::from_buffer((Real*)base_ptr, {num_particles}, {sizeof(Real)});
                case SPH::FieldType::Vector3:
                    return py::memoryview::from_buffer((Real*)base_ptr, {num_particles, 3}, {sizeof(Real) * 3, sizeof(Real)});
                case SPH::FieldType::UInt:
                    return py::memoryview::from_buffer((unsigned int*)base_ptr, {num_particles}, {sizeof(unsigned int)});
                default:
                    break;
            }
            return py::memoryview(py::buffer_info());
        })
        .def("getField", overload_cast_<const unsigned int>()(&SPH::FluidModel::getField))
		.def("getField", overload_cast_<const unsigned int>()(&SPH::FluidModel::getField))
		.def("getField", overload_cast_<const std::string&>()(&SPH::FluidModel::getField))
		.def("numberOfFields", &SPH::FluidModel::numberOfFields)
		.def("removeFieldByName", &SPH::FluidModel::removeFieldByName)
		.def("setNumActiveParticles", &SPH::FluidModel::setNumActiveParticles)
		.def("numberOfParticles", &SPH::FluidModel::numberOfParticles)
		.def("getEmitterSystem", &SPH::FluidModel::getEmitterSystem, py::return_value_policy::reference_internal)
		.def("reset", &SPH::FluidModel::reset)
		.def("performNeighborhoodSearchSort", &SPH::FluidModel::performNeighborhoodSearchSort)
		.def("initModel", &SPH::FluidModel::initModel)
		.def("numParticles", &SPH::FluidModel::numParticles)
		.def("numActiveParticles", &SPH::FluidModel::numActiveParticles)
		.def("getNumActiveParticles0", &SPH::FluidModel::getNumActiveParticles0)
		.def("setNumActiveParticles0", &SPH::FluidModel::setNumActiveParticles0)
		.def("emittedParticles", &SPH::FluidModel::emittedParticles)

		.def("getSurfaceTensionMethod", &SPH::FluidModel::getSurfaceTensionMethod)
		.def("setSurfaceTensionMethod", overload_cast_<const unsigned int>()(&SPH::FluidModel::setSurfaceTensionMethod))
		.def("setSurfaceTensionMethod", overload_cast_<const std::string&>()(&SPH::FluidModel::setSurfaceTensionMethod))
		.def("getViscosityMethod", &SPH::FluidModel::getViscosityMethod)
		.def("setViscosityMethod", overload_cast_<const unsigned int>()(&SPH::FluidModel::setViscosityMethod))
		.def("setViscosityMethod", overload_cast_<const std::string &>()(&SPH::FluidModel::setViscosityMethod))
		.def("getVorticityMethod", &SPH::FluidModel::getVorticityMethod)
		.def("setVorticityMethod", overload_cast_<const unsigned int>()(&SPH::FluidModel::setVorticityMethod))
		.def("setVorticityMethod", overload_cast_<const std::string &>()(&SPH::FluidModel::setVorticityMethod))
		.def("getDragMethod", &SPH::FluidModel::getDragMethod)
		.def("setDragMethod", overload_cast_<const unsigned int>()(&SPH::FluidModel::setDragMethod))
		.def("setDragMethod", overload_cast_<const std::string &>()(&SPH::FluidModel::setDragMethod))
		.def("getElasticityMethod", &SPH::FluidModel::getElasticityMethod)
		.def("setElasticityMethod", overload_cast_<const unsigned int>()(&SPH::FluidModel::setElasticityMethod))
		.def("setElasticityMethod", overload_cast_<const std::string &>()(&SPH::FluidModel::setElasticityMethod))

		.def("getSurfaceTensionBase", &SPH::FluidModel::getSurfaceTensionBase, py::return_value_policy::reference_internal)
		.def("getViscosityBase", &SPH::FluidModel::getViscosityBase, py::return_value_policy::reference_internal)
		.def("getVorticityBase", &SPH::FluidModel::getVorticityBase, py::return_value_policy::reference_internal)
		.def("getDragBase", &SPH::FluidModel::getDragBase, py::return_value_policy::reference_internal)
		.def("getElasticityBase", &SPH::FluidModel::getElasticityBase, py::return_value_policy::reference_internal)

		.def("setDragMethodChangedCallback", &SPH::FluidModel::setDragMethodChangedCallback)
		.def("setSurfaceMethodChangedCallback", &SPH::FluidModel::setSurfaceMethodChangedCallback)
		.def("setViscosityMethodChangedCallback", &SPH::FluidModel::setViscosityMethodChangedCallback)
		.def("setVorticityMethodChangedCallback", &SPH::FluidModel::setVorticityMethodChangedCallback)
		.def("setElasticityMethodChangedCallback", &SPH::FluidModel::setElasticityMethodChangedCallback)

		.def("computeSurfaceTension", &SPH::FluidModel::computeSurfaceTension)
		.def("computeViscosity", &SPH::FluidModel::computeViscosity)
		.def("computeVorticity", &SPH::FluidModel::computeVorticity)
		.def("computeDragForce", &SPH::FluidModel::computeDragForce)
		.def("computeElasticity", &SPH::FluidModel::computeElasticity)

		.def("saveState", &SPH::FluidModel::saveState)
		.def("loadState", &SPH::FluidModel::loadState)

		// .def("getPosition0", (Vector3r& (SPH::FluidModel::*)(const unsigned int))(&SPH::FluidModel::getPosition0)) // TODO: wont work by reference
		.def("getPosition0", (const Vector3r& (SPH::FluidModel::*)(const unsigned int)const)(&SPH::FluidModel::getPosition0))
		.def("setPosition0", &SPH::FluidModel::setPosition0)

		// .def("getPosition", (Vector3r& (SPH::FluidModel::*)(const unsigned int))(&SPH::FluidModel::getPosition)) // TODO: wont work by reference
		.def("getPosition", (const Vector3r& (SPH::FluidModel::*)(const unsigned int)const)(&SPH::FluidModel::getPosition))
		.def("setPosition", &SPH::FluidModel::setPosition)

		// .def("getVelocity", (Vector3r& (SPH::FluidModel::*)(const unsigned int))(&SPH::FluidModel::getVelocity)) // TODO: wont work by reference
		.def("getVelocity", (const Vector3r& (SPH::FluidModel::*)(const unsigned int)const)(&SPH::FluidModel::getVelocity))
		.def("setVelocity", &SPH::FluidModel::setVelocity)

		// .def("getVelocity0", (Vector3r& (SPH::FluidModel::*)(const unsigned int))(&SPH::FluidModel::getVelocity0)) // TODO: wont work by reference
		.def("getVelocity0", (const Vector3r& (SPH::FluidModel::*)(const unsigned int)const)(&SPH::FluidModel::getVelocity0))
		.def("setVelocity0", &SPH::FluidModel::setVelocity0)

		// .def("getAcceleration", (Vector3r& (SPH::FluidModel::*)(const unsigned int))(&SPH::FluidModel::getAcceleration)) // TODO: wont work by reference
		.def("getAcceleration", (const Vector3r& (SPH::FluidModel::*)(const unsigned int)const)(&SPH::FluidModel::getAcceleration))
		.def("setAcceleration", &SPH::FluidModel::setAcceleration)

		// .def("getMass", (Real& (SPH::FluidModel::*)(const unsigned int))(&SPH::FluidModel::getMass)) // TODO: wont work by reference
		.def("getMass", (const Real (SPH::FluidModel::*)(const unsigned int)const)(&SPH::FluidModel::getMass))
		.def("setMass", &SPH::FluidModel::setMass)

		// .def("getDensity", (Real& (SPH::FluidModel::*)(const unsigned int))(&SPH::FluidModel::getDensity)) // TODO: wont work by reference
		.def("getDensity", (const Real&(SPH::FluidModel::*)(const unsigned int)const)(&SPH::FluidModel::getDensity))
		.def("setMass", &SPH::FluidModel::setMass)

		// .def("getParticleId", (unsigned int& (SPH::FluidModel::*)(const unsigned int))(&SPH::FluidModel::getParticleId)) // TODO: wont work by reference
		.def("getParticleId", (const unsigned int&(SPH::FluidModel::*)(const unsigned int)const)(&SPH::FluidModel::getParticleId))

		// .def("getParticleState", (SPH::ParticleState& (SPH::FluidModel::*)(const unsigned int))(&SPH::FluidModel::getParticleState)) // TODO: wont work by reference
		.def("getParticleState", (const SPH::ParticleState&(SPH::FluidModel::*)(const unsigned int)const)(&SPH::FluidModel::getParticleState))
		.def("setParticleState", &SPH::FluidModel::setParticleState)

		// .def("getVolume", (Real& (SPH::FluidModel::*)(const unsigned int))(&SPH::FluidModel::getVolume)) // TODO: wont work by reference
		.def("getVolume", (const Real (SPH::FluidModel::*)(const unsigned int)const)(&SPH::FluidModel::getVolume));

		m_sub.def("makeVoidPointerFct", &makeVoidPointerFct, py::arg().noconvert());
}
