#include "common.h"

#include <SPlisHSPlasH/Common.h>
#include <SPlisHSPlasH/AnimationField.h>
#include <SPlisHSPlasH/AnimationFieldSystem.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "bind_pointer_vector.h"

namespace py = pybind11;

void AnimationFieldModule(py::module m) {
	py::class_<SPH::AnimationField>(m, "AnimationField")
		.def(py::init<>([](
			const std::string &particleFieldName,
			const Vector3r &pos, const Matrix3r & rotation, const Vector3r &scale,
			const std::string expression[3], const unsigned int type = 0) {
				return SPH::AnimationField(particleFieldName, pos, rotation, scale, expression, type); }))
		.def("setStartTime", &SPH::AnimationField::setStartTime)
		.def("setEndTime", &SPH::AnimationField::setEndTime)
		.def("step", &SPH::AnimationField::step);

	py::bind_pointer_vector<std::vector<SPH::AnimationField*>>(m, "AnimationFieldVector");

	py::class_<SPH::AnimationFieldSystem>(m, "AnimationFieldSystem")
		.def(py::init<>())
		.def("addAnimationField", &SPH::AnimationFieldSystem::addAnimationField)
		.def("numAnimationFields", &SPH::AnimationFieldSystem::numAnimationFields)
		.def("getAnimationFields", &SPH::AnimationFieldSystem::getAnimationFields, py::return_value_policy::reference_internal)
		.def("step", &SPH::AnimationFieldSystem::step)
		.def("reset", &SPH::AnimationFieldSystem::reset);
}