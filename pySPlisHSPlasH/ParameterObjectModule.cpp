//
// Created by sjeske on 1/22/20.
//
#include "common.h"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include <pybind11/native_enum.h> 

#include <SPlisHSPlasH/Common.h>
#include <ParameterObject.h>
#include <Parameter.h>

namespace py = pybind11;

void ParameterObjectModule(py::module m_sub){
    //auto m_sub = m.def_submodule("Common");
    auto param_obj = py::class_<GenParam::ParameterObject>(m_sub, "ParameterObject")
        // .def(py::init<>()) TODO: no constructor for now because this object does not need to be constructable
        .def("getType", &GenParam::ParameterObject::getType)
        .def("getName", &GenParam::ParameterObject::getName)
        .def("getLabel", &GenParam::ParameterObject::getLabel)
        .def("getDescription", &GenParam::ParameterObject::getDescription)
        .def("getGroup", &GenParam::ParameterObject::getGroup)
        .def("numParameters", &GenParam::ParameterObject::numParameters)
        .def("getEnumValues", [](GenParam::ParameterObject& obj, const unsigned int id) -> std::vector<GenParam::EnumParameter::EnumValue> {

            if (obj.getType(id) == GenParam::ParameterBase::DataTypes::ENUM)
            {
                auto param = dynamic_cast<GenParam::EnumParameter*>(obj.getParameter(id));
                return param->getEnumValues();
            }
            return {};
        })

        .def("getValueBool", &GenParam::ParameterObject::getValue<bool>)
        .def("getValueInt", &GenParam::ParameterObject::getValue<int>)
        .def("getValueUInt", &GenParam::ParameterObject::getValue<unsigned int>)
        .def("getValueFloat", &GenParam::ParameterObject::getValue<Real>)
        .def("getValueString", &GenParam::ParameterObject::getValue<std::string>)
        .def("getVec3ValueReal", [](GenParam::ParameterObject& obj, const unsigned int id) -> Vector3r { return Vector3r(obj.getVecValue<Real>(id)); })
        .def("getVec4ValueReal", [](GenParam::ParameterObject& obj, const unsigned int id) -> Vector4r { return Vector4r(obj.getVecValue<Real>(id)); })
        .def("getVec3ValueInt", [](GenParam::ParameterObject& obj, const unsigned int id) ->  Eigen::Matrix<int, 3, 1> { return  Eigen::Matrix<int, 3, 1>(obj.getVecValue<int>(id)); })
        .def("getVec3ValueUInt", [](GenParam::ParameterObject& obj, const unsigned int id) ->  Eigen::Matrix<unsigned int, 3, 1> { return  Eigen::Matrix<unsigned int, 3, 1>(obj.getVecValue<unsigned int>(id)); })

        .def("setValueBool", &GenParam::ParameterObject::setValue<bool>)
        .def("setValueInt", &GenParam::ParameterObject::setValue<int>)
        .def("setValueUInt", &GenParam::ParameterObject::setValue<unsigned int>)
        .def("setValueFloat", &GenParam::ParameterObject::setValue<Real>)
        .def("setValueString", &GenParam::ParameterObject::setValue<std::string>)
        .def("setVec3ValueReal", [](GenParam::ParameterObject& obj, const unsigned int id, Vector3r& val) { obj.setVecValue<Real>(id, val.data()); })
        .def("setVec4ValueReal", [](GenParam::ParameterObject& obj, const unsigned int id, Vector4r& val) { obj.setVecValue<Real>(id, val.data()); })
        .def("setVec3ValueUint", [](GenParam::ParameterObject& obj, const unsigned int id, Eigen::Matrix<unsigned int, 3, 1>& val) { obj.setVecValue<unsigned int>(id, val.data()); });

    py::class_<GenParam::EnumParameter::EnumValue>(m_sub, "EnumValue")
        .def_readonly("id", &GenParam::EnumParameter::EnumValue::id)
        .def_readonly("name", &GenParam::EnumParameter::EnumValue::name);

    // auto param_obj = py::class_<GenParam::ParameterBase>(m_sub, "ParameterBase")
    //     .def("getName", &GenParam::ParameterBase::getName)
    //     .def("setName", &GenParam::ParameterBase::setName)
    //     .def("getGroup", &GenParam::ParameterBase::getGroup)
    //     .def("setGroup", &GenParam::ParameterBase::setGroup)
    //     .def("getGroup", &GenParam::ParameterBase::getGroup)
    //     .def("setGroup", &GenParam::ParameterBase::setGroup)

    py::native_enum<GenParam::ParameterBase::DataTypes>(param_obj, "DataTypes", "enum.Enum")
        .value("BOOL", GenParam::ParameterBase::DataTypes::BOOL)
        .value("FLOAT", GenParam::ParameterBase::DataTypes::FLOAT)
        .value("DOUBLE", GenParam::ParameterBase::DataTypes::DOUBLE)
        .value("INT8", GenParam::ParameterBase::DataTypes::INT8)
        .value("INT16", GenParam::ParameterBase::DataTypes::INT16)
        .value("INT32", GenParam::ParameterBase::DataTypes::INT32)
        .value("UINT8", GenParam::ParameterBase::DataTypes::UINT8)
        .value("UINT16", GenParam::ParameterBase::DataTypes::UINT16)
        .value("UINT32", GenParam::ParameterBase::DataTypes::UINT32)
        .value("VEC_FLOAT", GenParam::ParameterBase::DataTypes::VEC_FLOAT)
        .value("VEC_DOUBLE", GenParam::ParameterBase::DataTypes::VEC_DOUBLE)
        .value("VEC_INT8",  GenParam::ParameterBase::DataTypes::VEC_INT8)
        .value("VEC_INT16", GenParam::ParameterBase::DataTypes::VEC_INT16)
        .value("VEC_INT32", GenParam::ParameterBase::DataTypes::VEC_INT32)
        .value("VEC_UINT8",  GenParam::ParameterBase::DataTypes::VEC_UINT8)
        .value("VEC_UINT16", GenParam::ParameterBase::DataTypes::VEC_UINT16)
        .value("VEC_UINT32", GenParam::ParameterBase::DataTypes::VEC_UINT32)
        .value("ENUM", GenParam::ParameterBase::DataTypes::ENUM)
        .value("STRING", GenParam::ParameterBase::DataTypes::STRING)
        .value("FUNCTION", GenParam::ParameterBase::DataTypes::FUNCTION)
        .export_values()
        .finalize();
}
