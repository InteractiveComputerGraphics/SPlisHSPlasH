//
// Created by stefan on 14.03.21.
//
#include <pybind11/pybind11.h>
#include <Partio.h>
#include <PartioAttribute.h>
#include <memory>

// Simple custom holder that works like unique_ptr
template <typename T>
class custom_ptr {
    T* impl;
public:
    custom_ptr(T* p) { impl = p; }
    T* get() const { return impl; }
    //T* release_ptr() { }
};

PYBIND11_DECLARE_HOLDER_TYPE(T, custom_ptr<T>, true);

namespace py = pybind11;

PYBIND11_MODULE(partio, m){
    m.def("read", [](const char* filename, const bool verbose){
            return Partio::read(filename, verbose);
    }, py::arg("filename"), py::arg("verbose")=true);
    m.def("readHeaders", [](const char* filename, const bool verbose){
            return Partio::readHeaders(filename, verbose);
    }, py::arg("filename"), py::arg("verbose")=true);
    m.def("write", [](const char* filename, const Partio::ParticlesData& obj, const bool forceCompressed, const bool verbose){
            Partio::write(filename, obj, forceCompressed, verbose);
    }, py::arg("filename"), py::arg("particlesData"), py::arg("forceCompressed")=false, py::arg("verbose")=true);
    m.def("create", &Partio::create);
    m.def("createInterleave", &Partio::createInterleave);
    m.def("cloneSchema", &Partio::cloneSchema);
    m.def("clone", &Partio::clone);

    py::class_<Partio::ParticlesInfo, custom_ptr<Partio::ParticlesInfo>>(m, "ParticlesInfo")
            .def("release", &Partio::ParticlesInfo::release)
            .def("numParticles", &Partio::ParticlesInfo::numParticles)
            .def("numAttributes", &Partio::ParticlesInfo::numAttributes)
            .def("numFixedAttributes", &Partio::ParticlesInfo::numFixedAttributes)
            .def("attributeInfo", [](Partio::ParticlesInfo& obj, const char * name){
                Partio::ParticleAttribute attr;
                obj.attributeInfo(name, attr);
                return attr;
            }, py::return_value_policy::copy)
            .def("fixedAttributeInfo", [](Partio::ParticlesInfo& obj, const char * name){
                Partio::FixedAttribute attr;
                obj.fixedAttributeInfo(name, attr);
                return attr;
            }, py::return_value_policy::copy)
            .def("attributeInfo", [](Partio::ParticlesInfo& obj, const int id){
                Partio::ParticleAttribute attr;
                obj.attributeInfo(id, attr);
                return attr;
            }, py::return_value_policy::copy)

            .def("fixedAttributeInfo", [](Partio::ParticlesInfo& obj, const int id){
                Partio::FixedAttribute attr;
                obj.fixedAttributeInfo(id, attr);
                return attr;
            }, py::return_value_policy::copy);

    py::class_<Partio::ParticlesData, Partio::ParticlesInfo, custom_ptr<Partio::ParticlesData>>(m, "ParticlesData")
            .def("data_buffer", [](Partio::ParticlesData & obj, Partio::ParticleAttribute& attr) -> py::memoryview{
                const int nparticles = obj.numParticles();
                const unsigned char* base_ptr = obj.data<unsigned char>(attr, 0);
                const size_t stride = obj.data<unsigned char>(attr, 1) - base_ptr;
                Partio::ParticleAttributeType t = attr.type;
                switch(t){
                    case Partio::ParticleAttributeType::VECTOR:
                        return py::memoryview::from_buffer(reinterpret_cast<const float*>(base_ptr), {nparticles, 3}, {stride, sizeof(float)});
                    case Partio::ParticleAttributeType::FLOAT:
                        return py::memoryview::from_buffer(reinterpret_cast<const float*>(base_ptr), {nparticles, 1}, {stride, sizeof(float)});
                    case Partio::ParticleAttributeType::INT:
                        return py::memoryview::from_buffer(reinterpret_cast<const int*>(base_ptr), {nparticles, 1}, {stride, sizeof(int)});
                    default: break;
                }

                return py::memoryview(py::buffer_info());
            });

    auto pdm = py::class_<Partio::ParticlesDataMutable, Partio::ParticlesData, custom_ptr<Partio::ParticlesDataMutable>>(m, "ParticlesDataMutable")
            .def("data_buffer_mutable", [](Partio::ParticlesDataMutable & obj, Partio::ParticleAttribute& attr) -> py::memoryview{
                const int nparticles = obj.numParticles();
                unsigned char* base_ptr = obj.dataWrite<unsigned char>(attr, 0);
                const size_t stride = obj.dataWrite<unsigned char>(attr, 1) - base_ptr;
                Partio::ParticleAttributeType t = attr.type;
                switch(t){
                    case Partio::ParticleAttributeType::VECTOR:
                        return py::memoryview::from_buffer(reinterpret_cast<float*>(base_ptr), {nparticles, 3}, {stride, sizeof(float)});
                    case Partio::ParticleAttributeType::FLOAT:
                        return py::memoryview::from_buffer(reinterpret_cast<float*>(base_ptr), {nparticles, 1}, {stride, sizeof(float)});
                    case Partio::ParticleAttributeType::INT:
                        return py::memoryview::from_buffer(reinterpret_cast<int*>(base_ptr), {nparticles, 1}, {stride, sizeof(int)});
                    default: break;
                }

                return py::memoryview(py::buffer_info());
            })
            .def("addAttribute", &Partio::ParticlesDataMutable::addAttribute)
            .def("addParticle", &Partio::ParticlesDataMutable::addParticle)
            .def("addParticles", &Partio::ParticlesDataMutable::addParticles);

    py::enum_<Partio::ParticleAttributeType>(m, "ParticleAttributeType")
            .value("NONE", Partio::ParticleAttributeType::NONE)
            .value("VECTOR", Partio::ParticleAttributeType::VECTOR)
            .value("FLOAT", Partio::ParticleAttributeType::FLOAT)
            .value("INT", Partio::ParticleAttributeType::INT)
            .value("INDEXEDSTR", Partio::ParticleAttributeType::INDEXEDSTR);

    py::class_<Partio::ParticleAttribute>(m, "ParticleAttribute")
            .def_readonly("type", &Partio::ParticleAttribute::type)
            .def_readonly("count", &Partio::ParticleAttribute::count)
            .def_readonly("name", &Partio::ParticleAttribute::name)
            .def_readonly("attributeIndex", &Partio::ParticleAttribute::attributeIndex);

    py::class_<Partio::FixedAttribute>(m, "FixedAttribute")
            .def_readonly("type", &Partio::FixedAttribute::type)
            .def_readonly("count", &Partio::FixedAttribute::count)
            .def_readonly("name", &Partio::FixedAttribute::name)
            .def_readonly("attributeIndex", &Partio::FixedAttribute::attributeIndex);
}