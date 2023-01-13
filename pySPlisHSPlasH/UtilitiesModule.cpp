//
// Created by sjeske on 1/22/20.
//
#include "common.h"

#include <vector>
#include <array>
#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/chrono.h>
#include <Utilities/BinaryFileReaderWriter.h>
#include <Utilities/Counting.h>
#include <Utilities/FileSystem.h>
#include <Utilities/OBJLoader.h>
#include <Utilities/PLYLoader.h>
#include <Utilities/PartioReaderWriter.h>
#include <SPlisHSPlasH/Utilities/GaussQuadrature.h>
#include <SPlisHSPlasH/Utilities/MathFunctions.h>
#include <SPlisHSPlasH/Utilities/MatrixFreeSolver.h>
#include <SPlisHSPlasH/Utilities/PoissonDiskSampling.h>
#include <SPlisHSPlasH/Utilities/RegularSampling2D.h>
#include <SPlisHSPlasH/Utilities/RegularTriangleSampling.h>
#include <SPlisHSPlasH/Utilities/SceneLoader.h>
#include <SPlisHSPlasH/Utilities/SDFFunctions.h>
#include <SPlisHSPlasH/Utilities/SimpleQuadrature.h>
#include <SPlisHSPlasH/Utilities/SurfaceSampling.h>
#include <SPlisHSPlasH/Utilities/VolumeSampling.h>
#include <SPlisHSPlasH/Utilities/WindingNumbers.h>
#include <SPlisHSPlasH/Utilities/MeshImport.h>

#include "bind_pointer_vector.h"


template<typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;


namespace py = pybind11;

using namespace pybind11::literals;

void UtilitiesModule(py::module m) {
    // Utilities submodule
    auto m_sub = m.def_submodule("Utilities");

    // ---------------------------------------
    // Binary File reader and writer
    // ---------------------------------------
    py::class_<SPH::BinaryFileWriter>(m_sub, "BinaryFileWriter")
            .def(py::init<>())
            .def("openFile", &SPH::BinaryFileWriter::openFile)
            .def("closeFile", &SPH::BinaryFileWriter::closeFile);

    py::class_<SPH::BinaryFileReader>(m_sub, "BinaryFileReader")
            .def(py::init<>())
            .def("openFile", &SPH::BinaryFileReader::openFile)
            .def("closeFile", &SPH::BinaryFileReader::closeFile);

    // ---------------------------------------
    // Counting class
    // ---------------------------------------
    py::class_<Utilities::AverageCount>(m_sub, "AverageCount")
            .def(py::init<>())
            .def_readwrite("sum", &Utilities::AverageCount::sum)
            .def_readwrite("numberOfCalls", &Utilities::AverageCount::numberOfCalls);

    py::class_<Utilities::Counting>(m_sub, "Counting")
            .def(py::init<>())
            .def_readwrite_static("m_averageCounts", &Utilities::Counting::m_averageCounts)
            .def_static("reset", &Utilities::Counting::reset)
            .def_static("increaseCounter", &Utilities::Counting::increaseCounter)
            .def_static("printAverageCounts", &Utilities::Counting::printAverageCounts)
            .def_static("printCounterSums", &Utilities::Counting::printCounterSums);

    // TODO: check if init counting need to be implemented
    m_sub.def("INCREASE_COUNTER", &Utilities::Counting::increaseCounter);

    // ---------------------------------------
    // File System utilities
    // ---------------------------------------
    py::class_<Utilities::FileSystem>(m_sub, "FileSystem")
            .def(py::init<>())
            .def_static("getFilePath", Utilities::FileSystem::getFilePath)
            .def_static("getFileName", Utilities::FileSystem::getFileName)
            .def_static("getFileNameWithExt", Utilities::FileSystem::getFileNameWithExt)
            .def_static("getFileExt", Utilities::FileSystem::getFileExt)
            .def_static("isRelativePath", Utilities::FileSystem::isRelativePath)
            .def_static("makeDir", Utilities::FileSystem::makeDir)
            .def_static("makeDirs", Utilities::FileSystem::makeDirs)
            .def_static("normalizePath", Utilities::FileSystem::normalizePath)
            .def_static("fileExists", Utilities::FileSystem::fileExists)
            .def_static("getProgramPath", Utilities::FileSystem::getProgramPath)
            .def_static("copyFile", Utilities::FileSystem::copyFile)
            .def_static("isFile", Utilities::FileSystem::isFile)
            .def_static("isDirectory", Utilities::FileSystem::isDirectory)
            .def_static("getFilesInDirectory", Utilities::FileSystem::getFilesInDirectory)
            .def_static("getFileMD5", Utilities::FileSystem::getFileMD5)
            .def_static("writeMD5File", Utilities::FileSystem::writeMD5File)
            .def_static("checkMD5", Utilities::FileSystem::checkMD5);

    // ---------------------------------------
    // Logger
    // ---------------------------------------
    py::enum_<Utilities::LogLevel>(m_sub, "LogLevel")
            .value("DEBUG", Utilities::LogLevel::DEBUG)
            .value("INFO", Utilities::LogLevel::INFO)
            .value("WARN", Utilities::LogLevel::WARN)
            .value("ERR", Utilities::LogLevel::ERR);

    py::class_<Utilities::ConsoleSink>(m_sub, "ConsoleSink")
            .def(py::init<>([](const Utilities::LogLevel minLevel) {
                return Utilities::ConsoleSink(minLevel);
            }))
            .def("write", &Utilities::ConsoleSink::write);

    // TODO: check if it is okay to use shared pointer and implement the actual logger functions
    py::class_<Utilities::FileSink, std::shared_ptr<Utilities::FileSink>>(m_sub, "FileSink")
            .def(py::init<>([](const Utilities::LogLevel minLevel, const std::string &fileName) {
                return std::make_shared<Utilities::FileSink>(minLevel, fileName);
            }))
            .def("write", &Utilities::FileSink::write);

    // ---------------------------------------
    // Object loader necessary vector types
    // ---------------------------------------
    py::bind_vector<std::vector<std::array<float, 3>>>(m_sub, "VecVec3f");
    py::bind_vector<std::vector<std::array<float, 2>>>(m_sub, "VecVec2f");
    py::bind_vector<std::vector<Utilities::MeshFaceIndices>>(m_sub, "VecMeshFaceIndices");

    // Todo: Bind missing attributes
    py::class_<Utilities::MeshFaceIndices>(m_sub, "MeshFaceIndices")
            .def(py::init<>())
            .def("__repr__", [](const Utilities::MeshFaceIndices &m) {
                std::stringstream s;
                s << "posIndices\n" << m.posIndices[0] << " " << m.posIndices[1] << " " << m.posIndices[2] << "\n";
                s << "texIndices\n" << m.texIndices[0] << " " << m.texIndices[1] << " " << m.texIndices[2] << "\n";
                s << "normalIndices\n" << m.normalIndices[0] << " " << m.normalIndices[1] << " " << m.normalIndices[2]
                  << "\n";
                return s.str();
            });

    py::class_<Utilities::OBJLoader>(m_sub, "OBJLoader")
            .def(py::init<>())
            .def_static("loadObj", &Utilities::OBJLoader::loadObj);

    py::class_<Utilities::PLYLoader>(m_sub, "PLYLoader")
        .def(py::init<>())
        .def_static("loadPly", &Utilities::PLYLoader::loadPly);

    // ---------------------------------------
    // Partio reader and writer
    // ---------------------------------------
    py::bind_vector<std::vector<Vector3r>>(m_sub, "VecVector3r", "std::vector<Vector3r> not to be confused with"
                                                                 "VecVec3f which ist std::vector<std::array<float,3>>");

    py::class_<Utilities::PartioReaderWriter>(m_sub, "PartioReaderWriter")
            .def(py::init<>())
            .def_static("readParticle", overload_cast_<const std::string &, const Vector3r &, const Matrix3r &,
                    const Real, std::vector<Vector3r> &>()(&Utilities::PartioReaderWriter::readParticles))
            .def_static("readParticle", overload_cast_<const std::string &, const Vector3r &, const Matrix3r &,
                    const Real, std::vector<Vector3r> &, std::vector<Vector3r> &>()(
                    &Utilities::PartioReaderWriter::readParticles))
            .def_static("readParticle", overload_cast_<const std::string &, const Vector3r &, const Matrix3r &,
                    const Real, std::vector<Vector3r> &, std::vector<Vector3r> &, Real &>()(
                    &Utilities::PartioReaderWriter::readParticles))
            .def_static("writeParticle", &Utilities::PartioReaderWriter::writeParticles);

    // ---------------------------------------
    // String tools
    // ---------------------------------------
    // TODO: String tools are omitted, because python is basically one big string tool


    // ---------------------------------------
    // System Info
    // ---------------------------------------
    // TODO: System info is also omitted. Will add if explicitly required by other classes

    // ---------------------------------------
    // Timing
    // ---------------------------------------
    // TODO: Timing and find a way for everything to be actually printed
    py::class_<Utilities::TimingHelper>(m_sub, "TimingHelper")
            .def(py::init<>())
            .def_readwrite("start", &Utilities::TimingHelper::start)
            .def_readwrite("name", &Utilities::TimingHelper::name);

    py::class_<Utilities::AverageTime>(m_sub, "AverageTime")
            .def(py::init<>())
            .def_readwrite("totalTime", &Utilities::AverageTime::totalTime)
            .def_readwrite("counter", &Utilities::AverageTime::counter)
            .def_readwrite("name", &Utilities::AverageTime::name);

    py::class_<Utilities::IDFactory>(m_sub, "IDFactory")
            .def(py::init<>())
            .def_static("getId", &Utilities::IDFactory::getId);

    py::class_<Utilities::Timing>(m_sub, "Timing")
            .def(py::init<>())
            .def_readwrite_static("m_dontPrintTimes", &Utilities::Timing::m_dontPrintTimes)
            .def_readwrite_static("m_startCounter", &Utilities::Timing::m_startCounter)
            .def_readwrite_static("m_stopCounter", &Utilities::Timing::m_stopCounter)
            .def_readwrite_static("m_timingStack", &Utilities::Timing::m_timingStack)
            .def_readwrite_static("m_averageTimes", &Utilities::Timing::m_averageTimes)
            .def_static("reset", &Utilities::Timing::reset)
            .def_static("startTiming", &Utilities::Timing::startTiming)
            .def_static("stopTiming", overload_cast_<bool>()(&Utilities::Timing::stopTiming))
            .def_static("stopTiming", overload_cast_<bool, int &>()(&Utilities::Timing::stopTiming))
            .def_static("printAverageTimes", &Utilities::Timing::printAverageTimes)
            .def_static("printTimeSums", &Utilities::Timing::printTimeSums);

    // ---------------------------------------
    // Gauss Quadrature
    // ---------------------------------------
    py::class_<SPH::GaussQuadrature>(m_sub, "GaussQuadrature")
            .def_static("integrate", &SPH::GaussQuadrature::integrate)
            .def_static("exportSamples", &SPH::GaussQuadrature::exportSamples);

    // ---------------------------------------
    // Math Functions
    // ---------------------------------------
    py::class_<SPH::MathFunctions>(m_sub, "MathFunctions")
            .def_static("extractRotation", &SPH::MathFunctions::extractRotation)
            .def_static("pseudoInverse", &SPH::MathFunctions::pseudoInverse)
            .def_static("svdWithInversionHandling", &SPH::MathFunctions::svdWithInversionHandling)
            .def_static("eigenDecomposition", &SPH::MathFunctions::eigenDecomposition)
            .def_static("jacobiRotate", &SPH::MathFunctions::jacobiRotate)
            .def_static("getOrthogonalVectors", &SPH::MathFunctions::getOrthogonalVectors);

    // ---------------------------------------
    // Matrix Replacement TODO: Pretty sure the user doesnt need the matrix replacement solved in the python API
    // ---------------------------------------
    py::class_<SPH::MatrixReplacement>(m_sub, "MatrixReplacement");

    // ---------------------------------------
    // Poisson Disk Sampling
    // ---------------------------------------
    py::class_<SPH::PoissonDiskSampling>(m_sub, "PoissonDiskSampling")
            .def(py::init<>())
            .def_static("floor", &SPH::PoissonDiskSampling::floor)
            .def("sampleMesh", &SPH::PoissonDiskSampling::sampleMesh);

    using CellPos = Eigen::Matrix<int, 3, 1, Eigen::DontAlign>;
    py::class_<SPH::PoissonDiskSampling::InitialPointInfo>(m_sub, "InitialPointInfo")
            .def(py::init<>())
            .def(py::init<CellPos, Vector3r, unsigned int>())
            .def_readwrite("cP", &SPH::PoissonDiskSampling::InitialPointInfo::cP)
            .def_readwrite("pos", &SPH::PoissonDiskSampling::InitialPointInfo::pos)
            .def_readwrite("ID", &SPH::PoissonDiskSampling::InitialPointInfo::ID);

    py::class_<SPH::PoissonDiskSampling::HashEntry>(m_sub, "HashEntry")
            .def(py::init<>())
            .def_readwrite("samples", &SPH::PoissonDiskSampling::HashEntry::samples)
            .def_readwrite("startIndex", &SPH::PoissonDiskSampling::HashEntry::startIndex);

    // ---------------------------------------
    // Regular Sampling 2D
    // ---------------------------------------
    py::class_<SPH::RegularSampling2D>(m_sub, "RegularSampling2D")
            .def_static("sampleMesh", &SPH::RegularSampling2D::sampleMesh);

    // ---------------------------------------
    // Scene Loader
    // ---------------------------------------
    py::class_<Utilities::SceneLoader>(m_sub, "SceneLoader")
            .def(py::init<>())
            .def("readScene", &Utilities::SceneLoader::readScene);

    auto m_sub_sub = m_sub.def_submodule("SceneLoaderStructs");
    using SceneInfo = Utilities::SceneLoader;

    py::class_<Utilities::BoundaryParameterObject, GenParam::ParameterObject>(m_sub_sub, "BoundaryData")
        .def(py::init<>())
        .def(py::init<std::string, std::string, Vector3r, Vector3r, Real, Vector3r, bool, bool,
                        Eigen::Matrix<float, 4, 1, Eigen::DontAlign>, std::string, bool, Real,
                        Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign>, unsigned int, bool>(),
                "samplesFile"_a = "", "meshFile"_a = "", "translation"_a = Vector3r::Zero(),
                "axis"_a = Vector3r(1,0,0), "angle"_a = 0.0, "scale"_a = Vector3r::Ones(),
                "isDynamic"_a = false, "isWall"_a = false,
                "color"_a = Eigen::Vector4f(1.f, 0.f, 0.f, 0.f), 
                "mapFile"_a = "", "mapInvert"_a = false, "mapThickness"_a = 0.0,
                "mapResolution"_a = Eigen::Matrix<unsigned int, 3, 1>(20, 20, 20),
                "samplingMode"_a = 0, "isAnimated"_a = false)
        .def_readwrite("samplesFile", &Utilities::BoundaryParameterObject::samplesFile)
        .def_readwrite("meshFile", &Utilities::BoundaryParameterObject::meshFile)
        .def_readwrite("translation", &Utilities::BoundaryParameterObject::translation)
        .def_readwrite("axis", &Utilities::BoundaryParameterObject::axis)
        .def_readwrite("angle", &Utilities::BoundaryParameterObject::angle)
        .def_readwrite("scale", &Utilities::BoundaryParameterObject::scale)
        .def_readwrite("dynamic", &Utilities::BoundaryParameterObject::dynamic)
        .def_readwrite("isWall", &Utilities::BoundaryParameterObject::isWall)
        .def_readwrite("color", &Utilities::BoundaryParameterObject::color)
        .def_readwrite("mapFile", &Utilities::BoundaryParameterObject::mapFile)
        .def_readwrite("mapInvert", &Utilities::BoundaryParameterObject::mapInvert)
        .def_readwrite("mapThickness", &Utilities::BoundaryParameterObject::mapThickness)
        .def_readwrite("mapResolution", &Utilities::BoundaryParameterObject::mapResolution)
        .def_readwrite("samplingMode", &Utilities::BoundaryParameterObject::samplingMode)
        .def_readwrite("isAnimated", &Utilities::BoundaryParameterObject::isAnimated);

    py::class_<Utilities::FluidModelParameterObject, GenParam::ParameterObject>(m_sub_sub, "FluidData")
        .def(py::init<>())
        .def(py::init<std::string, std::string, std::string, Vector3r, Vector3r, Real, Vector3r, Vector3r, Vector3r, unsigned char,
                bool, std::array<unsigned int, 3>>(),
                "id"_a="Fluid", "visMeshFile"_a = "", "samplesFile"_a, "translation"_a=Vector3r::Zero(), // TODO: samples file here has no default because it needs to be provided
                "axis"_a = Vector3r(1, 0, 0), "angle"_a = 0.0, "scale"_a=Vector3r::Ones(),
                "initialVelocity"_a=Vector3r::Zero(), "initialAngularVelocity"_a = Vector3r::Zero(), "mode"_a=0, "invert"_a=false,
                "resolutionSDF"_a=std::array<unsigned int, 3>({20, 20, 20}))
        .def_readwrite("id", &Utilities::FluidModelParameterObject::id)
        .def_readwrite("visMeshFile", &Utilities::FluidModelParameterObject::visMeshFile)
        .def_readwrite("samplesFile", &Utilities::FluidModelParameterObject::samplesFile)
        .def_readwrite("translation", &Utilities::FluidModelParameterObject::translation)
        .def_readwrite("axis", &Utilities::FluidModelParameterObject::axis)
        .def_readwrite("angle", &Utilities::FluidModelParameterObject::angle)
        .def_readwrite("scale", &Utilities::FluidModelParameterObject::scale)
        .def_readwrite("initialVelocity", &Utilities::FluidModelParameterObject::initialVelocity)
        .def_readwrite("initialAngularVelocity", &Utilities::FluidModelParameterObject::initialAngularVelocity)
        .def_readwrite("mode", &Utilities::FluidModelParameterObject::mode)
        .def_readwrite("invert", &Utilities::FluidModelParameterObject::invert)
        .def_readwrite("resolutionSDF", &Utilities::FluidModelParameterObject::resolutionSDF);

    py::class_<Utilities::FluidBlockParameterObject, GenParam::ParameterObject>(m_sub_sub, "FluidBlock")
        .def(py::init<>())
        .def(py::init<std::string, std::string, Vector3r, Vector3r, unsigned char, Vector3r, Vector3r, Vector3r, Vector3r>(),
                "id"_a="Fluid", "visMeshFile"_a="", "boxMin"_a=Vector3r::Zero(), "boxMax"_a = Vector3r::Ones(),
                "mode"_a=0, "translation"_a = Vector3r::Zero(), "scale"_a = Vector3r::Ones(), "initialVelocity"_a=Vector3r::Zero(), "initialVelocity"_a = Vector3r::Zero())
        .def_readwrite("id", &Utilities::FluidBlockParameterObject::id)
        .def_readwrite("visMeshFile", &Utilities::FluidBlockParameterObject::visMeshFile)
        .def_readwrite("boxMin", &Utilities::FluidBlockParameterObject::boxMin)
        .def_readwrite("boxMax", &Utilities::FluidBlockParameterObject::boxMax)
        .def_readwrite("mode", &Utilities::FluidBlockParameterObject::mode)
        .def_readwrite("translation", &Utilities::FluidBlockParameterObject::translation)
        .def_readwrite("scale", &Utilities::FluidBlockParameterObject::scale)
        .def_readwrite("initialVelocity", &Utilities::FluidBlockParameterObject::initialVelocity)
        .def_readwrite("initialAngularVelocity", &Utilities::FluidBlockParameterObject::initialAngularVelocity);

    py::class_<Utilities::EmitterParameterObject, GenParam::ParameterObject>(m_sub_sub, "EmitterData")
        .def(py::init<>())
        .def(py::init<std::string, unsigned int, unsigned int, Vector3r, Real, Vector3r, Real, Real, Real, unsigned int>(),
                "id"_a="Fluid", "width"_a=5, "height"_a=5, "x"_a=Vector3r::Zero(),
                "velocity"_a=1, "axis"_a = Vector3r(1, 0, 0), "angle"_a = 0.0, "emitStartTime"_a=0,
                "emitEndTime"_a=std::numeric_limits<Real>::max(), "type"_a=0)
        .def_readwrite("id", &Utilities::EmitterParameterObject::id)
        .def_readwrite("width", &Utilities::EmitterParameterObject::width)
        .def_readwrite("height", &Utilities::EmitterParameterObject::height)
        .def_readwrite("x", &Utilities::EmitterParameterObject::x)
        .def_readwrite("velocity", &Utilities::EmitterParameterObject::velocity)
        .def_readwrite("axis", &Utilities::EmitterParameterObject::axis)
        .def_readwrite("angle", &Utilities::EmitterParameterObject::angle)
        .def_readwrite("emitStartTime", &Utilities::EmitterParameterObject::emitStartTime)
        .def_readwrite("emitEndTime", &Utilities::EmitterParameterObject::emitEndTime)
        .def_readwrite("type", &Utilities::EmitterParameterObject::type);

    py::class_<Utilities::AnimationFieldParameterObject, GenParam::ParameterObject>(m_sub_sub, "AnimationFieldData")
        .def(py::init<>())
        .def(py::init<std::string, std::string, std::string, std::string, unsigned int, Vector3r, Vector3r, Real, Vector3r, Real, Real>(),
                "particleFieldName"_a="", "expressionX"_a="", "expressionY"_a="",
                "expressionZ"_a="", "shapeType"_a=0, "translation"_a=Vector3r::Zero(),
                "axis"_a = Vector3r(1, 0, 0), "angle"_a = 0.0, "scale"_a=Vector3r::Ones(), "startTime"_a=0,
                "endTime"_a=std::numeric_limits<Real>::max())
        .def_readwrite("particleFieldName", &Utilities::AnimationFieldParameterObject::particleFieldName)
                // .def_readwrite("expression", &Utilities::AnimationFieldParameterObject::expression) // TODO: bind this type manually
        .def_readwrite("shapeType", &Utilities::AnimationFieldParameterObject::shapeType)
        .def_readwrite("translation", &Utilities::AnimationFieldParameterObject::translation)
        .def_readwrite("axis", &Utilities::AnimationFieldParameterObject::axis)
        .def_readwrite("angle", &Utilities::AnimationFieldParameterObject::angle)
        .def_readwrite("scale", &Utilities::AnimationFieldParameterObject::scale)
        .def_readwrite("startTime", &Utilities::AnimationFieldParameterObject::startTime)
        .def_readwrite("endTime", &Utilities::AnimationFieldParameterObject::endTime);

    py::class_<Utilities::MaterialParameterObject, GenParam::ParameterObject>(m_sub_sub, "MaterialData")
        .def(py::init<>())
        .def(py::init<std::string, std::string, unsigned int, Real, Real, unsigned int, bool, Vector3r, Vector3r>(),
                "id"_a, "colorField"_a="velocity", "colorMapType"_a=1, "minVal"_a=0.0, "maxVal"_a=10.0, //TODO: an id has to be provided
                "maxEmitterParticles"_a=10000, "emitterReuseParticles"_a=false, "emitterBoxMin"_a=Vector3r(-1.0, -1.0, -1.0),
                "emitterBoxMax"_a=Vector3r(1.0, 1.0, 1.0))
        .def_readwrite("id", &Utilities::MaterialParameterObject::id)
        .def_readwrite("colorField", &Utilities::MaterialParameterObject::colorField)
        .def_readwrite("colorMapType", &Utilities::MaterialParameterObject::colorMapType)
        .def_readwrite("minVal", &Utilities::MaterialParameterObject::minVal)
        .def_readwrite("maxVal", &Utilities::MaterialParameterObject::maxVal)
        .def_readwrite("maxEmitterParticles", &Utilities::MaterialParameterObject::maxEmitterParticles)
        .def_readwrite("emitterReuseParticles", &Utilities::MaterialParameterObject::emitterReuseParticles)
        .def_readwrite("emitterBoxMin", &Utilities::MaterialParameterObject::emitterBoxMin)
        .def_readwrite("emitterBoxMax", &Utilities::MaterialParameterObject::emitterBoxMax);

    py::bind_pointer_vector<std::vector<Utilities::BoundaryParameterObject*>>(m_sub, "BoundaryParameterObjectVector");
    py::bind_pointer_vector<std::vector<Utilities::FluidModelParameterObject *>>(m_sub, "FluidModelParameterObjectVector");
    py::bind_pointer_vector<std::vector<Utilities::FluidBlockParameterObject *>>(m_sub, "FluidBlockParameterObjectVector");
    py::bind_pointer_vector<std::vector<Utilities::EmitterParameterObject *>>(m_sub, "EmitterParameterObjectVector");
    py::bind_pointer_vector<std::vector<Utilities::AnimationFieldParameterObject*>>(m_sub, "AnimationFieldParameterObjectVector");
    py::bind_pointer_vector<std::vector<Utilities::MaterialParameterObject *>>(m_sub, "MaterialParameterObjectVector");

    using SceneInfo = Utilities::SceneLoader;
    py::class_<SceneInfo::Scene>(m_sub_sub, "Scene")
        .def(py::init<>())
        .def_readwrite("boundaryModels", &SceneInfo::Scene::boundaryModels)
        .def_readwrite("fluidModels", &SceneInfo::Scene::fluidModels)
        .def_readwrite("fluidBlocks", &SceneInfo::Scene::fluidBlocks)
        .def_readwrite("emitters", &SceneInfo::Scene::emitters)
        .def_readwrite("animatedFields", &SceneInfo::Scene::animatedFields)
        .def_readwrite("materials", &SceneInfo::Scene::materials)
        .def_readwrite("particleRadius", &SceneInfo::Scene::particleRadius)
        .def_readwrite("sim2D", &SceneInfo::Scene::sim2D);

    // ---------------------------------------
    // SDF Functions TODO: implement discregrid
    // ---------------------------------------
    py::class_<Discregrid::CubicLagrangeDiscreteGrid>(m_sub, "DiscreteGrid");

    py::class_<Utilities::SDFFunctions>(m_sub, "SDFFunctions")
            .def_static("generateSDF", &Utilities::SDFFunctions::generateSDF)
            .def_static("computeBoundingBox", &Utilities::SDFFunctions::computeBoundingBox)
            .def_static("distance",
                        overload_cast_<Discregrid::CubicLagrangeDiscreteGrid *, const Vector3r &, const Real, Vector3r &, Vector3r &>()(
                                &Utilities::SDFFunctions::distance))
            .def_static("distance",
                        overload_cast_<Discregrid::CubicLagrangeDiscreteGrid *, const Vector3r &, const Real>()(
                                &Utilities::SDFFunctions::distance));

    // ---------------------------------------
    // Simple Quadrature
    // ---------------------------------------
    py::class_<SPH::SimpleQuadrature>(m_sub, "SimpleQuadrature")
            .def_readwrite_static("m_samplePoints", &SPH::SimpleQuadrature::m_samplePoints)
            .def_readwrite_static("m_volume", &SPH::SimpleQuadrature::m_volume)

            .def_static("determineSamplePointsInSphere", &SPH::SimpleQuadrature::determineSamplePointsInSphere)
            .def_static("determineSamplePointsInCircle", &SPH::SimpleQuadrature::determineSamplePointsInCircle)
            .def_static("integrate", &SPH::SimpleQuadrature::integrate);

    // ---------------------------------------
    // Surface Sampling Modes
    // ---------------------------------------
    py::enum_<SPH::SurfaceSamplingMode>(m_sub, "SurfaceSamplingMode")
            .value("PoissonDisk", SPH::SurfaceSamplingMode::PoissonDisk)
            .value("RegularTriangle", SPH::SurfaceSamplingMode::RegularTriangle)
            .value("Regular2D", SPH::SurfaceSamplingMode::Regular2D);

    // ---------------------------------------
    // Volume Sampling
    // ---------------------------------------
    py::class_<Utilities::VolumeSampling>(m_sub, "VolumeSampling")
            .def_static("sampleMesh", Utilities::VolumeSampling::sampleMesh);

    // ---------------------------------------
    // Winding Numbers
    // ---------------------------------------
    py::class_<Utilities::WindingNumbers>(m_sub, "WindingNumbers")
            .def_static("computeGeneralizedWindingNumber",
                        overload_cast_<const Vector3r &, const Vector3r &, const Vector3r &, const Vector3r &>()(
                                &Utilities::WindingNumbers::computeGeneralizedWindingNumber))
            .def_static("computeGeneralizedWindingNumber",
                        overload_cast_<const Vector3r &, const SPH::TriangleMesh &>()(
                                &Utilities::WindingNumbers::computeGeneralizedWindingNumber));

    // ---------------------------------------
    // Mesh Import
    // ---------------------------------------
    py::class_<SPH::MeshImport>(m_sub, "MeshImport")
        .def_static("importMesh", SPH::MeshImport::importMesh);
}
