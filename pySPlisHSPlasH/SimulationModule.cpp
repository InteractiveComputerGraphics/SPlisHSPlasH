//
// Created by sjeske on 1/22/20.
//
#include "common.h"

#include <SPlisHSPlasH/Simulation.h>
#include <Simulator/SimulatorBase.h>
#include <Simulator/BoundarySimulator.h>
#include <Simulator/StaticBoundarySimulator.h>
#include <Simulator/GUI/Simulator_GUI_Base.h>
#include <Simulator/SceneConfiguration.h>
#include <Simulator/Exporter/ExporterBase.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

#include <iostream>
#include <vector>
#include <string>

namespace py = pybind11;
using namespace pybind11::literals;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void py_init_simulator(SPH::SimulatorBase& obj,
                       std::string sceneFile = "data/Scenes/DoubleDamBreak.json", // TODO: change to empty default
                       std::string programName = "pySPlisHSPlasH",
                       bool useCache = true, 
                       std::string stateFile = "",
                       std::string outputDir = "",
                       bool initialPause = true,
                       bool useGui = true, 
					   float stopAt = -1.0f, 
					   std::string param="")
                       {
                           std::vector<std::string> argv;
                           argv.push_back(programName);
                           argv.push_back("--scene-file"); argv.push_back(sceneFile);
                           if(!useCache) argv.push_back("--no-cache");
                           argv.push_back("--stopAt");
                           argv.push_back(std::to_string(stopAt));
                           if (strcmp(param.c_str(), "") != 0) {
                               argv.push_back("--param");
                               argv.push_back(param);
                           }
                           if(strcmp(stateFile.c_str(), "") != 0) {
                                   argv.push_back("--state-file");
                                   argv.push_back(stateFile);
                           }
                           if(strcmp(outputDir.c_str(), "") != 0) {
                                   argv.push_back("--output-dir");
                                   argv.push_back(outputDir);
                           }
                           if(!initialPause) argv.push_back("--no-initial-pause");
                           if(!useGui) argv.push_back("--no-gui");
                           obj.init(argv, "pySPlisHSPlasH");
                       };

void SimulationModule(py::module m_sub){
    // ---------------------------------------
    // Enum Simulation Methods
    // ---------------------------------------
    py::enum_<SPH::SimulationMethods>(m_sub, "SimulationMethods")
            .value("WCSPH", SPH::SimulationMethods::WCSPH)
            .value("PCISPH", SPH::SimulationMethods::PCISPH)
            .value("PBF", SPH::SimulationMethods::PBF)
            .value("IISPH", SPH::SimulationMethods::IISPH)
            .value("DFSPH", SPH::SimulationMethods::DFSPH)
            .value("PF", SPH::SimulationMethods::PF)
            .value("NumSimulationMethods", SPH::SimulationMethods::NumSimulationMethods);

    // ---------------------------------------
    // Enum Boundary Handling Methods
    // ---------------------------------------
    py::enum_<SPH::BoundaryHandlingMethods>(m_sub, "BoundaryHandlingMethods")
            .value("Akinci2012", SPH::BoundaryHandlingMethods::Akinci2012)
            .value("Koschier2017", SPH::BoundaryHandlingMethods::Koschier2017)
            .value("Bender2019", SPH::BoundaryHandlingMethods::Bender2019)
            .value("NumSimulationMethods", SPH::BoundaryHandlingMethods::NumSimulationMethods);

    // ---------------------------------------
    // Class Simulation
    // ---------------------------------------
    py::class_<SPH::Simulation, GenParam::ParameterObject>(m_sub, "Simulation")
            .def_readwrite_static("SIM_2D", &SPH::Simulation::SIM_2D)
            .def_readwrite_static("PARTICLE_RADIUS", &SPH::Simulation::PARTICLE_RADIUS)
            .def_readwrite_static("GRAVITATION", &SPH::Simulation::GRAVITATION)
            .def_readwrite_static("CFL_METHOD", &SPH::Simulation::CFL_METHOD)
            .def_readwrite_static("CFL_FACTOR", &SPH::Simulation::CFL_FACTOR)
			.def_readwrite_static("CFL_MIN_TIMESTEPSIZE", &SPH::Simulation::CFL_MIN_TIMESTEPSIZE)
            .def_readwrite_static("CFL_MAX_TIMESTEPSIZE", &SPH::Simulation::CFL_MAX_TIMESTEPSIZE)
            .def_readwrite_static("ENABLE_Z_SORT", &SPH::Simulation::ENABLE_Z_SORT)

            .def_readwrite_static("KERNEL_METHOD", &SPH::Simulation::KERNEL_METHOD)
            .def_readwrite_static("GRAD_KERNEL_METHOD", &SPH::Simulation::GRAD_KERNEL_METHOD)
            .def_readwrite_static("ENUM_KERNEL_CUBIC", &SPH::Simulation::ENUM_KERNEL_CUBIC)
            .def_readwrite_static("ENUM_KERNEL_WENDLANDQUINTICC2", &SPH::Simulation::ENUM_KERNEL_WENDLANDQUINTICC2)
            .def_readwrite_static("ENUM_KERNEL_POLY6", &SPH::Simulation::ENUM_KERNEL_POLY6)
            .def_readwrite_static("ENUM_KERNEL_SPIKY", &SPH::Simulation::ENUM_KERNEL_SPIKY)
            .def_readwrite_static("ENUM_KERNEL_PRECOMPUTED_CUBIC", &SPH::Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC)
            .def_readwrite_static("ENUM_KERNEL_CUBIC_2D", &SPH::Simulation::ENUM_KERNEL_CUBIC_2D)
            .def_readwrite_static("ENUM_KERNEL_WENDLANDQUINTICC2_2D", &SPH::Simulation::ENUM_KERNEL_WENDLANDQUINTICC2_2D)
            .def_readwrite_static("ENUM_GRADKERNEL_CUBIC", &SPH::Simulation::ENUM_GRADKERNEL_CUBIC)
            .def_readwrite_static("ENUM_GRADKERNEL_WENDLANDQUINTICC2", &SPH::Simulation::ENUM_GRADKERNEL_WENDLANDQUINTICC2)
            .def_readwrite_static("ENUM_GRADKERNEL_POLY6", &SPH::Simulation::ENUM_GRADKERNEL_POLY6)
            .def_readwrite_static("ENUM_GRADKERNEL_SPIKY", &SPH::Simulation::ENUM_GRADKERNEL_SPIKY)
            .def_readwrite_static("ENUM_GRADKERNEL_PRECOMPUTED_CUBIC", &SPH::Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC)
            .def_readwrite_static("ENUM_GRADKERNEL_CUBIC_2D", &SPH::Simulation::ENUM_GRADKERNEL_CUBIC_2D)
            .def_readwrite_static("ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D", &SPH::Simulation::ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D)

            .def_readwrite_static("SIMULATION_METHOD", &SPH::Simulation::SIMULATION_METHOD)

            .def_readwrite_static("ENUM_CFL_NONE", &SPH::Simulation::ENUM_CFL_NONE)
            .def_readwrite_static("ENUM_CFL_STANDARD", &SPH::Simulation::ENUM_CFL_STANDARD)
            .def_readwrite_static("ENUM_CFL_ITER", &SPH::Simulation::ENUM_CFL_ITER)

            .def_readwrite_static("ENUM_SIMULATION_WCSPH", &SPH::Simulation::ENUM_SIMULATION_WCSPH)
            .def_readwrite_static("ENUM_SIMULATION_PCISPH", &SPH::Simulation::ENUM_SIMULATION_PCISPH)
            .def_readwrite_static("ENUM_SIMULATION_PBF", &SPH::Simulation::ENUM_SIMULATION_PBF)
            .def_readwrite_static("ENUM_SIMULATION_IISPH", &SPH::Simulation::ENUM_SIMULATION_IISPH)
            .def_readwrite_static("ENUM_SIMULATION_DFSPH", &SPH::Simulation::ENUM_SIMULATION_DFSPH)
            .def_readwrite_static("ENUM_SIMULATION_PF", &SPH::Simulation::ENUM_SIMULATION_PF)

            .def_readwrite_static("BOUNDARY_HANDLING_METHOD", &SPH::Simulation::BOUNDARY_HANDLING_METHOD)
            .def_readwrite_static("ENUM_AKINCI2012", &SPH::Simulation::ENUM_AKINCI2012)
            .def_readwrite_static("ENUM_KOSCHIER2017", &SPH::Simulation::ENUM_KOSCHIER2017)
            .def_readwrite_static("ENUM_BENDER2019", &SPH::Simulation::ENUM_BENDER2019)

            .def(py::init<>())
            .def("init", &SPH::Simulation::init)
            .def("reset", &SPH::Simulation::reset)
            .def_static("getCurrent", &SPH::Simulation::getCurrent, py::return_value_policy::reference)
            .def_static("setCurrent", &SPH::Simulation::setCurrent)
            .def_static("hasCurrent", &SPH::Simulation::hasCurrent)

            .def("addFluidModel", [](SPH::Simulation& current, const std::string &id, std::vector<Vector3r> fluidParticles, std::vector<Vector3r> fluidVelocities, const unsigned int nMaxEmitterParticles){
                if(fluidParticles.size() != fluidVelocities.size())
                    throw std::runtime_error("Sizes of position and velocity array must be equal");
                current.addFluidModel(id, fluidParticles.size(), fluidParticles.data(), fluidVelocities.data(), nMaxEmitterParticles);
            })
            .def("getFluidModel", &SPH::Simulation::getFluidModel, py::return_value_policy::reference_internal)
            .def("getFluidModelFromPointSet", &SPH::Simulation::getFluidModelFromPointSet, py::return_value_policy::reference_internal)
            .def("numberOfFluidModels", &SPH::Simulation::numberOfFluidModels)

            .def("addBoundaryModel", &SPH::Simulation::addBoundaryModel)
            .def("getBoundaryModel", &SPH::Simulation::getBoundaryModel, py::return_value_policy::reference_internal)
            .def("getBoundaryModelFromPointSet", &SPH::Simulation::getBoundaryModelFromPointSet, py::return_value_policy::reference_internal)
            .def("numberOfBoundaryModels", &SPH::Simulation::numberOfBoundaryModels)
            .def("updateBoundaryVolume", &SPH::Simulation::updateBoundaryVolume)

            .def("getAnimationFieldSystem", &SPH::Simulation::getAnimationFieldSystem, py::return_value_policy::reference_internal)

            .def("getBoundaryHandlingMethod", &SPH::Simulation::getBoundaryHandlingMethod)
            .def("setBoundaryHandlingMethod", &SPH::Simulation::setBoundaryHandlingMethod)

            .def("getKernel", &SPH::Simulation::getKernel)
            .def("setKernel", &SPH::Simulation::setKernel)
            .def("getGradKernel", &SPH::Simulation::getGradKernel)
            .def("setGradKernel", &SPH::Simulation::setGradKernel)

            .def("W_zero", &SPH::Simulation::W_zero)
            .def("W", &SPH::Simulation::W)
            .def("gradW", &SPH::Simulation::gradW)

            .def("getSimulationMethod", &SPH::Simulation::getSimulationMethod)
            .def("setSimulationMethod", &SPH::Simulation::setSimulationMethod)
            .def("setSimulationMethodChangedCallback", &SPH::Simulation::setSimulationMethodChangedCallback)
            .def("getTimeStep", &SPH::Simulation::getTimeStep, py::return_value_policy::reference_internal) // TODO: This returns abstract class pointer, figure out what to do with it
            .def("is2DSimulation", &SPH::Simulation::is2DSimulation)
            .def("zSortEnabled", &SPH::Simulation::zSortEnabled)
            .def("setParticleRadius", &SPH::Simulation::setParticleRadius)
            .def("getParticleRadius", &SPH::Simulation::getParticleRadius)
            .def("getSupportRadius", &SPH::Simulation::getSupportRadius)
            .def("updateTimeStepSize", &SPH::Simulation::updateTimeStepSize)
            .def("updateTimeStepSizeCFL", &SPH::Simulation::updateTimeStepSizeCFL)
            .def("performNeighborhoodSearch", &SPH::Simulation::performNeighborhoodSearch)
            .def("performNeighborhoodSearchSort", &SPH::Simulation::performNeighborhoodSearchSort)
            .def("computeNonPressureForces", &SPH::Simulation::computeNonPressureForces)
            .def("animateParticles", &SPH::Simulation::animateParticles)
            .def("emitParticles", &SPH::Simulation::emitParticles)
            .def("emittedParticles", &SPH::Simulation::emittedParticles)
            .def("getNeighborhoodSearch", &SPH::Simulation::getNeighborhoodSearch, py::return_value_policy::reference_internal)
            .def("saveState", &SPH::Simulation::saveState)
            .def("loadState", &SPH::Simulation::loadState)
            .def("addViscosityMethod", &SPH::Simulation::addViscosityMethod)
            .def("getViscosityMethods", &SPH::Simulation::getViscosityMethods)
            .def("addVorticityMethod", &SPH::Simulation::addVorticityMethod)
            .def("getVorticityMethods", &SPH::Simulation::getVorticityMethods)
            .def("addDragMethod", &SPH::Simulation::addDragMethod)
            .def("getDragMethods", &SPH::Simulation::getDragMethods)
            .def("numberOfPointSets", &SPH::Simulation::numberOfPointSets)
            .def("numberOfNeighbors", &SPH::Simulation::numberOfNeighbors)
            .def("getNeighbor", &SPH::Simulation::getNeighbor)
            .def("getNeighborList", &SPH::Simulation::getNeighborList); // TODO: Might not work because of array pointer



    // ---------------------------------------
    // EXEC SUBMODULE
    // ---------------------------------------
    m_sub = m_sub.def_submodule("Exec");

    // ---------------------------------------
    // Struct Exporter
    // ---------------------------------------
    py::class_<SPH::SimulatorBase::Exporter>(m_sub, "Exporter")
            .def_readwrite("key", &SPH::SimulatorBase::Exporter::m_key)
            .def_readwrite("name", &SPH::SimulatorBase::Exporter::m_name)
            .def_readwrite("decription", &SPH::SimulatorBase::Exporter::m_description)
            .def_readwrite("id", &SPH::SimulatorBase::Exporter::m_id)
            .def_readwrite("exporter", &SPH::SimulatorBase::Exporter::m_exporter);

    // ---------------------------------------
    // Simulator Base class
    // ---------------------------------------
    py::class_<SPH::SimulatorBase::SimulationMethod>(m_sub, "SimulationMethod")
            .def_readwrite("simulationMethod", &SPH::SimulatorBase::SimulationMethod::simulationMethod)
            .def_readwrite("simulation", &SPH::SimulatorBase::SimulationMethod::simulation)
            .def_readonly("model", &SPH::SimulatorBase::SimulationMethod::model); // TODO: this is public property but defined as readonly because of deleted assignment operator. check in future

    py::class_<SPH::SimulatorBase, GenParam::ParameterObject>(m_sub, "SimulatorBase")
            .def_readwrite_static("PAUSE", &SPH::SimulatorBase::PAUSE)
            .def_readwrite_static("PAUSE_AT", &SPH::SimulatorBase::PAUSE_AT)
            .def_readwrite_static("STOP_AT", &SPH::SimulatorBase::STOP_AT)
            .def_readwrite_static("NUM_STEPS_PER_RENDER", &SPH::SimulatorBase::NUM_STEPS_PER_RENDER)
            .def_readwrite_static("DATA_EXPORT_FPS", &SPH::SimulatorBase::DATA_EXPORT_FPS)
            .def_readwrite_static("PARTICLE_EXPORT_ATTRIBUTES", &SPH::SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES)
            .def_readwrite_static("STATE_EXPORT", &SPH::SimulatorBase::STATE_EXPORT)
            .def_readwrite_static("STATE_EXPORT_FPS", &SPH::SimulatorBase::STATE_EXPORT_FPS)
            .def_readwrite_static("RENDER_WALLS", &SPH::SimulatorBase::RENDER_WALLS)

            .def_readwrite_static("ENUM_WALLS_NONE", &SPH::SimulatorBase::ENUM_WALLS_NONE)
            .def_readwrite_static("ENUM_WALLS_PARTICLES_ALL", &SPH::SimulatorBase::ENUM_WALLS_PARTICLES_ALL)
            .def_readwrite_static("ENUM_WALLS_PARTICLES_NO_WALLS", &SPH::SimulatorBase::ENUM_WALLS_PARTICLES_NO_WALLS)
            .def_readwrite_static("ENUM_WALLS_GEOMETRY_ALL", &SPH::SimulatorBase::ENUM_WALLS_GEOMETRY_ALL)
            .def_readwrite_static("ENUM_WALLS_GEOMETRY_NO_WALLS", &SPH::SimulatorBase::ENUM_WALLS_GEOMETRY_NO_WALLS)

            .def(py::init<>())
            .def("run", &SPH::SimulatorBase::run)
            // .def("init", [](SPH::SimulatorBase& obj, std::vector<std::string> argv, std::string windowName){
            //     std::vector<const char *> cargv;
            //     cargv.reserve(argv.size());
            //     for (auto & elem : argv) {
            //         cargv.push_back(elem.c_str());
            //     }
            //     obj.init(argv.size(), const_cast<char**>(cargv.data()), windowName);
            // })
            .def("init", overload_cast_<std::vector<std::string>, const std::string&>()(&SPH::SimulatorBase::init))
            .def("init", &py_init_simulator, 
                            "sceneFile"_a = "data/Scenes/DoubleDamBreak.json",
                            "programName"_a = "pySPlisHSPlasH",
                            "useCache"_a = true,
                            "stateFile"_a = "",
                            "outputDir"_a = "",
                            "initialPause"_a = true,
                            "useGui"_a = true,
							"stopAt"_a = -1.0,
							"param"_a = "")
            .def("initSimulation", &SPH::SimulatorBase::initSimulation)
			.def("runSimulation", &SPH::SimulatorBase::runSimulation)
            .def("cleanup", &SPH::SimulatorBase::cleanup)

            .def("reset", &SPH::SimulatorBase::reset)
            .def("timeStep", &SPH::SimulatorBase::timeStep)
            .def("timeStepNoGUI", &SPH::SimulatorBase::timeStepNoGUI)

            .def_static("particleInfo", &SPH::SimulatorBase::particleInfo)
            .def("real2String", &SPH::SimulatorBase::real2String)
            .def("initDensityMap", &SPH::SimulatorBase::initDensityMap)
            .def("initVolumeMap", &SPH::SimulatorBase::initVolumeMap)

            .def("readParameters", &SPH::SimulatorBase::readParameters)
            .def("step", &SPH::SimulatorBase::step)

            .def("saveState", &SPH::SimulatorBase::saveState)
            // .def("loadStateDialog", &SPH::SimulatorBase::loadStateDialog)
            .def("loadState", &SPH::SimulatorBase::loadState)
            .def("writeFluidParticlesState", &SPH::SimulatorBase::writeFluidParticlesState)
            .def("readFluidParticlesState", &SPH::SimulatorBase::readFluidParticlesState)
            .def("writeBoundaryState", &SPH::SimulatorBase::writeBoundaryState)
            .def("readBoundaryState", &SPH::SimulatorBase::readBoundaryState)
            .def("writeParameterState", &SPH::SimulatorBase::writeParameterState)
            .def("readParameterState", &SPH::SimulatorBase::readParameterState)
            .def("writeParameterObjectState", &SPH::SimulatorBase::writeParameterObjectState)
            .def("readParameterObjectState", &SPH::SimulatorBase::readParameterObjectState)

            .def("updateBoundaryParticles", &SPH::SimulatorBase::updateBoundaryParticles)
            .def("updateDMVelocity", &SPH::SimulatorBase::updateDMVelocity)
            .def("updateVMVelocity", &SPH::SimulatorBase::updateVMVelocity)

            .def("getScalarField", &SPH::SimulatorBase::getScalarField)
            .def("updateScalarField", &SPH::SimulatorBase::updateScalarField)
            .def("determineMinMaxOfScalarField", &SPH::SimulatorBase::determineMinMaxOfScalarField)

            .def_static("loadObj", &SPH::SimulatorBase::loadObj)

            .def("getSceneLoader", &SPH::SimulatorBase::getSceneLoader, py::return_value_policy::reference_internal)

            .def("getExePath", &SPH::SimulatorBase::getExePath)
            
            .def("getUseParticleCaching", &SPH::SimulatorBase::getUseParticleCaching)
            .def("setUseParticleCaching", &SPH::SimulatorBase::setUseParticleCaching)
            .def("getUseGUI", &SPH::SimulatorBase::getUseGUI)
            .def("setUseGUI", &SPH::SimulatorBase::setUseGUI)

            .def("getColorField", &SPH::SimulatorBase::getColorField)
            .def("setColorField", &SPH::SimulatorBase::setColorField)

            .def("getColorMapType", &SPH::SimulatorBase::getColorMapType)
            .def("setColorMapType", &SPH::SimulatorBase::setColorMapType)
            .def("getRenderMaxValue", &SPH::SimulatorBase::getRenderMaxValue)
            .def("setRenderMaxValue", &SPH::SimulatorBase::setRenderMaxValue)
            .def("getRenderMinValue", &SPH::SimulatorBase::getRenderMinValue)
            .def("setRenderMinValue", &SPH::SimulatorBase::setRenderMinValue)
            .def("getOutputPath", &SPH::SimulatorBase::getOutputPath)

            .def("getStateFile", &SPH::SimulatorBase::getStateFile)
            .def("setStateFile", &SPH::SimulatorBase::setStateFile)

            .def("getBoundarySimulator", &SPH::SimulatorBase::getBoundarySimulator, py::return_value_policy::reference_internal)
            .def("setBoundarySimulator", &SPH::SimulatorBase::setBoundarySimulator)
            .def("getGui", &SPH::SimulatorBase::getGui, py::return_value_policy::reference_internal)
            .def("setGui", &SPH::SimulatorBase::setGui)
            .def("isStaticScene", &SPH::SimulatorBase::isStaticScene)

            .def("addParticleExporter", &SPH::SimulatorBase::addParticleExporter)
            .def("getParticleExporters", &SPH::SimulatorBase::getParticleExporters)
            .def("addRigidBodyExporter", &SPH::SimulatorBase::addRigidBodyExporter)
            .def("getRigidBodyExporters", &SPH::SimulatorBase::getRigidBodyExporters)

            .def("activateExporter", &SPH::SimulatorBase::activateExporter)

			.def("setTimeStepCB", &SPH::SimulatorBase::setTimeStepCB);

     // ---------------------------------------
    // SceneConfiguration
    // ---------------------------------------
    py::class_<SPH::SceneConfiguration>(m_sub, "SceneConfiguration")
            .def_static("getCurrent", &SPH::SceneConfiguration::getCurrent, py::return_value_policy::reference)
            .def_static("setCurrent", &SPH::SceneConfiguration::setCurrent)
            .def_static("hasCurrent", &SPH::SceneConfiguration::hasCurrent)
            .def("getSceneFile", &SPH::SceneConfiguration::getSceneFile)
            .def("getScene", &SPH::SceneConfiguration::getScene, py::return_value_policy::reference_internal);

    // ---------------------------------------
    // BoundarySimulator
    // ---------------------------------------
    py::class_<SPH::BoundarySimulator>(m_sub, "BoundarySimulator")
            .def(py::init<>())
            .def("init", &SPH::BoundarySimulator::init)
            .def("timeStep", &SPH::BoundarySimulator::timeStep)
            .def("initBoundaryData", &SPH::BoundarySimulator::initBoundaryData)
            .def("reset", &SPH::BoundarySimulator::reset)
            .def("updateBoundaryForces", &SPH::BoundarySimulator::updateBoundaryForces);

    // ---------------------------------------
    // StaticBoundarySimulator
    // ---------------------------------------
    py::class_<SPH::StaticBoundarySimulator, SPH::BoundarySimulator>(m_sub, "StaticBoundarySimulator")
            .def(py::init<SPH::SimulatorBase*>());

}

