#include <pybind11/pybind11.h>

// Put this here for now
#include "common.h"

#include <SPlisHSPlasH/Common.h>
#include <Utilities/Logger.h>
#include <Utilities/Timing.h>
#include <Utilities/Counting.h>
// INIT_LOGGING
// INIT_TIMING
// INIT_COUNTING
// TODO: Move this to simulator base class

namespace py = pybind11;

void AnimationFieldModule(py::module);
void ParameterObjectModule(py::module m);
void BoundaryModelModule(py::module);
void EmitterModule(py::module);
void FluidModelModule(py::module);
void GUIModule(py::module);
void ExporterModule(py::module);
void NonPressureForceBaseModule(py::module);
void UtilitiesModule(py::module);
void SimulationModule(py::module);
void SPHKernelsModule(py::module);
void RigidBodyModule(py::module m_sub);
void TimeModule(py::module);
void TriangleMeshModule(py::module);
void PBDModule(py::module);
void DFSPHModule(py::module);
void DragModule(py::module);
void ElasticityModule(py::module);
void IISPHModule(py::module);
void PBFModule(py::module);
void PCISPHModule(py::module);
void PFModule(py::module);
void SurfaceTensionModule(py::module);
void ViscosityModule(py::module);
void VorticityModule(py::module);
void WCSPHModule(py::module);

void ExtrasModule(py::module);


PYBIND11_MODULE(pysplishsplash, m) {
    ParameterObjectModule(m);
    SPHKernelsModule(m);
    AnimationFieldModule(m);
    UtilitiesModule(m);
    BoundaryModelModule(m);
    EmitterModule(m);
    FluidModelModule(m);
    NonPressureForceBaseModule(m);
    SimulationModule(m);
    RigidBodyModule(m);
    TimeModule(m);
    TriangleMeshModule(m);
    PBDModule(m);
    DFSPHModule(m);
    DragModule(m);
    ElasticityModule(m);
    IISPHModule(m);
    PBFModule(m);
    PCISPHModule(m);
    PFModule(m);
    SurfaceTensionModule(m);
    ViscosityModule(m);
    VorticityModule(m);
    WCSPHModule(m);
    GUIModule(m);
    ExporterModule(m);

    ExtrasModule(m);
}