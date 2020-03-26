//
// Created by Stefan on 13.02.2020.
//
#include "common.h"

#include <string>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void ExtrasModule(py::module m){
    auto m_sub = m.def_submodule("Extras");
    auto m_sub_sub = m_sub.def_submodule("Scenes");
    m_sub_sub.attr("AnimationFields_2D") = "data/Scenes/AnimationFields_2D.json";
    m_sub_sub.attr("BucklingModel_Bender2017") = "data/Scenes/BucklingModel_Bender2017.json";
    m_sub_sub.attr("BucklingModel_Peer2015") = "data/Scenes/BucklingModel_Peer2015.json";
    m_sub_sub.attr("BucklingModel_Peer2016") = "data/Scenes/BucklingModel_Peer2016.json";
    m_sub_sub.attr("BucklingModel_Takahashi2015") = "data/Scenes/BucklingModel_Takahashi2015.json";
    m_sub_sub.attr("BucklingModel_Weiler2018") = "data/Scenes/BucklingModel_Weiler2018.json";
    m_sub_sub.attr("Bunny_vs_Dragon") = "data/Scenes/Bunny_vs_Dragon.json";
    m_sub_sub.attr("CoilingModel_Bender2017") = "data/Scenes/CoilingModel_Bender2017.json";
    m_sub_sub.attr("CoilingModel_Peer2015") = "data/Scenes/CoilingModel_Peer2015.json";
    m_sub_sub.attr("CoilingModel_Peer2016") = "data/Scenes/CoilingModel_Peer2016.json";
    m_sub_sub.attr("CoilingModel_Takahashi2015") = "data/Scenes/CoilingModel_Takahashi2015.json";
    m_sub_sub.attr("CoilingModel_Weiler2018") = "data/Scenes/CoilingModel_Weiler2018.json";
    m_sub_sub.attr("DamBreakModel") = "data/Scenes/DamBreakModel.json";
    m_sub_sub.attr("DamBreakModelDragons") = "data/Scenes/DamBreakModelDragons.json";
    m_sub_sub.attr("DamBreakModel_2D") = "data/Scenes/DamBreakModel_2D.json";
    m_sub_sub.attr("DamBreakWithObjects") = "data/Scenes/DamBreakWithObjects.json";
    m_sub_sub.attr("DeformableModel") = "data/Scenes/DeformableModel.json";
    m_sub_sub.attr("DoubleDamBreak") = "data/Scenes/DoubleDamBreak.json";
    m_sub_sub.attr("DoubleDamBreakMultiPhase") = "data/Scenes/DoubleDamBreakMultiPhase.json";
    m_sub_sub.attr("DoubleDamBreakWithSphere") = "data/Scenes/DoubleDamBreakWithSphere.json";
    m_sub_sub.attr("DragTest") = "data/Scenes/DragTest.json";
    m_sub_sub.attr("Emitter") = "data/Scenes/Emitter.json";
    m_sub_sub.attr("Emitter_2D") = "data/Scenes/Emitter_2D.json";
    m_sub_sub.attr("Empty") = "data/Scenes/Empty.json";
    m_sub_sub.attr("GridModel_Akinci2012") = "data/Scenes/GridModel_Akinci2012.json";
    m_sub_sub.attr("GridModel_Bender2019") = "data/Scenes/GridModel_Bender2019.json";
    m_sub_sub.attr("GridModel_Koschier2017") = "data/Scenes/GridModel_Koschier2017.json";
    m_sub_sub.attr("MotorScene") = "data/Scenes/MotorScene.json";
    m_sub_sub.attr("MotorScene2") = "data/Scenes/MotorScene2.json";
    m_sub_sub.attr("MultiPhaseColoring") = "data/Scenes/MultiPhaseColoring.json";
    m_sub_sub.attr("Obstacle") = "data/Scenes/Obstacle.json";
    m_sub_sub.attr("Sampling_2D") = "data/Scenes/Sampling_2D.json";
    m_sub_sub.attr("ViscousBunny") = "data/Scenes/ViscousBunny.json";
}