#include "DynamicParameterSystem.h"
#include "FluidModel.h"
#include "Simulation.h"
#include "TimeManager.h"
#include "Utilities/Logger.h"
#include "WCSPH/TimeStepWCSPH.h"
#include "PF/TimeStepPF.h"
#include "ICSPH/TimeStepICSPH.h"
#include "Viscosity/Viscosity_Standard.h"
#include "Viscosity/Viscosity_Bender2017.h"
#include "Viscosity/Viscosity_Peer2015.h"
#include "Viscosity/Viscosity_Peer2016.h"
#include "Viscosity/Viscosity_Takahashi2015.h"
#include "Viscosity/Viscosity_Weiler2018.h"
#include "SurfaceTension/SurfaceTension_Becker2007.h"
#include "SurfaceTension/SurfaceTension_Akinci2013.h"
#include "SurfaceTension/SurfaceTension_He2014.h"
#include "SurfaceTension/SurfaceTension_Jeske2023.h"
#include "SurfaceTension/SurfaceTension_ZorillaRitter2020.h"
#include <algorithm>

using namespace SPH;

DynamicParameterSystem::DynamicParameterSystem() :
    m_t_double(0.0),
    m_dt_double(0.0)
{
    LOG_INFO << "Dynamic Parameters created"; //Debug
}

DynamicParameterSystem::~DynamicParameterSystem()
{
    reset();
}

void DynamicParameterSystem::step()
{
    const Real currentTime = TimeManager::getCurrent()->getTime();
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();

    m_t_double = currentTime;
    m_dt_double = dt;

    for (auto& schedule : m_schedules)
    {
        if (!schedule.active)
            continue;

        Real new_value;

        if (schedule.use_expression && schedule.compiled_expr.get() != nullptr)
        {
            // Evaluate
            new_value = static_cast<Real>(te_eval(schedule.compiled_expr.get()));
            if (std::isnan(new_value) || std::isinf(new_value)) {
                LOG_ERR << "Expression evaluation failed for " << schedule.parameter_name;
                schedule.active = false;
                continue;
            }
        }
        else
        {
            if (schedule.use_step_function) 
            {
                new_value = stepTimeline(schedule, currentTime);
                if (std::isnan(new_value))
                    continue;
            }  
            else
            {
                new_value = interpolateTimeline(schedule, currentTime);
                if (std::isnan(new_value))
                    continue;
            }
            if (currentTime >= schedule.timeline.back().first)
                schedule.active = false;
        }

        schedule.current_value = new_value;
        applyParameterValue(schedule.fluidId, schedule.parameter_name, new_value);
    }
}

void DynamicParameterSystem::reset()
{
    for (auto& schedule : m_schedules)
    {
        applyParameterValue(schedule.fluidId, schedule.parameter_name, schedule.default_value);
    }
    m_schedules.clear();
}

bool DynamicParameterSystem::addTimelineSchedule(const std::string& fluidId, const std::string& paramName, const std::vector<std::pair<Real, Real>>& timeline, Real defaultValue, bool stepFunction)
{
    ParameterSchedule schedule;
    schedule.fluidId = fluidId;
    schedule.parameter_name = paramName;
    schedule.timeline = timeline;
    schedule.default_value = defaultValue;
    schedule.current_value = defaultValue;
    schedule.use_expression = false;
    schedule.use_step_function = stepFunction;
    schedule.active = true;

    // Order timeline
    std::sort(schedule.timeline.begin(), schedule.timeline.end());

    m_schedules.push_back(std::move(schedule));

    LOG_INFO << "Added timeline schedule for " << paramName << " in Fluid " << fluidId
        << " with " << timeline.size() << " points"
        << (stepFunction ? " (step function)" : " (interpolated)");
    return true;
}

bool DynamicParameterSystem::addExpressionSchedule(const std::string& fluidId, const std::string& paramName, const std::string& expression, Real defaultValue)
{
    ParameterSchedule schedule;
    schedule.fluidId = fluidId;
    schedule.parameter_name = paramName;
    schedule.expression = expression;
    schedule.default_value = defaultValue;
    schedule.current_value = defaultValue;
    schedule.use_expression = true;
    schedule.active = true;

    if (!compileExpression(schedule))
    {
        LOG_ERR << "Failed to compile expression for " << paramName
            << ": " << expression;
        return false;
    }

    m_schedules.push_back(std::move(schedule)); // why the move here?

    LOG_INFO << "Added expression schedule for " << paramName << " in Fluid " << fluidId << ": " << expression;
    return true;
}

bool DynamicParameterSystem::compileExpression(ParameterSchedule& schedule)
{
    te_variable vars[] = {
        {"t", &m_t_double},
        {"dt", &m_dt_double}
    };
    const int numVars = 2;

    int err;
    te_expr* raw_expr = te_compile(schedule.expression.c_str(), vars, numVars, &err);

    if (raw_expr != nullptr && err == 0) {
        schedule.compiled_expr.reset(raw_expr);
    }

    return (schedule.compiled_expr != nullptr && err == 0);
}

Real DynamicParameterSystem::interpolateTimeline(ParameterSchedule& schedule, Real currentTime)
{
    const auto& timeline = schedule.timeline;

    if (timeline.empty())
        return schedule.default_value;

    // Before first
    if (currentTime < timeline.front().first)
    {
        schedule.current_timeline_index = -1;
        return std::numeric_limits<Real>::quiet_NaN();
    }

    // After last
    if (currentTime >= timeline.back().first)
    {
        return timeline.back().second;
    }

    int old_index = schedule.current_timeline_index;

    while (schedule.current_timeline_index < static_cast<int>(timeline.size()) - 1 && currentTime > timeline[schedule.current_timeline_index + 1].first)
    {
        schedule.current_timeline_index++;
    }

    // Logging
    if (schedule.current_timeline_index != old_index)
    {
        int i = schedule.current_timeline_index;
        LOG_INFO << schedule.parameter_name << " (Fluid: " << schedule.fluidId
            << ") will interpolate from " << timeline[i].second
            << " to " << timeline[i + 1].second
            << " between t=" << timeline[i].first
            << " and t=" << timeline[i + 1].first;
    }

    // Interpolate
    int i = schedule.current_timeline_index;
    const Real t1 = timeline[i].first;
    const Real v1 = timeline[i].second;
    const Real t2 = timeline[i + 1].first;
    const Real v2 = timeline[i + 1].second;

    const Real alpha = (currentTime - t1) / (t2 - t1);
    return v1 + alpha * (v2 - v1);
}

Real DynamicParameterSystem::stepTimeline(ParameterSchedule& schedule, Real currentTime)
{
    const auto& timeline = schedule.timeline;
    if (timeline.empty())
        return schedule.default_value;

    // Before first
    if (currentTime < timeline.front().first)
    {
        schedule.current_timeline_index = -1;
        return std::numeric_limits<Real>::quiet_NaN();
    }

    int old_index = schedule.current_timeline_index;

    while (schedule.current_timeline_index < static_cast<int>(timeline.size()) - 1 && currentTime >= timeline[schedule.current_timeline_index + 1].first)
    {
        schedule.current_timeline_index++;
    }

    // Logging
    if (schedule.current_timeline_index != old_index)
    {
        int i = schedule.current_timeline_index;
        LOG_INFO << schedule.parameter_name << " (Fluid: " << schedule.fluidId
            << ") changed to " << timeline[i].second
            << " at t=" << timeline[i].first;
    }

    return timeline[schedule.current_timeline_index].second;
}

void DynamicParameterSystem::applyParameterValue(const std::string& fluidId, const std::string& paramName, Real value)
{
    const Real currentTime = TimeManager::getCurrent()->getTime();
    Simulation* sim = Simulation::getCurrent();

    if (paramName == "stiffness")
    {
        if (sim->getSimulationMethod() == 0) // WCSPH
        {
            TimeStepWCSPH* timeStep = static_cast<TimeStepWCSPH*>(sim->getTimeStep());
            timeStep->setValue(TimeStepWCSPH::STIFFNESS, value);
            return;
        } 
        else if (sim->getSimulationMethod() == 5) // PF
        {
            TimeStepPF* timeStep = static_cast<TimeStepPF*>(sim->getTimeStep());
            timeStep->setValue(TimeStepPF::STIFFNESS, value);
            return;
        } else {
            LOG_WARN << "Not a valid simulation method.";
        }
    }

    if (paramName == "lambda")
    {
        if (sim->getSimulationMethod() == 6) // ICSPH
        {
            TimeStepICSPH* timeStep = static_cast<TimeStepICSPH*>(sim->getTimeStep());
            timeStep->setValue(TimeStepICSPH::LAMBDA, value);
            return;
        }
        else {
            LOG_WARN << "Not a valid simulation method.";
        }
    }

    for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
    {
        if (sim->getFluidModel(i)->getId() == fluidId)
        {
            FluidModel* model = sim->getFluidModel(i);
            if (paramName == "density0")
            {
                model->setDensity0(value);
            }
            else if (paramName == "viscosity")
            {
                if (model->getViscosityMethod() == 1) // Standard
                {
                    Viscosity_Standard* viscMethod = static_cast<Viscosity_Standard*>(model->getViscosityBase());
                    viscMethod->setValue(Viscosity_Standard::VISCOSITY_COEFFICIENT, value);
                }
                else if (model->getViscosityMethod() == 2) // Bender 2017
                {
                    Viscosity_Bender2017* viscMethod = static_cast<Viscosity_Bender2017*>(model->getViscosityBase());
                    viscMethod->setValue(Viscosity_Bender2017::VISCOSITY_COEFFICIENT, value);
                }
                else if (model->getViscosityMethod() == 3) // Peer 2015
                {
                    Viscosity_Peer2015* viscMethod = static_cast<Viscosity_Peer2015*>(model->getViscosityBase());
                    viscMethod->setValue(Viscosity_Peer2015::VISCOSITY_COEFFICIENT, value);
                }
                else if (model->getViscosityMethod() == 4) // Peer 2016
                {
                    Viscosity_Peer2016* viscMethod = static_cast<Viscosity_Peer2016*>(model->getViscosityBase());
                    viscMethod->setValue(Viscosity_Peer2016::VISCOSITY_COEFFICIENT, value);
                }
                else if (model->getViscosityMethod() == 5) // Takahashi 2015
                {
                    Viscosity_Takahashi2015* viscMethod = static_cast<Viscosity_Takahashi2015*>(model->getViscosityBase());
                    viscMethod->setValue(Viscosity_Takahashi2015::VISCOSITY_COEFFICIENT, value);
                }
                else if (model->getViscosityMethod() == 6) // Weiler 2018
                {
                    Viscosity_Weiler2018* viscMethod = static_cast<Viscosity_Weiler2018*>(model->getViscosityBase());
                    viscMethod->setValue(Viscosity_Weiler2018::VISCOSITY_COEFFICIENT, value);
                }
            }
            else if (paramName == "viscosityBoundary")
            {
                if (model->getViscosityMethod() == 1) // Standard
                {
                    Viscosity_Standard* viscMethod = static_cast<Viscosity_Standard*>(model->getViscosityBase());
                    viscMethod->setValue(Viscosity_Standard::VISCOSITY_COEFFICIENT_BOUNDARY, value);
                }
                else if (model->getViscosityMethod() == 6) // Weiler 2018
                {
                    Viscosity_Weiler2018* viscMethod = static_cast<Viscosity_Weiler2018*>(model->getViscosityBase());
                    viscMethod->setValue(Viscosity_Weiler2018::VISCOSITY_COEFFICIENT_BOUNDARY, value);
                }
            }
            else if (paramName == "surfaceTension")
            {
                if (model->getSurfaceTensionMethod() == 1) // Becker & Teschner 2007
                {
                    SurfaceTension_Becker2007* stMethod = static_cast<SurfaceTension_Becker2007*>(model->getSurfaceTensionBase());
                    stMethod->setValue(SurfaceTension_Becker2007::SURFACE_TENSION, value);
                }
                else if (model->getSurfaceTensionMethod() == 2) // Akinci 2013
                {
                    SurfaceTension_Akinci2013* stMethod = static_cast<SurfaceTension_Akinci2013*>(model->getSurfaceTensionBase());
                    stMethod->setValue(SurfaceTension_Akinci2013::SURFACE_TENSION, value);
                }
                else if (model->getSurfaceTensionMethod() == 3) // He 2014
                {
                    SurfaceTension_He2014* stMethod = static_cast<SurfaceTension_He2014*>(model->getSurfaceTensionBase());
                    stMethod->setValue(SurfaceTension_He2014::SURFACE_TENSION, value);
                }
                else if (model->getSurfaceTensionMethod() == 4) // Jeske 2023
                {
                    SurfaceTension_Jeske2023* stMethod = static_cast<SurfaceTension_Jeske2023*>(model->getSurfaceTensionBase());
                    stMethod->setValue(SurfaceTension_Jeske2023::SURFACE_TENSION, value);
                }
                else {
                    LOG_WARN << "Not a valid surface tension method.";
                }
            }
            else if (paramName == "surfaceTensionBoundary")
            {
                if (model->getSurfaceTensionMethod() == 1) // Becker & Teschner 2007
                {
                    SurfaceTension_Becker2007* stMethod = static_cast<SurfaceTension_Becker2007*>(model->getSurfaceTensionBase());
                    stMethod->setValue(SurfaceTension_Becker2007::SURFACE_TENSION_BOUNDARY, value);
                }
                else if (model->getSurfaceTensionMethod() == 2) // Akinci 2013
                {
                    SurfaceTension_Akinci2013* stMethod = static_cast<SurfaceTension_Akinci2013*>(model->getSurfaceTensionBase());
                    stMethod->setValue(SurfaceTension_Akinci2013::SURFACE_TENSION_BOUNDARY, value);
                }
                else if (model->getSurfaceTensionMethod() == 3) // He 2014
                {
                    SurfaceTension_He2014* stMethod = static_cast<SurfaceTension_He2014*>(model->getSurfaceTensionBase());
                    stMethod->setValue(SurfaceTension_He2014::SURFACE_TENSION_BOUNDARY, value);
                }
                else if (model->getSurfaceTensionMethod() == 4) // Jeske 2023
                {
                    SurfaceTension_Jeske2023* stMethod = static_cast<SurfaceTension_Jeske2023*>(model->getSurfaceTensionBase());
                    stMethod->setValue(SurfaceTension_Jeske2023::SURFACE_TENSION_BOUNDARY, value);
                }
                else {
                    LOG_WARN << "Not a valid surface tension method.";
                }
            }
            else if (paramName == "surfaceTensionViscosity")
            {
                if (model->getSurfaceTensionMethod() == 4) // Jeske 2023
                {
                    SurfaceTension_Jeske2023* stMethod = static_cast<SurfaceTension_Jeske2023*>(model->getSurfaceTensionBase());
                    stMethod->setValue(SurfaceTension_Jeske2023::VISCOSITY_COEFFICIENT, value);
                }
                else {
                    LOG_WARN << "Not a valid surface tension method.";
                }
            }
            else if (paramName == "surfaceTensionViscosityBoundary")
            {
                if (model->getSurfaceTensionMethod() == 4) // Jeske 2023
                {
                    SurfaceTension_Jeske2023* stMethod = static_cast<SurfaceTension_Jeske2023*>(model->getSurfaceTensionBase());
                    stMethod->setValue(SurfaceTension_Jeske2023::VISCOSITY_COEFFICIENT_BOUNDARY, value);
                }
                else {
                    LOG_WARN << "Not a valid surface tension method.";
                }
            }
            else
            {
                LOG_WARN << "Unknown parameter: " << paramName; 
            }
        }
    }
}
