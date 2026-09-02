#ifndef __DynamicParameterSystem_h__
#define __DynamicParameterSystem_h__

#include "Common.h"
#include "extern/tinyexpr/tinyexpr.h"
#include <vector>
#include <string>

namespace SPH
{
    //class FluidModel;
    //class ViscosityBase;

    struct TinyExprDeleter {
        void operator()(te_expr* expr) const {
            if (expr) {
                te_free(expr);
            }
        }
    };

    struct ParameterSchedule
    {
        std::string fluidId;
        std::string parameter_name;
        std::string expression;
        std::vector<std::pair<Real, Real>> timeline;
        int current_timeline_index;
        bool use_expression;
        bool use_step_function;
        std::unique_ptr<te_expr, TinyExprDeleter> compiled_expr;
        Real current_value;
        Real default_value;
        bool active;

        ParameterSchedule() : compiled_expr(nullptr), current_value(0.0), default_value(0.0), use_expression(false), active(true), current_timeline_index(-1) {}
        /*~ParameterSchedule() 
        {
            if (compiled_expr) 
            {
                te_free(compiled_expr);
                compiled_expr = nullptr;
            }
        }*/
    };

    class DynamicParameterSystem
    {
    private:
        std::vector<ParameterSchedule> m_schedules;

        double m_t_double;
        double m_dt_double;

    public:
        DynamicParameterSystem();
        virtual ~DynamicParameterSystem();

        void step();
        void reset();

        // Add schedules
        bool addTimelineSchedule(const std::string& fluidId, const std::string& paramName, const std::vector<std::pair<Real, Real>>& timeline,
            Real defaultValue = 1.0, bool stepFunction = false);

        bool addExpressionSchedule(const std::string& fluidId, const std::string& paramName,
            const std::string& expression, Real defaultValue = 1.0);

        size_t numberOfSchedules() const { return m_schedules.size(); }
        ParameterSchedule& getSchedule(const unsigned int i) { return m_schedules[i]; }
        const ParameterSchedule& getSchedule(const unsigned int i) const { return m_schedules[i]; }

    private:
        bool compileExpression(ParameterSchedule& schedule);
        Real interpolateTimeline(ParameterSchedule& schedule, Real currentTime);
        Real stepTimeline(ParameterSchedule& schedule, Real currentTime);
        void applyParameterValue(const std::string& fluidId, const std::string& paramName, Real value);
    };
}

#endif
