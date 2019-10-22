#include "AnimationField.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "FluidModel.h"
#include "Simulation.h"
#include "extern/tinyexpr/tinyexpr.h"
#include "Utilities/Logger.h"

using namespace SPH;

AnimationField::AnimationField(
	const std::string &particleFieldName,
	const Vector3r &pos, const Matrix3r & rotation, const Vector3r &scale,
	const std::string expression[3], const unsigned int type)
	: m_particleFieldName(particleFieldName)
	, m_x(pos)
	, m_rotation(rotation)
	, m_scale(scale)
	, m_type(type)
	, m_startTime(0)
	, m_endTime(REAL_MAX)
{
	for (int i = 0; i < 3; i++)
		m_expression[i] = expression[i];
}

AnimationField::~AnimationField(void)
{
}

void AnimationField::reset()
{
}

double getTime()
{
	return TimeManager::getCurrent()->getTime();
}

void AnimationField::step()
{
	Simulation *sim = Simulation::getCurrent();
	TimeManager *tm = TimeManager::getCurrent();
	const Real t = tm->getTime();
	const Real dt = tm->getTimeStepSize();

	if (t >= m_startTime && t <= m_endTime)
	{
		// animate particles		
		const unsigned int nModels = sim->numberOfFluidModels();
		for (unsigned int m = 0; m < nModels; m++)
		{
			FluidModel *fm = sim->getFluidModel(m);
			const unsigned int numParticles = fm->numActiveParticles();

			// find angular velocity field
			const FieldDescription *particleField = nullptr;
			for (unsigned int j = 0; j < fm->numberOfFields(); j++)
			{
				const FieldDescription &field = fm->getField(j);
				if (field.name == m_particleFieldName)
				{
					particleField = &field;
					break;
				}
			}

			if (particleField == nullptr)
				continue;

			#pragma omp parallel for schedule(static) default(shared)
			for (int i = 0; i < (int)numParticles; i++)
			{
				const Vector3r &xi = fm->getPosition(i);
				const Vector3r &vi = fm->getVelocity(i);
				if (inShape(m_type, xi, m_x, m_rotation, m_scale))
				{
					Eigen::Map<Vector3r> value((Real*) particleField->getFct(i));
 					te_variable vars[] = { {"t", &t}, {"dt", &dt}, 
										   {"x", &xi[0]}, {"y", &xi[1]}, {"z", &xi[2]},
										   {"vx", &vi[0]}, {"vy", &vi[1]}, {"vz", &vi[2]},
										   {"valuex", &value[0]}, {"valuey", &value[1]}, {"valuez", &value[2]},
										 };
 					const int numVars = 11;
 					int err;

					//////////////////////////////////////////////////////////////////////////
					// v_x
					//////////////////////////////////////////////////////////////////////////
					if (m_expression[0] != "")
					{
						te_expr *expr_vx = te_compile(m_expression[0].c_str(), vars, numVars, &err);
						if (expr_vx)
							value[0] = te_eval(expr_vx);
						te_free(expr_vx);

						if (err != 0)
							LOG_ERR << "Animation field: expression for x is wrong.";
					}

					//////////////////////////////////////////////////////////////////////////
					// v_y
					//////////////////////////////////////////////////////////////////////////
					if (m_expression[1] != "")
					{
						te_expr *expr_vy = te_compile(m_expression[1].c_str(), vars, numVars, &err);
						if (expr_vy)
							value[1] = te_eval(expr_vy);
						te_free(expr_vy);

						if (err != 0)
							LOG_ERR << "Animation field: expression for y is wrong.";
					}

					//////////////////////////////////////////////////////////////////////////
					// v_z
					//////////////////////////////////////////////////////////////////////////
					if (m_expression[2] != "")
					{
						te_expr *expr_vz = te_compile(m_expression[2].c_str(), vars, numVars, &err);
						if (expr_vz)
							value[2] = te_eval(expr_vz);
						te_free(expr_vz);

						if (err != 0)
							LOG_ERR << "Animation field: expression for z is wrong.";
					}
				}
			}
		}
	}
}

