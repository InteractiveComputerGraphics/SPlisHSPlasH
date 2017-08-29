#include "VorticityBase.h"

using namespace SPH;

VorticityBase::VorticityBase(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_vorticityCoeff = 0.01;
}

VorticityBase::~VorticityBase(void)
{
}



