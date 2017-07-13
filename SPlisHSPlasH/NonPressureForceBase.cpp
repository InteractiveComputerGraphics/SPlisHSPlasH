#include "NonPressureForceBase.h"

using namespace SPH;

NonPressureForceBase::NonPressureForceBase(FluidModel *model)
{
	m_model = model;
}

NonPressureForceBase::~NonPressureForceBase(void)
{
}



