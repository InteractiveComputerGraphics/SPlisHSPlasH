#include "ViscosityBase.h"

using namespace SPH;

ViscosityBase::ViscosityBase(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_viscosity = 0.02;
}

ViscosityBase::~ViscosityBase(void)
{
}



