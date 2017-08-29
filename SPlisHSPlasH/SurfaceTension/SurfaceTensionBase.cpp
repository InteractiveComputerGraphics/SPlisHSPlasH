#include "SurfaceTensionBase.h"

using namespace SPH;

SurfaceTensionBase::SurfaceTensionBase(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_surfaceTension = 0.05;
}

SurfaceTensionBase::~SurfaceTensionBase(void)
{
}



