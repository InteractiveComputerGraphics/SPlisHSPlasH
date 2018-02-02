#include "SurfaceTensionBase.h"

using namespace SPH;
using namespace GenParam;

int SurfaceTensionBase::SURFACE_TENSION_COEFFICIENT = -1;

SurfaceTensionBase::SurfaceTensionBase(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_surfaceTension = 0.05;
}

SurfaceTensionBase::~SurfaceTensionBase(void)
{
}

void SurfaceTensionBase::initParameters()
{
	NonPressureForceBase::initParameters();

	SURFACE_TENSION_COEFFICIENT = createNumericParameter("surfaceTension", "Surface tension coefficient", &m_surfaceTension);
	setGroup(SURFACE_TENSION_COEFFICIENT, "Surface tension");
	setDescription(SURFACE_TENSION_COEFFICIENT, "Coefficient for the surface tension computation");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(SURFACE_TENSION_COEFFICIENT));
	rparam->setMinValue(0.0);
}


