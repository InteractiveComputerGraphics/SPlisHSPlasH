#include "SurfaceTensionBase.h"

using namespace SPH;
using namespace GenParam;

int SurfaceTensionBase::SURFACE_TENSION = -1;
int SurfaceTensionBase::SURFACE_TENSION_BOUNDARY = -1;

SurfaceTensionBase::SurfaceTensionBase(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_surfaceTension = 0.05;
	m_surfaceTensionBoundary = 0.01;
}

SurfaceTensionBase::~SurfaceTensionBase(void)
{
}

void SurfaceTensionBase::initParameters()
{
	NonPressureForceBase::initParameters();

	SURFACE_TENSION = createNumericParameter("surfaceTension", "Surface tension coefficient", &m_surfaceTension);
	setGroup(SURFACE_TENSION, "Surface tension");
	setDescription(SURFACE_TENSION, "Coefficient for the surface tension computation");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(SURFACE_TENSION));
	rparam->setMinValue(0.0);

	SURFACE_TENSION_BOUNDARY = createNumericParameter("surfaceTensionBoundary", "Boundary surface tension coefficient", &m_surfaceTensionBoundary);
	setGroup(SURFACE_TENSION_BOUNDARY, "Surface tension");
	setDescription(SURFACE_TENSION_BOUNDARY, "Coefficient for the surface tension computation at the boundary");
	rparam = static_cast<RealParameter*>(getParameter(SURFACE_TENSION_BOUNDARY));
	rparam->setMinValue(0.0);
}


