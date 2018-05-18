#include "ViscosityBase.h"

using namespace SPH;
using namespace GenParam;

int ViscosityBase::VISCOSITY_COEFFICIENT = -1;


ViscosityBase::ViscosityBase(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_viscosity = 0.01;
}

ViscosityBase::~ViscosityBase(void)
{
}


void ViscosityBase::initParameters()
{
	NonPressureForceBase::initParameters();

	VISCOSITY_COEFFICIENT = createNumericParameter("viscosity", "Viscosity coefficient", &m_viscosity);
	setGroup(VISCOSITY_COEFFICIENT, "Viscosity");
	setDescription(VISCOSITY_COEFFICIENT, "Coefficient for the viscosity force computation");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(VISCOSITY_COEFFICIENT));
	rparam->setMinValue(0.0);
}


