#include "VorticityBase.h"

using namespace SPH;
using namespace GenParam;

int VorticityBase::VORTICITY_COEFFICIENT = -1;

VorticityBase::VorticityBase(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_vorticityCoeff = 0.01;
}

VorticityBase::~VorticityBase(void)
{
}

void VorticityBase::initParameters()
{
	NonPressureForceBase::initParameters();

	VORTICITY_COEFFICIENT = createNumericParameter("vorticity", "Vorticity transfer coefficient", &m_vorticityCoeff);
	setGroup(VORTICITY_COEFFICIENT, "Vorticity");
	setDescription(VORTICITY_COEFFICIENT, "Coefficient for the vorticity force computation");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(VORTICITY_COEFFICIENT));
	rparam->setMinValue(0.0);
}




