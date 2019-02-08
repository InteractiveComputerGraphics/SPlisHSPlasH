#include "ElasticityBase.h"

using namespace SPH;
using namespace GenParam;

int ElasticityBase::YOUNGS_MODULUS = -1;
int ElasticityBase::POISSON_RATIO = -1;


ElasticityBase::ElasticityBase(FluidModel *model) :
	NonPressureForceBase(model),
	m_youngsModulus(100000.0),
	m_poissonRatio(0.3)
{
}

ElasticityBase::~ElasticityBase(void)
{
}


void ElasticityBase::initParameters()
{
	NonPressureForceBase::initParameters();

	YOUNGS_MODULUS = createNumericParameter("youngsModulus", "Young`s modulus", &m_youngsModulus);
	setGroup(YOUNGS_MODULUS, "Elasticity");
	setDescription(YOUNGS_MODULUS, "Stiffness of the elastic material");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(YOUNGS_MODULUS));
	rparam->setMinValue(0.0);

	POISSON_RATIO = createNumericParameter("poissonsRatio", "Poisson`s ratio", &m_poissonRatio);
	setGroup(POISSON_RATIO, "Elasticity");
	setDescription(POISSON_RATIO, "Ratio of transversal expansion and axial compression");
	rparam = static_cast<RealParameter*>(getParameter(POISSON_RATIO));
	rparam->setMinValue(-1.0 + 1e-4);
	rparam->setMaxValue(0.5 - 1e-4);
}


