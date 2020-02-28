#include "ElasticityBase.h"

using namespace SPH;
using namespace GenParam;

int ElasticityBase::YOUNGS_MODULUS = -1;
int ElasticityBase::POISSON_RATIO = -1;


ElasticityBase::ElasticityBase(FluidModel *model) :
	NonPressureForceBase(model),
	m_youngsModulus(static_cast<Real>(100000.0)),
	m_poissonRatio(static_cast<Real>(0.3))
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
	rparam->setMinValue(static_cast<Real>(-1.0 + 1e-4));
	rparam->setMaxValue(static_cast<Real>(0.5 - 1e-4));
}


