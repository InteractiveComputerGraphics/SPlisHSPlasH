#include "ElasticityBase.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Utilities/Logger.h"

using namespace SPH;
using namespace GenParam;

int ElasticityBase::YOUNGS_MODULUS = -1;
int ElasticityBase::POISSON_RATIO = -1;
int ElasticityBase::FIXED_BOX_MIN = -1;
int ElasticityBase::FIXED_BOX_MAX = -1;



ElasticityBase::ElasticityBase(FluidModel *model) :
	NonPressureForceBase(model),
	m_youngsModulus(static_cast<Real>(100000.0)),
	m_poissonRatio(static_cast<Real>(0.3))
{
	m_fixedBoxMin.setZero();
	m_fixedBoxMax.setZero();
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

	ParameterBase::GetVecFunc<Real> getFct = [&]()-> Real* { return m_fixedBoxMin.data(); };
	ParameterBase::SetVecFunc<Real> setFct = [&](Real* val)
	{
		m_fixedBoxMin = Vector3r(val[0], val[1], val[2]);
		determineFixedParticles();
	};
	FIXED_BOX_MIN = createVectorParameter("fixedBoxMin", "Fixed box min", 3u, getFct, setFct);
	setGroup(FIXED_BOX_MIN, "Elasticity");
	setDescription(FIXED_BOX_MIN, "Minimum point of box of which contains the fixed particles.");
	getParameter(FIXED_BOX_MIN)->setReadOnly(true);


	ParameterBase::GetVecFunc<Real> getFct2 = [&]()-> Real* { return m_fixedBoxMax.data(); };
	ParameterBase::SetVecFunc<Real> setFct2 = [&](Real* val)
	{
		m_fixedBoxMax = Vector3r(val[0], val[1], val[2]);
		determineFixedParticles();
	};
	FIXED_BOX_MAX = createVectorParameter("fixedBoxMax", "Fixed box max", 3u, getFct2, setFct2);
	setGroup(FIXED_BOX_MAX, "Elasticity");
	setDescription(FIXED_BOX_MAX, "Maximum point of box of which contains the fixed particles.");
	getParameter(FIXED_BOX_MAX)->setReadOnly(true);
}

/** Mark all particles in the bounding box as fixed.
*/
void ElasticityBase::determineFixedParticles()
{
	const unsigned int numParticles = m_model->numActiveParticles();

	if (!m_fixedBoxMin.isZero() || !m_fixedBoxMax.isZero())
	{
		int counter = 0;
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& x = m_model->getPosition(i);
			if ((x[0] > m_fixedBoxMin[0]) && (x[1] > m_fixedBoxMin[1]) && (x[2] > m_fixedBoxMin[2]) &&
				(x[0] < m_fixedBoxMax[0]) && (x[1] < m_fixedBoxMax[1]) && (x[2] < m_fixedBoxMax[2]))
			{
				m_model->setParticleState(i, ParticleState::Fixed);
			}
			if (m_model->getParticleState(i) == ParticleState::Fixed)
				counter++;
		}
		LOG_INFO << "Fixed particles: " << counter;
	}
}


