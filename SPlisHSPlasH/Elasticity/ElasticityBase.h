#ifndef __ElasticityBase_h__
#define __ElasticityBase_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief Base class for all elasticity methods.
	*/
	class ElasticityBase : public NonPressureForceBase
	{
	protected:
		Real m_youngsModulus;
		Real m_poissonRatio;
		Vector3r m_fixedBoxMin;
		Vector3r m_fixedBoxMax;

		virtual void initParameters();
		void determineFixedParticles();

	public:
		static int YOUNGS_MODULUS;
		static int POISSON_RATIO;
		static int FIXED_BOX_MIN;
		static int FIXED_BOX_MAX;

		ElasticityBase(FluidModel *model);
		virtual ~ElasticityBase(void);
	};
}

#endif
