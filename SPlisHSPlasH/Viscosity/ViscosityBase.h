#ifndef __ViscosityBaseBase_h__
#define __ViscosityBaseBase_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief Base class for all viscosity methods.
	*/
	class ViscosityBase : public NonPressureForceBase
	{
	protected:
		Real m_viscosity;

		virtual void initParameters();

	public:
		static int VISCOSITY_COEFFICIENT;

		ViscosityBase(FluidModel *model);
		virtual ~ViscosityBase(void);
	};
}

#endif
