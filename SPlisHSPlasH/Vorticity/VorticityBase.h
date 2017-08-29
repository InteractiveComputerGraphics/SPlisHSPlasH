#ifndef __VorticityBaseBase_h__
#define __VorticityBaseBase_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief Base class for all vorticity methods.
	*/
	class VorticityBase : public NonPressureForceBase
	{
	protected:
		Real m_vorticityCoeff;

	public:
		VorticityBase(FluidModel *model);
		virtual ~VorticityBase(void);

		Real getVorticityCoeff() const { return m_vorticityCoeff; }
		void setVorticityCoeff(Real val) { m_vorticityCoeff = val; }
	};
}

#endif
