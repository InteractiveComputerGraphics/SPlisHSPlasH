#ifndef __SurfaceTensionBaseBase_h__
#define __SurfaceTensionBaseBase_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief Base class for all surface tension methods.
	*/
	class SurfaceTensionBase : public NonPressureForceBase
	{
	protected:
		Real m_surfaceTension;

	public:
		SurfaceTensionBase(FluidModel *model);
		virtual ~SurfaceTensionBase(void);

		Real getSurfaceTension() const { return m_surfaceTension; }
		void setSurfaceTension(Real val) { m_surfaceTension = val; }
	};
}

#endif
