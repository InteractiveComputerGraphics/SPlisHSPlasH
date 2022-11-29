#ifndef __XSPH_h__
#define __XSPH_h__

#include "Common.h"
#include "FluidModel.h"
#include "NonPressureForceBase.h"

namespace SPH
{
	/** \brief This class implements the XSPH method descibed by
	* J. J. Monaghan [Mon92].
	*
	* References:
	* - [Mon92] J. J. Monaghan. Smoothed Particle Hydrodynamics. Annual Review of Astronomy and Astrophysics, 1992, 30, 543-574. 
	* URL: https://www.annualreviews.org/doi/10.1146/annurev.aa.30.090192.002551
	*/
	class XSPH : public NonPressureForceBase
	{
	protected:
		Real m_fluidCoefficient;
		Real m_boundaryCoefficient;

		virtual void initParameters();

	public:
		static int FLUID_COEFFICIENT;
		static int BOUNDARY_COEFFICIENT;

		XSPH(FluidModel *model);
		virtual ~XSPH(void);

		virtual void step();
		virtual void reset();
	};
}

#endif
