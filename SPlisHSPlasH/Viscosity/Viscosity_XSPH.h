#ifndef __Viscosity_XSPH_h__
#define __Viscosity_XSPH_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ViscosityBase.h"

namespace SPH
{
	/** \brief This class implements the XSPH method descibed by
	* Schechter and Bridson \cite Schechter:2012.
	*/
	class Viscosity_XSPH : public ViscosityBase
	{
	public:
		Viscosity_XSPH(FluidModel *model);
		virtual ~Viscosity_XSPH(void);

		virtual void step();
		virtual void reset();
	};
}

#endif
