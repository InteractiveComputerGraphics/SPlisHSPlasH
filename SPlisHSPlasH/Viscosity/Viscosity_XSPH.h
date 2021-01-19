#ifndef __Viscosity_XSPH_h__
#define __Viscosity_XSPH_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ViscosityBase.h"

namespace SPH
{
	/** \brief This class implements the XSPH method descibed by
	* Schechter and Bridson [SB12].
	*
	* References:
	* - [SB12] Hagit Schechter and Robert Bridson. Ghost sph for animating water. ACM Trans. Graph., 31(4):61:1-61:8, July 2012. URL: http://doi.acm.org/10.1145/2185520.2185557
	*/
	class Viscosity_XSPH : public ViscosityBase
	{
	protected:
		Real m_boundaryViscosity;

		virtual void initParameters();

	public:
		static int VISCOSITY_COEFFICIENT_BOUNDARY;

		Viscosity_XSPH(FluidModel *model);
		virtual ~Viscosity_XSPH(void);

		virtual void step();
		virtual void reset();

		static NonPressureForceBase* creator(FluidModel* model) { return new Viscosity_XSPH(model); }
	};
}

#endif
