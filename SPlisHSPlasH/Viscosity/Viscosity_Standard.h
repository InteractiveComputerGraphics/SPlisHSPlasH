#ifndef __Viscosity_Standard_h__
#define __Viscosity_Standard_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ViscosityBase.h"

namespace SPH
{
	/** \brief This class implements the standard method for viscosity descibed e.g. by
	* Ihmsen et al. [IOS+14].\n\n
	* The method evaluates the term \f$\nu \nabla^2 \mathbf{v}\f$ and uses an approximation 
	* of the kernel Laplacian to improve the stability. This approximation is given in [IOS+14].
	*
	* References:
	* - [IOS+14] Markus Ihmsen, Jens Orthmann, Barbara Solenthaler, Andreas Kolb, and Matthias Teschner. SPH Fluids in Computer Graphics. In Sylvain Lefebvre and Michela Spagnuolo, editors, Eurographics 2014 - State of the Art Reports. The Eurographics Association, 2014. URL: http://dx.doi.org/10.2312/egst.20141034
	*/
	class Viscosity_Standard : public ViscosityBase
	{
	protected:
		Real m_boundaryViscosity;

		virtual void initParameters();

	public:
		static int VISCOSITY_COEFFICIENT_BOUNDARY;

		Viscosity_Standard(FluidModel *model);
		virtual ~Viscosity_Standard(void);

		virtual void step();
		virtual void reset();
	};
}

#endif
