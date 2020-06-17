#ifndef __DragForce_Gissler2017_h__
#define __DragForce_Gissler2017_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "DragBase.h"

namespace SPH
{
	/** \brief This class implements the drag force computation introduced
	* by Gissler et al. [GBP+17].
	*
	* References:
	* - [GPB+17] Christoph Gissler, Stefan Band, Andreas Peer, Markus Ihmsen, and Matthias Teschner. Approximate air-fluid interactions for SPH. In Virtual Reality Interactions and Physical Simulations, 1-10. April 2017. URL: http://dx.doi.org/10.2312/vriphys.20171081
	*/
	class DragForce_Gissler2017 : public DragBase
	{
	protected:
		const Real rho_a = static_cast<Real>(1.2041);
		const Real sigma = static_cast<Real>(0.0724);
		const Real mu_l = static_cast<Real>(0.00102);
		const Real C_F = static_cast<Real>(1.0/3.0);
		const Real C_k = static_cast<Real>(8.0);
		const Real C_d = static_cast<Real>(5.0);
		const Real C_b = static_cast<Real>(0.5);
		const Real mu_a = static_cast<Real>(0.00001845);

	public:
		DragForce_Gissler2017(FluidModel *model);
		virtual ~DragForce_Gissler2017(void);

		virtual void step();
		virtual void reset();
	};
}

#endif
