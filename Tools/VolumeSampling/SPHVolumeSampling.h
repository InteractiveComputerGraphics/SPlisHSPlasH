#pragma once

#include "SPHSamplingBase.h"
#include "SPlisHSPlasH/NeighborhoodSearch.h"


namespace SPH
{
	/** \brief This class implements the sampling method introduced
	* by Kugelstadt et al. [KBF+21].
	*
	* References:
	* - [KBF+21] Tassilo Kugelstadt, Jan Bender, José Antonio Fernández-Fernández, Stefan Rhys Jeske, Fabian Löschner, Andreas Longva.
	* Fast Corotated Elastic SPH Solids with Implicit Zero-Energy Mode Control. Proceedings of the ACM on Computer Graphics and
	* Interactive Techniques, 2021. URL: http://dx.doi.org/10.1145/3480142
	*/
	class SPHVolumeSampling : public SPHSamplingBase
	{
	protected:
		virtual void initSPHOptimization();
		virtual void step(Real &avg_density_error);	
		virtual void generateSamples();

	public:
		SPHVolumeSampling();
		~SPHVolumeSampling();

	};
}
