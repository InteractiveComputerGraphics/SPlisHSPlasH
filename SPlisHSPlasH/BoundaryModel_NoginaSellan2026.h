#ifndef __BoundaryModel_NoginaSellan2026_h__
#define __BoundaryModel_NoginaSellan2026_h__

#include "BoundaryModel_Koschier2017.h"


namespace SPH
{
	/** \brief The boundary model stores the information required for boundary handling
	* using the tube maps of Nogina and Sellán 2026 [NS26].
	*
	* A tube map is a density map in the sense of [KB17], i.e., it stores the signed distance
	* to the boundary (field 0) and the boundary density (field 1) on a cubic Lagrange grid.
	* In contrast to [KB17], the boundary density is not obtained by a numerical volume integral
	* but by a closed-form 1D integral of the projected cubic spline kernel along the normal
	* direction (optionally with a correction using the mean curvature (order 1) or the mean and Gaussian
	* curvature (order 2) of the surface). The map is generated in SimulatorBase::initTubeMap().
	*
	* Since tube maps are evaluated exactly like density maps, the model shares the interface of
	* BoundaryModel_Koschier2017 and is handled by all pressure and non-pressure solvers in the
	* same way (see Simulation::usesDensityMaps()).
	*
	* References:
	* - [KB17] Dan Koschier and Jan Bender. Density maps for improved SPH boundary handling. In ACM SIGGRAPH/Eurographics Symposium on Computer Animation, 1-10. July 2017. URL: http://dx.doi.org/10.1145/3099564.3099565
	* - [NS26] Daria Nogina and Silvia Sellán. Tube Maps: Fast SPH Boundary Handling in Tubular Coordinates. In SIGGRAPH 2026 Conference Papers. 2026.
	*/
	class BoundaryModel_NoginaSellan2026 : public BoundaryModel_Koschier2017
	{
		public:
			BoundaryModel_NoginaSellan2026() : BoundaryModel_Koschier2017() {}
			virtual ~BoundaryModel_NoginaSellan2026() {}
	};
}

#endif
