#ifndef __TubeMapDensity_h__
#define __TubeMapDensity_h__

namespace SPH
{
	/** \brief Closed-form boundary density of the tube maps of Nogina and Sellán 2026 [NS26].
	*
	* The boundary density at a point with signed distance d to the surface is computed by a 1D
	* integral along the normal direction t in [-h, h] of the cubic spline kernel projected onto the
	* tangent plane. The integrand is weighted by
	* - the linear fall-off Gamma(s) = 1 - factor * s / h, s = d + t, which smoothly cuts off the
	*   density at s = h / factor, and
	* - an area correction which accounts for the curvature of the surface. The order of the tube map
	*   determines the correction:
	*   - order 0: planar boundary, A(s) = 1
	*   - order 1: mean curvature H, A(s) = 1 + 0.5 H d - 0.5 H s
	*   - order 2: mean curvature H and Gaussian curvature K, A(s) = 1 + 0.5 H d - 0.5 H s + (2K + 0.5 H^2) s^2,
	*     together with a second kernel moment term.
	*   The curvatures are clamped to the range in which the expansion is valid (principal curvatures
	*   |k| <= 2/h, saddle-type Gaussian curvature K >= -1/(4h^2)) since discrete curvatures of fine or
	*   noisy meshes can be orders of magnitude larger.
	*
	* All polynomial integrals are evaluated exactly.
	*
	* References:
	* - [NS26] Daria Nogina and Silvia Sellán. Tube Maps: Fast SPH Boundary Handling in Tubular Coordinates. In SIGGRAPH 2026 Conference Papers. 2026.
	*/
	class TubeMapDensity
	{
	public:
		/** Compute the (unnormalized) boundary density.
		* @param d signed distance to the surface
		* @param H mean curvature at the closest surface point, H = -(k1 + k2), i.e., negative for convex surfaces (sphere of radius R: H = -2/R)
		* @param K Gaussian curvature at the closest surface point
		* @param h support radius of the cubic spline kernel
		* @param factor fall-off factor of Gamma
		* @param order order of the tube map (0, 1 or 2)
		*/
		static double boundaryDensity(const double d, const double H, const double K, const double h, const double factor, const unsigned int order);
	};
}

#endif
