#ifndef FOAMKERNEL_H
#define FOAMKERNEL_H

#define _USE_MATH_DEFINES
#include <math.h>
#include "SPlisHSPlasH/Common.h"

namespace SPH
{
	class FoamKernel
	{
	protected:
		static Real m_radius;
		static Real m_W_zero;
		static Real m_norm_fac;
	public:
		static Real getRadius() { return m_radius; }
		static void setRadius(Real val)
		{
			m_radius = val;
			//// The normalization factor is 1 divided by the integral over the kernel function.
			//// WolframAlpha code for 3D integral over kernel in spherical coordinates: int_0^(2pi) int_0^pi int_0^h ((1-r/h) r^2 sin(theta)) dr dtheta dphi
			//// r^2 sin(theta) is the Jacobian determinant, necessary for the transformation of the volume integral from Cartesian into spherical coordinates.
			//m_norm_fac = 1.0 / ((1.0 / 3.0) * M_PI * val * val * val);
			//m_W_zero = m_norm_fac;
			m_W_zero = 1.0;
		}

	public:
		static Real W(const Real r)
		{
			Real res = 0.0;
			const Real q = r / m_radius;
			if (q <= 1.0)
				//res = m_norm_fac * (1.0 - q);
				res = (static_cast<Real>(1.0) - q);
			return res;
		}

		static Real W(const Vector3r& r)
		{
			return W(r.norm());
		}

		static Real W_zero()
		{
			return m_W_zero;
		}
	};
}

#endif
