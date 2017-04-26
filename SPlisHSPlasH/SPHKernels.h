#ifndef SPHKERNELS_H
#define SPHKERNELS_H

#define _USE_MATH_DEFINES
#include <math.h>
#include "Common.h"
#include <algorithm>

namespace SPH
{
	/** \brief Cubic spline kernel.
	*/
	class CubicKernel
	{
	protected:
		static Real m_radius;
		static Real m_k;
		static Real m_l;
		static Real m_W_zero;
	public:
		static Real getRadius() { return m_radius; }
		static void setRadius(Real val)
		{
			m_radius = val;
			static const Real pi = static_cast<Real>(M_PI);

			const Real h3 = m_radius*m_radius*m_radius;
			m_k = 8.0 / (pi*h3);
			m_l = 48.0 / (pi*h3);
			m_W_zero = W(Vector3r::Zero());
		}

	public:
		static Real W(const Real r)
		{
			Real res = 0.0;
			const Real q = r / m_radius;
			if (q <= 1.0)
			{
				if (q <= 0.5)
				{
					const Real q2 = q*q;
					const Real q3 = q2*q;
					res = m_k * (6.0*q3 - 6.0*q2 + 1.0);
				}
				else
				{
					res = m_k * (2.0*pow(1.0 - q, 3));
				}
			}
			return res;
		}

		static Real W(const Vector3r &r)
		{
			return W(r.norm());
		}

		static Vector3r gradW(const Vector3r &r)
		{
			Vector3r res;
			const Real rl = r.norm();
			const Real q = rl / m_radius;
			if (q <= 1.0)
			{
				if (rl > 1.0e-6)
				{
					const Vector3r gradq = r * ((Real) 1.0 / (rl*m_radius));
					if (q <= 0.5)
					{
						res = m_l*q*((Real) 3.0*q - (Real) 2.0)*gradq;
					}
					else
					{
						const Real factor = 1.0 - q;
						res = m_l*(-factor*factor)*gradq;
					}
				}
			}
			else
				res.setZero();

			return res;
		}

		static Real W_zero()
		{
			return m_W_zero;
		}
	};

	/** \brief Poly6 kernel.
	*/
	class Poly6Kernel
	{
	protected:
		static Real m_radius;
		static Real m_k;
		static Real m_l;
		static Real m_m;
		static Real m_W_zero;
	public:
		static Real getRadius() { return m_radius; }
		static void setRadius(Real val)
		{
			m_radius = val;
			static const Real pi = static_cast<Real>(M_PI);
			m_k = 315.0 / (64.0*pi*pow(m_radius, 9));
			m_l = -945.0 / (32.0*pi*pow(m_radius, 9));
			m_m = m_l;
			m_W_zero = W(Vector3r::Zero());
		}

	public:

		/** 
		* W(r,h) = (315/(64 pi h^9))(h^2-|r|^2)^3
		*        = (315/(64 pi h^9))(h^2-r*r)^3
		*/
		static Real W(const Real r)
		{
			Real res = 0.0;
			const Real r2 = r*r;
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				res = pow(radius2 - r2, 3)*m_k;
			}
			return res;
		}

		static Real W(const Vector3r &r)
		{
			Real res = 0.0;
			const Real r2 = r.squaredNorm();
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				res = pow(radius2 - r2, 3)*m_k;
			}
			return res;
		}


		/**
		* grad(W(r,h)) = r(-945/(32 pi h^9))(h^2-|r|^2)^2
		*              = r(-945/(32 pi h^9))(h^2-r*r)^2
		*/
		static Vector3r gradW(const Vector3r &r)
		{
			Vector3r res;
			const Real r2 = r.squaredNorm();
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				Real tmp = radius2 - r2;
				res = m_l * tmp * tmp*r;
			}
			else
				res.setZero();

			return res;
		}

		/**
		* laplacian(W(r,h)) = (-945/(32 pi h^9))(h^2-|r|^2)(-7|r|^2+3h^2)
		*                   = (-945/(32 pi h^9))(h^2-r*r)(3 h^2-7 r*r)
		*/
		static Real laplacianW(const Vector3r &r)
		{
			Real res;
			const Real r2 = r.squaredNorm();
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				Real tmp = radius2 - r2;
				Real tmp2 = 3 * radius2 - 7 * r2;
				res = m_m * tmp  * tmp2;
			}
			else
				res = (Real) 0.;

			return res;
		}

		static Real W_zero()
		{
			return m_W_zero;
		}
	};

	/** \brief Spiky kernel.
	*/
	class SpikyKernel
	{
	protected:
		static Real m_radius;
		static Real m_k;
		static Real m_l;
		static Real m_W_zero;
	public:
		static Real getRadius() { return m_radius; }
		static void setRadius(Real val)
		{
			m_radius = val;
			const Real radius6 = pow(m_radius, 6);
			static const Real pi = static_cast<Real>(M_PI);
			m_k = 15.0 / (pi*radius6);
			m_l = -45.0 / (pi*radius6);
			m_W_zero = W(Vector3r::Zero());
		}

	public:

		/**
		* W(r,h) = 15/(pi*h^6) * (h-r)^3
		*/
		static Real W(const Real r)
		{
			Real res = 0.0;
			const Real r2 = r*r;
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				const Real hr3 = pow(m_radius - r, 3);
				res = m_k * hr3;
			}
			return res;
		}

		static Real W(const Vector3r &r)
		{
			Real res = 0.0;
			const Real r2 = r.squaredNorm();
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				const Real hr3 = pow(m_radius - sqrt(r2), 3);
				res = m_k * hr3;
			}
			return res;
		}


		/**
		* grad(W(r,h)) = -r(45/(pi*h^6) * (h-r)^2)
		*/
		static Vector3r gradW(const Vector3r &r)
		{
			Vector3r res;
			const Real r2 = r.squaredNorm();
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				const Real r_l = sqrt(r2);
				const Real hr = m_radius - r_l;
				const Real hr2 = hr*hr;
				res = m_l * hr2 * r * (static_cast<Real>(1.0) / r_l);
			}
			else
				res.setZero();

			return res;
		}

		static Real W_zero()
		{
			return m_W_zero;
		}
	};

	/** \brief Cohesion kernel used for the surface tension method of Akinci el al. \cite Akinci:2013.
	*/
	class CohesionKernel
	{
	protected:
		static Real m_radius;
		static Real m_k;
		static Real m_c;
		static Real m_W_zero;
	public:
		static Real getRadius() { return m_radius; }
		static void setRadius(Real val)
		{
			m_radius = val;
			static const Real pi = static_cast<Real>(M_PI);
			m_k = 32. / (pi*pow(m_radius, 9));
			m_c = pow(m_radius, 6) / 64.0;
			m_W_zero = W(Vector3r::Zero());
		}

	public:

		/**
		* W(r,h) = (32/(pi h^9))(h-r)^3*r^3					if h/2 < r <= h
		*          (32/(pi h^9))(2*(h-r)^3*r^3 - h^6/64		if 0 < r <= h/2
		*/
		static Real W(const Real r)
		{
			Real res = 0.0;
			const Real r2 = r*r;
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				const Real r1 = sqrt(r2);
				const Real r3 = r2*r1;
				if (r1 > 0.5*m_radius)
					res = m_k*pow(m_radius - r1, 3)*r3;
				else
					res = m_k*2.0*pow(m_radius - r1, 3)*r3 - m_c;

			}
			return res;
		}

		static Real W(const Vector3r &r)
		{
			Real res = 0.0;
			const Real r2 = r.squaredNorm();
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				const Real r1 = sqrt(r2);
				const Real r3 = r2*r1;
				if (r1 > 0.5*m_radius)
					res = m_k*pow(m_radius - r1, 3)*r3;
				else
					res = m_k*2.0*pow(m_radius - r1, 3)*r3 - m_c;

			}
			return res;
		}

		static Real W_zero()
		{
			return m_W_zero;
		}
	};

	/** \brief Adhesion kernel used for the surface tension method of Akinci el al. \cite Akinci:2013.
	*/
	class AdhesionKernel
	{
	protected:
		static Real m_radius;
		static Real m_k;
		static Real m_W_zero;
	public:
		static Real getRadius() { return m_radius; }
		static void setRadius(Real val)
		{
			m_radius = val;
			m_k = 0.007 / pow(m_radius, 3.25);
			m_W_zero = W(Vector3r::Zero());
		}

	public:

		/**
		* W(r,h) = (0.007/h^3.25)(-4r^2/h + 6r -2h)^0.25					if h/2 < r <= h
		*/
		static Real W(const Real r)
		{
			Real res = 0.0;
			const Real r2 = r*r;
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				const Real r = sqrt(r2);
				if (r > 0.5*m_radius)
					res = m_k*pow(-4.0*r2 / m_radius + 6.0*r - 2.0*m_radius, 0.25);
			}
			return res;
		}

		static Real W(const Vector3r &r)
		{
			Real res = 0.0;
			const Real r2 = r.squaredNorm();
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				const Real r = sqrt(r2);
				if (r > 0.5*m_radius)
					res = m_k*pow(-4.0*r2 / m_radius + 6.0*r - 2.0*m_radius, 0.25);
			}
			return res;
		}

		static Real W_zero()
		{
			return m_W_zero;
		}
	};


	/** \brief Precomputed kernel which is based on a lookup table as described by Bender and Koschier \cite Bender:2015, \cite Bender2017.
	*
	* The lookup tables can be used in combination with any kernel. 
	*/
	template<typename KernelType, unsigned int resolution = 1000u>
	class PrecomputedKernel
	{
	protected:
		static Real m_W[resolution];
		static Real m_gradW[resolution + 1];
		static Real m_radius;
		static Real m_invStepSize;
		static Real m_W_zero;
	public:
		static Real getRadius() { return m_radius; }
		static void setRadius(Real val)
		{
			m_radius = val;
			KernelType::setRadius(val);
			const Real stepSize = m_radius / (Real)resolution;
			m_invStepSize = 1.0 / stepSize;
			for (unsigned int i = 0; i < resolution; i++)
			{
				const Real posX = stepSize * (Real)i;		// Store kernel values in the middle of an interval
				m_W[i] = KernelType::W(posX);
				KernelType::setRadius(val);
				if (posX > 1.0e-9)
					m_gradW[i] = KernelType::gradW(Vector3r(posX, 0.0, 0.0))[0] / posX;
				else
					m_gradW[i] = 0.0;
			}
			m_gradW[resolution] = 0.0;
			m_W_zero = W(0.0);
		}

	public:
		static Real W(const Vector3r &r)
		{
			Real res = 0.0;
			const Real r2 = r.squaredNorm();
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				const Real r = sqrt(r2);
				const unsigned int pos = std::min<unsigned int>((unsigned int)(r * m_invStepSize), resolution);
				res = m_W[pos];
			}
			return res;
		}

		static Real W(const Real r)
		{
			Real res = 0.0;
			if (r <= m_radius)
			{
				const unsigned int pos = std::min<unsigned int>((unsigned int)(r * m_invStepSize), resolution);
				res = m_W[pos];
			}
			return res;
		}

		static Vector3r gradW(const Vector3r &r)
		{
			Vector3r res;
			const Real r2 = r.squaredNorm();
			const Real radius2 = m_radius*m_radius;
			if (r2 <= radius2)
			{
				const Real rl = sqrt(r2);
				const unsigned int pos = std::min<unsigned int>((unsigned int)(rl * m_invStepSize), resolution);
				res = m_gradW[pos] * r;
			}
			else
				res.setZero();

			return res;
		}

		static Real W_zero()
		{
			return m_W_zero;
		}
	};

	template<typename KernelType, unsigned int resolution>
	Real PrecomputedKernel<KernelType, resolution>::m_radius;
	template<typename KernelType, unsigned int resolution>
	Real PrecomputedKernel<KernelType, resolution>::m_W[resolution];
	template<typename KernelType, unsigned int resolution>
	Real PrecomputedKernel<KernelType, resolution>::m_gradW[resolution + 1];
	template<typename KernelType, unsigned int resolution>
	Real PrecomputedKernel<KernelType, resolution>::m_invStepSize;
	template<typename KernelType, unsigned int resolution>
	Real PrecomputedKernel<KernelType, resolution>::m_W_zero;
}

#endif
