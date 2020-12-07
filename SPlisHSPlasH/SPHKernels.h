#ifndef SPHKERNELS_H
#define SPHKERNELS_H

#define _USE_MATH_DEFINES
#include <math.h>
#include "Common.h"
#include <algorithm>
#ifdef USE_AVX
#include "SPlisHSPlasH/Utilities/AVX_math.h"
#endif

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
			const Real pi = static_cast<Real>(M_PI);

			const Real h3 = m_radius*m_radius*m_radius;
			m_k = static_cast<Real>(8.0) / (pi*h3);
			m_l = static_cast<Real>(48.0) / (pi*h3);
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
					res = m_k * (static_cast<Real>(6.0)*q3 - static_cast<Real>(6.0)*q2 + static_cast<Real>(1.0));
				}
				else
				{
					res = m_k * (static_cast<Real>(2.0)*pow(static_cast<Real>(1.0) - q, static_cast<Real>(3.0)));
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
			if ((rl > 1.0e-5) && (q <= 1.0))
			{
				const Vector3r gradq = r * (static_cast<Real>(1.0) / (rl*m_radius));
				if (q <= 0.5)
				{
					res = m_l*q*((Real) 3.0*q - static_cast<Real>(2.0))*gradq;
				}
				else
				{
					const Real factor = static_cast<Real>(1.0) - q;
					res = m_l*(-factor*factor)*gradq;
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
			const Real pi = static_cast<Real>(M_PI);
			m_k = static_cast<Real>(315.0) / (static_cast<Real>(64.0)*pi*pow(m_radius, static_cast<Real>(9.0)));
			m_l = -static_cast<Real>(945.0) / (static_cast<Real>(32.0)*pi*pow(m_radius, static_cast<Real>(9.0)));
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
				res = pow(radius2 - r2, static_cast<Real>(3.0))*m_k;
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
				res = pow(radius2 - r2, static_cast<Real>(3.0))*m_k;
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
			const Real radius6 = pow(m_radius, static_cast<Real>(6.0));
			const Real pi = static_cast<Real>(M_PI);
			m_k = static_cast<Real>(15.0) / (pi*radius6);
			m_l = -static_cast<Real>(45.0) / (pi*radius6);
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
				const Real hr3 = pow(m_radius - r, static_cast<Real>(3.0));
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
				const Real hr3 = pow(m_radius - sqrt(r2), static_cast<Real>(3.0));
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

	/** \brief quintic Wendland C2 kernel.
	*/
	class WendlandQuinticC2Kernel
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
			const Real pi = static_cast<Real>(M_PI);

			const Real h3 = m_radius*m_radius*m_radius;
			m_k = static_cast<Real>(21.0) / (static_cast<Real>(2.0)*pi*h3);
			m_l = -static_cast<Real>(210.0) / (pi*h3);
			m_W_zero = W(0.0);
		}

	public:
		static Real W(const Real r)
		{
			Real res = 0.0;
			const Real q = r / m_radius;
			if (q <= 1.0)
				res = m_k * pow(static_cast<Real>(1.0) - q, static_cast<Real>(4.0)) * (static_cast<Real>(4.0) * q + static_cast<Real>(1.0));
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
				const Vector3r gradq = r * (static_cast<Real>(1.0) / (rl*m_radius));
				res = m_l*gradq*pow(static_cast<Real>(1.0) - q, 3);
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

	/** \brief Cohesion kernel used for the surface tension method of Akinci el al. [ATT13].
	*
	* References:
	* - [AAT13] Nadir Akinci, Gizem Akinci, and Matthias Teschner. Versatile surface tension and adhesion for sph fluids. ACM Trans. Graph., 32(6):182:1-182:8, November 2013. URL: http://doi.acm.org/10.1145/2508363.2508395
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
			const Real pi = static_cast<Real>(M_PI);
			m_k = static_cast<Real>(32.0) / (pi*pow(m_radius, static_cast<Real>(9.0)));
			m_c = pow(m_radius, static_cast<Real>(6.0)) / static_cast<Real>(64.0);
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
					res = m_k*pow(m_radius - r1, static_cast<Real>(3.0))*r3;
				else
					res = m_k* static_cast<Real>(2.0)*pow(m_radius - r1, static_cast<Real>(3.0))*r3 - m_c;

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
					res = m_k*pow(m_radius - r1, static_cast<Real>(3.0))*r3;
				else
					res = m_k* static_cast<Real>(2.0)*pow(m_radius - r1, static_cast<Real>(3.0))*r3 - m_c;

			}
			return res;
		}

		static Real W_zero()
		{
			return m_W_zero;
		}
	};

	/** \brief Adhesion kernel used for the surface tension method of Akinci el al. [ATT13].
	*
	* References:
	* - [AAT13] Nadir Akinci, Gizem Akinci, and Matthias Teschner. Versatile surface tension and adhesion for sph fluids. ACM Trans. Graph., 32(6):182:1-182:8, November 2013. URL: http://doi.acm.org/10.1145/2508363.2508395
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
			m_k = static_cast<Real>(0.007) / pow(m_radius, static_cast<Real>(3.25));
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
				if (r > 0.5*m_radius)
					res = m_k*pow(-static_cast<Real>(4.0)*r2 / m_radius + static_cast<Real>(6.0)*r - static_cast<Real>(2.0)*m_radius, static_cast<Real>(0.25));
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
				const Real rl = sqrt(r2);
				if (rl > 0.5*m_radius)
					res = m_k*pow(-static_cast<Real>(4.0)*r2 / m_radius + static_cast<Real>(6.0)*rl - static_cast<Real>(2.0)*m_radius, static_cast<Real>(0.25));
			}
			return res;
		}

		static Real W_zero()
		{
			return m_W_zero;
		}
	};


	/** \brief Cubic spline kernel (2D).
	*/
	class CubicKernel2D
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
			const Real pi = static_cast<Real>(M_PI);

			const Real h2 = m_radius*m_radius;
			m_k = static_cast<Real>(40.0) / (static_cast<Real>(7.0) * (pi*h2));
			m_l = static_cast<Real>(240.0) / (static_cast<Real>(7.0) * (pi*h2));

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
					res = m_k * (static_cast<Real>(6.0)*q3 - static_cast<Real>(6.0)*q2 + static_cast<Real>(1.0));
				}
				else
				{
					res = m_k * (static_cast<Real>(2.0)*pow(static_cast<Real>(1.0) - q, static_cast<Real>(3.0)));
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
				if (rl > 1.0e-5)
				{
					const Vector3r gradq = r * (static_cast<Real>(1.0) / (rl*m_radius));
					if (q <= 0.5)
					{
						res = m_l*q*(static_cast<Real>(3.0)*q - static_cast<Real>(2.0))*gradq;
					}
					else
					{
						const Real factor = static_cast<Real>(1.0) - q;
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

	/** \brief Wendland Quintic C2 spline kernel (2D).
	*/
	class WendlandQuinticC2Kernel2D
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
			const Real pi = static_cast<Real>(M_PI);

			const Real h2 = m_radius*m_radius;
			m_k = static_cast<Real>(7.0) / (pi*h2);
			m_l = -static_cast<Real>(140.0) / (pi*h2);

			m_W_zero = W(Vector3r::Zero());
		}

	public:
		static Real W(const Real r)
		{
			Real res = 0.0;
			const Real q = r / m_radius;
			if (q <= 1.0)
				res = m_k * pow(static_cast<Real>(1.0) - q, static_cast<Real>(4.0)) * (static_cast<Real>(4.0) * q + static_cast<Real>(1.0));
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
				const Vector3r gradq = r * (static_cast<Real>(1.0) / (rl*m_radius));
				res = m_l*q*pow(static_cast<Real>(1.0) - q, static_cast<Real>(3.0))*gradq;
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


	/** \brief Precomputed kernel which is based on a lookup table as described by Bender and Koschier [BK15,BK17].
	*
	* The lookup tables can be used in combination with any kernel. 
	*
	* References:
	* - [BK15] Jan Bender and Dan Koschier. Divergence-free smoothed particle hydrodynamics. In ACM SIGGRAPH / Eurographics Symposium on Computer Animation, SCA '15, 147-155. New York, NY, USA, 2015. ACM. URL: http://doi.acm.org/10.1145/2786784.2786796
	* - [BK17] Jan Bender and Dan Koschier. Divergence-free SPH for incompressible and viscous fluids. IEEE Transactions on Visualization and Computer Graphics, 23(3):1193-1206, 2017. URL: http://dx.doi.org/10.1109/TVCG.2016.2578335
	*/
	template<typename KernelType, unsigned int resolution = 10000u>
	class PrecomputedKernel
	{
	protected:
		static Real m_W[resolution];
		static Real m_gradW[resolution + 1];
		static Real m_radius;
		static Real m_radius2;
		static Real m_invStepSize;
		static Real m_W_zero;
	public:
		static Real getRadius() { return m_radius; }
		static void setRadius(Real val)
		{
			m_radius = val;
			m_radius2 = m_radius * m_radius;
			KernelType::setRadius(val);
			const Real stepSize = m_radius / (Real)(resolution-1);
			m_invStepSize = static_cast<Real>(1.0) / stepSize;
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
			m_W_zero = W(static_cast<Real>(0));
		}

	public:
		static Real W(const Vector3r &r)
		{
			Real res = 0.0;
			const Real r2 = r.squaredNorm();
			if (r2 <= m_radius2)
			{
				const Real rl = sqrt(r2);
				const unsigned int pos = std::min<unsigned int>((unsigned int)(rl * m_invStepSize), resolution-2u);
				res = static_cast<Real>(0.5)*(m_W[pos]+ m_W[pos+1]);
			}
			return res;
		}

		static Real W(const Real r)
		{
			Real res = 0.0;
			if (r <= m_radius)
			{
				const unsigned int pos = std::min<unsigned int>((unsigned int)(r * m_invStepSize), resolution-2u);
				res = static_cast<Real>(0.5)*(m_W[pos] + m_W[pos + 1]);
			}
			return res;
		}

		static Vector3r gradW(const Vector3r &r)
		{
			Vector3r res;
			const Real rl = r.norm();
			if (rl <= m_radius)
			{
				//const Real rl = sqrt(r2);
				const unsigned int pos = std::min<unsigned int>(static_cast<unsigned int>(rl * m_invStepSize), resolution-2u);
				res = static_cast<Real>(0.5)*(m_gradW[pos] + m_gradW[pos + 1]) * r;
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

#ifdef USE_AVX
	/** \brief Cubic spline kernel.
	*/
	class CubicKernel_AVX
	{
	protected:
		static Real m_r;
		static Scalarf8 m_invRadius;
		static Scalarf8 m_invRadius2;
		static Scalarf8 m_k;
		static Scalarf8 m_l;
		static Real m_W_zero;
		static Scalarf8 m_zero;
		static Scalarf8 m_half;
		static Scalarf8 m_one;
		static Scalarf8 m_eps;


	public:
		static Real getRadius() { return m_r; }
		static void setRadius(Real val, bool is2D = false)
		{
			m_r = val;
			m_invRadius = Scalarf8(1.0f/val);
			m_invRadius2 = m_invRadius*m_invRadius;
			const Real pi = static_cast<Real>(M_PI);

			if (!is2D)
			{
				const Real h3 = m_r * m_r * m_r;
				m_k = Scalarf8(8.0f / static_cast<float>(pi*h3));
				m_l = Scalarf8(48.0f / static_cast<float>(pi*h3));
			}
			else
			{
				const Real h2 = m_r * m_r;
				m_k = static_cast<Real>(40.0) / (static_cast<Real>(7.0) * (pi*h2));
				m_l = static_cast<Real>(240.0) / (static_cast<Real>(7.0) * (pi*h2));
			}
			m_zero = Scalarf8(0.0f);
			m_half = Scalarf8(0.5f);
			m_one = Scalarf8(1.0f);
			m_eps = Scalarf8(1.0e-5f);
			Scalarf8 W_zero = W(m_zero);
			float tmp[8];
			W_zero.store(tmp);
			m_W_zero = tmp[0];
		}

	public:
		static Scalarf8 W(const Scalarf8 r)
		{
			Scalarf8 res;
			const Scalarf8 q = r * m_invRadius;

			const Scalarf8 v = m_one - q;

			// q <= 0.5
			const Scalarf8 res1 = m_k * (Scalarf8(-6.0f) * q*q * v + m_one);
			// 0.5 <= q <= 1
			const Scalarf8 res2 = (m_k * Scalarf8(2.0f)*(v*v*v));

			res = blend(q <= m_one, res2, m_zero);
			res = blend(q <= m_half, res1, res);
			
			return res;
		}

 		static Scalarf8 W(const Vector3f8 &r)
 		{
 			return W(r.norm());
 		}
 
 		static Vector3f8 gradW(const Vector3f8 &r)
 		{
			Scalarf8 res;
 			const Scalarf8 rl = r.norm();
 			const Scalarf8 q = rl * m_invRadius;

			// q <= 0.5
			const Scalarf8 res1 =  (m_l * m_invRadius2 * (Scalarf8(3.0f)*q - Scalarf8(2.0f)));

			// 0.5 <= q <= 1
			const Scalarf8 v = m_one - q;
			const Scalarf8 gradq = (m_invRadius / rl);
			const Scalarf8 res2 = gradq * (-m_l * (v*v));


			res = blend(q <= m_one, res2, m_zero);
			res = blend(q <= m_half, res1, res);
			res = blend(rl > m_eps, res, m_zero);

 			return r * res;
 		}

		static const Real& W_zero() 
		{
			return m_W_zero;
		}
	};

	/** \brief Poly6 kernel.
*/
	class Poly6Kernel_AVX
	{
	protected:
		static Real m_radius;
		static Real m_k;
		static Real m_l;
		static Scalarf8 m_radius_avx;
		static Real m_W_zero;
	public:
		static Real getRadius() { return m_radius; }
		static void setRadius(Real val)
		{
			m_radius = val;
			const Real pi = static_cast<Real>(M_PI);
 			m_k = static_cast<Real>(315.0) / (static_cast<Real>(64.0)*pi*pow(m_radius, static_cast<Real>(9)));
 			m_l = -static_cast<Real>(945.0) / (static_cast<Real>(32.0)*pi*pow(m_radius, static_cast<Real>(9)));
			m_radius_avx = Scalarf8(m_radius);
			Scalarf8 W_zero = W(Scalarf8(0.0f));
			float tmp[8];
			W_zero.store(tmp);
			m_W_zero = tmp[0];
		}

	public:

		/**
		* W(r,h) = (315/(64 pi h^9))(h^2-|r|^2)^3
		*        = (315/(64 pi h^9))(h^2-r*r)^3
		*/
		static Scalarf8 W(const Scalarf8 r)
		{
			Scalarf8 res;
			const Scalarf8 r2 = r * r;
			const Scalarf8 radius2 = m_radius_avx * m_radius_avx;
			const Scalarf8 t = (radius2 - r2);
			res = t*t*t*Scalarf8(m_k);
			return blend(r2 <= radius2, res, Scalarf8(0.0f));
		}

		static Scalarf8 W(const Vector3f8 &r)
		{
			Scalarf8 res;
			const Scalarf8 r2 = r.squaredNorm();
			const Scalarf8 radius2 = m_radius_avx * m_radius_avx;
			const Scalarf8 t = (radius2 - r2);
			res = t * t * t * Scalarf8(m_k);
			return blend(r2 <= radius2, res, Scalarf8(0.0f));
		}


		/**
		* grad(W(r,h)) = r(-945/(32 pi h^9))(h^2-|r|^2)^2
		*              = r(-945/(32 pi h^9))(h^2-r*r)^2
		*/
		static Vector3f8 gradW(const Vector3f8 &r)
		{
			Vector3f8 res;
			res.setZero();
			const Scalarf8 r2 = r.squaredNorm();
			const Scalarf8 radius2 = m_radius_avx * m_radius_avx;
			const Scalarf8 t = (radius2 - r2);

			const Vector3f8 res2 = r * (t*t*Scalarf8(m_l));
			return Vector3f8::blend(r2 <= radius2, res2, res);
		}

		//static Scalarf8 W_zero()
		//{
		//	return m_W_zero;
		//}

		static const Real& W_zero()
		{
			return m_W_zero;
		}
	};


	/** \brief Spiky kernel.
	*/
	class SpikyKernel_AVX
	{
	protected:
		static Real m_radius;
		static Real m_k;
		static Real m_l;
		static Scalarf8 m_radius_avx;
		static Scalarf8 m_W_zero;
		static Scalarf8 m_eps;
	public:
		static Real getRadius() { return m_radius; }
		static void setRadius(Real val)
		{
			m_radius = val;
			const Real radius6 = pow(m_radius, 6.0f);
 			const Real pi = static_cast<Real>(M_PI);
 			m_k = static_cast<Real>(15.0) / (pi*radius6);
 			m_l = -static_cast<Real>(45.0) / (pi*radius6);
			m_radius_avx = Scalarf8(m_radius);
 			m_W_zero = W(Scalarf8(0.0f));
 			m_eps = Scalarf8(1.0e-5f);
		}

	public:

		/**
		* W(r,h) = 15/(pi*h^6) * (h-r)^3
		*/
		static Scalarf8 W(const Scalarf8 r)
		{
			Scalarf8 res;
			const Scalarf8 t = (m_radius_avx - r);
			res = t*t*t*Scalarf8(m_k);
			return blend(r <= m_radius_avx, res, Scalarf8(0.0f));
		}

		static Scalarf8 W(const Vector3f8 &r)
		{
			return W(r.norm());
		}


		/**
		* grad(W(r,h)) = -r(45/(pi*h^6) * (h-r)^2)
		*/
		static Vector3f8 gradW(const Vector3f8 &r)
		{
			Vector3f8 res;
			res.setZero();
			const Scalarf8 r2 = r.squaredNorm();
			const Scalarf8 radius2 = m_radius_avx * m_radius_avx;
			const Scalarf8 r_l = r2.sqrt();
			const Scalarf8 t = (m_radius_avx - r_l);

			Vector3f8 res2 = r * (t*t*Scalarf8(m_l)/r_l);

			res2 = Vector3f8::blend(r_l > m_eps, res2, res);
			res = Vector3f8::blend(r2 <= radius2, res2, res);
			return res;
		}

		static Scalarf8 W_zero()
		{
			return m_W_zero;
		}
	};

#endif

	template<typename KernelType, unsigned int resolution>
	Real PrecomputedKernel<KernelType, resolution>::m_radius;
	template<typename KernelType, unsigned int resolution>
	Real PrecomputedKernel<KernelType, resolution>::m_radius2;
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
