#include "TubeMapDensity.h"
#include <algorithm>
#include <array>
#include <cmath>

using namespace SPH;

namespace
{
	/** Exact integral of P(t) * (E0 + E1 t + E2 t^2 + E3 t^3) over [t1, t2] for a polynomial P
	* of degree 5 given in ascending powers.
	*/
	double integratePolyTimesCubic(const std::array<double, 6> &P,
		const double E0, const double E1, const double E2, const double E3,
		const double t1, const double t2)
	{
		auto primitive = [&](const double t) -> double
		{
			double res = 0.0;
			// powers t^(i+1), ..., t^(i+4)
			double t1pow = t;
			double t2pow = t1pow * t;
			double t3pow = t2pow * t;
			double t4pow = t3pow * t;
			for (int i = 0; i < 6; ++i)
			{
				const double pi = P[i];
				if (pi != 0.0)
				{
					res += pi * (E0 * t1pow / (i + 1.0)
							   + E1 * t2pow / (i + 2.0)
							   + E2 * t3pow / (i + 3.0)
							   + E3 * t4pow / (i + 4.0));
				}
				t1pow = t2pow;
				t2pow = t3pow;
				t3pow = t4pow;
				t4pow = t4pow * t;
			}
			return res;
		};
		return primitive(t2) - primitive(t1);
	}

	/** Exact integral of B(t) * (g0 + g1 t) over [t1, t2] for a polynomial B of degree 7
	* given in ascending powers.
	*/
	double integratePolyTimesLinear(const std::array<double, 8> &B,
		const double g0, const double g1,
		const double t1, const double t2)
	{
		auto primitive = [&](const double t) -> double
		{
			double res = 0.0;
			// powers t^(i+1), t^(i+2)
			double t1pow = t;
			double t2pow = t1pow * t;
			for (int i = 0; i < 8; ++i)
			{
				const double bi = B[i];
				if (bi != 0.0)
				{
					res += bi * (g0 * t1pow / (i + 1.0)
							   + g1 * t2pow / (i + 2.0));
				}
				t1pow = t2pow;
				t2pow = t2pow * t;
			}
			return res;
		};
		return primitive(t2) - primitive(t1);
	}

	// Cubic spline kernel projected onto a plane as polynomial in q = |t|/h (ascending powers).
	// Inner: 0 <= q <= 0.5, outer: 0.5 < q <= 1
	const std::array<double, 6> P_SLICE_IN  = { 0.0875, 0.0, -0.5, 0.0,  1.5, -1.2 };
	const std::array<double, 6> P_SLICE_OUT = { 0.1,    0.0, -1.0, 2.0, -1.5,  0.4 };

	// Second moment of the projected kernel as polynomial in q = |t|/h (ascending powers).
	// Inner: 0 <= q <= 0.5, outer: 0.5 < q <= 1
	const std::array<double, 8> B_SLICE_IN  = { -31.0, 0.0, 196.0, 0.0, -560.0, 0.0, 1120.0, -768.0 };
	const std::array<double, 8> B_SLICE_OUT = { -32.0, 0.0, 224.0, 0.0, -1120.0, 1792.0, -1120.0, 256.0 };

	// Limits for the curvature correction (h: support radius):
	// - principal curvatures |k1|, |k2| <= MAX_PRINCIPAL_CURVATURE / h. Spheres and cylinders
	//   give non-negative densities up to this value.
	// - saddle-type Gaussian curvature K >= -MAX_SADDLE_CURVATURE / h^2. Stronger saddles
	//   (|k1| = |k2| > 0.5/h) lead to negative boundary densities.
	const double MAX_PRINCIPAL_CURVATURE = 2.0;
	const double MAX_SADDLE_CURVATURE = 0.25;

	/** Clamp the curvatures to the range in which the expansion of the tube map is valid.
	* H = -(k1 + k2) and K = k1 * k2 are the mean and Gaussian curvature in terms of the
	* principal curvatures (the bounds are symmetric, so the sign of H does not matter). Discrete curvatures of fine or noisy meshes often violate these
	* bounds or are even inconsistent (K > H^2/4, i.e., no real principal curvatures).
	* Values within the bounds are not modified.
	*/
	void clampCurvature(const double h, double &H, double &K)
	{
		const double kMax = MAX_PRINCIPAL_CURVATURE / h;

		// real principal curvatures require K <= H^2/4
		if (K > 0.25 * H * H)
			K = 0.25 * H * H;

		// limit the principal curvatures
		const double disc = std::sqrt(std::max(0.25 * H * H - K, 0.0));
		const double k1 = 0.5 * H + disc;
		const double k2 = 0.5 * H - disc;
		if ((std::abs(k1) > kMax) || (std::abs(k2) > kMax))
		{
			const double k1c = std::max(-kMax, std::min(kMax, k1));
			const double k2c = std::max(-kMax, std::min(kMax, k2));
			H = k1c + k2c;
			K = k1c * k2c;
		}

		// limit saddle-type curvature
		const double KMin = -MAX_SADDLE_CURVATURE / (h * h);
		if (K < KMin)
			K = KMin;
	}
}

double TubeMapDensity::boundaryDensity(const double d, const double H_, const double K_, const double h, const double factor, const unsigned int order)
{
	// integration bounds: kernel support and cut-off of Gamma
	const double s_gamma = h / factor;
	const double t_cutoff = s_gamma - d;
	const double t_start = -h;
	const double t_end = std::min(h, t_cutoff);
	if (t_end <= t_start)
		return 0.0;

	// order 0: planar boundary, order 1: mean curvature only
	double H = (order >= 1) ? H_ : 0.0;
	double K = (order >= 2) ? K_ : 0.0;
	if (order >= 1)
		clampCurvature(h, H, K);

	// Gamma(d+t) = g0 + g1*t
	const double g0 = 1.0 - (factor * d) / h;
	const double g1 = -factor / h;

	// A(s) = c0 + c1*s + c2*s^2 with s = d + t
	const double c0 = 1.0 + 0.5 * H * d;
	const double c1 = -0.5 * H;
	const double c2 = (order >= 2) ? 2.0 * K + 0.5 * H * H : 0.0;

	// A(d+t) = A0 + A1*t + A2*t^2
	const double A0 = c0 + c1 * d + c2 * d * d;
	const double A1 = c1 + 2.0 * c2 * d;
	const double A2 = c2;

	// Gamma(d+t) * A(d+t) = E0 + E1*t + E2*t^2 + E3*t^3
	const double E0 = g0 * A0;
	const double E1 = g0 * A1 + g1 * A0;
	const double E2 = g0 * A2 + g1 * A1;
	const double E3 = g1 * A2;

	// coefficient of the second moment term (second order only)
	const double bcoef = (order >= 2) ? (3.0 / 64.0) * H * H - (1.0 / 16.0) * K : 0.0;

	// normalization of the projected kernel and scaling of the second moment
	const double K_A = 16.0 / h;
	const double scale_B = h / 35.0;

	// piecewise integration between the knots of the cubic spline
	const double knots[5] = { -h, -0.5*h, 0.0, 0.5*h, h };
	double total = 0.0;
	for (int seg = 0; seg < 4; ++seg)
	{
		const double a = std::max(t_start, knots[seg]);
		const double b = std::min(t_end, knots[seg + 1]);
		if (b <= a)
			continue;

		const double mid = 0.5 * (a + b);
		const bool inner = std::abs(mid) < 0.5 * h;
		const std::array<double, 6> &P_ref = inner ? P_SLICE_IN : P_SLICE_OUT;
		const std::array<double, 8> &B_ref = inner ? B_SLICE_IN : B_SLICE_OUT;

		// q^k = |t|^k / h^k: odd powers change sign for t < 0
		const double sign = (mid < 0.0) ? -1.0 : 1.0;

		std::array<double, 6> P_t;
		double h_pow = 1.0;
		for (int k = 0; k < 6; ++k)
		{
			const double sgnk = (k % 2 != 0) ? sign : 1.0;
			P_t[k] = K_A * sgnk * P_ref[k] / h_pow;
			h_pow *= h;
		}

		std::array<double, 8> B_t;
		h_pow = 1.0;
		for (int k = 0; k < 8; ++k)
		{
			const double sgnk = (k % 2 != 0) ? sign : 1.0;
			B_t[k] = scale_B * sgnk * B_ref[k] / h_pow;
			h_pow *= h;
		}

		total += integratePolyTimesCubic(P_t, E0, E1, E2, E3, a, b);
		if (order >= 2)
			total += bcoef * integratePolyTimesLinear(B_t, g0, g1, a, b);
	}
	return total;
}
