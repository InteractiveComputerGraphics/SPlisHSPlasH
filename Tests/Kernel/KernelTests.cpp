#include "SPlisHSPlasH/Common.h"

// Let Catch provide main():
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "SPlisHSPlasH/SPHKernels.h"

using namespace SPH;

#ifdef USE_DOUBLE
	double eps = 1.0e-5;
#else
	float eps = 1.0e-4f;
#endif 

TEMPLATE_TEST_CASE("3D Kernel is normalized, positive", "", CubicKernel, Poly6Kernel, SpikyKernel, WendlandQuinticC2Kernel, PrecomputedKernel<CubicKernel>)
{
	const Real supportRadius = static_cast<Real>(4.0*0.025);
	TestType::setRadius(supportRadius);
	const unsigned int numberOfSteps = 50;
	const Real stepSize = static_cast<Real>(2.0)*supportRadius / (Real)(numberOfSteps - 1);
	Vector3r xi;
	xi.setZero();
	Real sum = 0.0;
	Vector3r sumV = Vector3r::Zero();
	bool positive = true;
	Real V = pow(stepSize, 3);
	for (unsigned int i = 0; i < numberOfSteps; i++)
	{
		for (unsigned int j = 0; j < numberOfSteps; j++)
		{
			for (unsigned int k = 0; k < numberOfSteps; k++)
			{
				const Vector3r xj(-supportRadius + i*stepSize, -supportRadius + j*stepSize, -supportRadius + k*stepSize);
				const Real W = TestType::W(xi - xj);
				sum += W *V;
				sumV += TestType::gradW(xi - xj) * V;
				if (W < -eps)
					positive = false;
			}
		}
	}
	REQUIRE(fabs(sum - 1.0) < eps);
	REQUIRE(sumV.norm() < eps);
	REQUIRE(positive);
}

TEMPLATE_TEST_CASE("2D Kernel is normalized, positive", "", CubicKernel2D, WendlandQuinticC2Kernel2D, PrecomputedKernel<CubicKernel2D>)
{
	const Real supportRadius = static_cast<Real>(4.0*0.025);
	TestType::setRadius(supportRadius);
	const unsigned int numberOfSteps = 50;
	const Real stepSize = static_cast<Real>(2.0)*supportRadius / (Real)(numberOfSteps - 1);
	Vector3r xi;
	xi.setZero();
	Real sum = 0.0;
	Vector3r sumV = Vector3r::Zero();
	bool positive = true;
	Real V = pow(stepSize, 2);
	for (unsigned int i = 0; i < numberOfSteps; i++)
	{
		for (unsigned int j = 0; j < numberOfSteps; j++)
		{
			const Vector3r xj(-supportRadius + i*stepSize, -supportRadius + j*stepSize, 0.0);
			const Real W = TestType::W(xi - xj);
			sum += W * V;
			sumV += TestType::gradW(xi - xj) * V;
			if (W < -eps)
				positive = false;
		}
	}
	REQUIRE(fabs(sum - 1.0) < eps);
	REQUIRE(sumV.norm() < eps);
	REQUIRE(positive);
}

#ifdef USE_AVX

TEMPLATE_TEST_CASE("AVX Kernel is normalized", "", CubicKernel_AVX, Poly6Kernel_AVX, SpikyKernel_AVX)
{
	const Real supportRadius = static_cast<Real>(4.0*0.025);
	TestType::setRadius(supportRadius);
	const unsigned int numberOfSteps = 50;
	const Real stepSize = static_cast<Real>(2.0)*supportRadius / (Real)(numberOfSteps - 1);
	Vector3f8 xi_avx;
	xi_avx.setZero();
	Scalarf8 sum(0.0);
	Vector3f8 sumV;
	sumV.setZero();

 	for (unsigned int i = 0; i < numberOfSteps; i++)
 	{
 		for (unsigned int j = 0; j < numberOfSteps; j++)
 		{
 			for (unsigned int k = 0; k < numberOfSteps; k += 8)
 			{
 				const unsigned int count = std::min(numberOfSteps - k, 8u);
 				const Scalarf8 V = convert_zero(stepSize*stepSize*stepSize, count);
 				Vector3r xj[8];
 				for (unsigned int l = 0; l < 8; l++)
 				{
 					if (l < count)
 						xj[l] = Vector3r(-supportRadius + i * stepSize, -supportRadius + j * stepSize, -supportRadius + (k+l) * stepSize);
 					else
 						xj[l].setZero();
  				}
 				Vector3f8 xj_avx(xj);
 				const Scalarf8 W = TestType::W(xi_avx - xj_avx);
 				sum += W * V;
				sumV += TestType::gradW(xi_avx - xj_avx) * V;
 			}
 		}
 	}
	Real res = sum.reduce();
	Vector3r res2;
	res2[0] = sumV.x().reduce();
	res2[1] = sumV.y().reduce();
	res2[2] = sumV.z().reduce();
	REQUIRE(fabs(res - 1.0) < eps);
	REQUIRE(res2.norm() < eps);
}

TEST_CASE("AVX Kernel is compared with corresponding kernel", "[CubicKernel]")
{
	const Real supportRadius = static_cast<Real>(4.0*0.025);
	CubicKernel::setRadius(supportRadius);
	CubicKernel_AVX::setRadius(supportRadius);
	const unsigned int numberOfSteps = 50;
	const Real stepSize = static_cast<Real>(2.0)*supportRadius / (Real)(numberOfSteps - 1);

	bool chk = true;
	bool chk2 = true;
	for (unsigned int i = 0; i < numberOfSteps; i++)
	{
		for (unsigned int j = 0; j < numberOfSteps; j++)
		{
			for (unsigned int k = 0; k < numberOfSteps; k += 8)
			{
				const unsigned int count = std::min(numberOfSteps - k, 8u);
				Vector3r xj[8];
				for (unsigned int l = 0; l < 8; l++)
				{
					if (l < count)
						xj[l] = Vector3r(-supportRadius + i * stepSize, -supportRadius + j * stepSize, -supportRadius + (k + l) * stepSize);
					else
						xj[l].setZero();
				}
				Vector3f8 xj_avx(xj);
				const Scalarf8 W = CubicKernel_AVX::W(xj_avx);

				// Kernel
				Real W_values[8];
				W.store(W_values);
				for (unsigned int l = 0; l < count; l++)
				{
					if (fabs(W_values[l] - CubicKernel::W(xj[l])) > 1.0e-3)
					{
						chk = false;
						//std::cout << fabs(W_values[l] - CubicKernel::W(xj[l])) << "\n";
					}
				}

				// Gradient of kernel
				const Vector3f8 gradW = CubicKernel_AVX::gradW(xj_avx);
				Vector3r gradW_values[8];
				gradW.store(gradW_values);
				for (unsigned int l = 0; l < count; l++)
				{
					Vector3r gradW_ = CubicKernel::gradW(xj[l]);
					if ((gradW_values[l] - gradW_).norm() > 3.0e-2)
					{
						chk2 = false;
						//std::cout << (gradW_values[l] - gradW_).norm()  << "\n";
					}
				}
			}
		}
	}

	// small steps
	bool chk3 = true;
	bool chk4 = true;
	const Real width = static_cast<Real>(1e-5);
	const unsigned int numberOfSteps2 = 200;
	const Real stepSize2 = width / (Real)(numberOfSteps2 - 1);
	for (unsigned int k = 0; k < numberOfSteps2; k += 8)
	{
		const unsigned int count = std::min(numberOfSteps2 - k, 8u);
		Vector3r xj[8];
		for (unsigned int l = 0; l < 8; l++)
		{
			if (l < count)
				xj[l] = Vector3r(0.0, 0.0, -width + (k + l) * stepSize2);
			else
				xj[l].setZero();
		}
		Vector3f8 xj_avx(xj);
		const Scalarf8 W = CubicKernel_AVX::W(xj_avx);

		// Kernel
		Real W_values[8];
		W.store(W_values);
		for (unsigned int l = 0; l < count; l++)
		{
			if (fabs(W_values[l] - CubicKernel::W(xj[l])) > 1.0e-3)
			{
				chk3 = false;
				//std::cout << fabs(W_values[l] - CubicKernel::W(xj[l])) << "\n";
			}
		}

		// Gradient of kernel
		const Vector3f8 gradW = CubicKernel_AVX::gradW(xj_avx);
		Vector3r gradW_values[8];
		gradW.store(gradW_values);
		for (unsigned int l = 0; l < count; l++)
		{
			Vector3r gradW_ = CubicKernel::gradW(xj[l]);
			if ((gradW_values[l] - gradW_).norm() > 3.0e-2)
			{
				chk4 = false;
				//std::cout << (gradW_values[l] - gradW_).norm() << "\n";
			}
		}
	}
	REQUIRE(chk);
	REQUIRE(chk2);
	REQUIRE(chk3);
	REQUIRE(chk4);
}

#endif