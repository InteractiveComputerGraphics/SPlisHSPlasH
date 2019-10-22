#include "SPlisHSPlasH/Common.h"

// Let Catch provide main():
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "SPlisHSPlasH/SPHKernels.h"

using namespace SPH;

TEMPLATE_TEST_CASE("3D Kernel is normalized, positive", "", CubicKernel, Poly6Kernel, SpikyKernel, WendlandQuinticC2Kernel, PrecomputedKernel<CubicKernel>)
{
	const Real supportRadius = 4.0*0.025;
	TestType::setRadius(supportRadius);
	const unsigned int numberOfSteps = 50;
	const Real stepSize = 2.0*supportRadius / (Real)(numberOfSteps - 1);
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
				if (W < -1.0e-5)
					positive = false;
			}
		}
	}
	REQUIRE(fabs(sum - 1.0) < 1.0e-5);
	REQUIRE(sumV.norm() < 1.0e-5);
	REQUIRE(positive);
}

TEMPLATE_TEST_CASE("2D Kernel is normalized, positive", "", CubicKernel2D, WendlandQuinticC2Kernel2D, PrecomputedKernel<CubicKernel2D>)
{
	const Real supportRadius = 4.0*0.025;
	TestType::setRadius(supportRadius);
	const unsigned int numberOfSteps = 50;
	const Real stepSize = 2.0*supportRadius / (Real)(numberOfSteps - 1);
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
			if (W < -1.0e-5)
				positive = false;
		}
	}
	REQUIRE(fabs(sum - 1.0) < 1.0e-5);
	REQUIRE(sumV.norm() < 1.0e-5);
	REQUIRE(positive);
}