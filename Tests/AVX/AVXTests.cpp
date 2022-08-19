#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/Utilities/AVX_math.h"
#define _USE_MATH_DEFINES
#include <math.h>

// Let Catch provide main():
#define CATCH_CONFIG_MAIN

#include "catch.hpp"

using namespace SPH;

#ifdef USE_DOUBLE
	double eps = 1.0e-5;
#else
	float eps = 1.0e-4f;
#endif 

TEST_CASE("AVX quaternion tests", "[Quaternion8f]")
{
	// set some rotation matrices and convert them to quaternions
	Matrix3r R[8];
	Quaternionr q[8];
	for (int i = 0; i < 8; i++)
	{
		Real angle = static_cast<Real>((0.1 * i) * M_PI);
		R[i] = AngleAxisr(angle, Vector3r(2, 1, 0.25).normalized());
		q[i] = Quaternionr(R[i]);
	}
	// convert to avx quaternion
	Quaternion8f q_avx;
	q_avx.set(q);

	// convert everything back and check equality
	Matrix3r R2[8];
	Quaternionr q2[8];
	for (int i = 0; i < 8; i++)
	{
		R2[i] = q[i].toRotationMatrix();
		REQUIRE(R[i].isApprox(R2[i], eps));
	}

	// convert back from avx quaternion
	q_avx.store(q2);
	Vector3f8 col1_avx, col2_avx, col3_avx;	//columns of the rotation matrix
	q_avx.toRotationMatrix(col1_avx, col2_avx, col3_avx);
	Vector3r col1[8], col2[8], col3[8];
	col1_avx.store(col1);
	col2_avx.store(col2);
	col3_avx.store(col3);
	for (int i = 0; i < 8; i++)
	{
		REQUIRE(q[i].isApprox(q2[i], eps));

		R2[i].col(0) = col1[i];
		R2[i].col(1) = col2[i];
		R2[i].col(2) = col3[i];
		REQUIRE(R[i].isApprox(R2[i], eps));
	}

	REQUIRE(true);
}

