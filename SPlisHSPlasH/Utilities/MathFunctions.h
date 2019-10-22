#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include "SPlisHSPlasH/Common.h"

// ------------------------------------------------------------------------------------
namespace SPH
{
	class MathFunctions
	{
	public:
		/** Implementation of the paper: \n
		 * Matthias Müller, Jan Bender, Nuttapong Chentanez and Miles Macklin, 
		 * "A Robust Method to Extract the Rotational Part of Deformations", 
		 * ACM SIGGRAPH Motion in Games, 2016
		 */
		static void extractRotation(const Matrix3r &A, Quaternionr &q, const unsigned int maxIter);

		static void pseudoInverse(const Matrix3r &a, Matrix3r &res);
		static void svdWithInversionHandling(const Matrix3r &A, Vector3r &sigma, Matrix3r &U, Matrix3r &VT);
		static void eigenDecomposition(const Matrix3r &A, Matrix3r &eigenVecs, Vector3r &eigenVals);
		static void jacobiRotate(Matrix3r &A, Matrix3r &R, int p, int q);

		/** Returns two orthogonal vectors to vec which are also orthogonal to each other.
		*/
		static void getOrthogonalVectors(const Vector3r &vec, Vector3r &x, Vector3r &y);
	};
}

#endif