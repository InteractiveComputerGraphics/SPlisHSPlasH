#include "MathFunctions.h"
#include <cfloat>

using namespace SPH;


// ----------------------------------------------------------------------------------------------
void MathFunctions::extractRotation(const Matrix3r &A, Quaternionr &q,	const unsigned int maxIter)
{
	for (unsigned int iter = 0; iter < maxIter; iter++)
	{
		Matrix3r R = q.matrix();
		Vector3r omega = (R.col(0).cross(A.col(0)) + R.col(1).cross(A.col(1)) + R.col(2).cross(A.col(2))) * 
			(1.0 / fabs(R.col(0).dot(A.col(0)) + R.col(1).dot(A.col(1)) + R.col(2).dot(A.col(2)) + 1.0e-9));
		Real w = omega.norm();
		if (w < 1.0e-9)
			break;
		q = Quaternionr(AngleAxisr(w, (1.0 / w) * omega)) *	q;
		q.normalize();
	}
}

void MathFunctions::pseudoInverse(const Matrix3r &a, Matrix3r &res)
{
	const Real epsilon = std::numeric_limits<Real>::epsilon();
	const Eigen::JacobiSVD<Matrix3r> svd(a, Eigen::ComputeFullU | Eigen::ComputeFullV);
	const Real tolerance = epsilon * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs()(0);
	res = svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}


/** Perform a singular value decomposition of matrix A: A = U * sigma * V^T.
* This function returns two proper rotation matrices U and V^T which do not
* contain a reflection. Reflections are corrected by the inversion handling
* proposed by Irving et al. 2004.
*/
void MathFunctions::svdWithInversionHandling(const Matrix3r &A, Vector3r &sigma, Matrix3r &U, Matrix3r &VT)
{

	Matrix3r AT_A, V;
	AT_A = A.transpose() * A;

	Vector3r S;

	// Eigen decomposition of A^T * A
	eigenDecomposition(AT_A, V, S);

	// Detect if V is a reflection .
	// Make a rotation out of it by multiplying one column with -1.
	const Real detV = V.determinant();
	if (detV < 0.0)
	{
		Real minLambda = REAL_MAX;
		unsigned char pos = 0;
		for (unsigned char l = 0; l < 3; l++)
		{
			if (S[l] < minLambda)
			{
				pos = l;
				minLambda = S[l];
			}
		}
		V(0, pos) = -V(0, pos);
		V(1, pos) = -V(1, pos);
		V(2, pos) = -V(2, pos);
	}

	if (S[0] < 0.0) S[0] = 0.0;		// safety for sqrt
	if (S[1] < 0.0) S[1] = 0.0;
	if (S[2] < 0.0) S[2] = 0.0;

	sigma[0] = sqrt(S[0]);
	sigma[1] = sqrt(S[1]);
	sigma[2] = sqrt(S[2]);

	VT = V.transpose();

	//
	// Check for values of hatF near zero
	//
	unsigned char chk = 0;
	unsigned char pos = 0;
	for (unsigned char l = 0; l < 3; l++)
	{
		if (fabs(sigma[l]) < 1.0e-4)
		{
			pos = l;
			chk++;
		}
	}

	if (chk > 0)
	{
		if (chk > 1)
		{
			U.setIdentity();
		}
		else
		{
			U = A * V;
			for (unsigned char l = 0; l < 3; l++)
			{
				if (l != pos)
				{
					for (unsigned char m = 0; m < 3; m++)
					{
						U(m, l) *= static_cast<Real>(1.0) / sigma[l];
					}
				}
			}

			Vector3r v[2];
			unsigned char index = 0;
			for (unsigned char l = 0; l < 3; l++)
			{
				if (l != pos)
				{
					v[index++] = Vector3r(U(0, l), U(1, l), U(2, l));
				}
			}
			Vector3r vec = v[0].cross(v[1]);
			vec.normalize();
			U(0, pos) = vec[0];
			U(1, pos) = vec[1];
			U(2, pos) = vec[2];
		}
	}
	else
	{
		Vector3r sigmaInv(static_cast<Real>(1.0) / sigma[0], static_cast<Real>(1.0) / sigma[1], static_cast<Real>(1.0) / sigma[2]);
		U = A * V;
		for (unsigned char l = 0; l < 3; l++)
		{
			for (unsigned char m = 0; m < 3; m++)
			{
				U(m, l) *= sigmaInv[l];
			}
		}
	}

	const Real detU = U.determinant();

	// U is a reflection => inversion
	if (detU < 0.0)
	{
		//std::cout << "Inversion!\n";
		Real minLambda = REAL_MAX;
		unsigned char pos = 0;
		for (unsigned char l = 0; l < 3; l++)
		{
			if (sigma[l] < minLambda)
			{
				pos = l;
				minLambda = sigma[l];
			}
		}

		// invert values of smallest singular value
		sigma[pos] = -sigma[pos];
		U(0, pos) = -U(0, pos);
		U(1, pos) = -U(1, pos);
		U(2, pos) = -U(2, pos);
	}
}

// ----------------------------------------------------------------------------------------------
void MathFunctions::eigenDecomposition(const Matrix3r &A, Matrix3r &eigenVecs, Vector3r &eigenVals)
{
	const int numJacobiIterations = 10;
	const Real epsilon = static_cast<Real>(1e-15);

	Matrix3r D = A;

	// only for symmetric matrices!
	eigenVecs.setIdentity();	// unit matrix
	int iter = 0;
	while (iter < numJacobiIterations) {	// 3 off diagonal elements
											// find off diagonal element with maximum modulus
		int p, q;
		Real a, max;
		max = fabs(D(0, 1));
		p = 0; q = 1;
		a = fabs(D(0, 2));
		if (a > max) { p = 0; q = 2; max = a; }
		a = fabs(D(1, 2));
		if (a > max) { p = 1; q = 2; max = a; }
		// all small enough -> done
		if (max < epsilon) break;
		// rotate matrix with respect to that element
		jacobiRotate(D, eigenVecs, p, q);
		iter++;
	}
	eigenVals[0] = D(0, 0);
	eigenVals[1] = D(1, 1);
	eigenVals[2] = D(2, 2);
}

// ----------------------------------------------------------------------------------------------
void MathFunctions::jacobiRotate(Matrix3r &A, Matrix3r &R, int p, int q)
{
	// rotates A through phi in pq-plane to set A(p,q) = 0
	// rotation stored in R whose columns are eigenvectors of A
	if (A(p, q) == 0.0)
		return;

	Real d = (A(p, p) - A(q, q)) / (static_cast<Real>(2.0)*A(p, q));
	Real t = static_cast<Real>(1.0) / (fabs(d) + sqrt(d*d + static_cast<Real>(1.0)));
	if (d < 0.0) t = -t;
	Real c = static_cast<Real>(1.0) / sqrt(t*t + 1);
	Real s = t*c;
	A(p, p) += t*A(p, q);
	A(q, q) -= t*A(p, q);
	A(p, q) = A(q, p) = 0.0;
	// transform A
	int k;
	for (k = 0; k < 3; k++) {
		if (k != p && k != q) {
			Real Akp = c*A(k, p) + s*A(k, q);
			Real Akq = -s*A(k, p) + c*A(k, q);
			A(k, p) = A(p, k) = Akp;
			A(k, q) = A(q, k) = Akq;
		}
	}
	// store rotation in R
	for (k = 0; k < 3; k++) {
		Real Rkp = c*R(k, p) + s*R(k, q);
		Real Rkq = -s*R(k, p) + c*R(k, q);
		R(k, p) = Rkp;
		R(k, q) = Rkq;
	}
}

// ----------------------------------------------------------------------------------------------
void MathFunctions::getOrthogonalVectors(const Vector3r &vec, Vector3r &x, Vector3r &y)
{
	// Get plane vectors x, y
	Vector3r v(1, 0, 0);

	// Check, if v has same direction as vec
	if (fabs(v.dot(vec)) > 0.999)
		v = Vector3r(0, 1, 0);

	x = vec.cross(v);
	y = vec.cross(x);
	x.normalize();
	y.normalize();
}

