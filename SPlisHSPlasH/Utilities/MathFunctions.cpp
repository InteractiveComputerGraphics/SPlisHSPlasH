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

// ----------------------------------------------------------------------------------------------
void MathFunctions::APD_Newton(const Matrix3r& F, Quaternionr& q)
{
	//one iteration is sufficient for plausible results
	for (int it = 0; it < 1; it++)
	{
		//transform quaternion to rotation matrix
		Matrix3r R;
		R = q.matrix();

		//columns of B = RT * F
		Vector3r B0 = R.transpose() * F.col(0);
		Vector3r B1 = R.transpose() * F.col(1);
		Vector3r B2 = R.transpose() * F.col(2);

		Vector3r gradient(B2[1] - B1[2], B0[2] - B2[0], B1[0] - B0[1]);

		//compute Hessian, use the fact that it is symmetric
		const Real h00 = B1[1] + B2[2];
		const Real h11 = B0[0] + B2[2];
		const Real h22 = B0[0] + B1[1];
		const Real h01 = static_cast<Real>(-0.5) * (B1[0] + B0[1]);
		const Real h02 = static_cast<Real>(-0.5) * (B2[0] + B0[2]);
		const Real h12 = static_cast<Real>(-0.5) * (B2[1] + B1[2]);

		const Real detH = static_cast<Real>(-1.0) * h02 * h02 * h11 + static_cast<Real>(2.0) * h01 * h02 * h12 - h00 * h12 * h12 - h01 * h01 * h22 + h00 * h11 * h22;

		Vector3r omega;
		//compute symmetric inverse
		const Real factor = static_cast<Real>(-0.25) / detH;
		omega[0] = (h11 * h22 - h12 * h12) * gradient[0]
			+ (h02 * h12 - h01 * h22) * gradient[1]
			+ (h01 * h12 - h02 * h11) * gradient[2];
		omega[0] *= factor;

		omega[1] = (h02 * h12 - h01 * h22) * gradient[0]
			+ (h00 * h22 - h02 * h02) * gradient[1]
			+ (h01 * h02 - h00 * h12) * gradient[2];
		omega[1] *= factor;

		omega[2] = (h01 * h12 - h02 * h11) * gradient[0]
			+ (h01 * h02 - h00 * h12) * gradient[1]
			+ (h00 * h11 - h01 * h01) * gradient[2];
		omega[2] *= factor;

		//if det(H) = 0 use gradient descent, never happened in our tests, could also be removed 
		if (fabs(detH) < static_cast<Real>(1.0e-9))
			omega = -gradient;

		//instead of clamping just use gradient descent. also works fine and does not require the norm
		if (omega.dot(gradient) > 0.0)
			omega = gradient * static_cast<Real>(-0.125);

		const Real l_omega2 = omega.squaredNorm();
		const Real w = (static_cast<Real>(1.0) - l_omega2) / (static_cast<Real>(1.0) + l_omega2);
		const Vector3r vec = omega * (static_cast<Real>(2.0) / (static_cast<Real>(1.0) + l_omega2));
		q = q * Quaternionr(w, vec.x(), vec.y(), vec.z());		//no normalization needed because the Cayley map returs a unit quaternion
	}
}

// ----------------------------------------------------------------------------------------------
// iARAP helper: quartic solver (Lin et al. 2022)
// ----------------------------------------------------------------------------------------------
namespace {

const double M_2PI = 6.28318530717958647692;
const double quartic_eps = 1e-14;

// Solve cubic equation x^3 + a*x^2 + b*x + c = 0
unsigned int solveP3(double* x, double a, double b, double c)
{
	double a2 = a * a;
	double q = (a2 - 3 * b) / 9;
	double r = (a * (2 * a2 - 9 * b) + 27 * c) / 54;
	double r2 = r * r;
	double q3 = q * q * q;
	double A, B;
	if (r2 < q3)
	{
		double t = r / std::sqrt(q3);
		if (t < -1) t = -1;
		if (t > 1) t = 1;
		t = std::acos(t);
		a /= 3;
		q = -2 * std::sqrt(q);
		x[0] = q * std::cos(t / 3) - a;
		x[1] = q * std::cos((t + M_2PI) / 3) - a;
		x[2] = q * std::cos((t - M_2PI) / 3) - a;
		return 3;
	}
	else
	{
		A = -std::pow(std::fabs(r) + std::sqrt(r2 - q3), 1.0 / 3.0);
		if (r < 0) A = -A;
		B = (A == 0 ? 0 : q / A);
		a /= 3;
		x[0] = (A + B) - a;
		x[1] = -0.5 * (A + B) - a;
		x[2] = 0.5 * std::sqrt(3.0) * (A - B);
		if (std::fabs(x[2]) < quartic_eps) { x[2] = x[1]; return 2; }
		return 1;
	}
}

// Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d = 0
void solveQuartic(double a, double b, double c, double d, double roots[4])
{
	double a3 = -b;
	double b3 = a * c - 4.0 * d;
	double c3 = -a * a * d - c * c + 4.0 * b * d;

	double x3[3];
	unsigned int iZeroes = solveP3(x3, a3, b3, c3);

	double y = x3[0];
	if (iZeroes != 1)
	{
		if (std::fabs(x3[1]) > std::fabs(y)) y = x3[1];
		if (std::fabs(x3[2]) > std::fabs(y)) y = x3[2];
	}

	double q1, q2, p1, p2, D, sqD;

	D = y * y - 4 * d;
	if (std::fabs(D) < quartic_eps)
	{
		q1 = q2 = y * 0.5;
		D = a * a - 4 * (b - y);
		if (std::fabs(D) < quartic_eps)
			p1 = p2 = a * 0.5;
		else
		{
			sqD = std::sqrt(D);
			p1 = (a + sqD) * 0.5;
			p2 = (a - sqD) * 0.5;
		}
	}
	else
	{
		sqD = std::sqrt(D);
		q1 = (y + sqD) * 0.5;
		q2 = (y - sqD) * 0.5;
		p1 = (a * q1 - c) / (q1 - q2);
		p2 = (c - a * q2) / (q1 - q2);
	}

	D = p1 * p1 - 4 * q1;
	if (D < 0.0)
	{
		roots[0] = -p1 * 0.5;
		roots[1] = -p1 * 0.5;
	}
	else
	{
		sqD = std::sqrt(D);
		roots[0] = (-p1 + sqD) * 0.5;
		roots[1] = (-p1 - sqD) * 0.5;
	}

	D = p2 * p2 - 4 * q2;
	if (D < 0.0)
	{
		roots[2] = -p2 * 0.5;
		roots[3] = -p2 * 0.5;
	}
	else
	{
		sqD = std::sqrt(D);
		roots[2] = (-p2 + sqD) * 0.5;
		roots[3] = (-p2 - sqD) * 0.5;
	}
}

} // anonymous namespace

// ----------------------------------------------------------------------------------------------
void MathFunctions::iARAP(const Matrix3r& F, Matrix3r& R)
{
	// Compute Cauchy-Green invariants
	const Real I1 = F.squaredNorm();  // ||F||²
	const Matrix3r FtF = F.transpose() * F;
	const Real I2 = FtF.squaredNorm();  // ||F^T F||²
	const Real J = F.determinant();

	// Quartic polynomial: t⁴ - 2·I₁·t² - 8·J·t + (2·I₂ - I₁²)
	double roots[4];
	solveQuartic(0.0, -2.0 * I1, -8.0 * J, I1 * I1 - 2.0 * (I1 * I1 - I2), roots);

	// Find trace term f (largest root = σ₁ + σ₂ + σ₃)
	double f = roots[0];
	for (int k = 1; k < 4; k++)
		if (roots[k] > f) f = roots[k];

	// Compute R = df1*g1 + df2*g2 + dfJ*gJ
	const Real denom = 4.0 * f * f * f - 4.0 * I1 * f - 8.0 * J;

	if (std::fabs(denom) < static_cast<Real>(1e-12))
	{
		// Degenerate case: fall back to SVD polar decomposition
		Eigen::JacobiSVD<Matrix3r> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		R = svd.matrixU() * svd.matrixV().transpose();
		if (R.determinant() < 0)
			R.col(2) = -R.col(2);
		return;
	}

	const Real df1 = (2.0 * f * f + 2.0 * I1) / denom;
	const Real df2 = -2.0 / denom;
	const Real dfJ = (8.0 * f) / denom;

	// g1 = 2F, g2 = 4F·F^T·F, gJ = cof(F)
	Matrix3r g1 = 2.0 * F;
	Matrix3r g2 = 4.0 * F * FtF;
	Matrix3r gJ;
	gJ.col(0) = F.col(1).cross(F.col(2));
	gJ.col(1) = F.col(2).cross(F.col(0));
	gJ.col(2) = F.col(0).cross(F.col(1));

	R = df1 * g1 + df2 * g2 + dfJ * gJ;

	// Validate: if R is not a proper rotation, fall back to SVD
	const Real detR = R.determinant();
	const Real orthoErr = (R.transpose() * R - Matrix3r::Identity()).squaredNorm();
	if (std::fabs(detR - 1.0) > static_cast<Real>(1e-4) || orthoErr > static_cast<Real>(1e-4))
	{
		Eigen::JacobiSVD<Matrix3r> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		R = svd.matrixU() * svd.matrixV().transpose();
		if (R.determinant() < 0)
			R.col(2) = -R.col(2);
	}
}

// ----------------------------------------------------------------------------------------------
// Analytical eigenvector computation for 3x3 symmetric matrix (from iARAP reference)
// ----------------------------------------------------------------------------------------------
namespace {

void computeEigenvector0(Real a00, Real a01, Real a02, Real a11, Real a12, Real a22,
                         Real eval0, Vector3r& evec0)
{
	Vector3r row0(a00 - eval0, a01, a02);
	Vector3r row1(a01, a11 - eval0, a12);
	Vector3r row2(a02, a12, a22 - eval0);
	Vector3r r0xr1 = row0.cross(row1);
	Vector3r r0xr2 = row0.cross(row2);
	Vector3r r1xr2 = row1.cross(row2);
	Real d0 = r0xr1.squaredNorm();
	Real d1 = r0xr2.squaredNorm();
	Real d2 = r1xr2.squaredNorm();

	Real dmax = d0;
	int imax = 0;
	if (d1 > dmax) { dmax = d1; imax = 1; }
	if (d2 > dmax) { imax = 2; }

	if (imax == 0)
		evec0 = r0xr1 / std::sqrt(d0);
	else if (imax == 1)
		evec0 = r0xr2 / std::sqrt(d1);
	else
		evec0 = r1xr2 / std::sqrt(d2);
}

void computeOrthogonalComplement(const Vector3r& W, Vector3r& U, Vector3r& V)
{
	Real invLength;
	if (std::fabs(W[0]) > std::fabs(W[1]))
	{
		invLength = static_cast<Real>(1.0) / std::sqrt(W[0] * W[0] + W[2] * W[2]);
		U = Vector3r(-W[2] * invLength, static_cast<Real>(0.0), W[0] * invLength);
	}
	else
	{
		invLength = static_cast<Real>(1.0) / std::sqrt(W[1] * W[1] + W[2] * W[2]);
		U = Vector3r(static_cast<Real>(0.0), W[2] * invLength, -W[1] * invLength);
	}
	V = W.cross(U);
}

void computeEigenvector1(Real a00, Real a01, Real a02, Real a11, Real a12, Real a22,
                         const Vector3r& evec0, Real eval1, Vector3r& evec1)
{
	Vector3r U, V;
	computeOrthogonalComplement(evec0, U, V);

	Vector3r AU(a00 * U[0] + a01 * U[1] + a02 * U[2],
	            a01 * U[0] + a11 * U[1] + a12 * U[2],
	            a02 * U[0] + a12 * U[1] + a22 * U[2]);

	Vector3r AV(a00 * V[0] + a01 * V[1] + a02 * V[2],
	            a01 * V[0] + a11 * V[1] + a12 * V[2],
	            a02 * V[0] + a12 * V[1] + a22 * V[2]);

	Real m00 = U.dot(AU) - eval1;
	Real m01 = U.dot(AV);
	Real m11 = V.dot(AV) - eval1;

	Real absM00 = std::fabs(m00);
	Real absM01 = std::fabs(m01);
	Real absM11 = std::fabs(m11);

	if (absM00 >= absM11)
	{
		Real maxAbsComp = std::max(absM00, absM01);
		if (maxAbsComp > static_cast<Real>(0.0))
		{
			if (absM00 >= absM01)
			{
				m01 /= m00;
				m00 = static_cast<Real>(1.0) / std::sqrt(static_cast<Real>(1.0) + m01 * m01);
				m01 *= m00;
			}
			else
			{
				m00 /= m01;
				m01 = static_cast<Real>(1.0) / std::sqrt(static_cast<Real>(1.0) + m00 * m00);
				m00 *= m01;
			}
			evec1 = m01 * U - m00 * V;
		}
		else
			evec1 = U;
	}
	else
	{
		Real maxAbsComp = std::max(absM11, absM01);
		if (maxAbsComp > static_cast<Real>(0.0))
		{
			if (absM11 >= absM01)
			{
				m01 /= m11;
				m11 = static_cast<Real>(1.0) / std::sqrt(static_cast<Real>(1.0) + m01 * m01);
				m01 *= m11;
			}
			else
			{
				m11 /= m01;
				m01 = static_cast<Real>(1.0) / std::sqrt(static_cast<Real>(1.0) + m11 * m11);
				m11 *= m01;
			}
			evec1 = m11 * U - m01 * V;
		}
		else
			evec1 = U;
	}
}

} // anonymous namespace

// ----------------------------------------------------------------------------------------------
void MathFunctions::ARAP_eigenvalues(const Matrix3r& F, Vector3r& sigma)
{
	// Compute Cauchy-Green invariants
	const Real I1 = F.squaredNorm();
	const Matrix3r FtF = F.transpose() * F;
	const Real I2 = FtF.squaredNorm();
	const Real J = F.determinant();

	// Solve quartic: t^4 - 2*I1*t^2 - 8*J*t + (2*I2 - I1^2) = 0
	double roots[4];
	solveQuartic(0.0, -2.0 * I1, -8.0 * J, I1 * I1 - 2.0 * (I1 * I1 - I2), roots);

	// Find largest root x4 = σ₁+σ₂+σ₃
	int idx_x4 = 0;
	for (int k = 1; k < 4; k++)
		if (roots[k] > roots[idx_x4]) idx_x4 = k;

	const double x4 = roots[idx_x4];

	// Get the other three roots
	double others[3];
	int oi = 0;
	for (int k = 0; k < 4; k++)
		if (k != idx_x4) others[oi++] = roots[k];

	// Singular values: sigₖ = (xₖ + x4) / 2
	Real sigs[3];
	sigs[0] = static_cast<Real>((others[0] + x4) * 0.5);
	sigs[1] = static_cast<Real>((others[1] + x4) * 0.5);
	sigs[2] = static_cast<Real>((others[2] + x4) * 0.5);

	// Sort ascending
	if (sigs[0] > sigs[1]) std::swap(sigs[0], sigs[1]);
	if (sigs[0] > sigs[2]) std::swap(sigs[0], sigs[2]);
	if (sigs[1] > sigs[2]) std::swap(sigs[1], sigs[2]);

	sigma[0] = sigs[0];  // smallest
	sigma[1] = sigs[1];  // middle
	sigma[2] = sigs[2];  // largest
}

// ----------------------------------------------------------------------------------------------
void MathFunctions::ARAP_decomposition(const Matrix3r& F, Matrix3r& R, Vector3r& sigma, Matrix3r& V)
{
	// Compute Cauchy-Green invariants
	const Real I1 = F.squaredNorm();
	const Matrix3r FtF = F.transpose() * F;
	const Real I2 = FtF.squaredNorm();
	const Real J = F.determinant();

	// Solve quartic: t^4 - 2*I1*t^2 - 8*J*t + (2*I2 - I1^2) = 0
	double roots[4];
	solveQuartic(0.0, -2.0 * I1, -8.0 * J, I1 * I1 - 2.0 * (I1 * I1 - I2), roots);

	// Find largest root x4 = σ₁+σ₂+σ₃ = f
	int idx_x4 = 0;
	for (int k = 1; k < 4; k++)
		if (roots[k] > roots[idx_x4]) idx_x4 = k;

	const double f = roots[idx_x4];

	// Get singular values from quartic roots
	double others[3];
	int oi = 0;
	for (int k = 0; k < 4; k++)
		if (k != idx_x4) others[oi++] = roots[k];

	Real sigs[3];
	sigs[0] = static_cast<Real>((others[0] + f) * 0.5);
	sigs[1] = static_cast<Real>((others[1] + f) * 0.5);
	sigs[2] = static_cast<Real>((others[2] + f) * 0.5);

	// Sort ascending
	if (sigs[0] > sigs[1]) std::swap(sigs[0], sigs[1]);
	if (sigs[0] > sigs[2]) std::swap(sigs[0], sigs[2]);
	if (sigs[1] > sigs[2]) std::swap(sigs[1], sigs[2]);

	sigma[0] = sigs[0];  // sig1 = smallest
	sigma[1] = sigs[1];  // sig2 = middle
	sigma[2] = sigs[2];  // sig3 = largest

	// Compute R via iARAP formula
	const Real denom = 4.0 * f * f * f - 4.0 * I1 * f - 8.0 * J;
	if (std::fabs(denom) < static_cast<Real>(1e-12))
	{
		R.setIdentity();
		V.setIdentity();
		return;
	}

	const Real df1 = (2.0 * f * f + 2.0 * I1) / denom;
	const Real df2 = -2.0 / denom;
	const Real dfJ = (8.0 * f) / denom;

	Matrix3r g1 = 2.0 * F;
	Matrix3r g2 = 4.0 * F * FtF;
	Matrix3r gJ;
	gJ.col(0) = F.col(1).cross(F.col(2));
	gJ.col(1) = F.col(2).cross(F.col(0));
	gJ.col(2) = F.col(0).cross(F.col(1));

	R = df1 * g1 + df2 * g2 + dfJ * gJ;

	// Compute S = R^T F (symmetric)
	Matrix3r S = R.transpose() * F;

	// Compute eigenvectors of S analytically using sigma as eigenvalues
	// Following iARAP reference: V columns are eigenvectors, ordered by sigma
	const Real norm = S(0, 1) * S(0, 1) + S(0, 2) * S(0, 2) + S(1, 2) * S(1, 2);

	if (norm > static_cast<Real>(1e-14))
	{
		Vector3r V0, V1, V2;
		const Real q = (S(0, 0) + S(1, 1) + S(2, 2)) / static_cast<Real>(3.0);
		const Real b00 = S(0, 0) - q;
		const Real b11 = S(1, 1) - q;
		const Real b22 = S(2, 2) - q;
		const Real p = std::sqrt((b00 * b00 + b11 * b11 + b22 * b22 + norm * static_cast<Real>(2.0)) / static_cast<Real>(6.0));
		const Real c00 = b11 * b22 - S(1, 2) * S(1, 2);
		const Real c01 = S(0, 1) * b22 - S(1, 2) * S(0, 2);
		const Real c02 = S(0, 1) * S(1, 2) - b11 * S(0, 2);
		const Real det = (b00 * c00 - S(0, 1) * c01 + S(0, 2) * c02) / (p * p * p);
		Real halfDet = det * static_cast<Real>(0.5);
		halfDet = std::min(std::max(halfDet, static_cast<Real>(-1.0)), static_cast<Real>(1.0));

		if (halfDet >= static_cast<Real>(0.0))
		{
			// sig3 is largest eigenvalue
			computeEigenvector0(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), sigma[2], V2);
			computeEigenvector1(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), V2, sigma[1], V1);
			V0 = V1.cross(V2);
		}
		else
		{
			// sig1 is smallest eigenvalue
			computeEigenvector0(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), sigma[0], V0);
			computeEigenvector1(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), V0, sigma[1], V1);
			V2 = V0.cross(V1);
		}

		// V columns ordered by sigma: V.col(0) -> sig1 (smallest), V.col(2) -> sig3 (largest)
		V.col(0) = V0;
		V.col(1) = V1;
		V.col(2) = V2;
	}
	else
	{
		// S is diagonal
		V.setIdentity();
	}
}
