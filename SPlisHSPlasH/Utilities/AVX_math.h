#ifndef AVX_MATH_H
#define AVX_MATH_H

#include <ctime>
#include <iomanip>
#include <vector>
#include <memory>
#include <immintrin.h>

#include "SPlisHSPlasH/Common.h"


// ----------------------------------------------------------------------------------------------
//vector of 8 float values to represent 8 scalars
class Scalarf8
{
public:
	__m256 v; 

	Scalarf8() {}

	Scalarf8(float f) {	v = _mm256_set1_ps(f); }

// 	Scalarf8(float f0, float f1, float f2, float f3, float f4, float f5, float f6, float f7) {
// 		v = _mm256_setr_ps(f0, f1, f2, f3, f4, f5, f6, f7);
// 	}

	Scalarf8(Real f0, Real f1, Real f2, Real f3, Real f4, Real f5, Real f6, Real f7) 
	{
		v = _mm256_setr_ps((float)f0, (float)f1, (float)f2, (float)f3, 
						   (float)f4, (float)f5, (float)f6, (float)f7);
	}

	Scalarf8(float const * p) 
	{
		v = _mm256_loadu_ps(p);
	}

	Scalarf8(__m256 const & x) {
		v = x;
	}
	
	inline void setZero() { v = _mm256_setzero_ps(); }

	Scalarf8 & operator = (__m256 const & x) {
		v = x;
		return *this;
	}

	inline Scalarf8 sqrt() const 
	{
		return _mm256_sqrt_ps(v);
	}

	inline Scalarf8 rsqrt() const
	{
		return _mm256_rsqrt_ps(v);
	}

	Scalarf8 & load(float const * p) {
		v = _mm256_loadu_ps(p);
		return *this;
	}
	
	void store(float * p) const {
		_mm256_storeu_ps(p, v);
	}

	inline float reduce() const {
		/* ( x3+x7, x2+x6, x1+x5, x0+x4 ) */
		const __m128 x128 = _mm_add_ps(_mm256_extractf128_ps(v, 1), _mm256_castps256_ps128(v));
		/* ( -, -, x1+x3+x5+x7, x0+x2+x4+x6 ) */
		const __m128 x64 = _mm_add_ps(x128, _mm_movehl_ps(x128, x128));
		/* ( -, -, -, x0+x1+x2+x3+x4+x5+x6+x7 ) */
		const __m128 x32 = _mm_add_ss(x64, _mm_shuffle_ps(x64, x64, 0x55));
		/* Conversion to float is a no-op on x86-64 */
		return _mm_cvtss_f32(x32);
	}
};

inline Scalarf8 operator - (Scalarf8& a) {
	return _mm256_sub_ps(_mm256_set1_ps(0.0), a.v);
}


static inline Scalarf8 multiplyAndAdd(const Scalarf8& a, const Scalarf8& b, const Scalarf8& c)
{
	return _mm256_fmadd_ps(a.v, b.v, c.v);
}

static inline Scalarf8 multiplyAndSubtract(const Scalarf8& a, const Scalarf8& b, const Scalarf8& c)
{
	return _mm256_fmsub_ps(a.v, b.v, c.v);
}

static inline Scalarf8 operator + (Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_add_ps(a.v, b.v);
}

static inline Scalarf8 & operator += (Scalarf8 & a, Scalarf8 const & b) {
	a.v = _mm256_add_ps(a.v, b.v);
	return a;
}

static inline Scalarf8 operator - (Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_sub_ps(a.v, b.v);
}

static inline Scalarf8 & operator -= (Scalarf8 & a, Scalarf8 const & b) {
	a.v = _mm256_sub_ps(a.v, b.v);
	return a;
}

static inline Scalarf8 operator * (Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_mul_ps(a.v, b.v);
}

static inline Scalarf8 & operator *= (Scalarf8 & a, Scalarf8 const & b) {
	a.v = _mm256_mul_ps(a.v, b.v);
	return a;
}

static inline Scalarf8 operator / (Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_div_ps(a.v, b.v);
}

static inline Scalarf8 & operator /= (Scalarf8 & a, Scalarf8 const & b) {
	a.v = _mm256_div_ps(a.v, b.v);
	return a;
}

static inline Scalarf8 operator == (Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_cmp_ps(a.v, b.v, 0);
}

static inline Scalarf8 operator != (Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_cmp_ps(a.v, b.v, 4);
}

static inline Scalarf8 operator < (Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_cmp_ps(a.v, b.v, 1);
}

static inline Scalarf8 operator <= (Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_cmp_ps(a.v, b.v, 2);
}

static inline Scalarf8 operator > (Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_cmp_ps(b.v, a.v, 1);
}

static inline Scalarf8 operator >= (Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_cmp_ps(b.v, a.v, 2);
}

static inline Scalarf8 max(Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_max_ps(a.v, b.v);
}

template <int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
static inline __m256 constant8f() {
	static const union {
		int     i[8];
		__m256  ymm;
	} u = { { i0,i1,i2,i3,i4,i5,i6,i7 } };
	return u.ymm;
}

static inline Scalarf8 abs(Scalarf8 const & a) {
	__m256 mask = constant8f<0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF>();
	return _mm256_and_ps(a.v, mask);
}

//does the same as for (int i = 0; i < 8; i++) result[i] = c[i] ? a[i] : b[i];
//the elemets in c must be either 0 (false) or 0xFFFFFFFF (true)
static inline Scalarf8 blend(Scalarf8 const & c, Scalarf8 const & a, Scalarf8 const & b) {
	return _mm256_blendv_ps(b.v, a.v, c.v);
}

static inline Scalarf8 convert_zero(const unsigned int *idx, const Real *x, const unsigned char count = 8u)
{
	Scalarf8 v;
	switch (count)
	{
	case 1u:
		v.v = _mm256_setr_ps(x[idx[0]], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f); break;
	case 2u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f); break;
	case 3u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f); break;
	case 4u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], 0.0f, 0.0f, 0.0f, 0.0f); break;
	case 5u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], x[idx[4]], 0.0f, 0.0f, 0.0f); break;
	case 6u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], x[idx[4]], x[idx[5]], 0.0f, 0.0f); break;
	case 7u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], x[idx[4]], x[idx[5]], x[idx[6]], 0.0f); break;
	case 8u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], x[idx[4]], x[idx[5]], x[idx[6]], x[idx[7]]); break;
	}
	return v;
}

static inline Scalarf8 convert_zero(const Real x, const unsigned char count = 8u)
{
	Scalarf8 v;
	switch (count)
	{
	case 1u:
		v.v = _mm256_setr_ps(x, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f); break;
	case 2u:
		v.v = _mm256_setr_ps(x, x, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f); break;
	case 3u:
		v.v = _mm256_setr_ps(x, x, x, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f); break;
	case 4u:
		v.v = _mm256_setr_ps(x, x, x, x, 0.0f, 0.0f, 0.0f, 0.0f); break;
	case 5u:
		v.v = _mm256_setr_ps(x, x, x, x, x, 0.0f, 0.0f, 0.0f); break;
	case 6u:
		v.v = _mm256_setr_ps(x, x, x, x, x, x, 0.0f, 0.0f); break;
	case 7u:
		v.v = _mm256_setr_ps(x, x, x, x, x, x, x, 0.0f); break;
	case 8u:
		v.v = _mm256_setr_ps(x, x, x, x, x, x, x, x); break;
	}
	return v;
}

static inline Scalarf8 convert_one(const unsigned int *idx, const Real *x, const unsigned char count = 8u)
{
	Scalarf8 v;
	switch (count)
	{
	case 1u:
		v.v = _mm256_setr_ps(x[idx[0]], 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f); break;
	case 2u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f); break;
	case 3u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], 1.0f, 1.0f, 1.0f, 1.0f, 1.0f); break;
	case 4u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], 1.0f, 1.0f, 1.0f, 1.0f); break;
	case 5u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], x[idx[4]], 1.0f, 1.0f, 1.0f); break;
	case 6u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], x[idx[4]], x[idx[5]], 1.0f, 1.0f); break;
	case 7u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], x[idx[4]], x[idx[5]], x[idx[6]], 1.0f); break;
	case 8u:
		v.v = _mm256_setr_ps(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], x[idx[4]], x[idx[5]], x[idx[6]], x[idx[7]]); break;
	}
	return v;
}

// ----------------------------------------------------------------------------------------------
//3 dimensional vector of Scalar8f to represent 8 3d vectors
class Vector3f8
{
public:

	Scalarf8 v[3];

	Vector3f8() {}
	Vector3f8(const bool ) { v[0] = Scalarf8(0.0f); v[1] = Scalarf8(0.0f); v[2] = Scalarf8(0.0f); }
	Vector3f8(const Scalarf8 &x, const Scalarf8 &y, const Scalarf8 &z) { v[0] = x; v[1] = y; v[2] = z; }
	Vector3f8(const Scalarf8 &x) { v[0] = v[1] = v[2] = x; }
	Vector3f8(const Vector3f &x) 
	{
		v[0].v = _mm256_set1_ps(x[0]); 
		v[1].v = _mm256_set1_ps(x[1]);
		v[2].v = _mm256_set1_ps(x[2]);
	}
	Vector3f8(const Vector3f &v0, const Vector3f &v1, const Vector3f &v2, const Vector3f &v3, const Vector3f &v4, const Vector3f &v5, const Vector3f &v6, const Vector3f &v7)
	{
		Scalarf8 x(v0[0], v1[0], v2[0], v3[0], v4[0], v5[0], v6[0], v7[0]);
		Scalarf8 y(v0[1], v1[1], v2[1], v3[1], v4[1], v5[1], v6[1], v7[1]);
		Scalarf8 z(v0[2], v1[2], v2[2], v3[2], v4[2], v5[2], v6[2], v7[2]);
		v[0] = x; v[1] = y; v[2] = z;
	}
	Vector3f8(Vector3f const *x)
	{
		Scalarf8 vx(x[0][0], x[1][0], x[2][0], x[3][0], x[4][0], x[5][0], x[6][0], x[7][0]);
		Scalarf8 vy(x[0][1], x[1][1], x[2][1], x[3][1], x[4][1], x[5][1], x[6][1], x[7][1]);
		Scalarf8 vz(x[0][2], x[1][2], x[2][2], x[3][2], x[4][2], x[5][2], x[6][2], x[7][2]);
		v[0] = vx; v[1] = vy; v[2] = vz;
	}

	inline void setZero() { v[0].v = _mm256_setzero_ps(); v[1].v = _mm256_setzero_ps(); v[2].v = _mm256_setzero_ps();
	}

	inline Scalarf8& operator [] (int i) { return v[i]; }
	inline const Scalarf8& operator [] (int i) const { return v[i]; }

	inline Scalarf8& x() { return v[0]; }
	inline Scalarf8& y() { return v[1]; }
	inline Scalarf8& z() { return v[2]; }

	inline const Scalarf8& x() const { return v[0]; }
	inline const Scalarf8& y() const { return v[1]; }
	inline const Scalarf8& z() const { return v[2]; }
	
	inline Scalarf8 dot(const Vector3f8& a) const {
		Scalarf8 res; 
		//res.v = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(v[0].v, a.v[0].v), _mm256_mul_ps(v[1].v, a.v[1].v)), _mm256_mul_ps(v[2].v, a.v[2].v));
		res.v = _mm256_fmadd_ps(v[0].v, a.v[0].v, _mm256_fmadd_ps(v[1].v, a.v[1].v, _mm256_mul_ps(v[2].v, a.v[2].v)));
		return res;
	}

	//dot product
	inline Scalarf8 operator * (const Vector3f8& a) const {
		Scalarf8 res;
		//res.v = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(v[0].v, a.v[0].v), _mm256_mul_ps(v[1].v, a.v[1].v)), _mm256_mul_ps(v[2].v, a.v[2].v));
		res.v = _mm256_fmadd_ps(v[0].v, a.v[0].v, _mm256_fmadd_ps(v[1].v, a.v[1].v, _mm256_mul_ps(v[2].v, a.v[2].v)));
		return res;
	}

	inline void cross(const Vector3f8& a, const Vector3f8& b) {
		v[0] = a.v[1] * b.v[2] - a.v[2] * b.v[1];
		v[1] = a.v[2] * b.v[0] - a.v[0] * b.v[2];
		v[2] = a.v[0] * b.v[1] - a.v[1] * b.v[0];
	}

	//cross product
	inline const Vector3f8 operator % (const Vector3f8& a) const {
		return Vector3f8(v[1] * a.v[2] - v[2] * a.v[1],
			v[2] * a.v[0] - v[0] * a.v[2],
			v[0] * a.v[1] - v[1] * a.v[0]);
	}

	inline Vector3f8& operator *= (const Scalarf8 &s) {
		v[0] *= s;
		v[1] *= s;
		v[2] *= s;
		return *this;
	}

	inline const Vector3f8 operator / (const Scalarf8 &s) const {
		return Vector3f8(v[0] / s, v[1] / s, v[2] / s);
	}

	inline Vector3f8& operator /= (const Scalarf8 &s) {
		v[0] = v[0] / s;
		v[1] = v[1] / s;
		v[2] = v[2] / s;
		return *this;
	}

	inline const Vector3f8 operator - () const {
		return Vector3f8(Scalarf8(-1.0) * v[0], Scalarf8(-1.0) * v[1], Scalarf8(-1.0) * v[2]);
	}

	inline Scalarf8 squaredNorm() const {
		//return _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(v[0].v, v[0].v), _mm256_mul_ps(v[1].v, v[1].v)), _mm256_mul_ps(v[2].v, v[2].v));
		return _mm256_fmadd_ps(v[0].v, v[0].v, _mm256_fmadd_ps(v[1].v, v[1].v, _mm256_mul_ps(v[2].v, v[2].v)));
	}

	inline Scalarf8 norm() const
	{
		//return _mm256_sqrt_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(v[0].v, v[0].v), _mm256_mul_ps(v[1].v, v[1].v)), _mm256_mul_ps(v[2].v, v[2].v)));
		return _mm256_sqrt_ps(_mm256_fmadd_ps(v[0].v, v[0].v, _mm256_fmadd_ps(v[1].v, v[1].v, _mm256_mul_ps(v[2].v, v[2].v))));
	}

	inline void normalize()
	{
		const Scalarf8 norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
		v[0] = v[0] / norm;
		v[1] = v[1] / norm;
		v[2] = v[2] / norm;
		*this = blend(norm >= Scalarf8(1e-5f), *this, Vector3f8(true));
	}

	//does the same as for (int i = 0; i < 8; i++) result[i] = c[i] ? a[i] : b[i];
	//the elements in c must be either 0 (false) or 0xFFFFFFFF (true)
	static inline Vector3f8 blend(Scalarf8 const & c, Vector3f8 const & a, Vector3f8 const & b) {
		Vector3f8 result;
		result.x() = _mm256_blendv_ps(b.x().v, a.x().v, c.v);
		result.y() = _mm256_blendv_ps(b.y().v, a.y().v, c.v);
		result.z() = _mm256_blendv_ps(b.z().v, a.z().v, c.v);
		return result;
	}

	inline void store(std::vector<Vector3r>& Vf) const
	{
		for (int i = 0; i < 3; i++)
		{
			float val[8];
			v[i].store(val);
			for (int k = 0; k < 8; k++)
				Vf[k][i] = val[k];
		}
	}

	inline void store(Vector3r* Vf) const
	{
		for (int i = 0; i < 3; i++)
		{
			float val[8];
			v[i].store(val);
			for (int k = 0; k < 8; k++)
				Vf[k][i] = val[k];
		}
	}

	inline Vector3r reduce() const {
		Vector3r res;
		res[0] = v[0].reduce();
		res[1] = v[1].reduce();
		res[2] = v[2].reduce();
		return res;
	}
};

inline Vector3f8 operator + (Vector3f8 const& a, Vector3f8 const& b) {
	Vector3f8 res;
	res.v[0].v = _mm256_add_ps(a[0].v, b[0].v);
	res.v[1].v = _mm256_add_ps(a[1].v, b[1].v);
	res.v[2].v = _mm256_add_ps(a[2].v, b[2].v);
	return res;
}

inline Vector3f8 operator - (Vector3f8 const& a, Vector3f8 const& b) {
	Vector3f8 res;
	res.v[0].v = _mm256_sub_ps(a[0].v, b[0].v);
	res.v[1].v = _mm256_sub_ps(a[1].v, b[1].v);
	res.v[2].v = _mm256_sub_ps(a[2].v, b[2].v);
	return res;
}

inline Vector3f8& operator += (Vector3f8& a, Vector3f8 const& b) {
	a[0].v = _mm256_add_ps(a[0].v, b[0].v);
	a[1].v = _mm256_add_ps(a[1].v, b[1].v);
	a[2].v = _mm256_add_ps(a[2].v, b[2].v);
	return a;
}

inline Vector3f8& operator -= (Vector3f8& a, Vector3f8 const& b) {
	a[0].v = _mm256_sub_ps(a[0].v, b[0].v);
	a[1].v = _mm256_sub_ps(a[1].v, b[1].v);
	a[2].v = _mm256_sub_ps(a[2].v, b[2].v);
	return a;
}

inline Vector3f8 operator * (Vector3f8 const& a, const Scalarf8& s) {
	Vector3f8 res;
	res.v[0].v = _mm256_mul_ps(a[0].v, s.v);
	res.v[1].v = _mm256_mul_ps(a[1].v, s.v);
	res.v[2].v = _mm256_mul_ps(a[2].v, s.v);
	return res;
}

inline Vector3f8 convertVec_zero(const unsigned int* idx, const Real* v, const unsigned char count = 8u)
{
	Vector3f8 x;
	switch (count)
	{
	case 1u:
		x.v[0].v = _mm256_setr_ps(v[3*idx[0]+0], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[3*idx[0]+1], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[3*idx[0]+2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 2u:
		x.v[0].v = _mm256_setr_ps(v[3*idx[0]+0], v[3*idx[1]+0], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[3*idx[0]+1], v[3*idx[1]+1], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[3*idx[0]+2], v[3*idx[1]+2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 3u:
		x.v[0].v = _mm256_setr_ps(v[3*idx[0]+0], v[3*idx[1]+0], v[3*idx[2]+0], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[3*idx[0]+1], v[3*idx[1]+1], v[3*idx[2]+1], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[3*idx[0]+2], v[3*idx[1]+2], v[3*idx[2]+2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 4u:
		x.v[0].v = _mm256_setr_ps(v[3*idx[0]+0], v[3*idx[1]+0], v[3*idx[2]+0], v[3*idx[3]+0], 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[3*idx[0]+1], v[3*idx[1]+1], v[3*idx[2]+1], v[3*idx[3]+1], 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[3*idx[0]+2], v[3*idx[1]+2], v[3*idx[2]+2], v[3*idx[3]+2], 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 5u:
		x.v[0].v = _mm256_setr_ps(v[3*idx[0]+0], v[3*idx[1]+0], v[3*idx[2]+0], v[3*idx[3]+0], v[3*idx[4]+0], 0.0f, 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[3*idx[0]+1], v[3*idx[1]+1], v[3*idx[2]+1], v[3*idx[3]+1], v[3*idx[4]+1], 0.0f, 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[3*idx[0]+2], v[3*idx[1]+2], v[3*idx[2]+2], v[3*idx[3]+2], v[3*idx[4]+2], 0.0f, 0.0f, 0.0f);
		break;
	case 6u:
		x.v[0].v = _mm256_setr_ps(v[3*idx[0]+0], v[3*idx[1]+0], v[3*idx[2]+0], v[3*idx[3]+0], v[3*idx[4]+0], v[3*idx[5]+0], 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[3*idx[0]+1], v[3*idx[1]+1], v[3*idx[2]+1], v[3*idx[3]+1], v[3*idx[4]+1], v[3*idx[5]+1], 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[3*idx[0]+2], v[3*idx[1]+2], v[3*idx[2]+2], v[3*idx[3]+2], v[3*idx[4]+2], v[3*idx[5]+2], 0.0f, 0.0f);
		break;
	case 7u:
		x.v[0].v = _mm256_setr_ps(v[3*idx[0]+0], v[3*idx[1]+0], v[3*idx[2]+0], v[3*idx[3]+0], v[3*idx[4]+0], v[3*idx[5]+0], v[3*idx[6]+0], 0.0f);
		x.v[1].v = _mm256_setr_ps(v[3*idx[0]+1], v[3*idx[1]+1], v[3*idx[2]+1], v[3*idx[3]+1], v[3*idx[4]+1], v[3*idx[5]+1], v[3*idx[6]+1], 0.0f);
		x.v[2].v = _mm256_setr_ps(v[3*idx[0]+2], v[3*idx[1]+2], v[3*idx[2]+2], v[3*idx[3]+2], v[3*idx[4]+2], v[3*idx[5]+2], v[3*idx[6]+2], 0.0f);
		break;
	case 8u:
		x.v[0].v = _mm256_setr_ps(v[3*idx[0]+0], v[3*idx[1]+0], v[3*idx[2]+0], v[3*idx[3]+0], v[3*idx[4]+0], v[3*idx[5]+0], v[3*idx[6]+0], v[3*idx[7]+0]);
		x.v[1].v = _mm256_setr_ps(v[3*idx[0]+1], v[3*idx[1]+1], v[3*idx[2]+1], v[3*idx[3]+1], v[3*idx[4]+1], v[3*idx[5]+1], v[3*idx[6]+1], v[3*idx[7]+1]);
		x.v[2].v = _mm256_setr_ps(v[3*idx[0]+2], v[3*idx[1]+2], v[3*idx[2]+2], v[3*idx[3]+2], v[3*idx[4]+2], v[3*idx[5]+2], v[3*idx[6]+2], v[3*idx[7]+2]);
	}
	return x;
}

static inline Vector3f8 convertVec_zero(const unsigned int *idx, const Vector3r *v, const unsigned char count = 8u)
{
	Vector3f8 x;
	switch (count)
	{
	case 1u:
		x.v[0].v = _mm256_setr_ps(v[idx[0]][0], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[idx[0]][1], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[idx[0]][2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 2u:
		x.v[0].v = _mm256_setr_ps(v[idx[0]][0], v[idx[1]][0], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[idx[0]][1], v[idx[1]][1], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[idx[0]][2], v[idx[1]][2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 3u:
		x.v[0].v = _mm256_setr_ps(v[idx[0]][0], v[idx[1]][0], v[idx[2]][0], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[idx[0]][1], v[idx[1]][1], v[idx[2]][1], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[idx[0]][2], v[idx[1]][2], v[idx[2]][2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 4u:
		x.v[0].v = _mm256_setr_ps(v[idx[0]][0], v[idx[1]][0], v[idx[2]][0], v[idx[3]][0], 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[idx[0]][1], v[idx[1]][1], v[idx[2]][1], v[idx[3]][1], 0.0f, 0.0f, 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[idx[0]][2], v[idx[1]][2], v[idx[2]][2], v[idx[3]][2], 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 5u:
		x.v[0].v = _mm256_setr_ps(v[idx[0]][0], v[idx[1]][0], v[idx[2]][0], v[idx[3]][0], v[idx[4]][0], 0.0f, 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[idx[0]][1], v[idx[1]][1], v[idx[2]][1], v[idx[3]][1], v[idx[4]][1], 0.0f, 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[idx[0]][2], v[idx[1]][2], v[idx[2]][2], v[idx[3]][2], v[idx[4]][2], 0.0f, 0.0f, 0.0f);
		break;
	case 6u:
		x.v[0].v = _mm256_setr_ps(v[idx[0]][0], v[idx[1]][0], v[idx[2]][0], v[idx[3]][0], v[idx[4]][0], v[idx[5]][0], 0.0f, 0.0f);
		x.v[1].v = _mm256_setr_ps(v[idx[0]][1], v[idx[1]][1], v[idx[2]][1], v[idx[3]][1], v[idx[4]][1], v[idx[5]][1], 0.0f, 0.0f);
		x.v[2].v = _mm256_setr_ps(v[idx[0]][2], v[idx[1]][2], v[idx[2]][2], v[idx[3]][2], v[idx[4]][2], v[idx[5]][2], 0.0f, 0.0f);
		break;
	case 7u:
		x.v[0].v = _mm256_setr_ps(v[idx[0]][0], v[idx[1]][0], v[idx[2]][0], v[idx[3]][0], v[idx[4]][0], v[idx[5]][0], v[idx[6]][0], 0.0f);
		x.v[1].v = _mm256_setr_ps(v[idx[0]][1], v[idx[1]][1], v[idx[2]][1], v[idx[3]][1], v[idx[4]][1], v[idx[5]][1], v[idx[6]][1], 0.0f);
		x.v[2].v = _mm256_setr_ps(v[idx[0]][2], v[idx[1]][2], v[idx[2]][2], v[idx[3]][2], v[idx[4]][2], v[idx[5]][2], v[idx[6]][2], 0.0f);
		break;
	case 8u:
		x.v[0].v = _mm256_setr_ps(v[idx[0]][0], v[idx[1]][0], v[idx[2]][0], v[idx[3]][0], v[idx[4]][0], v[idx[5]][0], v[idx[6]][0], v[idx[7]][0]);
		x.v[1].v = _mm256_setr_ps(v[idx[0]][1], v[idx[1]][1], v[idx[2]][1], v[idx[3]][1], v[idx[4]][1], v[idx[5]][1], v[idx[6]][1], v[idx[7]][1]);
		x.v[2].v = _mm256_setr_ps(v[idx[0]][2], v[idx[1]][2], v[idx[2]][2], v[idx[3]][2], v[idx[4]][2], v[idx[5]][2], v[idx[6]][2], v[idx[7]][2]);
	}
	return x;
}

// ----------------------------------------------------------------------------------------------
//3x3 dimensional matrix of Scalar8f to represent 8 3x3 matrices
class Matrix3f8
{
public:
	Scalarf8 m[3][3];

	Matrix3f8() {  }

	Matrix3f8(const Matrix3f &x)
	{
		const Vector3f8 m1(x.col(0));
		const Vector3f8 m2(x.col(1));
		const Vector3f8 m3(x.col(2));
		m[0][0] = m1.x();
		m[1][0] = m1.y();
		m[2][0] = m1.z();

		m[0][1] = m2.x();
		m[1][1] = m2.y();
		m[2][1] = m2.z();

		m[0][2] = m3.x();
		m[1][2] = m3.y();
		m[2][2] = m3.z();
	}

	//constructor to create matrix from 3 column vectors
	Matrix3f8(const Vector3f8& m1, const Vector3f8& m2, const Vector3f8& m3)
	{
		m[0][0] = m1.x();
		m[1][0] = m1.y();
		m[2][0] = m1.z();

		m[0][1] = m2.x();
		m[1][1] = m2.y();
		m[2][1] = m2.z();

		m[0][2] = m3.x();
		m[1][2] = m3.y();
		m[2][2] = m3.z();
	}

	inline void setZero() 
	{ 
		m[0][0].v = _mm256_setzero_ps();
		m[1][0].v = _mm256_setzero_ps();
		m[2][0].v = _mm256_setzero_ps();

		m[0][1].v = _mm256_setzero_ps();
		m[1][1].v = _mm256_setzero_ps();
		m[2][1].v = _mm256_setzero_ps();

		m[0][2].v = _mm256_setzero_ps();
		m[1][2].v = _mm256_setzero_ps();
		m[2][2].v = _mm256_setzero_ps();
	}

	inline Scalarf8& operator()(int i, int j) { return m[i][j]; }

	inline void setCol(int i, const Vector3f8& v)
	{
		m[0][i] = v.x();
		m[1][i] = v.y();
		m[2][i] = v.z();
	}

	inline void setCol(int i, const Scalarf8& x, const Scalarf8& y, const Scalarf8& z)
	{
		m[0][i] = x;
		m[1][i] = y;
		m[2][i] = z;
	}

	inline Matrix3f8 operator * (const Scalarf8 &b) const
	{
		Matrix3f8 A;

		A.m[0][0] = m[0][0] * b;
		A.m[0][1] = m[0][1] * b;
		A.m[0][2] = m[0][2] * b;

		A.m[1][0] = m[1][0] * b;
		A.m[1][1] = m[1][1] * b;
		A.m[1][2] = m[1][2] * b;

		A.m[2][0] = m[2][0] * b;
		A.m[2][1] = m[2][1] * b;
		A.m[2][2] = m[2][2] * b;

		return A;
	}

	inline Vector3f8 operator * (const Vector3f8 &b) const
	{
		Vector3f8 A;

		A.v[0] = m[0][0] * b.v[0] + m[0][1] * b.v[1] + m[0][2] * b.v[2];
		A.v[1] = m[1][0] * b.v[0] + m[1][1] * b.v[1] + m[1][2] * b.v[2];
		A.v[2] = m[2][0] * b.v[0] + m[2][1] * b.v[1] + m[2][2] * b.v[2];

		return A;
	}

	inline Matrix3f8 operator * (const Matrix3f8 &b) const
	{
		Matrix3f8 A;

		A.m[0][0] = m[0][0] * b.m[0][0] + m[0][1] * b.m[1][0] + m[0][2] * b.m[2][0];
		A.m[0][1] = m[0][0] * b.m[0][1] + m[0][1] * b.m[1][1] + m[0][2] * b.m[2][1];
		A.m[0][2] = m[0][0] * b.m[0][2] + m[0][1] * b.m[1][2] + m[0][2] * b.m[2][2];

		A.m[1][0] = m[1][0] * b.m[0][0] + m[1][1] * b.m[1][0] + m[1][2] * b.m[2][0];
		A.m[1][1] = m[1][0] * b.m[0][1] + m[1][1] * b.m[1][1] + m[1][2] * b.m[2][1];
		A.m[1][2] = m[1][0] * b.m[0][2] + m[1][1] * b.m[1][2] + m[1][2] * b.m[2][2];

		A.m[2][0] = m[2][0] * b.m[0][0] + m[2][1] * b.m[1][0] + m[2][2] * b.m[2][0];
		A.m[2][1] = m[2][0] * b.m[0][1] + m[2][1] * b.m[1][1] + m[2][2] * b.m[2][1];
		A.m[2][2] = m[2][0] * b.m[0][2] + m[2][1] * b.m[1][2] + m[2][2] * b.m[2][2];

		return A;
	}

	inline Matrix3f8& operator += (const Matrix3f8& a) {
		m[0][0] += a.m[0][0];
		m[1][0] += a.m[1][0];
		m[2][0] += a.m[2][0];

		m[0][1] += a.m[0][1];
		m[1][1] += a.m[1][1];
		m[2][1] += a.m[2][1];

		m[0][2] += a.m[0][2];
		m[1][2] += a.m[1][2];
		m[2][2] += a.m[2][2];
		return *this;
	}

	inline Matrix3f8 transpose() const
	{
		Matrix3f8 A;
		A.m[0][0] = m[0][0]; A.m[0][1] = m[1][0]; A.m[0][2] = m[2][0];
		A.m[1][0] = m[0][1]; A.m[1][1] = m[1][1]; A.m[1][2] = m[2][1];
		A.m[2][0] = m[0][2]; A.m[2][1] = m[1][2]; A.m[2][2] = m[2][2];

		return A;
	}

	inline Scalarf8 determinant() const
	{
		return  m[0][1] * m[1][2] * m[2][0] - m[0][2] * m[1][1] * m[2][0] + m[0][2] * m[1][0] * m[2][1] 
			  - m[0][0] * m[1][2] * m[2][1] - m[0][1] * m[1][0] * m[2][2] + m[0][0] * m[1][1] * m[2][2];
	}

	inline void store(std::vector<Matrix3r>& Mf) const
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				float val[8];
				m[i][j].store(val);
				for (int k = 0; k < 8; k++)
					Mf[k](i, j) = val[k];
			}
		}
	}

	inline Matrix3r reduce() const {
		Matrix3r A;
		A(0, 0) = m[0][0].reduce();
		A(0, 1) = m[0][1].reduce();
		A(0, 2) = m[0][2].reduce();

		A(1, 0) = m[1][0].reduce();
		A(1, 1) = m[1][1].reduce();
		A(1, 2) = m[1][2].reduce();

		A(2, 0) = m[2][0].reduce();
		A(2, 1) = m[2][1].reduce();
		A(2, 2) = m[2][2].reduce();
		return A;
	}
};

inline Matrix3f8& operator -= (Matrix3f8& a, Matrix3f8 const& b) {
	a.m[0][0].v = _mm256_sub_ps(a.m[0][0].v, b.m[0][0].v);
	a.m[0][1].v = _mm256_sub_ps(a.m[0][1].v, b.m[0][1].v);
	a.m[0][2].v = _mm256_sub_ps(a.m[0][2].v, b.m[0][2].v);

	a.m[1][0].v = _mm256_sub_ps(a.m[1][0].v, b.m[1][0].v);
	a.m[1][1].v = _mm256_sub_ps(a.m[1][1].v, b.m[1][1].v);
	a.m[1][2].v = _mm256_sub_ps(a.m[1][2].v, b.m[1][2].v);

	a.m[2][0].v = _mm256_sub_ps(a.m[2][0].v, b.m[2][0].v);
	a.m[2][1].v = _mm256_sub_ps(a.m[2][1].v, b.m[2][1].v);
	a.m[2][2].v = _mm256_sub_ps(a.m[2][2].v, b.m[2][2].v);

	return a;
}

static inline Matrix3f8 convertMat_zero(const unsigned int* idx, const Matrix3r* v, const unsigned char count = 8u)
{
	Matrix3f8 x;
	switch (count)
	{
	case 1u:
		x.m[0][0] = _mm256_setr_ps(v[idx[0]](0, 0), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[0][1] = _mm256_setr_ps(v[idx[0]](0, 1), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[0][2] = _mm256_setr_ps(v[idx[0]](0, 2), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][0] = _mm256_setr_ps(v[idx[0]](1, 0), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][1] = _mm256_setr_ps(v[idx[0]](1, 1), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][2] = _mm256_setr_ps(v[idx[0]](1, 2), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][0] = _mm256_setr_ps(v[idx[0]](2, 0), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][1] = _mm256_setr_ps(v[idx[0]](2, 1), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][2] = _mm256_setr_ps(v[idx[0]](2, 2), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 2u:
		x.m[0][0] = _mm256_setr_ps(v[idx[0]](0, 0), v[idx[1]](0, 0), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[0][1] = _mm256_setr_ps(v[idx[0]](0, 1), v[idx[1]](0, 1), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[0][2] = _mm256_setr_ps(v[idx[0]](0, 2), v[idx[1]](0, 2), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][0] = _mm256_setr_ps(v[idx[0]](1, 0), v[idx[1]](1, 0), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][1] = _mm256_setr_ps(v[idx[0]](1, 1), v[idx[1]](1, 1), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][2] = _mm256_setr_ps(v[idx[0]](1, 2), v[idx[1]](1, 2), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][0] = _mm256_setr_ps(v[idx[0]](2, 0), v[idx[1]](2, 0), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][1] = _mm256_setr_ps(v[idx[0]](2, 1), v[idx[1]](2, 1), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][2] = _mm256_setr_ps(v[idx[0]](2, 2), v[idx[1]](2, 2), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 3u:
		x.m[0][0] = _mm256_setr_ps(v[idx[0]](0, 0), v[idx[1]](0, 0), v[idx[2]](0, 0), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[0][1] = _mm256_setr_ps(v[idx[0]](0, 1), v[idx[1]](0, 1), v[idx[2]](0, 1), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[0][2] = _mm256_setr_ps(v[idx[0]](0, 2), v[idx[1]](0, 2), v[idx[2]](0, 2), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][0] = _mm256_setr_ps(v[idx[0]](1, 0), v[idx[1]](1, 0), v[idx[2]](1, 0), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][1] = _mm256_setr_ps(v[idx[0]](1, 1), v[idx[1]](1, 1), v[idx[2]](1, 1), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][2] = _mm256_setr_ps(v[idx[0]](1, 2), v[idx[1]](1, 2), v[idx[2]](1, 2), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][0] = _mm256_setr_ps(v[idx[0]](2, 0), v[idx[1]](2, 0), v[idx[2]](2, 0), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][1] = _mm256_setr_ps(v[idx[0]](2, 1), v[idx[1]](2, 1), v[idx[2]](2, 1), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][2] = _mm256_setr_ps(v[idx[0]](2, 2), v[idx[1]](2, 2), v[idx[2]](2, 2), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 4u:
		x.m[0][0] = _mm256_setr_ps(v[idx[0]](0, 0), v[idx[1]](0, 0), v[idx[2]](0, 0), v[idx[3]](0, 0), 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[0][1] = _mm256_setr_ps(v[idx[0]](0, 1), v[idx[1]](0, 1), v[idx[2]](0, 1), v[idx[3]](0, 1), 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[0][2] = _mm256_setr_ps(v[idx[0]](0, 2), v[idx[1]](0, 2), v[idx[2]](0, 2), v[idx[3]](0, 2), 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][0] = _mm256_setr_ps(v[idx[0]](1, 0), v[idx[1]](1, 0), v[idx[2]](1, 0), v[idx[3]](1, 0), 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][1] = _mm256_setr_ps(v[idx[0]](1, 1), v[idx[1]](1, 1), v[idx[2]](1, 1), v[idx[3]](1, 1), 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[1][2] = _mm256_setr_ps(v[idx[0]](1, 2), v[idx[1]](1, 2), v[idx[2]](1, 2), v[idx[3]](1, 2), 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][0] = _mm256_setr_ps(v[idx[0]](2, 0), v[idx[1]](2, 0), v[idx[2]](2, 0), v[idx[3]](2, 0), 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][1] = _mm256_setr_ps(v[idx[0]](2, 1), v[idx[1]](2, 1), v[idx[2]](2, 1), v[idx[3]](2, 1), 0.0f, 0.0f, 0.0f, 0.0f);
		x.m[2][2] = _mm256_setr_ps(v[idx[0]](2, 2), v[idx[1]](2, 2), v[idx[2]](2, 2), v[idx[3]](2, 2), 0.0f, 0.0f, 0.0f, 0.0f);
		break;
	case 5u:
		x.m[0][0] = _mm256_setr_ps(v[idx[0]](0, 0), v[idx[1]](0, 0), v[idx[2]](0, 0), v[idx[3]](0, 0), v[idx[4]](0, 0), 0.0f, 0.0f, 0.0f);
		x.m[0][1] = _mm256_setr_ps(v[idx[0]](0, 1), v[idx[1]](0, 1), v[idx[2]](0, 1), v[idx[3]](0, 1), v[idx[4]](0, 1), 0.0f, 0.0f, 0.0f);
		x.m[0][2] = _mm256_setr_ps(v[idx[0]](0, 2), v[idx[1]](0, 2), v[idx[2]](0, 2), v[idx[3]](0, 2), v[idx[4]](0, 2), 0.0f, 0.0f, 0.0f);
		x.m[1][0] = _mm256_setr_ps(v[idx[0]](1, 0), v[idx[1]](1, 0), v[idx[2]](1, 0), v[idx[3]](1, 0), v[idx[4]](1, 0), 0.0f, 0.0f, 0.0f);
		x.m[1][1] = _mm256_setr_ps(v[idx[0]](1, 1), v[idx[1]](1, 1), v[idx[2]](1, 1), v[idx[3]](1, 1), v[idx[4]](1, 1), 0.0f, 0.0f, 0.0f);
		x.m[1][2] = _mm256_setr_ps(v[idx[0]](1, 2), v[idx[1]](1, 2), v[idx[2]](1, 2), v[idx[3]](1, 2), v[idx[4]](1, 2), 0.0f, 0.0f, 0.0f);
		x.m[2][0] = _mm256_setr_ps(v[idx[0]](2, 0), v[idx[1]](2, 0), v[idx[2]](2, 0), v[idx[3]](2, 0), v[idx[4]](2, 0), 0.0f, 0.0f, 0.0f);
		x.m[2][1] = _mm256_setr_ps(v[idx[0]](2, 1), v[idx[1]](2, 1), v[idx[2]](2, 1), v[idx[3]](2, 1), v[idx[4]](2, 1), 0.0f, 0.0f, 0.0f);
		x.m[2][2] = _mm256_setr_ps(v[idx[0]](2, 2), v[idx[1]](2, 2), v[idx[2]](2, 2), v[idx[3]](2, 2), v[idx[4]](2, 2), 0.0f, 0.0f, 0.0f);
		break;
	case 6u:
		x.m[0][0] = _mm256_setr_ps(v[idx[0]](0, 0), v[idx[1]](0, 0), v[idx[2]](0, 0), v[idx[3]](0, 0), v[idx[4]](0, 0), v[idx[5]](0, 0), 0.0f, 0.0f);
		x.m[0][1] = _mm256_setr_ps(v[idx[0]](0, 1), v[idx[1]](0, 1), v[idx[2]](0, 1), v[idx[3]](0, 1), v[idx[4]](0, 1), v[idx[5]](0, 1), 0.0f, 0.0f);
		x.m[0][2] = _mm256_setr_ps(v[idx[0]](0, 2), v[idx[1]](0, 2), v[idx[2]](0, 2), v[idx[3]](0, 2), v[idx[4]](0, 2), v[idx[5]](0, 2), 0.0f, 0.0f);
		x.m[1][0] = _mm256_setr_ps(v[idx[0]](1, 0), v[idx[1]](1, 0), v[idx[2]](1, 0), v[idx[3]](1, 0), v[idx[4]](1, 0), v[idx[5]](1, 0), 0.0f, 0.0f);
		x.m[1][1] = _mm256_setr_ps(v[idx[0]](1, 1), v[idx[1]](1, 1), v[idx[2]](1, 1), v[idx[3]](1, 1), v[idx[4]](1, 1), v[idx[5]](1, 1), 0.0f, 0.0f);
		x.m[1][2] = _mm256_setr_ps(v[idx[0]](1, 2), v[idx[1]](1, 2), v[idx[2]](1, 2), v[idx[3]](1, 2), v[idx[4]](1, 2), v[idx[5]](1, 2), 0.0f, 0.0f);
		x.m[2][0] = _mm256_setr_ps(v[idx[0]](2, 0), v[idx[1]](2, 0), v[idx[2]](2, 0), v[idx[3]](2, 0), v[idx[4]](2, 0), v[idx[5]](2, 0), 0.0f, 0.0f);
		x.m[2][1] = _mm256_setr_ps(v[idx[0]](2, 1), v[idx[1]](2, 1), v[idx[2]](2, 1), v[idx[3]](2, 1), v[idx[4]](2, 1), v[idx[5]](2, 1), 0.0f, 0.0f);
		x.m[2][2] = _mm256_setr_ps(v[idx[0]](2, 2), v[idx[1]](2, 2), v[idx[2]](2, 2), v[idx[3]](2, 2), v[idx[4]](2, 2), v[idx[5]](2, 2), 0.0f, 0.0f);
		break;
	case 7u:
		x.m[0][0] = _mm256_setr_ps(v[idx[0]](0, 0), v[idx[1]](0, 0), v[idx[2]](0, 0), v[idx[3]](0, 0), v[idx[4]](0, 0), v[idx[5]](0, 0), v[idx[6]](0, 0), 0.0f);
		x.m[0][1] = _mm256_setr_ps(v[idx[0]](0, 1), v[idx[1]](0, 1), v[idx[2]](0, 1), v[idx[3]](0, 1), v[idx[4]](0, 1), v[idx[5]](0, 1), v[idx[6]](0, 1), 0.0f);
		x.m[0][2] = _mm256_setr_ps(v[idx[0]](0, 2), v[idx[1]](0, 2), v[idx[2]](0, 2), v[idx[3]](0, 2), v[idx[4]](0, 2), v[idx[5]](0, 2), v[idx[6]](0, 2), 0.0f);
		x.m[1][0] = _mm256_setr_ps(v[idx[0]](1, 0), v[idx[1]](1, 0), v[idx[2]](1, 0), v[idx[3]](1, 0), v[idx[4]](1, 0), v[idx[5]](1, 0), v[idx[6]](1, 0), 0.0f);
		x.m[1][1] = _mm256_setr_ps(v[idx[0]](1, 1), v[idx[1]](1, 1), v[idx[2]](1, 1), v[idx[3]](1, 1), v[idx[4]](1, 1), v[idx[5]](1, 1), v[idx[6]](1, 1), 0.0f);
		x.m[1][2] = _mm256_setr_ps(v[idx[0]](1, 2), v[idx[1]](1, 2), v[idx[2]](1, 2), v[idx[3]](1, 2), v[idx[4]](1, 2), v[idx[5]](1, 2), v[idx[6]](1, 2), 0.0f);
		x.m[2][0] = _mm256_setr_ps(v[idx[0]](2, 0), v[idx[1]](2, 0), v[idx[2]](2, 0), v[idx[3]](2, 0), v[idx[4]](2, 0), v[idx[5]](2, 0), v[idx[6]](2, 0), 0.0f);
		x.m[2][1] = _mm256_setr_ps(v[idx[0]](2, 1), v[idx[1]](2, 1), v[idx[2]](2, 1), v[idx[3]](2, 1), v[idx[4]](2, 1), v[idx[5]](2, 1), v[idx[6]](2, 1), 0.0f);
		x.m[2][2] = _mm256_setr_ps(v[idx[0]](2, 2), v[idx[1]](2, 2), v[idx[2]](2, 2), v[idx[3]](2, 2), v[idx[4]](2, 2), v[idx[5]](2, 2), v[idx[6]](2, 2), 0.0f);
		break;
	case 8u:
		x.m[0][0] = _mm256_setr_ps(v[idx[0]](0, 0), v[idx[1]](0, 0), v[idx[2]](0, 0), v[idx[3]](0, 0), v[idx[4]](0, 0), v[idx[5]](0, 0), v[idx[6]](0, 0), v[idx[7]](0, 0));
		x.m[0][1] = _mm256_setr_ps(v[idx[0]](0, 1), v[idx[1]](0, 1), v[idx[2]](0, 1), v[idx[3]](0, 1), v[idx[4]](0, 1), v[idx[5]](0, 1), v[idx[6]](0, 1), v[idx[7]](0, 1));
		x.m[0][2] = _mm256_setr_ps(v[idx[0]](0, 2), v[idx[1]](0, 2), v[idx[2]](0, 2), v[idx[3]](0, 2), v[idx[4]](0, 2), v[idx[5]](0, 2), v[idx[6]](0, 2), v[idx[7]](0, 2));
		x.m[1][0] = _mm256_setr_ps(v[idx[0]](1, 0), v[idx[1]](1, 0), v[idx[2]](1, 0), v[idx[3]](1, 0), v[idx[4]](1, 0), v[idx[5]](1, 0), v[idx[6]](1, 0), v[idx[7]](1, 0));
		x.m[1][1] = _mm256_setr_ps(v[idx[0]](1, 1), v[idx[1]](1, 1), v[idx[2]](1, 1), v[idx[3]](1, 1), v[idx[4]](1, 1), v[idx[5]](1, 1), v[idx[6]](1, 1), v[idx[7]](1, 1));
		x.m[1][2] = _mm256_setr_ps(v[idx[0]](1, 2), v[idx[1]](1, 2), v[idx[2]](1, 2), v[idx[3]](1, 2), v[idx[4]](1, 2), v[idx[5]](1, 2), v[idx[6]](1, 2), v[idx[7]](1, 2));
		x.m[2][0] = _mm256_setr_ps(v[idx[0]](2, 0), v[idx[1]](2, 0), v[idx[2]](2, 0), v[idx[3]](2, 0), v[idx[4]](2, 0), v[idx[5]](2, 0), v[idx[6]](2, 0), v[idx[7]](2, 0));
		x.m[2][1] = _mm256_setr_ps(v[idx[0]](2, 1), v[idx[1]](2, 1), v[idx[2]](2, 1), v[idx[3]](2, 1), v[idx[4]](2, 1), v[idx[5]](2, 1), v[idx[6]](2, 1), v[idx[7]](2, 1));
		x.m[2][2] = _mm256_setr_ps(v[idx[0]](2, 2), v[idx[1]](2, 2), v[idx[2]](2, 2), v[idx[3]](2, 2), v[idx[4]](2, 2), v[idx[5]](2, 2), v[idx[6]](2, 2), v[idx[7]](2, 2));
	}
	return x;
}


inline void dyadicProduct(const Vector3f8 &a, const Vector3f8 &b, Matrix3f8 &res)
{
	res.m[0][0] = _mm256_mul_ps(a.x().v, b.x().v); res.m[0][1] = _mm256_mul_ps(a.x().v, b.y().v); res.m[0][2] = _mm256_mul_ps(a.x().v, b.z().v);
	res.m[1][0] = _mm256_mul_ps(a.y().v, b.x().v); res.m[1][1] = _mm256_mul_ps(a.y().v, b.y().v); res.m[1][2] = _mm256_mul_ps(a.y().v, b.z().v);
	res.m[2][0] = _mm256_mul_ps(a.z().v, b.x().v); res.m[2][1] = _mm256_mul_ps(a.z().v, b.y().v); res.m[2][2] = _mm256_mul_ps(a.z().v, b.z().v);
}

// ----------------------------------------------------------------------------------------------
//4 dimensional vector of Scalar8f to represent 8 quaternions
class Quaternion8f 
{
public:

	Scalarf8  q[4];

	inline Quaternion8f() { q[0] = 0.0; q[1] = 0.0; q[2] = 0.0; q[3] = 1.0; }

	inline Quaternion8f(Scalarf8 x, Scalarf8 y, Scalarf8 z, Scalarf8 w) {
		q[0] = x; q[1] = y; q[2] = z; q[3] = w;
	}

	inline Quaternion8f(Vector3f8& v) {
		q[0] = v[0]; q[1] = v[1]; q[2] = v[2]; q[3] = 0.0;
	}

	inline Scalarf8 & operator [] (int i) { return q[i]; }
	inline Scalarf8   operator [] (int i) const { return q[i]; }

	inline Scalarf8 & x() { return q[0]; }
	inline Scalarf8 & y() { return q[1]; }
	inline Scalarf8 & z() { return q[2]; }
	inline Scalarf8 & w() { return q[3]; }

	inline Scalarf8 x() const { return q[0]; }
	inline Scalarf8 y() const { return q[1]; }
	inline Scalarf8 z() const { return q[2]; }
	inline Scalarf8 w() const { return q[3]; }

	inline const Quaternion8f operator*(const Quaternion8f& a) const {
		return
			Quaternion8f(q[3] * a.q[0] + q[0] * a.q[3] + q[1] * a.q[2] - q[2] * a.q[1],
				q[3] * a.q[1] - q[0] * a.q[2] + q[1] * a.q[3] + q[2] * a.q[0],
				q[3] * a.q[2] + q[0] * a.q[1] - q[1] * a.q[0] + q[2] * a.q[3],
				q[3] * a.q[3] - q[0] * a.q[0] - q[1] * a.q[1] - q[2] * a.q[2]);
	}

	inline void toRotationMatrix(Matrix3f8& R)
	{
		const Scalarf8 tx = Scalarf8(2.0) * q[0];
		const Scalarf8 ty = Scalarf8(2.0) * q[1];
		const Scalarf8 tz = Scalarf8(2.0) * q[2];
		const Scalarf8 twx = tx*q[3];
		const Scalarf8 twy = ty*q[3];
		const Scalarf8 twz = tz*q[3];
		const Scalarf8 txx = tx*q[0];
		const Scalarf8 txy = ty*q[0];
		const Scalarf8 txz = tz*q[0];
		const Scalarf8 tyy = ty*q[1];
		const Scalarf8 tyz = tz*q[1];
		const Scalarf8 tzz = tz*q[2];

	    R.m[0][0] = Scalarf8(1.0) - (tyy + tzz);
		R.m[0][1] = txy - twz;
		R.m[0][2] = txz + twy;
		R.m[1][0] = txy + twz;
		R.m[1][1] = Scalarf8(1.0) - (txx + tzz);
		R.m[1][2] = tyz - twx;
		R.m[2][0] = txz - twy;
		R.m[2][1] = tyz + twx;
		R.m[2][2] = Scalarf8(1.0) - (txx + tyy);
	}

	inline void toRotationMatrix(Vector3f8& R1, Vector3f8& R2, Vector3f8& R3)
	{
		const Scalarf8 tx = Scalarf8(2.0) * q[0];
		const Scalarf8 ty = Scalarf8(2.0) * q[1];
		const Scalarf8 tz = Scalarf8(2.0) * q[2];
		const Scalarf8 twx = tx*q[3];
		const Scalarf8 twy = ty*q[3];
		const Scalarf8 twz = tz*q[3];
		const Scalarf8 txx = tx*q[0];
		const Scalarf8 txy = ty*q[0];
		const Scalarf8 txz = tz*q[0];
		const Scalarf8 tyy = ty*q[1];
		const Scalarf8 tyz = tz*q[1];
		const Scalarf8 tzz = tz*q[2];

		R1[0] = Scalarf8(1.0) - (tyy + tzz);
		R2[0] = txy - twz;
		R3[0] = txz + twy;
		R1[1] = txy + twz;
		R2[1] = Scalarf8(1.0) - (txx + tzz);
		R3[1] = tyz - twx;
		R1[2] = txz - twy;
		R2[2] = tyz + twx;
		R3[2] = Scalarf8(1.0) - (txx + tyy);
	}

	inline void store(std::vector<Quaternionr>& qf) const
	{
		float x[8], y[8], z[8], w[8];
		q[0].store(x);
		q[1].store(y);
		q[2].store(z);
		q[3].store(w);

		for (int i = 0; i < 8; i++)
		{
			qf[i].x() = x[i];
			qf[i].y() = y[i];
			qf[i].z() = z[i];
			qf[i].w() = w[i];
		}
	}

	inline void store(Quaternionr *qf) const
	{
		float x[8], y[8], z[8], w[8];
		q[0].store(x);
		q[1].store(y);
		q[2].store(z);
		q[3].store(w);

		for (int i = 0; i < 8; i++)
		{
			qf[i].x() = x[i];
			qf[i].y() = y[i];
			qf[i].z() = z[i];
			qf[i].w() = w[i];
		}
	}

	inline void set(const Quaternionr *qf)
	{
		float x[8], y[8], z[8], w[8];
		for(int i=0; i<8; i++)
		{
			x[i] = static_cast<float>(qf[i].x());
			y[i] = static_cast<float>(qf[i].y());
			z[i] = static_cast<float>(qf[i].z());
			w[i] = static_cast<float>(qf[i].w());
		}
		Scalarf8 s;
		s.load(x);
		q[0] = s;
		s.load(y);
		q[1] = s; 
		s.load(z);
		q[2] = s; 
		s.load(w);
		q[3] = s;
	}
};

// ----------------------------------------------------------------------------------------------
//alligned allocator so that vectorized types can be used in std containers
//from: https://stackoverflow.com/questions/8456236/how-is-a-vectors-data-aligned
template <typename T, std::size_t N = 32>
class AlignmentAllocator {
public:
	typedef T value_type;
	typedef std::size_t size_type;
	typedef std::ptrdiff_t difference_type;

	typedef T * pointer;
	typedef const T * const_pointer;

	typedef T & reference;
	typedef const T & const_reference;

public:
	inline AlignmentAllocator() throw () { }

	template <typename T2>
	inline AlignmentAllocator(const AlignmentAllocator<T2, N> &) throw () { }

	inline ~AlignmentAllocator() throw () { }

	inline pointer adress(reference r) {
		return &r;
	}

	inline const_pointer adress(const_reference r) const {
		return &r;
	}

	inline pointer allocate(size_type n) {
#ifdef _WIN32
		return (pointer)_aligned_malloc(n * sizeof(value_type), N);
#elif __linux__
		// NB! Argument order is opposite from MSVC/Windows
		return (pointer) aligned_alloc(N, n * sizeof(value_type));
#else
#error "Unknown platform"
#endif
	}

	inline void deallocate(pointer p, size_type) {
#ifdef _WIN32
		_aligned_free(p);
#elif __linux__
		free(p);
#else
#error "Unknown platform"
#endif
	}

	inline void construct(pointer p, const value_type & wert) {
		new (p) value_type(wert);
	}

	inline void destroy(pointer p) {
		p->~value_type();
	}

	inline size_type max_size() const throw () {
		return size_type(-1) / sizeof(value_type);
	}

	template <typename T2>
	struct rebind {
		typedef AlignmentAllocator<T2, N> other;
	};

	bool operator!=(const AlignmentAllocator<T, N>& other) const {
		return !(*this == other);
	}

	// Returns true if and only if storage allocated from *this
	// can be deallocated from other, and vice versa.
	// Always returns true for stateless allocators.
	bool operator==(const AlignmentAllocator<T, N>& other) const {
		return true;
	}
};

namespace SPH
{
	class MathFunctions_AVX
	{
	public:
		/** computes the APD of 8 deformation gradients. (Alg. 3 from the paper: Kugelstadt et al. "Fast Corotated FEM using Operator Splitting", CGF 2018)
		*/
		static void APD_Newton_AVX(const Vector3f8& F1, const Vector3f8& F2, const Vector3f8& F3, Quaternion8f& q)
		{
			//one iteration is sufficient for plausible results
			for (int it = 0; it < 1; it++)
			{
				//transform quaternion to rotation matrix
				Matrix3f8 R;
				q.toRotationMatrix(R);

				//columns of B = RT * F
				Vector3f8 B0 = R.transpose() * F1;
				Vector3f8 B1 = R.transpose() * F2;
				Vector3f8 B2 = R.transpose() * F3;

				Vector3f8 gradient(B2[1] - B1[2], B0[2] - B2[0], B1[0] - B0[1]);

				//compute Hessian, use the fact that it is symmetric
				Scalarf8 h00 = B1[1] + B2[2];
				Scalarf8 h11 = B0[0] + B2[2];
				Scalarf8 h22 = B0[0] + B1[1];
				Scalarf8 h01 = Scalarf8(-0.5) * (B1[0] + B0[1]);
				Scalarf8 h02 = Scalarf8(-0.5) * (B2[0] + B0[2]);
				Scalarf8 h12 = Scalarf8(-0.5) * (B2[1] + B1[2]);

				Scalarf8 detH = Scalarf8(-1.0) * h02 * h02 * h11 + Scalarf8(2.0) * h01 * h02 * h12 - h00 * h12 * h12 - h01 * h01 * h22 + h00 * h11 * h22;

				Vector3f8 omega;
				//compute symmetric inverse
				const Scalarf8 factor = Scalarf8(-0.25) / detH;
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

				omega = Vector3f8::blend(abs(detH) < 1.0e-9f, gradient * Scalarf8(-1.0), omega);	//if det(H) = 0 use gradient descent, never happened in our tests, could also be removed 

				//instead of clamping just use gradient descent. also works fine and does not require the norm
				//Scalarf8 useGD = blend(omega * gradient > Scalarf8(0.0), Scalarf8(1.0), Scalarf8(-1.0));
				//omega = Vector3f8::blend(useGD > Scalarf8(0.0), gradient * Scalarf8(-0.125), omega);
				omega = Vector3f8::blend(omega * gradient > Scalarf8(0.0), gradient * Scalarf8(-0.125), omega);

				Scalarf8 l_omega2 = omega.squaredNorm();
				const Scalarf8 w = (1.0 - l_omega2) / (1.0 + l_omega2);
				const Vector3f8 vec = omega * (2.0 / (1.0 + l_omega2));
				q = q * Quaternion8f(vec.x(), vec.y(), vec.z(), w);		//no normalization needed because the Cayley map returs a unit quaternion
			}
		}
	};
}

#endif
