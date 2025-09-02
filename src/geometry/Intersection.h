#pragma once

#include "intrinsic/intrinsics.h"
#include "geometry_lib.h"

#define _div _mm_div_ps
#define _mul _mm_mul_ps
#define _add _mm_add_ps
#define _sub _mm_sub_ps
#define _dp _mm_dp_ps /// dot product
#define _abs(p) _mm_andnot_ps(_mm_set1_ps(-0.f), p)

/// set __m128 vector as SSE::Vector structure
#define _set(v) _mm_setr_ps(v.coordinates[0], v.coordinates[1], v.coordinates[2], 0.0)

/// cross product
#define _xp(a, b) _sub(_mul(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,0,2,1)),\
    _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,1,0,2))),\
    _mul(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,1,0,2)),\
    _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,0,2,1))))

/// length of the __m128 vector
#define _len(v) _mm_cvtss_f32(_mm_sqrt_ss(_dp(v, v, 0x71)));


#define EPS_PROJECTION		0.0000174532836589830883577820272085
//#define EPS_PROJECTION		0.00174532836589830883577820272085
#define EPS_LAYONLINE		0.05

const float EPS_INTERSECT = 0.08;
//const float EPS_INTERSECT = 0.0008;
//const float EPS_MERGE = 0.08;
const float EPS_MERGE = 0.001;
const float EPS_INSIDE = -0.06;

inline bool is_inside_i(__m128 x, __m128 p1, __m128 p2, __m128 normal)
{
	__m128 m_eps = _mm_set_ss(EPS_INSIDE);

	__m128 p1_p2 = _mm_sub_ps(p2, p1);
	__m128 p1_x = _mm_sub_ps(x, p1);
	__m128 dir = _mm_dp_ps(_cross_product(p1_p2, p1_x), normal, MASK_1LOW);

	return _mm_ucomigt_ss(dir, m_eps);
}

inline bool is_layOnLine_i(__m128 _x, __m128 _a, __m128 _b)
{
	__m128 ab = _mm_sub_ps(_a, _b);
	__m128 ax = _mm_sub_ps(_a, _x);
	__m128 bx = _mm_sub_ps(_b, _x);

	__m128 sqr_len_ab = _mm_dp_ps(ab, ab, MASK_1LOW);
	__m128 sqr_len_ax = _mm_dp_ps(ax, ax, MASK_1LOW);
	__m128 sqr_len_bx = _mm_dp_ps(bx, bx, MASK_1LOW);

	return (sqr_len_ax[0] + sqr_len_bx[0] < sqr_len_ab[0] + sqr_len_ab[0]*EPS_LAYONLINE);
}

// OPT: try to return 'ok' insted 'x'
inline __m128 intersect_i(__m128 _a1, __m128 _a2, __m128 _b1, __m128 _b2,
						  __m128 _normal_to_facet, bool &ok)
{
	__m128 _v_a = _mm_sub_ps(_a2, _a1);
	__m128 _v_b = _mm_sub_ps(_b2, _b1);

	// normal of new plane
	__m128 _normal_to_line = _cross_product(_v_b, _normal_to_facet);

	// normalize normal
	__m128 _normal_n = _normalize(_normal_to_line);

	// intersection vector and new plane
	__m128 _dp0 = _mm_dp_ps(_v_a, _normal_n, MASK_FULL);

	__m128 _sign_mask = _mm_set1_ps(-0.f);
	__m128 _abs_dp = _mm_andnot_ps(_sign_mask, _dp0);

	if (_abs_dp[0] < EPS_INTERSECT)
	{
		ok = false;
		return _dp0;
	}

	__m128 _dp1 = _mm_dp_ps(_a1, _normal_n, MASK_FULL);
	__m128 m_d_param = _mm_dp_ps(_b1, _normal_n, MASK_FULL);

    __m128 _a = _mm_sub_ps(_dp1, m_d_param);
    __m128 _t = _mm_div_ps(_a, _dp0);

    __m128 _m = _mm_mul_ps(_t, _v_a);

	ok = true;
    return _mm_sub_ps(_a1, _m);
}

/**
 * @brief Intersects two vectors laid on the same plane
 * @param _a1 point in first vector
 * @param _b1 point in second vector
 * @param _v_a first vector
 * @param _v_b second vector
 * @param _normal_to_facet normal to plane
 * @param ok true if vectors are not parallel
 * @return intersection point
 */
inline __m128 intersect_iv(__m128 _a1, __m128 _b1, __m128 _v_a, __m128 _v_b,
						   __m128 _normal_to_facet, bool &ok)
{
	// normal of new plane // OPT: try to do buffer for other variables from this
	__m128 _normal_to_line = _cross_product(_v_b, _normal_to_facet);
	__m128 _normal_n = _normalize(_normal_to_line);

	// intersection vector and new plane
	__m128 _dp0 = _mm_dp_ps(_v_a, _normal_n, MASK_FULL);

	__m128 _sign_mask = _mm_set1_ps(-0.f);
	__m128 _abs_dp = _mm_andnot_ps(_sign_mask, _dp0);

	if (_abs_dp[0] < EPS_INTERSECT)
	{
		ok = false;
		return _dp0;
	}

	__m128 _dp1 = _mm_dp_ps(_a1, _normal_n, MASK_FULL);
	__m128 m_d_param = _mm_dp_ps(_b1, _normal_n, MASK_FULL);

    __m128 _a = _mm_sub_ps(_dp1, m_d_param);
    __m128 _t = _mm_div_ps(_a, _dp0);

    __m128 _m = _mm_mul_ps(_t, _v_a);

	ok = true;
    return _mm_sub_ps(_a1, _m);
}

void computeIntersection(const Point3f &s, const Point3f &e,
						 const Point3f &p1, const Point3f &p2, const Point3f &normal,
						 Point3f &x);
