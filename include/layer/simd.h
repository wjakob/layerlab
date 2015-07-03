/*
    simd.h -- Useful declarations for SIMD optimized code

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/common.h>

#if defined(_MSC_VER)
# include <intrin.h>
#else
# include <immintrin.h>
#endif

#define _mm256_set_ss(value)    _mm256_insertf128_ps(_mm256_setzero_ps(), _mm_set_ss(value), 0)
#define _mm256_set_sd(value)    _mm256_insertf128_pd(_mm256_setzero_pd(), _mm_set_sd(value), 0)
#define _mm256_splat0_pd(value) _mm256_permute2f128_pd(_mm256_shuffle_pd(value, value, 0x0), value, 0x00)
#define _mm256_splat1_pd(value) _mm256_permute2f128_pd(_mm256_shuffle_pd(value, value, 0x3), value, 0x00)
#define _mm256_splat2_pd(value) _mm256_permute2f128_pd(_mm256_shuffle_pd(value, value, 0x0), value, 0x11)
#define _mm256_splat3_pd(value) _mm256_permute2f128_pd(_mm256_shuffle_pd(value, value, 0xC), value, 0x11)

#if !defined(MM_ALIGN16)
# if defined(__GNUC__)
#  define MM_ALIGN16 __attribute__ ((aligned (16)))
# else
#  define MM_ALIGN16 __declspec(align(16))
# endif
#endif

#if !defined(MM_ALIGN32)
# if defined(__GNUC__)
#  define MM_ALIGN32 __attribute__ ((aligned (32)))
# else
#  define MM_ALIGN32 __declspec(align(32))
# endif
#endif

NAMESPACE_BEGIN(layer)
NAMESPACE_BEGIN(simd)

/// Allocate an aligned region of memory
void * malloc(size_t size);

/// Free an aligned region of memory
void free(void *ptr);

#if defined(__AVX__)
/// Perform four simultaneous double precision horizontal additions using AVX
inline void hadd(__m256d a, __m256d b, __m256d c, __m256d d, double *target) {
    /* See http://stackoverflow.com/questions/10833234/4-horizontal-double-precision-sums-in-one-go-with-avx */
    __m256d sumab = _mm256_hadd_pd(a, b);
    __m256d sumcd = _mm256_hadd_pd(c, d);
    __m256d blend = _mm256_blend_pd(sumab, sumcd, 0x0C);
    __m256d perm = _mm256_permute2f128_pd(sumab, sumcd, 0x21);
    __m256d result =  _mm256_add_pd(perm, blend);
    _mm256_store_pd(target, result);
}

/// Perform two simultaneous double precision horizontal additions using AVX
inline void hadd(__m256d a, __m256d b, double *target) {
    /* See http://stackoverflow.com/questions/9775538/fastest-way-to-do-horizontal-vector-sum-with-avx-instructions */
    __m256d sum = _mm256_hadd_pd(a, b);
    __m128d sum_high = _mm256_extractf128_pd(sum, 1);
    __m128d result = _mm_add_pd(sum_high, _mm256_castpd256_pd128(sum));
    _mm_store_pd(target, result);
}

/// Perform a double precision horizontal addition using AVX
inline double hadd(__m256d x) {
    double MM_ALIGN16 result[2];
    hadd(x, x, result);
    return result[0];
}
#endif

NAMESPACE_END(simd)
NAMESPACE_END(layer)
