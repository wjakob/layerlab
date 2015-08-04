/*
    fourier.h -- Functions for sampling and evaluating Fourier series

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/color.h>
#include <layer/math.h>

#if defined(__AVX__)
#define FOURIER_SCALAR     0
#define FOURIER_VECTORIZED 1
#define LANE_WIDTH         8
#else
#define FOURIER_SCALAR     1
#define FOURIER_VECTORIZED 0
#define LANE_WIDTH         1
#endif

#define LANE_SIZE_BYTES (LANE_WIDTH * 4)

NAMESPACE_BEGIN(layer)

/**
 * \brief Evaluate an even Fourier series (i.e. containing
 * only cosine terms).
 *
 * \param coeffs
 *     Coefficient storage
 * \param nCoeffs
 *     Denotes the size of \c coeffs
 * \param phi
 *     Angle for which the series should be evaluated
 * \return
 *     The value of the series for this angle
 * \remark
 *    In the Python API, the \c nCoeffs parameter is automatically
 *    computed from the length of the input arrays and thus not needed.
 */
Float evalFourier(const float *coeffs, size_t nCoeffs, Float phi);

/**
 * \brief Simultaneously evaluate <em>three</em> even Fourier series
 * corresponding to the color channels (Y, R, B) and return
 * a spectral power distribution
 *
 * \param coeffs
 *     Coefficient storage
 * \param nCoeffs
 *     Denotes the size of \c coeffs
 * \param phi
 *     Angle for which the series should be evaluated
 * \return
 *     The value of the series for this angle
 * \remark
 *    In the Python API, the \c nCoeffs parameter is automatically
 *    computed from the length of the input arrays and thus not needed.
 */
Color3 evalFourier3(float *const coeffs[3], size_t nCoeffs, Float phi);

/**
 * \brief Sample a angle from an even Fourier series
 * (i.e. containing only cosine terms).
 *
 * This is done by importance sampling a uniform piecewise constant
 * approximation with 2^nRecursions elements, which is constructed
 * on the fly in a recursive manner.
 *
 * \param coeffs
 *     Coefficient storage
 * \param nCoeffs
 *     Denotes the size of \c coeffs
 * \param recip
 *     Precomputed array of integer reciprocals, i.e. <tt>inf, 1, 1/2, 1/3,
 *     1/4, ..</tt> of size <tt>nCoeffs-1</tt>. This is used to reduce
 *     integer-to-FP pipeline stalls and division latencies at runtime.
 * \param sample
 *     A uniformly distributed random sample in the interval <tt>[0,1]</tt>
 * \param[out] pdf
 *     This parameter is used to return the probability density of
 *     the sampling scheme evaluated at the generated point on the
 *     underlying piecewise constant approximation (on [0, \pi])
 * \param[out] phi
 *     Used to return the sampled angle (on [0, \pi])
 * \return
 *     The importance weight (i.e. the value of the Fourier series
 *     divided by \c pdf)
 */
Float sampleFourier(const float *coeffs, const Float *recip, size_t nCoeffs,
                    Float sample, Float &pdf, Float &phi);

/**
 * \brief Sample a angle from <em>three</em> Fourier series
 * corresponding to the color channels (Y, R, B)
 *
 * This is done by importance sampling a uniform piecewise constant
 * approximation with respect to the luminance, which is constructed
 * on the fly in a recursive manner.
 *
 * \param coeffs
 *     Coefficient storage
 * \param nCoeffs
 *     Denotes the size of \c coeffs
 * \param recip
 *     Precomputed array of integer reciprocals, i.e. <tt>inf, 1, 1/2, 1/3,
 *     1/4, ..</tt> of size <tt>nCoeffs-1</tt>. This is used to reduce
 *     integer-to-FP pipeline stalls and division latencies at runtime.
 * \param sample
 *     A uniformly distributed random sample in the interval <tt>[0,1]</tt>
 * \param[out] pdf
 *     This parameter is used to return the probability density of
 *     the sampling scheme evaluated at the generated point on the
 *     underlying piecewise constant approximation (on [0, \pi])
 * \param[out] phi
 *     Used to return the sampled angle (on [0, \pi])
 * \return
 *     The importance weight (i.e. the value of the Fourier series
 *     divided by \c pdf)
 */
Color3 sampleFourier3(float *const coeffs[3], const double *recip, int nCoeffs,
                      Float sample, Float &pdf, Float &phi);

/**
 * \brief Evaluate the probability density of the sampling scheme
 * implemented by \ref sampleFourier() at the position \c phi.
 *
 * \param coeffs
 *     Coefficient storage
 * \param nCoeffs
 *     Denotes the size of \c coeffs
 * \return
 *     The continuous probability density on [0, \pi]
 */
inline Float pdfFourier(const float *coeffs, int nCoeffs, Float phi) {
    return evalFourier(coeffs, nCoeffs, phi) * math::InvTwoPi / (Float)coeffs[0];
}

/**
 * \brief Computes the Fourier series of a product of even Fourier series
 * using discrete convolution.
 *
 * The input series are assumed to have \c ka and \c kb coefficients, and the
 * output must have room for <tt>ka+kb-1</tt> coefficients.
 *
 * \remark a
 *    First input array of Fourier coefficients
 * \remark ka
 *    Size of the first input array
 * \remark b
 *    Second input array of Fourier coefficients
 * \remark kb
 *    Size of the second input array
 * \remark c
 *    Pointer into the output array
 * \remark
 *    In the Python API, the \c ka and \c kb parameters are automatically
 *    computed from the lengths of the input arrays and thus not needed.
 *    The \c parameter is also removed (the result array is directly returned).
 */
void convolveFourier(const Float *a, size_t ka, const Float *b, size_t kb, Float *c);

/**
 * \brief Compute a Fourier series of the given even function by integrating it
 * against the basis functions using Filon quadrature
 *
 * Filon quadrature works by constructing a piecewise quadratic interpolant of
 * the original function. The Fourier basis functions are then integrated
 * against this representation, which has an analytic solution. This avoids all
 * of the problems of traditional quadrature schemes involving highly
 * oscillatory integrals. It is thus possible to obtain accurate coefficients
 * even for high orders.
 *
 * \param f
 *    Function to be integrated
 * \param[out] coeffs
 *    Output buffer used to store the computed coefficients. The function adds
 *    the computed coefficients to the buffer rather than overwriting the
 *    existing contents.
 * \param nCoeffs
 *    Desired number of coefficients
 * \param nEvals
 *    Desired resolution of the piecewise quadratic interpolant
 * \param a
 *    Start of the integration, can optionally be set to values other than
 *    zero. Note that the Fourier basis functions are not orthogonal anymore in
 *    this case.
 * \param b
 *    End of the integration, can be set to values other than pi. Note that the
 *    Fourier basis functions are not orthogonal anymore in this case.
 * \remark
 *    In the Python API, the \c coeffs array is directly returned.
 */
void filonIntegrate(const std::function<Float(Float)> &f, Float *coeffs,
                    size_t nCoeffs, int nEvals, Float a = 0, Float b = math::Pi);

#if FOURIER_SCALAR == 1

#define fourier_aligned_alloca(size) \
    __align_helper(alloca(size), size)

static inline float *__align_helper(void *ptr, size_t size) {
    memset(ptr, 0, size);
    return (float *) ptr;
}

#else

/**
 * \brief Helper functions for allocating temporary stack memory for
 * Fourier coefficients.
 *
 * This macro works like alloca(), except that it also ensures that
 * the resulting buffer is properly aligned so that it can be used
 * with the SSE/AVX-vectorized evaluation and sampling routines.
 *
 * Given an allocation, we require memory for
 *
 *  1 float at alignment -4
 *  N quadruplets at alignment 0
 *
 * SSE: 12 bytes in addition to make sure that
 *      this alignment can be established
 * AVX: 28 bytes in addition to make sure that
 *      this alignment can be established
 *
 * i.e.
 *
 * SSE:  size = 4 + 16*((size-4 + 12) / 16) + 12;
 * SSE:  size = 4 + 32*((size-4 + 28) / 32) + 28;
 *
 * and to align:
 *
 * SSE:  buffer += 12 - buffer mod 16
 * AVX:  buffer += 28 - buffer mod 32
 */

#define fourier_aligned_alloca(size) \
    __align_helper(alloca(LANE_SIZE_BYTES*(((LANE_SIZE_BYTES-2*sizeof(float))+size) / LANE_SIZE_BYTES) + LANE_SIZE_BYTES), \
                          LANE_SIZE_BYTES*(((LANE_SIZE_BYTES-2*sizeof(float))+size) / LANE_SIZE_BYTES) + LANE_SIZE_BYTES)

namespace {
    static inline float *__align_helper(void *ptr, size_t size) {
        memset(ptr, 0, size);
        return (float *) ((uint8_t *) ptr + (LANE_SIZE_BYTES - sizeof(float) - ((uintptr_t) ptr) % LANE_SIZE_BYTES));
    }
};

#endif

NAMESPACE_END(layer)
