/*
    microfacet.h -- Microfacet BSDF evaluation routines

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/common.h>
#include <vector>
#include <complex>

NAMESPACE_BEGIN(layer)

/**
 * \brief Smith's 1D shadowing masking term for the Beckmann microfacet
 * distribution
 *
 * \param v
 *    Incident direction
 * \param m
 *    Microsurface normal
 * \param alpha
 *    Beckmann roughness parameter
 */
extern Float smithG1(const Vector3 &v, const Vector3 &m, Float alpha);

/**
 * \brief Evaluate the Beckmann distribution-based microfacet BSDF by
 * Walter et al. using the mu_i, mu_o, phi_d parameterization
 *
 * \param mu_i
 *    Incident zenith angle cosine
 * \param mu_o
 *    Exitant zenith angle cosine
 * \param eta
 *    Relative index of refraction (complex)
 * \param alpha
 *    Beckmann roughness parameter
 * \param phi_d
 *    Azimuthal difference angle
 */
extern Float microfacet(Float mu_o, Float mu_i, std::complex<Float> eta,
                        Float alpha, Float phi_d);

/**
 * \brief Evaluate the Beckmann distribution-based microfacet BSDF by
 * Walter et al. using the mu_i, mu_o, phi_d parameterization. This
 * version leaves out the exponential term
 *
 * \param mu_i
 *    Incident zenith angle cosine
 * \param mu_o
 *    Exitant zenith angle cosine
 * \param eta
 *    Relative index of refraction (complex)
 * \param k
 *    Absorption coefficient
 * \param alpha
 *    Beckmann roughness parameter
 * \param phi_d
 *    Azimuthal difference angle
 */
extern Float microfacetNoExp(Float mu_o, Float mu_i, std::complex<Float> eta,
                             Float alpha, Float phi_d);

/**
 * \brief Return Fourier series coefficients for an exponential
 * of a cosine, specifically the expression "exp(A+B*cos(phi))"
 *
 * \param
 *    A The 'A' coefficient in above expression
 * \param A
 *    The 'B' coefficient in above expression
 * \param relerr
 *    Relative error goal
 */
extern void expCosFourierSeries(Float A, Float B, Float relerr,
                                std::vector<Float> &coeffs);

/**
 * \brief Compute a Fourier series of the Beckmann-distribution based
 * microfacet BSDF by Walter et al. (covers both the dielectric and
 * conducting case). This version leaves out the exponential term
 *
 * \param mu_i
 *    Incident zenith angle cosine
 * \param mu_o
 *    Exitant zenith angle cosine
 * \param eta
 *    Relative index of refraction (complex)
 * \param alpha
 *    Beckmann roughness parameter
 * \param n
 *    Indicates a desired number of Fourier coefficients (usually 8-12 are
 *    plenty -- in cases where the function contains high frequencies, it will
 *    automatically increase <tt>n</tt>).
 * \param phiMax
 *    The implementation minimizes the fitting error on the interval [0, phiMax].
 *    If in doubt, set phiMax = math::Pi
 * \param result
 *    Storage for the generated Fourier coefficients
 */
extern void microfacetNoExpFourierSeries(Float mu_o, Float mu_i,
                                         std::complex<Float> eta, Float alpha,
                                         int n, Float phiMax,
                                         std::vector<Float> &result);

/**
 * \brief Compute a Fourier series of the Beckmann-distribution based
 * microfacet BSDF by Walter et al. (covers both the dielectric and
 * conducting case)
 *
 * \param mu_i
 *    Incident zenith angle cosine
 * \param mu_o
 *    Exitant zenith angle cosine
 * \param eta
 *    Relative index of refraction (complex)
 * \param alpha
 *    Beckmann roughness parameter
 * \param n
 *    Indicates a desired maximum number of Fourier coefficients. The
 *    implementation will blur out higher Frequency content to try to
 *    achieve this number.
 * \param relerr
 *    A relative error threshold after which series terms can safely
 *    be truncated
 * \param result
 *    Storage for the generated Fourier coefficients
 */
extern void microfacetFourierSeries(Float mu_o, Float mu_i,
                                    std::complex<Float> eta, Float alpha, int n,
                                    Float relerr, std::vector<Float> &result);

NAMESPACE_END(layer)
