/*
    hg.h -- Henyey-Greenstein model evaluation routines

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/math.h>
#include <vector>

NAMESPACE_BEGIN(layer)

/**
 * \brief Evaluate the HG model using the mu_i, mu_o, phi_d parameterization
 *
 * \param mu_i
 *    Incident zenith angle cosine
 * \param mu_o
 *    Exitant zenith angle cosine
 * \param phi_d
 *    Azimuthal difference angle
 * \param g
 *    Anisotropy parameter
 */
Float hg(Float mu_o, Float mu_i, Float g, Float phi_d);

/**
 * \brief Compute a Fourier series of the HG model
 *
 * This function first finds the 0-th and 1st-order Fourier coefficients using
 * elliptic integrals.
 *
 * The other coefficients can then be determined much more efficiently; the
 * approach here is based on the idea that the ratio of adjacent coefficients
 * eventually converges to a constant value. Using a 2nd-order Taylor
 * expansion, we can obtain a fairly accurate estimate of this ratio somewhere
 * "in the middle" (i.e. for large $n$, but well before the aforementioned
 * convergence).
 *
 * Using a backwards recurrence scheme, we can then determine all previous
 * ratios and, thereby (using the explicitly computed first values), all
 * Fourier coefficients.
 *
 * This approach is based on the article
 *
 * "A Recurrence Formula For Computing Fourier Components of the
 *  Henyey-Greenstein Phase Function" by E.G. Yanovitskij
 *
 * Journal of Quantitative Spectroscopy & Radiative Transfer, 57, no 1. 1977
 *
 * \param mu_i
 *    Incident zenith angle cosine
 * \param mu_o
 *    Exitant zenith angle cosine
 * \param g
 *    Anisotropy parameter
 * \param kmax
 *    Indicates a desired maximum number of Fourier coefficients. The
 *    implementation will blur out higher Frequency content to try to
 *    achieve this number.
 * \param result
 *    Storage for the generated Fourier coefficients
 */
void hgFourierSeries(Float mu_o, Float mu_i, Float g, int kmax,
					 Float relerr, std::vector<Float> &result);

NAMESPACE_END(layer)
