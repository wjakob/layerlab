/*
    hg.cpp -- Henyey-Greenstein model evaluation routines

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#include <layer/microfacet.h>
#include <layer/frame.h>
#include <layer/math.h>
#include <layer/fourier.h>
#include <layer/fresnel.h>
#include <Eigen/SVD>
#include <numeric>

NAMESPACE_BEGIN(layer)

static int expcosCoefficientCount(Float B, Float relerr) {
    Float prod = 1, invB = 1.0f / B;
    if (B == 0)
        return 1;

    for (int i=0; ; ++i) {
        prod /= 1 + i * invB;

        if (prod < relerr)
            return i+1;
    }
}

static Float modBesselRatio(Float B, Float k) {
    const Float eps = std::numeric_limits<Float>::epsilon(),
                invTwoB = 2.0f / B;

    Float i = (Float) k,
           D = 1 / (invTwoB * i++),
           Cd = D, C = Cd;

    while (std::abs(Cd) > eps * std::abs(C)) {
        Float coeff = invTwoB * i++;
        D = 1 / (D + coeff);
        Cd *= coeff*D - 1;
        C += Cd;
    }

    return C;
}

void expCosFourierSeries(Float A, Float B, Float relerr, std::vector<Float> &coeffs) {
    /* Determine the required number of coefficients and allocate memory */
    int n = expcosCoefficientCount(B, relerr);
    coeffs.resize(n);

    /* Determine the last ratio and work downwards */
    coeffs[n-1] = modBesselRatio(B, n - 1);
    for (int i=n-2; i>0; --i)
        coeffs[i] = B / (2*i + B*coeffs[i+1]);

    /* Evaluate the exponentially scaled I0 and correct scaling */
    coeffs[0] = math::i0e(B) * std::exp(A+B);

    /* Apply the ratios & factor of two upwards */
    Float prod = 2*coeffs[0];
    for (int i=1; i<n; ++i) {
        prod *= coeffs[i];
        if (std::abs(prod) < coeffs[0] * relerr) {
            coeffs.erase(coeffs.begin() + i, coeffs.end());
            break;
        }
        coeffs[i] = prod;
    }
}

Float smithG1(const Vector3 &v, const Vector3 &m, Float alpha) {
    const Float tanTheta = std::abs(Frame::tanTheta(v));

    /* Can't see the back side from the front and vice versa */
    if (v.dot(m) * Frame::cosTheta(v) <= 0)
        return 0.0f;

    Float a = 1.0f / (alpha * tanTheta);
    if (a < 1.6f) {
        /* Use a fast and accurate (<0.35% rel. error) rational
           approximation to the shadowing-masking function */
        const Float aSqr = a * a;
        return (3.535f * a + 2.181f * aSqr)
             / (1.0f + 2.276f * a + 2.577f * aSqr);
    }

    return 1.f;
}

Float microfacet(Float mu_o, Float mu_i, std::complex<Float> eta_,
                 Float alpha, Float phi_d) {
    Float sinThetaI = math::safe_sqrt(1-mu_i*mu_i),
          sinThetaO = math::safe_sqrt(1-mu_o*mu_o),
          cosPhi = std::cos(phi_d),
          sinPhi = std::sin(phi_d);

    Vector wi(-sinThetaI, 0, -mu_i);
    Vector wo(sinThetaO*cosPhi, sinThetaO*sinPhi, mu_o);
    bool reflect = -mu_i*mu_o > 0;

    if (mu_o == 0 || mu_i == 0)
        return 0.f;
    
    bool conductor = eta_.imag() != 0.0f;
    if (conductor && !reflect)
        return 0.0f;
    std::complex<Float> eta =
        (-mu_i > 0 || conductor) ? eta_ : std::complex<Float>(1) / eta_;

    Vector H = (wi + wo * (reflect ? 1.0f : eta.real())).normalized();
    H *= math::signum(Frame::cosTheta(H));

    Float cosThetaH2 = Frame::cosTheta2(H),
          exponent = -Frame::tanTheta2(H) / (alpha*alpha),
          D = std::exp(exponent) / (math::Pi * alpha*alpha * cosThetaH2*cosThetaH2),
          F = !conductor ? fresnelDielectric(wi.dot(H), eta_.real())
                         : fresnelConductor(std::abs(wi.dot(H)), eta),
          G = smithG1(wi, H, alpha) * smithG1(wo, H, alpha);

    if (reflect) {
        return F * D * G / (4.0f * std::abs(mu_i*mu_o));
    } else {
        Float sqrtDenom = wi.dot(H) + eta.real() * wo.dot(H);

        return std::abs(((1 - F) * D * G * eta.real() * eta.real() * wi.dot(H)
            * wo.dot(H)) / (mu_i*mu_o * sqrtDenom * sqrtDenom));
    }
}

Float microfacetNoExp(Float mu_o, Float mu_i, std::complex<Float> eta_,
                      Float alpha, Float phi_d) {
    Float sinThetaI = math::safe_sqrt(1-mu_i*mu_i),
          sinThetaO = math::safe_sqrt(1-mu_o*mu_o),
          cosPhi = std::cos(phi_d),
          sinPhi = std::sin(phi_d);

    Vector wi(-sinThetaI, 0, -mu_i);
    Vector wo(sinThetaO*cosPhi, sinThetaO*sinPhi, mu_o);
    bool reflect = -mu_i*mu_o > 0;

    if (mu_o == 0 || mu_i == 0)
        return 0.f;

    bool conductor = eta_.imag() != 0.0f;
    if (conductor && !reflect)
        return 0.0f;
    std::complex<Float> eta =
        (-mu_i > 0 || conductor) ? eta_ : std::complex<Float>(1) / eta_;

    Vector H = (wi + wo * (reflect ? 1.0f : eta.real())).normalized();
    H *= math::signum(Frame::cosTheta(H));

    Float cosThetaH2 = Frame::cosTheta2(H),
          D = (Float) 1 / (math::Pi * alpha*alpha * cosThetaH2*cosThetaH2),
          F = !conductor ? fresnelDielectric(wi.dot(H), eta_.real())
                         : fresnelConductor(std::abs(wi.dot(H)), eta),
          G = smithG1(wi, H, alpha) * smithG1(wo, H, alpha);

    if (reflect) {
        return F * D * G / (4.0f * std::abs(mu_i*mu_o));
    } else {
        Float sqrtDenom = wi.dot(H) + eta.real() * wo.dot(H);

        return std::abs(((1 - F) * D * G * eta.real() * eta.real() * wi.dot(H)
            * wo.dot(H)) / (mu_i*mu_o * sqrtDenom * sqrtDenom));
    }
}

static Float Bmax(size_t n, Float relerr) {
    if (relerr >= 1e-1f)
        return 0.1662f*std::pow((Float) n, (Float) 2.05039);
    else if (relerr >= 1e-2f)
        return 0.0818f*std::pow((Float) n, (Float) 2.04982);
    else if (relerr >= 1e-3f)
        return 0.0538f*std::pow((Float) n, (Float) 2.05001);
    else if (relerr >= 1e-4f)
        return 0.0406f*std::pow((Float) n, (Float) 2.04686);
    else if (relerr >= 1e-5f)
        return 0.0337f*std::pow((Float) n, (Float) 2.03865);
    else if (relerr >= 1e-6f)
        return 0.0299f*std::pow((Float) n, (Float) 2.02628);
    else
        throw std::runtime_error("Bmax(): unknown relative error bound!");
}

void microfacetNoExpFourierSeries(Float mu_o, Float mu_i, std::complex<Float> eta_,
                                  Float alpha, size_t n, Float phiMax,
                                  std::vector<Float> &result) {

    bool reflect = -mu_i * mu_o > 0;

    Float sinMu2 = math::safe_sqrt((1 - mu_i * mu_i) * (1 - mu_o * mu_o)),
          phiCritical = 0.0f;

    bool conductor = (eta_.imag() != 0.0f);
    std::complex<Float> eta =
        (-mu_i > 0 || conductor) ? eta_ : std::complex<Float>(1) / eta_;

    if (reflect) {
        if (!conductor)
            phiCritical = math::safe_acos((2*eta.real()*eta.real()-mu_i*mu_o-1)/sinMu2);
    } else if (!reflect) {
        if (conductor)
            throw std::runtime_error("lowfreqFourierSeries(): encountered refraction case for a conductor");
        Float etaDenser = (eta.real() > 1 ? eta.real() : 1 / eta.real());
        phiCritical = math::safe_acos((1 - etaDenser * mu_i * mu_o) /
                                      (etaDenser * sinMu2));
    }

    if (!conductor && phiCritical > math::Epsilon &&
        phiCritical < math::Pi - math::Epsilon &&
        phiCritical < phiMax - math::Epsilon) {
        /* Uh oh, some high frequency content leaked in the generally low frequency part.
           Increase the number of coefficients so that we can capture it. Fortunately, this
           happens very rarely. */
        n = std::max(n, (size_t) 100);
    }

    VectorX coeffs(n);
    coeffs.setZero();
    std::function<Float(Float)> integrand = std::bind(
        &microfacetNoExp, mu_o, mu_i, eta_, alpha, std::placeholders::_1);

    const int nEvals = 200;
    if (reflect) {
        if (phiCritical > math::Epsilon && phiCritical < phiMax-math::Epsilon) {
            filonIntegrate(integrand, coeffs.data(), n, nEvals, 0, phiCritical);
            filonIntegrate(integrand, coeffs.data(), n, nEvals, phiCritical, phiMax);
        } else {
            filonIntegrate(integrand, coeffs.data(), n, nEvals, 0, phiMax);
        }
    } else {
        filonIntegrate(integrand, coeffs.data(), n, nEvals, 0,
                       std::min(phiCritical, phiMax));
    }

    if (phiMax < math::Pi - math::Epsilon) {
        /* Precompute some sines and cosines */
        VectorX cosPhi(n), sinPhi(n);
        for (size_t i=0; i<n; ++i) {
            sinPhi[i] = std::sin(i*phiMax);
            cosPhi[i] = std::cos(i*phiMax);
        }

        /* The fit only occurs on a subset [0, phiMax], where the Fourier
           Fourier basis functions are not orthogonal anymore! The following
           then does a change of basis to proper Fourier coefficients. */
        MatrixX A(n, n);

        for (MatrixX::Index i=0; i < (MatrixX::Index) n; ++i) {
            for (MatrixX::Index j=0; j <= (MatrixX::Index) i; ++j) {
                if (i != j) {
                    A(i, j) = A(j, i) = (i * cosPhi[j] * sinPhi[i] -
                                         j * cosPhi[i] * sinPhi[j]) /
                                        (i * i - j * j);
                } else if (i != 0) {
                    A(i, i) = (std::sin(2 * i * phiMax) + 2 * i * phiMax) / (4 * i);
                } else {
                    A(i, i) = phiMax;
                }
            }
        }

        auto svd = A.bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
        const MatrixX &U = svd.matrixU();
        const MatrixX &V = svd.matrixV();
        const VectorX &sigma = svd.singularValues();

        if (sigma[0] == 0) {
            result.clear();
            result.push_back(0);
            return;
        }

        VectorX temp = VectorX::Zero(n);
        coeffs[0] *= math::Pi;
        coeffs.tail(n-1) *= 0.5 * math::Pi;
        for (size_t i=0; i<n; ++i) {
            if (sigma[i] < 1e-9f * sigma[0])
                break;
            temp += V.col(i) * U.col(i).dot(coeffs) / sigma[i];
        }
        coeffs = temp;
    }

    result.resize(coeffs.size());
    memcpy(result.data(), coeffs.data(), sizeof(Float) * coeffs.size());
}

void microfacetFourierSeries(Float mu_o, Float mu_i, std::complex<Float> eta_,
                             Float alpha, size_t n, Float relerr,
                             std::vector<Float> &result) {
    bool reflect = -mu_i * mu_o > 0;

    /* Compute the 'A' and 'B' constants, as well as the critical azimuth */
    Float A, B;
    Float sinMu2 = math::safe_sqrt((1 - mu_i * mu_i) * (1 - mu_o * mu_o));

    bool conductor = (eta_.imag() != 0.0f);
    std::complex<Float> eta =
        (-mu_i > 0 || conductor) ? eta_ : std::complex<Float>(1) / eta_;

    if (reflect) {
        Float temp = 1.0f / (alpha * (mu_i - mu_o));
        A = (mu_i * mu_i + mu_o * mu_o - 2) * temp * temp;
        B = 2 * sinMu2 * temp * temp;
    } else {
        if (conductor) {
            /* No refraction in conductors */
            result.clear();
            result.push_back(0.0f);
            return;
        } else {
            Float temp = 1.0f / (alpha * (mu_i - eta.real() * mu_o));
            A = (mu_i * mu_i - 1 + eta.real() * eta.real() * (mu_o * mu_o - 1)) * temp * temp;
            B = 2 * eta.real() * sinMu2 * temp * temp;
        }
    }

    /* Minor optimization: don't even bother computing the Fourier series
       if the contribution to the scattering model is miniscule */
    if (math::i0e(B) * std::exp(A+B) < 1e-10) {
        result.clear();
        result.push_back(0.0f);
        return;
    }

    Float B_max = Bmax(n, relerr);
    if (B > B_max) {
        A = A + B - B_max + std::log(math::i0e(B) / math::i0e(B_max));
        B = B_max;
    }

    std::vector<Float> lowfreq_coeffs, expcos_coeffs;

    /* Compute Fourier coefficients of the exponential term */
    expCosFourierSeries(A, B, relerr, expcos_coeffs);

    /* Compute Fourier coefficients of the low-frequency term
       Only fit in the region where the result is actually
       going to make some sort of difference given the convolution
       with expcos_coeffs */
    Float phiMax = math::safe_acos(1 + std::log(relerr) / B);

    microfacetNoExpFourierSeries(mu_o, mu_i, eta_, alpha,
                                 12, phiMax, lowfreq_coeffs);

    /* Perform discrete circular convolution of the two series */
    result.resize(lowfreq_coeffs.size() + expcos_coeffs.size() - 1);

    convolveFourier(lowfreq_coeffs.data(), lowfreq_coeffs.size(),
        expcos_coeffs.data(), expcos_coeffs.size(), result.data());

    /* Truncate the series if error bounds are satisfied */
    for (size_t i=0; i<result.size(); ++i) {
        assert(std::isfinite(result[i]));
        if (result[i] == 0 || std::abs(result[i]) < result[0] * relerr) {
            result.erase(result.begin() + i, result.end());
            break;
        }
    }
}

NAMESPACE_END(layer)
