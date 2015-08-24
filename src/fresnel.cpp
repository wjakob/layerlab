/*
    fresnel.cpp -- Fresnel coefficients for dielectrics and conductors

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#include <layer/fresnel.h>
#include <layer/math.h>
#include <iostream>

NAMESPACE_BEGIN(layer)

Float fresnelDielectric(Float cosThetaI_, Float &cosThetaT_, Float eta) {
    if (eta == 1) {
        cosThetaT_ = -cosThetaI_;
        return 0.0f;
    }

    /* Using Snell's law, calculate the squared sine of the
       angle between the normal and the transmitted ray */
    Float scale = (cosThetaI_ > 0) ? 1/eta : eta,
          cosThetaTSqr = 1 - (1-cosThetaI_*cosThetaI_) * (scale*scale);

    /* Check for total internal reflection */
    if (cosThetaTSqr <= 0.0f) {
        cosThetaT_ = 0.0f;
        return 1.0f;
    }

    /* Find the absolute cosines of the incident/transmitted rays */
    Float cosThetaI = std::abs(cosThetaI_);
    Float cosThetaT = std::sqrt(cosThetaTSqr);

    Float Rs = (cosThetaI - eta * cosThetaT)
             / (cosThetaI + eta * cosThetaT);
    Float Rp = (eta * cosThetaI - cosThetaT)
             / (eta * cosThetaI + cosThetaT);

    cosThetaT_ = (cosThetaI_ > 0) ? -cosThetaT : cosThetaT;

    /* No polarization -- return the unpolarized reflectance */
    return 0.5f * (Rs * Rs + Rp * Rp);
}

Float fresnelConductor(Float cosThetaI, std::complex<Float> eta_) {
    /* Modified from "Optics" by K.D. Moeller, University Science Books, 1988 */
    Float eta = eta_.real(), k = eta_.imag();

    Float cosThetaI2 = cosThetaI*cosThetaI,
          sinThetaI2 = 1-cosThetaI2,
          sinThetaI4 = sinThetaI2*sinThetaI2;

    Float temp1 = eta*eta - k*k - sinThetaI2,
          a2pb2 = math::safe_sqrt(temp1*temp1 + 4*k*k*eta*eta),
          a     = math::safe_sqrt(0.5f * (a2pb2 + temp1));

    Float term1 = a2pb2 + cosThetaI2,
          term2 = 2*a*cosThetaI;

    Float Rs2 = (term1 - term2) / (term1 + term2);

    Float term3 = a2pb2*cosThetaI2 + sinThetaI4,
          term4 = term2*sinThetaI2;

    Float Rp2 = Rs2 * (term3 - term4) / (term3 + term4);

    return 0.5f * (Rp2 + Rs2);
}

Float fresnelDielectricIntegral(Float eta) {
    /* Fast mode: the following code approximates the
     * diffuse Frensel reflectance for the eta<1 and
     * eta>1 cases. An evalution of the accuracy led
     * to the following scheme, which cherry-picks
     * fits from two papers where they are best.
     */
    if (eta < 1) {
        /* Fit by Egan and Hilgeman (1973). Works
           reasonably well for "normal" IOR values (<2).

           Max rel. error in 1.0 - 1.5 : 0.1%
           Max rel. error in 1.5 - 2   : 0.6%
           Max rel. error in 2.0 - 5   : 9.5%
        */
        return -1.4399f * (eta * eta)
              + 0.7099f * eta
              + 0.6681f
              + 0.0636f / eta;
    } else {
        /* Fit by d'Eon and Irving (2011)
         *
         * Maintains a good accuracy even for
         * unrealistic IOR values.
         *
         * Max rel. error in 1.0 - 2.0   : 0.1%
         * Max rel. error in 2.0 - 10.0  : 0.2%
         */
        Float invEta = 1.0f / eta,
              invEta2 = invEta*invEta,
              invEta3 = invEta2*invEta,
              invEta4 = invEta3*invEta,
              invEta5 = invEta4*invEta;

        return 0.919317f - 3.4793f * invEta
             + 6.75335f * invEta2
             - 7.80989f * invEta3
             + 4.98554f * invEta4
             - 1.36881f * invEta5;
    }
}

Float fresnelConductorIntegral(std::complex<Float> eta) {
    /* 10 point Gauss-Lobatto rule */

    Float nodes[] =   { 0.                    , 0.04023304591677057118,
                        0.13061306744724748841, 0.26103752509477773369,
                        0.41736052116680649737, 0.58263947883319344712,
                        0.73896247490522226631, 0.8693869325527525671 ,
                        0.95976695408322942882, 1.                     };
    Float weights[] = { 0.01111111111111111154, 0.0666529954255350443,
                        0.11244467103156326193, 0.14602134183984186167,
                        0.16376988059194869107, 0.16376988059194869107,
                        0.14602134183984186167, 0.11244467103156326193,
                        0.0666529954255350443,  0.01111111111111111154 };

    Float value = 0;
    for (int i=0; i<10; ++i)
        value += fresnelConductor(std::sqrt(nodes[i]), eta) * weights[i];

    return value;
}

NAMESPACE_END(layer)
