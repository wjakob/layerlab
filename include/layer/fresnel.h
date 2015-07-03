/*
    fresnel.h -- Fresnel coefficients for dielectrics and conductors

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/common.h>
#include <complex>

NAMESPACE_BEGIN(layer)

/**
 * \brief Calculates the unpolarized Fresnel reflection coefficient
 * at a planar interface between two dielectrics
 *
 * This function also computes the transmitted direction and returns it using
 * the \c cosThetaT argument. When encountering total internal reflection, it
 * sets <tt>cosThetaT=0</tt> and returns the value 1.
 *
 * When <tt>cosThetaI < 0</tt>, the function computes the Fresnel reflectance
 * from the \a internal boundary, which is equivalent to calling the function
 * with arguments <tt>fresnelDielectric(abs(cosThetaI), cosThetaT, 1/eta)</tt>.
 *
 * \remark When accessed from Python, this function has the signature
 * "<tt>F, cosThetaT = fresnelDielectric(cosThetaI, eta)</tt>".
 *
 * \param cosThetaI
 *      Cosine of the angle between the normal and the incident ray
 *      (may be negative)
 * \param cosThetaT
 *      Argument used to return the cosine of the angle between the normal
 *      and the transmitted ray, will have the opposite sign of \c cosThetaI
 * \param eta
 *      Relative refractive index
 */
extern Float fresnelDielectric(Float cosThetaI, Float &cosThetaT, Float eta);

/**
 * \brief Calculates the unpolarized Fresnel reflection coefficient
 * at a planar interface between two dielectrics
 *
 * This is just a convenience wrapper function around the other \c
 * fresnelDielectric function, which does not return the transmitted direction
 * cosine in case it is not needed by the application.
 *
 * \param cosThetaI
 *      Cosine of the angle between the normal and the incident ray
 * \param eta
 *      Relative refractive index
 */
inline Float fresnelDielectric(Float cosThetaI, Float eta) {
    Float unused;
    return fresnelDielectric(cosThetaI, unused, eta);
}

/**
 * \brief Calculates the unpolarized Fresnel reflection coefficient
 * at a planar interface having a complex-valued relative index of
 * refraction
 *
 * \remark The name of this function is a slight misnomer, since it supports
 * the general case of a complex-valued relative index of refraction (rather
 * than being restricted to conductors)
 *
 * \param cosThetaI
 *      Cosine of the angle between the normal and the incident ray
 * \param eta
 *    Relative index of refraction (complex)
 */
extern Float fresnelConductor(Float cosThetaI, std::complex<Float> eta);

/**
 * \brief Calculates the diffuse unpolarized Fresnel reflectance of
 * a dielectric material (sometimes referred to as "Fdr").
 *
 * This value quantifies what fraction of diffuse incident illumination
 * will, on average, be reflected at a dielectric material boundary
 *
 * \param eta
 *      Relative refraction coefficient
 */
extern Float fresnelDielectricIntegral(Float eta);

/**
 * \brief Calculates the diffuse unpolarized Fresnel reflectance of
 * a conductor
 *
 * This value quantifies what fraction of diffuse incident illumination
 * will, on average, be reflected at a conductive material boundary
 *
 * \param eta
 *      Relative refractive index (real component)
 * \param k
 *      Relative refractive index (imaginary component)
 */
extern Float fresnelConductorIntegral(std::complex<Float> eta);

NAMESPACE_END(layer)
