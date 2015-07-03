/*
    color.h -- Data types for representing RGB colors

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/common.h>
#include <Eigen/Core>

NAMESPACE_BEGIN(layer)

/// Data type for representing linearized RGB color values
template <typename Scalar> struct TColor3 : public Eigen::Array<Scalar, 3, 1> {
public:
    enum {
        Dimension = 3
    };

    typedef Eigen::Array<Scalar, 3, 1> Base;

    /// Initialize the color vector with a uniform value
    TColor3(Scalar value = (Scalar) 0) : Base(value, value, value) { }

    /// Initialize the color vector with specific per-channel values
    TColor3(Scalar r, Scalar g, Scalar b) : Base(r, g, b) { }

    /// Assign a color from a dense Eigen expression template
    template <typename Derived> TColor3(const Eigen::DenseBase<Derived>& p)
        : Base(p) { }

    /// Assign a color from a dense Eigen expression template
    template <typename Derived> TColor3 &operator=(const Eigen::DenseBase<Derived>& p) {
        this->Base::operator=(p);
        return *this;
    }

    /// Return a reference to the red channel
    Scalar &r() { return Base::x(); }
    /// Return a reference to the red channel (const version)
    const Scalar &r() const { return Base::x(); }

    /// Return a reference to the green channel
    Scalar &g() { return Base::y(); }
    /// Return a reference to the green channel (const version)
    const Scalar &g() const { return Base::y(); }

    /// Return a reference to the blue channel
    Scalar &b() { return Base::z(); }
    /// Return a reference to the blue channel (const version)
    const Scalar &b() const { return Base::z(); }

    /// Clamp to the positive range
    TColor3 clamp() const {
        return TColor3(std::max(r(), (Scalar) 0),
                       std::max(g(), (Scalar) 0),
                       std::max(b(), (Scalar) 0));
    }

    /// Check if the color vector contains a NaN/Inf/negative value
    bool isValid() const {
        for (int i = 0; i < 3; ++i) {
            Scalar value = Base::coeff(i);
            if (value < 0 || !std::isfinite(value))
                return false;
        }
        return true;
    }

    /// Convert from sRGB to linear RGB
    TColor3 toLinearRGB() const {
        TColor3 result;

        for (int i = 0; i < 3; ++i) {
            Scalar value = Base::coeff(i);

            if (value <= (Scalar) 0.04045)
                result[i] = value * ((Scalar) 1.0 / (Scalar) 12.92);
            else
                result[i] = std::pow((value + (Scalar) 0.055) *
                                     (1 / (Scalar) 1.055), (Scalar) 2.4);
        }

        return result;
    }

    /// Convert from linear RGB to sRGB
    TColor3 toSRGB() const {
        TColor3 result;
        for (int i = 0; i < 3; ++i) {
            Scalar value = Base::coeff(i);
            if (value <= (Scalar) 0.0031308)
                result[i] = (Scalar) 12.92 * value;
            else
                result[i] = (Scalar) 1.055 * std::pow(value,
                        (1 / (Scalar) 2.4)) - (Scalar) 0.055;
        }
        return result;
    }

    /// Return the associated luminance
    Scalar getLuminance() const {
        return Base::coeff(0) * (Scalar) 0.212671 +
               Base::coeff(1) * (Scalar) 0.715160 +
               Base::coeff(2) * (Scalar) 0.072169;
    }

    /// Stream operator
    friend std::ostream& operator<<(std::ostream& os, const TColor3 &c) {
        os << c.transpose(); return os;
    }
};

typedef TColor3<float>  Color3f;
typedef TColor3<double> Color3d;
typedef TColor3<Float>  Color3;

NAMESPACE_END(layer)
