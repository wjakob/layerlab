/*
    vector.h -- This file contains templates and specializations for
    2/3D points, vectors, and normals over different underlying data types.
    These all transform differently under homogeneous coordinate
    transformations and are thus implemented as separate types.

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/common.h>
#include <Eigen/Core>
#include <Eigen/src/Geometry/OrthoMethods.h> // for cross products

NAMESPACE_BEGIN(layer)

/* Dynamic vectors */
typedef Eigen::Matrix<Float,  Eigen::Dynamic, 1> VectorX;
typedef Eigen::Matrix<float,  Eigen::Dynamic, 1> VectorXf;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

/* Dynamic matrices */
typedef Eigen::Matrix<Float,  Eigen::Dynamic, Eigen::Dynamic> MatrixX;
typedef Eigen::Matrix<float,  Eigen::Dynamic, Eigen::Dynamic> MatrixXf;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;

/// Generic N-dimensional vector data structure based on Eigen::Matrix
template <typename Scalar, int _Dimension> struct TVector : public Eigen::Matrix<Scalar, _Dimension, 1> {
public:
    enum {
        Dimension = _Dimension
    };

    typedef TVector<Scalar, Dimension>          VectorType;
    typedef TPoint<Scalar, Dimension>           PointType;
    typedef Eigen::Matrix<Scalar, Dimension, 1> Base;

    /// Create a new vector with constant component values
    TVector(Scalar value = (Scalar) 0) { Base::setConstant(value); }

    /// Create a new 2D vector (type error if \c Dimension != 2)
    TVector(Scalar x, Scalar y) : Base(x, y) { }

    /// Create a new 3D vector (type error if \c Dimension != 3)
    TVector(Scalar x, Scalar y, Scalar z) : Base(x, y, z) { }

    /// Create a new 4D vector (type error if \c Dimension != 4)
    TVector(Scalar x, Scalar y, Scalar z, Scalar w) : Base(x, y, z, w) { }

    /// Construct a vector from a dense Eigen expression template
    template <typename Derived> TVector(const Eigen::DenseBase<Derived>& p)
        : Base(p) { }

    /// Assign a vector from a dense Eigen expression template
    template <typename Derived> TVector &operator=(const Eigen::DenseBase<Derived>& p) {
        this->Base::operator=(p); return *this;
    }

    /// Stream operator
    friend std::ostream& operator<<(std::ostream& os, const TVector& v) {
        os << v.transpose(); return os;
    }
};

/// Generic N-dimensional point data structure based on Eigen::Matrix
template <typename Scalar, int _Dimension> struct TPoint : public Eigen::Matrix<Scalar, _Dimension, 1> {
public:
    enum {
        Dimension = _Dimension
    };

    typedef TVector<Scalar, Dimension>          VectorType;
    typedef TPoint<Scalar, Dimension>           PointType;
    typedef Eigen::Matrix<Scalar, Dimension, 1> Base;

    /// Create a new point with constant component values
    TPoint(Scalar value = (Scalar) 0) { Base::setConstant(value); }

    /// Create a new 2D point (type error if \c Dimension != 2)
    TPoint(Scalar x, Scalar y) : Base(x, y) { }

    /// Create a new 3D point (type error if \c Dimension != 3)
    TPoint(Scalar x, Scalar y, Scalar z) : Base(x, y, z) { }

    /// Create a new 4D point (type error if \c Dimension != 4)
    TPoint(Scalar x, Scalar y, Scalar z, Scalar w) : Base(x, y, z, w) { }

    /// Construct a point from a dense Eigen expression template
    template <typename Derived> TPoint(const Eigen::DenseBase<Derived>& p)
        : Base(p) { }

    /// Assign a point from a dense Eigen expression template
    template <typename Derived> TPoint &operator=(const Eigen::DenseBase<Derived>& p) {
        this->Base::operator=(p); return *this;
    }

    /// Stream operator
    friend std::ostream& operator<<(std::ostream& os, const TPoint& v) {
        os << v.transpose(); return os;
    }
};

/// 3-dimensional surface normal representation
template <typename Scalar> struct TNormal3 : public TVector<Scalar, 3> {
public:
    enum {
        Dimension = 3
    };

    typedef TVector<Scalar, Dimension>          VectorType;
    typedef TPoint<Scalar, Dimension>           PointType;
    typedef VectorType                          Base;

    /// Create a new normal with constant component values
    TNormal3(Scalar value = (Scalar) 0) { Base::setConstant(value); }

    /// Create a new 3D normal
    TNormal3(Scalar x, Scalar y, Scalar z) : Base(x, y, z) { }

    /// Construct a normal from a dense Eigen expression template
    template <typename Derived> TNormal3(const Eigen::DenseBase<Derived>& p)
        : Base(p) { }

    /// Assign a normal from a dense Eigen expression template
    template <typename Derived> TNormal3 &operator=(const Eigen::DenseBase<Derived>& p) {
        this->Base::operator=(p);
        return *this;
    }

    /// Stream operator
    friend std::ostream& operator<<(std::ostream& os, const TNormal3& v) {
        os << v.transpose(); return os;
    }
};

/**
 * \brief Given the unit vector n, find an orthonormal basis {s, t, n}
 *
 * Based on
 * "Building an Orthonormal Basis from a 3D Unit Vector Without Normalization"
 * by Jeppe Revall Frisvad
 * in "Journal of Graphics Tools" 16(3), pp. 151-159, August 2012.
 */
template <typename Scalar>
void coordinateSystem(const TVector<Scalar, 3> &n, TVector<Scalar, 3> &s,
                      TVector<Scalar, 3> &t) {
    if (n.z() < (Scalar) -0.9999999) { // Handle the singularity
        s = TVector<Scalar, 3>(0, -1, 0);
        t = TVector<Scalar, 3>(-1, 0, 0);
        return;
    }
    const Scalar a = (Scalar) 1 / ((Scalar) 1 + n.z());
    const Scalar b = -n.x() * n.y() * a;
    s = TVector<Scalar, 3>((Scalar) 1 - n.x() * n.x() * a, b, -n.x());
    t = TVector<Scalar, 3>(b, (Scalar) 1 - n.y() * n.y() * a, -n.y());
}

NAMESPACE_END(layer)
