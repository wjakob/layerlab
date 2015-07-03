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

#include <layer/vector.h>
#include <layer/math.h>

NAMESPACE_BEGIN(layer)

/**
 * \brief Stores a three-dimensional orthonormal coordinate frame
 *
 * This class converts between different cartesian coordinate systems and
 * efficiently computes certain derived quantities based on spherical
 * coordinates (e.g. \ref cosTheta(), \ref tanTheta(), ..).
 */
template <typename Scalar> struct TFrame3 {
    typedef TVector<Scalar, 3> VectorType;
    typedef TNormal3<Scalar>   NormalType;

    /// First tangent
    VectorType s;
    /// Second tangent
    VectorType t;
    /// Normal direction
    NormalType n;

    /// Default constructor -- performs no initialization!
    TFrame3() { }

    /// Construct a new coordinate frame from a single vector
    TFrame3(const VectorType &n) : n(n) {
        coordinateSystem(n, s, t);
    }

    /// Construct a frame from the given orthonormal vectors
    TFrame3(const VectorType &s, const VectorType &t, const VectorType &n)
        : s(s), t(t), n(n) { }

    /// Convert from world coordinates to local coordinates
    VectorType toLocal(const VectorType &v) const {
        return VectorType(
            v.dot(s),
            v.dot(t),
            v.dot(n)
        );
    }

    /// Convert from local coordinates to world coordinates
    VectorType toWorld(const VectorType &v) const {
        return s * v.x() + t * v.y() + n * v.z();
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared sine of the angle between the normal and v */
    static Scalar sinTheta2(const VectorType &v) {
        return (Scalar) 1 - v.z() * v.z();
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared cosine of the angle between the normal and v */
    static Scalar cosTheta2(const VectorType &v) {
        return v.z() * v.z();
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared tangent of the angle between the normal and v */
    static Scalar tanTheta2(const VectorType &v) {
        return sinTheta2(v) / cosTheta2(v);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the sine of the angle between the normal and v */
    static Scalar sinTheta(const VectorType &v) {
        return math::safe_sqrt(sinTheta2(v));
    }

    /** \brief Assuming that the given direction is in the local coordinate 
     * system, return the cosine of the angle between the normal and v */
    static Scalar cosTheta(const VectorType &v) {
        return v.z();
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the tangent of the angle between the normal and v */
    static Scalar tanTheta(const VectorType &v) {
        return sinTheta(v) / cosTheta(v);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared sine of the phi parameter in spherical
     * coordinates */
    static Scalar sinPhi2(const VectorType &v) {
        Scalar sinTheta2 = TFrame3::sinTheta2(v);
        if (sinTheta2 == 0)
            return 0;
        return math::clamp(v.y() * v.y() / sinTheta2, (Scalar) 0, (Scalar) 1);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared cosine of the phi parameter in spherical
     * coordinates */
    static Scalar cosPhi2(const VectorType &v) {
        Scalar sinTheta2 = TFrame3::sinTheta2(v);
        if (sinTheta2 == 0)
            return 0;
        return math::clamp(v.x() * v.x() / sinTheta2, (Scalar) 0, (Scalar) 1);
    }

    /** \brief Assuming that the given direction is in the local coordinate 
     * system, return the sine of the phi parameter in spherical coordinates */
    static Scalar sinPhi(const VectorType &v) {
        Scalar sinTheta = TFrame3::sinTheta(v);
        if (sinTheta == 0)
            return 0;
        return math::clamp(v.y() / sinTheta, (Scalar) -1, (Scalar) 1);
    }

    /** \brief Assuming that the given direction is in the local coordinate 
     * system, return the cosine of the phi parameter in spherical coordinates */
    static Scalar cosPhi(const VectorType &v) {
        Scalar sinTheta = TFrame3::sinTheta(v);
        if (sinTheta == 0)
            return 0;
        return math::clamp(v.x() / sinTheta, (Scalar) -1, (Scalar) 1);
    }

    /// Equality test
    bool operator==(const TFrame3 &frame) const {
        return frame.s == s && frame.t == t && frame.n == n;
    }

    /// Inequality test
    bool operator!=(const TFrame3 &frame) const {
        return !operator==(frame);
    }
 
    /// Stream operator
    friend std::ostream& operator<<(std::ostream& os, const TFrame3& f) {
        os << "Frame[" << std::endl
           << "  s = " << f.s << "," << std::endl
           << "  t = " << f.t << "," << std::endl
           << "  n = " << f.n << std::endl
           << "]";
        return os;
    }
};

NAMESPACE_END(layer)
