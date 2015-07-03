/*
    common.h -- Basic macros and type definitions

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#if !defined(NAMESPACE_BEGIN)
#define NAMESPACE_BEGIN(name) namespace name {
#endif
#if !defined(NAMESPACE_END)
#define NAMESPACE_END(name) }
#endif

#define EIGEN_DONT_PARALLELIZE
#define EIGEN_DEFAULT_IO_FORMAT \
    Eigen::IOFormat(7, 0, ", ", ";\n", "", "", "[", "]")

#include <cstdint>
#include <cstddef>
#include <filesystem/fwd.h>

NAMESPACE_BEGIN(layer)

namespace fs = ::filesystem;

typedef double Float;

/* Forward declarations */
template <typename Scalar, int Dimension>  struct TVector;
template <typename Scalar, int Dimension>  struct TPoint;
template <typename Scalar>                 struct TNormal3;
template <typename Scalar>                 struct TFrame3;

/* Lots of aliases for various dimensions and data types */
typedef TVector<float, 1>       Vector1f;
typedef TVector<float, 2>       Vector2f;
typedef TVector<float, 3>       Vector3f;
typedef TVector<float, 4>       Vector4f;
typedef TVector<double, 1>      Vector1d;
typedef TVector<double, 2>      Vector2d;
typedef TVector<double, 3>      Vector3d;
typedef TVector<double, 4>      Vector4d;
typedef TVector<Float, 1>       Vector1;
typedef TVector<Float, 2>       Vector2;
typedef TVector<Float, 3>       Vector3;
typedef TVector<Float, 4>       Vector4;
typedef TVector<int32_t, 1>     Vector1i;
typedef TVector<int32_t, 2>     Vector2i;
typedef TVector<int32_t, 3>     Vector3i;
typedef TVector<int32_t, 4>     Vector4i;
typedef TVector<uint32_t, 1>    Vector1u;
typedef TVector<uint32_t, 2>    Vector2u;
typedef TVector<uint32_t, 3>    Vector3u;
typedef TVector<uint32_t, 4>    Vector4u;
typedef TVector<size_t, 1>      Vector1s;
typedef TVector<size_t, 2>      Vector2s;
typedef TVector<size_t, 3>      Vector3s;
typedef TVector<size_t, 4>      Vector4s;
typedef TPoint<float, 1>        Point1f;
typedef TPoint<float, 2>        Point2f;
typedef TPoint<float, 3>        Point3f;
typedef TPoint<float, 4>        Point4f;
typedef TPoint<double, 1>       Point1d;
typedef TPoint<double, 2>       Point2d;
typedef TPoint<double, 3>       Point3d;
typedef TPoint<double, 4>       Point4d;
typedef TPoint<Float, 1>        Point1;
typedef TPoint<Float, 2>        Point2;
typedef TPoint<Float, 3>        Point3;
typedef TPoint<Float, 4>        Point4;
typedef TPoint<int, 1>          Point1i;
typedef TPoint<int, 2>          Point2i;
typedef TPoint<int, 3>          Point3i;
typedef TPoint<int, 4>          Point4i;
typedef TNormal3<Float>         Normal3;
typedef TNormal3<float>         Normal3f;
typedef TNormal3<double>        Normal3d;
typedef TFrame3<float>          Frame3f;
typedef TFrame3<double>         Frame3d;
typedef TFrame3<Float>          Frame3;
typedef Normal3                 Normal;
typedef Vector3                 Vector;
typedef Point3                  Point;
typedef Frame3                  Frame;

NAMESPACE_END(layer)
