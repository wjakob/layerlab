/*
    spline.h -- Spline evaluation and sampling routines

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/math.h>
#include <layer/vector.h>
#include <iostream>

NAMESPACE_BEGIN(layer)
NAMESPACE_BEGIN(spline)

#define GET_SPLINE_UNIFORM(idx) \
    Scalar f0 = values[idx], \
           f1 = values[idx+1], \
           d0, d1; \
    if (idx > 0) \
        d0 = (Scalar) 0.5 * (f1 - values[idx-1]); \
    else \
        d0 = f1 - f0; \
    \
    if (idx + 2 < size) \
        d1 = (Scalar) 0.5 * (values[idx+2] - f0); \
    else \
        d1 = f1 - f0;

#define GET_SPLINE_NONUNIFORM(idx) \
    Scalar f0    = values[idx], \
           f1    = values[idx+1], \
           x0    = nodes[idx], \
           x1    = nodes[idx+1], \
           width = x1 - x0, \
           d0, d1; \
    if (idx > 0) \
        d0 = width * (f1 - values[idx-1]) / (x1 - nodes[idx-1]); \
    else \
        d0 = f1 - f0; \
    \
    if (idx + 2 < size) \
        d1 = width * (values[idx+2] - f0) / (nodes[idx+2] - x0); \
    else \
        d1 = f1 - f0;

// -----------------------------------------------------------------------
//! @{ \name Functions for evaluating and sampling cubic Catmull-Rom splines
// -----------------------------------------------------------------------


/**
 * \brief Compute the definite integral and derivative of a cubic spline that
 * is parameterized by the function values and derivatives at the endpoints
 * of the interval <tt>[0, 1]</tt>.
 *
 * \param f0
 *      The function value at the left position
 * \param f1
 *      The function value at the right position
 * \param d0
 *      The function derivative at the left position
 * \param d1
 *      The function derivative at the right position
 * \param t
 *      The parameter variable
 * \return
 *      The interpolated function value at \c t
 */
template <typename Scalar>
Scalar evalSpline(Scalar f0, Scalar f1, Scalar d0, Scalar d1, Scalar t) {
    Scalar t2 = t*t, t3 = t2*t;
    return ( 2*t3 - 3*t2 + 1) * f0 +
           (-2*t3 + 3*t2)     * f1 +
           (   t3 - 2*t2 + t) * d0 +
           (   t3 - t2)       * d1;
}

/**
 * \brief Compute the value and derivative of a cubic spline that is
 * parameterized by the function values and derivatives of the
 * interval <tt>[0, 1]</tt>.
 *  
 * \param f0
 *      The function value at the left position
 * \param f1
 *      The function value at the right position
 * \param d0
 *      The function derivative at the left position
 * \param d1
 *      The function derivative at the right position
 * \param t
 *      The parameter variable
 * \return
 *      The interpolated function value and
 *      its derivative at \c t
 */
template <typename Scalar>
std::pair<Scalar, Scalar> evalSplineD(Scalar f0, Scalar f1, Scalar d0,
                                      Scalar d1, Scalar t) {
    Scalar t2 = t*t, t3 = t2*t;
    return std::make_pair(
        /* Function value */
        ( 2*t3 - 3*t2 + 1) * f0 +
        (-2*t3 + 3*t2)     * f1 +
        (   t3 - 2*t2 + t) * d0 +
        (   t3 - t2)       * d1,

        /* Derivative */
        ( 6*t2 - 6*t)      * f0 +
        (-6*t2 + 6*t)      * f1 +
        ( 3*t2 - 4*t + 1)  * d0 +
        ( 3*t2 - 2*t)      * d1
    );
}

/**
 * \brief Compute the definite integral and value of a cubic spline
 * that is parameterized by the function values and derivatives of 
 * the interval <tt>[0, 1]</tt>.
 *  
 * \param f0
 *      The function value at the left position
 * \param f1
 *      The function value at the right position
 * \param d0
 *      The function derivative at the left position
 * \param d1
 *      The function derivative at the right position
 * \return
 *      The definite integral and the interpolated
 *      function value at \c t
 */
template <typename Scalar>
std::pair<Scalar, Scalar> evalSplineI(Scalar f0, Scalar f1, Scalar d0,
                                      Scalar d1, Scalar t) {
    Scalar t2 = t*t, t3 = t2*t, t4 = t2*t2;
    const Scalar H = (Scalar) 0.5;
    const Scalar T = (Scalar) (1.0 / 3.0);
    const Scalar Q = (Scalar) 0.25;

    return std::make_pair(
        /* Definite integral */
        ( H*t4 - t3 + t)         * f0 +
        (-H*t4 + t3)             * f1 +
        ( Q*t4 - 2*T*t3 + H*t2)  * d0 +
        ( Q*t4 - T*t3)           * d1,

        /* Function value */
        ( 2*t3 - 3*t2 + 1)       * f0 +
        (-2*t3 + 3*t2)           * f1 +
        (   t3 - 2*t2 + t)       * d0 +
        (   t3 - t2)             * d1
    );
}

/**
 * \brief Evaluate a cubic spline interpolant of a \a uniformly sampled 1D function
 *
 * The implementation relies on Catmull-Rom splines, i.e. it uses finite
 * differences to approximate the derivatives at the endpoints of each spline
 * segment.
 *
 * \param min
 *      Position of the first node
 * \param max
 *      Position of the last node
 * \param values
 *      Array containing \c size regularly spaced evaluations in the range [\c
 *      min, \c max] of the approximated function.
 * \param size
 *      Denotes the size of the \c values array
 * \param x
 *      Evaluation point
 * \param extrapolate
 *      Extrapolate values when \c x is out of range? (default: \c false)
 * \remark
 *      The Python API lacks the \c size parameter, which is inferred
 *      automatically from the size of the input array.
 * \remark
 *      The Python API provides a vectorized version which evaluates
 *      the function for many arguments \c x.
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>
 *      and \c x lies outside of [\c min, \c max]
 */
template <typename Scalar>
Scalar eval1D(Scalar min, Scalar max, const Scalar *values,
              size_t size, Scalar x, bool extrapolate = false) {
    /* Give up when given an out-of-range or NaN argument */
    if (!(x >= min && x <= max) && !extrapolate)
        return (Scalar) 0;

    /* Transform 'x' so that nodes lie at integer positions */
    Scalar t = ((x - min) * (size - 1)) / (max - min);

    /* Find the index of the left node in the queried subinterval */
    size_t idx = std::max((size_t) 0, std::min((size_t) t, size - 2));

    GET_SPLINE_UNIFORM(idx);

    /* Compute the relative position within the interval */
    t -= (Scalar) idx;

    return evalSpline(f0, f1, d0, d1, t);
}

/**
 * \brief Evaluate a cubic spline interpolant of a \a non-uniformly sampled 1D function
 *
 * The implementation relies on Catmull-Rom splines, i.e. it uses finite
 * differences to approximate the derivatives at the endpoints of each spline
 * segment.
 *
 * \param nodes
 *      Array containing \c size non-uniformly spaced values denoting positions
 *      the where the function to be interpolated was evaluated. They must be
 *      provided in \a increasing order.
 * \param values
 *      Array containing function evaluations matched to the entries of \c
 *      nodes.
 * \param size
 *      Denotes the size of the \c nodes and \c values array
 * \param x
 *      Evaluation point
 * \param extrapolate
 *      Extrapolate values when \c x is out of range? (default: \c false)
 * \remark
 *      The Python API lacks the \c size parameter, which is inferred
 *      automatically from the size of the input array
 * \remark
 *      The Python API provides a vectorized version which evaluates
 *      the function for many arguments \c x.
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>
 *      and \c x lies outside of \a [\c min, \c max]
 */
template <typename Scalar>
Scalar eval1D(const Scalar *nodes, const Scalar *values, size_t size,
              Scalar x, bool extrapolate = false) {
    /* Give up when given an out-of-range or NaN argument */
    if (!(x >= nodes[0] && x <= nodes[size-1]) && !extrapolate)
        return (Scalar) 0;

    /* Find the index of the left node in the queried subinterval */
    size_t idx = math::findInterval(size, 
        [&](size_t i) { return nodes[i] <= x; }
    );

    GET_SPLINE_NONUNIFORM(idx);

    /* Compute the relative position within the interval */
    Scalar t = (x - nodes[idx]) / width;

    return evalSpline(f0, f1, d0, d1, t);
}

/**
 * \brief Computes a prefix sum of integrals over segments of a \a uniformly
 * sampled 1D Catmull-Rom spline interpolant
 *
 * This is useful for sampling spline segments as part of an importance
 * sampling scheme (in conjunction with \ref sample1D)
 *
 * \param min
 *      Position of the first node
 * \param max
 *      Position of the last node
  * \param values
 *      Array containing \c size regularly spaced evaluations in the range
 *      [\c min, \c max] of the approximated function.
 * \param size
 *      Denotes the size of the \c values array
 * \param[out] out
 *      An array with \c size entries, which will be used to store the
 *      prefix sum
 * \remark
 *      The Python API lacks the \c size and \c out parameters. The former 
 *      is inferred automatically from the size of the input array, and \c out
 *      is returned as a list.
 */
template <typename Scalar>
void integrate1D(Scalar min, Scalar max, const Scalar *values,
                 size_t size, Scalar *out) {
    const Scalar width = (max - min) / (size - 1);
    Scalar sum = 0;
    out[0] = 0;
    for (size_t idx = 0; idx < size - 1; ++idx) {
        GET_SPLINE_UNIFORM(idx);

        sum += ((d0 - d1) * (Scalar) (1.0 / 12.0) +
                (f0 + f1) * (Scalar) 0.5) * width;

        out[idx + 1] = sum;
    }
}

/**
 * \brief Computes a prefix sum of integrals over segments of a \a non-uniformly
 * sampled 1D Catmull-Rom spline interpolant
 *
 * This is useful for sampling spline segments as part of an importance
 * sampling scheme (in conjunction with \ref sample1D)
 *
 * \param nodes
 *      Array containing \c size non-uniformly spaced values denoting positions
 *      the where the function to be interpolated was evaluated. They must be
 *      provided in \a increasing order.
 * \param values
 *      Array containing function evaluations matched to the entries of
 *      \c nodes.
 * \param size
 *      Denotes the size of the \c values array
 * \param[out] out
 *      An array with \c size entries, which will be used to store the
 *      prefix sum
 * \remark
 *      The Python API lacks the \c size and \c out parameters. The former 
 *      is inferred automatically from the size of the input array, and \c out
 *      is returned as a list.
 */
template <typename Scalar>
void integrate1D(const Scalar *nodes, const Scalar *values,
                 size_t size, Scalar *out) {
    Scalar sum = 0;
    out[0] = 0;
    for (size_t idx = 0; idx < size - 1; ++idx) {
        GET_SPLINE_NONUNIFORM(idx);

        sum += ((d0 - d1) * (Scalar) (1.0 / 12.0) +
                (f0 + f1) * (Scalar) 0.5) * width;

        out[idx + 1] = sum;
    }
}

/**
 * \brief Invert a cubic spline interpolant of a \a uniformly sampled 1D function.
 * The spline interpolant must be <em>monotonically increasing</em>.
 *
 * \param min
 *      Position of the first node
 * \param max
 *      Position of the last node
 * \param values
 *      Array containing \c size regularly spaced evaluations in the range
 *      [\c min, \c max] of the approximated function.
 * \param size
 *      Denotes the size of the \c values array
 * \param y
 *      Input parameter for the inversion
 * \return
 *      The spline parameter \c t such that <tt>eval1D(..., t)=y</tt>
 */
template <typename Scalar>
Scalar invert1D(Scalar min, Scalar max, const Scalar *values, size_t size, Scalar y) {
    /* Give up when given an out-of-range or NaN argument */
    if (!(y > values[0]))
        return min;
    else if (!(y < values[size-1]))
        return max;

    /* Map y to a spline interval by searching through the
       'values' array (which is assumed to be monotonic) */
    size_t idx = math::findInterval(size,
        [&](size_t i) { return values[i] <= y; }
    );

    const Scalar width = (max - min) / (size - 1);
    GET_SPLINE_UNIFORM(idx);

    /* Invert the spline interpolant using Newton-Bisection */
    Scalar a = 0, b = 1, t = (Scalar) 0.5;
    Scalar value, deriv;

    while (true) {
        /* Fall back to a bisection step when t is out of bounds */
        if (!(t >= a && t <= b))
            t = (Scalar) 0.5 * (a + b);

        /* Evaluate the spline and its derivative */
        std::tie(value, deriv)
            = evalSplineD(f0, f1, d0, d1, t);
        value -= y;

        /* Stop the iteration if converged */
        if (std::abs(value) < (Scalar) 1e-6 || b-a < (Scalar) 1e-6)
            break;

        /* Update the bisection bounds */
        if (value > 0)
            b = t;
        else
            a = t;

        /* Perform a Newton step */
        t -= value / deriv;
    }

    return min + (idx + t) * width;
}

/**
 * \brief Invert a cubic spline interpolant of a \a non-uniformly sampled 1D function.
 * The spline interpolant must be <em>monotonically increasing</em>.
 *
 * \param nodes
 *      Array containing \c size non-uniformly spaced values denoting positions
 *      the where the function to be interpolated was evaluated. They must be
 *      provided in \a increasing order.
 * \param values
 *      Array containing function evaluations matched to the entries of
 *      \c nodes.
 * \param size
 *      Denotes the size of the \c values array
 * \param y
 *      Input parameter for the inversion
 * \return
 *      The spline parameter \c t such that <tt>eval1D(..., t)=y</tt>
 */
template <typename Scalar>
Scalar invert1D(const Scalar *nodes, const Scalar *values, size_t size, Scalar y) {
    /* Give up when given an out-of-range or NaN argument */
    if (!(y > values[0]))
        return nodes[0];
    else if (!(y < values[size-1]))
        return nodes[size-1];

    /* Map y to a spline interval by searching through the
       'values' array (which is assumed to be monotonic) */
    size_t idx = math::findInterval(size,
        [&](size_t i) { return values[i] <= y; }
    );

    GET_SPLINE_NONUNIFORM(idx);

    /* Invert the spline interpolant using Newton-Bisection */
    Scalar a = 0, b = 1, t = (Scalar) 0.5;

    while (true) {
        /* Fall back to a bisection step when t is out of bounds */
        if (!(t >= a && t <= b))
            t = (Scalar) 0.5 * (a + b);

        /* Evaluate the spline and its derivative */
        Scalar value, deriv;
        std::tie(value, deriv)
            = evalSplineD(f0, f1, d0, d1, t);
        value -= y;

        /* Stop the iteration if converged */
        if (std::abs(value) < (Scalar) 1e-6 || b-a < (Scalar) 1e-6)
            break;

        /* Update the bisection bounds */
        if (value > 0)
            b = t;
        else
            a = t;

        /* Perform a Newton step */
        t -= value / deriv;
    }

    return x0 + t*width;
}

/**
 * \brief Importance sample a segment of a \a uniformly sampled 1D Catmull-Rom
 * spline interpolant
 *
 * \param min
 *      Position of the first node
 * \param max
 *      Position of the last node
 * \param values
 *      Array containing \c size regularly spaced evaluations in the range [\c
 *      min, \c max] of the approximated function.
 * \param cdf
 *      Array containing a cumulative distribution function computed by \ref
 *      integrate1D().
 * \param size
 *      Denotes the size of the \c values array
 * \param sample
 *      A uniformly distributed random sample in the interval <tt>[0,1]</tt>
 * \param[out] fval
 *      If set to a non-null pointer, this argument will be used to return
 *      the value of the spline at the sampled position
 * \param[out] pdf
 *      If set to a non-null pointer, this argument will be used to return
 *      the probability density at the sampled position (which only differs
 *      from \c fval when the function does not integrate to 1)
 * \remark
 *      The Python API lacks the \c size, \c fval and \c pdf parameters. The
 *      first is automatically inferred from the size of the input array, and
 *      \c fval and \c pdf are returned as the second and third element of the
 *      return value, which is now a tuple.
 * \return
 *      The sampled position
 */
template <typename Scalar>
Scalar sample1D(Scalar min, Scalar max, const Scalar *values, const Scalar *cdf,
                size_t size, Scalar sample, Scalar *fval, Scalar *pdf) {
    /* Scale by the definite integral of the function (in case
       it is not normalized) */
    Scalar last = cdf[size-1];
    sample *= last;

    /* Map y to a spline interval by searching through the
       monotonic 'cdf' array */
    size_t idx = math::findInterval(size,
        [&](size_t i) { return cdf[i] <= sample; }
    );

    const Scalar width = (max - min) / (size - 1);
    GET_SPLINE_UNIFORM(idx);

    // Re-scale the sample after having choosen the interval
    sample = (sample - cdf[idx]) / width;

    /* Importance sample linear interpolant as initial guess for 't'*/
    Scalar t;
    if (f0 != f1)
        t = (f0 - math::safe_sqrt(f0 * f0 + 2 * sample * (f1 - f0))) / (f0 - f1);
    else
        t = sample / f0;

    Scalar a = 0, b = 1, value, deriv;
    while (true) {
        /* Fall back to a bisection step when t is out of bounds */
        if (!(t >= a && t <= b))
            t = 0.5f * (a + b);

        /* Evaluate the definite integral and its derivative
           (i.e. the spline) */
        std::tie(value, deriv)
            = evalSplineI(f0, f1, d0, d1, t);
        value -= sample;

        /* Stop the iteration if converged */
        if (std::abs(value) < (Scalar) 1e-6 || b-a < (Scalar) 1e-6)
            break;

        /* Update the bisection bounds */
        if (value > 0)
            b = t;
        else
            a = t;

        // Perform a Newton step
        t -= value / deriv;
    }

    /* Return the value and PDF if requested */
    if (fval)
        *fval = deriv;
    if (pdf)
        *pdf = deriv / last;

    return min + (idx + t) * width;
}

/**
 * \brief Importance sample a segment of a \a non-uniformly sampled 1D Catmull-Rom
 * spline interpolant
 *
 * \param nodes
 *      Array containing \c size non-uniformly spaced values denoting positions
 *      the where the function to be interpolated was evaluated. They must be
 *      provided in \a increasing order.
 * \param values
 *      Array containing function evaluations matched to
 *      the entries of \c nodes.
 * \param cdf
 *      Array containing a cumulative distribution function computed by \ref
 *      integrate1D().
 * \param size
 *      Denotes the size of the \c values array
 * \param sample
 *      A uniformly distributed random sample in the interval <tt>[0,1]</tt>
 * \param[out] fval
 *      If set to a non-null pointer, this argument will be used to return
 *      the value of the spline at the sampled position
 * \param[out] pdf
 *      If set to a non-null pointer, this argument will be used to return
 *      the probability density at the sampled position (which only differs
 *      from \c fval when the function does not integrate to 1)
 * \remark
 *      The Python API lacks the \c size, \c fval and \c pdf parameters. The
 *      first is automatically inferred from the size of the input array, and
 *      \c fval and \c pdf are returned as the second and third element of the
 *      return value, which is now a tuple.
 * \return
 *      The sampled position
 */
template <typename Scalar>
Scalar sample1D(const Scalar *nodes, const Scalar *values, const Scalar *cdf,
                size_t size, Scalar sample, Scalar *fval, Scalar *pdf) {
    /* Scale by the definite integral of the function (in case
       it is not normalized) */
    Scalar last = cdf[size-1];
    sample *= last;

    /* Map y to a spline interval by searching through the
       monotonic 'cdf' array */
    size_t idx = math::findInterval(size,
        [&](size_t i) { return cdf[i] <= sample; }
    );

    GET_SPLINE_NONUNIFORM(idx);

    // Re-scale the sample after having choosen the interval
    sample = (sample - cdf[idx]) / width;

    /* Importance sample linear interpolant as initial guess for 't'*/
    Scalar t;
    if (f0 != f1)
        t = (f0 - math::safe_sqrt(f0 * f0 + 2 * sample * (f1 - f0))) / (f0 - f1);
    else
        t = sample / f0;

    Scalar a = 0, b = 1, value, deriv;
    while (true) {
        /* Fall back to a bisection step when t is out of bounds */
        if (!(t >= a && t <= b))
            t = 0.5f * (a + b);

        /* Evaluate the definite integral and its derivative
           (i.e. the spline) */
        std::tie(value, deriv)
            = evalSplineI(f0, f1, d0, d1, t);
        value -= sample;

        /* Stop the iteration if converged */
        if (std::abs(value) < (Scalar) 1e-6 || b-a < (Scalar) 1e-6)
            break;

        /* Update the bisection bounds */
        if (value > 0)
            b = t;
        else
            a = t;

        // Perform a Newton step
        t -= value / deriv;
    }

    /* Return the value and PDF if requested */
    if (fval)
        *fval = deriv;
    if (pdf)
        *pdf = deriv / last;

    return x0 + width * t;
}

/**
 * \brief Compute weights to perform a spline-interpolated lookup on a
 * \a uniformly sampled 1D function.
 *
 * The implementation relies on Catmull-Rom splines, i.e. it uses finite
 * differences to approximate the derivatives at the endpoints of each spline
 * segment. The resulting weights are identical those internally used by \ref
 * sample1D().
 *
 * \param min
 *      Position of the first node
 * \param max
 *      Position of the last node
 * \param size
 *      Denotes the number of function samples
 * \param x
 *      Evaluation point
 * \param[out] weights
 *      Pointer to a weight array of size 4 that will be populated
 * \param[out] offset
 *      Offset into the function samples associated with weights[0]
 * \param extrapolate
 *      Extrapolate values when \c x is out of range? (default: \c false)
 * \remark
 *      In the Python API, the \c offset and \c weights parameters are returned
 *      as the second and third elements of a triple.
 * \return
 *      \c true on success and \c false when <tt>extrapolate=false</tt>
 *      and \c x lies outside of [\c min, \c max]
 */
template <typename Scalar>
bool evalSplineWeights(Scalar min, Scalar max, size_t size, Scalar x, ssize_t &offset,
                       Scalar *weights, bool extrapolate = false) {
    /* Give up when given an out-of-range or NaN argument */
    if (!(x >= min && x <= max) && !extrapolate)
        return false;

    /* Transform 'x' so that nodes lie at integer positions */
    Scalar t = (x - min) * (size - 1) / (max - min);

    /* Find the index of the left node in the queried subinterval */
    size_t idx = std::max((size_t) 0, std::min((size_t) t, size - 2));

    /* Compute the relative position within the interval */
    t -= (Scalar) idx;
    Scalar t2 = t * t, t3 = t2 * t;

    /* Function value weights */
    weights[0] = 0.f;
    weights[1] =  2 * t3 - 3 * t2 + 1;
    weights[2] = -2 * t3 + 3 * t2;
    weights[3] = 0.f;
    offset = (ssize_t) idx - 1;

    /* Turn derivative weights into node weights using
       an appropriate chosen finite differences stencil */
    Float d0 = t3 - 2*t2 + t, d1 = t3 - t2;
    if (idx > 0) {
        weights[2] += d0 * 0.5f;
        weights[0] -= d0 * 0.5f;
    } else {
        weights[2] += d0;
        weights[1] -= d0;
    }

    if (idx + 2 < size) {
        weights[3] += d1 * 0.5f;
        weights[1] -= d1 * 0.5f;
    } else {
        weights[2] += d1;
        weights[1] -= d1;
    }

    return true;
}

/**
 * \brief Compute weights to perform a spline-interpolated lookup on a
 * \a non-uniformly sampled 1D function.
 *
 * The implementation relies on Catmull-Rom splines, i.e. it uses finite
 * differences to approximate the derivatives at the endpoints of each spline
 * segment. The resulting weights are identical those internally used by \ref
 * sample1D().
 *
 * \param nodes
 *      Array containing \c size non-uniformly spaced values denoting positions
 *      the where the function to be interpolated was evaluated. They must be
 *      provided in \a increasing order.
 * \param size
 *      Denotes the size of the \c nodes array
 * \param x
 *      Evaluation point
 * \param[out] weights
 *      Pointer to a weight array of size 4 that will be populated
 * \param[out] offset
 *      Offset into the function samples associated with weights[0]
 * \param extrapolate
 *      Extrapolate values when \c x is out of range? (default: \c false)
 * \remark
 *      The Python API lacks the \c size parameter, which is inferred
 *      automatically from the size of the input array. The \c offset
 *      and \c weights parameters are returned as the second and third
 *      elements of a triple.
 * \return
 *      \c true on success and \c false when <tt>extrapolate=false</tt>
 *      and \c x lies outside of [\c min, \c max]
 */
template <typename Scalar>
bool evalSplineWeights(const Scalar *nodes, size_t size, Scalar x, ssize_t &offset,
                       Scalar *weights, bool extrapolate = false) {
    /* Give up when given an out-of-range or NaN argument */
    if (!(x >= nodes[0] && x <= nodes[size-1]) && !extrapolate)
        return false;

    /* Find the index of the left node in the queried subinterval */
    size_t idx = math::findInterval(size, 
        [&](size_t i) { return nodes[i] <= x; }
    );

    Scalar x0 = nodes[idx],
           x1 = nodes[idx+1],
           width = x1 - x0;

    /* Compute the relative position within the interval and powers of 't' */
    Scalar t  = (x - x0) / width,
           t2 = t * t,
           t3 = t2 * t;

    /* Function value weights */
    weights[0] = 0.f;
    weights[1] = 2*t3 - 3*t2 + 1;
    weights[2] = -2*t3 + 3*t2;
    weights[3] = 0.f;

    offset = (ssize_t) idx - 1;
    
    /* Turn derivative weights into node weights using
       an appropriate chosen finite differences stencil */
    float d0 = t3 - 2*t2 + t, d1 = t3 - t2;
    if (idx > 0) {
        float factor = width / (nodes[idx+1]-nodes[idx-1]);
        weights[2] += d0 * factor;
        weights[0] -= d0 * factor;
    } else {
        weights[2] += d0;
        weights[1] -= d0;
    }

    if (idx + 2 < size) {
        float factor = width / (nodes[idx+2]-nodes[idx]);
        weights[3] += d1 * factor;
        weights[1] -= d1 * factor;
    } else {
        weights[2] += d1;
        weights[1] -= d1;
    }

    return true;
}

/**
 * \brief Evaluate a cubic spline interpolant of a uniformly sampled 2D function
 *
 * This implementation relies on a tensor product of Catmull-Rom splines, i.e.
 * it uses finite differences to approximate the derivatives for each dimension
 * at the endpoints of spline patches.
 *
 * \param nodes1
 *      Arrays containing \c size1 non-uniformly spaced values denoting
 *      positions the where the function to be interpolated was evaluated 
 *      on the \c X axis (in increasing order)
 * \param size1
 *      Denotes the size of the \c nodes1 array
 * \param nodes
 *      Arrays containing \c size2 non-uniformly spaced values denoting
 *      positions the where the function to be interpolated was evaluated 
 *      on the \c Y axis (in increasing order)
 * \param size2
 *      Denotes the size of the \c nodes2 array
 * \param values
 *      A 2D floating point array of <tt>size1*size2</tt> cells containing
 *      irregularly spaced evaluations of the function to be interpolated.
 *      Consecutive entries of this array correspond to increments in the \c X
 *      coordinate.
 * \param x
 *      \c X coordinate of the evaluation point
 * \param y
 *      \c Y coordinate of the evaluation point
 * \param extrapolate
 *      Extrapolate values when \c p is out of range? (default: \c false)
 * \remark
 *      The Python API lacks the \c size1 and \c size2 parameters, which are
 *      inferred automatically from the size of the input arrays.
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt> and
 *      <tt>(x,y)</tt> lies outside of the node range
 */
template <typename Scalar>
Scalar eval2D(const Scalar *nodes1, size_t size1, const Scalar *nodes2,
              size_t size2, const Scalar *values, Float x, Float y,
              bool extrapolate = false) {
    Scalar weights[2][4];
    ssize_t offset[2];

    /* Compute interpolation weights separately for each dimension */
    if (!evalSplineWeights(nodes1, size1, x, offset[0], weights[0], extrapolate) ||
        !evalSplineWeights(nodes2, size2, y, offset[1], weights[1], extrapolate))
        return 0.f;

    Scalar result = 0.f;
    for (int yi=0; yi<=3; ++yi) {
        Scalar wy = weights[1][yi];
        for (int xi=0; xi<=3; ++xi) {
            Scalar wxy = weights[0][xi] * wy;
            if (wxy == 0)
                continue;

            ssize_t pos = (offset[1] + y) * size1 + offset[0] + x;
            result += values[pos] * wxy;
        }
    }
    return result;
}

#undef SPLINE_DERIVATIVES_UNIFORM
#undef SPLINE_DERIVATIVES_NONUNIFORM

// -----------------------------------------------------------------------
/*! @} */

NAMESPACE_END(spline)
NAMESPACE_END(layer)
