/*
    quad.h -- Functions for numerical quadrature

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/math.h>
#include <cassert>
#include <stdexcept>
#include <iostream>

NAMESPACE_BEGIN(layer)
NAMESPACE_BEGIN(quad)

/**
 * \brief Computes the nodes and weights of a Gauss-Legendre quadrature
 * (aka "Gaussian quadrature") rule with the given number of evaluations.
 *
 * Integration is over the interval \f$[-1, 1]\f$. Gauss-Legendre quadrature
 * maximizes the order of exactly integrable polynomials achieves this up to
 * degree \f$2n-1\f$ (where \f$n\f$ is the number of function evaluations).
 *
 * This method is numerically well-behaved until about \f$n=200\f$
 * and then becomes progressively less accurate. It is generally not a
 * good idea to go much higher---in any case, a composite or
 * adaptive integration scheme will be superior for large \f$n\f$.
 *
 * \param n
 *     Desired number of evalution points
 * \param[out] nodes
 *     Length-\c n array used to store the nodes of the quadrature rule
 * \param[out] nodes
 *     Length-\c n array used to store the weights of the quadrature rule
 * \remark
 *     In the Python API, the \c nodes and \c weights field are returned
 *     as a tuple
 */
template <typename Scalar>
void gaussLegendre(int n, Scalar *nodes, Scalar *weights) {
    if (n-- < 1)
        throw std::runtime_error("gaussLegendre(): n must be >= 1");

    if (n == 0) {
        nodes[0] = 0;
        weights[0] = 2;
    } else if (n == 1) {
        nodes[0] = (Scalar) -std::sqrt(1.0 / 3.0);
        nodes[1] = -nodes[0];
        weights[0] = weights[1] = 1;
    }

    int m = (n+1)/2;
    for (int i=0; i<m; ++i) {
        /* Initial guess for this root using that of a Chebyshev polynomial */
        double x = -std::cos((double) (2*i + 1) / (double) (2*n + 2) * math::Pi_d);
        int it = 0;

        while (true) {
            if (++it > 20)
                throw std::runtime_error(
                    "gaussLobatto(" + std::to_string(n) +
                        "): did not converge after 20 iterations!");

            /* Search for the interior roots of P_{n+1}(x) using Newton's method. */
            std::pair<double, double> L = math::legendre_pd(n+1, x);
            double step = L.first / L.second;
            x -= step;

            if (std::abs(step) <= 4 * std::abs(x) * std::numeric_limits<double>::epsilon())
                break;
        }

        std::pair<double, double> L = math::legendre_pd(n+1, x);
        weights[i] = weights[n - i] =
            (Scalar)(2 / ((1 - x * x) * (L.second * L.second)));
        nodes[i] = (Scalar) x;
        nodes[n - i] = (Scalar) -x;
        assert(i == 0 || x > nodes[i-1]);
    }

    if ((n % 2) == 0) {
        std::pair<double, double> L = math::legendre_pd(n+1, 0.0);
        weights[n/2] = (double) (2 / (L.second*L.second));
        nodes[n/2] = 0;
    }
}

/**
 * \brief Computes the nodes and weights of a Gauss-Lobatto quadrature
 * rule with the given number of evaluations.
 *
 * Integration is over the interval \f$[-1, 1]\f$. Gauss-Lobatto quadrature
 * is preferable to Gauss-Legendre quadrature whenever the endpoints of the
 * integration domain should explicitly be included. It maximizes the order
 * of exactly integrable polynomials subject to this constraint and achieves
 * this up to degree \f$2n-3\f$ (where \f$n\f$ is the number of function
 * evaluations).
 *
 * This method is numerically well-behaved until about \f$n=200\f$
 * and then becomes progressively less accurate. It is generally not a
 * good idea to go much higher---in any case, a composite or
 * adaptive integration scheme will be superior for large \f$n\f$.
 *
 * \param n
 *     Desired number of evalution points
 * \param[out] nodes
 *     Length-\c n array used to store the nodes of the quadrature rule
 * \param[out] nodes
 *     Length-\c n array used to store the weights of the quadrature rule
 * \remark
 *     In the Python API, the \c nodes and \c weights field are returned
 *     as a tuple
 */
template <typename Scalar>
void gaussLobatto(int n, Scalar *nodes, Scalar *weights) {
    if (n-- < 2)
        throw std::runtime_error("gaussLobatto(): n must be >= 1");

    nodes[0] = -1;
    nodes[n] =  1;
    weights[0] = weights[n] = 2 / (Scalar) (n * (n+1));

    int m = (n+1)/2;
    for (int i=1; i<m; ++i) {
        /* Initial guess for this root -- see "On the Legendre-Gauss-Lobatto Points
           and Weights" by Seymor V. Parter, Journal of Sci. Comp., Vol. 14, 4, 1999 */

        double x = -std::cos((i + 0.25) * math::Pi_d / n - 3/(8*n*math::Pi_d * (i + 0.25)));
        int it = 0;

        while (true) {
            if (++it > 20)
                throw std::runtime_error(
                    "gaussLobatto(" + std::to_string(n) +
                        "): did not converge after 20 iterations!");

            /* Search for the interior roots of P_n'(x) using Newton's method. The same
               roots are also shared by P_{n+1}-P_{n-1}, which is nicer to evaluate. */

            std::pair<double, double> Q = math::legendre_pd_diff(n, x);
            double step = Q.first / Q.second;
            x -= step;

            if (std::abs(step) <= 4 * std::abs(x) * std::numeric_limits<double>::epsilon())
                break;
        }

        double Ln = math::legendre_p(n, x);
        weights[i] = weights[n - i] = (Scalar) (2 / ((n * (n + 1)) * Ln * Ln));
        nodes[i] = (Scalar) x;
        nodes[n - i] = (Scalar) -x;
        assert(x > nodes[i-1]);
    }

    if ((n % 2) == 0) {
        double Ln = math::legendre_p(n, 0.0);
        weights[n / 2] = (Scalar) (2 / ((n * (n + 1)) * Ln * Ln));
        nodes[n/2] = 0;
    }
}

/**
 * \brief Computes the nodes and weights of a composite Simpson quadrature
 * rule with the given number of evaluations.
 *
 * Integration is over the interval \f$[-1, 1]\f$, which will be split into
 * \f$(n-1) / 2\f$ sub-intervals with overlapping endpoints. A 3-point
 * Simpson rule is applied per interval, which is exact for polynomials of
 * degree three or less.
 *
 * \param n
 *     Desired number of evalution points. Must be an odd number bigger than 3.
 * \param[out] nodes
 *     Length-\c n array used to store the nodes of the quadrature rule
 * \param[out] nodes
 *     Length-\c n array used to store the weights of the quadrature rule
 * \remark
 *     In the Python API, the \c nodes and \c weights field are returned
 *     as a tuple
 */
template <typename Scalar>
void compositeSimpson(int n, Scalar *nodes, Scalar *weights) {
    if (n % 2 != 1 || n < 3)
        throw std::runtime_error("compositeSimpson(): n must be >= 3 and odd");

    n = (n - 1) / 2;

    Scalar h      = (Scalar) 2 / (Scalar) (2 * n),
           weight = h * (Scalar) (1.0 / 3.0);

    for (int i = 0; i < n; ++i) {
        Float x = -1 + h * (2*i);
        nodes[2*i]     = x;
        nodes[2*i+1]   = x+h;
        weights[2*i]   = (i == 0 ? 1 : 2) * weight;
        weights[2*i+1] = 4 * weight;
    }
 
    nodes[2*n] = 1;
    weights[2*n] = weight;

}

/**
 * \brief Computes the nodes and weights of a composite Simpson 3/8 quadrature
 * rule with the given number of evaluations.
 *
 * Integration is over the interval \f$[-1, 1]\f$, which will be split into
 * \f$(n-1) / 3\f$ sub-intervals with overlapping endpoints. A 4-point
 * Simpson rule is applied per interval, which is exact for polynomials of
 * degree four or less.
 *
 * \param n
 *     Desired number of evalution points. Must be an odd number bigger than 3.
 * \param[out] nodes
 *     Length-\c n array used to store the nodes of the quadrature rule
 * \param[out] nodes
 *     Length-\c n array used to store the weights of the quadrature rule
 * \remark
 *     In the Python API, the \c nodes and \c weights field are returned
 *     as a tuple
 */
template <typename Scalar>
void compositeSimpson38(int n, Scalar *nodes, Scalar *weights) {
    if ((n-1) % 3 != 0 || n < 4)
        throw std::runtime_error("compositeSimpson38(): n-1 must be divisible by 3");

    n = (n - 1) / 3;

    Scalar h      = (Scalar) 2 / (Scalar) (3 * n),
           weight = h * (Scalar) (3.0 / 8.0);

    for (int i = 0; i < n; ++i) {
        Float x = -1 + h * (3*i);
        nodes[3*i]     = x;
        nodes[3*i+1]   = x+h;
        nodes[3*i+2]   = x+2*h;
        weights[3*i]   = (i == 0 ? 1 : 2) * weight;
        weights[3*i+1] = 3 * weight;
        weights[3*i+2] = 3 * weight;
    }
 
    nodes[3*n] = 1;
    weights[3*n] = weight;
}

//! @}
// -----------------------------------------------------------------------

NAMESPACE_END(quad)
NAMESPACE_END(layer)
