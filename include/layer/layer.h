/*
    layer.h -- Microfacet BSDF evaluation routines

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/vector.h>
#include <Eigen/SparseCore>

NAMESPACE_BEGIN(layer)

#define ERROR_GOAL     1e-3f
#define DROP_THRESHOLD 1e-9f

typedef Eigen::SparseMatrix<Float, Eigen::ColMajor> MatrixS;

/**
 * \brief Helper class, which stores one Fourier order of a
 * layer scattering function
 */
struct LayerMode {
public:
    /**
     * \brief Create a new layer mode data structure
     *
     * \param size
     *     Resolution of the underlying matrix, must be a multiple of 2
     */
    LayerMode(size_t size = 0) {
        if (size % 2 == 1)
            throw std::runtime_error("LayerMode(): 'size' must be a multiple of 2");
        MatrixS::Index hSize = (MatrixS::Index) size / 2;
        reflectionTop.resize(hSize, hSize);
        reflectionBottom.resize(hSize, hSize);
        transmissionTopBottom.resize(hSize, hSize);
        transmissionBottomTop.resize(hSize, hSize);
    }

    /// Return the underlying matrix resolution (i.e. the number of discretizations in \mu)
    size_t resolution() const {
        return reflectionTop.rows() * 2;
    }

    /// Reset scattering matrix to that of a clear non-scattering layer
    void clear() {
        reflectionTop.setZero();
        reflectionBottom.setZero();
        transmissionTopBottom.setIdentity();
        transmissionBottomTop.setIdentity();
    }

    /// Reverse the layer
    void reverse() {
        reflectionTop.swap(reflectionBottom);
        transmissionTopBottom.swap(transmissionBottomTop);
    }

    /// Is this layer represented by a diagonal matrix?
    bool isDiagonal() const {
        int n = reflectionTop.rows();
        return
            reflectionTop.nonZeros() == 0 &&
            reflectionBottom.nonZeros() == 0 &&
            transmissionTopBottom.nonZeros() == n &&
            transmissionBottomTop.nonZeros() == n;
    }

    /// Access a matrix entry
    Float coeff(MatrixS::Index i, MatrixS::Index j) const {
        int n = reflectionTop.rows();

        if (i < n && j < n)
            return transmissionBottomTop.coeff(i, j);
        else if (i >= n && j >= n)
            return transmissionTopBottom.coeff(i-n, j-n);
        else if (i < n && j >= n)
            return reflectionTop.coeff(i, j-n);
        else if (i >= n && j < n)
            return reflectionBottom.coeff(i-n, j);
        else
            throw std::runtime_error("LayerMode::coeff(): out of bounds!");
    }

    /// Return the number of nonzero coefficients
    size_t nonZeros() const {
        return
            reflectionTop.nonZeros() +
            reflectionBottom.nonZeros() +
            transmissionTopBottom.nonZeros() +
            transmissionBottomTop.nonZeros();
    }

    void addScaled(const LayerMode &layer, Float scale) {
        reflectionTop         += layer.reflectionTop         * scale;
        reflectionBottom      += layer.reflectionBottom      * scale;
        transmissionTopBottom += layer.transmissionTopBottom * scale;
        transmissionBottomTop += layer.transmissionBottomTop * scale;
    }

    /// Return a human-readable summary
    std::string toString() const;
public:
    /// Reflection matrix on the top interface
    MatrixS reflectionTop;
    /// Reflection matrix on the bottom interface
    MatrixS reflectionBottom;
    /// Transmission matrix going from the top to the bottom interface
    MatrixS transmissionTopBottom;
    /// Transmission matrix going from the bottom to the top interface
    MatrixS transmissionBottomTop;
};


/**
 * \brief Discretized layer reflection model
 *
 * Describes the linear response to illumination that is incident along a
 * chosen set of zenith angles. The azimuthal dependence is modeled as
 * an even real Fourier transform. Each Fourier order is stored in a special
 * \c LayerMode data structure.
 */
class Layer {
public:
    /**
     * \brief Create a new layer with the given discretization in zenith angles
     * \param nodes
     *    A vector with the zenith angle cosines of the chosen discretization
     * \param weights
     *    Associated weights for each angle. Usually, the 'nodes' and 'weights'
     *    are generated using some kind of quadrature rule (e.g. Gauss-Legendre,
     *    Gauss-Lobatto, etc.)
     * \param nFourierOrders
     *    Specifies the number of coefficients to use for the azimuthal Fourier
     *    expansion (default: 1)
     */
    Layer(const VectorX &nodes, const VectorX &weights, size_t nFourierOrders = 1);

    /// Return the number of Fourier orders
    size_t fourierOrders() const { return m_modes.size(); }

    /// Return the number of nodes (i.e. the number of discretizations in \mu)
    size_t resolution() const { return m_nodes.size(); }

    /**
     * \brief Initialize the layer with a diffuse base layer
     *
     * \param albedo
     *    The layer's diffuse reflectance (in <tt>[0, 1]</tt>)
     */
    void setDiffuse(Float albedo);

    /**
     * \brief Initialize the layer with an isotropic phase function
     *
     * \param albedo
     *    The layer's single scattering albedo reflectance (in <tt>[0, 1]</tt>)
     */
    void setIsotropic(Float albedo);

    /**
     * \brief Initialize the layer with a Henyey-Greenstein phase function
     *
     * \param albedo
     *    The layer's single scattering albedo reflectance (in <tt>[0, 1]</tt>)
     * \param g
     *    The layer's HG anisotropy parameter(in <tt>[-1, 1]</tt>)
     */
    void setHenyeyGreenstein(Float albedo, Float g);

    /**
     * \brief Initialize the layer with a von Mises-Fisher phase function
     *
     * \param albedo
     *    The layer's single scattering albedo reflectance (in <tt>[0, 1]</tt>)
     * \param kappa
     *    The layer's kappa anisotropy parameter
     */
    void setVonMisesFisher(Float albedo, Float kappa);

    /**
     * \brief Initialize the layer with a microfacet model (dielectric or conductor)
     *
     * \param eta
     *    Relative index of refraction (complex)
     * \param alpha
     *    Beckmann roughness coefficient
     * \param conserveEnergy
     *    Correct for energy loss due to multiple scattering?
     *    Default: \c true
     * \param fourierOrders
     *    Number of fourier orders that should be used internally
     *    in the computation. Defaults to the value returned by 
     *    \ref fourierOrders()
     */
    void setMicrofacet(std::complex<Float> eta, Float alpha,
                       bool conserveEnergy = true,
                       size_t fourierOrders = 0);

    /**
     * \brief Combine two layers using the adding equations
     *
     * \param l1
     *    Input layer 1
     * \param l2
     *    Input layer 2
     * \param output
     *    Used to return the resulting layer
     * \param homogeneous
     *    When both layers are homogenous, (i.e. if their two sides are
     *    indistinguishable, this flag should be set to \c true to get a
     *    speed-up). Default: \c false
     * \remark
     *    In the Python API, the \c output parameter is directly returned
     */
    static void add(const Layer &l1, const Layer &l2, Layer &output,
                    bool homogeneous = false);

    /**
     * \brief Initialize the layer with a Matusik-style BRDF data file
     *
     * \param path
     *    Filename of the BRDF data file
     * \param channel
     *    Color channel to extract in <tt>[0..2]</tt>
     * \param fourierOrders
     *    Number of fourier orders that should be used internally
     *    in the computation. Defaults to the value returned by 
     *    \ref fourierOrders()
     */
    void setMatusik(const fs::path &path, int channel, int fourierOrders = -1);

    /**
     * \brief Append a layer above the current one
     *
     * This is just a convenience wrapper around \ref Layer::add()
     *
     * \param l
     *    The layer to be appended
     * \param homogeneous
     *    When the layers are homogenous, (i.e. if their two sides are
     *    indistinguishable, this flag should be set to \c true to get a
     *    speed-up). Default: \c false
     */
    void addToTop(const Layer &l, bool homogeneous = false) {
        Layer::add(l, *this, *this, homogeneous);
    }

    /**
     * \brief Append a layer below the current one
     *
     * This is just a convenience wrapper around \ref Layer::add()
     *
     * \param l
     *    The layer to be appended
     * \param homogeneous
     *    When the layers are homogenous, (i.e. if their two sides are
     *    indistinguishable, this flag should be set to \c true to get a
     *    speed-up). Default: \c false
     */
    void addToBottom(const Layer &l, bool homogeneous = false) {
        Layer::add(*this, l, *this, homogeneous);
    }

    /// Solve for the transport matrix of a layer with the given optical thickness (using Adding-Doubling)
    void expand(Float tau);

    /// Return the used integration weights
    const VectorX &weights() const { return m_weights; }

    /// Return the used integration nodes
    const VectorX &nodes() const { return m_nodes; }

    /// Reset scattering matrix to that of a clear non-scattering layer
    void clear();

    /// Reverse the layer
    void reverse();

    /// Look up a mode of the azimuthal Fourier expansion
    LayerMode &operator[](size_t i) { return m_modes[i]; }

    /// Look up a mode of the azimuthal Fourier expansion (const version)
    const LayerMode &operator[](size_t i) const { return m_modes[i]; }

    /**
     * \brief Return a dense representation of a mode's scattering matrix
     *
     * All integration weights and cosine foreshortening factors are
     * removed so that the coefficients can be interpreted as sampled
     * BSDF values
     */
    MatrixX matrix(size_t mode = 0) const;

    /// Evaluate the BSDF for a given pair of zenith angle cosines
    Float eval(Float mu_o, Float mu_i, Float phi_d = 0) const;

    /// Write the layer coefficients to a sparse file
    void write(const std::string &filename);

    /// Return a human-readable summary
    std::string toString() const;

protected:
    /// Helper struct for sparse matrix construction
    struct Quartet {
        uint32_t l, o, i;
        Float value;

        Quartet(uint32_t l, uint32_t o, uint32_t i, Float value)
         : l(l), o(o), i(i), value(value) { }
        Quartet(size_t l, size_t o, size_t i, Float value)
         : l((uint32_t) l), o((uint32_t) o), i((uint32_t) i), value(value) { }
    };

    /// Initialize from a list of quartets
    void setQuartets(const std::vector<Quartet> &quartets);
protected:
    /// Storage for all of the Fourier modes
    std::vector<LayerMode> m_modes;

    /// Integration nodes
    VectorX m_nodes;

    /// Integration weights
    VectorX m_weights;
};

/// Heuristic to guess a suitable number of parameters (Microfacet model)
extern std::pair<int, int> parameterHeuristicMicrofacet(Float alpha, std::complex<Float> &eta);

/// Heuristic to guess a suitable number of parameters (HG model)
extern std::pair<int, int> parameterHeuristicHG(Float g);

NAMESPACE_END(layer)
