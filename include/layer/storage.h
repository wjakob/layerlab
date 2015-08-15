/*
    layer.h -- Sparse data structure for storing layers

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/layer.h>
#include <layer/color.h>
#include <layer/mmap.h>
#include <filesystem/path.h>

#if defined(_MSC_VER)
#pragma warning(disable: 4200) //  warning C4200 : nonstandard extension used : zero - sized array in struct / union
#endif

NAMESPACE_BEGIN(layer)

/*
 * \brief Storage class for isotropic BSDFs
 *
 * This class implements sparse storage support for isotropic BSDFs which are
 * point-sampled as a function of the incident and exitant zenith angles and
 * expanded into Fourier coefficients as a function of the azimuthal difference
 * angle.
 */
class BSDFStorage {
public:
    typedef uint32_t OffsetType;

    /// Map an existing BSDF storage file into memory
    BSDFStorage(const fs::path &filename, bool readOnly = true);

    /// Return the number of Fourier coefficients
    size_t maxOrder() const { return (size_t) m_header->nMaxOrder; }

    /// Return the number of color channels
    size_t channelCount() const { return (size_t) m_header->nChannels; }

    /// Return the resolution of the discretization in \mu_i and \mu_o
    size_t nodeCount() const { return (size_t) m_header->nNodes; }

    /// Return the number of basis functions stored in this file (usually just 1)
    size_t basisCount() const { return (size_t) m_header->nBases; }

    /// Return the number of model parameters
    size_t parameterCount() const { return (size_t) m_header->nParameters; }

    /// Return the number of samples associated with parameter \c i
    size_t parameterSampleCount(size_t i) const { return (size_t) m_paramSampleCounts[i]; }

    /// Return the sample positions associated with parameter \c i
    const float *parameterSamplePositions(size_t i) const { return m_paramSamplePositionsNested[i]; }

    /// Does this file store coefficients for the harmonic extrapolation-based model?
    bool extrapolated() const;

    /// Return the size of the underlying representation in bytes
    size_t size() const;

    /// Return metadata attached to the BSDF file (if any)
    const std::string &metadata() const { return m_metadata; }

    /// Return the relative index of refraction
    Float eta() const { return (Float) m_header->eta; }

    /// Set the relative index of refraction
    void setEta(Float eta) { m_header->eta = (float) eta; }

    /// Return the Beckmann-equivalent roughness (0: bottom, 1: top surface)
    Float alpha(int index) const { assert(index >= 0 && index <= 1); return (Float) m_header->alpha[index]; }

    /// Return the Beckmann-equivalent roughness (0: bottom, 1: top surface)
    void setAlpha(int index, Float alpha) { assert(index >= 0 && index <= 1); m_header->alpha[index] = (float) alpha; }

    /// Return the nodes of the underlying discretization in \mu_i and \mu_o
    const float *getNodes() const { return m_nodes; }

    /// Return a pointer to the coefficients of the CDF associated with the incident angle \c i
    float *cdf(size_t i) { return m_cdfMu + i*nodeCount()*basisCount(); }

    /// Return a pointer to the coefficients of the CDF associated with the incident angle \c i
    const float *cdf(int i) const { return m_cdfMu + i*nodeCount()*basisCount(); }

    /// Evaluate the model for the given values of \mu_i, \mu_o, and \phi_d
    Color3 eval(Float mu_i, Float mu_o, Float phi_d, const float *basisCoeffs = NULL) const;

    /// Evaluate the model for the given values of \mu_i, \mu_o, and \phi_d
    Float pdf(Float mu_i, Float mu_o, Float phi_d, const float *basisCoeffs = NULL) const;

    /// Importance sample the model
    Color3 sample(Float mu_i, Float &mu_o, Float &phi_d, Float &pdf,
                  const Point2 &sample,
                  const float *basisCoeffs = NULL) const;

    /// For debugging: return a Fourier series for the given parameters
    void interpolateSeries(Float mu_i, Float mu_o, int basis, int channel, float *coeffs) const;

    /// Forcefully release all resources
    void close() { delete m_mmap; m_mmap = NULL; m_header = NULL; m_coeffs = NULL; m_cdfMu = NULL; m_nodes = NULL; }

    /// Return a string representation
    std::string toString() const;

    /// Create a BSDF storage file from a Layer data structure (monochromatic)
    static BSDFStorage *fromLayer(const fs::path &filename, const Layer *layer,
            bool extrapolate = false, const std::string &metadata = "") {
        const Layer *layers[1] = { layer };
        return BSDFStorage::fromLayerGeneral(filename, layers, 1, 1, 0, NULL, NULL, extrapolate, metadata);
    }

    /// Create a BSDF storage file from three Layer data structures (RGB)
    static BSDFStorage *fromLayerRGB(const fs::path &filename, const Layer *layerR,
            const Layer *layerG, const Layer *layerB, bool extrapolate = false, const std::string &metadata = "") {
        const Layer *layers[3] = { layerR, layerG, layerB };
        return BSDFStorage::fromLayerGeneral(filename, layers, 3, 1, 0, NULL, NULL, extrapolate, metadata);
    }

    /// Create a BSDF storage file from three Layer data structures (most general interface)
    static BSDFStorage *fromLayerGeneral(const fs::path &filename,
            const Layer **layers, size_t nChannels, size_t nBases = 1, size_t nParameters = 0,
            const size_t *paramSampleCounts = NULL, const float **paramSamplePositions = NULL,
            bool extrapolate = false, const std::string &metadata = "");

    std::string stats() const;

    /// Virtual destructor
    virtual ~BSDFStorage();
protected:
    struct Header {
        uint8_t identifier[7];     // Set to 'SCATFUN'
        uint8_t version;           // Currently version is 1
        uint32_t flags;            // 0x01: file contains a BSDF, 0x02: uses harmonic extrapolation
        uint32_t nNodes;           // Number of samples in the elevational discretization

        uint32_t nCoeffs;          // Total number of Fourier series coefficients stored in the file
        uint32_t nMaxOrder;        // Coeff. count for the longest series occuring in the file
        uint32_t nChannels;        // Number of color channels (usually 1 or 3)
        uint32_t nBases;           // Number of BSDF basis functions (relevant for texturing)

        uint32_t nMetadataBytes;   // Size of descriptive metadata that follows the BSDF data
        uint32_t nParameters;      // Number of textured material parameters
        uint32_t nParameterValues; // Total number of BSDF samples for all textured parameters
        float eta;                 // Relative IOR through the material (eta(bottom) / eta(top))

        float alpha[2];            // Beckmann-equiv. roughness on the top (0) and bottom (1) side
        float unused[2];           // Unused fields to pad the header to 64 bytes

        float data[0];             // BSDF data starts here
    };

    /// Create a new BSDF storage file for the given amount of coefficients etc
    BSDFStorage(const fs::path &filename, size_t nNodes, size_t nChannels,
            size_t nMaxOrder, size_t nCoeffs, size_t nBases = 1,
            size_t nParameters = 0, const size_t *paramSampleCounts = NULL,
            const float **paramSamplePositions = NULL, bool extrapolate = false,
            const std::string &metadata = "");

    /// Return a posize_ter to the underlying sparse offset table
    OffsetType *offsetTable(size_t o = 0, size_t i = 0)
        { return m_offsetTable + 2*(o + i * nodeCount()); }

    /// Return a posize_ter to the underlying sparse offset table (const version)
    const OffsetType *offsetTable(size_t o = 0, size_t i = 0) const
        { return m_offsetTable + 2*(o + i * nodeCount()); }

    /// Return the sparse data offset of the given incident and exitant angle pair
    const float *coeff(size_t o, size_t i) const { return m_coeffs + offsetTable(o, i)[0]; }

    /// Return the sparse data offset of the given incident and exitant angle pair
    float *coeff(size_t o, size_t i) { return m_coeffs + offsetTable(o, i)[0]; }

    /// Return the sparse data offset and size of the given incident and exitant angle pair
    const float *coeff(size_t o, size_t i, size_t basis, size_t channel) const {
        const OffsetType *offsetPtr = offsetTable(o, i);
        OffsetType offset = offsetPtr[0], size = offsetPtr[1];
        return m_coeffs + offset + basis * size + basisCount()*size*channel;
    }

    /// Return the sparse data offset and size of the given incident and exitant angle pair
    float *coeff(size_t o, size_t i, size_t basis, size_t channel) {
        const OffsetType *offsetPtr = offsetTable(o, i);
        OffsetType offset = offsetPtr[0], size = offsetPtr[1];
        return m_coeffs + offset + basis * size + basisCount()*size*channel;
    }

    /// Return the sparse data size of the given incident and exitant angle pair
    OffsetType coeffCount(size_t o, size_t i) const {
        return offsetTable(o, i)[1];
    }

    /// Return the sparse data offset and size of the given incident and exitant angle pair
    std::pair<const float *, OffsetType> coeffAndCount(size_t o, size_t i) const {
        const OffsetType *offset = offsetTable(o, i);
        return std::make_pair(m_coeffs + offset[0], offset[1]);
    }

    /// Return the sparse data offset and size of the given incident and exitant angle pair
    std::pair<float *, OffsetType> coeffAndCount(size_t o, size_t i) {
        const OffsetType *offset = offsetTable(o, i);
        return std::make_pair(m_coeffs + offset[0], offset[1]);
    }

    /// Return the sparse data offset and size of the given incident and exitant angle pair
    std::pair<const float *, OffsetType> coeffAndCount(size_t o, size_t i, size_t basis) const {
        const OffsetType *offsetPtr = offsetTable(o, i);
        OffsetType offset = offsetPtr[0], size = offsetPtr[1];
        return std::make_pair(m_coeffs + offset + basis * size, size);
    }

    /// Return the sparse data offset and size of the given incident and exitant angle pair
    std::pair<float *, OffsetType> coeffAndCount(size_t o, size_t i, size_t basis) {
        const OffsetType *offsetPtr = offsetTable(o, i);
        OffsetType offset = offsetPtr[0], size = offsetPtr[1];
        return std::make_pair(m_coeffs + offset + basis * size, size);
    }

    /// Return the sparse data offset and size of the given incident and exitant angle pair
    std::pair<const float *, OffsetType> coeffAndCount(size_t o, size_t i, size_t basis, size_t channel) const {
        const OffsetType *offsetPtr = offsetTable(o, i);
        OffsetType offset = offsetPtr[0], size = offsetPtr[1];
        return std::make_pair(m_coeffs + offset + basis * size + basisCount()*size*channel, size);
    }

    /// Return the sparse data offset and size of the given incident and exitant angle pair
    std::pair<float *, OffsetType> coeffAndCount(size_t o, size_t i, size_t basis, size_t channel) {
        const OffsetType *offsetPtr = offsetTable(o, i);
        OffsetType offset = offsetPtr[0], size = offsetPtr[1];
        return std::make_pair(m_coeffs + offset + basis * size + basisCount()*size*channel, size);
    }

    /// Evaluate the discrete CDF that is used to sample a zenith angle spline segment
    float evalLatitudinalCDF(size_t knotOffset, float *knotWeights, size_t index, const float *basisCoeffs) const {
        const size_t n = nodeCount(), m = basisCount();
        const float *cdf = m_cdfMu + (knotOffset*n + index) * m;
        float result = 0;

        for (size_t i=0; i<4; ++i) {
            float weight = knotWeights[i];
            if (weight != 0)
                for (size_t basis=0; basis<m; ++basis)
                    result += cdf[i*n*m+basis] * weight * basisCoeffs[basis];
        }

        return result;
    }

    /// Evaluate the zeroeth-order Fourier coefficient given size_terpolation weights
    float evalLatitudinalAverage(size_t knotOffset, float *knotWeights, size_t index, const float *basisCoeffs) const {
        float result = 0.0f;
        for (size_t i=0; i<4; ++i) {
            float weight = knotWeights[i];
            if (weight == 0)
                continue;
            std::pair<const float *, OffsetType> coeffAndCount =
                this->coeffAndCount(index, knotOffset + i);
            const OffsetType count = coeffAndCount.second;
            if (!count)
                continue;

            const float *coeff = coeffAndCount.first;
            for (size_t basis=0; basis<basisCount(); ++basis)
                result += weight * basisCoeffs[basis] * coeff[basis*count];
        }

        return result;
    }
protected:
    MemoryMappedFile *m_mmap;
    Header *m_header;
    float *m_nodes;
    float *m_cdfMu;
    OffsetType *m_offsetTable;
    float *m_coeffs;
    double *m_reciprocals;
    uint32_t *m_paramSampleCounts;
    fs::path m_filename;
    float *m_paramSamplePositions;
    float **m_paramSamplePositionsNested;
    std::string m_metadata;
};

NAMESPACE_END(layer)
