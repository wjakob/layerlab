/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <layer/storage.h>
#include <layer/spline.h>
#include <layer/fourier.h>
#include <layer/log.h>
#include <layer/simd.h>
#include <tbb/tbb.h>
#include <atomic>

#if defined(__MSVC__)
# include <intrin.h>
#else
# include <immintrin.h>
#endif

NAMESPACE_BEGIN(layer)

#define BSDF_STORAGE_HEADER_ID          "SCATFUN"
#define BSDF_STORAGE_VERSION            1
#define BSDF_STORAGE_FLAGS_EXTRAPOLATED 2
#define BSDF_STORAGE_FLAGS_BSDF         1
#define BSDF_STORAGE_HEADER_SIZE        64

static const float __basisCoeffsDefault[3] = { 1.0, 1.0, 1.0 };

BSDFStorage::BSDFStorage(const fs::path &filename, size_t nNodes, size_t nChannels,
            size_t nMaxOrder, size_t nCoeffs, size_t nBases, size_t nParameters,
            const size_t *paramSampleCounts, const float **paramSamplePositions, bool extrapolate,
            const std::string &metadata) : m_header(nullptr), m_reciprocals(nullptr),
            m_filename(filename), m_paramSamplePositionsNested(nullptr) {

    if (nChannels != 1 && nChannels != 3)
        Error("Only 1 and 3-channel files are supported!");

    if (extrapolate && nMaxOrder != 3)
        Error("Only three Fourier orders should be specified "
            "for the extrapolated storage format!");

    size_t nBasesPred = 1, nParameterValues = 0;

    for (size_t i=0; i<nParameters; ++i) {
        nParameterValues += paramSampleCounts[i];
        nBasesPred *= paramSampleCounts[i];
    }

    if (nBasesPred != nBases)
        Error("BSDFStorage::BSDFStorage(): provided an invalid number of basis functions");

    size_t size = BSDF_STORAGE_HEADER_SIZE + // Header
        sizeof(float)*nNodes +               // Node locations
        sizeof(uint32_t)*nParameters +       // Parameter sample counts
        sizeof(float)*nParameterValues +     // Parameter sample positions
        sizeof(float)*nNodes*nNodes*nBases + // CDF in \mu
        sizeof(OffsetType)*nNodes*nNodes*2 + // Offset + size table
        sizeof(float)*nCoeffs +              // Fourier coefficients
        metadata.size();                     // Metadata

    size_t uncompressedSize = size - sizeof(float)*nCoeffs
        + nNodes*nNodes*nChannels*nBases*nMaxOrder*sizeof(float);

    Log("Creating sparse BSDF storage file \"%s\":", filename.str());
    Log("  Discretizations in mu  : %d", nNodes);
    if (!extrapolate)
        Log("  Max. Fourier orders    : %d", nMaxOrder);
    else
        Log("  Harmonic extrapolation : yes");
    Log("  Color channels         : %d", nChannels);
    Log("  Textured parameters    : %d", nParameters);
    Log("  Basis functions        : %d", nBases);
    Log("  Uncompressed size      : %s", memString(uncompressedSize));
    Log("  Actual size            : %s (reduced to %.2f%%)", memString(size),
            100 * size / (Float) uncompressedSize);

    m_mmap = new MemoryMappedFile(filename, size);
    m_header = (Header *) m_mmap->data();

    const char *id = BSDF_STORAGE_HEADER_ID;

    const size_t len = strlen(BSDF_STORAGE_HEADER_ID);
    for (size_t i=0; i<len; ++i)
        m_header->identifier[i] = id[i];

    m_header->version = BSDF_STORAGE_VERSION;
    m_header->flags = 0;
    m_header->flags |= BSDF_STORAGE_FLAGS_BSDF;
    if (extrapolate)
        m_header->flags |= BSDF_STORAGE_FLAGS_EXTRAPOLATED;
    m_header->nNodes = (uint32_t) nNodes;
    m_header->nParameters = (uint16_t) nParameters;
    m_header->nMaxOrder = (uint32_t) nMaxOrder;
    m_header->nChannels = (uint32_t) nChannels;
    m_header->nBases = (uint32_t) nBases;
    m_header->nParameterValues = (uint16_t) nParameterValues;
    m_header->nCoeffs = (uint32_t) nCoeffs;
    m_header->nMetadataBytes = (uint32_t) metadata.size();
    m_header->eta = 1.0f; // default

    m_nodes = m_header->data;
    m_paramSampleCounts = (uint32_t *) (m_nodes + nNodes);
    m_paramSamplePositions = (float *) (m_paramSampleCounts + nParameters);
    m_cdfMu = m_paramSamplePositions + nParameterValues;
    m_offsetTable = (OffsetType *) (m_cdfMu + nNodes*nNodes*nBases);
    m_coeffs = (float *) (m_offsetTable + nNodes*nNodes*2);

    size_t idx = 0;
    m_paramSamplePositionsNested = new float*[nParameters];
    for (size_t i=0; i<nParameters; ++i) {
        m_paramSampleCounts[i] = (uint32_t) paramSampleCounts[i];
        m_paramSamplePositionsNested[i] = m_paramSamplePositions + idx;
        for (size_t j = 0; j<m_paramSampleCounts[i]; ++j)
            m_paramSamplePositions[idx++] = (float) paramSamplePositions[i][j];
    }

    memcpy(m_coeffs + nCoeffs, metadata.c_str(), metadata.size());
    m_metadata = metadata;

    int extra = LANE_WIDTH + LANE_WIDTH-1;
    m_reciprocals = (double *) simd::malloc((nMaxOrder+extra) * sizeof(double));
    memset(m_reciprocals, 0, sizeof(double) * (nMaxOrder+extra));
    m_reciprocals += LANE_WIDTH-1;
    for (uint32_t i=0; i<nMaxOrder+LANE_WIDTH; ++i)
        m_reciprocals[i] = 1.0 / (double) i;
}

BSDFStorage::BSDFStorage(const fs::path &filename, bool readOnly)
        : m_header(nullptr), m_reciprocals(nullptr), m_filename(filename),
          m_paramSamplePositionsNested(nullptr) {
    static_assert(sizeof(Header) == BSDF_STORAGE_HEADER_SIZE, "Header size mismatch!");

    m_mmap = new MemoryMappedFile(filename, readOnly);
    if (m_mmap->size() < sizeof(Header))
        Error("BSDF storage file \"%s\" has a truncated header!",
              filename.str());

    m_header = (Header *) m_mmap->data();
    const char *id = BSDF_STORAGE_HEADER_ID;
    const size_t len = strlen(BSDF_STORAGE_HEADER_ID);
    if (memcmp(id, m_header->identifier, len) != 0)
        Error("BSDF storage file \"%s\" has a corrupt header!",
              filename.str().c_str());

    size_t
        nNodes = m_header->nNodes,
        nMaxOrder = m_header->nMaxOrder,
        nChannels = m_header->nChannels,
        nBases = m_header->nBases,
        nParameters = m_header->nParameters,
        nCoeffs = m_header->nCoeffs,
        nParameterValues = m_header->nParameterValues,
        nMetadataBytes = m_header->nMetadataBytes;

    size_t size = BSDF_STORAGE_HEADER_SIZE + // Header
        sizeof(float)*nNodes +               // Node locations
        sizeof(uint32_t)*nParameters +       // Parameter sample counts
        sizeof(float)*nParameterValues +     // Parameter sample positions
        sizeof(float)*nNodes*nNodes*nBases + // CDF in \mu
        sizeof(OffsetType)*nNodes*nNodes*2 + // Offset + size table
        sizeof(float)*nCoeffs +              // Fourier coeff
        nMetadataBytes;                      // Metadata

    size_t uncompressedSize = size - sizeof(float)*nCoeffs
        + nNodes*nNodes*nChannels*nBases*nMaxOrder*sizeof(float);

    if (m_mmap->size() != size)
        Error("BSDF storage file \"%s\" has an invalid size! (it"
            " is potentially truncated)", filename.str());

    Log("Mapped sparse BSDF storage file \"%s\" into memory:", filename.str());
    Log("  Discretizations in mu  : %d", nNodes);
    if (!extrapolated())
        Log("  Max. Fourier orders    : %d", nMaxOrder);
    else
        Log("  Harmonic extrapolation : yes");
    Log("  Color channels         : %d", nChannels);
    Log("  Textured parameters    : %d", nParameters);
    Log("  Basis functions        : %d", nBases);
    Log("  Uncompressed size      : %s", memString(uncompressedSize));
    Log("  Actual size            : %s (reduced to %.2f%%)", memString(size),
            100 * size / (Float) uncompressedSize);

    m_nodes = m_header->data;
    m_paramSampleCounts = (uint32_t *) (m_nodes + nNodes);
    m_paramSamplePositions = (float *) (m_paramSampleCounts + nParameters);
    m_cdfMu = m_paramSamplePositions + nParameterValues;
    m_offsetTable = (OffsetType *) (m_cdfMu + nNodes*nNodes*nBases);
    m_coeffs = (float *) (m_offsetTable + nNodes*nNodes*2);

    size_t idx = 0;
    m_paramSamplePositionsNested = new float*[nParameters];
    for (size_t i=0; i<nParameters; ++i) {
        m_paramSamplePositionsNested[i] = m_paramSamplePositions + idx;
        idx += m_paramSampleCounts[i];
    }

    m_metadata.resize(nMetadataBytes);
    memcpy(&m_metadata[0], m_coeffs + nCoeffs, nMetadataBytes);

    int extra = LANE_WIDTH + LANE_WIDTH-1;
    m_reciprocals = (double *) simd::malloc((nMaxOrder+extra) * sizeof(double));
    memset(m_reciprocals, 0, sizeof(double) * (nMaxOrder+extra));
    m_reciprocals += LANE_WIDTH-1;
    for (uint32_t i=0; i<nMaxOrder+LANE_WIDTH; ++i)
        m_reciprocals[i] = 1.0 / (double) i;
}

BSDFStorage::~BSDFStorage() {
    if (m_mmap)
        delete m_mmap;
    if (m_reciprocals)
        simd::free(m_reciprocals - (LANE_WIDTH-1));
    if (m_paramSamplePositionsNested)
        delete[] m_paramSamplePositionsNested;
}


Float integrateCubicInterp1DN(size_t idx, const Float *nodes, const Float *values, size_t size) {
	Float f0       = values[idx],
	      f1       = values[idx+1],
	      width    = nodes[idx+1] - nodes[idx],
	      d0, d1;

	/* Approximate the derivatives */
	if (idx > 0)
		d0 = width * (f1 - values[idx-1]) / (nodes[idx+1] - nodes[idx-1]);
	else
		d0 = f1 - f0;

	if (idx + 2 < size)
		d1 = width * (values[idx+2] - f0) / (nodes[idx+2] - nodes[idx]);
	else
		d1 = f1 - f0;

	return ((d0-d1) * (Float) (1.0 / 12.0) + (f0+f1) * 0.5f) * width;
}


BSDFStorage *BSDFStorage::fromLayerGeneral(const fs::path &filename,
            const Layer **layers, size_t nChannels, size_t nBases, size_t nParameters,
            const size_t *paramSampleCounts, const float **paramSamplePositions,
            bool extrapolate, const std::string &metadata) {
    const Layer &layer = *layers[0];
    size_t n = layer.resolution(), h = n/2;

    /* Insert an explicit mu=0 node to simplify the evaluation / sampling code */
    VectorX nodes = (VectorX(n + 2) 
            << layer.nodes().head(h).reverse(),
               0, 0, layer.nodes().tail(h)).finished();

    Log("BSDFStorage::fromLayerGeneral(): merging %d layers into \"%s\" "
        "- analyzing sparsity pattern..", nBases * nChannels, filename.str());

    /* Do a huge matrix transpose */
    size_t maxCoeffs = nodes.size()*nodes.size()*nBases*nChannels*
        (extrapolate ? 3 : layer.fourierOrders());
    size_t nNodes = (size_t) nodes.size();
    std::atomic<size_t> m(0);

    OffsetType *offsetTable = new OffsetType[nodes.size()*nodes.size()*2];

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, nNodes),
        [&](const tbb::blocked_range<size_t> &range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                for (size_t o = 0; o < nNodes; ++o) {
                    MatrixS::Index ip, op;
                    size_t offset = (o + i * nNodes) * 2;

                    if (i == h || i == h+1 || o == h || o == h+1) {
                        offsetTable[offset + 0] = 0;
                        offsetTable[offset + 1] = 0;
                        continue;
                    }

                    ip = (MatrixS::Index) (i < h ? (h-i-1) : (i-2));
                    op = (MatrixS::Index) (o < h ? (h-o-1) : (o-2));

                    size_t nCoeffs = 0;
                    for (size_t basis=0; basis<nBases; ++basis) {
                        for (size_t ch=0; ch<nChannels; ++ch) {
                            size_t sparseSize = 0;
                            float ref = std::abs((float) (*layers[basis*nChannels+ch])[0].coeff(op, ip));
                            float partialSum = 0;
                            if (ref != 0) {
                                sparseSize = (size_t) layer.fourierOrders();
                                for (size_t j=(size_t) layer.fourierOrders()-1; j>=1; --j) {
                                    float value = (float) (*layers[basis*nChannels+ch])[j].coeff(op, ip);
                                    partialSum += std::abs(value);
                                    if (partialSum <= ref * ERROR_GOAL || value == 0)
                                        sparseSize = j;
                                }
                                nCoeffs = std::max(nCoeffs, sparseSize);
                            }
                        }
                    }
                    if (extrapolate && nCoeffs > 0)
                        nCoeffs = 3;

                    do {
                        size_t m_cur = m;
                        if (nCoeffs <= m_cur || m.compare_exchange_strong(m_cur, nCoeffs))
                            break;
                    } while(true);

                    offsetTable[offset + 0] = 0;
                    offsetTable[offset + 1] = (OffsetType) nCoeffs;
                }
            }
        });

    if (extrapolate)
        m = 3;

    /* Compute the offset table */
    size_t totalCoeffs = 0;
    for (size_t i=0; i<nNodes*nNodes; ++i) {
        offsetTable[2*i] = (OffsetType) totalCoeffs;
        totalCoeffs += offsetTable[2*i + 1] * nBases*nChannels;
    }

    Log("Done. Number of coeff: %d" " / %d" ", sparsity=%.2f%%",
        totalCoeffs, maxCoeffs, 100 * (Float) totalCoeffs / (Float) maxCoeffs);

    BSDFStorage *storage = new BSDFStorage(filename, nNodes, nChannels, m,
        totalCoeffs, nBases, nParameters, paramSampleCounts, paramSamplePositions,
        extrapolate, metadata);

    Log("Copying data into sparse BSDF file ..");
    for (size_t i=0; i<nNodes; ++i)
        storage->m_nodes[i] = (float) nodes[i];

    memcpy(storage->offsetTable(), offsetTable,
        nNodes*nNodes*2*sizeof(OffsetType));

    /* Do a huge matrix transpose */
    for (size_t i=0; i<nNodes; ++i) {
        for (size_t o=0; o<nNodes; ++o) {
            std::pair<float *, OffsetType> coeffsAndCount = storage->coeffAndCount(o, i);
            float *coeffs = coeffsAndCount.first;
            OffsetType size = coeffsAndCount.second;

            MatrixS::Index ip, op;

            if (i == h || o == h) {
                assert(size == 0);
                continue;
            }

            ip = (MatrixS::Index) (i < h ? (h-i-1) : (i-2));
            op = (MatrixS::Index) (o < h ? (h-o-1) : (o-2));

            float weight = (float) std::abs(nodes[o] / (math::Pi * nodes[i] * layer.weights()[ip]));


            if (nChannels == 1) {
                for (size_t basis=0; basis<nBases; ++basis) {
                    for (OffsetType j=0; j<size; ++j) {
                        float value = (float) (*layers[basis])[j].coeff(op, ip) *
                            weight * (j == 0 ? 0.5f : 1.0f);
                        if (!std::isfinite(value))
                            Warn("Encountered invalid data: %f", value);

                        *coeffs++ = value;
                    }
                }
            } else if (nChannels == 3) {
                float *coeffsY = coeffs;
                float *coeffsR = coeffsY + size*nBases;
                float *coeffsB = coeffsR + size*nBases;

                for (size_t basis=0; basis<nBases; ++basis) {
                    for (OffsetType j=0; j<size; ++j) {
                        float weight2 = weight * (j == 0 ? 0.5f : 1.0f);
                        float R = (float) (*layers[basis*nChannels+0])[j].coeff(op, ip) * weight2;
                        float G = (float) (*layers[basis*nChannels+1])[j].coeff(op, ip) * weight2;
                        float B = (float) (*layers[basis*nChannels+2])[j].coeff(op, ip) * weight2;

                        float Y = R * 0.212671f + G * 0.715160f + B * 0.072169f;
                        if (!std::isfinite(Y))
                            Warn("Encountered invalid data: %f", Y);

                        *coeffsY++ = Y; *coeffsR++ = R; *coeffsB++ = B;
                    }
                }
            }
        }
    }

    Log("Computing cumulative distributions for importance sampling ..");

    /* Create an importance sampling CDF */
    VectorX splineValues(nNodes), cdf(nNodes);
    for (size_t i = 0; i < nNodes; ++i) {
        for (size_t basis = 0; basis < nBases; ++basis) {
            for (size_t o = 0; o < nNodes; ++o) {
                auto c = storage->coeffAndCount(o, i, basis, 0);
                splineValues[o] = c.second > 0 ? *c.first : 0;
            }
            spline::integrate1D(nodes.data(), splineValues.data(), nNodes, cdf.data());
            float *cdf_ptr = storage->cdf(i) + basis;
            for (size_t k=0; k<nNodes; ++k)
                cdf_ptr[k * nBases] = (float) cdf[k];
        }
    }

#if 0
    if (extrapolate) {
        Log("Performing harmonic extrapolation ..");
        SAssert(totalCoeffs % 3 == 0);
        for (size_t i=0; i<totalCoeffs; i += 3)
            HarmonicExtrapolation::transform(storage->m_coeffs + i, storage->m_coeffs + i);
    }
#endif
    Log("BSDFStorage::fromLayerGeneral(): Done.");

    return storage;
}

size_t BSDFStorage::size() const {
    if (!m_mmap)
        return 0;
    return m_mmap->size();
}

bool BSDFStorage::extrapolated() const {
    return m_header->flags & BSDF_STORAGE_FLAGS_EXTRAPOLATED;
}

Color3 BSDFStorage::eval(Float mu_i, Float mu_o, Float phi_d, const float *basisCoeffs) const {
    if (!basisCoeffs) {
        assert(basisCount() == 1);
        basisCoeffs = __basisCoeffsDefault;
    }

    ssize_t knotOffsetO, knotOffsetI;
    float knotWeightsO[4], knotWeightsI[4];

    spline::evalSplineWeights(m_nodes, m_header->nNodes, (float) mu_o, knotOffsetO, knotWeightsO);
    spline::evalSplineWeights(m_nodes, m_header->nNodes, (float) mu_i, knotOffsetI, knotWeightsI);

    size_t nChannels = channelCount(), nBases = basisCount();
    OffsetType nCoeffs = 0;

    float *coeffs[3];
    for (size_t i=0; i<nChannels; ++i)
        coeffs[i] = fourier_aligned_alloca(maxOrder() * sizeof(float));

    for (int i=0; i<4; ++i) {
        for (int o=0; o<4; ++o) {
            float weight = knotWeightsO[o] * knotWeightsI[i];
            if (weight == 0)
                continue;

            std::pair<const float *, OffsetType> coeffAndCount =
                this->coeffAndCount(knotOffsetO + o, knotOffsetI + i);

            const float *source = coeffAndCount.first;
            OffsetType count = coeffAndCount.second;

            if (count == 0)
                continue;

            nCoeffs = std::max(nCoeffs, count);

            for (size_t channel=0; channel<nChannels; ++channel) {
                for (size_t basis=0; basis<nBases; ++basis) {
                    float interpWeight = weight * basisCoeffs[channel*nBases+basis];
                    if (interpWeight == 0) {
                        source += count;
                        continue;
                    }
                    float *target = coeffs[channel];
                    OffsetType remainder = count;

                    #if MTS_FOURIER_VECTORIZED == 1
                        /* Copy first (unaligned) element using scalar arithmetic */
                        *target++ += *source++ * interpWeight; --remainder;

                        /* Copy as many elements as possible using AVX */
                        __m256 weight_vec = _mm256_set1_ps(interpWeight);
                        OffsetType simdCount = remainder & ~7;
                        for (OffsetType k=0; k<simdCount; k += 8)
                            _mm256_store_ps(target+k, _mm256_add_ps(_mm256_load_ps(target+k),
                                _mm256_mul_ps(_mm256_loadu_ps(source+k), weight_vec)));

                        source += simdCount; target += simdCount; remainder -= simdCount;
                    #endif

                    for (OffsetType k=0; k<remainder; ++k)
                        *target++ += *source++ * interpWeight;
                }
            }
        }
    }

    Color3 result;
    if (nCoeffs == 0 || coeffs[0][0] == 0.0f) {
        result = Color3(0.0f);
#if 0
    } else if (m_header->flags & BSDF_STORAGE_FLAGS_EXTRAPOLATED) {
        float phi_d_sp = (float) phi_d;

        for (size_t ch=0; ch<nChannels; ++ch) {
            coeffs[ch][0] = std::max(0.0f, coeffs[ch][0]);
            coeffs[ch][1] = std::max(0.0f, std::min(1.0f, coeffs[ch][1]));
            coeffs[ch][2] = std::max(1e-6f,coeffs[ch][2]);
        }

        if (nChannels == 1) {
            result = Color3((Float) HarmonicExtrapolation::eval(coeffs[0], phi_d_sp));
        } else {
            Float Y = HarmonicExtrapolation::eval(coeffs[0], phi_d_sp);
            Float R = HarmonicExtrapolation::eval(coeffs[1], phi_d_sp);
            Float B = HarmonicExtrapolation::eval(coeffs[2], phi_d_sp);
            Float G = 1.39829f*Y - 0.100913f*B - 0.297375f*R;
            result.fromLinearRGB(R, G, B);
        }
#endif
    } else if (nChannels == 1) {
        result = Color3(std::max((Float) 0, (Float) evalFourier(coeffs[0], nCoeffs, phi_d)));
    } else {
        result = evalFourier3(coeffs, nCoeffs, phi_d);
    }
    result.clamp();

    return result;
}

Float BSDFStorage::pdf(Float mu_i, Float mu_o, Float phi_d, const float *basisCoeffs) const {
    if (!basisCoeffs) {
        assert(basisCount() == 1);
        basisCoeffs = __basisCoeffsDefault;
    }

    ssize_t knotOffsetO, knotOffsetI;
    float knotWeightsO[4], knotWeightsI[4];

    spline::evalSplineWeights(m_nodes, m_header->nNodes, (float) mu_o, knotOffsetO, knotWeightsO);
    spline::evalSplineWeights(m_nodes, m_header->nNodes, (float) mu_i, knotOffsetI, knotWeightsI);

    size_t nBases = basisCount();
    OffsetType nCoeffs = 0;

    float *coeffs = fourier_aligned_alloca(maxOrder() * sizeof(float));

    for (int i=0; i<4; ++i) {
        for (int o=0; o<4; ++o) {
            float weight = knotWeightsO[o] * knotWeightsI[i];
            if (weight == 0)
                continue;

            std::pair<const float *, OffsetType> coeffAndCount = 
                this->coeffAndCount(knotOffsetO + o, knotOffsetI + i);

            const float *source = coeffAndCount.first;
            OffsetType count = coeffAndCount.second;

            if (count == 0)
                continue;

            nCoeffs = std::max(nCoeffs, count);

            for (size_t basis=0; basis<nBases; ++basis) {
                float interpWeight = weight * basisCoeffs[basis];
                if (interpWeight == 0) {
                    source += count;
                    continue;
                }
                float *target = coeffs;
                OffsetType remainder = count;

                #if MTS_FOURIER_VECTORIZED == 1
                    /* Copy first (unaligned) element using scalar arithmetic */
                    *target++ += *source++ * interpWeight; --remainder;

                    /* Copy as many elements as possible using AVX */
                    __m256 weight_vec = _mm256_set1_ps(interpWeight);
                    OffsetType simdCount = remainder & ~7;
                    for (OffsetType k=0; k<simdCount; k += 8)
                        _mm256_store_ps(target+k, _mm256_add_ps(_mm256_load_ps(target+k),
                            _mm256_mul_ps(_mm256_loadu_ps(source+k), weight_vec)));

                    source += simdCount; target += simdCount; remainder -= simdCount;
                #endif

                for (OffsetType k=0; k<remainder; ++k)
                    *target++ += *source++ * interpWeight;
            }
        }
    }

    Float pdfMu = coeffs[0] / evalLatitudinalCDF(knotOffsetI,
            knotWeightsI, nodeCount()-1, basisCoeffs);

    if (nCoeffs == 0 || coeffs[0] == 0.0f) {
        return 0.0f;
#if 0
    } else if (m_header->flags & BSDF_STORAGE_FLAGS_EXTRAPOLATED) {
        coeffs[0] = std::max(0.0f, coeffs[0]);
        coeffs[1] = std::max(0.0f, std::min(1.0f, coeffs[1]));
        coeffs[2] = std::max(1e-6f,coeffs[2]);

        return HarmonicExtrapolation::pdf(coeffs, (float) phi_d) * pdfMu;
#endif
    } else {
        return std::max((Float) 0, pdfFourier(coeffs, nCoeffs, phi_d) * pdfMu);
    }
}

#if 0
Color3 BSDFStorage::sample(Float mu_i, Float &mu_o, Float &phi_d,
        Float &pdf, const Point2 &sample, const float *basisCoeffs) const {
    if (!basisCoeffs) {
        assert(basisCount() == 1);
        basisCoeffs = __basisCoeffsDefault;
    }

    size_t knotOffsetI, n = nodeCount();
    float knotWeightsI[4];

    /* Lookup spline nodes and weights for mu_i */
    spline::evalSplineWeights(m_nodes, m_header->nNodes, (float) mu_i, knotOffsetI, knotWeightsI);

    /* Account for energy loss */
    float normalization = evalLatitudinalCDF(knotOffsetI, knotWeightsI, n-1, basisCoeffs);
    float sample_y = (float) sample.y * normalization;

    /* Binary search for the spline segment containing the outgoing angle */
    size_t first = 0, len = n;
    while (len > 0) {
        ssize_t half = len >> 1, middle = first + half;
        if (evalLatitudinalCDF(knotOffsetI, knotWeightsI, middle, basisCoeffs) < sample_y) {
            first = middle + 1;
            len -= half + 1;
        } else {
            len = half;
        }
    }

    size_t index = std::min(n-2, std::max((size_t) 0, first-1));

    /* The spline segment to be sampled has been chosen. Determine the
       values of its nodes and then use the inversion method to sample
       an exact position within the segment */
    float cdf0  = evalLatitudinalCDF(knotOffsetI, knotWeightsI, index, basisCoeffs),
          cdf1  = evalLatitudinalCDF(knotOffsetI, knotWeightsI, index+1, basisCoeffs),
          f0    = evalLatitudinalAverage(knotOffsetI, knotWeightsI, index, basisCoeffs),
          f1    = evalLatitudinalAverage(knotOffsetI, knotWeightsI, index+1, basisCoeffs),
          width = m_nodes[index+1] - m_nodes[index],
          d0, d1;

    /* Catmull-Rom spline: approximate the derivatives at the endpoints
       using finite differences */
    if (index > 0) {
        d0 = width / (m_nodes[index+1] - m_nodes[index-1]) *
            (f1 - evalLatitudinalAverage(knotOffsetI, knotWeightsI, index-1, basisCoeffs));
    } else {
        d0 = f1 - f0;
    }

    if (index + 2 < n) {
        d1 = width / (m_nodes[index+2] - m_nodes[index]) *
            (evalLatitudinalAverage(knotOffsetI, knotWeightsI, index+2, basisCoeffs) - f0);
    } else {
        d1 = f1 - f0;
    }

    /* Bracketing interval and starting guess */
    float a = 0, c = 1, b;

    b          = (sample_y-cdf0) / (cdf1-cdf0);
    sample_y   = (sample_y-cdf0) / width;

    if (f0 != f1) /* Importance sample linear interpolant */
        b = (f0-math::safe_sqrt(f0*f0 + b * (f1*f1-f0*f0))) / (f0-f1);

    /* Invert CDF using Newton-Bisection */
    int it = 0;
    while (true) {
        if (!(b >= a && b <= c))
            b = 0.5f * (a + c);

        /* CDF and PDF in Horner form */
        float value = b*(f0 + b*(.5f*d0 + b*((1.0f/3.0f) * (-2*d0-d1)
            + f1 - f0 + b*(0.25f*(d0 + d1) + 0.5f * (f0 - f1))))) - sample_y;
        float deriv = f0 + b*(d0 + b*(-2*d0 - d1 + 3*(f1-f0) + b*(d0 + d1 + 2*(f0 - f1))));

        if (std::abs(value) < 1e-6f * deriv || ++it > 10) {
            mu_o = m_nodes[index] + width*b;
            break;
        }

        if (value > 0)
            c = b;
        else
            a = b;

        b -= value / deriv;
    }
    /* Outgoing zenith angle has been sampled -- interpolate
       Fourier coeff and sample the series */
    ssize_t knotOffsetO;
    float knotWeightsO[4];
    spline::evalSplineWeights(m_nodes, m_header->nNodes, (float) mu_o, knotOffsetO, knotWeightsO);

    size_t nChannels = channelCount(), nBases = basisCount();
    OffsetType nCoeffs = 0;

    float *coeffs[3];
    for (size_t i=0; i<nChannels; ++i)
        coeffs[i] = fourier_aligned_alloca(maxOrder() * sizeof(float));

    for (int i=0; i<4; ++i) {
        for (int o=0; o<4; ++o) {
            float weight = knotWeightsO[o] * knotWeightsI[i];
            if (weight == 0)
                continue;

            std::pair<const float *, OffsetType> coeffAndCount = coeffAndCount(knotOffsetO + o, knotOffsetI + i);

            const float *source = coeffAndCount.first;
            OffsetType count = coeffAndCount.second;

            if (count == 0)
                continue;

            nCoeffs = std::max(nCoeffs, count);

            for (size_t channel=0; channel<nChannels; ++channel) {
                for (size_t basis=0; basis<nBases; ++basis) {
                    float interpWeight = weight * basisCoeffs[channel*nBases+basis];
                    if (interpWeight == 0) {
                        source += count;
                        continue;
                    }
                    float *target = coeffs[channel];
                    OffsetType remainder = count;

                    #if MTS_FOURIER_VECTORIZED == 1
                        /* Copy first (unaligned) element using scalar arithmetic */
                        *target++ += *source++ * interpWeight; --remainder;

                        /* Copy as many elements as possible using AVX */
                        __m256 weight_vec = _mm256_set1_ps(interpWeight);
                        OffsetType simdCount = remainder & ~7;
                        for (OffsetType k=0; k<simdCount; k += 8)
                            _mm256_store_ps(target+k, _mm256_add_ps(_mm256_load_ps(target+k),
                                _mm256_mul_ps(_mm256_loadu_ps(source+k), weight_vec)));

                        source += simdCount; target += simdCount; remainder -= simdCount;
                    #endif

                    for (OffsetType k=0; k<remainder; ++k)
                        *target++ += *source++ * interpWeight;
                }
            }
        }
    }

    Float pdfMu = coeffs[0][0] / normalization;
    Float pdfPhi = 0;
    Color3 weight;

    if (coeffs[0][0] == 0) {
        weight = Color3(0.0f);
    } else if (m_header->flags & BSDF_STORAGE_FLAGS_EXTRAPOLATED) {
        for (size_t ch=0; ch<nChannels; ++ch) {
            coeffs[ch][0] = std::max(0.0f, coeffs[ch][0]);
            coeffs[ch][1] = std::max(0.0f, std::min(1.0f, coeffs[ch][1]));
            coeffs[ch][2] = std::max(1e-6f,coeffs[ch][2]);
        }

        Float phiWeight = HarmonicExtrapolation::sample(coeffs[0], phi_d, sample.x);
        float phi_d_sp = (float) phi_d;
        pdfPhi = HarmonicExtrapolation::pdf(coeffs[0], phi_d_sp);

        if (nChannels == 1) {
            weight = Color3(phiWeight * 2 * M_PI / pdfMu);
        } else {
            Float Y = HarmonicExtrapolation::eval(coeffs[0], phi_d_sp);
            Float R = HarmonicExtrapolation::eval(coeffs[1], phi_d_sp);
            Float B = HarmonicExtrapolation::eval(coeffs[2], phi_d_sp);
            Float G = 1.39829f*Y - 0.100913f*B - 0.297375f*R;
            weight.fromLinearRGB(R, G, B);
            weight /= pdfPhi * pdfMu;
        }
    } else if (nChannels == 1) {
        weight = Color3(std::max((Float) 0.0f, sampleFourier(coeffs[0], m_reciprocals,
            nCoeffs, (float) sample.x, pdfPhi, phi_d) / pdfMu));
    } else {
        weight = sampleFourier3(coeffs, m_reciprocals, nCoeffs,
            (float) sample.x, pdfPhi, phi_d) / pdfMu;
    }
    weight.clamp();

    pdf = std::max((Float) 0, pdfPhi * pdfMu);

    #if 1
        if (!std::isfinite(phi_d) || !std::isfinite(weight.getLuminance())) {
            cout << "Coeffs: ";
            for (size_t i=0; i<nCoeffs; ++i) {
                cout << coeffs[0][i] << ", ";
            }
            cout << "Sample=" << sample.toString() << endl;
            cout << "phi_d=" << phi_d << endl;
            cout << "pdfMu=" << pdfMu << endl;
            cout << "weight=" << weight.toString() << endl;
            cout << endl;
            Error("Internal error while sampling: phi_d: %f, weight: %s",
                phi_d, weight.toString().c_str());
        }
    #endif

    return weight;
}
#endif

void BSDFStorage::interpolateSeries(Float mu_i, Float mu_o, int basis, int channel, float *coeffs) const {
    ssize_t knotOffsetO, knotOffsetI;
    float knotWeightsO[4], knotWeightsI[4];

    spline::evalSplineWeights(m_nodes, m_header->nNodes, (float) mu_o, knotOffsetO, knotWeightsO);
    spline::evalSplineWeights(m_nodes, m_header->nNodes, (float) mu_i, knotOffsetI, knotWeightsI);

    memset(coeffs, 0, sizeof(float) * maxOrder());

    for (int i=0; i<4; ++i) {
        for (int o=0; o<4; ++o) {
            float weight = knotWeightsO[o] * knotWeightsI[i];
            if (weight == 0)
                continue;

            std::pair<const float *, OffsetType> coeffAndCount =
                this->coeffAndCount(knotOffsetO + o, knotOffsetI + i, basis, channel);

            const float *source = coeffAndCount.first;
            OffsetType count = coeffAndCount.second;
            float *target = coeffs;

            for (OffsetType k=0; k<count; ++k)
                *target++ += *source++ * weight;
        }
    }
}

std::string BSDFStorage::stats() const {
    std::ostringstream oss;

    size_t
        nNodes = m_header->nNodes,
        nMaxOrder = m_header->nMaxOrder,
        nChannels = m_header->nChannels,
        nBases = m_header->nBases,
        nParameters = m_header->nParameters,
        nCoeffs = m_header->nCoeffs,
        nParameterValues = m_header->nParameterValues,
        nMetadataBytes = m_header->nMetadataBytes;

    size_t size = BSDF_STORAGE_HEADER_SIZE + // Header
        sizeof(float)*nNodes +               // Node locations
        sizeof(uint32_t)*nParameters +       // Parameter sample counts
        sizeof(float)*nParameterValues +     // Parameter sample positions
        sizeof(float)*nNodes*nNodes*nBases + // CDF in \mu
        sizeof(OffsetType)*nNodes*nNodes*2 + // Offset + size table
        sizeof(float)*nCoeffs +              // Fourier coeff
        nMetadataBytes;                      // Metadata

    size_t uncompressedSize = size - sizeof(float)*nCoeffs
        + nNodes*nNodes*nChannels*nBases*nMaxOrder*sizeof(float);

    oss.precision(2);
    oss << " Discretizations in mu  : " << nNodes << std::endl;
    if (!extrapolated())
        oss << " Max. Fourier orders    : " << nMaxOrder << std::endl;
    else
        oss <<  " Harmonic extrapolation : yes" << std::endl;
    oss << " Color channels         : " << nChannels << std::endl;
    oss << " Textured parameters    : " << nParameters << std::endl;
    oss << " Basis functions        : " << nBases << std::endl;
    oss << " Uncompressed size      : " << memString(uncompressedSize) << std::endl;
    oss << " Actual size            : " << memString(size);
    oss << " (reduced to " << (100 * size / (Float) uncompressedSize) << "%)";
    return oss.str();
}

std::string BSDFStorage::toString() const {
    std::ostringstream oss;
    oss << "BSDFStorage[" << std::endl
        << "  mmap = " << m_mmap->toString() << "," << std::endl;
    if (m_header) {
        oss << "  nNodes = " << m_header->nNodes << "," << std::endl
            << "  nMaxOrder = " << m_header->nMaxOrder << "," << std::endl
            << "  nChannels = " << m_header->nChannels << "," << std::endl
            << "  nBases = " << m_header->nBases << "," << std::endl
            << "  eta = " << m_header->eta << std::endl;
    }
    oss << "]";
    return oss.str();
}

NAMESPACE_END(layer)
