#include <layer/layer.h>
#include <layer/math.h>
#include <layer/hg.h>
#include <layer/microfacet.h>
#include <layer/fresnel.h>
#include <layer/spline.h>
#include <layer/log.h>
#include <filesystem/path.h>
#include <Eigen/SparseLU>
#include <Eigen/Geometry>
#include <tbb/tbb.h>

#if defined(HAVE_FFTW)
    #include <fftw3.h>
#endif

NAMESPACE_BEGIN(layer)

namespace {
    template <typename VectorType> MatrixS sparseDiagonal(const VectorType &vec) {
        MatrixS result(vec.size(), vec.size());
        for (MatrixS::Index i = 0; i < vec.size(); ++i)
            result.insert(i, i) = vec[i];
        result.makeCompressed();
        return result;
    }

    void sparsify(const MatrixX &dense, MatrixS &sparse) {
        sparse.setZero();

        for (MatrixX::Index j = 0; j < dense.cols(); ++j) {
            for (MatrixX::Index i = 0; i < dense.rows(); ++i) {
                Float value = dense.coeff(i, j);
                if (value != 0)
                    sparse.insert(i, j) = value;
            }
        }
        sparse.makeCompressed();
    }

    void scaleColumns(LayerMode &mode, const VectorX &d) {
        if ((size_t) d.size() != mode.resolution())
            throw std::runtime_error("scaleColumns(): size mismatch!");
        MatrixS scale = sparseDiagonal(d.head(d.size()/2));
        mode.transmissionBottomTop = mode.transmissionBottomTop * scale;
        mode.reflectionBottom = mode.reflectionBottom * scale;
        scale = sparseDiagonal(d.tail(d.size()/2));
        mode.transmissionTopBottom = mode.transmissionTopBottom * scale;
        mode.reflectionTop = mode.reflectionTop * scale;
    }

    void applySurfaceIntegrationWeights(Layer &layer) {
        for (size_t l=0; l<layer.fourierOrders(); ++l)
            scaleColumns(layer[l], layer.weights().cwiseProduct(
                                       layer.nodes().cwiseAbs()) *
                                       math::Pi * (l == 0 ? 2 : 1));
    }

    void applyMediumIntegrationWeights(Layer &layer) {
        for (size_t l=0; l<layer.fourierOrders(); ++l)
            scaleColumns(layer[l], layer.weights() * math::Pi * (l == 0 ? 2 : 1));
    }
};

Layer::Layer(const VectorX &nodes, const VectorX &weights, size_t nFourierOrders)
    : m_modes(nFourierOrders, LayerMode(nodes.size())), m_nodes(nodes), m_weights(weights) {
    if (nodes.size() < 2)
        throw std::runtime_error("Need at least 2 integration nodes!");
    else if (nodes.size() % 2 == 1)
        throw std::runtime_error("The number of integration nodes must be even!");
    for (int i=0; i<nodes.size(); ++i)
        if (nodes[i] == 0)
            throw std::runtime_error("The set of integrations includes mu=0 -- this is not allowed.");

    if (nodes[0] < nodes[1]) {
        size_t n = (size_t) nodes.size();
        /* Order integration weights so that they are usable for adding-doubling */
        m_weights.head(n/2).reverseInPlace();
        m_nodes.head(n/2).reverseInPlace();
    }
}

void Layer::setQuartets(const std::vector<Quartet> &quartets) {
    std::vector<std::vector<Eigen::Triplet<Float>>>
        tripletsTbt(fourierOrders()),
        tripletsTtb(fourierOrders()),
        tripletsRb(fourierOrders()),
        tripletsRt(fourierOrders());

    size_t approxSize = quartets.size() / (4 * fourierOrders());

    for (size_t i=0; i<fourierOrders(); ++i) {
        tripletsTbt[i].reserve(approxSize);
        tripletsTtb[i].reserve(approxSize);
        tripletsRb[i].reserve(approxSize);
        tripletsRt[i].reserve(approxSize);
    }

    size_t n = resolution() / 2;
    for (auto const &quartet: quartets) {
        typedef MatrixS::Index Index;
        if (quartet.o < n && quartet.i < n)
            tripletsTbt[quartet.l].emplace_back(Index(quartet.o), Index(quartet.i), quartet.value);
        else if (quartet.o >= n && quartet.i >= n)
            tripletsTtb[quartet.l].emplace_back(Index(quartet.o-n), Index(quartet.i-n), quartet.value);
        else if (quartet.o <n && quartet.i >= n)
            tripletsRt[quartet.l].emplace_back(Index(quartet.o), Index(quartet.i-n), quartet.value);
        else if (quartet.o >=n && quartet.i < n)
            tripletsRb[quartet.l].emplace_back(Index(quartet.o-n), Index(quartet.i), quartet.value);
        else
            throw std::runtime_error("Layer::setFromQuartets(): internal error!");
    }

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, fourierOrders(), 1),
        [&](const tbb::blocked_range<size_t> &range) {
            for (size_t l = range.begin(); l < range.end(); ++l) {
                m_modes[l].reflectionTop        .setFromTriplets(tripletsRt [l].begin(), tripletsRt [l].end());
                m_modes[l].reflectionBottom     .setFromTriplets(tripletsRb [l].begin(), tripletsRb [l].end());
                m_modes[l].transmissionTopBottom.setFromTriplets(tripletsTtb[l].begin(), tripletsTtb[l].end());
                m_modes[l].transmissionBottomTop.setFromTriplets(tripletsTbt[l].begin(), tripletsTbt[l].end());
            }
        }
    );
}

void Layer::setHenyeyGreenstein(Float albedo, Float g) {
    std::vector<Quartet> quartets;
    tbb::spin_mutex mutex;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, resolution()),
        [&](const tbb::blocked_range<size_t> &range) {
            std::vector<Quartet> quartetsLocal;
            quartetsLocal.reserve(fourierOrders() * resolution());
            std::vector<Float> result;
            for (size_t i = range.begin(); i < range.end(); ++i) {
                for (size_t o = 0; o <= i; ++o) {
                    hgFourierSeries(m_nodes[o], m_nodes[i], g, (int) fourierOrders(),
                                    ERROR_GOAL, result);
                    for (size_t l=0; l<std::min(fourierOrders(), result.size()); ++l) {
                        quartetsLocal.emplace_back(l, o, i, result[l] * albedo);
                        if (i != o)
                            quartetsLocal.emplace_back(l, i, o, result[l] * albedo);
                    }
                }
                tbb::spin_mutex::scoped_lock lock(mutex);
                quartets.insert(quartets.end(), quartetsLocal.begin(), quartetsLocal.end());
            }
        }
    );
    setQuartets(quartets);
    applyMediumIntegrationWeights(*this);
}

void Layer::setVonMisesFisher(Float albedo, Float kappa) {
    std::vector<Quartet> quartets;
    tbb::spin_mutex mutex;

    Float scale;
    if (kappa == 0)
        scale = albedo / (4 * math::Pi);
    else
        scale = albedo * kappa / (4 * math::Pi * std::sinh(kappa));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, resolution()),
        [&](const tbb::blocked_range<size_t> &range) {
            std::vector<Quartet> quartetsLocal;
            quartetsLocal.reserve(fourierOrders() * resolution());
            std::vector<Float> result;
            for (size_t i = range.begin(); i < range.end(); ++i) {
                Float mu_i = m_nodes[i];
                for (size_t o = 0; o <= i; ++o) {
                    Float mu_o = m_nodes[o];
                    Float A = kappa * mu_i * mu_o;
                    Float B = kappa * math::safe_sqrt((1 - mu_i * mu_i) *
                                                      (1 - mu_o * mu_o));
                    expCosFourierSeries(A, B, ERROR_GOAL, result);

                    for (size_t l=0; l<std::min(fourierOrders(), result.size()); ++l) {
                        quartetsLocal.emplace_back(l, o, i, result[l] * scale);
                        if (i != o)
                            quartetsLocal.emplace_back(l, i, o, result[l] * scale);
                    }
                }
                tbb::spin_mutex::scoped_lock lock(mutex);
                quartets.insert(quartets.end(), quartetsLocal.begin(), quartetsLocal.end());
            }
        }
    );

    setQuartets(quartets);
    applyMediumIntegrationWeights(*this);
}

void Layer::setDiffuse(Float albedo) {
    std::vector<Quartet> quartets;
    quartets.reserve(resolution() * resolution() / 2);

    size_t n = resolution(), h = n/2;
    for (size_t i=0; i<n; ++i) {
        for (size_t o=0; o<n; ++o) {
            if ((i < h && o >= h) || (o < h && i >= h))
                quartets.emplace_back(0, o, i, albedo * math::InvPi);
        }
    }

    setQuartets(quartets);
    applySurfaceIntegrationWeights(*this);
}

void Layer::setIsotropic(Float albedo) {
    std::vector<Quartet> quartets;
    quartets.reserve(resolution() * resolution());

    size_t n = resolution();
    for (size_t i=0; i<n; ++i)
        for (size_t o=0; o<n; ++o)
            quartets.emplace_back(0, o, i, albedo * math::InvFourPi);

    setQuartets(quartets);
    applyMediumIntegrationWeights(*this);
}

void Layer::setMicrofacet(std::complex<Float> eta, Float alpha, bool conserveEnergy,
                          size_t fourierOrdersTarget) {
    size_t n = resolution(), h = n/2;
    std::vector<Quartet> quartets;
    tbb::spin_mutex mutex;

    fourierOrdersTarget = std::max(fourierOrdersTarget, fourierOrders());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, resolution()),
        [&](const tbb::blocked_range<size_t> &range) {
            std::vector<Quartet> quartetsLocal;
            quartetsLocal.reserve(fourierOrdersTarget * resolution());
            std::vector<Float> result;
            for (size_t i = range.begin(); i < range.end(); ++i) {
                for (size_t o=0; o<n; ++o) {
                    /* Sign flip due to different convention (depth values
                     * increase opposite to the normal direction) */
                    microfacetFourierSeries(-m_nodes[o], -m_nodes[i], eta,
                                            alpha, fourierOrdersTarget,
                                            ERROR_GOAL, result);

                    for (size_t l=0; l<std::min(fourierOrders(), result.size()); ++l)
                        quartetsLocal.emplace_back(l, o, i, result[l]);
                }
            }
            tbb::spin_mutex::scoped_lock lock(mutex);
            quartets.insert(quartets.end(), quartetsLocal.begin(), quartetsLocal.end());
        }
    );

    setQuartets(quartets);

    /* Add a pseudo-diffuse term to capture lost energy */
    if (conserveEnergy && eta.imag() == 0) {
        /* Case 1: Dielectrics */
        VectorX W = m_weights.tail(h).cwiseProduct(m_nodes.tail(h)) * 2 * math::Pi;
        LayerMode &l = m_modes[0];

        VectorX Mb  = (W.asDiagonal() * MatrixX(l.reflectionBottom)).colwise().sum();
        VectorX Mt  = (W.asDiagonal() * MatrixX(l.reflectionTop)).colwise().sum();
        VectorX Mtb = (W.asDiagonal() * MatrixX(l.transmissionTopBottom)).colwise().sum();
        VectorX Mbt = (W.asDiagonal() * MatrixX(l.transmissionBottomTop)).colwise().sum();

        /* Determine how much energy we'd like to put into the transmission component
           (proportional to the current reflection/reflaction split) */
        VectorX Atb = (VectorX::Ones(h) - Mt - Mtb).cwiseProduct(Mtb.cwiseQuotient(Mt + Mtb));
        VectorX Abt = (VectorX::Ones(h) - Mb - Mbt).cwiseProduct(Mbt.cwiseQuotient(Mb + Mbt));
        Atb = Atb.cwiseMax(VectorX::Zero(h));
        Abt = Abt.cwiseMax(VectorX::Zero(h));

        /* Create a correction matrix which contains as much of the desired
           energy as possible, while maintaining symmetry and energy conservation */
        MatrixX Ctb = Abt*Atb.transpose() / std::max(W.dot(Abt), W.dot(Atb) / (eta.real()*eta.real()));
        MatrixX Cbt = Ctb.transpose() / (eta.real()*eta.real());

        sparsify(MatrixX(l.transmissionTopBottom) + Ctb, l.transmissionTopBottom);
        sparsify(MatrixX(l.transmissionBottomTop) + Cbt, l.transmissionBottomTop);

        /* Update missing energy terms */
        Mtb = (W.asDiagonal() * MatrixX(l.transmissionTopBottom)).colwise().sum();
        Mbt = (W.asDiagonal() * MatrixX(l.transmissionBottomTop)).colwise().sum();

        /* Put the rest of the missing energy into the reflection component */
        VectorX At = VectorX::Ones(h) - Mt - Mtb;
        VectorX Ab = VectorX::Ones(h) - Mb - Mbt;
        At = At.cwiseMax(VectorX::Zero(h));
        Ab = Ab.cwiseMax(VectorX::Zero(h));

        sparsify(MatrixX(l.reflectionTop)    + At*At.transpose() / W.dot(At), l.reflectionTop);
        sparsify(MatrixX(l.reflectionBottom) + Ab*Ab.transpose() / W.dot(Ab), l.reflectionBottom);
    } else if (conserveEnergy && eta.imag() != 0) {
        /* Case 2: Conductors */
        VectorX W = m_weights.tail(h).cwiseProduct(m_nodes.tail(h)) * 2 * math::Pi;

        /* Compute a reference matrix for a material *without* Fresnel effects */
        MatrixX refMatrix(n, n);

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, resolution()),
            [&](const tbb::blocked_range<size_t> &range) {
                std::vector<Float> result;
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    /* Parallel loop over 'i' */
                    for (size_t o = 0; o < n; ++o) {
                        microfacetFourierSeries(
                            -m_nodes[o], -m_nodes[i], std::complex<Float>(0.0f, 1.0f), alpha,
                            fourierOrdersTarget, ERROR_GOAL, result);
                        refMatrix(o, i) = result.size() > 0 ? result[0] : 0.0f;
                    }
                }
            }
        );

        MatrixX reflectionTopRef = MatrixX(refMatrix).block(0, h, h, h);

        VectorX Mt = VectorX::Ones(h) - (W.asDiagonal() * reflectionTopRef).colwise().sum().transpose();
        Mt = Mt.cwiseMax(VectorX::Zero(h));

        Float F = fresnelConductorIntegral(eta);
        Float E = 1 - W.dot(Mt) * math::InvPi;

        Float factor = F*E / (1-F*(1-E));

        MatrixX C = Mt * Mt.transpose() * (factor/ W.dot(Mt));

        sparsify(MatrixX(m_modes[0].reflectionTop) + C, m_modes[0].reflectionTop);
        sparsify(MatrixX(m_modes[0].reflectionBottom) + C, m_modes[0].reflectionBottom);
    }

    applySurfaceIntegrationWeights(*this);
}

void Layer::add(const Layer &layer1, const Layer &layer2, Layer &output, bool homogeneous) {
    if (output.resolution() != layer1.resolution() &&
        output.resolution() != layer2.resolution() &&
        output.fourierOrders() != layer1.fourierOrders() &&
        output.fourierOrders() != layer2.fourierOrders())
        throw std::runtime_error("Layer::addLayer(): incompatible sizes!");

    size_t n = output.resolution() / 2;
    MatrixS I((int) n, (int) n);
    I.setIdentity();

    /* Special case: it is possible to save quite a bit of computation when we
       know that both layers are homogeneous and of the same type */
    if (homogeneous) {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, layer1.fourierOrders(), 1),
            [&](const tbb::blocked_range<size_t> &range) {
                MatrixS Rb, Ttb;
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const LayerMode &l1 = layer1[i], &l2 = layer2[i];
                    LayerMode &lo = output[i];

                    /* Gain for downward radiation */
                    Eigen::SparseLU<MatrixS, Eigen::AMDOrdering<int>> G_tb;
                    G_tb.compute(I - l1.reflectionBottom * l2.reflectionTop);

                    /* Transmission at the bottom due to illumination at the top */
                    MatrixS result = G_tb.solve(l1.transmissionTopBottom);
                    Ttb = l2.transmissionTopBottom * result;

                    /* Reflection at the bottom */
                    MatrixS temp = l1.reflectionBottom * l2.transmissionBottomTop;
                    result = G_tb.solve(temp);
                    Rb = l2.reflectionBottom + l2.transmissionTopBottom * result;

                    #if defined(DROP_THRESHOLD)
                        Ttb.prune((Float) 1, (Float) DROP_THRESHOLD);
                        Rb.prune((Float) 1, (Float) DROP_THRESHOLD);
                    #endif

                    lo.transmissionTopBottom = Ttb;
                    lo.transmissionBottomTop = Ttb;
                    lo.reflectionTop = Rb;
                    lo.reflectionBottom = Rb;
                }
            }
        );
    } else {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, layer1.fourierOrders(), 1),
            [&](const tbb::blocked_range<size_t> &range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const LayerMode &l1 = layer1[i], &l2 = layer2[i];
                    LayerMode &lo = output[i];

                    /* Gain for downward radiation */
                    Eigen::SparseLU<MatrixS, Eigen::AMDOrdering<int>> G_tb;
                    G_tb.compute(I - l1.reflectionBottom * l2.reflectionTop);

                    /* Gain for upward radiation */
                    Eigen::SparseLU<MatrixS, Eigen::AMDOrdering<int>> G_bt;
                    G_bt.compute(I - l2.reflectionTop * l1.reflectionBottom);

                    /* Transmission at the bottom due to illumination at the top */
                    MatrixS result = G_tb.solve(l1.transmissionTopBottom);
                    MatrixS Ttb = l2.transmissionTopBottom * result;

                    /* Reflection at the bottom */
                    MatrixS temp = l1.reflectionBottom * l2.transmissionBottomTop;
                    result = G_tb.solve(temp);
                    MatrixS Rb = l2.reflectionBottom + l2.transmissionTopBottom * result;

                    /* Transmission at the top due to illumination at the bottom */
                    result = G_bt.solve(l2.transmissionBottomTop);
                    MatrixS Tbt = l1.transmissionBottomTop * result;

                    /* Reflection at the top */
                    temp = l2.reflectionTop * l1.transmissionTopBottom;
                    result = G_bt.solve(temp);
                    MatrixS Rt = l1.reflectionTop + l1.transmissionBottomTop * result;

                    #if defined(DROP_THRESHOLD)
                        Ttb.prune((Float) 1, (Float) DROP_THRESHOLD);
                        Tbt.prune((Float) 1, (Float) DROP_THRESHOLD);
                        Rb.prune((Float) 1, (Float) DROP_THRESHOLD);
                        Rt.prune((Float) 1, (Float) DROP_THRESHOLD);
                    #endif

                    lo.transmissionTopBottom = Ttb;
                    lo.transmissionBottomTop = Tbt;
                    lo.reflectionTop = Rt;
                    lo.reflectionBottom = Rb;
                }
            }
        );
    }
}

#if !defined(HAVE_FFTW)
void Layer::setMatusik(const fs::path &, int, int) {
    throw std::runtime_error("setMatusik(): You need to recompile with support for FFTW!");
}
#else
void Layer::setMatusik(const fs::path &path, int ch, int order) {
    if (ch < 0 || ch >= 3)
        throw std::runtime_error("Channel must be between 1 and 3");
    double scale[3] = { 1.0/1500.0, 1.15/1500.0, 1.66/1500.0 };

    FILE *f = fopen(path.str().c_str(), "rb");
    if (f == nullptr)
        throw std::runtime_error("I/O error: could not open file " + path.str());

    struct {
        int res_theta_h;
        int res_theta_d;
        int res_phi_d;
    } header;

    if (fread(&header, sizeof(int), 3, f) != 3)
        throw std::runtime_error("I/O error while loading header");

    Log("Loading Matusik-style BRDF data file \"%s\" (%ix%ix%i)",
        path.str(), header.res_theta_h, header.res_theta_d, header.res_phi_d);

    size_t nValues = 3 * header.res_theta_h * header.res_theta_d * header.res_phi_d;

    double *storage = new double[nValues];
    if (fread(storage, sizeof(double), nValues, f) != nValues)
        throw std::runtime_error("I/O error while loading file contents");

    fftw_plan_with_nthreads(1);

    order = std::max(order, (int) fourierOrders());
    int fftSize = order * 4;
    fftw_plan plan = fftw_plan_r2r_1d(fftSize, NULL, NULL, FFTW_REDFT00, FFTW_ESTIMATE);

    size_t n = resolution(), nEntries = n*n;
    tbb::spin_mutex mutex;
    std::vector<Quartet> quartets;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, nEntries, 1),
        [&](const tbb::blocked_range<size_t> &range) {
            for (size_t entry = range.begin(); entry < range.end(); ++entry) {
                int i = entry / n, o = entry % n;

                Float cosThetaI = m_nodes[i],
                      sinThetaI = std::sqrt(1-cosThetaI*cosThetaI),
                      cosThetaO = m_nodes[o],
                      sinThetaO = std::sqrt(1-cosThetaO*cosThetaO);

                Vector wi(sinThetaI, 0, cosThetaI);

                if (cosThetaI * cosThetaO > 0 || cosThetaI < 0)
                    continue;

                double *data = (double *) fftw_malloc(fftSize * sizeof(double));

                for (int j=0; j<fftSize; ++j) {
                    Float phi_d = M_PI * j/(Float) (fftSize-1),
                          cosPhi = std::cos(phi_d),
                          sinPhi = std::sin(phi_d);
                    Vector wo(-sinThetaO*cosPhi, -sinThetaO*sinPhi, -cosThetaO);

                    Vector half = (wi + wo).normalized();

                    Float theta_half = std::acos(half.z());
                    Float phi_half = std::atan2(half.y(), half.x());

                    Vector diff =
                          Eigen::AngleAxis<Float>(-theta_half, Vector(0, 1, 0)) *
                         (Eigen::AngleAxis<Float>(-phi_half, Vector(0, 0, 1)) * wi);

                    int theta_half_idx = std::min(std::max(0, (int) std::sqrt(
                        ((theta_half / (M_PI/2.0))*header.res_theta_h) * header.res_theta_h)), header.res_theta_h-1);

                    Float theta_diff = std::acos(diff.z());
                    Float phi_diff = std::atan2(diff.y(), diff.x());

                    if (phi_diff < 0)
                        phi_diff += M_PI;

                    int phi_diff_idx = std::min(std::max(0, int(phi_diff / M_PI * header.res_phi_d)), header.res_phi_d - 1);

                    int theta_diff_idx = std::min(std::max(0, int(theta_diff / (M_PI * 0.5) * header.res_theta_d)),
                        header.res_theta_d - 1);

                    int ind = phi_diff_idx +
                        theta_diff_idx * header.res_phi_d +
                        theta_half_idx * header.res_phi_d * header.res_theta_d;

                    data[j] = storage[ind+ch*header.res_theta_h*header.res_theta_d*header.res_phi_d] * scale[ch] * 2; /// XXX too dark?
                }
                double *spectrum = (double *) fftw_malloc(fftSize * sizeof(double));
                fftw_execute_r2r(plan, data, spectrum);

                for (int j=0; j<fftSize; ++j)
                    spectrum[j] /= (double) (fftSize-1);
                spectrum[0] /= 2;
                spectrum[fftSize-1] /= 2;

                double ref = std::abs(spectrum[0]);
                size_t sparseSize = 0;
                double partialSum = 0;
                if (ref != 0) {
                    sparseSize = fourierOrders();
                    for (size_t j= fourierOrders()-1; j>=1; --j) {
                        double value = (float) spectrum[j];
                        partialSum += std::abs(value);
                        if (partialSum <= ref * ERROR_GOAL)
                            sparseSize = j;
                    }
                }

                tbb::spin_mutex::scoped_lock lock(mutex);
                for (size_t l=0; l<sparseSize; ++l)
                    quartets.push_back(Quartet(l, o, i, (Float) spectrum[l]));
                fftw_free(spectrum);
                fftw_free(data);
            }
        }
    );

    fftw_destroy_plan(plan);
    delete[] storage;

    setQuartets(quartets);
    applySurfaceIntegrationWeights(*this);
}
#endif

void Layer::reverse() {
    for (auto &m: m_modes)
        m.reverse();
}

void Layer::clear() {
    for (auto &m: m_modes)
        m.clear();
}

void Layer::expand(Float target_tau) {
    /* Heuristic for choosing the initial width of a layer based on
       "Discrete Space Theory of Radiative Transfer" by Grant and Hunt
       Proc. R. Soc. London 1969 */
    Float tau = std::min(m_nodes.cwiseAbs().minCoeff() * 2, (Float) std::pow(2.0, -15.0));

    size_t doublings = (size_t) std::ceil(std::log(target_tau / tau) / std::log(2.0));
    tau = target_tau * std::pow((Float) 2.0f, -(Float) doublings);

    size_t n = resolution() / 2;

    MatrixS I((int) n, (int) n);
    I.setIdentity();

    MatrixS rowScale = sparseDiagonal(m_nodes.tail(n).cwiseInverse() * tau);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, fourierOrders(), 1),
        [&](const tbb::blocked_range<size_t> &range) {
            for (size_t i = range.begin(); i < range.end(); ++i) {
                LayerMode &mode = m_modes[i];
                MatrixS Rt = rowScale * mode.reflectionTop;
                MatrixS Ttb = I + rowScale * (mode.transmissionTopBottom - I);

                mode.reflectionTop = Rt;
                mode.reflectionBottom = Rt;
                mode.transmissionTopBottom = Ttb;
                mode.transmissionBottomTop = Ttb;
            }
        }
    );

    for (size_t i = 0; i < (size_t) doublings; ++i)
        addToTop(*this, true);
}

std::string Layer::toString() const {
    size_t nz = 0, nz_max = resolution() * resolution() * fourierOrders();
    for (auto const &mode: m_modes)
        nz += mode.nonZeros();
    std::ostringstream oss;
    oss.precision(2);
    oss << "Layer[resolution=" << resolution() << "x" << resolution()
        << ", fourierOrders=" << fourierOrders() << ", nonZeros=" << nz << "/"
        << nz_max << " (" << ((Float) nz / (Float) nz_max * 100) << "%)]";
    return oss.str();
}

MatrixX Layer::matrix(size_t l) const {
    if (l >= m_modes.size())
        throw std::runtime_error("Layer::matrix(): out of bounds!");
    size_t n = resolution(), h = n / 2;
    MatrixX M(n, n);

    M.topLeftCorner(h, h)     = m_modes[l].transmissionBottomTop;
    M.bottomRightCorner(h, h) = m_modes[l].transmissionTopBottom;
    M.topRightCorner(h, h)    = m_modes[l].reflectionTop;
    M.bottomLeftCorner(h, h)  = m_modes[l].reflectionBottom;

    for (size_t i = 0; i < n; ++i)
        M.col(i) /= m_weights[i] * std::abs(m_nodes[i]) * (l == 0 ? 2.f : 1.f) * math::Pi;

    M.block(0, 0, h, n) = M.block(0, 0, h, n).colwise().reverse().eval();
    M.block(0, 0, h, n) = M.block(0, 0, n, h).rowwise().reverse().eval();

    return M;
}

Float Layer::eval(Float mu_o, Float mu_i, Float phi_d) const {
    int n = m_nodes.size(), h = n / 2;
    ssize_t offset_o, offset_i;
    Float weights_o[4], weights_i[4];

    if (mu_o < 0 || (mu_o == 0 && mu_i > 0)) {
        spline::evalSplineWeights(m_nodes.data() + h, h, -mu_o, offset_o, weights_o, true);
    } else {
        spline::evalSplineWeights(m_nodes.data() + h, h, mu_o, offset_o, weights_o, true);
        offset_o += h;
    }

    if (mu_i < 0 || (mu_i == 0 && mu_o > 0)) {
        spline::evalSplineWeights(m_nodes.data() + h, h, -mu_i, offset_i, weights_i, true);
    } else {
        spline::evalSplineWeights(m_nodes.data() + h, h, mu_i, offset_i, weights_i, true);
        offset_i += h;
    }

    Float result = 0;
    for (size_t l=0; l<fourierOrders(); ++l) {
        Float sum = 0;
        for (int o = 0; o < 4; ++o) {
            for (int i = 0; i < 4; ++i) {
                Float weight = weights_o[o] * weights_i[i];
                if (weight == 0)
                    continue;
                weight /= 
                    std::abs(m_nodes[offset_i + i]) *
                    m_weights[offset_i + i] *
                    (l == 0 ? 2 : 1);

                sum += m_modes[l].coeff(offset_o + o, offset_i + i) * weight;
            }
        }
        result += sum * std::cos(phi_d * l);
    }

    return std::max((Float) 0.f, result / math::Pi);
}

std::string LayerMode::toString() const {
    size_t nz = nonZeros(), nz_max = resolution() * resolution();
    std::ostringstream oss;
    oss.precision(2);
    oss << "LayerMode[resolution=" << resolution() << "x" << resolution()
        << ", nonZeros=" << nz << "/" << nz_max << " ("
        << ((Float) nz / (Float) nz_max * 100) << "%)]";
    return oss.str();
}

std::pair<int, int> parameterHeuristicMicrofacet(Float alpha, std::complex<Float> &eta) {
    alpha = std::min(alpha, (Float) 1);
    if (eta.real() < 1 && eta.imag() == 0)
        eta = std::complex<Float>(1.f) / eta;

    static const Float c[][9] = {
        /* IOR    A_n      B_n     C_n       D_n      A_m      B_m      C_m      D_m                                 */
        {  0.0, 35.275,  14.136,  29.287,  1.8765,   39.814,  88.992, -98.998,  39.261  },  /* Generic conductor     */
        {  1.1, 256.47, -73.180,  99.807,  37.383,  110.782,  57.576,  94.725,  14.001  },  /* Dielectric, eta = 1.1 */
        {  1.3, 100.264, 28.187,  64.425,  14.850,   45.809,  17.785, -7.8543,  12.892  },  /* Dielectric, eta = 1.3 */
        {  1.5, 74.176,  27.470,  42.454,  9.6437,   31.700,  44.896, -45.016,  19.643  },  /* Dielectric, eta = 1.5 */
        {  1.7, 80.098,  17.016,  50.656,  7.2798,   46.549,  58.592, -73.585,  25.473  },  /* Dielectric, eta = 1.7 */
    };

    int i0 = 0, i1 = 0;

    if (eta.imag() == 0) { /* Dielectric case */
        for (int i=1; i<4; ++i) {
            if (eta.real() >= c[i][0] && eta.real() <= c[i+1][0]) {
                if (std::abs(eta.real()-c[i][0]) < 0.05f) {
                    i1 = i0 = i;
                } else if (std::abs(eta-c[i+1][0]) < 0.05f) {
                    i0 = i1 = i+1;
                } else {
                    i0 = i; i1 = i+1;
                }
            }
        }

        if (!i0)
            throw std::runtime_error("Index of refraction is out of bounds (must be between 1.1 and 1.7)!");
    }

    Float n0 = std::max(c[i0][1] + c[i0][2]*std::pow(std::log(alpha), (Float) 4)*alpha, c[i0][3]+c[i0][4]*std::pow(alpha, (Float) -1.2f));
    Float n1 = std::max(c[i1][1] + c[i1][2]*std::pow(std::log(alpha), (Float) 4)*alpha, c[i1][3]+c[i1][4]*std::pow(alpha, (Float) -1.2f));
    Float m0 = std::max(c[i0][5] + c[i0][6]*std::pow(std::log(alpha), (Float) 4)*alpha, c[i0][7]+c[i0][8]*std::pow(alpha, (Float) -1.2f));
    Float m1 = std::max(c[i1][5] + c[i1][6]*std::pow(std::log(alpha), (Float) 4)*alpha, c[i1][7]+c[i1][8]*std::pow(alpha, (Float) -1.2f));

    int n_i = (int) std::ceil(std::max(n0, n1));
    int m_i = (int) std::ceil(std::max(m0, m1));

    if (n_i % 2 == 1)
        n_i += 1;

    return std::make_pair(n_i, m_i);
}

std::pair<int, int> parameterHeuristicHG(Float g) {
    g = std::abs(g);
    Float m = 5.4f/(1.0f - g) - 1.3f;
    Float n = 8.6f/(1.0f - g) - 0.2f;
    int n_i = (int) std::ceil(n), m_i = (int) std::ceil(m);
    if (n_i % 2 == 1)
        n_i += 1;
    return std::make_pair(n_i, m_i);
}

NAMESPACE_END(layer)
