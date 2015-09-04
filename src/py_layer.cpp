#include <layer/hg.h>
#include <layer/microfacet.h>
#include <layer/vector.h>
#include <layer/layer.h>
#include <layer/storage.h>

#include "python.h"

void python_export_layer(py::module &m) {
    /* hg.h bindings */
    m.def("hg", py::vectorize(&hg), D(hg), py::arg("mu_o"), py::arg("mu_i"), py::arg("g"),
          py::arg("phi_d"));

    m.def("hgFourierSeries", [](Float mu_o, Float mu_i, Float g, int kmax, Float relerr) {
        std::vector<Float> result;
        hgFourierSeries(mu_o, mu_i, g, kmax, relerr, result);
        return result;
    }, D(hgFourierSeries), py::arg("mu_o"), py::arg("mu_i"), py::arg("g"), py::arg("kmax"), py::arg("relerr"));

    /* microfacet.h bindings */
    m.def("smithG1", &smithG1, D(smithG1), py::arg("v"), py::arg("m"), py::arg("alpha"));
    m.def("microfacetNoExp", py::vectorize(&microfacetNoExp), D(microfacetNoExp), py::arg("mu_o"),
          py::arg("mu_i"), py::arg("eta"), py::arg("alpha"), py::arg("phi_d"));

    m.def("microfacet", py::vectorize(&microfacet), D(microfacet), py::arg("mu_o"),
          py::arg("mu_i"), py::arg("eta"), py::arg("alpha"), py::arg("phi_d"));

    m.def("microfacetNoExpFourierSeries",
          [](Float mu_o, Float mu_i, std::complex<Float> eta, Float alpha,
             int n, Float phiMax) {
          std::vector<Float> result;
          microfacetNoExpFourierSeries(mu_o, mu_i, eta, alpha, n, phiMax, result);
          return result;
    },  D(microfacetNoExpFourierSeries), py::arg("mu_o"), py::arg("mu_i"),
        py::arg("eta"), py::arg("alpha"), py::arg("n"), py::arg("phiMax"));

    m.def("microfacetFourierSeries",
          [](Float mu_o, Float mu_i, std::complex<Float> eta, Float alpha,
             int n, Float relerr) {
          std::vector<Float> result;
          microfacetFourierSeries(mu_o, mu_i, eta, alpha, n, relerr, result);
          return result;
    },  D(microfacetFourierSeries), py::arg("mu_o"), py::arg("mu_i"),
        py::arg("eta"), py::arg("alpha"), py::arg("n"), py::arg("relerr"));

    m.def("expCosFourierSeries", [](Float A, Float B, Float relerr) {
        std::vector<Float> result;
        expCosFourierSeries(A, B, relerr, result);
        return result;
    }, D(expCosFourierSeries), py::arg("A"), py::arg("B"), py::arg("relerr"));

    py::class_<LayerMode>(m, "LayerMode")
        .def(py::init<size_t>(), D(LayerMode, LayerMode))
        .def(py::init<>())
        .def("reverse", &LayerMode::reverse, D(LayerMode, reverse))
        .def("clear", &LayerMode::clear, D(LayerMode, clear))
        .def("nonZeros", &LayerMode::nonZeros, D(LayerMode, nonZeros))
        .def("__repr__", &LayerMode::toString, D(LayerMode, toString))
        .def_property_readonly("reflectionTop", [](const LayerMode &m) -> MatrixX { return m.reflectionTop; }, D(LayerMode, reflectionTop))
        .def_property_readonly("reflectionBottom", [](const LayerMode &m) -> MatrixX { return m.reflectionBottom; }, D(LayerMode, reflectionBottom))
        .def_property_readonly("transmissionTopBottom", [](const LayerMode &m) -> MatrixX { return m.transmissionTopBottom; }, D(LayerMode, transmissionTopBottom))
        .def_property_readonly("transmissionBottomTop", [](const LayerMode &m) -> MatrixX { return m.transmissionBottomTop; }, D(LayerMode, transmissionBottomTop));

    py::class_<Layer>(m, "Layer", D(Layer))
        .def(py::init<const VectorX &, const VectorX &, int>(), py::arg("nodes"), py::arg("weights"), py::arg("nFourierOrders") = 1)
        .def(py::init<const Layer&>())
        .def("reverse", &Layer::reverse, D(Layer, reverse))
        .def("clear", &Layer::clear, D(Layer, clear))
        .def("setDiffuse", &Layer::setDiffuse, D(Layer, setDiffuse), py::arg("albedo"))
        .def("setHenyeyGreenstein", &Layer::setHenyeyGreenstein, D(Layer, setHenyeyGreenstein), py::arg("albedo"), py::arg("g"))
        .def("setIsotropic", &Layer::setIsotropic, D(Layer, setIsotropic), py::arg("albedo"))
        .def("setVonMisesFisher", &Layer::setVonMisesFisher, D(Layer, setVonMisesFisher), py::arg("albedo"), py::arg("kappa"))
        .def("setMatusik", &Layer::setMatusik, D(Layer, setMatusik), py::arg("filename"), py::arg("channel"), py::arg("fourierOrders") = 0)
        .def("setMicrofacet", &Layer::setMicrofacet, D(Layer, setMicrofacet),
            py::arg("eta"), py::arg("alpha"), py::arg("conserveEnergy") = false,
            py::arg("fourierOrders") = 0)
        .def("__repr__", &Layer::toString, D(Layer, toString))
        .def_static("add", [](const Layer &l1, const Layer &l2, bool homogeneous) { Layer l3(l1.nodes(), l1.weights()); Layer::add(l1, l2, l3, homogeneous); return l3; }, D(Layer, add))
        .def_static("add", [](const Layer &l1, const Layer &l2) { Layer l3(l1.nodes(), l1.weights()); Layer::add(l1, l2, l3); return l3; })
        .def("addToTop", [](Layer &l1, const Layer &l2, bool homogeneous) { l1.addToTop(l2, homogeneous); }, D(Layer, addToTop))
        .def("addToTop", [](Layer &l1, const Layer &l2) { l1.addToTop(l2); })
        .def("addToBottom", [](Layer &l1, const Layer &l2, bool homogeneous) { l1.addToBottom(l2, homogeneous); }, D(Layer, addToBottom))
        .def("addToBottom", [](Layer &l1, const Layer &l2) { l1.addToBottom(l2); })
        .def("expand", &Layer::expand, D(Layer, expand))
        .def("eval", [](const Layer &l, py::array_dtype<Float> mu_o, py::array_dtype<Float> mu_i, py::array_dtype<Float> phi_d) {
            return py::vectorize([&l](Float mu_o, Float mu_i, Float phi_d) { return l.eval(mu_o, mu_i, phi_d); })(mu_o, mu_i, phi_d);
        }, py::arg("mu_o"), py::arg("mu_i"), py::arg("phi_d") = 0)
        .def("matrix", &Layer::matrix, py::arg("matrix") = 0)
        .def("__getitem__", [](Layer &m, size_t i) -> LayerMode& {
            if (i >= m.fourierOrders())
                throw py::index_error();
            return m[i];
         }, D(Layer, operator_array), py::return_value_policy::reference_internal)
        .def_property_readonly("resolution", &Layer::resolution, D(Layer, resolution))
        .def_property_readonly("fourierOrders", &Layer::fourierOrders, D(Layer, fourierOrders))
        .def_property_readonly("weights", &Layer::weights, D(Layer, weights))
        .def_property_readonly("nodes", &Layer::nodes, D(Layer, nodes));

    py::class_<BSDFStorage>(m, "BSDFStorage", D(BSDFStorage))
        .def(py::init<const fs::path &, bool>())
        .def(py::init<const fs::path &>())
        .def("close", &BSDFStorage::close, D(BSDFStorage, close))
        .def_property_readonly("maxOrder", &BSDFStorage::maxOrder, D(BSDFStorage, maxOrder))
        .def_property_readonly("channelCount", &BSDFStorage::channelCount, D(BSDFStorage, channelCount))
        .def_property_readonly("nodeCount", &BSDFStorage::nodeCount, D(BSDFStorage, nodeCount))
        .def_property_readonly("basisCount", &BSDFStorage::basisCount, D(BSDFStorage, basisCount))
        .def_property_readonly("parameterCount", &BSDFStorage::parameterCount, D(BSDFStorage, parameterCount))
        .def_property_readonly("size", &BSDFStorage::size, D(BSDFStorage, size))
        .def_property_readonly("metadata", &BSDFStorage::metadata, D(BSDFStorage, metadata))
        .def_property_readonly("eta", &BSDFStorage::eta, D(BSDFStorage, eta))
        .def_property_readonly("extrapolated", &BSDFStorage::extrapolated, D(BSDFStorage, extrapolated))
        .def("alpha", [](const BSDFStorage &m, size_t i) {
                 if (i >= 2)
                     throw py::index_error();
                 return m.alpha((int) i);
             }, D(BSDFStorage, alpha))
        .def("setAlpha", [](BSDFStorage &m, size_t i, float v) {
                 if (i >= 2)
                     throw py::index_error();
                 m.setAlpha((int) i, v);
             }, D(BSDFStorage, setAlpha))
        .def("parameterSampleCount", [](const BSDFStorage &m, size_t i) {
                 if (i >= m.parameterCount())
                     throw py::index_error();
                 return m.parameterSampleCount(i);
             }, D(BSDFStorage, parameterSampleCount))
        .def("parameterSamplePositions", [](const BSDFStorage &m, size_t i) {
                 if (i >= m.parameterCount())
                     throw py::index_error();
                 py::list list;
                 for (size_t j=0; j<m.parameterSampleCount(i); ++j)
                     list.append(py::float_(m.parameterSamplePositions(i)[j]));
                return list;
             }, D(BSDFStorage, parameterSampleCount))
        .def_static("fromLayer", [](const fs::path &path, const Layer *layer) {
                return BSDFStorage::fromLayer(path, layer);
            }, D(BSDFStorage, fromLayer))
        .def_static("fromLayerRGB", [](const fs::path &path, const Layer *layerR, const Layer *layerG, const Layer *layerB) {
                return BSDFStorage::fromLayerRGB(path, layerR, layerG, layerB);
            }, D(BSDFStorage, fromLayer))
        .def("__repr__", &BSDFStorage::toString);

    m.def("parameterHeuristicMicrofacet", parameterHeuristicMicrofacet,
            py::arg("alpha"), py::arg("eta"), D(parameterHeuristicMicrofacet));

    m.def("parameterHeuristicHG", parameterHeuristicHG,
            py::arg("g"), D(parameterHeuristicHG));
}
