#include <layer/quad.h>
#include <layer/vector.h>
#include "python.h"

void python_export_quad(py::module &m_) {
    /* quad.h bindings */
    py::module m = m_.def_submodule("quad", "Functions for numerical quadrature");
    
    m.def("gaussLegendre", [](int n) {
        std::vector<Float> nodes(n), weights(n);
        quad::gaussLegendre(n, nodes.data(), weights.data());
        return std::make_pair(py::array(n, nodes.data()), py::array(n, weights.data()));
    }, D(quad, gaussLegendre));

    m.def("gaussLobatto", [](int n) {
        std::vector<Float> nodes(n), weights(n);
        quad::gaussLobatto(n, nodes.data(), weights.data());
        return std::make_pair(py::array(n, nodes.data()), py::array(n, weights.data()));
    }, D(quad, gaussLobatto));

    m.def("compositeSimpson", [](int n) {
        std::vector<Float> nodes(n), weights(n);
        quad::compositeSimpson(n, nodes.data(), weights.data());
        return std::make_pair(py::array(n, nodes.data()), py::array(n, weights.data()));
    }, D(quad, compositeSimpson));
    
    m.def("compositeSimpson38", [](int n) {
        std::vector<Float> nodes(n), weights(n);
        quad::compositeSimpson38(n, nodes.data(), weights.data());
        return std::make_pair(py::array(n, nodes.data()), py::array(n, weights.data()));
    }, D(quad, compositeSimpson38));
}
