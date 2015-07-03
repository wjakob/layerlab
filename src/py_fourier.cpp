#include <layer/fourier.h>
#include <pybind/functional.h>

#include "python.h"

using namespace layer;
namespace py = pybind;

void python_export_fourier(py::module &m_) {
    py::module m= m_.def_submodule("fourier", "Functions for sampling and evaluating Fourier series");

    /* fourier.h bindings */
    m.def("filonIntegrate", [](const std::function<Float(Float)> &f, int nCoeffs, int nEvals, Float a, Float b) {
        std::vector<Float> coeffs(nCoeffs, 0);
        filonIntegrate(f, &coeffs[0], nCoeffs, nEvals, a, b);
        return coeffs;
    }, D(filonIntegrate), py::arg("f"), py::arg("nCoeffs"), py::arg("nEvals"), py::arg("a") = 0, py::arg("b") = math::Pi);

    m.def("convolveFourier", [](const std::vector<Float> &a, const std::vector<Float> &b) {
        std::vector<Float> c(a.size() + b.size() - 1);
        convolveFourier(&a[0], (int) a.size(), &b[0], (int) b.size(), &c[0]);
        return c;
    }, D(convolveFourier), py::arg("a"), py::arg("b"));

    m.def("evalFourier", [](const std::vector<float> &coeffs, py::array_dtype<Float> phi) {
        float *temp = fourier_aligned_alloca(coeffs.size() * sizeof(float));
        memcpy(temp, &coeffs[0], sizeof(float) * coeffs.size());
        return py::vectorize([&](Float phi) { return evalFourier(temp, coeffs.size(), phi); })(phi);
    }, D(evalFourier), py::arg("coeffs"), py::arg("phi"));

    m.def("evalFourier3", [](const std::vector<std::vector<float>> &coeffs, Float phi) {
        if (coeffs.size() != 3)
            throw std::runtime_error("Incompatible input");
        size_t size =
            std::max({ coeffs[0].size(), coeffs[1].size(), coeffs[0].size() });
        float *temp[3];
        for (int i=0; i<3; ++i) {
            temp[i] = fourier_aligned_alloca(size * sizeof(float));
            memcpy(temp[i], &coeffs[i][0], sizeof(float) * coeffs[i].size());
        }
        return evalFourier3(temp, coeffs[0].size(), phi);
    }, D(evalFourier3), py::arg("coeffs"), py::arg("phi"));
}
