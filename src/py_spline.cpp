#include <layer/spline.h>
#include <layer/vector.h>
#include "python.h"

#define D(...) DOC(layer, __VA_ARGS__ )

using namespace layer::spline;

void python_export_spline(py::module &m) {
    /* spline.h bindings */
    py::module spline = m.def_submodule(
        "spline", "Functions for evaluating and sampling Catmull-Rom splines");

    spline.def("evalSpline",   evalSpline<Float>, D(spline, evalSpline));
    spline.def("evalSplineD", evalSplineD<Float>, D(spline, evalSplineD));
    spline.def("evalSplineI", evalSplineI<Float>, D(spline, evalSplineI));

    spline.def("eval1D", [](Float min, Float max, const VectorX &values,
                            Float x, bool extrapolate) {
        return eval1D(min, max, values.data(), values.size(), x, extrapolate);
    }, D(spline, eval1D));

    spline.def("eval1D", [](Float min, Float max, const VectorX &values,
                            const MatrixX &x, bool extrapolate) -> MatrixX {
        MatrixX result(x.rows(), x.cols());
        for (int i=0; i<x.rows(); ++i)
            for (int j=0; j<x.cols(); ++j)
                result(i, j) = eval1D(min, max, values.data(), values.size(), x(i, j), extrapolate);
        return result;
    }, D(spline, eval1D));

    spline.def("eval1D", [](VectorX &nodes, const VectorX &values,
                            Float x, bool extrapolate) {
        if (nodes.size() != values.size())
            throw std::runtime_error("'nodes' and 'values' must have a matching size!");
        return eval1D(nodes.data(), values.data(), nodes.size(), x, extrapolate);
    }, D(spline, eval1D, 2));

    spline.def("eval1D", [](VectorX &nodes, const VectorX &values,
                            const MatrixX &x, bool extrapolate) -> MatrixX {
        if (nodes.size() != values.size())
            throw std::runtime_error("'nodes' and 'values' must have a matching size!");
        MatrixX result(x.rows(), x.cols());
        for (int i=0; i<x.rows(); ++i)
            for (int j=0; j<x.cols(); ++j)
                result(i, j) = eval1D(nodes.data(), values.data(), nodes.size(), x(i, j), extrapolate);
        return result;
    }, D(spline, eval1D, 2));

    spline.def("invert1D", [](Float min, Float max, const VectorX &values,
                              Float y) {
        return invert1D(min, max, values.data(), values.size(), y);
    }, D(spline, invert1D));

    spline.def("invert1D", [](const VectorX &nodes, const VectorX &values,
                              Float y) {
        if (nodes.size() != values.size())
            throw std::runtime_error("'nodes' and 'values' must have a matching size!");
        return invert1D(nodes.data(), values.data(), values.size(), y);
    }, D(spline, invert1D, 2));

    spline.def("integrate1D", [](Float min, Float max, const VectorX &values) {
        VectorX result(values.size());
        integrate1D(min, max, values.data(), values.size(), result.data());
        return result;
    }, D(spline, integrate1D));

    spline.def("integrate1D", [](const VectorX &nodes, const VectorX &values) {
        if (nodes.size() != values.size())
            throw std::runtime_error("'nodes' and 'values' must have a matching size!");
        std::vector<Float> result(values.size());
        integrate1D(nodes.data(), values.data(), values.size(), result.data());
        return result;
    }, D(spline, integrate1D, 2));

    spline.def("sample1D", [](Float min, Float max, const VectorX &values,
                              const VectorX &cdf, Float sample) {
        if (values.size() != cdf.size())
            throw std::runtime_error("'values' and 'cdf' must have a matching size!");
        Float pos, fval, pdf;
        pos = sample1D(min, max, values.data(), cdf.data(), values.size(),
                       sample, &fval, &pdf);
        return std::make_tuple(pos, fval, pdf);
    }, D(spline, sample1D));

    spline.def("sample1D", [](const VectorX &nodes, const VectorX &values,
                              const VectorX &cdf, Float sample) {
        if (nodes.size() != values.size() || nodes.size() != cdf.size())
            throw std::runtime_error("'nodes', 'values', and 'cdf' must have a matching size!");
        Float pos, fval, pdf;
        pos = sample1D(nodes.data(), values.data(), cdf.data(), nodes.size(), sample, &fval, &pdf);
        return std::make_tuple(pos, fval, pdf);
    }, D(spline, sample1D, 2));

    spline.def("evalSplineWeights", [](Float min, Float max, size_t size, Float x, bool extrapolate) {
        Float weights[4] = { 0, 0, 0, 0 };
        ssize_t offset = 0;
        bool success = evalSplineWeights(min, max, size, x, offset, weights, extrapolate);
        return std::make_tuple(
            success, offset, std::make_tuple(weights[0], weights[1], weights[2], weights[3])
        );
    }, D(spline, evalSplineWeights));

    spline.def("evalSplineWeights", [](const VectorX &nodes, Float x, bool extrapolate) {
        Float weights[4] = { 0, 0, 0, 0};
        ssize_t offset = 0;
        bool success = evalSplineWeights(nodes.data(), nodes.size(), x, offset, weights, extrapolate);
        return std::make_tuple(
            success, offset, std::make_tuple(weights[0], weights[1], weights[2], weights[3])
        );
    }, D(spline, evalSplineWeights, 2));

    spline.def("eval2D", [](const VectorX &nodes1, const VectorX &nodes2,
                            const MatrixX &values, Float x, Float y,
                            bool extrapolate) {
        if (values.rows() != nodes1.size() || values.cols() != nodes2.size())
            throw std::runtime_error("'nodes' and 'values' must have a matching size!");

        return eval2D(nodes1.data(), nodes1.size(), nodes2.data(), nodes2.size(), values.data(), y, x, extrapolate);
    }, D(spline, eval2D));

    spline.def("eval2D", [](const VectorX &nodes1, const VectorX &nodes2,
                            const MatrixX &values, const MatrixX &x, const MatrixX &y,
                            bool extrapolate) {
        if (values.rows() != nodes1.size() || values.cols() != nodes2.size())
            throw std::runtime_error("'nodes' and 'values' must have a matching size!");
        if (x.rows() != nodes1.size() || x.cols() != y.size())
            throw std::runtime_error("'x' and 'y' must have a matching size!");

        MatrixX result(x.rows(), x.cols());
        for (int i=0; i<x.rows(); ++i)
            for (int j=0; j<x.cols(); ++j)
                result(i, j) = eval2D(
                    nodes1.data(), nodes1.size(), nodes2.data(), nodes2.size(),
                    values.data(), y(i, j), x(i, j), extrapolate);
        return result;

    }, D(spline, eval2D));
}
