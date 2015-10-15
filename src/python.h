#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <layer/common.h>
#include <layer/vector.h>
#include "py_doc.h"

#define D(...) DOC(layer, __VA_ARGS__ )

#define PYTHON_DECLARE(name) \
    extern void python_export_##name(py::module &)
#define PYTHON_IMPORT(name) \
    python_export_##name(m)

using namespace layer;
namespace py = pybind11;
