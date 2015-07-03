#pragma once

#include <pybind/pybind.h>
#include <pybind/operators.h>
#include <pybind/complex.h>
#include <pybind/numpy.h>
#include <pybind/stl.h>
#include <layer/common.h>
#include <layer/vector.h>
#include "py_doc.h"

#define D(...) DOC(layer, __VA_ARGS__ )

#define PYTHON_DECLARE(name) \
    extern void python_export_##name(py::module &)
#define PYTHON_IMPORT(name) \
    python_export_##name(m)

using namespace layer;
namespace py = pybind;
