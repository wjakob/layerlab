/*
    src/python.cpp -- Layer lab Python bindings

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#include <layer/common.h>
#include "python.h"

PYTHON_DECLARE(math);
PYTHON_DECLARE(vector);
PYTHON_DECLARE(spline);
PYTHON_DECLARE(fourier);
PYTHON_DECLARE(filesystem);
PYTHON_DECLARE(layer);
PYTHON_DECLARE(fresnel);
PYTHON_DECLARE(quad);

PYBIND11_PLUGIN(layerlab) {
    py::module m("layerlab", "Layer lab Python plugin");

    PYTHON_IMPORT(math);
    PYTHON_IMPORT(vector);
    PYTHON_IMPORT(spline);
    PYTHON_IMPORT(fourier);
    PYTHON_IMPORT(filesystem);
    PYTHON_IMPORT(layer);
    PYTHON_IMPORT(fresnel);
    PYTHON_IMPORT(quad);

    return m.ptr();
}
