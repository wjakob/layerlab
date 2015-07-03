import sys, os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

from layerlab import *
import numpy as np
from numpy.linalg import norm

values = np.array([1.0, 5.0, 3.0, 6.0])
nodes = np.array([5.0, 20/3, 25/3, 10.0])
eps = 1e-5

@pytest.mark.parametrize("input, expected1, expected2", 
    [ (4, 0, 0.328), (5, 1, 1), (6, 3.832, 3.832), (10, 6, 6), (11, 0, 6.36) ])
def test_eval1D_1(input, expected1, expected2):
    assert abs(spline.eval1D(nodes, values, input, False) - expected1) < eps
    assert abs(spline.eval1D(nodes, values, input, True ) - expected2) < eps
    assert abs(spline.eval1D(5, 10, values, input, False) - expected1) < eps
    assert abs(spline.eval1D(5, 10, values, input, True ) - expected2) < eps

def test_eval1D_2():
    input = np.array([4.0, 5.0, 6.0, 10.0, 11.0])
    expected1 = np.array([0, 1, 3.832, 6, 0])
    expected2 = np.array([0.328, 1, 3.832, 6, 6.36])
    assert norm(spline.eval1D(5, 10, values, input, False) - expected1) < eps
    assert norm(spline.eval1D(5, 10, values, input, True ) - expected2) < eps
    assert norm(spline.eval1D(nodes, values, input, False) - expected1) < eps
    assert norm(spline.eval1D(nodes, values, input, True ) - expected2) < eps

@pytest.mark.parametrize("input, expected1, expected2", 
    [ (4, 0, 0.328), (5, 1, 1), (6, 3.832, 3.832), (10, 6, 6), (11, 0, 6.36) ])
def test_eval1D_3(input, expected1, expected2):
    results = [
        spline.evalSplineWeights(nodes, input, False),
        spline.evalSplineWeights(5, 10, len(values), input, False),
        spline.evalSplineWeights(nodes, input, True),
        spline.evalSplineWeights(5, 10, len(values), input, True)
    ]
    for index, result in enumerate(results):
        value = 0
        for j in range(4):
            if result[2][j] != 0:
                value += values[result[1] + j] * result[2][j]
        if index < 2:
            assert abs(value - expected1) < eps
        else:
            assert abs(value - expected2) < eps

@pytest.mark.parametrize("input, expected", 
    [ (2, 4), (4.0, 4.644584273224155), (5.5, 5.3522) ])
def test_invert1D(input, expected):
    x = np.array([3.0, 4, 5, 6])
    f = np.array([1.0, 2, 5, 6])
    assert abs(spline.invert1D(x, f, input) - expected) < eps
    assert abs(spline.invert1D(3, 6, f, input) - expected) < eps

def test_integrate1D():
    expected = np.array([0, 5.416666666666669, 12.152777777777782, 19.305555555555557])
    assert norm(spline.integrate1D(5, 10, values) - expected) < eps
    assert norm(spline.integrate1D(nodes, values) - expected) < eps

@pytest.mark.parametrize("input, expected", 
    [ (0, 5), (0.124213, 6), (0.5, 7.57735), (1, 10) ])
def test_sample1D(input, expected):
    cdf = spline.integrate1D(nodes, values)
    assert abs(spline.sample1D(5, 10, values, cdf, input)[0] - expected) < eps
    assert abs(spline.sample1D(nodes, values, cdf, input)[0] - expected) < eps
