#!/usr/bin/env python

"""
runs tests
"""
import pytest
from pytest import approx
from datetime import datetime

import flipchem

def test_read_geophys():
    import flipchem

    expected = (71.0,74.0962962962963,[10.,4.,12.,2.,3.,18.,15.625])

    f107, f107a, ap = flipchem.read_geophys(datetime(2017,1,4,2))

    assert(expected[0] == approx(f107,nan_ok=True))
    assert(expected[1] == approx(f107a,nan_ok=True))
    for i in range(len(expected[2])):        
        assert(expected[2][i] == approx(ap[i],nan_ok=True))

if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])

