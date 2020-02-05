#!/usr/bin/env python

"""
runs tests
"""
import pytest
from pytest import approx
from datetime import datetime

import flipchem


def test_Flipchem():
    from flipchem import Flipchem

    expected = (19.588119506835938,93.07436627254928,-22.65814170352014,
                459733062500.0,29630238281.25,10509367187.5,35145942.68798828,
                92176940.91796875,1110899.625,190275.46875,1)

    date = datetime(2017,1,4,2)
    fc = Flipchem(date)

    glat = 74.72955
    glon = -94.90576
    alt = 190.0
    ne = 5.0e11
    te = ti = 600.
    outputs = fc.get_point(glat,glon,alt,ne,te,ti)
    print(outputs)
    for i in range(len(expected)):        
        assert(expected[i] == approx(outputs[i],nan_ok=True))


def test_MSIS():
    from flipchem import MSIS

    expected = (48448644.72026611,32174722574.090504,84293171066.5116,
                14685672500.55534,335726315.83510375,5.575271108703837e-12,
                2278649.8769347477,1861179.1671634344,1.2966025498855705e-29,
                713.9961689976568,465.49210257127993)

    date = datetime(2017,1,4,2)
    glat = 74.72955
    glon = -94.90576
    alt = 130.0

    msis = MSIS(date)
    outputs = msis.get_point(glat,glon,alt)
    for i in range(len(expected)):        
        assert(expected[i] == approx(outputs[i],nan_ok=True))


def test_read_geophys():
    import flipchem

    expected = (71.0,74.0962962962963,[10.,4.,12.,2.,3.,18.,15.625])

    f107, f107a, ap = flipchem.read_geophys(datetime(2017,1,4,2))

    assert(expected[0] == approx(f107,nan_ok=True))
    assert(expected[1] == approx(f107a,nan_ok=True))
    for i in range(len(expected[2])):        
        assert(expected[2][i] == approx(ap[i],nan_ok=True))


def test_compute_ion_neutral_collfreq():
    from flipchem import compute_ion_neutral_collfreq

    expected = [88.1420437471394,88.85557877009926,90.46249148037654,50.89153579892685,48.42139426329943]

    neutral_densities = (2278649.8769347477,48448644.72026611,1861179.1671634344,32174722574.090504,84293171066.5116,14685672500.55534)
    ion_masses = [14.0,16.0,28.0,30.0,32.0]
    Ti = 1000.0
    Tn = 465.49210257127993
    nu_in = list()  
    for mass in ion_masses: 
        nu_in.append(compute_ion_neutral_collfreq(neutral_densities, Tn, mass, Ti)) 

    for i in range(len(expected)):        
        assert(expected[i] == approx(nu_in[i],nan_ok=True))


def test_compute_electron_neutral_collfreq():
    from flipchem import compute_electron_neutral_collfreq
    
    expected = 2052.242918871202

    neutral_densities = (2278649.8769347477,48448644.72026611,1861179.1671634344,32174722574.090504,84293171066.5116,14685672500.55534)
    Te = 1000.0
    nu_en = compute_electron_neutral_collfreq(neutral_densities, Te) 

    assert(expected == approx(nu_en,nan_ok=True))


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])

