#!/usr/bin/env python

"""
runs tests
"""
import pytest
from pytest import approx
from datetime import datetime

import flipchem


def test_MSIS():
    from flipchem import MSIS

    expected = (2278649876934.748, 48448644720266.11, 1861179167163.4343,
                3.2174722574090504e+16, 8.42931710665116e+16,
                1.468567250055534e+16, 335726315835103.75,
                5.575271108703837e-06, 1.2966025498855703e-23, 
                713.9961689976568, 465.49210257127993)


    date = datetime(2017,1,4,2)
    glat = 74.72955
    glon = -94.90576
    alt = 130.0

    msis = MSIS(date)
    outputs = msis.get_point(glat,glon,alt)
    for i in range(len(expected)):        
        assert(expected[i] == approx(outputs[i],nan_ok=True))

def test_compute_ion_neutral_collfreq():
    from flipchem import compute_ion_neutral_collfreq

    expected = [88.6752126659757,89.34495926234734,90.79723670237881,
                51.21018235127066,48.7255261408254]

    neutral_densities = (2278649876934.748, 48448644720266.11, 1861179167163.4343,
                         3.2174722574090504e+16, 8.42931710665116e+16,
                         1.468567250055534e+16, 335726315835103.75,
                         5.575271108703837e-06, 1.2966025498855703e-23, 
                         713.9961689976568, 465.49210257127993)
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

    neutral_densities = (2278649876934.748, 48448644720266.11, 1861179167163.4343,
                         3.2174722574090504e+16, 8.42931710665116e+16,
                         1.468567250055534e+16, 335726315835103.75,
                         5.575271108703837e-06, 1.2966025498855703e-23, 
                         713.9961689976568, 465.49210257127993)
    Te = 1000.0
    nu_en = compute_electron_neutral_collfreq(neutral_densities, Te) 

    assert(expected == approx(nu_en,nan_ok=True))


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])

