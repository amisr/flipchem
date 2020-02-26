flipchem
========
.. image:: https://travis-ci.com/amisr/flipchem.svg?branch=master
    :target: https://travis-ci.com/amisr/flipchem

Overview
--------
`flipchem` provides a python wrapper of the flipchem ionospheric photochemistry model developed by Phil Richards [Richards2011]_. Specifically, this code wraps the version of flipchem that was used for [Richards2010]_. The model requires NRLMSIS-00 neutral density and the f107, f107a, and AP geophysical parameters, so both of these have been packaged with `flipchem`. NRLMSIS-00 is provided by wrapping the C version of the code written by Dominik Brodowski, which is based on the original Fortran version of the model [Picone2002]_.

Use Case
--------

The package is intended to be used for fitting Incoherent Scatter Radar (ISR) data, which involves fitting multiple gates of data at fixed time intervals. Also, the geophysical indicies are available with approximately daily time resolution. This is why the indicies are read when the `Flipchem` class is initialized and the `get_point*` functions do not accept date as an argument. The intended usage is:

1. initialize the Flipchem class for one time interval
2. call one of the get_point* functions for each of the range gates in the one time interval
3. initialize a new Flipchem class for the next time interval
4. goto 2

This facilitates usage of the multiprocessing library when processing ISR data, if each time interval is treated as a pool of asynchronous jobs.

Quick Start
-----------

Since this package is a wrapper around both C and Fortran source code, one must have a C and Fortran compiler installed and one must also install `numpy <https://numpy.readthedocs.io/en/latest/>`_ before attempting installation. `Numpy` is used during installation of `flipchem` to compile the Fortran via f2py [Oliphant2006]_.

Then, installation of `flipchem` can be accomplished using `pip`::

    pip install git+https://github.com/amisr/flipchem.git@v2020.2.0

And to make a profile of ion densities, one can try this::


    from datetime import datetime
    import numpy as np
    import flipchem
    date = datetime(2017,1,4,18)
    fc = flipchem.Flipchem(date)
    
    glat = 74.72955
    glon = -94.90576
    alts = np.linspace(90,350,num=100)
    nes = 5.0e11*np.exp(1-(alts-250)/70-np.exp(-(alts-250)/70))
    tes = 2500*(1-np.exp(-(alts-90)/70))+300
    tis = 1750*(1-np.exp(-(alts-90)/70))+300
    
    Op = np.zeros((alts.shape))
    O2p = np.zeros((alts.shape))
    NOp = np.zeros((alts.shape))
    N2p = np.zeros((alts.shape))
    Np = np.zeros((alts.shape))
    successes = np.zeros((alts.shape))
    
    for i,(alt,ne,te,ti) in enumerate(zip(alts,nes,tes,tis)):
        outputs = fc.get_point(glat,glon,alt,ne,te,ti)
        Op[i] = outputs[3]
        O2p[i] = outputs[4]
        NOp[i] = outputs[5]
        N2p[i] = outputs[6]
        Np[i] = outputs[7]
        successes[i] = outputs[-1]
    
    fig = pyplot.figure(figsize=(15,10))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(nes,alts,'--',label='$N_e$',lw=4)
    ax1.plot(Op,alts,label='$O^+$')
    ax1.plot(O2p,alts,label='$O_2^+$')
    ax1.plot(NOp,alts,label='$NO^+$')
    ax1.plot(N2p,alts,label='$N_2^+$')
    ax1.plot(Np,alts,label='$N^+$')
    ax1.set_xscale('log')
    ax1.set_xlim([1e9,1e12])
    ax1.set_ylabel('Altitude (km)')
    ax1.set_xlabel('Number Density (m$^{-3}$)')
    l = ax1.legend()
    
    ax2.plot(tes,alts,label='$T_e$')
    ax2.plot(tis,alts,label='$T_i$')
    ax2.set_xlabel('Temperature (K)')
    l = ax2.legend()

Documentation
-------------

You can finde more detailed documentation, including installation guides for Windows 10, macOS, and Linux here: https://flipchem.readthedocs.io

References
----------

.. [Oliphant2006] Oliphant, T. E. (2006). A guide to NumPy (Vol. 1). Trelgol Publishing USA.
.. [Picone2002] Picone, J. M., Hedin, A. E., Drob, D. P., and Aikin, A. C. (2002). NRLMSISE‐00 empirical model of the atmosphere: Statistical comparisons and scientific issues, J. Geophys. Res., 107(A12), 1468, doi:10.1029/2002JA009430. 
.. [Richards2010] Richards, P. G., Bilitza, D., and Voglozin, D. (2010), Ion density calculator (IDC): A new efficient model of ionospheric ion densities, Radio Sci., 45, RS5007, doi:10.1029/2009RS004332.
.. [Richards2011] Richards, P. G. (2011). Reexamination of ionospheric photochemistry, J. Geophys. Res., 116, A08307, doi:10.1029/2011JA016613.
