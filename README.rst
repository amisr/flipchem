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

Installation
------------

Since this package is a wrapper around both C and Fortran source code, one must install `numpy <https://numpy.readthedocs.io/en/latest/>`_ before attempting installation. `Numpy` is used during installation of `flipchem` to compile the Fortran [Oliphant2006]_. 

Installation of `flipchem` can be accomplished using `pip`::

    pip install git+https://github.com/amisr/flipchem.git@v2020.1.1


Example Usage
-------------

First, it is very important to update the geophysical parameters files. There is an archive of these files available at `<https://amisr.com/geophys_params/>`_. They are organized by year, so once updated, you only need to update the current year. To update all years::

    import flipchem
    flipchem.update_geophys(year='all')


and to update a specific year, say 1991::

    import flipchem
    flipchem.update_geophys(year=1991)


Get the O+, O2+, NO+, N2+, and N+ ion densities at 130km above the RISR-N ISR::

    from datetime import datetime
    from flipchem import Flipchem

    date = datetime(2017,1,4,2)
    fc = Flipchem(date)

    glat = 74.72955
    glon = -94.90576
    alt = 130.0
    ne = 5.0e11
    te = ti = 600.
    outputs = fc.get_point(glat,glon,alt,ne,te,ti)
    lthrs,sza,dec,Op,O2p,NOp,N2p,Np,NNO,N2D,ITERS = outputs


Or, we can get the ion fractions instead by replacing the last 2 lines with::

    outputs = fc.get_point_fractions(glat,glon,alt,ne,te,ti)
    lthrs,sza,dec,Op,O2p,NOp,N2p,Np,NNO,N2D,ITERS = outputs


One can also access MSIS directly::

    from flipchem import MSIS
    from datetime import datetime

    date = datetime(2017,1,4,2)
    glat = 74.72955
    glon = -94.90576
    alt = 130.0

    msis = MSIS(date)
    outputs = msis.get_point(glat,glon,alt)
    He,O,N2,O2,Ar,Mass,H,N,AnomO,Texo,Tn = outputs

And one can also access the f107, f107a, and AP::

    from datetime import datetime
    import flipchem

    f107, f107a, ap = flipchem.read_geophys(datetime(2017,1,4,2))


And there is code available with the MSIS wrapper that provides ion-neutral and electron-neutral collision frequencies::

    from flipchem import MSIS
    from flipchem import compute_ion_neutral_collfreq
    from flipchem import compute_electron_neutral_collfreq
    from datetime import datetime

    date = datetime(2017,1,4,2)
    glat = 74.72955
    glon = -94.90576
    alt = 130.0

    msis = MSIS(date)
    outputs = msis.get_point(glat,glon,alt)
    He,O,N2,O2,Ar,Mass,H,N,AnomO,Texo,Tn = outputs
    
    # N+, O+, N2+, NO+, O2+
    ion_masses = [14.0,16.0,28.0,30.0,32.0]
    Te = Ti = 1000.0
    nu_in = list()
    neutral_densities = (H,He,N,O,N2,O2)
    for mass in ion_masses:
        nu_in.append(compute_ion_neutral_collfreq(neutral_densities, Tn, mass, Ti))
    nu_en = compute_electron_neutral_collfreq(neutral_densities, Te)

Example Notebook
----------------

`Here you can find an example notebook that shows how to get altitude profiles of ion densities <https://nbviewer.jupyter.org/github/amisr/flipchem/blob/v2020.1.1/notebooks/usage_examples.ipynb>`_. Is there an example missing that you would like to see? Feel free to suggest one!

.. [Oliphant2006] Oliphant, T. E. (2006). A guide to NumPy (Vol. 1). Trelgol Publishing USA.
.. [Picone2002] Picone, J. M., Hedin, A. E., Drob, D. P., and Aikin, A. C. (2002). NRLMSISE‐00 empirical model of the atmosphere: Statistical comparisons and scientific issues, J. Geophys. Res., 107(A12), 1468, doi:10.1029/2002JA009430. 
.. [Richards2010] Richards, P. G., Bilitza, D., and Voglozin, D. (2010), Ion density calculator (IDC): A new efficient model of ionospheric ion densities, Radio Sci., 45, RS5007, doi:10.1029/2009RS004332.
.. [Richards2011] Richards, P. G. (2011). Reexamination of ionospheric photochemistry, J. Geophys. Res., 116, A08307, doi:10.1029/2011JA016613.