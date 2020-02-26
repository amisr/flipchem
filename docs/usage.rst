Some Example Usage
==================


geophys params
--------------

First, it is very important to update the geophysical parameters files. There is an archive of these files available at `<https://amisr.com/geophys_params/>`_. They are organized by year, so once updated, you only need to update the current year. To update all years::

    import flipchem
    flipchem.update_geophys(year='all')


and to update a specific year, say 1991::

    import flipchem
    flipchem.update_geophys(year=1991)

And one can read the f107, f107a, and AP::

    from datetime import datetime
    import flipchem

    f107, f107a, ap = flipchem.read_geophys(datetime(2017,1,4,2))


flipchem ion densities
----------------------

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


Or, we can get the ion fractions instead by replacing setting the fractions keyword to True::

    outputs = fc.get_point(glat,glon,alt,ne,te,ti,fractions=True)
    lthrs,sza,dec,Op,O2p,NOp,N2p,Np,NNO,N2D,ITERS = outputs


calling nrlmsise-00
-------------------

One can also access MSIS directly::

    from flipchem import MSIS
    from datetime import datetime

    date = datetime(2017,1,4,2)
    glat = 74.72955
    glon = -94.90576
    alt = 130.0

    msis = MSIS(date)
    outputs = msis.get_point(glat,glon,alt)
    H,He,N,O,N2,O2,Ar,Mass,AnomO,Texo,Tn = outputs


ion-neutral and electron-neutral collisions
-------------------------------------------

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
    H,He,N,O,N2,O2,Ar,Mass,AnomO,Texo,Tn = outputs
    
    # N+, O+, N2+, NO+, O2+
    ion_masses = [14.0,16.0,28.0,30.0,32.0]
    Te = Ti = 1000.0
    nu_in = list()
    neutral_densities = (H,He,N,O,N2,O2)
    for mass in ion_masses:
        nu_in.append(compute_ion_neutral_collfreq(neutral_densities, Tn, mass, Ti))
    nu_en = compute_electron_neutral_collfreq(neutral_densities, Te)

Example Jupyter Notebook
-------------------------

Do you prefer working in jupyter notebooks? `Here you can find an example notebook that shows how to get altitude profiles of ion densities <https://nbviewer.jupyter.org/github/amisr/flipchem/blob/v2020.2.0/notebooks/usage_examples.ipynb>`_.

Is there an example missing that you would like to see? Feel free to suggest one!