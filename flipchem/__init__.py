# -*- coding: utf-8 -*-
from __future__ import absolute_import

from . import ext
from .geophys import read_geophys, update_geophys
from .msis import MSIS
from .msis import compute_ion_neutral_collfreq, compute_electron_neutral_collfreq
from .msis import compute_electron_neutral_collfreq
from .flipchem import Flipchem

__version__ = '2020.1.0'
__doc__ = """
This package contains the a python wrapper of the flipchem ionospheric
photochemistry model developed by Phil Richards:

Richards, P. G. (2011), Reexamination of ionospheric photochemistry,
J. Geophys. Res., 116, A08307, doi:10.1029/2011JA016613.

The package is intended to be used for fitting Incoherent Scatter Radar
data, which involves fitting range gates at fixed time intervals. Also,
running flipchem requires both a neutral atmosphere model and geophysical
indicies. This is why the indicies are read when the Flipchem class is
initialized and the 'get_point*' functions do not accept date. The
intended usage is:
    1) initialize the Flipchem class for one time interval
    2) call one of the get_point* functions for each of the range gates
       in the one time interval
    3) initialize a new Flipchem class for the next time interval

This usage case facilitates usage of the multiprocessing library, where
each time interval is treated as a pool of asynchronous jobs.

"""