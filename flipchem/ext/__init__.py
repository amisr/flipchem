# -*- coding: utf-8 -*-
import os

# if on readthedocs, don't even import these extensions
if os.environ.get('READTHEDOCS', None) == 'True':
    chemion = None
    getltsza = None
    gtd7 = None
    _f_flipchem = None
    _c_msis = None
else:
    from ._f_flipchem import chemion, getltsza
    from ._c_msis import gtd7