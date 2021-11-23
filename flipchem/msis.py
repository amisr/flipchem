# -*- coding: utf-8 -*-
# Based on code written by Mike Nicolls in 2007?
# Modified and documented: Ashton S. Reimer 2020

from __future__ import division, absolute_import, print_function
from datetime import datetime

import numpy as np
from flipchem import read_geophys
from flipchem.ext import _c_msis as _msis


def compute_ion_neutral_collfreq(densities, Tn, mi, Ti=None):
    """This code calculates the elastic and resonant ion-neutral collision
    frequencies following Chapter 4 of [1]_.

    Parameters
    ==========
    densities : array_like
        An array of neutral densities in this order: H, He, N, O, N2, O2, Ar
        with units of number per cubic meter.
    Tn : float
        The mean neutral temperature in Kelvin
    mi : integer
        The ion mass in amu
    Ti : float, optional
        The ion temperature in Kelvin

    Returns
    =======
    nu_in : float
        The the total ion-neutral collision frequency summed over the
        collisions between the input ion and H, He, N, O, N2, O2, Ar

    References
    ==========

    .. [1] Schunk, R., & Nagy, A. (2000). Ionospheres: Physics, Plasma
           Physics, and Chemistry (Cambridge Atmospheric and Space Science
           Series). Cambridge: Cambridge University Press. 99-054707
           ISBN: 0 521 60770 1
    .. [2] Gaiser C, Fellmuth B. Polarizability of Helium, Neon, and Argon: 
           New Perspectives for Gas Metrology. Phys Rev Lett. 2018 Mar 23;
           120(12):123203. doi: 10.1103/PhysRevLett.120.123203. 
            
    """
    # set the ion temperature if None
    if Ti is None:
        Ti = Tn

    # Resonant collision frequencies (mi==##) come from
    # table 4.5 page 99 of Schunk and Nagy Ionospheres text book 2000
    # Non-resonant collisions come from equation 4.88 on page 83 and
    # table 4.1. There's some unit conversion needed, but basically
    # the equations come from:
    #    n_in = 2.21*pi*n_n*sqrt(gamma_n*e**2*m_n)/sqrt(m_i*(m_i+m_n))
    # which reduces to (after unit conversion to SI):
    #    n_in = 2.58790619679528e-3*n_n*sqrt(gamma_n*m_n)/sqrt(m_i*(m_i+m_n))
    # where gamma_n is in units of 1e-24 cm^3 and the masses are in amu
    # then, taking the values for gamma in table 4.1 and the amu for the 
    # neutrals, we can get an equation that is:
    #    n_in = const * n_n / sqrt(m_i*(m_i+m_n))
    # for each ion-neutral pair. Here's the constants for H, He, N, O, N2, O2,
    # and Ar:
    Hconst  = 2.118436361550408e-15
    Heconst = 4.796838917e-15 # From 10.1103/PhysRevLett.120.123203
    Nconst  = 1.0293931179488746e-14
    Oconst  = 9.0841307137989e-15
    N2const = 1.8168261427597804e-14
    O2const = 1.8518806819759527e-14
    Arconst = 4.291797148e-14 # From 10.1103/PhysRevLett.120.123203

    # define some other constants
    Tr = (Tn + Ti) / 2.0
    sqrtTr = np.sqrt(Tr)
    log10Tr = np.log10(Tr)

    # Now calculate the total ion-neutral collision frequency
    nu_in = 0.0

    # H, no resonant because we don't include H+ in ISR fitting
    nu_in += densities[0] * Hconst / np.sqrt(mi * (mi + 1.0))
    # He, no resonant because we don't include He+ in ISR fitting
    nu_in += densities[1] * Heconst / np.sqrt(mi * (mi + 4.0))
    # Ar, no resonant because we don't include Ar+ in ISR fitting
    nu_in += densities[6] * Arconst / np.sqrt(mi * (mi + 40.0))
    # N
    if mi == 14.0:
        nu_in += 1.0e-6 * densities[2] * 3.83e-11 * sqrtTr * (1.0 - 0.063 * log10Tr)**2.0
    else:
        nu_in += densities[2]* Nconst / np.sqrt(mi * (mi + 14.0))
    # O
    if mi == 16.0:
        nu_in += 1.0e-6 * densities[3] * 3.67e-11 * sqrtTr * (1.0 - 0.064 * log10Tr)**2.0
    else:
        nu_in += densities[3] * Oconst / np.sqrt(mi * (mi + 16.0))
    # N2
    if mi == 28.0:
        nu_in += 1.0e-6 * densities[4] * 5.14e-11 * sqrtTr * (1.0 - 0.069 * log10Tr)**2.0
    else:
        nu_in += densities[4] * N2const / np.sqrt(mi * (mi + 28.0))
    # 02
    if mi == 32.0 and Tr > 800.0:
        nu_in += 1.0e-6 * densities[5] * 2.59e-11 * sqrtTr * (1.0 - 0.073 * log10Tr)**2.0
    else:
        nu_in += densities[5] * O2const / np.sqrt(mi * (mi + 32.0))
    
    return nu_in


def compute_electron_neutral_collfreq(densities, Te):
    """
    This code calculates the elastic electron-neutral collision frequencies
    following Chapter 4 of [3]_.

    Parameters
    ==========
    densities : array_like
        An array of neutral densities in this order: H, He, N, O, N2, O2 with
        units of number per cubic meter.
        See Notes for comment about Nitrogen.
    Te : float
        The electron temperature in Kelvin

    Returns
    =======
    nu_en : float
        The the total electron-neutral collision frequency summed over the
        collisions between electrons and H, He, O, N2, O2.

    Notes
    =====
    The current code DOES NOT include an electron-Nitrogen collision frequency.

    References
    ==========

    .. [3] Schunk, R., & Nagy, A. (2000). Ionospheres: Physics, Plasma
           Physics, and Chemistry (Cambridge Atmospheric and Space Science
           Series). Cambridge: Cambridge University Press. 99-054707
           ISBN: 0 521 60770 1

    """

    # Elastic electron-neutral collision frequencies only
    # Table 4.6 page 99 of Schunk and Nagy Ionospheres text book 2000
    sqrtTe = np.sqrt(Te)
    nu_en = 0.0
    nu_en += 1e-6 * densities[0] * 4.5e-9 * (1.0 - 1.35e-4 * Te) * sqrtTe # H
    nu_en += 1e-6 * densities[1] * 4.6e-10 * sqrtTe # He
    # no N in Schunk and Nagy!
    nu_en += 1e-6 * densities[3] * 8.9e-11 * (1.0 + 5.7e-4 * Te) * sqrtTe # O
    nu_en += 1e-6 * densities[4] * 2.33e-11 * (1.0 - 1.2e-4 * Te) * Te # N2
    nu_en += 1e-6 * densities[5] * 1.82e-10 * (1.0 + 3.6e-2 * sqrtTe) * sqrtTe # O2
    
    return nu_en


class MSIS():
    """A python wrapper to the NRLMSISE-00 C library version of the code written
    by Dominik Brodowski, which is based on the original Fortran version of
    the model. See [3]_.

    Parameters
    ==========
    date : :class:`datetime.datetime`
        Date and time for which to evaluate MSIS.

    Attributes
    ==========
    f107 : float
        The F10.7 solar flux for the previous day.
    f107a : float
        The 81 day average F10.7 solar flux.
    ap : array_like
        An array of AP index values. See :func:`read_geophys`.

    Notes
    =====
    This code automatically grabs the f10.7, f10.7a, and AP, are
    required to run MSIS, from a local cache. See :func:`read_geophys`
    for more details.

    Examples
    ========
    ::

        from flipchem import MSIS
        from datetime import datetime

        date = datetime(2017,1,4,2)
        alt = 100
        glat = 60
        glon = -70

        msis = MSIS(date)
        outputs = msis.get_point(glat,glon,alt)
        H,He,N,O,N2,O2,Ar,Mass,AnomO,Texo,Tn = outputs

    References
    ==========

    .. [3] Picone, J. M., Hedin, A. E., Drob, D. P., and Aikin, A. C. (2002).
           NRLMSISE‐00 empirical model of the atmosphere: Statistical comparisons
           and scientific issues, J. Geophys. Res., 107(A12), 1468,
           doi:10.1029/2002JA009430. 

    """

    def __init__(self,date):
        self.date = date

        output = read_geophys(self.date)
        self.f107 = output[0]
        self.f107a = output[1]
        self.ap = output[2]

        self._inputs = _msis.new_nrlmsise_input()
        self._flags = _msis.new_nrlmsise_flags()
        self._aparray = _msis.new_ap_array()
        self._outputs = _msis.new_nrlmsise_output()
        self._doublearray2 = _msis.new_doubleArray(2)
        self._doublearray7 = _msis.new_doubleArray(7)
        self._doublearray9 = _msis.new_doubleArray(9)
        self._flags_sw = _msis.new_doubleArray(24)
        self._flags_swc = _msis.new_doubleArray(24)
        self._flags_switches = _msis.new_intArray(24)


    def __del__(self):
        _msis.delete_nrlmsise_input(self._inputs)
        _msis.delete_nrlmsise_flags(self._flags)
        _msis.delete_ap_array(self._aparray)
        _msis.delete_nrlmsise_output(self._outputs)
        _msis.delete_doubleArray(self._doublearray2)
        _msis.delete_doubleArray(self._doublearray7)
        _msis.delete_doubleArray(self._doublearray9)
        _msis.delete_doubleArray(self._flags_sw)
        _msis.delete_doubleArray(self._flags_swc)
        _msis.delete_intArray(self._flags_switches)


    def get_point(self,glat,glon,alt):
        """Evaluates the flipchem model for the input geodetic coordinates

        Parameters
        ==========
        glat : float
            Geodetic Latitude
        glon : float
            Geodetic Longitude
        alt : float
            Altitude above the Geodetic surface of the Earth

        Returns
        =======
        H : float
            Hydrogen number density in number per cubic meter
        He : float
            Helium number density in number per cubic meter
        N : float
            Nitrogen number density in number per cubic meter
        O : float
            Atomic Oxygen number density in number per cubic meter
        N2 : float
            Diatomic Nitrogen number density in number per cubic meter
        O2 : float
            Diatomic Oxygen number density in number per cubic meter
        Ar : float
            Argon number density in number per cubic meter
        Mass : float
            Total mass in kilograms
        AnomO : float
            Anomlous Oxygen number density in number per cubic meter
        Texo : float
            Exosphere temperature in Kelvin
        Tn : float
            Mean neutral temperature in Kelvin

        """
        # extract year and doy from given datetime
        year = self.date.year
        doy = int((self.date - datetime(year,1,1)).total_seconds()/86400) + 1
        hrut = self.date.hour + self.date.minute/60. + self.date.second/3600.0

        # run msis
        densities, temps = self._call_library(year,doy,hrut,alt,glat,glon,self.ap,
                                              self.f107a,self.f107,-1)

        He = densities[0] * 1e6
        O  = densities[1] * 1e6
        N2 = densities[2] * 1e6
        O2 = densities[3] * 1e6
        Ar = densities[4] * 1e6
        Mass = densities[5] * 1e6
        H  = densities[6] * 1e6
        N  = densities[7] * 1e6
        AnomO = densities[8] * 1e6
        Texo = temps[0]
        Tn = temps[1]
            
        return H,He,N,O,N2,O2,Ar,Mass,AnomO,Texo,Tn


    def _call_library(self,yr,doy,hrUT,alt,glat,glong,ap_array,f107a,f107,ap=-1):
        
        # input structure
        _msis.nrlmsise_input_year_set(self._inputs,int(yr))
        _msis.nrlmsise_input_doy_set(self._inputs,int(doy))
        _msis.nrlmsise_input_sec_set(self._inputs,_msis.PyFloat_AsDouble(hrUT * 3600))
        _msis.nrlmsise_input_alt_set(self._inputs,_msis.PyFloat_AsDouble(alt))
        _msis.nrlmsise_input_g_lat_set(self._inputs,_msis.PyFloat_AsDouble(glat))
        _msis.nrlmsise_input_g_long_set(self._inputs,_msis.PyFloat_AsDouble(glong))
        _msis.nrlmsise_input_lst_set(self._inputs,_msis.PyFloat_AsDouble(hrUT + glong / 15.0))
        _msis.nrlmsise_input_f107A_set(self._inputs,_msis.PyFloat_AsDouble(f107a))
        _msis.nrlmsise_input_f107_set(self._inputs,_msis.PyFloat_AsDouble(f107))
        _msis.nrlmsise_input_ap_set(self._inputs,_msis.PyFloat_AsDouble(ap))

        # ap_array
        for i in range(7):
            _msis.doubleArray___setitem__(self._doublearray7,i,float(ap_array[i]))
        _msis.ap_array_a_set(self._aparray,self._doublearray7)
        _msis.nrlmsise_input_ap_a_set(self._inputs,self._aparray)

        # set the flags in the flags structure
        sw = [0,1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
        for i in range(24):
            _msis.doubleArray___setitem__(self._flags_sw,i,0)
            _msis.doubleArray___setitem__(self._flags_swc,i,0)
            _msis.intArray___setitem__(self._flags_switches,i,sw[i])
        _msis.nrlmsise_flags_sw_set(self._flags,self._flags_sw)
        _msis.nrlmsise_flags_swc_set(self._flags,self._flags_swc)
        _msis.nrlmsise_flags_switches_set(self._flags,self._flags_switches)
        
        # call MSIS
        _msis.gtd7(self._inputs,self._flags,self._outputs)

        # get the outputs!
        densities = list()
        temp = _msis.nrlmsise_output_d_get(self._outputs)
        d = _msis.doubleArray_frompointer(temp)
        for i in range(9):
            temp = _msis.doubleArray___getitem__(d,i)
            densities.append(temp)
        densities = np.array(densities)

        temperatures = list()
        temp = _msis.nrlmsise_output_t_get(self._outputs)
        t = _msis.doubleArray_frompointer(temp)
        for i in range(2):
            temp = _msis.doubleArray___getitem__(t,i)
            temperatures.append(temp)
        temperatures = np.array(temperatures)

        return densities, temperatures
        
