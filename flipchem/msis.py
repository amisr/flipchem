
"""
    Based on code written by Mike Nicolls in 2007?

    Modified and documented: Ashton Reimer 2019

"""
from __future__ import division, absolute_import, print_function
from datetime import datetime

import numpy as np
from flipchem import read_geophys
import flipchem.ext._c_msis as _msis


def compute_ion_neutral_collfreq(d, Tn, mi, Ti):
    """
    This code calculates the elastic and resonant ion-neutral collision
    frequencies and the elastic electron-neutral collision frequencies
    following Chapter 4 of "Ionospheres" by Schunk and Nagy:

    Schunk, R., & Nagy, A. (2000). Ionospheres: Physics, Plasma
    Physics, and Chemistry (Cambridge Atmospheric and Space Science
    Series). Cambridge: Cambridge University Press.
    99-054707 ISBN: 0 521 60770 1"""
    if Ti is None:
        Ti = Tn

    # Resonant collision frequencies (mi==##) come from
    # table 4.5 page 99 of Schunk and Nagy Ionospheres text book 2000
    # Non-resonant collisions come from equation 4.88 on page 83 and
    # table 4.1. There's some unit conversion needed, but basically
    # the equations come from:
    #    n_in = 2.21*pi*n_n*sqrt(gamma_n*e**2*m_n)/sqrt(m_i*(m_i+m_n))
    # which reduces to (after unit conversion to SI):
    #    n_in = 2.5880819319480624e-3*n_n*sqrt(gamma_n*m_n)/sqrt(m_i*(m_i+m_n))
    # where gamma_n is in units of 1e-24 cm^3 and the masses are in amu
    # then, taking the values for gamma in table 4.1 and the amu for the 
    # neutrals, we can get an equation that is:
    #    n_in = const * n_n / sqrt(m_i*(m_i+m_n))
    # for each ion. Here's the constants for H+, He+, N+, O+, N2+, and O2+:
    Hconst  = 2.118436361550408e-15
    Heconst = 2.372016271579909e-15 # assuming 100% He-4
    Nconst  = 1.0293931179488746e-14
    Oconst  = 9.0841307137989e-15
    N2const = 1.8168261427597804e-14
    O2const = 1.8518806819759527e-14

    # define some other constants
    Tr = (Tn + Ti) / 2.0
    sqrtTr = np.sqrt(Tr)
    log10Tr = np.log10(Tr)

    # Now calculate the total ion-neutral collision frequency
    nu_in = 0.0

    # H, no resonant because we don't include H+ in ISR fitting
    nu_in += 1.0e6 * d[0] * Hconst / np.sqrt(mi * (mi + 16.0))
    # He, no resonant because we don't include He+ in ISR fitting
    nu_in += 1.0e6 * d[1] * Heconst / np.sqrt(mi * (mi + 16.0))
    # N
    if mi == 14.0:
        nu_in += d[2] * 3.83e-11 * sqrtTr * (1.0 - 0.063 * log10Tr)**2.0
    else:
        nu_in += 1.0e6 * d[2]* Nconst / np.sqrt(mi * (mi + 16.0))
    # O
    if mi == 16.0:
        nu_in += d[3] * 3.67e-11 * sqrtTr * (1.0 - 0.064 * log10Tr)**2.0
    else:
        nu_in += 1.0e6 * d[3] * Oconst / np.sqrt(mi * (mi + 16.0))
    # N2
    if mi == 28.0:
        nu_in += d[4] * 5.14e-11 * sqrtTr * (1.0 - 0.069 * log10Tr)**2.0
    else:
        nu_in += 1.0e6 * d[4] * N2const / np.sqrt(mi * (mi + 28.0))
    # 02
    if mi == 32.0 and Tr > 800.0:
        nu_in += d[5] * 2.59e-11 * sqrtTr * (1.0 - 0.073 * log10Tr)**2.0
    else:
        nu_in += 1.0e6 * d[5] * O2const / np.sqrt(mi * (mi + 32.0))
    
    return nu_in


def compute_electron_neutral_collfreq(d, Te):
    """
    This code calculates the elastic electron-neutral collision frequencies
    following Chapter 4 of "Ionospheres" by Schunk and Nagy:

    Schunk, R., & Nagy, A. (2000). Ionospheres: Physics, Plasma
    Physics, and Chemistry (Cambridge Atmospheric and Space Science
    Series). Cambridge: Cambridge University Press.
    99-054707 ISBN: 0 521 60770 1"""

    # Elastic electron-neutral collision frequencies only
    # Table 4.6 page 99 of Schunk and Nagy Ionospheres text book 2000
    sqrtTe = np.sqrt(Te)
    nu_en = 0.0
    nu_en += d[0] * 4.5e-9 * (1.0 - 1.35e-4 * Te) * sqrtTe # H
    nu_en += d[1] * 4.6e-10 * sqrtTe # He
    # no N in Schunk and Nagy!
    nu_en += d[3] * 8.9e-11 * (1.0 + 5.7e-4 * Te) * sqrtTe # O
    nu_en += d[4] * 2.33e-11 * (1.0 - 1.2e-4 * Te) * Te # N2
    nu_en += d[5] * 1.82e-10 * (1.0 + 3.6e-2 * sqrtTe) * sqrtTe # O2
    
    return nu_en


class MSIS():
    """
    A python wrapper to the NRLMSIS-00 C library version of the code written
    by Dominik Brodowski, which is based on the original Fortran version of
    the model. See:

    Picone, J. M., Hedin, A. E., Drob, D. P., and Aikin, A. C. (2002).
    NRLMSISE‐00 empirical model of the atmosphere: Statistical comparisons
    and scientific issues, J. Geophys. Res., 107(A12), 1468,
    doi:10.1029/2002JA009430. 

    This code automatically grabs the f10.7, f10.7a, and AP, are
    required to run MSIS, from a local cache. See read_geophys
    for more details.

    Example Usage:

        from flipchem import MSIS
        from datetime import datetime

        date = datetime(2017,1,4,2)
        alt = 100
        glat = 60
        glon = -70

        msis = MSIS(date)
        outputs = msis.get_point(glat,glon,alt)
        He,O,N2,O2,Ar,Mass,H,N,AnomO,Texo,Tn = outputs

    """

    def __init__(self,date):
        self.date = date

        output = read_geophys(self.date)
        self.f107 = output[0]
        self.f107a = output[1]
        self.ap = output[2]

        self.inputs = _msis.new_nrlmsise_input()
        self.flags = _msis.new_nrlmsise_flags()
        self.aparray = _msis.new_ap_array()
        self.outputs = _msis.new_nrlmsise_output()
        self.doublearray2 = _msis.new_doubleArray(2)
        self.doublearray7 = _msis.new_doubleArray(7)
        self.doublearray9 = _msis.new_doubleArray(9)
        self.flags_sw = _msis.new_doubleArray(24)
        self.flags_swc = _msis.new_doubleArray(24)
        self.flags_switches = _msis.new_intArray(24)


    def __del__(self):
        _msis.delete_nrlmsise_input(self.inputs)
        _msis.delete_nrlmsise_flags(self.flags)
        _msis.delete_ap_array(self.aparray)
        _msis.delete_nrlmsise_output(self.outputs)
        _msis.delete_doubleArray(self.doublearray2)
        _msis.delete_doubleArray(self.doublearray7)
        _msis.delete_doubleArray(self.doublearray9)
        _msis.delete_doubleArray(self.flags_sw)
        _msis.delete_doubleArray(self.flags_swc)
        _msis.delete_intArray(self.flags_switches)


    def get_point(self,glat,glon,alt):

        # extract year and doy from given datetime
        year = self.date.year
        doy = int((self.date - datetime(year,1,1)).total_seconds()/86400) + 1
        hrut = self.date.hour + self.date.minute/60. + self.date.second/3600.0

        # run msis
        densities, temps = self._call_library(year,doy,hrut,alt,glat,glon,self.ap,
                                              self.f107a,self.f107,-1)

        He = densities[0]
        O  = densities[1]
        N2 = densities[2]
        O2 = densities[3]
        Ar = densities[4]
        Mass = densities[5]
        H  = densities[6]
        N  = densities[7]
        AnomO = densities[8]
        Texo = temps[0]
        Tn = temps[1]
            
        return He,O,N2,O2,Ar,Mass,H,N,AnomO,Texo,Tn


    def _call_library(self,yr,doy,hrUT,alt,glat,glong,ap_array,f107a,f107,ap=-1):
        
        # input structure
        _msis.nrlmsise_input_year_set(self.inputs,int(yr))
        _msis.nrlmsise_input_doy_set(self.inputs,int(doy))
        _msis.nrlmsise_input_sec_set(self.inputs,_msis.PyFloat_AsDouble(hrUT * 3600))
        _msis.nrlmsise_input_alt_set(self.inputs,_msis.PyFloat_AsDouble(alt))
        _msis.nrlmsise_input_g_lat_set(self.inputs,_msis.PyFloat_AsDouble(glat))
        _msis.nrlmsise_input_g_long_set(self.inputs,_msis.PyFloat_AsDouble(glong))
        _msis.nrlmsise_input_lst_set(self.inputs,_msis.PyFloat_AsDouble(hrUT + glong / 15.0))
        _msis.nrlmsise_input_f107A_set(self.inputs,_msis.PyFloat_AsDouble(f107a))
        _msis.nrlmsise_input_f107_set(self.inputs,_msis.PyFloat_AsDouble(f107))
        _msis.nrlmsise_input_ap_set(self.inputs,_msis.PyFloat_AsDouble(ap))

        # ap_array
        for i in range(7):
            _msis.doubleArray___setitem__(self.doublearray7,i,float(ap_array[i]))
        _msis.ap_array_a_set(self.aparray,self.doublearray7)
        _msis.nrlmsise_input_ap_a_set(self.inputs,self.aparray)

        # set the flags in the flags structure
        sw = [0,1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
        for i in range(24):
            _msis.doubleArray___setitem__(self.flags_sw,i,0)
            _msis.doubleArray___setitem__(self.flags_swc,i,0)
            _msis.intArray___setitem__(self.flags_switches,i,sw[i])
        _msis.nrlmsise_flags_sw_set(self.flags,self.flags_sw)
        _msis.nrlmsise_flags_swc_set(self.flags,self.flags_swc)
        _msis.nrlmsise_flags_switches_set(self.flags,self.flags_switches)
        
        # call MSIS
        _msis.gtd7(self.inputs,self.flags,self.outputs)

        # get the outputs!
        densities = list()
        temp = _msis.nrlmsise_output_d_get(self.outputs)
        d = _msis.doubleArray_frompointer(temp)
        for i in range(9):
            temp = _msis.doubleArray___getitem__(d,i)
            densities.append(temp)
        densities = np.array(densities)

        temperatures = list()
        temp = _msis.nrlmsise_output_t_get(self.outputs)
        t = _msis.doubleArray_frompointer(temp)
        for i in range(2):
            temp = _msis.doubleArray___getitem__(t,i)
            temperatures.append(temp)
        temperatures = np.array(temperatures)

        return densities, temperatures
        
