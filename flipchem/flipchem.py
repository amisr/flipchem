# -*- coding: utf-8 -*-
# Based on code written by Mike Nicolls in 2007?
# Modified and documented: Ashton S. Reimer 2020

from __future__ import division, absolute_import, print_function

from datetime import datetime

import numpy as np
from flipchem import read_geophys, MSIS
from flipchem.ext import chemion as _chemion
from flipchem.ext import getltsza as _getltsza

class Flipchem():
    """A python wrapper to the flipchem ionospheric photochemistry model
    developed by Phil Richards [1]_. Specifically, this code wraps the version
    of flipchem that was used for [2]_. All ion densities, except for O+ are
    calculated from chemical equilibrium.

    Parameters
    ==========
    date : :class:`datetime.datetime`
        Date and time for which to evaluate the flipchem model.
    altop : float, optional
        Altitude above which ion fractions are set to 100% O+. Only used by
        the `get_point_fractions` function.

    Attributes
    ==========
    f107 : float
        The F10.7 solar flux for the previous day.
    f107a : float
        The 81 day average F10.7 solar flux.
    ap : array_like
        An array of AP index values. See `read_geophys`.

    Notes
    =====
    flipchem requires an input neutral atmosphere and the geophysical
    indicies: f107 and f107a. NRLMSIS-00 is used to provide a neutral
    atmosphere (see the :class:`flipchem.MSIS`). See `read_geophys` for more
    details about geophysical indicies.

    Examples
    ========
    from datetime import datetime
    from flipchem import Flipchem

    date = datetime(2017,1,4,2)
    fc = Flipchem(date)

    glat = 60
    glon = -70
    alt = 300
    ne = 5.0e11
    te = ti = 500.
    outputs = fc.get_point_fractions(glat,glon,alt,ne,te,ti)
    LTHRS,SZAD,DEC,OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D,ITERS = outputs

    References
    ==========

    .. [1] Richards, P. G. (2011), Reexamination of ionospheric photochemistry,
           J. Geophys. Res., 116, A08307, doi:10.1029/2011JA016613.

    .. [2] Richards, P. G., Bilitza, D., and Voglozin, D. (2010), Ion density
           calculator (IDC): A new efficient model of ionospheric ion densities,
           Radio Sci., 45, RS5007, doi:10.1029/2009RS004332. 


    """

    def __init__(self,date,altop=1000.0):
        
        self.date = date
        self.altop = altop

        output = read_geophys(self.date)
        self.f107 = output[0]
        self.f107a = output[1]
        self.ap = output[2]


    def get_point(self,glat,glon,alt,ne,te,ti,user_no=-1.0,user_oplus=-1.0,
                  msis_outputs=None,fractions=False):
        """Evaluates the flipchem model for the input geodetic coordinates

        Parameters
        ==========
        glat : float
            Geodetic Latitude
        glon : float
            Geodetic Longitude
        alt : float
            Altitude above the Geodetic surface of the Earth
        ne : float
            Electron density in number per cubic meter
        te : float
            Electron temperature in Kelvin
        ti : float
            Ion temperature in Kelvin
        user_no : float, optional
            If positive, replaces the NO density in number per cubic meter determined
            from chemical equilibrium with user input density.
        user_oplus : float, optional
            If positive, used as the O+ density in number per cubic meter. If negative
            the O+ number density is determined from photochemistry.
        msis_outputs : tuple, optional
            A tuple of the outputs of `msis.get_point`
        fractions : bool, optional
            A flag for whether the returned ion densities should be instead returned
            as fractions.

        Returns
        =======
        lthrs : float
            Local time in decimal hours
        sza : float
            The solar zenith angle in degrees
        dec : float
            The solar declination angle in degrees
        oxplus : float
            The O+ density in number per cubic meter
        o2plus : float
            The O2+ density in number per cubic meter
        noplus : float
            The NO+ density in number per cubic meter
        n2plus : float
            The N2+ density in number per cubic meter
        nplus : float
            The N+ density in number per cubic meter
        nno : float
            The NO density in number per cubic meter
        n2d : float
            The N2(D) density in number per cubic meter
        iters : integer
            The number of iterations required before chemical equilibrium was achieved

        """

        if alt > self.altop:
            OXPLUS = ne
            O2PLUS = 0.0
            NOPLUS = 0.0
            N2PLUS = 0.0
            NPLUS = 0.0
            lthrs,szad,dec = self._call_getltsza(date,lat,lon)

        else:

            # call msis to get neutral densities if needed
            if msis_outputs is None:
                msis = MSIS(self.date)
                msis_outputs = msis.get_point(glat,glon,alt)

            H,He,N,O,N2,O2,Ar,Mass,AnomO,Texo,Tn = msis_outputs

            # now call flipchem
            flip_outputs = self._call_flip(self.date,glat,glon,alt,self.ap,self.f107a,
                                           self.f107,te,ti,Tn,O,O2,N2,He,H,N,ne,
                                           user_no=user_no,user_oplus=user_oplus)

            lthrs,szad,dec,OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D,ITERS = flip_outputs

            # scale into SI units (m^-3)
            OXPLUS = OXPLUS * 1.0e6
            O2PLUS = O2PLUS * 1.0e6
            NOPLUS = NOPLUS * 1.0e6
            N2PLUS = N2PLUS * 1.0e6
            NPLUS = NPLUS * 1.0e6

        # return fractions if requested
        if fractions:
            OXPLUS = OXPLUS / ne
            O2PLUS = O2PLUS / ne
            NOPLUS = NOPLUS / ne
            N2PLUS = N2PLUS / ne
            NPLUS = NPLUS / ne

        return lthrs,szad,dec,OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D,ITERS

    def _call_flip(self,date,lat,lon,alt,ap,f107a,f107,te,ti,tn,OXN,O2N,N2N,HEN,HN,N4S,NE,user_no=-1.0,user_oplus=-1.0):
        """Hidden method used to call chemion the fortran subroutine.

        Not intented to be used by users.
        """
        # get local time and solar zenith angle/declination
        lthrs,szad,decd = self._call_getltsza(date,lat,lon)

        #chemion flags, see FLIP-CHEM.f for details
        #file output off
        jprint = 0

        # needs neutral densities in per cubic centimeter
        OXN *= 1.0e-6
        O2N *= 1.0e-6
        N2N *= 1.0e-6
        HEN *= 1.0e-6
        HN *= 1.0e-6
        N4S *= 1.0e-6

        # As per paragraph [13] of Richards 2011, input N4S density needs to
        # be halved if it comes from NRLMSIS-00.
        # all density units need to be in cm^-3
        outputs = _chemion(jprint,alt,f107,f107a,te,ti,tn,OXN,O2N,N2N,HEN,HN,user_no,0.5*N4S,NE*1.0e-6,user_oplus,szad)
        OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D,iters = outputs
            
        return lthrs,szad,decd,OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D,iters

    def _call_getltsza(self,date,lat,lon):
        """ Hidden method used to call the getltsza fortran subroutine.

        Not intented to be used by users.
        """
        # extract year and doy from given datetime
        year = date.year
        doy = int((date - datetime(year,1,1)).total_seconds()/86400) + 1
        curtime = date.hour + date.minute/60. + date.second/3600.0
        utsecs = curtime * 3600.0
        yyyyddd = year*1000 + doy
        
        lthrs,sza,dec = _getltsza(yyyyddd,utsecs,lat,lon)

        # convert to degrees
        szad = sza*180.0/np.pi
        decd = dec*180.0/np.pi

        return lthrs,szad,decd