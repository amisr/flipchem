# -*- coding: utf-8 -*-
"""
    Based on code written by Mike Nicolls in 2007?

    Modified and documented: Ashton Reimer 2019

"""
from __future__ import division, absolute_import, print_function
from datetime import datetime

import numpy as np
from flipchem import read_geophys, MSIS
from flipchem.ext import chemion as _chemion
from flipchem.ext import getltsza as _getltsza

class Flipchem():
    """
    A python wrapper to the flipchem ionospheric photochemistry model
    developed by Phil Richards:

    Richards, P. G. (2011), Reexamination of ionospheric photochemistry,
    J. Geophys. Res., 116, A08307, doi:10.1029/2011JA016613. 

    flipchem requires an input neutral atmosphere and the geophysical
    indicies: f107, f107a, and AP. NRLMSIS-00 is used to provide a 
    neutral atmosphere (see the flipchem.MSIS class). See read_geophys
    for more details about geophysical indicies.


    Example Usage:

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
        LTHRS,SZAD,DEC,OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D,success = outputs

    """

    def __init__(self,date,altop=1000.0):
        
        self.date = date
        self.altop = altop

        output = read_geophys(self.date)
        self.f107 = output[0]
        self.f107a = output[1]
        self.ap = output[2]


    def get_point(self,glat,glon,alt,ne,te,ti,user_no=-1.0,user_oplus=-1.0,msis_outputs=None):

        if msis_outputs is None:
            msis = MSIS(self.date)
            msis_outputs = msis.get_point(glat,glon,alt)

        HEdens,Odens,N2dens,O2dens,ARdens,MassDens,Hdens,Ndens,AnomOdens,Texo,Tn = msis_outputs

        # now call flipchem
        flip_outputs = self.call_flip(self.date,glat,glon,alt,self.ap,self.f107a,
                                      self.f107,te,ti,Tn,Odens,O2dens,N2dens,HEdens,
                                      Ndens,ne,user_no=user_no,user_oplus=user_oplus)

        LTHRS,SZAD,DEC,OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D = flip_outputs[:-1]
        success = bool(flip_outputs[-1])

        # scale into SI units (m^-3)
        OXPLUS = OXPLUS * 1.0e6
        O2PLUS = O2PLUS * 1.0e6
        NOPLUS = NOPLUS * 1.0e6
        N2PLUS = N2PLUS * 1.0e6
        NPLUS = NPLUS * 1.0e6

        return LTHRS,SZAD,DEC,OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D,success

    def get_point_fractions(self,glat,glon,alt,ne,te,ti,user_no=-1.0,user_oplus=-1.0,msis_outputs=None,minval=0.001,maxval=1.0,altop=300.0):
        flip_outputs = self.get_point(glat,glon,alt,ne,te,ti,user_no=-1.0,user_oplus=-1.0,msis_outputs=msis_outputs)
        LTHRS,SZAD,DEC,OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D,success = flip_outputs

        # compute fractions
        if alt > altop:
            OXPLUS = 1.0
            O2PLUS = 0.0
            NOPLUS = 0.0
            N2PLUS = 0.0
            NPLUS = 0.0
        else:
            OXPLUS = OXPLUS / ne
            O2PLUS = O2PLUS / ne
            NOPLUS = NOPLUS / ne
            N2PLUS = N2PLUS / ne
            NPLUS = NPLUS / ne
        
        # check ranges
        if OXPLUS < minval:
            OXPLUS = 0.0
        elif OXPLUS > maxval:
            OXPLUS = 1.0
        if O2PLUS < minval:
            O2PLUS = 0.0
        elif O2PLUS > maxval:
            O2PLUS = 1.0
        if NOPLUS < minval:
            NOPLUS = 0.0
        elif NOPLUS > maxval:
            NOPLUS = 1.0
        if N2PLUS < minval:
            N2PLUS = 0.0
        elif N2PLUS > maxval:
            N2PLUS = 1.0
        if NPLUS < minval:
            NPLUS = 0.0
        elif NPLUS > maxval:
            NPLUS = 1.0

        return LTHRS,SZAD,DEC,OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D,success

    def call_flip(self,date,lat,lon,alt,ap,f107a,f107,te,ti,tn,OXN,O2N,N2N,HEN,N4S,NE,user_no=-1.0,user_oplus=-1.0):

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

        #chemion flags, see FLIP-CHEM.f for details
        #file output off
        jprint = 0

        # As per paragraph [13] of Richards 2011, input N4S density needs to
        # be halved if it comes from NRLMSIS-00.
        # all density units need to be in cm^-3
        outputs = _chemion(jprint,alt,ap,f107,f107a,te,ti,tn,OXN,O2N,N2N,HEN,
                           user_no,0.5*N4S,NE*1.0e-6,user_oplus,szad)

        OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D,success = outputs
            
        return lthrs,szad,decd,OXPLUS,O2PLUS,NOPLUS,N2PLUS,NPLUS,NNO,N2D,success
