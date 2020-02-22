# -*- coding: utf-8 -*-
"""
    Based on code written by Mike Nicolls in 2007?

    Modified and documented: Ashton Reimer 2019

    Compatible with python2.7 and python3

"""
from __future__ import division, absolute_import, print_function

import os
import numpy as np
from datetime import datetime


try:
    # python 3
    from urllib.request import urlopen
except:
    # python 2
    from urllib import urlopen

GEOPHYSDIR = os.path.join(os.path.abspath(os.path.dirname(__file__)),'dat')

# Reads geophysical indicies from an ASCII file
# Daily updated "geophys_params" files are provided
# at: https://amisr.com/geophys_params/
def read_geophys(date):
    """Parses a directory of files, expecting one file per year. See `Notes~
    for file structure.

    Parameters
    ==========
    date : :class:`datetime.datetime`
        Date and time for which to read geophysical parameters

    Returns
    =======

    f107 : float
        F10.7 for the previous day
    f107a : float
        An 81 day average of the F10.7 values
    ap : array_like
        An array of AP values where each index in the array corresponds to::
        0 : daily AP
        1 : 3 hr AP index for current time
        2 : 3 hr AP index for 3 hrs before current time
        3 : 3 hr AP index for 6 hrs before current time
        4 : 3 hr AP index for 9 hrs before current time
        5 : Average of eight 3 hr AP indicies from 12 to 33 hrs 
                prior to current time
        6 : Average of eight 3 hr AP indicies from 36 to 57 hrs 
                prior to current time 
    
    Notes
    =====
    Each file in the directory of files is expected to have the following
    structure:

    1901012529 81013271720 3 0 7 97  4  5 12  6  7  2  0  3  50.21---069.50
    1901022529 9 0 0 0 0 3 3 0 0  7  0  0  0  0  2  2  0  0  00.00---072.70
    19010325291010 0 0 3 3 0 0 3 20  4  0  0  2  2  0  0  2  10.00---070.20
    etc.
    


    Usage
    =====

    from datetime import datetime
    import flipchem

    f107, f107a, ap = flipchem.read_geophys(datetime(2017,1,4,2))


    """
    # extract year and doy from given datetime
    year = date.year
    doy = int((date - datetime(year,1,1)).total_seconds()/86400) + 1
    curtime = date.hour + date.minute/60. + date.second/3600.0

    # If no geophysical parameters for current year, try previous year
    if os.path.exists(os.path.join(GEOPHYSDIR,str(year))) == False:
        year -= 1
    if os.path.exists(os.path.join(GEOPHYSDIR,str(year))) == False:
        raise IOError('Geophys param directory %s/%s does not exist.' % (GEOPHYSDIR, str(year)))

    # open and read the file for the previous, current year, and next year
    year = str(year)
    with open(os.path.join(GEOPHYSDIR,year)) as f:
        lines = f.readlines()
    if os.path.exists(os.path.join(GEOPHYSDIR,str(int(year)-1))) == True:
        with open(os.path.join(GEOPHYSDIR,str(int(year)-1))) as f:
            lines_py = f.readlines()
    else:
        lines_py = list()

    if os.path.exists(os.path.join(GEOPHYSDIR,str(int(year) + 1))) == True:
        f = open(os.path.join(GEOPHYSDIR,str(int(year) + 1)))
        lines_ny = f.readlines()
    else:
        lines_ny = list()
    
    # if we don't have data up to the current doy, we'll settle for the latest data
    if len(lines) < doy:
        doy = len(lines)
    
    # F107d - previous day
    if doy == 1:
        try:
            F107D = float(lines_py[-1][65:70])
        except: # if can't get prev day, just get current day
            F107D = float(lines[doy - 1][65:70])
    else:
        F107D = float(lines[doy - 2][65:70])
    
    # AP
    AP = np.zeros((7))
    AP[0] = float(lines[doy - 1][55:58]) # daily, for current day 
    TINDEX = int(curtime / 3.0)
    tlines = lines
    tind = doy - 1
    for i in range(4):
        AP[i + 1] = float(tlines[tind][31+TINDEX*3:34+TINDEX*3]) # for current time 
        if TINDEX > 0:
            TINDEX -= TINDEX
        else:
            if tind == 0:
                tlines = lines_py
                tind = len(tlines) - 1
            else:
                tind -= 1
            TINDEX = 7
    for i in range(8):
        AP[5] += float(tlines[tind][31+TINDEX*3:34+TINDEX*3]) # for current time 
        if TINDEX > 0:
            TINDEX -= TINDEX - 1
        else:
            if tind == 0:
                tlines = lines_py
                tind = len(tlines) - 1
            else:
                tind = tind - 1
            TINDEX = 7
    AP[5] /= 8.0
    for i in range(8):
        AP[6] += float(tlines[tind][31+TINDEX*3:34+TINDEX*3]) # for current time 
        if TINDEX > 0:
            TINDEX -= 1
        else:
            if tind == 0:
                tlines = lines_py
                tind = len(tlines) - 1
            else:
                tind -= 1
            TINDEX = 7
    AP[6] /= 8.0    
        
    F107A = 0
    imin = doy - 1 - 40
    imax = doy + 40
    
    if imax > len(lines) and len(lines_ny) == 0:
        imin = imin - (imax - len(lines)) + len(lines_ny)
        imax = imax - (imax - len(lines)) + len(lines_ny)
    
    # F107A
    lines2 = []
    for aa in range(imin,imax):
        try:
            if aa < 0:
                lines2.append(lines_py[len(lines_py) + aa])
            elif aa >= len(lines):        
                lines2.append(lines_ny[aa-len(lines)])
            else:
                lines2.append(lines[aa])
        except:
            pass
    
    F107A = 0.0
    for aa in range(len(lines2)):
        try:
            F107A = F107A + float(lines2[aa][65:70])
        except:
            try:
                F107A = F107A + float(lines2[aa + 1][65:70])
            except:
                F107A = F107A + float(lines2[aa - 1][65:70])            
    F107A = F107A / len(lines2)

    return F107D, F107A, AP


# helper function for downloading data for a specific year
def _download_year(year,base_url=None):

    url = os.path.join(base_url,str(year))
    try:
        resp = urlopen(url)
        content = resp.read()
    except Exception as e:
        text = 'Problem encountered trying to download geophysical'
        text += 'parameters from %s\nException details:\n%s' % (url,str(e))
        raise Exception(text)

    return content.decode('utf-8')


# update geophysical parameter files
def update_geophys(year=None,base_url='https://amisr.com/geophys_params/'):
    """This function downloads geophysical parameter files from the default
    url: %s

    Parameters
    ==========
    year : int, optional
        The year for which to download a geophysical parameter file. If
        year is None, all files are downloaded.

    base_url : str, optional
        The URL at which the geophysical parameter files are hosted.

    Notes
    =====

    The files are updated on a daily basis, so it is important to update
    the local copy of these files on a regular basis.

    Usage
    =====

    import flipchem
    flipchem.update_geophys(2020)


    """ % (base_url)

    # current year
    cur_year = datetime.now().year

    if year is None:
        years = [cur_year]
    elif year == 'all':
        years = range(1947,cur_year+1)
    else:
        years = [year]

    # iterate through the years we need to download
    print("Updating geophysical parameters files:")
    for yr in years:
        print("   ...doing %s" % str(yr))
        # download the data for the year
        content = _download_year(yr,base_url=base_url)
        # write it to file
        save_path = os.path.join(GEOPHYSDIR,str(yr))
        with open(save_path,'w') as f:
            f.write(content)
    print("Done!")




