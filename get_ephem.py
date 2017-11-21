#!/usr/bin/env python
"""
  Methods for querying the JPL Horizons database. 

  Instructions for keywords and options available here:
    ftp://ssd.jpl.nasa.gov/pub/ssd/horizons_batch_example.long

  v0: M. Adamkovics
  v1: K. de Kleer
  v2: 2017-06-15 E. Molter
            added naif_lookup
            adapted to fit into whats_up.py
  v3: C. Moeckel 
            Updated the docstrings 
"""

from urllib.request import urlopen
import numpy as np
from time import strftime, gmtime, time
from datetime import datetime,timedelta
import sys, os

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def naif_lookup(target):
    """Find the NAIF target number for the planet 

    For a given planet, (enter planet as a string), it will determine 
    the corresponding NAIF idenfier. Information found at 
    naif_id_table.txt
    
    Parameters
    -------
    target : [-] str
        [-] Name of target 

    Returns
    ----------
    target : [-] int
        [-] NAIF identifier for target



    Example
    -------
    >>> get_ephem.naif_lookup('Jupiter')

    References
    ------------
    https://ssd.jpl.nasa.gov/horizons.cgi

    Notes
    -------
    10/25/2017, CM, Initial Commit
    """    


    target = target.upper().strip(', \n')
    with open('naif_id_table.txt','r') as f:
        for line in f:
            l = line.split(',')
            if is_number(target):
                #print(l[0].strip(', \'\n') == target)
                if l[0].strip(', \'\n') == target:
                    code = target
            else:
                if l[1].strip(', \'\n') == target:
                    code = l[0].strip(', \n')
        try:
            code
        except:
            print('WARN: NAIF code not in lookup table. If code fails, ensure target can be queried in Horizons.')
            code = target
                            
    if len(code) == 7: #minor body
        if code[0] == '2': #asteroid
            return code[1:]+';'
        elif code[0] =='1': #comet
            sys.exit('ERROR: Comets cannot be looked up by NAIF ID; Horizons generates multiple codes for every comet. Try your search in the Horizons Web tool, select the target body you want, and then copy the exact string into this code and it *may* work.')
        return code
    return code

def get_ephemerides(code, tstart, tend, nstep, obs_code = '500') :
    """ Ephemeris for a celestial minor body for the given time interval

    
    Obtain ephemeris data from JPL Horizons webservers on the target 
    body. 
    
    Parameters
    -------
    code : [-] str
        [-] NAIF identifier for the planet 
    tstart : [-] str
        ['yyyy-mm-dd hh:mm'] Starting time of observations 
    tend : [-] str
        ['yyyy-mm-dd hh:mm'] Ending time of observations   
    tstep : [-] str
        [-] Number of ephermis

    Returns
    ----------
    data : [-] str
        [-] The following information
    observatory_coords : [3x1] float
        [latxlongxkm] NAIF observatory coordinates 
    radii : [3x1] flaot
        [km] Radius of target body 


    0  UTdate UTtime
    1  Empty 
    2  Empty 
    3  RA (J2000) (hh:mm:ss.ff)
    4  DEC (J2000) (hh:mm:ss.ff)
    5  d(RA)/dt*cosD (arcsec^2/h) (cosine of the declination)
    6  d(DEC)/dt 
    7  Azimuth (deg) , topocentric 
    8  Elevation (deg)
    9  Airmass (-)
    10 Extinction 
    11 APmag (magnitude)
    12 Surface Brightness (magnitude per arcsec^2)
    13 Satellite angular separ/vis. (arcsec) The angle between the
        center of target object and the center of the primary body it 
        revolves around,
    14 Empty 
    15 Target angular diameter (arcsec)
    16 Observer sub-longitude, (deg) planetodetic
    17 Observer sub-latitude, (deg)
    18 Sun sub-longitude (deg), planetodetic
    19 Sun sub-latitude (deg), planetodetic
    20 Sub-Sun position angle (deg)
    21 Sub-Sun distance (arcsec)
    22 North Pole position angle  (deg)  
    23 North pole distance(arcsec) Distance from the sub observer point 
        to the north pole 
    24 Heliocentric range (AU)
    25 Heliocentric range-rate (km/s)
    26 Observer range (AU)
    27 Observer range-rate (km/s)
    28 One-way (down-leg) light-time (minutes)
    29 Sun-Target-Observer ~PHASE angle (deg)
    30 North pole RA (deg)
    31 North pole DEC (deg)


    Keywords
    ----------
    obs_code : [-] str 
        [-] NAIF observatory code obtained from 
        http://www.minorplanetcenter.net/iau/lists/ObsCodesF.html
        500 -> Geocentric observer
        -5 : VLA

    Example
    -------
    Ios ('501', information at 2017-02-17 08:24
    >>> [data,R] = (get_ephem.get_ephemerides('599','2017-02-17 08:24',
                '2017-02-17 08:25','1',-5))[0,2]
    AD = data[0][14] # (") Ang-Diam

    References
    ------------
    https://ssd.jpl.nasa.gov/horizons.cgi

    Notes
    -------
    10/25/2017, CM, update Commit
    """  

    tstart_obj = datetime.strptime(tstart,'%Y-%m-%d %H:%M')
    tstart_UT = datetime.strftime(tstart_obj,"'%Y-%m-%d %H:%M'")
    tend_obj = datetime.strptime(tend,'%Y-%m-%d %H:%M')
    tend_UT = datetime.strftime(tend_obj,"'%Y-%m-%d %H:%M'")


    if nstep == '0': 
        nstep = '1' 
        print('nstep changed to 1')

    http = "http://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1"
    make_ephem = "&MAKE_EPHEM='YES'&TABLE_TYPE='OBSERVER'"
    command    = "&COMMAND=" + (code)
    center     = "&CENTER="+(obs_code)  #500 is geocentric
    t_start    = "&START_TIME=" + tstart_UT
    t_stop     = "&STOP_TIME=" + tend_UT
    t_step     = "&STEP_SIZE='" + nstep + "'"
    quantities = "&QUANTITIES='1,3,4,8,9,12,13,14,15,16,17,19,20,21,24,32'"
    csv        = "&CSV_FORMAT='YES'"

    # 1,2,4,9,10,13,14,19-21,24,29,32

    url = http+make_ephem+command+center+t_start+t_stop+t_step+quantities+csv
    ephem = urlopen( url ).readlines()

    inephem = False
    data = []
    for i in range(len(ephem)) :
        if type(ephem[i]) != str:
            ephem[i] = ephem[i].decode('UTF-8')
        if inephem == False:
            if ephem[i].startswith('Target radii '):
                radii = ephem[i].split(':')[1] 
                radii = radii.split('km')[0]
                R = np.zeros(3) # Radius (km)
                for j in range(3):
                    R[j] = np.float(radii.split('x')[j])
                R = R[[0,2,1]]
            # get observatory lat, lon, alt for later
            if ephem[i].startswith('Center geodetic'):
                l = ephem[i].split(':')[1]
                latlonalt = l.split()[0]
                [lon,lat,alt] = [float(s.strip(', \n')) for s in latlonalt.split(',')]
                observatory_coords = [lat,lon,alt]
            if ephem[i].startswith('$$SOE') :
                inephem=True
                #data = [ephem[i+1].decode('UTF-8').split(',')]
        elif inephem == True:
            if ephem[i].startswith("$$EOE") :
                inephem=False
            else:
                data.append(ephem[i].split(','))
    try:
        out = np.asarray(data)[:,:-1]
        return out, observatory_coords, R
    except:
        sys.exit('ERROR: Ephemeris data not found. Check that the target has valid ephemeris data for the specified time range.')




def read_ephem_line(arr):
    '''Helper to ephemeris.__init__. Converts ephemeris data to float, putting np.nan for "n.a."'''
    arr_out = []
    for s in arr:
        if s.strip() == 'n.a.':
            arr_out.append(np.nan)
        else:
            arr_out.append(float(s))
    return np.asarray(arr_out) 



### Class ###

class Ephemeris():
    '''Functions relevant to an ephemeris for a single target'''
    
    def __init__(self,target,obs_code,tstart,tend,stepsize):
        '''Immediately run get_ephemerides, then set a bunch of 
        class variables corresponding to different information found in
        the ephemeris.
        '''
        self.target = target
        self.obs_code = obs_code
        self.ephem, self.observatory_coords = get_ephemerides(naif_lookup(self.target),self.obs_code,tstart,tend,stepsize)
        self.times = self.ephem[:,0]
        self.sun = [s.strip() for s in self.ephem[:,1]]
        self.moon = [s.strip() for s in self.ephem[:,2]]
        self.ra = self.ephem[:,3] #dms
        self.dec = self.ephem[:,4] #dms
        self.dra = np.asarray([float(s) for s in self.ephem[:,5]]) #arcsec hr-1 
        self.ddec = np.asarray([float(s) for s in self.ephem[:,6]]) #arcsec hr-1
        self.azimuth = np.asarray([float(s) for s in self.ephem[:,7]]) #degrees, North = 0 = 360
        self.elevation = np.asarray([float(s) for s in self.ephem[:,8]]) #degrees, above horizon
        self.airmass = read_ephem_line(self.ephem[:,9]) 
        self.extinction = read_ephem_line(self.ephem[:,10]) #magnitudes; currently not used
        self.vmag = read_ephem_line(self.ephem[:,11]) 
        self.sbrt = read_ephem_line(self.ephem[:,12]) #currently not used
        self.ang_sep = read_ephem_line(self.ephem[:,13]) #arcsec
        self.visibility = [s.strip(' ') for s in self.ephem[:,14]]
        self.ang_diam = read_ephem_line(self.ephem[:,15]) #arcsec
        
        ##planet orientation information only used for plot_planet_system
        self.ob_lon = read_ephem_line(self.ephem[:,16]) #degrees, positive to west
        self.ob_lat = read_ephem_line(self.ephem[:,17]) #degrees
        self.np_ang = read_ephem_line(self.ephem[:,18]) #degrees
        self.np_dist = read_ephem_line(self.ephem[:,19]) #arcsec