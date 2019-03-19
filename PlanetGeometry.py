"""Initialize python repository for pyPR geometry.

Notes
-----
11/17/17,CM, Initial Commit
"""

import sys, os
from time import strftime, gmtime, time
from datetime import datetime,timedelta

import numpy as np
import scipy.ndimage 
import math 

from astropy import constants as cst
from astropy.io import fits

import matplotlib.pyplot as plt
import importlib
import time

# Debug 
import warnings
import IPython




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

    naif = 'naif_id_table.txt'

    try:
        open(naif)
    except FileNotFoundError: 
        try: 
            naif = 'pypr/' + naif
            open(naif)
        except:  
            print(os.getcwd())
            sys.exit('Wrong folder, make sure you are a folder above pyPR') 

    target = target.upper().strip(', \n')
    with open(naif)  as f:
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

def get_ephemerides(code, tstart, tend , nstep , obs_code = '500') :
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
        [km] Radius of target body {Equator, meridian, pole}


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

    from urllib.request import urlopen
    import urllib 

    # IF no end time is provided, simulate batch mode 
    try: 
        tstart_obj = datetime.strptime(tstart,'%Y-%m-%d %H:%M:%S.%f')
        tstart_UT = datetime.strftime(tstart_obj,"'%Y-%m-%d %H:%M:%S.%f'")
        tend_obj = datetime.strptime(tend,'%Y-%m-%d %H:%M:%S.%f')
        tend_UT = datetime.strftime(tend_obj,"'%Y-%m-%d %H:%M:%S.%f'")
    except ValueError: 
        try: 
            tstart_obj = datetime.strptime(tstart,'%Y-%m-%d %H:%M:%S')
            tstart_UT = datetime.strftime(tstart_obj,"'%Y-%m-%d %H:%M:%S'")
            tend_obj = datetime.strptime(tend,'%Y-%m-%d %H:%M:%S')
            tend_UT = datetime.strftime(tend_obj,"'%Y-%m-%d %H:%M:%S'")
        except ValueError: 
            tstart_obj = datetime.strptime(tstart,'%Y-%m-%d %H:%M')
            tstart_UT = datetime.strftime(tstart_obj,"'%Y-%m-%d %H:%M'")
            tend_obj = datetime.strptime(tend,'%Y-%m-%d %H:%M')
            tend_UT = datetime.strftime(tend_obj,"'%Y-%m-%d %H:%M'")

    if tstart_UT > tend_UT: 
        sys.exit('Start time is larger than end time.')

    if nstep == '0': 
        nstep = '1' 
        print('nstep changed to 1')

    http       = "https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1"
    make_ephem = "&MAKE_EPHEM='YES'&TABLE_TYPE='OBSERVER'"
    command    = "&COMMAND=" + (code)
    center     = "&CENTER="+(obs_code)  #500 is geocentric
    t_start    = "&START_TIME=" + tstart_UT
    t_stop     = "&STOP_TIME=" + tend_UT
    t_step     = "&STEP_SIZE='" + nstep + "'"
    quantities = "&QUANTITIES='1,3,4,8,9,12,13,14,15,16,17,19,20,21,24,32'"
    csv        = "&CSV_FORMAT='YES'"

    # 1,2,4,9,10,13,14,19-21,24,29,32

    url = urllib.parse.quote(http+make_ephem+command+center+t_start+t_stop+t_step+quantities+csv, safe=':/?=&,')    
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


def ephem_uncertainties(target,orange):
    """Time averaged unceratinties in the ephemeris. 
    
    Based on a paper by Will Folkner, at JPL, the time-averaged 
    uncertainties are return. For more details, the paper should be 
    consulted. As the time average is dependent on various factors, 
    such as range, spacecraft tracking availablity etc 
    
    Parameters
    ----------
    target : str
        [-] Name of the planet 
    ra : float
        [deg] Right ascension 
    dec : float
        [deg] Declination 
    orange : float
        [m] Range to the planet 
    
    Keyword Arguments
    ----------

    
    Returns
    -------
    dra : float
        [deg] Right ascension uncertainty
    ddec : float
        [deg] Declination uncertainty
    dorange : float
        [m] Range to the planet uncertainty

    
    Warnings
    -------
    
    
    Example
    -------
    
    >>>function(arg1, arg2)
    out1
    
    References
    ------------
    Folkner, W. M. "Uncertainties in the JPL planetary ephemeris." 
    Proceedings of the Journees. 2010. 
    https://syrte.obspm.fr/journees2010/pdf/Folkner.pdf

    Todo
    ----- 
    
    Notes
    -------
    12/01/17, CM, Initial commit
    """

    # Values as obtained from the paper 
    dra_paper = np.array([1,0.4,5,20,3,250,1000,1000])*1e3
    ddec_paper = np.array([1,0.2,2,20,5,100,400,1000])*1e3
    dorange_paper = np.array([0.5,0.1,1,5,0.25,150,500,2000])*1e3

    if target.casefold().strip() == 'Mercury'.casefold():
        pl = 0
    elif target.casefold().strip() == 'Venus'.casefold():
        pl = 1
    elif target.casefold().strip() == 'Mars'.casefold():
        pl = 2
    elif target.casefold().strip() == 'Jupiter'.casefold():
        pl = 3
    elif target.casefold().strip() == 'Saturn'.casefold():
        pl = 4
    elif target.casefold().strip() == 'Uranus'.casefold():
        pl = 5
    elif target.casefold().strip() == 'Neptune'.casefold():
        pl = 6
    elif target.casefold().strip() == 'Pluto'.casefold():
        pl = 7
    else:
        warnings.warn('Body not found. No uncertainties prescribed. Evrything should still work')
        dra,ddec,dorange = 0,0,0
        return [dra,ddec,dorange]

    dra = np.degrees(dra_paper[pl]/orange) # Right ascension uncertainty
    ddec = np.degrees(ddec_paper[pl]/orange) # Declination uncertainty 
    dorange = dorange_paper[pl] # (m) Range uncertainty
        
    return  [dra,ddec,dorange]



def read_ephem_line(arr):
    '''Helper to ephemeris.__init__. Converts ephemeris data to float, putting np.nan for "n.a."'''
    arr_out = []
    for s in arr:
        if s.strip() == 'n.a.':
            arr_out.append(np.nan)
        else:
            arr_out.append(float(s))
    return np.asarray(arr_out) 




# Make an oblate spheroid 
def triaxialellipsoid(R,x,y): 
    """Create an oblate spheroid 

    Create a 3D ellipsoid with a,b,c referring to x,y,z

    # Note that this is a difference convention. Z is in this case chosen 
    to come out of the plane. a,b,c and refer to the standard ellipsoid axis 
    defined by x,a out of plane, z,c upwards, y,b completes the frame

      y,c   ^ 
            |
            |
            |
             ------> x,b
           / 
          /
         /
      z,c 

    Parameters
    -------
    R [3x1] float 
        [m] a,b,c three semi major axis of the system 
    x [1] float 
        [m] x axis on which the the system will be build  
    y [1] float 
        [m] y axis 
    


    Returns
    ----------
    x,y,z [1xN] float 
        [m] Three dimensional distribution  

    Keywords
    ----------


    Example
    -------

    References
    ------------
    Triaxial ellipsoid http://mathworld.wolfram.com/Ellipsoid.html

    Notes
    -------
    10/24/2017, CM, Initial Commit
    """

    # x = np.linspace(-1,1,N)*a
    # y = np.linspace(-1,1,N)*b

    a,b,c = R

    # warnings.filterwarnings("ignore")
    xv, yv = np.meshgrid(x, y) 

    # Create the surface 
    surf = 1. - xv**2/b**2 - yv**2/c**2
    surf[surf<0] = 0.0
    zv= a*np.sqrt(surf)

    # remove nan outside of the planet 
    zv[~np.isfinite(zv)] = 0
    
    return(xv,yv,zv)





# Read in the Jupiter brightness map 
def brightnessmap(R,r_pla, x,y,T_dab,p, conv = 0.01): 
    """Create a brightness map of Jupiter

    Based on x,y and grid, create a brightness map of Jupiter. Iterative 
    solution that solves for the brightness distribution 

    Parameters
    -------
    x : [1xN] float
        [R] Normlized x - coordinate  
    y : [1xN] float
        [R] Normlized y - coordinate  
    T_dab: [1] float
        Disk averaged brightness temperature 
    R : [3x1] float 
        [pix] Radius of the planet [equator,polar,equator]
    p : [1x2] float 
        [EW,NS] Limb darkening coefficient 



    Returns
    ----------
    T : [NxN] float 
        [K or Jy] Brightness distribution on the sky 

    Keywords
    ----------
    Jy : Boolean
        If true, output in Jansky 
        If false, output in K 

    Example
    -------
    >>> R = np.array([72e6,66e6,72e6])
    >>> R = R/R[0]
    >>> N = 200 
    >>> T_dab = 350
    >>> p = np.array([0.08,0.08])
    BrightnessMap(R,N,T_dab,p)

    References
    ------------

    Notes
    -------
    10/24/2017, CM, Initial Commit
    """
    R = R/R[0]*r_pla
    # Call the three dimension surface and return value in pixels 
    (xv,yv,zv) = triaxialellipsoid(R,x,y)
    # Avoid singularities 
    zv[zv==0] = float('nan')
    # Obtain the emission angle 
    th = np.arccos(zv/np.sqrt(xv**2+yv**2+zv**2))
    

    # Where the value is nan, the model should be zero 
    th[~np.isfinite(th)] = np.pi/2.
    zv[~np.isfinite(zv)] = 0.

    # Iterate to make the model agree to the disk averaged 
    # brightness temperature 
    pixels = np.ones_like(zv) 
    cond = 1 
    T_scale = T_dab

    # Build a uniform disk model 
    if (p == [0,0]).all() or (p == 0).all() : 
        T = np.ones_like(zv)*T_dab
        T[zv==0] = 0
        return T
    # Build a limb-darkened disk model 
    else: 
        # IPython.embed()
        while cond: 
            # Avoiding numeric issues 
            cos = np.cos(th)
            cos[cos<1e-5] = 0
            # Exponent p, should always be smaller than p max 
            pexp = (p[0]+(np.abs(yv)/R[2])*(p[1]-p[0]))
            pexp[pexp > np.max(p)] = np.max(p)
            pexp[pexp < np.min(p)] = np.min(p)
            T = T_scale*cos**(pexp)
            # IPython.embed()
            T_model = (np.sum(T)/np.sum(pixels[zv>0.0])) 
            if np.abs(T_dab - T_model) < conv: 
                cond = 0 
            else: 
                T_scale = T_scale + (T_dab - T_model)
                # print('Disk averaged brighntess temperature\n',
                #     'Target T: {:3.2f}'.format(T_dab), \
                #     'Computed T: {:3.2f}'.format(T_model), \
                #     'Peak T: {:3.2f} '.format(T_scale))

        print('Disk averaged brighntess temperature\n',
            'Target T: {:3.2f}'.format(T_dab), \
            'Computed T: {:3.2f}'.format(T_model), \
            'Peak T: {:3.2f} '.format(T_scale))

        T[zv<1e-5] = 0
        return T

def emissionangle(lat_c, lat_d, ob_lat_d, orange, R): 
    """Calculate the zonally averaged emission angle based on the 
    goedetic and geocentric latitudes 

    Based on x,y and grid, create a brightness map of Jupiter. Iterative 
    solution that solves for the brightness distribution 

    Parameters
    -------
 



    Returns
    ----------


    Keywords
    ----------


    Example
    -------
    lat_c =(np.arange(-np.pi/2,np.pi/2,np.pi/179))
    ob_lat_d = np.radians(-3.27) 
    R = np.array([71492., 71492., 66854.])
    orange = 8.0491221e+11  # [m] 5.4 AU
    # Estimate the distance to the center as: 
    r_l = R[2]+(R[0]-R[2])*np.cos(ob_lat_d)
    f = 1 - R[2]/R[0] 
    lat_d = geoc2geod( orange+r_l, lat_c, f )[0]   
    alpha, lat = emissionangle(lat_c, lat_d, ob_lat_d, orange, R)

    References
    ------------

    Notes
    -------
    10/24/2017, CM, Initial Commit
    """

    # [rad] Convert observer geodetic latitude to geocentric latitude
    ob_lat_c = geod2geoc(ob_lat_d, orange, R[0], 1-R[2]/R[0]) 


    # Find the angle between the Center-Observer and the point to the viewing point 
    lat_p = np.radians(np.arange(-89,90,1)) 
    lat = lat_p+ob_lat_c


    # [m] Find local radius for a given 
    r_s = R[0]*R[2]/(np.sqrt((R[2]*np.cos(lat_p))**2 + (R[0]*np.sin(lat_p))**2))

    # [m] Calculate the surface intercept distance from observer using cosine law
    s = np.sqrt(r_s**2 + orange**2 - 2*r_s*orange*np.cos(np.radians(ob_lat_c)))

    # [deg] Calculate the angle beta between local tangent and the local radius 
    beta = np.abs(lat_d - lat_c)

    # Calculate the angle of the Triangle, Center-Observer-Surface intercept 
    gam = np.pi/2 - np.arcsin(orange/s*np.sin(np.abs(lat_p))) - beta

    # Calculate the emission angle 
    alpha = np.pi/2 - gam 

    return alpha, lat


def zonal_residual_model(R, r_pla, x_pix, y_pix, p, ob_lat_d, orange, T_res, th_res, conv = 0.01, plotting = False ): 
    """Create a brightness map with structure of Jupiter 

    Based on x,y and grid, create a brightness map of Jupiter. Iterative 
    solution that solves for the brightness distribution 

    Parameters
    -------


    R : [3x1] float
        [R] Principal axis of the ellipsoid [equator,polar,equator] 
        Rotated into correct view 
    r_pla : [1] 
        [pix] size of the planet in pixel  
    x_pix: [N] int
        [pix] x_axis for the model (redundant, could be inferred from r_pla)
    y_pix: [N] int
        [pix] y_axis for the model  (redundant, could be inferred from r_pla)
    p : [1x2] float 
        [EW,NS] Limb darkening coefficient 
    ob_lat_d: [1] 
        [rad] Sub observer latitude 
    orange: [1] 
        [m] Observer range 
    T_res : [1xN] float
        [K] Brightness temperature residual, zonally averaged, limb-darkened 
    th_res : [1xN] float
        [rad] Corresponding latitudes for the residual temperatures 
    

    Returns
    ----------
    T : [NxN] float 
        [K or Jy] Brightness distribution on the sky 

    Keywords
    ----------
    conv = 0.01 : float
        Converge criteria for limb darkening disk computation 
        A high value will create a disk very close to the peak temperature. 
        A low value will iterate until the disk averaged temperature matches the required value 

    Example
    -------
    >>> R = np.array([72e6,66e6,72e6])
    >>> R = R/R[0]
    >>> N = 200 
    >>> T_dab = 350
    >>> p = np.array([0.08,0.08])
    BrightnessMap(R,N,T_dab,p)

    References
    ------------

    Notes
    -------
    10/24/2017, CM, Initial Commit
    """

    
    from mpl_toolkits.basemap import Basemap, addcyclic
    import scipy.interpolate

    # Correct for the limb darkening given that the residuals are on 
    # top of a limb darkened disk. 
    # ------------------------------------------------------------------
    


    # Obtain emission angle. 
    lat_c =(np.arange(-np.pi/2,np.pi/2,np.pi/179))
    flat = 1 - R[2]/R[0] 
    lat_d = geoc2geod( orange + R[2]+(R[0]-R[2])*np.cos(ob_lat_d), lat_c, flat )[0]  
    alpha, lat_emission = emissionangle(lat_c, lat_d, (ob_lat_d), orange, R)

    # Input residual striructure in limb darkened. 

    # Pad residual with 0 to 

    f = scipy.interpolate.interp1d(np.degrees(th_res), T_res, bounds_error=False,fill_value=np.mean(T_res))
    lats = np.arange(-89,90,2)
    lons =  np.arange(0,360,2)
    T = f(lats) - np.mean(f(lats)) # Residual structure after and levelled to zero 

    # Detrend the residuals from emission angle by diving by cos(theta)**p
    f = scipy.interpolate.interp1d(np.degrees(lat_emission), alpha, bounds_error=False, fill_value=0)
    T_db = T/np.cos(f(lats))**p[0] 


    # Create a map due zonally averaged values 
    T_map = np.transpose(np.tile(T_db,(180,1)))

    # Bug in Basemaps
    try: 
        T_map, lons = addcyclic(T_map, lons)
    except: 
        print('Basemap bug has not yet been removed')

    # create Basemap instance for orthographic projection.
    m = Basemap(projection='ortho',lat_0=ob_lat_d
                ,lon_0=90)

    # compute map projection coordinates for lat/lon grid.
    x, y = m(*np.meshgrid(lons,lats))
    lonpt, latpt = m(x,y,inverse=True)


    # make filled contour plot.
    if plotting: 
        plt.figure()
        cs = m.contourf(x,y,T_map,30,cmap=plt.cm.jet)
        plt.colorbar(cs)
        m.drawmapboundary() # draw a line around the map region
        m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
        m.drawmeridians(np.arange(0.,420.,60.),labels=[0,0,0,1]) # draw meridians
        plt.title('Orthographic Projection') # add a title
        plt.show() 

    
    # Projected axis of ellipsod R[1], and R[2]

    # Call the three dimension surface 
    (xv,yv,zv) = triaxialellipsoid(R/R[0]*r_pla,x_pix,y_pix)
    # Avoid singularities 
    zv[zv==0] = float('nan')
    # Obtain the emission angle 
    th = np.arccos(zv/np.sqrt(xv**2+yv**2+zv**2))
    
    # Where the value is nan, the model should be zero 
    th[~np.isfinite(th)] = np.pi/2.
    zv[~np.isfinite(zv)] = 0.

    # Interpolate the structure onto the same grid as the limb darkened disk 
    # Note that the ld model is in pixel, so we covert them here  
    x_temp = x.flatten()/m.rmajor # [R_E] units of RE
    x_t = x_temp[x_temp<2.0]*r_pla # [pix] Multiply with R_E
    y_temp = y.flatten()/m.rmajor
    y_t = y_temp[x_temp<2.0]*R[2]/R[0]*r_pla # [pix] Multiply with R_P 
    T_temp = T_map.flatten()
    T_t = T_temp[x_temp<2.0]

    # Center the array so that 0,0 refers to the center of the planet 
    x_t -= r_pla 
    y_t -= R[2]/R[0]*r_pla
    p_in=(np.stack((x_t,y_t))).transpose()

    # Interpolate on the output 
    T_resmap = scipy.interpolate.griddata(p_in,T_t,(xv.T, yv.T), fill_value = 0, method='cubic')

    # plotting the residuals 
    if plotting: 
        fig, axs = fig, axs = plt.subplots(1, 1)
        cs = plt.contourf(xv.T, yv.T,T_resmap,30,cmap=plt.cm.jet)
        plt.title('Zonal residual NOT limb darkened')
        plt.colorbar(cs)
        axs.set_aspect('equal', 'box')
        plt.show()

    # Element wise multiplication of theta with T_resmap to limb-darken the residuals 
                # Avoiding numeric issues 
    cos = np.cos(th)
    cos[cos<1e-5] = 0
    # Exponent p, should always be smaller than p max 
    pexp = (p[0]+(np.abs(yv)*r_pla/R[2])*(p[1]-p[0]))
    pexp[pexp > np.max(p)] = np.max(p)
    pexp[pexp < np.min(p)] = np.min(p)
    T_res_ld = T_resmap*cos**(pexp)

    # 
    # T_res_ld = np.multiply(T_resmap,np.cos(th)**np.mean(p)) 

    # Normalize so that the structure has zero flux 
    pixels = np.ones_like(zv) # Pixels that contain flux 
    T_res_ld -= (np.nansum(T_res_ld)/np.sum(pixels[zv>0.0])) 
    if plotting: 
        fig, axs = fig, axs = plt.subplots(1, 1)
        cs = plt.contourf(xv.T, yv.T,T_res_ld,30,cmap=plt.cm.jet)
        plt.title('Zonal residual limb darkened')
        axs.set_aspect('equal', 'box')
        plt.colorbar(cs)
        plt.show()

    # Subtraction forces off-planet residuals below zero 
    T_res_ld[zv<0.1] = 0 
    
    return T_res_ld.T


# Conversions 

def jy2tb(flux, solidangle, freq):
    """Convert integrated flux to disk averaged brightness temperature 

    The integrated flux can be obtained by looking at Amp vs UVdistance 
    and interpolating to zero spacing alternatively, you can integrate a 
    map given in Jy/beam. The easiest way to do is by converting to 
    Jy/pixel (using the NRAO formula https://science.nrao.edu/facilities/vla/proposing/TBconv)
     and then integrating the resulting map 

    Parameters
    -------
    flux: [1] float 
        [Jy] integrated 
    solidangle: [1] float
        [strad] Solid angle of the planet. Calculated when calling the 
                ephemeris 
    freq: [1] float
        [Hz] Center frequency


    Returns
    ----------
    T : [1] float 
        [K] Disk averaged brightness temperature 


    Keywords
    ----------


    Example
    -------
    >>> jy2tb(13.431,Pl.sa,220) 
    92.45

    References
    ------------
    http://starlink.eao.hawaii.edu/docs/sun213.htx/sun213se9.html 


    Notes
    -------
    7/31/2018, CM, Initial Commit
    """   

    Jy = 1e-26
    lam = cst.c.value/freq

    T_da = flux*Jy*lam**2/(2*cst.k_B.value)/solidangle
    

    return T_da


def tb2jy(T,nu,pixscale):
    """Convert brightness temperature to pixel 



    Parameters
    -------
    T : [NxN] float 
        [K] Brightness model 
    nu : [1] float 
        [Hz] Frequency of obseravations
    pixscale : [1] float 
        [arc/pixel] Pixel scale, how many arcsec is one pixel 



    Returns
    ----------
    

    Keywords
    ----------


    Example
    -------
    >>> tb2jy

    References
    ------------

    Notes
    -------
    10/25/2017, CM, Initial Commit
    """   

    c = cst.c.value # (m/s) Speed of light 
    kb = cst.k_B.value # (m^2kg s^-2K^-1)
    jy = 1e-26 

    # Convert into flux scale 
    # flux = T*2.*kb*nu**2/c**2/jy*np.pi*arcsec2rad(pixscale)**2/(4*np.log(2))
    flux = T*2.*kb*nu**2/c**2/jy*arcsec2rad(pixscale)**2

    return flux


def arcsec2rad(arcsec): 
    ''' Convert arcseconds to radians


    Parameters
    -------

    arcsec : [1] float 
        [arcsec] value in arcsencs 



    Returns
    ----------
    rad : [1]f float 
        [rad] Value in radians 

    Keywords
    ----------


    Example
    -------
    >>> arcsec2rad

    References
    ------------

    Notes
    -------
    10/24/2017, CM, Initial Commit
    '''
    return arcsec/3600.*np.pi/180.


def lldvis(mu,asize,lam, ): 
    '''Computes the visibility function of a linear limbdarkened disk 
    I(theta) = I(0)*(1 - mu*(1-cos(theta))) 


    
    
    Parameters
    -------
    asize : [1] int
        [arcsecond] Size of the planet in arcseconds    


    Returns
    ----------
    firstnull: [] flpat 
        [] The first null in wavenumbers 
    

    Keywords
    ----------


    Example
    -------
    mu = 0.01; asize = arcsec2rad(50); lam = 0.3/8  
    Vsq,B = ldvis(mu,asize,lam)
    uvwave = B/lam
    plt.plot(uvwave,Vsq)

    References
    ------------
    Hanburry Brown - The effect of Limb darkened disk 
    Baines - Testing limb darkening laws using NPOI observatons  

    Notes
    -------
    08/06/2018, CM, Initial Commit
    ''' 




    from scipy.special import jv

    B = np.arange(1,1000)
    x = np.pi*B*asize/lam 

    Vsq = ((1.0-mu)/2. + mu/3.)**(-2)*((1.0-mu)*jv(1,x)/x + mu*(np.pi/2)**0.5*jv(1.5,x)/x**1.5)**2 

    return Vsq, B 


def cldvis(k,amp,lam,B,a): 
    '''Computes the visibility function of a cos(theta)**p limb-darkened disk   

    
    
    Parameters
    -------
    k : [1] float
        [-] Limb darkening parameter to be plotted 
    amp : [1] float
        [Jy] zero baseline flux amplitude 
    lam :[1] float 
        [m] Wavelength of the observation 
    B : [N] float
        [m] Baseline array      
    a : [1] float 
        [arcseconds] size of the planet (diameter) 

    Returns
    ----------
    V_abs: [N] flpat 
        [Jy] Visibility function 
    

    Keywords
    ----------


    Example
    -------
    k = 0.1; amp = 13.5; lam = 0.3/224; B = np.arange(0,500); a = np.radians(2.2/3600) 
    V = ldvis2(k,amp,lam,B,a)
    plt.plot(B,V,label='limb darkened: k = {:2.2f}'.format(k)) 
    plt.legend()
    plt.show()

    References
    ------------
    http://www.iue.tuwien.ac.at/phd/minixhofer/node59.html

    Notes
    -------
    08/08/2018, CM, Initial Commit
    '''

    from scipy.special import jv  
    from scipy import integrate 
    import matplotlib.pyplot as plt  

    rho = B/2./lam

    # val = np.abs((2*jv(1,2*np.pi*rho*a)/(2*np.pi*rho*a)))
    integral_val = np.zeros_like(rho)
    for i in range(len(rho)): 
        f = lambda r: r*np.sqrt(1-r**2/a**2)**k*jv(0,2*np.pi*rho[i]*r) 
        integral_val[i] = np.abs(2*np.pi*scipy.integrate.quad(f,0,a)[0] )


    # Normalize and scale to  
    integral_val = amp*integral_val/integral_val[0]


    return integral_val


def fringesize(asize):
    '''Approximate the bessel's function first null for a uniform disk  

    
    
    Parameters
    -------
    asize : [1] int
        [arcsecond] Size of the planet in arcseconds    


    Returns
    ----------
    firstnull: [] flpat 
        [] The first null in wavenumbers 
    gen

    Keywords
    ----------


    Example
    -------
    >>> 

    References
    ------------
    Basic Radio Interferometry – Geometry
    Rick Perley, NRAO/Socorro

    Notes
    -------
    08/06/2018, CM, Initial Commit
    '''      

    rad2arc = 1./arcsec2rad(1) # [] arceseconds per radians 
    print('The first null can be found at {:2.2}'.format(str(rad2arc/asize))) 

    return rad2arc/asize



def axisforcasamodel(imsize, planetsize, radius=False): 
    """Obtain the axis for casa, based on the input supplied to tclean 

    Beamsize parameters can be obtained from: 
    https://science.nrao.edu/facilities/vla/docs/manuals/
    oss/performance/resolution
    
    Parameters
    -------
    imsize : [1] int
        [-] Size of the pixels   
    planetsize : [1] int
        [arcsec] size of planet in pixel



    Returns
    ----------
    

    Keywords
    ----------
    radius: [1] Boolean 
        [R] Return the axis in radii  

    Example
    -------
    >>> axisforcasamodel(420,39,10,3.1)
    References
    ------------

    Notes
    -------
    10/24/2017, CM, Initial Commit
    """    
    

    # Radius of Jupiter in pixels 

    r_pla = planetsize/2

    # Compute the padding left and right 
    x = np.linspace(-imsize/2.,imsize/2.,imsize)
    y = x 


    if radius: 
        x,y = x/r_pla,y/r_pla

    return x,y  


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def geo2eci(r,lon,lat,):
    """Conversion from geocentric coords to eci coords. 

    Parameters
    ----------
    r : float
        [m] Radius from center 
    lon : float
        [rad] Longitude
    lat : float
        [rad] Latitude


    Keyword Arguments
    ----------


    Returns
    -------
    r : [3x1] float
        [m,m,m] Position in rectengular coordinates

    Warnings
    -------
    Assuming spherical Earth rather than WGS Spheroid. To Be Added

    Example
    -------

    References
    ------------

    Notes
    -------
    08/30/18, SB, Initial commit, Geocentric only
    """
    


    # Convert latitude to spherical polar angle
    th = np.pi/2. - lat

    r_eci = r*np.array([
                        np.sin(th)*np.cos(lon),
                        np.sin(th)*np.sin(lon),
                        np.cos(th)
                       ])
    return r_eci

def geoc2geod( r, lat_c, f ): 
    """Conversion from geocentric to geodetic. '
    
    
    Planetocentric: Defined with respect to a sphere 
        Longitude: Positive eastwards 
        Latitude: Positive northward 
    Planetodetic: Defined with respect to the local tangent plane
        Longitude: Positive eastwards   
        Latitude: Positive northward 
    Planetographic: Defined with respect to the local tangent plane, 
    and longitude is increasing with time for the observer 
        Longitude: Positive WESTwards   
        Latitude: Positive northward 
    
    Source: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials
    /pdf/individual_docs/17_frames_and_coordinate_systems.pdf (slide 23-25)


    Parameters
    ----------
    r : float
        [m] Radius from center 
    lat_c : float
        [rad] Geocentric Latitude
    f: float 
        [-] Flattening parameter (1-R_p/R_e) or (R_e-R_p)/R_e


    Keyword Arguments
    ----------


    Returns
    -------
    lat_d : float
        [rad] Geodetic Latitude
    h: float 
        [m] Height above the surface  

    Warnings
    -------
    Assuming spherical Earth rather than WGS Spheroid. To Be Added

    Example
    -------
    import pyPR.PlanetGeometry as PG
    lat_c = np.radians(48.66594); f = 1/298.257223563
    r = 6366453 # 400 m above Paris
    np.degrees(PG.geoc2geod( r, lat_c, f ))[0] # 48.8563

    

    References
    ------------
    https://www.mathworks.com/help/aeroblks/geocentrictogeodeticlatitude.html 

    Notes
    -------
    12/11/18, CM, Initial commit, Geocentric only
    """
    # 
    # Find the coordinates of the surface intersection point  
    x_s = (1-f)*r/(np.sqrt(np.tan(lat_c)**2 + (1 - f)**2))
    y_s = np.sqrt(r**2 - x_s**2)*(1-f) 

    # Find the angle between the surface intersect and the point where 
    # the surface tanget crosses the equator 
    lat_s = np.arctan(np.tan(lat_c)/(1-f)**2) 

    # Check for sign due atan ambiguity  
    lat_s = (np.sign(lat_c)/np.sign(lat_s))*lat_s

    # Calculate the local radius for a given ellipsoid 
    r_s = x_s/np.cos(lat_c) 

    # Distance from surface intercept to observer 
    l = r - r_s 

    # Difference in latitude at surface intercept 
    dlat_s = lat_s - lat_c 

    # Calculate the mean altitude of the observer above the surface 
    h = l*np.cos(dlat_s)

    # Curvate of the surface at the intercept to calculate the local tangent 
    rho_s = r*(1-f)**2/(1-(2*f-f**2)*np.sin(lat_s)**2)**(3/2)

    # Difference between the two latitude 
    dlat = np.arctan(l*np.sin(dlat_s)/(rho_s + h))

    # Calculate the geodetic latitude 
    lat_d = lat_s - dlat

    return lat_d, h

def geod2geoc( lat_d, h, r_e,  f ): 
    """Conversion from geodetic to geocentric. '
    

    
    Planetocentric: Defined with respect to a sphere 
        Longitude: Positive eastwards 
        Latitude: Positive northward 
    Planetodetic: Defined with respect to the local tangent plane
        Longitude: Positive eastwards   
        Latitude: Positive northward 
    Planetographic: Defined with respect to the local tangent plane, 
    and longitude is increasing with time for the observer 
        Longitude: Positive WESTwards   
        Latitude: Positive northward 
    
    Source: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials
    /pdf/individual_docs/17_frames_and_coordinate_systems.pdf (slide 23-25)


    Parameters
    ----------
    h : float
        [m] height of the observer 
    r_e : float
        [m] equatorial radius of the planet 
    lat_d : float
        [rad] Geodetic Latitude
    flat: float 
        [-] Flattening parameter (1-R_p/R_e) or (R_e-R_p)/R_e


    Keyword Arguments
    ----------


    Returns
    -------
    lat_c : float
        [rad] Geocentric Latitude

    Warnings
    -------
    Assuming spherical Earth rather than WGS Spheroid. To Be Added

    Example
    -------
    import pyPR.PlanetGeometry as PG 
    # Earth WGS - 84, Paris 
    f = 1/298.257223563; lat_d = np.radians(48.8567); r_e = 6378.137e3; h = 400
    np.degrees(PG.geod2geoc(h, r_e, lat_d, f)) # 48.665943820540754



    References
    ------------
    https://www.mathworks.com/help/aeroblks/geodetictogeocentriclatitude.html

    Notes
    -------
    12/11/18, CM, Initial commit, Geocentric only
    """

    # Calculate the geocentric latitude at the surface intercept 
    lat_s = np.arctan((1-f)**2*np.tan(lat_d))
    
    lat_c = (np.arctan((h*np.sin(lat_d) + r_e*np.sin(lat_s))/
                      (h*np.cos(lat_d) + r_e*np.cos(lat_s))))

    return lat_c


def rotateprincipalaxis(R, ob_lat): 
    """Obtain the principal axis of an ellopsoid seen at an angle. 

    Applicable for a 2D body, such as prolate spheroid

    When an ellipsoid is seen off the equator, the principal axis change. 
    A planet seen face on, is squished due to rotational flattening, 
    and the principal axis are R_equator, and R_pole. With increasing 
    observer latitude, the projection becomes more cirular, until we 
    observe the planet above the pole as a perfect sphere. 

    Input maybe be in any unit! 

    
    Parameters
    -------
    R : [3x1] float
        [m] Principal rotation axis of the planet (equator, pole, equator)
    obs_lat : [1] float
        [deg] Observer's latitude 



    Returns
    ----------
    R_rot: [3x1] float 
        [m]: Rotated principal axis 

    Keywords
    ----------
 

    Example
    -------
    >>> import pyPR.PlanetGeometry as PG
    >>> PG.rotateprincipalaxis(np.array([71492., 71492., 66854.]),np.radians(-3.27))
    References
    ------------
    https://www.geometrictools.com/Documentation/PerspectiveProjectionEllipsoid.pdf 

    Notes
    -------
    12/10/2018 CM, Initial Commit
    """    
    
   
    # Rotate around the y axis to correct for the sub latitude 
    # Note that the longitude should not matter, for an rotational ellipsoid.
    # (pixel) Polar radius for triaxial ellipsoid

    a = (R[1]*np.cos(np.radians(ob_lat))**2 
            + R[2]*np.sin(np.radians(ob_lat))**2 )
    # (pixel) Equatorial radius for triaxial ellipsoid 
    b = (R[2]*np.cos(np.radians(ob_lat))**2 
            + R[1]*np.sin(np.radians(ob_lat))**2)

    return [a,b]

def rotateprincipalaxis_3D(R, ob_lat_d, ob_lon, ob_range): 
    """Obtain the principal axis of an ellopsoid seen at an angle. 
    Applicable for a 3D body 

    When an ellipsoid is seen off the equator, the principal axis change. 
    A planet seen face on, is squished due to rotational flattening, 
    and the principal axis are R_equator, and R_pole. With increasing 
    observer latitude, the projection becomes more cirular, until we 
    observe the planet above the pole as a perfect sphere. 

    Input maybe be in any unit! 

    
    Parameters
    -------
    R : [3x1] float
        [m] Principal rotation axis of the planet (equator, pole, equator)
    obs_range : [1] float
        [m] Observer's distance from surface 
            Best obtained by calculating the distance of light travel time 
            from NAIF P.lt*c

    Remember that ephemeris is in degree 

    Returns
    ----------
    R_rot: [3x1] float 
        [m]: Rotated principal axis 

    Keywords
    ----------
    

    Example
    -------
    # Example by running ephem on some Jupiter data 
    R = np.array([71492., 71492., 66854.]) 
    ob_range = 8.0491221e+11
    ob_lat_d = np.radians(-3.267)
    ob_lon = np.radians(73.12)
    axis , center = rotateprincipalaxis_3D(R, ob_lat_d, ob_lon, ob_range)

    Warnings
    -------
    It might be avisable to scale the distances.  Test first unscaled 

    References
    ------------
    https://www.geometrictools.com/Documentation/PerspectiveProjectionEllipsoid.pdf
    https://math.stackexchange.com/questions/2243974/perspective-projection-of-ellipsoid-to-ellipse-solving-for-standard-form-ellips 
    
    Notes
    -------
    12/10/2018 CM, Initial Commit
    """ 

    import numpy.linalg as npl

    # Given the ob_lat, ob_long, and the range we can compute the observer 
    # Horizons: Obrv range, obs sub-lon, obs sub_lat, f from radius 
    # pypr: lt*c, ob_lon, ob_lat

    # There is a problem for negative latitudes 
    # ob_lat_d = np.abs(ob_lat_d)

    # Flattening parameter of the planet 
    f = (R[0]-R[2])/R[0] 
    
    # [rad] Convert into geocentric latitude based on altitude, radius and latitude
    ob_lat_c = geod2geoc( ob_lat_d, ob_range, R[0], f )  

    # [m] Obtain local radius 
    r_s = R[0]*R[2]/(np.sqrt((R[2]*np.cos(ob_lat_c))**2 + (R[0]*np.sin(ob_lat_c))**2))
 
    # Obtain the coordinates of the observer in an Jovian-centered frame 
    r_J = geo2eci(ob_range+r_s,ob_lon,ob_lat_c).reshape(3,1)
    # Normalized vector to define projection plane 
    n = -r_J/npl.norm(r_J) 

    # A plane is defined by a point and a normal (Use the y intercept with the body 
    # Plane not going through origin, but through point x_s 
    x_s = np.array([0,0,0]) 
    lam = np.dot(n.T,x_s)

    # Ellipsoid is defined as x.T*A*x + b*x + c = 0 

    # Set up the quadratic form of the ellipsoid equation 
    # http://mathworld.wolfram.com/QuadraticSurface.html
    A = np.diag([1./R[0]**2,1./R[1]**2,1./R[2]**2])
    c = 1 

    # Obtaining a repeated real root for the intersection of the ray 
    # eminating from the viewing point assures that we are the tangent 
    # point. Two distinct roots mean we are crossing the spheroid 

    M = (2*A@r_J)@(2*A@r_J).T - 4*(r_J.T@A@r_J + c)*A


    # Use temporary variables in the view plane to convert to a standard ellipse 
    # r_J_hat is normal vector of the viewing plane. In this case the viewing plane 
    # is aligned 
    # u,v, and n build a orthonormal basis 
    # make a random vector, find the perpendicular part to n and scale it to 1 
    u = np.random.randn(3)
    u -= np.reshape((u.dot(n) * n / npl.norm(n)**2),(3,))
    u /= npl.norm(u)
    u = np.reshape(u,(3,1))
    n = np.reshape(n,(3,1))
    v = np.cross(n,u,axis=0)
    v /= npl.norm(v)


    # Finding the parametric equation for the projected ellipse 
    k0 = u.T@M@u
    k1 = u.T@M@v
    k2 = v.T@M@v
    k3 = 2*(lam - np.dot(n.T,r_J))*u.T@M@n
    k4 = 2*(lam - np.dot(n.T,r_J))*v.T@M@n
    k5 = (lam - np.dot(n.T,r_J))**2*n.T@M@n

    # Conver to standard quadritc form for an ellipsoid 
    P = (np.array([[np.asscalar(k0),np.asscalar(k1)],
               [np.asscalar(k1),np.asscalar(k2)]])) 
    sign = 1 
    if np.all(npl.eig(P)[0]<0): 
        sign = -1 

    eva, eve = npl.eig(sign*P) 
    Rot = eve
    
    D = np.diag([eva[0],eva[1]])
    q = sign*np.array([np.asscalar(k3),np.asscalar(k4)])
    r = sign*np.asscalar(k5)


    beta = Rot@q 
    phi = np.abs((beta[0]**2/(4*D[0,0]) + beta[1]**2/(4*D[1,1]) - r)/(D[0,0]*D[1,1])) 

    center = np.array([ -beta[0]/(2*D[0,0]), -beta[1]/(2*D[1,1])])
    axis = np.array([np.sqrt(D[1,1]*phi),np.sqrt(D[0,0]*phi)])
    axis = np.sort(axis)[::-1]

    return axis,center

def planetmodel(nu,T_dab,p, R, imsize,planet,pixscale, Jansky = True, uniform=False, conv = 0.01): 
    # Planet is in pixel already 

    t0 = time.time()  # start time
    (x,y,) = axisforcasamodel(imsize, planet/pixscale)
    t1 = time.time() # end time
    print('1',t1 - t0,'s')  


    t0 = time.time()  # start time


    if uniform :     
        # Call the three dimension surface 
        t0 = time.time()  # start time
        #IPython.embed()
        (xv,yv,zv) = triaxialellipsoid(R/R[0]*planet/2/pixscale,x,y)
        t1 = time.time() # end time
        print('2',t1 - t0,'s') 
        # Avoid singularities 
        zv[zv==0] = float('nan')
        zv[~np.isfinite(zv)] = 0.
        T = np.ones_like(zv)*T_dab
        T[zv==0] = 0
        model = T
        # the code to time goes here
    
    else :
        print(planet/2/pixscale, planet, pixscale)
        model = brightnessmap(R,planet/2,x,y,T_dab,p,conv)

 

    if Jansky: 
        model = tb2jy(model, nu, pixscale)
    return model

def calcimsize(planet, pixscale):
    """Calculate the optimal image size for Casa 

    Casa works better if the imsize is a multiple of 2,3,5,7 
    All of them multiplied is 210. Function gives the closest multiple 
    to 210. Returns the imsize so that Jupiter is at least three times 
    in the frame. 
    
    Parameters
    -------
    planet : [1] int
        [arcsec] Size of the planet in pixels  
    pixscale : [1] float 
        [arc/pixel] Pixel scale, how many arcsec is one pixel 


    Returns
    ----------
    imsize : [1] int
        [-] Recommended size in pixels for Casa 

    Keywords
    ----------


    Example
    -------
    >>> calcimsize(210,0.15)

    References
    ------------

    Notes
    -------
    10/26/2017, CM, Initial Commit
    """

    base = 210 # (pixel) minimum size 
 
    PPlant = planet/pixscale # Planet in pixel 

    imsize = np.ceil(2.0*planet/pixscale/base)*base 

    return np.int(imsize)

def hms2deg(hms): 
        """Convert hh:mm:ss.ss to deg 
        
        Parameters
        -------
        self : [1] str
            [] 'hh mm ss.ss' 

        Returns
        ----------
        deg : [1] float
            [deg] Angle converted into degrees

        Keywords
        ----------


        Example
        -------
        >>> ra = ra.hms2deg()

        References
        ------------

        Notes
        -------
        11/18/2017, CM, Initial Commit
        """
        

        if len((hms.rsplit(' '))) > 1:
            split = (hms.rsplit(' '))
        elif len((hms.rsplit(':'))) > 1:
            split = (hms.rsplit(':'))
        elif len((hms.rsplit('.'))) > 1:
            split = (hms.rsplit('.'))
        else : 
            print('Did not recognize format')


        angle = ((float(split[0]) + 
                float(split[1])/60 + 
                float(split[2])/3600)*360/24)

        return angle

def dms2deg(dms): 
        """Convert degree mm:ss.ss to deg 
        
        Parameters
        -------
        self : [1] str
            [] 'hh mm ss.ss' 

        Returns
        ----------
        deg : [1] float
            [deg] Angle converted into degrees

        Keywords
        ----------


        Example
        -------
        >>> ra = ra.dms2deg()

        References
        ------------

        Notes
        -------
        11/18/2017, CM, Initial Commit
        """

        if len((dms.rsplit(' '))) > 1:
            split = (dms.rsplit(' '))
        elif len((dms.rsplit(':'))) > 1:
            split = (dms.rsplit(':'))
        elif len((dms.rsplit('.'))) > 1:
            split = (dms.rsplit('.'))
        else : 
            print('Did not recognize format')
        
       
        if (np.sign(float(split[0]))) == 0: 
            if split[0][0] == '-': 
                sign = -1.
            else: sign = +1. 
        else: sign =  np.sign(float(split[0])) 

        if len(split) == 4 : 
            angle = (sign*
                (np.abs(float(split[0])) 
                 + float(split[1])/60 
                 + (float(split[2])+float(split[3])*10**(-(len(split[3]))))/3600))  
        else: 
            angle = (sign*
                (np.abs(float(split[0])) 
                 + float(split[1])/60 
                 + (float(split[2]))/3600))  

        return angle

def deg2hms(angle): 
        """Convert angel in deg to hour format  
        Used for right ascension 
        
        Parameters
        -------
        deg : [1] float
            [deg] Angle converted into degrees

        Returns
        ----------
        self : [1] str
            [] 'hh mm ss.ss' 

        Keywords
        ----------


        Example
        -------
        >>> ra = ra.h2deg()

        References
        ------------

        Notes
        -------
        11/18/2017, CM, Initial Commit
        """
        h = np.fix(angle/360*24.) 
        m = np.fix((angle/360*24. - h )*60)
        s = np.fix((((angle/360*24. - h )*60) - m)*60)
        f = (((((angle/360*24. - h )*60) - m)*60) - s)*100

        string = '{:n}h {:n}m {:n}.{:n}s'.format(h,m,s,f) 

        return string

def deg2dms(angle): 
        '''Convert angel in deg to hour format . Used for right ascension
        
        Parameters
        -------
        deg : [1] float
            [deg] Angle converted into degrees

        Returns
        ----------
        self : [1] str
            [] 'hh mm ss.ss' 

        Keywords
        ----------


        Example
        -------
        >>> ra = ra.h2deg()

        References
        ------------

        Notes
        -------
        11/18/2017, CM, Initial Commit
        '''

        sign = np.sign(angle)
        angle = np.abs(angle) 
        d = np.fix(angle)
        dd = (angle - d)
        m = np.fix(dd*60)
        dm = dd*60 - m
        s = np.fix(dm*60)
        ds = dm*60 - s 
        f = round((ds*100000))


        string = '{:n}d {:n}m {:n}.{:n}s'.format(sign*d,m,s,f) 

        return string

def c_to_g(lat,R_e,R_p): 
    '''Convert planetocentric latitud eto planetographic latitude 
    
            Parameters
            -------
            lat : [1] float
                [rad] latitude 

            f : [1] float
                [-] (a-b)/a

            Returns
            ----------
            self : [1] str
                [] 'hh mm ss.ss' 

            Keywords
            ----------


            Example
            -------
            >>> 

            References
            ------------
            2016 - Implications of MAVEN’s planetographic coordinate
            system for comparisons to other recent Mars orbital
            missions

            Notes
            -------
            11/18/2017, CM, Initial Commit
    '''


    threshold = np.radians(89)
    lat[np.where(np.abs(lat)>=threshold)] = np.nan
    f = (R_e - R_p)/R_e
    return np.arctan(np.tan(lat)*(1-f)**2)
    
def zonalstructure(fitsfile, f_interp ,th_interp = np.array([]), residual = False ):
    '''Read out the zonal structure. 

    Current version reads out the residuals from Imke de Pater, 2014 
    VLA scans for Jupiter.  

    Parameters
    -------
    fitsfile : [1] string
        [] Name of the input fits file 
    f_interp : [1] float 
        [Hz] Required frequency  

    Returns
    ----------
    T : [1] K
        [deg] Disk averaged brightness temperature 
    p : [1] 
        [-] limb darkening coefficient 

    Keywords
    ----------
    th_interp: [Nx1] 
        [rad] Planetographic Latitude to be interpolated on 
    residual : [1] boolean
        [] Only return the residual structure after subtracting the ld disk 


    Example
    -------
    import pyPR.PlanetGeometry as PG  
    fitsfile = 'Data/VLA_ZonalScans_2014.fits'
    nu = 9e9 # [Hz] Frequency 
    T_res,th = zonalstructure(fitsfile, nu)

    References
    ------------
    2018, Imke de Pater, Jupiter’s ammonia distribution derived from VLA maps at 3–37 GHz

    Notes
    -------
    12/23/2018, CM, Initial Commit
    ''' 


    hdul = fits.open('pyPR/'+fitsfile)     
  
  
    # [deg] Emission angle  
    th_center = hdul[0].header['CRPIX1']
    dth       = hdul[0].header['CDELT1']    
    nth       = hdul[0].header['NAXIS2'] # Axis appear flipped 

    # Frequency      
    f_center  = hdul[0].header['CRPIX2'] 
    df        = hdul[0].header['CDELT2'] 
    nf        = hdul[0].header['NAXIS1'] 



    # Build the axis 
    th_l = th_center - (nth-1)/2*dth 
    th_u = th_center + (nth-1)/2*dth     
    th = np.linspace(th_l,th_u,nth,)
    
    f_l = f_center - (nf-1 )/2*df 
    f_u = f_center + (nf-1 )/2*df     
    f = np.linspace(f_l,f_u,nf,)

    # Read the temperature data 
    T_g = hdul[0].data # [Hz] frequency [30] 
    

    # Interpolation  
    T_map = scipy.interpolate.interp2d(f,th,T_g, kind='cubic') 

    # Interpolate the residuals onto the required frequency and angle grid 
    
    if th_interp.size == 0: 
        th_interp = np.radians(th)
    
    # Input is in degrees!
    T_interp = T_map(f_interp,np.degrees(th_interp)).flatten()

    # Adding the disk back in if required 
    if residual: 
        T,p = interpmodelparams(f_interp, printoutput = False)
        T_disk = T*np.cos((th_interp))**p    
        T_interp -=  T_disk 

    return T_interp, th_interp 


def interpmodelparams(f, planetname = 'jupiter', units = 'hz', printoutput = False): 
    '''Read out disk averaged brightness temp from Imke de Pater, 2016, 
    peering below the clouds 

    Parameters
    -------
    center : [1] float
        [Hz] center frequency or wavelength 

    Returns
    ----------
    T : [1] K
        [deg] Disk averaged brightness temperature 
    p : [1] 
        [-] limb darkening coefficient 

    Keywords
    ----------


    Example
    -------
    import pyPR.PlanetGeometry as PG  
    PG.interpmodelparams(10e9)

    References
    ------------

    Notes
    -------
    5/7/2018, CM, Initial Commit
    '''

    if planetname.casefold().strip() != 'jupiter': 
        sys.exit('Functionality only available for Jupiter')


    temperature = (np.array(
           [[0.06321, 144.92740],
            [0.06526, 148.91039],
            [0.06791, 152.53045],
            [0.07067, 156.15050],
            [0.07472, 159.40699],
            [0.07838, 162.66410],
            [0.08221, 164.47194],
            [0.08487, 167.00566],
            [0.08761, 168.81475],
            [0.09263, 170.98429],
            [0.09486, 171.70705],
            [0.09871, 171.34161],
            [0.10436, 173.87346],
            [0.11122, 177.12932],
            [0.12234, 179.29573],
            [0.13246, 179.65180],
            [0.14002, 178.92279],
            [0.14340, 180.00787],
            [0.14687, 181.09295],
            [0.15404, 182.17616],
            [0.16155, 181.81009],
            [0.17078, 181.08108],
            [0.18054, 180.71439],
            [0.18936, 181.79760],
            [0.19703, 182.88143],
            [0.20666, 183.96463],
            [0.22552, 184.32008],
            [0.23280, 184.31758],
            [0.24031, 184.31508],
            [0.25607, 183.94776],
            [0.27503, 183.57982],
            [0.29306, 183.03135],
            [0.30980, 182.12118],
            [0.33010, 181.02923],
            [0.35173, 179.57496],
            [0.36887, 178.48426],
            [0.39305, 177.39231],
            [0.41879, 175.57573],
            [0.43920, 173.76039],
            [0.45697, 173.03263],
            [0.47545, 171.94255],
            [0.50659, 169.76365],
            [0.53977, 167.58475],
            [0.57969, 165.04290],
            [0.62752, 162.13811],
            [0.66332, 159.59752],
            [0.70675, 157.05629],
            [0.76507, 154.15151],
            [0.86853, 149.06906],
            [0.95517, 144.35144],
            [1.07577, 139.26962],
            [1.16445, 134.55324],
            [1.21148, 132.01389],
            [1.25051, 130.56212],
            [1.32203, 131.64470],
            [1.36481, 133.81611],
            [1.43155, 136.34859],
            [1.51356, 139.60507],
            [1.58757, 142.13755],
            [1.67845, 144.30708],
            [1.76047, 146.11492],
            [1.89115, 149.73248],
            [2.03153, 153.35003],
            [2.14792, 156.60652],
            [2.32585, 161.31040],
            [2.51852, 166.01429],
            [2.83780, 172.52664],
            [3.14727, 179.76487],
            [3.63227, 189.89852],
            [4.29336, 201.84189],
            [4.95525, 213.42481],
            [5.95141, 227.54083],
            [7.09139, 242.01980],
            [8.38255, 255.41243],
            [10.0684, 271.34004],
            [12.5845, 290.16307],
            [15.1161, 307.17764],
            [18.6000, 331.43668],
            [22.1692, 353.16200],
            [25.5936, 371.26664],
            [28.6192, 385.75061],
            [30.2636, 392.99259],] ))

    interpT = scipy.interpolate.interp1d(cst.c.value/(temperature[:,0]/100),temperature[:,1])

    # limb darkening 
    # GHz, limb darkening coefficient 
    limbd = (np.array(
            [[4.52  ,0.16 ], 
             [5.49  ,0.16 ],
             [6.5   ,0.16 ],
             [7.5   ,0.16 ], 
             [8.5   ,0.16 ],
             [9.52  ,0.16 ],
             [10.46 ,0.16 ],
             [11.46 ,0.16 ],
             [13.18 ,0.08 ],
             [14.21 ,0.08 ],
             [15.18 ,0.08 ],
             [16.21 ,0.08 ],
             [17.38 ,0.06 ],
             [25.00 ,0.06 ],]))


    intperld = scipy.interpolate.interp1d(limbd[:,0]*1e9,limbd[:,1])


    # Convert to units of cm 
    if units.casefold().strip() == 'hz': 
        frequency = f 
    elif units.casefold().strip() == 'm': 
        frequency = cst.c.value/f 
    elif units.casefold().strip() == 'cm':
        frequency = cst.c.value/(f/100) 
    else: 
        sys.exit('Units not recognized.')

    # Disk averaged brightness temperature 
    # Limb darking coefficient 
    T = interpT(frequency)
    p = intperld(frequency)

    if printoutput:
        print('The disk averaged brightness temperature is: {:3.2f}'.format(float(T)))
        print('The limb darkening coeffient is: {:3.3f}'.format(float(p)))


    return np.array([T,p])




# Class definitions 

def mwe(planetname = 'Jupiter',tobs='2017-02-02 06:32:49',nu=22e9,T = 132.7,
        p = np.array([0.075,0.065]), beamsize = 0.7): 
    Pl = Planet(planetname)  # Initiate the Planet
    Pl.ephemeris(tobs)      # Querry the Planet's ephemeris
    Pl.initmodel('M')   # Initiate the model 
    Pl.M.gen_casa(nu,T,p, beamsize = beamsize)  # Generate the model 
    Pl.M.plot(Pl.M.data) # Plot the model

class Planet:

    def __init__(self,target): 
        self.name = target 

    def ephemeris(self, tstart, tend=None, nstep=1, obs_code = '-5'):
        '''Immediately run get_ephemerides, then set a bunch of 
        class variables corresponding to different information found in
        the ephemeris.
        '''

        # Be default increment time by one step to provide single epehemeris

        if tend is None:
            try:
                tend_obj = datetime.strptime(tstart,'%Y-%m-%d %H:%M:%S.%f')
                tend_obj = tend_obj + timedelta(0,1)
                tend = datetime.strftime(tend_obj,'%Y-%m-%d %H:%M:%S.%f')
            except ValueError: 
                try: 
                    tend_obj = datetime.strptime(tstart,'%Y-%m-%d %H:%M:%S')
                    tend_obj = tend_obj + timedelta(0,1)
                    tend = datetime.strftime(tend_obj,'%Y-%m-%d %H:%M:%S')
                except ValueError: 
                    tend_obj = datetime.strptime(tstart,'%Y-%m-%d %H:%M')
                    tend_obj = tend_obj + timedelta(0,0,0,0,1)
                    tend = datetime.strftime(tend_obj,'%Y-%m-%d %H:%M')
               
                
        



        # Determine, the number of steps 
        intv = np.linspace(0,nstep-1,nstep,dtype = int).tolist()

      
        self.target = self.name
        self.obs_code = obs_code
        self.ephem, self.observatory_coords, self.radius = (
            get_ephemerides(naif_lookup(self.target),
            tstart,tend,str(nstep),self.obs_code))
        self.ellipticity = np.sqrt((self.radius[0]**2-self.radius[1]**2)/self.radius[0]**2)
        self.flattening = 1. - np.sqrt(1. - self.ellipticity**2)
        self.time = self.ephem[intv,0] # UTdate UTtime
        self.sun = [s.strip() for s in self.ephem[intv,1]]
        self.moon = [s.strip() for s in self.ephem[intv,2]]
        self.ra = hms2deg(self.ephem[intv,3][0]) # RA (J2000) (hh:mm:ss.ff) converted to degree
        self.dec = dms2deg(self.ephem[intv,4][0]) # DEC (J2000) (hh:mm:ss.ff) converted to degree
        # (arcsec^2/h) (cosine of the declination) arcsec hr-1 
        self.dra = np.asarray([float(s) for s in self.ephem[intv,5]]) 
        #  d(DEC)/dt 
        self.ddec = np.asarray([float(s) for s in self.ephem[intv,6]]) 
        try: 
            # (degrees) azimuth , North = 0 = 360
            self.azimuth = np.asarray([float(s) for s in self.ephem[intv,7]])    
        except ValueError: 
            self.azimuth = np.asarray([np.nan for s in self.ephem[intv,7]])
        try: 
            # (degrees) elevation degrees, above horizon   
            self.elevation = np.asarray([float(s) for s in self.ephem[intv,8]])
        except ValueError: 
            self.elevation = np.asarray([np.nan for s in self.ephem[intv,8]])
        try:
            # Airmass 
            self.airmass = read_ephem_line(self.ephem[intv,9]) 
        except ValueError: 
            self.airmass = np.nan
        try:    
            # Extinction
            self.extinction = read_ephem_line(self.ephem[intv,10]) 
        except ValueError: 
            # Extinction
            self.extinction = np.nan 

        # Visual magntidue 
        self.apmag = read_ephem_line(self.ephem[intv,11]) 
        # Surface brightness
        self.sbrt = read_ephem_line(self.ephem[intv,12]) 
        # (arcsec) Satellite angular separ/vis. (arcsec) The angle between the
        # center of target object and the center of the primary body it 
        # revolves around,
        self.ang_sep = read_ephem_line(self.ephem[intv,13]) 
        # Visibility 
        self.visibility = [s.strip(' ') for s in self.ephem[intv,14]]
        # Target angular diameter (arcsec)
        self.ang_diam = read_ephem_line(self.ephem[intv,15]) #arcsec
        # Planet orientation information 
        # (degree) Observer sub-longitude, planetodetic, positive to west
        self.ob_lon = read_ephem_line(self.ephem[intv,16])  
        # (degree) Observer sub-latitude, planetodetic,
        self.ob_lat = read_ephem_line(self.ephem[intv,17]) 
        # (degree) sun sub-longitude, planetodetic, positive to west
        self.sun_lon = read_ephem_line(self.ephem[intv,18])  
        # (degree) sun sub-latitude, planetodetic,
        self.sun_lat = read_ephem_line(self.ephem[intv,19]) #degrees
        # (degree) Sub-Sun position angle 
        self.ssun_dis = read_ephem_line(self.ephem[intv,20])  
        # (arcsec) Sub-Sun distance 
        self.ssun_ang = read_ephem_line(self.ephem[intv,21]) 
        # (degrees) North pole position angle 
        self.np_ang = read_ephem_line(self.ephem[intv,22])  
        # (arcsec) North pole distance(arcsec) Distance from the sub 
        # observer point to the north pole 
        self.np_dis = read_ephem_line(self.ephem[intv,23]) 
        # (AU) Heliocentric range (AU)
        self.hrange = read_ephem_line(self.ephem[intv,24]) 
        # (km/s) Heliocentric range-rate (km/s)
        self.hrrate = read_ephem_line(self.ephem[intv,25])
        # (AU) Observer range (AU)
        self.orange = read_ephem_line(self.ephem[intv,26]) 
        # (km/s) Observer range-rate (km/s)
        self.hrrate = read_ephem_line(self.ephem[intv,27])
        # (min) One-way (down-leg) light-time
        self.lt = read_ephem_line(self.ephem[intv,28])
        # (degree) Sun-Target-Observer ~PHASE angle 
        self.ph_ang = read_ephem_line(self.ephem[intv,29])
        # (degree) North pole right asencsion 
        self.np_ra = read_ephem_line(self.ephem[intv,30])
        # (degree) North pole declination 
        self.np_dec = read_ephem_line(self.ephem[intv,31])




        # Equatorial radius for triaxial ellipsoid 
        self.principalaxis = np.copy(self.radius)
        self.principalaxis[1:3] = (rotateprincipalaxis_3D(
                            self.radius, np.radians(self.ob_lat), 
                            np.radians(self.ob_lon), self.orange*cst.au.value))[0] 
        


        # [sterradians] solid angle extended by the disc as seen by earth  
        # (http://starlink.eao.hawaii.edu/docs/sun213.htx/sun213se9.html)
        self.sa = np.pi*(np.sqrt(self.principalaxis[1]*self.principalaxis[2]*1e6)/self.orange/cst.au.value)**2



    def initmodel(self,modelname):
        try:
            self.radius
            setattr(Planet,modelname,Model(self.name, planet = self))
            print('The model can be accessed under: '\
             '{:s}.{:s}.'.format(self.name,modelname))
        except AttributeError:
            print('First querry the Planet\'s properties with epemeris'\
                  ' or define the radius: <Planet>.radius = np.array([R,R,R])')

    def model(self):
        print('Accessing Model')
        self.Model(self.name)


    
class Model: 
    """Class containing the required information to create, modify and 
        export a model.
        
        Parameters
        -------
        nu: [1] float
            [Hz] Frequency of observations 
        T : [1] float
            [K] Disk averaged brightness temperature

        Returns
        ----------
        deg : [1] float
            [deg] Angle converted into degrees

        Keywords
        ----------
        planet querries information required for the model 

        Example
        -------


        References
        ------------

        Notes
        -------
        11/18/2017, CM, Initial Commit
        """

    def __init__(self, name, planet = None):
        self.name = name 
        if planet is None: 
            print('Manual mode: Execute gen_casa_input to set input parameter ') 
            print('For automatic mode: Initialize with <Planet>.initmodel')
        else: 
            # For the model 
            self.ang_diam = planet.ang_diam 
            self.radius = planet.radius
            self.ob_lat = planet.ob_lat
            self.ob_lon = planet.ob_lon
            self.orange = planet.orange
            self.np_ang = planet.np_ang 
            # For the export 
            self.ra = planet.ra
            self.dec = planet.dec
            self.orange = planet.orange
            self.time = planet.time[0]

    def gen_casa_input(self, ang_diam, radius, ob_lat, np_ang,beamsize):
        ''' Setting the input parameters for gen_casa manually ''' 

        self.ang_diam = ang_diam 
        self.radius   = radius
        self.ob_lat   = ob_lat
        self.np_ang   = np_ang 

    
    def gen_casa(self, nu, T, p, beamsize = 0.7, psfsampling=5, 
                 Jansky = True, setimsize = False, conv = 0.01, 
                 T_res=np.zeros(100),
                 th_res=np.arange(-np.pi/2,np.pi/2,np.pi/100), plotting=False ):
        """Generate a model for casa 
        
        Parameters
        ---------- 

        nu : [1] float
            [Hz] Frequency 
        T : [1] float
            [K] Target disk averaged temperature 
        p : [1x2] float 
            [EW,NS] Limb darkening coefficient 
        ob_lat_d: [1] 
            [rad] Sub observer latitude 
        orange: [1] 
            [m] Observer range 
        
        th_res : [1xN] float
            [rad] Corresponding latitudes for the residual temperatures 
        


        Returns
        ----------


        Keywords
        ----------
        beamsize : [1] float 
            [arcsecs] Beam size to be used 
        psfsamplin : [1] float
            [-] Number of pixels per beam

        Jansky : [1] boolean
            [-] Output in Jy/pixel instead of brightness temperature
   
        conv : [1] float
            [-] Convergence criteria for iteration to match disk averaged temperature 
                High conv will create disk with the center temperature 
        T_res: [Nx1] float 
            [K] Zonal temperature profile for creating a zonally average model 
        th_res: [Nx1] float 
            [rad] Corresponding planetographic latitude 



        Example
        -------


        References
        ------------

        Notes
        -------
        11/18/2017, CM, Initial Commit
        12/15/2018, CM, added zonal averaged model capabilities 
        """

        try:
            self.ang_diam 
            self.radius   
            self.ob_lat   
            self.np_ang   
        except AttributeError:
            sys.exit('Attributes for gen_casa are missing.'\
                   'Execute gen_casa_input or initialize a planet')


        if np.size(p) == 1:
            p = np.array([p,p])

        # Assign a model type for documentation 
        self.modeltype = 'pyPR-CASA'
        self.Jansky = Jansky 
        self.pixscale = beamsize/psfsampling 

        # Assign the parameter 
        self.obsfrequency = nu 
        self.T_da = T 
        self.limbdarkening = p


        # Radius of the planet in pixel 
        self.r_pla = self.ang_diam/self.pixscale/2
        self.planetsize = 2*self.r_pla

        # Normalize the axis for the triaxial ellipsoid and convert to pixel 


        # Rotate around x axis to correct for the sub latitude 
        # (pixel) [equator, polar, equator]

        self.principalaxis = np.copy(self.radius)
        # print(np.radians(self.ob_lat),np.radians(self.ob_lon),self.radius)
        self.principalaxis[1:3] = (rotateprincipalaxis_3D(self.radius, 
                                                     np.radians(self.ob_lat), 
                                                     np.radians(self.ob_lon), 
                                                     self.orange*cst.au.value))[0]
        
        R_temp  = self.principalaxis/self.principalaxis[0]*self.r_pla

        if setimsize: 
            self.imsize = setimsize
        else: 
            self.imsize = calcimsize(self.ang_diam,self.pixscale) # (pixel)
        print('Use/check the following parameters for your casa deconvolution:')
        print('Imsize: ', self.imsize) 
        print('Cell : ', self.pixscale)
      
        rotangle = -(self.np_ang)
        # I can probably short cut this ToDo 
        # model = (planetmodel(self.obsfrequency, self.T_da, 
        #         self.limbdarkening, self.principalaxis, self.imsize, self.planetsize, 
        #         self.pixscale, Jansky, conv = conv))
        # self.data = scipy.ndimage.rotate(model,rotangle,order=0,reshape = False)

        # Create axis of the model  
         # self.ang_diam/self.pixscale is 
        rotangle = -(self.np_ang)
        (x,y,) = axisforcasamodel(self.imsize, self.planetsize/self.pixscale)
        ld_model = (brightnessmap(self.principalaxis, self.r_pla, x, y, 
                        self.T_da, self.limbdarkening, conv = 0.01))
        self.Tdata = scipy.ndimage.rotate(ld_model,rotangle,order=0,reshape = False)
        self.Bdata = scipy.ndimage.rotate(tb2jy(ld_model, nu, self.pixscale),rotangle,order=0,reshape = False)
        # To agree with previous version 
        self.data = self.Bdata

        if np.all(T_res != np.zeros(100)): 
            zonal_model = (zonal_residual_model(self.principalaxis, self.r_pla, 
                                             x, y, 
                                             self.limbdarkening, 
                                             np.radians(self.ob_lat), self.orange*cst.au.value, 
                                             T_res, th_res, 
                                             conv = 0.01, plotting = plotting))
            self.Tzonaldata = scipy.ndimage.rotate(zonal_model,rotangle,order=0,reshape = False)
            self.Bzonaldata = scipy.ndimage.rotate(tb2jy(zonal_model, nu, self.pixscale),rotangle,order=0,reshape = False)
    



    def gen_gaussian_structure(self,n_a,sig_xi,sig_yi,scalefactor, spacingfactor = 1, inc_negative=True): 
        """Generate random structure on Jupiter 

        This function will place Gaussian structure on Jupiter, where 
        the sigma parameters determine the shape the Gaussian structure
        and the scale factor determines the magnitude of the signal. The 
        spaceing factor determines how tightly spaced the individual 
        signals are. For tight spacing you get an overlap of signals. 

        
        Parameters
        -------
        n_a: [1] float
            [-] Number of annuli that contain structure
        sig_xi [1] int
            [-] Size of the enhancement structure (x-axis) 
        sig_yi [1] int
            [-] Size of the enhancement structure (y-axis) 
        scalefactor [1] 
            [-] Magnitude of the enhancement. It can be defined with 
                respect to the background structure. Such as 10% of the 
                maxium on the disk  

        Returns
        ----------
        deg : [NxN] float
            [-] Structure on top of the background 

        Keywords
        ----------
        spacingfactor : [1] float
            [-] Modify how tightly spaced the enhancement is. Factor of 
                1 very tight, 100 very widely spaced. 


        Example
        -------
        beam = 0.08
        sig_x = 1
        sig_y = 3

        Jupiter.initmodel('Model') 
        nu = 22e9; T = 132.7; p = np.array([0.08,0.08]); 
        Jupiter.Model.gen_casa(nu,T,p,beamsize = beam,)
        scalefactor = 0.1*np.max(Jupiter.Model.data)
        Jupiter.Model.gen_structure(n_a = 7, sig_xi=sig_x, sig_yi=sig_y, scalefactor = scalefactor)

        References
        ------------

        Notes
        -------
        6/6/2018, CM, Added negative residuals 
        7/6/2018, CM, Added negative residuals 
        """
        from scipy.stats import multivariate_normal
        from scipy.interpolate import interp2d
        # Remove the tilt of the planet, and project the planet on the 
        # sky with y axis pointing up
        rotangle = (self.np_ang)
        self.data = scipy.ndimage.rotate(self.data,rotangle,order=0,reshape = False)

        # Find the planet in pixels 
        r_pix = self.planetsize/2 # Radius in pixel 

        # n_a = 5 # [-] Number of annuli on the planet 
        for i in range(1,n_a):
            # Compute the local radius in the spiral 
            r_l = r_pix/n_a*(n_a-i)

            # Build a mesh based on the expected size outwards to 3 sigma 
            sig_x, sig_y = sig_xi+r_pix/100/n_a*(i-1)**2, sig_yi + r_pix/100/n_a*(i-1)**2 # size of enhancement in pixel 
            print(sig_x,sig_y)
            # 
            var = multivariate_normal(mean=[0,0], cov=[[sig_x**2,0],[0,sig_y**2]])


            # Decrease the scalefactor to avoid heavy concentration in the middle of the planet 
            # scalefactor = scalefactor

            # Calculate the number of steps through the unit circle 
            dth_m = (spacingfactor*sig_x)/(2*np.pi*r_l)

            th = np.zeros(1) #(i-1)*np.pi/4
            i2 = 0 
            # step through the unit circle 
            while th[-1] < 2*np.pi: #+(i-1)*np.pi/4: 
                i2 += 1
                # Convert the base center to x and y 
                th = np.append(th,th[-1]+dth_m*i2)
                x_l = r_l*np.cos(th[-1]) 
                y_l = r_l*np.sin(th[-1]) 
                
                
                # Mesh to interpolate on (plus on upper end to make it odd)
                x_m, y_m = np.mgrid[-3*sig_x:4*sig_x:1, -3*sig_y:4*sig_y:1]
                pos = np.empty(x_m.shape + (2,))
                pos[:, :, 0] = x_m; pos[:, :, 1] = y_m 

                # Mesh were to add the enhancement (enhanecment grid, plus planet grid, plus image grid) 
                x_i = x_m + np.round(x_l) + self.imsize/2
                y_i = y_m + np.round(y_l) + self.imsize/2
                # IPython.embed()
                # Lower bounds for slicing 
                x_ilb = x_i[0,0].astype(int)  
                y_ilb = y_i[0,0].astype(int)
                if inc_negative:
                    sign = [-1,1][np.random.randint(-1,2)]
                else: 
                    sign = 1

                self.data[x_ilb:x_ilb+x_m.shape[0],y_ilb:y_ilb+x_m.shape[1]] +=  \
                var.pdf(pos)/np.max(var.pdf(pos))*scalefactor*sign  



        # Add the tilt of the planet, and project the planet on the sky with y axis pointing up
        rotangle = -(self.np_ang)
        self.data = scipy.ndimage.rotate(self.data,rotangle,order=0,reshape = False) 


    def gen_banded_structure(self,scalefactor,width = 1, spacingfactor = 1,inc_negative=True): 
        """Generate random structure on Jupiter 

        This function will place banded structure on the planet, with 
        the size of the bands increasing away from the equator. The 
        various parameter allow to control the width, the spacingfactor 
        and the brightness.  

        
        Parameters
        -------  
        scalefactor [1] 
            [-] Magnitude of the enhancement. It can be defined with 
                respect to the background structure. Such as 10% of the 
                maxium on the disk  

        Returns
        ----------
        deg : [NxN] float
            [-] Structure on top of the background 

        Keywords
        ----------
        spacingfactor : [1] float
            [-] Modify how tightly spaced the enhancement is. Factor of 
                1 very tight, 100 very widely spaced.

        width: [1]
            [-] Width of the bands in pixel  
        inc_negative: [Bol]
            [-] Include also negative residuals 


        Example
        -------
        beam = 0.08
        spacingfactor = 2
        width = 1 

        Jupiter.initmodel('Model') 
        nu = 22e9; T = 132.7; p = np.array([0.08,0.08]); 
        Jupiter.Model.gen_casa(nu,T,p,beamsize = beam,)
        scalefactor = 0.1*np.max(Jupiter.Model.data)
        (J.M.gen_banded_structure(width = width, 
        scalefactor = scalefactor,spacingfactor = spacingfactor))

        References
        ------------

        Notes
        -------
        7/6/2018, CM, Initial Commit
        """

        # Remove the tilt of the planet, and project the planet on the 
        # sky with y axis pointing up
        rotangle = (self.np_ang)
        self.data = scipy.ndimage.rotate(self.data,rotangle,order=0,reshape = False)

        # Find the planet in pixels 
        r_pix = self.planetsize/2 # Radius in pixel 

        # Step through the negative y axis
        i = int(self.imsize/2) - int(r_pix) 


        # Minimum flux in the picture 
        f_min = np.min(self.data[np.where(self.data>0)])

        while i < int(self.imsize/2) + int(r_pix): 

            if inc_negative:
                sign = [-1,1][np.random.randint(-1,2)]
            else: 
                sign = 1

            i_range = int(width**(5/3))
            print(i_range)
            c_pix = self.data[i+i_range,int(self.imsize/2)]

            if c_pix > f_min: 
                for j in range(i,i+i_range):
                    self.data[j,np.where(self.data[j,:] > 
                        f_min)] += scalefactor*sign

            i += spacingfactor*i_range
            width += 1

        # Add the tilt of the planet, and project the planet on the sky with y axis pointing up
        rotangle = -(self.np_ang)
        self.data = scipy.ndimage.rotate(self.data,rotangle,order=0,reshape = False) 

    def combine_structure(self, background, bands, gaussians, bands_regions=[0,0,1,1], gaussians_regions=[1,0,0,0] ): 
        """Combine banded and gaussian structure on Jupiter 

        This function combines the banded and the gaussian structures. 
        The regions keyword allows to determine where the structures 
        are supposed to be placed. Overlapping structures are allowed.  
        
        Parameters
        -------
        background: [NxN] float
            [-] Limb-darkened model 
        bands [NxN] float
            [-] Banded structure without the background 
        gaussians [NxN] float
            [-] Gaussian structure without the background 
        

        Returns
        ----------
        deg : [NxN] float
            [-] Structure on top of the background 

        Keywords
        ----------
        bands_regions : [1x4] float
            [-] Regions where the enhancement can be found 
                Distributed in quartiles corresponding to the four 
                regions  
                
                Regions: 
                  2  1
                  3  4

        Example
        -------
        J.M.combine_structure(background,bands,gaussians,bands_regions=[0,0,1,1], gaussians_regions=[1,0,0,1])

        References
        ------------

        Notes
        -------
        7/6/2018, CM, Initial Commit
        """

        # Remove the tilt of the planet, and project the planet on the 
        # sky with y axis pointing up
        rotangle = (self.np_ang)
        background = scipy.ndimage.rotate(background,rotangle,order=0,reshape = False)
        bands = scipy.ndimage.rotate(bands,rotangle,order=0,reshape = False)
        gaussians = scipy.ndimage.rotate(gaussians,rotangle,order=0,reshape = False)

        # Find the planet in pixels 
        r_pix = self.planetsize/2 # Radius in pixel 
        model = background 

        if bands_regions[0]: 
            model[int(self.imsize/2):,int(self.imsize/2):] += bands[int(self.imsize/2):,int(self.imsize/2):]
        if bands_regions[1]: 
            model[int(self.imsize/2):,0:int(self.imsize/2)] += bands[int(self.imsize/2):,0:int(self.imsize/2)]
        if bands_regions[2]: 
            model[0:int(self.imsize/2),0:int(self.imsize/2)] += bands[0:int(self.imsize/2),0:int(self.imsize/2)]
        if bands_regions[3]: 
            model[0:int(self.imsize/2),int(self.imsize/2):] += bands[0:int(self.imsize/2),int(self.imsize/2):] 

        if gaussians_regions[0]: 
            model[int(self.imsize/2):,int(self.imsize/2):] += gaussians[int(self.imsize/2):,int(self.imsize/2):]
        if gaussians_regions[1]: 
            model[int(self.imsize/2):,0:int(self.imsize/2)] += gaussians[int(self.imsize/2):,0:int(self.imsize/2)]
        if gaussians_regions[2]: 
            model[0:int(self.imsize/2),0:int(self.imsize/2)] += gaussians[0:int(self.imsize/2),0:int(self.imsize/2)]
        if gaussians_regions[3]: 
            model[0:int(self.imsize/2),int(self.imsize/2):] += gaussians[0:int(self.imsize/2),int(self.imsize/2):]


        # model[0:int(imsize/2),0:int(imsize/2)] += gaussians[0:int(imsize/2),0:int(imsize/2)]
        # # Adding the banded structure 
        # model[:,int(imsize/2):-1] += bands[:,int(imsize/2):-1] 


        self.data = model 
        rotangle = -(self.np_ang)
        self.data = scipy.ndimage.rotate(self.data,rotangle,order=0,reshape = False) 



    def gen_uniform(self, nu, T, beamsize = 0.7, psfsampling=5, 
                 Jansky = True, setimsize = False, ):


        """Generate a uniform disk model for casa 
        
        Parameters
        -------
        beamsize: [1] float
            [arcs] Used to estimate the size of the image and the pixel
                   to arcsecond converstion

        Returns
        ----------
        deg : [1] float
            [deg] Angle converted into degrees

        Keywords
        ----------
        psfsamplin : [1] float
            [-] Number of pixels per beam
        Jansky : [1] float
            [K] Disk averaged brightness temperature
        T : [1] float
            [K] Disk averaged brightness temperature


        Example
        -------


        References
        ------------

        Notes
        -------
        11/18/2017, CMimpo, Initial Commit
        """

        try:
            self.ang_diam 
            self.radius   
            self.ob_lat   
            self.np_ang   
        except AttributeError:
            sys.exit('Attributes for gen_casa are missing.'\
                   'Execute gen_casa_input or initialize a planet')


        # Assign a model type for documentation 
        self.modeltype = 'pyPR-CASA'
        self.Jansky = Jansky 
        self.pixscale = beamsize/psfsampling 

        # Assign the parameter 
        self.obsfrequency = nu 
        self.T_da = T 



        # Radius of the planet in pixel 
        r_pla = self.ang_diam/self.pixscale/2
        self.planetsize = 2*r_pla

        # Normalize the axis for the triaxial ellipsoid and convert to pixel
        R = self.radius/self.radius[0]*r_pla
        # Rotate around x axis to correct for the sub latitude 
        # (pixel) Polar radius for triaxial ellipsoid
        R[1] = (R[1]*np.cos(np.radians(self.ob_lat))**2 
                + R[2]*np.sin(np.radians(self.ob_lat))**2 )
        # (pixel) Equatorial radius for triaxial ellipsoid 
        R[2] = (R[2]*np.cos(np.radians(self.ob_lat))**2 
                + R[1]*np.sin(np.radians(self.ob_lat))**2)
        
        if setimsize: 
            self.imsize = setimsize
        else: 
            self.imsize = calcimsize(self.ang_diam,self.pixscale) # (pixel)
        print('Use/check the following parameters for your casa deconvolution:')
        print('Imsize: ', self.imsize) 
        print('Cell : ', self.pixscale)


        self.limbdarkening = 0 

        model = (planetmodel(self.obsfrequency, self.T_da, 
                self.limbdarkening, R, self.imsize, self.planetsize, 
                self.pixscale, Jansky,uniform = True))
 
        # First iteration rotation. Needs automization 
        rotangle = -(self.np_ang)
        self.data = scipy.ndimage.rotate(model,rotangle,order=0,reshape = False)

        # Store pixelscale correctly for brightness model in brightness temperature 
        self.pixscale = self.ang_diam/self.planetsize




    def gen_general(self,nu,T,p,radius,imsize,planetsize,rotangle=0,):
        # rotangle = self.np_ang
        self.modeltype = 'pyPR-general'
        self.Jansky = False 
        self.imsize = imsize
        self.planetsize = planetsize
        self.pixscale = 0 # Not possible without distance information 
        # Radius of the planet in pixel 
        # Normalize the axis for the triaxial ellipsoid and convert to pixel
        R = radius/radius[0]*planetsize
        # Rotate around x axis to correct for the sub latitude  
        Jansky = False # Jansky requires a beam size        
        self.data = (planetmodel(self.obsfrequency, self.T_da, 
                self.limbdarkening, R, self.imsize, self.planetsize, 
                1., Jansky=False))
        # First iteration rotation. Needs automization 
        self.data = (scipy.ndimage.rotate(self.model,
                rotangle,order=0,reshape = False))

    def export(self, importfitsname, exportname, exportdata = None,units = 'Jy/pixel'): 
        # Import header information from a CASA image 
        print('When changing the size of image, make sure to load a new header.')

        fname = importfitsname +'.fits'
        outfile = exportname +'.fits'
        print(fname)
        # # Export 

        im = fits.open(fname,ignore_missing_end=True)
        hdu_out = im
        if exportdata is None: 
            hdu_out[0].data = self.Bdata
        else: 
            hdu_out[0].data = exportdata

        # hdu_out[0].header['BUNIT'] = units 
        # hdu_out[0].header['BMIN'] = np.float(self.pixscale)/3600 #set beam size equal to one pixel so uvsub doesnt get confused
        # hdu_out[0].header['BMAJ'] = np.float(self.pixscale)/3600


        hdu_out[0].header['NAXIS'] = 4
        hdu_out[0].header['NAXIS1'] = hdu_out[0].data.shape[0] 
        hdu_out[0].header['NAXIS2'] = hdu_out[0].data.shape[1]
        hdu_out[0].header['NAXIS3'] = 1
        hdu_out[0].header['NAXIS4'] = 1
        hdu_out[0].header['BUNIT'] = units 
        hdu_out[0].header['BMIN'] = np.float(self.pixscale)/3600 #set beam size equal to one pixel so uvsub doesnt get confused
        hdu_out[0].header['BMAJ'] = np.float(self.pixscale)/3600
        hdu_out[0].header['CRPIX1'] =hdu_out[0].data.shape[0]/2
        hdu_out[0].header['CDELT1']  =  -1*np.sign(hdu_out[0].header['CRVAL1'])*self.pixscale/360  
        hdu_out[0].header['CRPIX2'] =hdu_out[0].data.shape[1]/2
        hdu_out[0].header['CDELT2']  =  -1*np.sign(hdu_out[0].header['CRVAL2'])*self.pixscale/3600  


        hdu_out[0].writeto(outfile, overwrite=True)
        print('Model written to ', outfile, '\n\n')



    def exportasfits_input(self,ra,dec,time ):
        ''' Setting the input parameters for exportasfits manually ''' 

        self.ra  = ra
        self.dec = dec

    def exportasfits(self,data, exportname = 'output', units = 'Jy/pixel',ephemeris = True, header = True): 
        """ Import header infromation from CASA data 
    
        Extended description of the function.
        
        Parameters
        ----------
        arg1 : int
            [Unit] Description of arg1
        pixscale : int
            [arcseconds] Size of one pixel in arcseconds 
        
        Keyword Arguments
        ----------
        pt1: type
           [Unit] Description, Default: true
        
        Returns
        -------
        out1 :
            [Unit] Description of return value
        out2 :
            [Unit] Description of return value
        
        Warnings
        -------
        
        
        Example
        -------
        text
        
        >>>function(arg1, arg2)
        out1
        
        References
        ------------
        2-16 - Geiser - Representations of spectral coordinates in FITS
        https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html
        
        Todo
        ----- 
        
        Notes
        -------
        mm/dd/yy, Initials of Author, Short description of update
        """
    
        if ephemeris:
            try:
                self.ra
                self.dec 
                self.pixscale

                [dra,ddec] = (ephem_uncertainties(
                              self.name,self.orange*cst.au.value)[0:2])
            except AttributeError:
                sys.exit('Missing information for the export',\
                        ' First run exportasfits_input')

        now = datetime.now()

        hdu = fits.PrimaryHDU(data)
        hdulist = fits.HDUList([hdu])


        if header:
            hdulist[0].header['SIMPLE']  =  True #Standard FITS
            hdulist[0].header['BITPIX']  =  -32 #Floating point (32 bit)
            hdulist[0].header['NAXIS']   =   4                                                  
            hdulist[0].header['NAXIS1']  =   self.imsize                                                  
            hdulist[0].header['NAXIS2']  =   self.imsize                                                   
            hdulist[0].header['NAXIS3']  =                    1                                                  
            hdulist[0].header['NAXIS4']  =                    1                                                  
            hdulist[0].header['EXTEND']  =   True   #?                                                
            hdulist[0].header['BSCALE']  =   1.000000000000E+00 #PHYSICAL = PIXEL*BSCALE + BZERO                 
            hdulist[0].header['BZERO']   =   0.000000000000E+00                                              
            hdulist[0].header['BMAJ']    =   self.pixscale/3600                                                  
            hdulist[0].header['BMIN']    =   self.pixscale/3600                                                   
            hdulist[0].header['BPA']     =   45. # Position angle 
            hdulist[0].header['BTYPE']   = 'Intensity'                                                           
            hdulist[0].header['OBJECT']  = self.name.upper()                                                                                                                            
            hdulist[0].header['BUNIT']   = units #Brightness (pixel) unit                     
            hdulist[0].header['RADESYS'] = 'ICRS    ' 
            if ephemeris:                                                            
                hdulist[0].header['LONPOLE'] =   1.800000000000E+02                                                  
                hdulist[0].header['LATPOLE'] =   self.dec


            # hdulist[0].header['PC1_1']   =   1.000000000000E+00                                                  
            # hdulist[0].header['PC2_1']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC3_1']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC4_1']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC1_2']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC2_2']   =   1.000000000000E+00                                                  
            # hdulist[0].header['PC3_2']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC4_2']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC1_3']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC2_3']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC3_3']   =   1.000000000000E+00                                                  
            # hdulist[0].header['PC4_3']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC1_4']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC2_4']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC3_4']   =   0.000000000000E+00                                                  
            # hdulist[0].header['PC4_4']   =   1.000000000000E+00                                                  
            
            if ephemeris: 
                hdulist[0].header['CTYPE1']  = 'RA---SIN'                                                            
                hdulist[0].header['CRVAL1']  =   self.ra #'{:02E}'.format(self.ra)                                                
                hdulist[0].header['CDELT1']  =  -1*np.sign(self.ra)*self.pixscale/3600                                                 
                hdulist[0].header['CRPIX1']  =  np.ceil(self.imsize/2) # Reference pixel                                                
                hdulist[0].header['CUNIT1']  = 'deg     '   

                hdulist[0].header['CTYPE2']  = 'DEC--SIN'                                                            
                hdulist[0].header['CRVAL2']  =  self.dec #'{:02E}'.format(self.dec)                                              
                hdulist[0].header['CDELT2']  =  -1*np.sign(self.dec)*self.pixscale/3600                                                   
                hdulist[0].header['CRPIX2']  =   np.ceil(self.imsize/2)                                                  
                hdulist[0].header['CUNIT2']  = 'deg     '
          
            hdulist[0].header['CTYPE3']  = 'FREQ    '                                                            
            hdulist[0].header['CRVAL3']  =   self.obsfrequency                                                 
            hdulist[0].header['CDELT3']  =   1.000000000000E+00   # self.obsfrequency/2.75 # Spectral sampling, should be 1?                                                   
            hdulist[0].header['CRPIX3']  =   1.000000000000E+00   
            hdulist[0].header['CUNIT3']  = 'Hz '

            hdulist[0].header['CTYPE4']  = 'STOKES  '                                                            
            hdulist[0].header['CRVAL4']  =   1.000000000000E+00                                                  
            hdulist[0].header['CDELT4']  =   1.000000000000E+00                                                  
            hdulist[0].header['CRPIX4']  =   1.000000000000E+00                                                  
            hdulist[0].header['CUNIT4']  = '        '  
            # hdulist[0].header['CTYPE5']  = 'T-Daverage    '                                                            
            # hdulist[0].header['CRVAL5']  =  self.T_da
            # hdulist[0].header['CUNIT5']  = 'K     '
            
            # try: 
            #     if self.limbdarkening[0] != self.limbdarkening[1]:
            #         hdulist[0].header['CTYPE6']  = 'Limb darkening coefficient EW '                                                            
            #         hdulist[0].header['CRVAL6']  =  self.limbdarkening[0].tolist() 
            #         hdulist[0].header['CTYPE7']  = 'Limb darkening coefficient NS '                                                            
            #         hdulist[0].header['CRVAL7']  =  self.limbdarkening[1].tolist() 
            #     else: 
            #         hdulist[0].header['CTYPE6']  = 'Limb darkening coefficient '                                                            
            #         hdulist[0].header['CRVAL6']  =  self.limbdarkening[0].tolist() 
            # except: 
            #     hdulist[0].header['CTYPE6']  = 'Limb darkening coefficient '                                                            
            #     hdulist[0].header['CRVAL6']  =  self.limbdarkening



            hdulist[0].header['PV2_1']   =   0.000000000000E+00                                                  
            hdulist[0].header['PV2_2']   =   0.000000000000E+00                                                  
            hdulist[0].header['RESTFRQ'] =   self.obsfrequency #Rest Frequency (Hz)                             
            hdulist[0].header['SPECSYS'] = 'LSRK    '           #Spectral reference frame    
            # http://tdc-www.harvard.edu/wcstools/aips27.pdf Page 1-2 ellaborates on the relation ship                     
            # hdulist[0].header['ALTRVAL'] =  -0.000000000000E+00 #Alternate frequency reference value             
            # hdulist[0].header['ALTRPIX'] =   1.000000000000E+00 #Alternate frequency reference pixel             
            hdulist[0].header['VELREF']  =                  257                  
            #1 LSR, 2 HEL, 3 OBS, +256 Radiocasacore non-standard usage: 4 LSD, 5 GEO, 6 SOU, 7 GAL                 
            if ephemeris: 
                hdulist[0].header['TELESCOP']= 'Model  '                                                            
                hdulist[0].header['OBSERVER']= 'C. Moeckel'  
                modeltime = datetime.strptime(self.time.strip(),'%Y-%b-%d %H:%M:%S.%f')                                      
                hdulist[0].header['DATE-OBS']=  modeltime.strftime("%Y-%m-%dT%H:%M:%S.%f")                                         
                hdulist[0].header['TIMESYS'] = 'UTC     '                                                            
                hdulist[0].header['OBSRA']   =   self.ra                                                  
                hdulist[0].header['OBSDEC']  =   self.dec                                                 
                hdulist[0].header['OBSGEO-X']=  -1.601156673287E+06    # VLA                                              
                hdulist[0].header['OBSGEO-Y']=  -5.041988986066E+06    # VLA                                                
                hdulist[0].header['OBSGEO-Z']=   3.554879236821E+06    # VLA                                              
                hdulist[0].header['OBJECT']  = self.name.upper()                                                             
                hdulist[0].header['TELESCOP']= 'Model    '                                                            
                hdulist[0].header['INSTRUME']= 'Model    '                                                            
                hdulist[0].header['DISTANCE']=   0.000000000000E+00                                                  
            
            hdulist[0].header['DATE']    = now.strftime("%Y-%m-%dT%H:%M:%S.%f") #Date FITS file was written              
            hdulist[0].header['ORIGIN']  = 'pyPR'
        try: 
            outfile = (pathdata+exportname+'.fits')
        except NameError: 
            outfile = (exportname + '.fits')
        # # Export
        hdulist.writeto(outfile, overwrite=True)
        print('Model written to ', outfile)


    def shift_image(self, dx,dy):
        # shift the data by a number of pixels 
        # dx positive to the right (positive RA) 
        # dy positive downwards  (negative RA) 
        self.Bdata = (scipy.ndimage.interpolation.shift(self.Bdata,
            [dy,dx],mode = 'wrap')) 
        self.Bdata[np.abs(self.Bdata)<1e-10] = 0.0

    def shift_center(self, d_p_ra,d_p_dec):
        # shift the right asencsion and declination for units in pixel 
        self.ra = float(self.ra + d_p_ra*self.pixscale/3600)
        self.dec = float(self.dec + d_p_dec*self.pixscale/3600)

    def shift_header(self, d_p_ra,d_p_dec):
        # shift the right asencsion and declination for units in hms and dms
        self.ra = float(self.ra +hms2deg(d_p_ra))
        self.dec = float(self.dec +dms2deg(d_p_dec))

    def set_center(self, p_ra,p_dec):
        # shift the right asencsion and declination for units in hms and dms
        self.ra = float(hms2deg(p_ra))
        self.dec = float(dms2deg(p_dec))

    # def maskforplanet(self,nu,T,p,beamsize,psfsampling=5,): 
    def maskforplanet(self,scalefactor=1.2,export = False,rotangle=0): 
        # rotangle = self.np_ang
        try:
            self.planetsize 
            self.imsize    
        except AttributeError:
            print('Attributes for maskforplanet are missing.'\
                   'Create a model first')
            return None


        planetsize = self.planetsize*scalefactor

        R = self.radius/self.radius[0]*planetsize/2.
        (x,y,) = axisforcasamodel(self.imsize, planetsize)
        (xv,yv,zv) = triaxialellipsoid(R,x,y)

        # # Avoid singularities 
        zv[zv>1] = 1
        zv[zv<1] = 0

        # Rotate the mask 
        self.planetmask = scipy.ndimage.rotate(zv,rotangle,order=0,reshape = False)

        

    def maskforsynchrotron(self,scalefactor=2.5,rotangle=0):
        ''' Create a mask for the Sychrotron

        # rotangle = self.np_ang
        '''
        planetsize = self.planetsize*1.1
        # Half the planet size 
        Phalf = (np.int(planetsize/2.))
        # Potentially a problem with int 
        # Half the image sie 
        Nhalf = (np.int(self.imsize/2.))
        x = np.linspace(-np.pi/2.,np.pi/2.,Phalf*2.)
        y = np.zeros(np.int(self.imsize))
        # set the middle values to non zero a cosine function 
       
        y[Nhalf-Phalf:Nhalf+Phalf] = scalefactor*np.cos(x)**0.75*Phalf

        # xv,yv = (np.meshgrid(np.linspace(0,self.imsize-1,self.imsize),
        #         np.linspace(0,self.imsize-1,self.imsize)))

        mask = np.zeros((self.imsize,self.imsize)) 

        # IPython.embed()
        for i in range(np.int(self.imsize/2.)):
            indices = (np.ones(self.imsize)*(Nhalf-i)-y<0).tolist()
            mask[i,indices] = 1 

            mask[Nhalf:,:] = np.flipud(mask[0:Nhalf])

        self.synchrotronmask = mask.T 
        # Rotate into place based on the north pole angle 

        self.synchrotronmask = (scipy.ndimage.rotate(mask.T,rotangle,
                order=0,reshape = False))   


 

    def plot(self,data,title=''):
        fig,ax = plt.subplots()
        try: 
            x = (np.linspace(-self.imsize/2.,self.imsize/2,self.imsize)
                *self.pixscale)
            y = x
            C = plt.contourf(x,y,data,30)
            plt.title(title)
            plt.xlabel('X [arcseconds]')
            plt.ylabel('Y [arcseconds]')
        except AttributeError:
            C = plt.contourf(data,30)
            plt.title(title)
            plt.xlabel('X [pix]')
            plt.ylabel('Y [pix]')
        try: 
            plt.title('Brightness model: T = {:2.1f} K'.format(self.T_da))
        except AttributeError:
            sys.exit('No temperature defined. Initialize the model parameter')  
   
        cbar = plt.colorbar(C)
        if self.Jansky: 
            cbar.set_label('(Jy/beam)')
        else: 
            cbar.set_label('(K)')
        ax.set_aspect('equal', 'datalim') 
        plt.ion()
        # Add a small reference system for RA and DEC
        xref = self.imsize*0.3*self.pixscale 
        yref = -self.imsize*0.4*self.pixscale
        vl = self.imsize*0.075*self.pixscale # vector length 
        ax.arrow(xref, yref, vl, 0, head_width=1, head_length=3, fc='k', ec='k')
        ax.text(xref,yref-0.65*vl,'RA: {:2.1f}'.format(float(self.ra)))
        ax.text(xref-3*vl,yref+0.75*vl,'DEC: {:2.1f}'.format(float(self.dec)))

        ax.arrow(xref, yref,  0, vl, head_width=1, head_length=3, fc='k', ec='k')
        plt.show() 