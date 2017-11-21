#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Jupiter brightness  


External : numpy, astropy, ephem, time, warnings

"""

import numpy as np
import scipy.ndimage 
import math 
from astropy import constants as cst
import warnings
import PlanetGeometry as pg 
import IPython 
from astropy.io import fits
import matplotlib.pyplot as plt 
import get_ephem
import importlib 
# %matplotlib

pathdata = '/Users/chris/Documents/Research/VLA/VLA2017/'

# Make an oblate spheroid 
def triaxialellipsoid(a,b,c,x,y): 
    """Create an oblate spheroid 

    Create a 3D ellipsoid with a,b,c referring to x,y,z

      y,b   ^ 
            |
            |
            |
             ------> x,a
           / 
          /
         /
      z,c 

    Parameters
    -------
    a [1] float 
        [m] Semi  axis 
    b [1] float 
        [m] Semi  axis
    c [1] float 
        [m] Semi axi 
    x [1] float 
        [m] x axis 
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

    # warnings.filterwarnings("ignore")
    xv, yv = np.meshgrid(x, y) 

    # Create the surface 
    surf = 1. - xv**2/a**2 - yv**2/b**2
    surf[surf<0] = 0.0
    zv= c*np.sqrt(surf)

    # remove nan outside of the planet 
    zv[~np.isfinite(zv)] = 0
    
    return(xv,yv,zv)





# Read in the Jupiter brightness map 
def brightnessmap(R,x,y,T_dab,p): 
    """Create a brightness map of Jupiter

    Based on x,y and grid, create a brightness map of Jupiter. Iterative 
    solution that solves for the brightness distribution 

    Parameters
    -------
    x : [1xN] float
        [R] Normlized x - coordinate  
    y : [1xN] float
        [R] Normlized y - coordinate  
    T_d: [1] float
        Disk averaged brightness temperature 
    R : [3x1] float 
        Radius of the planet [equator,polar,equator]
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
    pg.BrightnessMap(R,N,T_dab,p)

    References
    ------------

    Notes
    -------
    10/24/2017, CM, Initial Commit
    """

    # Call the three dimension surface 
    (xv,yv,zv) = pg.triaxialellipsoid(R[0],R[1],R[2],x,y)

    # Avoid singularities 
    zv[zv==0] = float('nan')
    # Obtain the emission angle 
    th = np.arccos(zv/np.sqrt(xv**2+yv**2+zv**2))


    # Where the value is nan, the model should be zero 
    th[~np.isfinite(th)] = np.pi/2.
    zv[~np.isfinite(zv)] = 0

    # Iterate to make the model agree to the disk averaged 
    # brightness temperature 
    pixels = np.ones_like(zv) 
    cond = 1 
    T_scale = T_dab
    # IPython.embed()
    while cond: 
        cos = np.cos(th)
        cos[cos<1e-5] = 0
        T = T_scale*cos**(p[0]+(np.abs(yv)/R[1])*(p[1]-p[0]))
        T_model = (np.sum(T)/np.sum(pixels[zv>0.0]))       
        if np.abs(T_dab - T_model) < 1: 
            cond = 0 
        else: 
            T_scale = T_scale + (T_dab - T_model)
            print(T_dab,T_model,T_scale)

    return T


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
    >>> pg.tbtojy
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
    # flux = T*2.*kb*nu**2/c**2/jy*np.pi*pg.arcsec2rad(pixscale)**2/(4*np.log(2))
    flux = T*2.*kb*nu**2/c**2/jy*pg.arcsec2rad(pixscale)**2

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
    >>> pg.arcsec2rad

    References
    ------------

    Notes
    -------
    10/24/2017, CM, Initial Commit
    '''
    return arcsec/3600.*np.pi/180.

def axisforcasamodel(imsize, planetsize, pixscale, radius=False): 
    """Obtain the axis for casa, based on the input supplied to tclean 

    Beamsize parameters can be obtained from: 
    https://science.nrao.edu/facilities/vla/docs/manuals/
    oss/performance/resolution
    
    Parameters
    -------
    imsize : [1] int
        [-] Size of the pixels   
    planetsize : [1] int
        [arcsec] size of planet in arcsec
    pixscale : [1] float 
        [arc/pixel] Pixel scale, how many arcsec is one pixel 



    Returns
    ----------
    

    Keywords
    ----------
    radius: [1] Boolean 
        [R] Return the axis in radii  

    Example
    -------
    >>> pg.axisforcasamodel(420,39,10,3.1)
    References
    ------------

    Notes
    -------
    10/24/2017, CM, Initial Commit
    """    
    

    # Radius of Jupiter in pixels 

    r_pla = planetsize/pixscale/2

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



def planetmodel(imsize,planet,pixscale,nu,T_dab,R,p): 


    (x,y,) = ( pg.axisforcasamodel(imsize, planet, pixscale))
    T_model = pg.brightnessmap(R,x,y,T_dab,p)

    # fig,ax = plt.subplots()
    # C = plt.contourf(T_model)
    # cbar = plt.colorbar(C)
    # plt.show()

    Jy_model = pg.tb2jy(T_model,nu,pixscale)

    return Jy_model

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
    >>> pg.calcimsize(210,0.15)

    References
    ------------

    Notes
    -------
    10/26/2017, CM, Initial Commit
    """

    base = 210 # (pixel) minimum size 
 
    PPlant = planet/pixscale # Planet in pixel 

    imsize = np.ceil(2.0*planet/pixscale/base)*base 

    return imsize 

# main loop 
if __name__ == '__main__':

    savemodel = True 

    code = '599'
    tstart = '2017-02-02 08:24'
    tend = '2017-02-02 08:25'
    nstep = '1' 
    obsrvtry = '-5' # VLA

    [data,obscoord,R] = (get_ephem.get_ephemerides(code, tstart, 
                        tend, nstep, obsrvtry))


    print(R)

    # Rotate around the y axis to find the correct tilt 


    # Jupiter 
    jupiter = float(data[0][15]) # (arcseconds) Apparent size of the target
    beamsize = 0.8  # (arcseconds) Size of the beam obtained from casa test run 
    psfsampling = 5 # (-) Number of points per beam
    pixscale = beamsize/psfsampling 
    # Radius of the planet in pixel 
    r_pla = jupiter/pixscale/2


    # Normalize the axis for the triaxial ellipsoid and convert to pixel
    R = R/R[0]*r_pla
    print(R)
    # Rotate around x axis to correct for the sub latitude 
    obsslat = float(data[0][17]) # Sub-observer latitude 

    imsize = pg.calcimsize(jupiter,pixscale) # (pixel)
    print('Use the following parameters for your casa deconvolution:\n ')
    print('When changing the size of image, make sure to load a new header.')
    print('Imsize: ', imsize) 
    print('Cell : ', pixscale)
     

    # Brightness model parameterss
    nu = 22e9 
    T_dab = 132
    p = np.array([0.075,0.065])


    Jy_jupiter = pg.planetmodel(imsize,jupiter,pixscale,nu,T_dab,R,p) 
    # First iteration rotation. Needs automization 
    rotangle = -(np.float(data[0][22]))
    Jy_jupiter = scipy.ndimage.rotate(Jy_jupiter,rotangle,order=0,reshape = False)


    fig,ax = plt.subplots()
    C = plt.contourf(Jy_jupiter)
    cbar = plt.colorbar(C)
    plt.draw()

    # # Uranus 
    # imsize = 420 # (pixel)
    # uranus = 3.6609 # (arcseconds) Size of Jupiter
    # psfsampling = 10 # (-) Number of points per beam 
    # beamsize = 0.2 # (arcseconds) Size of the beam 

    # # Brightness model parameterss
    # nu = 10e9
    # T_dab = 175
    # p = np.array([0.2,0.2])
    # R = np.array([25559,24973,25559])*1.
    # R = R/R[0]
    # # Rotate by 97 deg 
    # R = np.abs(np.dot(pg.rotation_matrix([0,1,0],np.radians(0)), R))

    # Jy_uranus = pg.Planetmodel(imsize,uranus,psfsampling,
    #             beamsize,nu,T_dab,R,p) 



    if savemodel: 
        # Import header infromation from CASA data 
        fname = pathdata+'dirty_S630.fits'
        outfile = pathdata+'ldmodel_T'+str(np.int(T_dab))+'_S'+str(np.int(imsize))+'.fits'
        # # Export 

        im = fits.open(fname,ignore_missing_end=True)
        hdu_out = im
        hdu_out[0].data = Jy_jupiter
        hdu_out[0].header['BUNIT'] = 'Jy/pixel ' 
        hdu_out[0].header['BMIN'] = pixscale/3600 #set beam size equal to one pixel so uvsub doesnt get confused
        hdu_out[0].header['BMAJ'] = pixscale/3600

        hdu_out[0].writeto(outfile, overwrite=True)
        print('Model written to ', outfile)
    else:
        print('Model is not saved')
