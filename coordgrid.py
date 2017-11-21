#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from .get_ephem import get_ephemerides, naif_lookup
from nirc2_reduce.image import Image
from nirc2_reduce.phot import nearest_stand_filt
from datetime import datetime, timedelta
import warnings
from skimage import feature
from image_registration.chi2_shifts import chi2_shift
from image_registration.fft_tools.shift import shiftnd, shift2d
from pyproj import Proj
from scipy.interpolate import interp2d
#from .fit_gaussian import fitgaussian
from astropy.modeling import models, fitting


def lat_lon(x,y,ob_lon,ob_lat,pixscale_km,np_ang,req,rpol):
    '''Find latitude and longitude on planet given x,y pixel locations and
    planet equatorial and polar radius'''
    np_ang = -np_ang
    x1 = pixscale_km*(np.cos(np.radians(np_ang))*x - np.sin(np.radians(np_ang))*y)
    y1 = pixscale_km*(np.sin(np.radians(np_ang))*x + np.cos(np.radians(np_ang))*y)
    olrad = np.radians(ob_lat)
    
    #set up quadratic equation for ellipsoid
    r2 = (req/rpol)**2
    a = 1 + r2*(np.tan(olrad))**2 #second order
    b = 2*y1*r2*np.sin(olrad) / (np.cos(olrad)**2) #first order
    c = x1**2 + r2*y1**2 / (np.cos(olrad))**2 - req**2 #constant

    radical = b**2 - 4*a*c
    #will equal nan outside planet since radical < 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") #suppresses error for taking sqrt nan
        x3s1=(-b+np.sqrt(radical))/(2*a)
        x3s2=(-b-np.sqrt(radical))/(2*a)
    z3s1=(y1+x3s1*np.sin(olrad))/np.cos(olrad)
    z3s2=(y1+x3s2*np.sin(olrad))/np.cos(olrad)
    odotr1=x3s1*np.cos(olrad)+z3s1*np.sin(olrad)
    odotr2=x3s2*np.cos(olrad)+z3s2*np.sin(olrad)
    #the two solutions are front and rear intersections with planet
    #only want front intersection
    
    #tricky way of putting all the positive solutions into one array
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") #suppresses error for taking < nan
        odotr2[odotr2 < 0] = np.nan
        x3s2[odotr2 < 0] = np.nan
        z3s2[odotr2 < 0] = np.nan
        odotr1[odotr1 < 0] = odotr2[odotr1 < 0]
        x3s1[odotr1 < 0] = x3s2[odotr1 < 0]
        z3s1[odotr1 < 0] = z3s2[odotr1 < 0]
    
    odotr,x3,z3 = odotr1,x3s1,z3s1
    y3 = x1
    r = np.sqrt(x3**2 + y3**2 + z3**2)
    
    #lon_e = np.degrees(np.arctan(y3/x3)) + ob_lon
    lon_e = np.degrees(np.arctan2(x3,y3)-np.pi/2) + ob_lon 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") #suppresses error for taking < nan
        lon_e[lon_e < 0] += 360
    lat_c = np.degrees(np.arcsin(z3/r))
    lat_g = np.degrees(np.arctan(r2*np.tan(np.radians(lat_c))))
    #plt.imshow(lon_e, origin = 'lower left')
    #plt.show()
    return lat_g, lat_c, lon_e

def surface_normal(lat_g, lon_e, ob_lon):
    '''Returns the normal vector to the surface of the planet.
    Take dot product with sub-obs or sub-sun vector to find cosine of emission angle'''
    nx = np.cos(np.radians(lat_g))*np.cos(np.radians(lon_e-ob_lon))
    ny = np.cos(np.radians(lat_g))*np.sin(np.radians(lon_e-ob_lon))
    nz = np.sin(np.radians(lat_g))
    return np.asarray([nx,ny,nz])

def emission_angle(ob_lat, surf_n):
    '''Return the cosine of the emission angle of surface wrt observer'''
    ob = np.asarray([np.cos(np.radians(ob_lat)),0,np.sin(np.radians(ob_lat))])
    return np.dot(surf_n.T,ob)
    
def get_filt_info(filt):
    '''Helper to I/F. Will find flux of sun in given filter'''
    with open('/Users/emolter/Python/nirc2_reduce/filter_passbands/sun_fluxes.txt','r') as f:
        f.readline() #header
        for line in f:
            l = line.split(',')
            fl = l[0].strip(', \n')
            if fl == filt:
                wl = float(l[1].strip(', \n'))
                sun_mag = float(l[2].strip(', \n'))
                return wl, sun_mag
                
def airmass_correction(air_t, air_c, filt):
    '''Helper to I/F. Computes correction factor to photometry based on airmass.
       Multiply 
       air_t is airmass of science target.
       air_c is airmass of calibrator.
       filt options are j, h, k, l, m. Use nearest one for narrowband filters'''
    cdict = {'j': 0.102,
             'h': 0.059,
             'k': 0.088,
             'l': 0.093,
             'm': 0.220} #from https://www2.keck.hawaii.edu/realpublic/inst/nirc/exts.html
    if filt == None:
        return 1.0
    tau = cdict[filt]
    factor = np.exp(tau*air_t)/np.exp(tau*air_c)
    return factor

class CoordGrid:
    
    def __init__(self, infile, lead_string = None, pixscale = 0.009942, req = 24764, rpol = 24341, scope = 'keck'):
        '''Pull ephemeris data, calculate lat and lon. Pixscale in arcsec, req and rpol in km'''
        
        self.im = Image(infile)
        self.req = req
        self.rpol = rpol
        if scope == 'keck':
            self.data = self.im.data
        elif scope == 'vla':
            self.data = self.im.data[0,0,:,:]
        self.pixscale_arcsec = pixscale
        
        #pull and reformat header info
        targ = self.im.header['OBJECT'].split('_')[0]
        targ = targ.split(' ')[0]
        self.target = targ
        date = self.im.header['DATE-OBS']
        if scope == 'vla':
            expstart = date.split('T')[1]
            date = date.split('T')[0]
        elif scope == 'keck':
            expstart = self.im.header['EXPSTART']
        imsize_x = self.data.shape[0]
        imsize_y = self.data.shape[1]
        tstart = datetime.strptime(date+' '+expstart[:5],'%Y-%m-%d %H:%M')
        tend = tstart + timedelta(minutes=1)
        tstart = datetime.strftime(tstart, '%Y-%m-%d %H:%M')
        tend = datetime.strftime(tend, '%Y-%m-%d %H:%M')
        
        #pull ephemeris data
        naif = naif_lookup(targ)
        if scope == 'keck':
            obscode = 568
        elif scope == 'vla':
            obscode = -5
        ephem = get_ephemerides(naif, obscode, tstart, tend, '1 minutes')[0][0] #just the row for start time
        ephem = [val.strip(' ') for val in ephem]
        time = ephem[0]
        ra, dec = ephem[3], ephem[4]
        dra, ddec = float(ephem[5]), float(ephem[6])
        az, el = float(ephem[7]), float(ephem[8])
        self.airmass, extinction = float(ephem[9]), float(ephem[10])
        apmag, sbrt = float(ephem[11]), float(ephem[12])
        self.ang_diam = float(ephem[15])
        self.ob_lon, self.ob_lat = float(ephem[16]), float(ephem[17])
        self.sun_lon, self.sun_lat = float(ephem[18]), float(ephem[19])
        self.np_ang, self.np_dist = float(ephem[20]), float(ephem[21])
        self.sun_dist = float(ephem[22])*1.496e8 #from AU to km
        self.dist = float(ephem[24])*1.496e8 #from AU to km        
        
        self.pixscale_km = self.dist*np.radians(pixscale/3600)
        avg_circumference = 2*np.pi*((req + rpol)/2.0)
        self.deg_per_px = self.pixscale_km * (1/avg_circumference) * 360 #approximate conversion between degrees and pixels at sub-observer point
    
        if lead_string != None:
            #if you already did the edge detection and centering and are loading centered image
            self.lat_g = Image(lead_string+'_latg.fits').data
            self.lon_e = Image(lead_string+'_lone.fits').data
            self.err_x = Image(lead_string+'_errx.fits').data
            self.err_y = Image(lead_string+'_erry.fits').data
            self.centered = self.im.data
            self.model_planet = np.nan_to_num(self.lat_g * 0.0 + 1.0)
        
        else:
            xcen, ycen = int(imsize_x/2), int(imsize_y/2) #pixels at center of planet
            xx = np.arange(imsize_x) - xcen
            yy = np.arange(imsize_y) - ycen
            x,y = np.meshgrid(xx,yy)
            self.lat_g, self.lat_c, self.lon_e = lat_lon(x,y,self.ob_lon,self.ob_lat,self.pixscale_km,self.np_ang,req,rpol)
    
        self.surf_n = surface_normal(self.lat_g, self.lon_e, self.ob_lon)
        self.mu = emission_angle(self.ob_lat, self.surf_n)

    def ioverf(self, filt, flux_per, stand_airmass):
        '''Compute I/F ratio given an image in cts s-1 and a conversion between
        cts s-1 and erg s-1 cm-2 um-1 sr-1'''
        wl, sun_flux_earth = get_filt_info(filt)
        sun_flux = sun_flux_earth * (1/np.pi)*(1.496e8/self.sun_dist)**2 #factor of pi because F = pi*B
        sr_per_px = np.radians(self.pixscale_arcsec/3600)**2
        sun_flux_density = sun_flux * sr_per_px # from erg s-1 cm-2 um-1 sr-1 to erg s-1 cm-2 um-1 px-1
        print('Sun flux density ', sun_flux_density)
        
        #photometry correction
        airmass_filt = nearest_stand_filt[filt]
        air_corr = airmass_correction(self.airmass, stand_airmass, airmass_filt)
        print('Airmass correction ',air_corr)
        self.data = self.data * flux_per * air_corr / sun_flux_density
        if hasattr(self, 'centered'):
            self.centered = self.centered * flux_per * air_corr / sun_flux_density

    def edge_detect(self, low_thresh = 0.01, high_thresh = 0.05, sigma = 5, plot = True):
        '''Uses skimage canny algorithm to find edges of planet, correlates
        that with edges of model, '''
        self.model_planet = np.nan_to_num(self.lat_g * 0.0 + 1.0)
        edges = feature.canny(self.data/np.max(self.data), sigma=sigma, low_threshold = low_thresh, high_threshold = high_thresh)
        model_edges = feature.canny(self.model_planet, sigma=sigma, low_threshold = low_thresh, high_threshold = high_thresh)
    
        [dx,dy,dxerr,dyerr] = chi2_shift(model_edges,edges)
        self.x_shift = -dx #need if we want to shift another filter the same amount
        self.y_shift = -dy
        print('Pixel shift X, Y = ', self.x_shift, self.y_shift)
        
        #error in position on surface is approximately equal to projected error times emission angle - in reality there is some asymmetry that would be important if error bars are large
        #These are lat/lon error in the x-hat and y-hat directions; their magnitude is correct for lat-lon space but their direction is not
        self.err_x = dxerr/(self.deg_per_px*self.mu)
        self.err_y = dyerr/(self.deg_per_px*self.mu)
        print('    Lat/Lon error at sub-obs point in x-hat direction = '+str(dxerr/self.deg_per_px))
        print('    Lat/Lon error at sub-obs point in y-hat direction = '+str(dyerr/self.deg_per_px))   
        
        self.centered = shift2d(self.data,-1*dx,-1*dy)
        self.edges = shift2d(edges,-1*dx,-1*dy)

        if plot:
            fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(10, 5))
            
            ax0.imshow(self.data, origin = 'lower left')
            ax0.set_title('Image')
            
            ax1.imshow(edges, origin = 'lower left')
            ax1.set_title('Canny filter, $\sigma=$%d'%sigma)
            
            ax2.imshow(self.edges, origin = 'lower left', alpha = 0.5)
            ax2.imshow(model_edges, origin = 'lower left', alpha = 0.5)
            ax2.set_title('Overlay model and data')
            
            plt.show()
    
    
    def manual_shift(self,dx,dy):
        self.centered = shift2d(self.data,dx,dy)
        
    def plot_latlon(self):
        '''Make pretty plot of lat_g and lon_e overlaid on planet'''
        fig, (ax0, ax1) = plt.subplots(1,2, figsize = (12,6))
        
        #little circle around planet - now does not depend on self.edges existing
        planetedge = np.copy(self.lat_g)
        nans = np.isnan(planetedge)
        planetedge[np.invert(nans)] = 100
        planetedge[nans] = 0
        
        #latitudes
        ax0.imshow(self.centered, origin = 'lower left')
        levels_lat = np.arange(-90,105,15)
        label_levels_lat = np.arange(-90,60,30)
        ctr_lat = ax0.contour(self.lat_g, levels_lat, colors='white', linewidths=2)
        ax0.clabel(ctr_lat, label_levels_lat, inline=1, inline_spacing = 2, fontsize=16, fmt='%d')
        ax0.contour(planetedge, colors = 'white', linewidths = 1)
        #ax0.set_title('Latitudes', fontsize = 18)
        ax0.get_xaxis().set_ticks([])
        ax0.axes.get_yaxis().set_ticks([])
        
        #longitudes
        ax1.imshow(self.centered, origin = 'lower left')
        #hack here to avoid discontinuity in contours - split longs in half
        with warnings.catch_warnings():
            warnings.simplefilter("ignore") #suppresses error for taking < nan
            lon_e1 = np.copy(self.lon_e)
            lon_e1[lon_e1 >= 180] = np.nan
            lon_e2 = np.copy(self.lon_e)
            lon_e2[lon_e2 < 180] = np.nan
        
        levels_lon = range(0,360,30)
        levels_lon_hack = [1] + list(levels_lon[1:]) #make contour at zero actually 1 - otherwise won't plot it since it's at the edge
        ctr_lon1 = ax1.contour(lon_e1, levels_lon_hack, colors='white', linewidths=2)
        ctr_lon2 = ax1.contour(lon_e2, levels_lon_hack, colors='white', linewidths=2)
                
        fmt = {}
        vals = np.arange(0,360,30)
        for l, v in zip(levels_lon_hack, vals):
            fmt[l] = str(int(v)) #make it so the labels say the right things despite hack
        ax1.clabel(ctr_lon1, levels_lon_hack, fmt = fmt, inline=1, inline_spacing = 2, fontsize=16)
        ax1.clabel(ctr_lon2, levels_lon_hack, fmt = fmt, inline=1, inline_spacing = 2, fontsize=16)
        ax1.contour(planetedge, colors = 'white', linewidths = 1)
        #ax1.set_title('Longitudes', fontsize = 18)
        ax1.get_xaxis().set_ticks([])
        ax1.axes.get_yaxis().set_ticks([])        
                
        plt.tight_layout()
        plt.savefig('lat_lon_overlay.png')
        plt.show()
        
    def write(self, lead_string):
        '''Tertiary data products'''
        hdulist_out = self.im.hdulist
        #centered data
        hdulist_out[0].header['OBJECT'] = self.target+'_CENTERED'
        hdulist_out[0].data = self.centered
        hdulist_out[0].writeto(lead_string + '_centered.fits', overwrite=True)
        #latitudes
        hdulist_out[0].header['OBJECT'] = self.target+'_LATITUDES'
        hdulist_out[0].data = self.lat_g
        hdulist_out[0].writeto(lead_string + '_latg.fits', overwrite=True)
        #longitudes
        hdulist_out[0].header['OBJECT'] = self.target+'_LONGITUDES'
        hdulist_out[0].data = self.lon_e
        hdulist_out[0].writeto(lead_string + '_lone.fits', overwrite=True)
        #errors only exist if edge_detect was run. if manual shift, just ignore
        try:
            #error in x*mu
            hdulist_out[0].header['OBJECT'] = self.target+'_XERR'
            hdulist_out[0].data = self.err_x
            hdulist_out[0].writeto(lead_string + '_errx.fits', overwrite=True)
            #error in y*mu
            hdulist_out[0].header['OBJECT'] = self.target+'_YERR'
            hdulist_out[0].data = self.err_y
            hdulist_out[0].writeto(lead_string + '_erry.fits', overwrite=True)
        except:
            pass

    def bootstrap_func(self, order = 2):
        '''Takes a navigated image, plots flux as function of emission angle,
        fits (to nth order) to the minimum flux vs emission angle curve.
        returns the fit coefficients to IoverF(mu)'''
    
        onplanet = np.copy(self.centered)
        onplanet *= self.model_planet
    
        vals = onplanet.flatten()
        mus = self.mu.flatten()
    
        vals_2 = vals[vals > 0] #remove zeros on outside of image
        mus_2 = mus[vals > 0]
        
        np.savetxt('flux_vs_mu.txt', np.asarray([mus_2, vals_2]))
    
        '''It looks like the minimum value at each emission angle follows a fairly
        regular distribution. Try to make a function fit the bottom of it.'''
        bins = np.arange(0,1.01,0.01)
        x = np.digitize(mus_2, bins)
        mins = [np.min(vals_2[np.where(x == i)]) if len(vals_2[np.where(x == i)]) > 0 else 0.0 for i in range(bins.shape[0])]
        mins = np.asarray(mins)
        bins = bins[np.where(mins > 0)]
        mins = mins[np.where(mins > 0)]
        
        z = np.polyfit(bins,mins, order)
        func = np.poly1d(z)
    
        plt.semilogy(mus_2, vals_2, linestyle = '', marker = '.', color = 'k', markersize = 1)
        plt.semilogy(bins, func(bins), color = 'r')
        #plt.semilogy(bins, z[0]*bins**2 + z[1]*bins*1.1 + z[2])
        plt.xlabel(r'Emission angle $\mu$')
        plt.ylabel('I/F')
        plt.show()
    
        print('Polynomial parameters ... + A2*x**2 + A1*x + A2 = ',z)
        return z
        
    def locate_feature(self):
        plt.imshow(self.centered, origin = 'lower left')
        plt.show()
        print('Define a box around the feature you want to track. Note x,y are reversed in image due to weird Python indexing!')
        pix_l = input('Enter lower left pixel x,y separated by a comma: ')
        pix_u = input('Enter upper right pixel x,y separated by a comma: ')
        
        p0x, p0y = int(pix_l.split(',')[0].strip(', \n')),int(pix_l.split(',')[1].strip(', \n'))
        p1x, p1y = int(pix_u.split(',')[0].strip(', \n')),int(pix_u.split(',')[1].strip(', \n'))
        region = self.centered[p0x:p1x,p0y:p1y]        

        #Brightest spot in feature
        maxloc = np.where(self.centered == np.max(region))
        maxlat, maxlon = self.lat_g[maxloc], self.lon_e[maxloc]
        print('----------')
        print('Brightest spot in feature is at: ')
        print('    X, Y pixel = ', maxloc)
        print('    '+str(maxlat[0])+' latitude, '+str(maxlon[0])+' longitude')
        
        #Gaussian fit
        A0 = np.sum(region)
        x0, y0 = (p1x - p0x)/2, (p1y - p0y)/2
        x_std0, y_std0 = region.shape[0]/2, region.shape[1]/2
        theta0 = 0.0
        g_init = models.Gaussian2D(A0, x0, y0, x_std0, y_std0, theta0)
        a, b = np.mgrid[:region.shape[0], :region.shape[1]]
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, a, b, region)

        #determine useful parameters
        xf, yf = g.x_mean + p0x, g.y_mean + p0y
        latf, lonf = self.lat_g[int(round(xf)),int(round(yf))], self.lon_e[int(round(xf)), int(round(yf))]
        
        #replace Gaussian onto grid we had before
        eval_g = g(a, b)
        g_overlay = np.zeros(self.centered.shape)
        g_overlay[p0x:p1x,p0y:p1y] = eval_g
        
        #estimate lat and lon errors
        frac = 0.5 #can probably constrain the center much more, but depends on morphology of storm over time
        fracmax = np.where(g_overlay/np.max(g_overlay) > frac)
        inlat, inlon = self.lat_g[fracmax], self.lon_e[fracmax]
        minlat, maxlat = np.min(inlat), np.max(inlat)
        minlon, maxlon = np.min(inlon), np.max(inlon) #wrapping issues still present!
        
        lat_errl, lat_erru = np.abs(latf - minlat), np.abs(maxlat - latf) 
        lon_errl, lon_erru = np.abs(lonf - minlon), np.abs(maxlon - lonf)
        
        #print everything useful
        print('----------')
        print('Gaussian fit parameters: ')
        print('    X, Y pixel = ', xf, yf)
        print('    sigma_x, sigma_y (pixels) = ', g.x_stddev + 0, g.y_stddev + 0)
        print('    '+str(latf)+' latitude, '+str(lonf)+' longitude')
        #print('    Latitude bounds: '+str(minlat)+', '+str(maxlat))
        #print('    Longitude bounds: '+str(minlon)+', '+str(maxlon))
        print('    Latitude lower, upper error: '+str(lat_errl)+', '+str(lat_erru))
        print('    Longitude lower, upper error: '+str(lon_errl)+', '+str(lon_erru))
        print('Remember: X is vertical and Y is horizontal!')
        
        print('----------')
        print('Uncertainty in position from edge detection: ')
        print('    Lat/Lon error in x-hat direction = '+str(self.err_x[int(round(xf)),int(round(yf))]))
        print('    Lat/Lon error in y-hat direction = '+str(self.err_y[int(round(xf)),int(round(yf))]))
        
        fig, (ax0) = plt.subplots(1,1, figsize = (8,5))
        ax0.imshow(region, origin = 'lower left')
        levels = [0.5, 0.999]
        cs = ax0.contour(eval_g/np.max(eval_g), levels, colors = 'white')
        plt.show()

        '''Error is more or less related to width of Gaussian because that is
        your uncertainty in where the center of the storm is. There doesnt seem
        to be a better way of doing things for a highly non-Gaussian extended
        structure. Error does not arise from random noise in data so there is no
        point estimating the error in the Gaussian fit. 
        The error will be uneven in +/- due to geometry, especially
        near the limb - should account for that. Will add this error in quadrature
        to the error in fitting the latitude-longitude grid, which I also need 
        to calculate and will also depend on geometry. Need to think about
        whether to include geometry in both places or just once'''
        
            
'''        
    def deproject(self):
        ''''''
        projector = Proj(proj='eqc', lon_0 = 0, lat_ts = 0, a = 24764000.0, b = 24341000.0) #Equidistant Cylindrical (Plate Caree)
        x,y = projector(self.lon_e,self.lat_g,errcheck = False) #no error checking lets invalid values return 1.e30
        x[np.where(x > 1e29)] = np.nan
        y[np.where(y > 1e29)] = np.nan
        x = x.flatten()
        y = y.flatten()
        z = self.centered.flatten()
        interp = interp2d(x,y,z, kind = 'linear')
        print(interp(50000,100000))
        #yind = y.argsort()
        #projected = self.centered[xind]
        #print(projected.shape)
        #
        #plt.imshow(projected, origin = 'lower left')
        #plt.show()

    '''