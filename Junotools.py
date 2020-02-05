""" Dealing with Juno data 

Notes
-----
05/05/19,CM, Initial Commit
"""

import glob 

import pandas as pd
import numpy as np
import astropy as ap 
import astropy.time as apt

from importlib import reload  

from datetime import datetime, timedelta
import matplotlib.pyplot as plt 
from matplotlib import cm  
cmap = cm.plasma 

import IPython

import warnings, sys

warnings.simplefilter(action = "ignore", category = RuntimeWarning)


# If you want to reload 
    # import pyPR.Junotools as jt
    # reload(jt) 


# Functions to be written 
# Look for PJ products and download them 

def PJFileLocation(path, number): 

    file_I = sorted(glob.glob(path + 'PJ{:d}/PDS/MWR{:02d}RI*'.format(number,number)))
    file_G = sorted(glob.glob(path + 'PJ{:d}/PDS/MWR{:02d}RG*'.format(number,number)))

    # The data might not be read in 

    return file_I,file_G

def ReadPDSProducts(file_I,file_G): 
    """ Read the data products for Juno and return one combined data frame 

    
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
      
    import pyPR.Junotools as jt
    path = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/PJ3/'
    file_I = [path+'MWR03RI2016346160000_R30000_V01.csv',path+'MWR03RI2016346170000_R30000_V01.csv',path+'MWR03RI2016346180000_R30000_V01.csv']
    file_G = [path+'MWR03RG2016346160000_R30000_V01.csv',path+'MWR03RG2016346170000_R30000_V01.csv',path+'MWR03RG2016346180000_R30000_V01.csv']


    pj2 = jt.readPDSproducts(file_I,file_G)

    References
    ------------
    

    Todo
    ----- 
    
    Notes
    -------
    mm/dd/yy, Initials of Author, Short description of update
    """
    
    if  np.size(file_I) !=  np.size(file_G):
        sys.exit('Make sure you have the same number of Intrument and Geo files.')

    if np.size(file_I) < 1: 
        sys.exit('No files found to read in! Check directory!')
        return 
    
    print('Reading calibrated data from PDS')

    # Read in the temperature data 
    instrument = pd.DataFrame() 
    for I,G in zip(file_I,file_G): 
        temp_i = pd.read_csv(I) 
        temp_g = pd.read_csv(G)  
        # Join spacecraft amd instrument data 
        temp = pd.concat([temp_g,temp_i],axis=1,sort=False) 

        instrument = pd.concat([instrument,temp],ignore_index=True)

    return instrument 

# Class for the Channel  
class C:   

    def __init__(self,number): 
        self.channel = number 


class PJ:

    def __init__(self,number): 
        self.PJnumber = number 



    def readdata(self,pathJ,synchrotron= True):

        # Define the 6 channels 
        self.C1 = C(1)
        self.C2 = C(2)
        self.C3 = C(3)
        self.C4 = C(4)
        self.C5 = C(5)
        self.C6 = C(6)

        file_I,file_G = PJFileLocation(pathJ,self.PJnumber)
        pj = ReadPDSProducts(file_I,file_G)


        self.pj = pj 
       # pj['t_utc_doy'].apply(datetime.strptime, args=['%Y-%jT%H:%M:%S.%f'])


        print('Processing calibrated data')

        self.doy_obs    = pj['t_utc_doy']
        self.time       = pj['t_utc_doy'].iloc[:,1].apply(datetime.strptime, args=['%Y-%jT%H:%M:%S.%f']) 
        # Not if throws an error that it should be str and not an float, make sure there is an equivalent GRDR and IRDR for each hour
        self.jd         = self.time.apply(lambda x: x.to_julian_date()).values 
        self.t          =  (self.jd - np.floor(self.jd))*86400 
        self.indices    = (np.linspace(0,len(self.t),len(self.t)+1)) 
        self.indices.astype(int)   
        self.ndata      = len(self.t  )

        # Find the time of perijove pass 
        pj_idx          = pj['range_JnJc'].idxmin() 
        self.jd_pj      = apt.Time(datetime.strptime( pj['t_utc_doy'].values[pj_idx][0],'%Y-%jT%H:%M:%S.%f')).jd     
        self.t_pj       = (self.jd - self.jd_pj) *86400
        self.rotation   = self.t_pj/30

        self.target     = 'Jupiter'
        self.range      = pj['range_JnJc'].values   
        self.r_j        = np.array([71492., 71492., 66854.])
        self.jca        = np.degrees(np.arcsin(self.r_j[0]/self.range)) # Jupiter Central Angle 
        self.beamoffset = 2

        self.ob_lat_c   = pj['PC_lat_JsJnJc'].values # System 3 - Right handed 
        # Does not exist
        # Use PG to convert 
        self.ob_lat_g   = geoc2geod(self.ob_lat_c)

        self.ob_lon     = pj['PC_lon_JsJnJc'].values # System 3 - Right handed 
        self.ob_lon_lh  = np.mod(360 - self.ob_lon,360) # System 3 - Left handed 

        # # Channel 1
        chn = 1
        print(f'Channel {chn}')
        # Antenna information 
        self.C1.hpbw  = 20.6 # deg 
        self.C1.f     = 0.6e9  # [Hz] Center frequency  
        self.C1.ld    = 0.45
        # Boresight information
        self.C1.lat_c = pj['PC_lat_JsB1'].values
        self.C1.lat_g = pj['PG_lat_JsB1'].values
        self.C1.lon   = pj['PC_lon_JsB1'].values
        # Information on the beam location 
        self.C1.alpha_boresight = np.degrees(pj['angle_JnJc_B1'].values) # Angle between boresight and antenna direciton 
        self.C1.planet          = self.C1.alpha_boresight + self.beamoffset*self.C1.hpbw < self.jca   
        self.C1.indices_planet  = self.indices[np.where(self.C1.planet==True)].astype(int)   
        self.C1.sky             = self.C1.alpha_boresight > self.jca + self.beamoffset*self.C1.hpbw  
        self.C1.indices_sky     = self.indices[np.where(self.C1.sky==True)].astype(int)   
        self.C1.skyandplanet    = (self.C1.sky  == False) & (self.C1.planet == False) 
        self.C1.indices_sandp   = self.indices[np.where(self.C1.skyandplanet==True)].astype(int) 
        # Antenna temperatures 
        self.C1.T_a   = np.mean([pj['R1_1TA'].values, pj['R1_2TA'].values],axis=0) 
        self.C1.T_a_lg= pj['R1_1TA'].values 
        self.C1.T_a_hg= pj['R1_2TA'].values 
        self.C1.eangle= pj['emssn_angle_JsB1'].values 
        # Compute the uncleaned zonal average 
        self.C1.T_n, self.C1.p, self.C1.lat_g_i = self.zonalaverage(1)

        if synchrotron: 
            # Synchrotron  
            self.C1.bs_lat_J2000    = np.degrees(np.arctan(pj['S3RH_z_B1'].values/np.sqrt(pj['S3RH_x_B1'].values **2+pj['S3RH_y_B1'].values **2)))
            self.C1.bs_lon_J2000    = np.degrees(np.arctan2(pj['S3RH_y_B1'].values  ,pj['S3RH_x_B1'].values)) 
            self.C1.bs_lat_mag      = np.degrees(np.arctan(pj['JMag_z_B1'].values/np.sqrt(pj['JMag_x_B1'].values **2+pj['JMag_y_B1'].values **2)))
            self.C1.bs_lon_mag      = np.degrees(np.arctan2(pj['JMag_y_B1'].values  ,pj['JMag_x_B1'].values)) 

            # convert to 0 - 180 (south pole = 0) & 0 - 360 
            self.C1.bs_lat_J2000    += 90 
            self.C1.bs_lat_J2000[np.where((self.C1.bs_lat_J2000>179.5)==True)] = 179   
            self.C1.bs_lon_J2000    =  np.mod(self.C1.bs_lon_J2000 ,360)
            self.C1.bs_lon_J2000[np.where((self.C1.bs_lon_J2000>359.5)==True)] = 0   
            self.C1.bs_lat_mag      += 90
            self.C1.bs_lat_mag[np.where((self.C1.bs_lat_mag>179.5)==True)] = 179   
            self.C1.bs_lon_mag      = np.mod(self.C1.bs_lon_mag ,360)      
            self.C1.bs_lon_mag[np.where((self.C1.bs_lon_mag>359.5)==True)] = 0 

            self.C1.sync_mag        = [[[] for i in range(360)] for i in range(180)]
            self.C1.sync_J2000      = [[[] for i in range(360)] for i in range(180)]
            for i in range(self.ndata): 
                if self.C1.sky[i] == True: 
                    self.C1.sync_mag[int(np.round(self.C1.bs_lat_mag[i]))][int(np.round(self.C1.bs_lon_mag[i]))].append(self.C1.T_a[i])
                    self.C1.sync_J2000[int(np.round(self.C1.bs_lat_J2000[i]))][int(np.round(self.C1.bs_lon_J2000[i]))].append(self.C1.T_a[i])

            self.C1.sync_m = np.zeros((180,360))
            self.C1.sync_J = np.zeros((180,360))

            for i1 in range(180):
                for i2 in range(360):  
                    if len(self.C1.sync_mag[i1][i2]) < 1:
                        self.C1.sync_m[i1,i2] = 0 
                    else: 
                         self.C1.sync_m[i1,i2] = np.mean(self.C1.sync_mag[i1][i2])

                    if len(self.C1.sync_J2000[i1][i2]) < 1:
                        self.C1.sync_J[i1,i2] = 0 
                    else: 
                         self.C1.sync_J[i1,i2] = np.mean(self.C1.sync_J2000[i1][i2])




        #Channel 2 
        chn = 2
        print(f'Channel {chn}')
        # Antenna information 
        self.C2.hpbw  = 21.0 # deg 
        self.C2.f     = 1.248e9  # [Hz] Center frequency  
        self.C2.ld    = 0.2
        # Boresight information
        self.C2.lat_c = pj['PC_lat_JsB2'].values
        self.C2.lat_g = pj['PG_lat_JsB2'].values
        self.C2.lon   = pj['PC_lon_JsB2'].values
        # Information on the beam location 
        self.C2.alpha_boresight = np.degrees(pj['angle_JnJc_B2'].values) # Angle between boresight and antenna direciton 
        self.C2.planet          = self.C2.alpha_boresight + self.beamoffset*self.C2.hpbw < self.jca   
        self.C2.indices_planet  = self.indices[np.where(self.C2.planet==True)].astype(int)    
        self.C2.sky             = self.C2.alpha_boresight > self.jca + self.beamoffset*self.C2.hpbw  
        self.C2.indices_sky     = self.indices[np.where(self.C1.sky==True)].astype(int)   
        self.C2.skyandplanet    = (self.C2.sky  == False) & (self.C2.planet == False) 
        self.C2.indices_sandp   = self.indices[np.where(self.C2.skyandplanet==True)].astype(int) 
        # Antenna temperatures 
        self.C2.T_a   = np.mean([pj['R2_1TA'].values, pj['R2_2TA'].values],axis=0) 
        self.C2.T_a_lg= pj['R2_1TA'].values 
        self.C2.T_a_hg= pj['R2_2TA'].values 
        self.C2.eangle= pj['emssn_angle_JsB2'].values 
        # Compute the uncleaned zonal average 
        self.C2.T_n, self.C2.p, self.C2.lat_g_i = self.zonalaverage(2)

        # Synchrotron
        if synchrotron: 
            self.C2.bs_lat_J2000       = np.degrees(np.arctan(pj['S3RH_z_B2'].values/np.sqrt(pj['S3RH_x_B2'].values **2+pj['S3RH_y_B2'].values **2)))
            self.C2.bs_lon_J2000       = np.degrees(np.arctan2(pj['S3RH_y_B2'].values  ,pj['S3RH_x_B2'].values)) 
            self.C2.bs_lat_mag         = np.degrees(np.arctan(pj['JMag_z_B2'].values/np.sqrt(pj['JMag_x_B2'].values **2+pj['JMag_y_B2'].values **2)))
            self.C2.bs_lon_mag         = np.degrees(np.arctan2(pj['JMag_y_B2'].values  ,pj['JMag_x_B2'].values)) 
            # convert to 0 - 180 (south pole = 0) & 0 - 360 
            self.C2.bs_lat_J2000    += 90 
            self.C2.bs_lat_J2000[np.where((self.C2.bs_lat_J2000>179.5)==True)] = 179   
            self.C2.bs_lon_J2000    =  np.mod(self.C2.bs_lon_J2000 ,360)
            self.C2.bs_lon_J2000[np.where((self.C2.bs_lon_J2000>359.5)==True)] = 0   
            self.C2.bs_lat_mag      += 90
            self.C2.bs_lat_mag[np.where((self.C2.bs_lat_mag>179.5)==True)] = 179   
            self.C2.bs_lon_mag      = np.mod(self.C2.bs_lon_mag ,360)      
            self.C2.bs_lon_mag[np.where((self.C2.bs_lon_mag>359.5)==True)] = 0 

            self.C2.sync_mag        = [[[] for i in range(360)] for i in range(180)]
            self.C2.sync_J2000      = [[[] for i in range(360)] for i in range(180)]

            for i in range(self.ndata): 
                if self.C2.sky[i] == True: 
                    self.C2.sync_mag[int(np.round(self.C2.bs_lat_mag[i]))][int(np.round(self.C2.bs_lon_mag[i]))].append(self.C2.T_a[i])
                    self.C2.sync_J2000[int(np.round(self.C2.bs_lat_J2000[i]))][int(np.round(self.C2.bs_lon_J2000[i]))].append(self.C2.T_a[i])

            self.C2.sync_m = np.zeros((180,360))
            self.C2.sync_J = np.zeros((180,360))

            for i1 in range(180):
                for i2 in range(360):  
                    if len(self.C2.sync_mag[i1][i2]) < 1:
                        self.C2.sync_m[i1,i2] = 0 
                    else: 
                         self.C2.sync_m[i1,i2] = np.mean(self.C2.sync_mag[i1][i2])

                    if len(self.C2.sync_J2000[i1][i2]) < 1:
                        self.C2.sync_J[i1,i2] = 0 
                    else: 
                         self.C2.sync_J[i1,i2] = np.mean(self.C2.sync_J2000[i1][i2])


        #Channel 3 
        chn = 3
        print(f'Channel {chn}')
        # Antenna information 
        self.C3.hpbw  = 12.1 #deg 
        self.C3.f     = 2.597e9  # [Hz] Center frequency  
        self.C3.ld    = 0.16

        # Boresight information
        self.C3.lat_c = pj['PC_lat_JsB3'].values
        self.C3.lat_g = pj['PG_lat_JsB3'].values
        self.C3.lon   = pj['PC_lon_JsB3'].values
        # Information on the beam location 
        self.C3.alpha_boresight = np.degrees(pj['angle_JnJc_B3'].values) # Angle between boresight and antenna direciton 
        self.C3.planet          = self.C3.alpha_boresight + self.beamoffset*self.C3.hpbw < self.jca   
        self.C3.indices_planet  = self.indices[np.where(self.C3.planet==True)].astype(int)    
        self.C3.sky             = self.C3.alpha_boresight > self.jca + self.beamoffset*self.C3.hpbw  
        self.C3.indices_sky     = self.indices[np.where(self.C1.sky==True)].astype(int)   
        self.C3.skyandplanet    = (self.C3.sky  == False) & (self.C3.planet == False) 
        self.C3.indices_sandp   = self.indices[np.where(self.C3.skyandplanet==True)].astype(int) 
        # Antenna temperatures
        self.C3.T_a   = pj['R3TA'].values 
        self.C3.eangle= pj['emssn_angle_JsB3'].values 
        # Compute the uncleaned zonal average 
        self.C3.T_n, self.C3.p, self.C3.lat_g_i = self.zonalaverage(3)

        # Synchrotron  
        if synchrotron: 
            self.C3.bs_lat_J2000       = np.degrees(np.arctan(pj['S3RH_z_B3'].values/np.sqrt(pj['S3RH_x_B3'].values **2+pj['S3RH_y_B3'].values **2)))
            self.C3.bs_lon_J2000       = np.degrees(np.arctan2(pj['S3RH_y_B3'].values  ,pj['S3RH_x_B3'].values)) 
            self.C3.bs_lat_mag         = np.degrees(np.arctan(pj['JMag_z_B3'].values/np.sqrt(pj['JMag_x_B3'].values **2+pj['JMag_y_B3'].values **2)))
            self.C3.bs_lon_mag         = np.degrees(np.arctan2(pj['JMag_y_B3'].values  ,pj['JMag_x_B3'].values)) 
            self.C3.sync_mag        = [[[] for i in range(360)] for i in range(180)]
            self.C3.sync_J2000      = [[[] for i in range(360)] for i in range(180)]

            # convert to 0 - 180 (south pole = 0) & 0 - 360 
            self.C3.bs_lat_J2000    += 90 
            self.C3.bs_lat_J2000[np.where((self.C3.bs_lat_J2000>179.5)==True)] = 179   
            self.C3.bs_lon_J2000    =  np.mod(self.C3.bs_lon_J2000 ,360)
            self.C3.bs_lon_J2000[np.where((self.C3.bs_lon_J2000>359.5)==True)] = 0   
            self.C3.bs_lat_mag      += 90
            self.C3.bs_lat_mag[np.where((self.C3.bs_lat_mag>179.5)==True)] = 179   
            self.C3.bs_lon_mag      = np.mod(self.C3.bs_lon_mag ,360)      
            self.C3.bs_lon_mag[np.where((self.C3.bs_lon_mag>359.5)==True)] = 0

            for i in range(self.ndata): 
                if self.C3.sky[i] == True: 
                    self.C3.sync_mag[int(np.round(self.C3.bs_lat_mag[i]))][int(np.round(self.C3.bs_lon_mag[i]))].append(self.C3.T_a[i])
                    self.C3.sync_J2000[int(np.round(self.C3.bs_lat_J2000[i]))][int(np.round(self.C3.bs_lon_J2000[i]))].append(self.C3.T_a[i])

            self.C3.sync_m = np.zeros((180,360))
            self.C3.sync_J = np.zeros((180,360))

            for i1 in range(180):
                for i2 in range(360):  
                    if len(self.C3.sync_mag[i1][i2]) < 1:
                        self.C3.sync_m[i1,i2] = 0 
                    else: 
                         self.C3.sync_m[i1,i2] = np.mean(self.C3.sync_mag[i1][i2])

                    if len(self.C3.sync_J2000[i1][i2]) < 1:
                        self.C3.sync_J[i1,i2] = 0 
                    else: 
                         self.C3.sync_J[i1,i2] = np.mean(self.C3.sync_J2000[i1][i2])

        #Channel 4 
        chn = 4
        print(f'Channel {chn}')
        # Antenna information 
        self.C4.hpbw  = 12.1 # deg 
        self.C4.f     = 5.215e9  # [Hz] Center frequency  
        self.C4.ld    = 0.16
        # Boresight information
        self.C4.lat_c = pj['PC_lat_JsB4'].values
        self.C4.lat_g = pj['PG_lat_JsB4'].values
        self.C4.lon   = pj['PC_lon_JsB4'].values
        # Information on the beam location 
        self.C4.alpha_boresight = np.degrees(pj['angle_JnJc_B4'].values) # Angle between boresight and antenna direciton 
        self.C4.planet          = self.C4.alpha_boresight + self.beamoffset*self.C4.hpbw < self.jca   
        self.C4.indices_planet  = self.indices[np.where(self.C4.planet==True)].astype(int)    
        self.C4.sky             = self.C4.alpha_boresight > self.jca + self.beamoffset*self.C4.hpbw  
        self.C4.indices_sky     = self.indices[np.where(self.C1.sky==True)].astype(int)   
        self.C4.skyandplanet    = (self.C4.sky  == False) & (self.C4.planet == False) 
        self.C4.indices_sandp   = self.indices[np.where(self.C4.skyandplanet==True)].astype(int) 
        # Antenna temperatures
        self.C4.T_a   = pj['R4TA'].values 
        self.C4.eangle= pj['emssn_angle_JsB4'].values 
        # Compute the uncleaned zonal average 
        self.C4.T_n, self.C4.p, self.C4.lat_g_i = self.zonalaverage(4)

        # Synchrotron 
        if synchrotron: 
            self.C4.bs_lat_J2000       = np.degrees(np.arctan(pj['S3RH_z_B4'].values/np.sqrt(pj['S3RH_x_B4'].values **2+pj['S3RH_y_B4'].values **2)))
            self.C4.bs_lon_J2000       = np.degrees(np.arctan2(pj['S3RH_y_B4'].values  ,pj['S3RH_x_B4'].values)) 
            self.C4.bs_lat_mag         = np.degrees(np.arctan(pj['JMag_z_B4'].values/np.sqrt(pj['JMag_x_B4'].values **2+pj['JMag_y_B4'].values **2)))
            self.C4.bs_lon_mag         = np.degrees(np.arctan2(pj['JMag_y_B4'].values  ,pj['JMag_x_B4'].values)) 
            self.C4.sync_mag        = [[[] for i in range(360)] for i in range(180)]
            self.C4.sync_J2000      = [[[] for i in range(360)] for i in range(180)]
            # convert to 0 - 180 (south pole = 0) & 0 - 360 
            self.C4.bs_lat_J2000    += 90 
            self.C4.bs_lat_J2000[np.where((self.C4.bs_lat_J2000>179.5)==True)] = 179   
            self.C4.bs_lon_J2000    =  np.mod(self.C4.bs_lon_J2000 ,360)
            self.C4.bs_lon_J2000[np.where((self.C4.bs_lon_J2000>359.5)==True)] = 0   
            self.C4.bs_lat_mag      += 90
            self.C4.bs_lat_mag[np.where((self.C4.bs_lat_mag>179.5)==True)] = 179   
            self.C4.bs_lon_mag      = np.mod(self.C4.bs_lon_mag ,360)      
            self.C4.bs_lon_mag[np.where((self.C4.bs_lon_mag>359.5)==True)] = 0 

            for i in range(self.ndata):
                if self.C4.sky[i] == True:  
                    self.C4.sync_mag[int(np.round(self.C4.bs_lat_mag[i]))][int(np.round(self.C4.bs_lon_mag[i]))].append(self.C4.T_a[i])
                    self.C4.sync_J2000[int(np.round(self.C4.bs_lat_J2000[i]))][int(np.round(self.C4.bs_lon_J2000[i]))].append(self.C4.T_a[i])

            self.C4.sync_m = np.zeros((180,360))
            self.C4.sync_J = np.zeros((180,360))

            for i1 in range(180):
                for i2 in range(360):  
                    if len(self.C4.sync_mag[i1][i2]) < 1:
                        self.C4.sync_m[i1,i2] = 0 
                    else: 
                         self.C4.sync_m[i1,i2] = np.mean(self.C4.sync_mag[i1][i2])

                    if len(self.C4.sync_J2000[i1][i2]) < 1:
                        self.C4.sync_J[i1,i2] = 0 
                    else: 
                         self.C4.sync_J[i1,i2] = np.mean(self.C4.sync_J2000[i1][i2])


        #Channel 5
        chn = 5
        print(f'Channel {chn}') 
        # Antenna information 
        self.C5.hpbw  = 12.0 # deg 
        self.C5.f     = 10.004e9  # [Hz] Center frequency  
        self.C5.ld    = 0.16
        # Boresight information
        self.C5.lat_c = pj['PC_lat_JsB5'].values
        self.C5.lat_g = pj['PG_lat_JsB5'].values
        self.C5.lon   = pj['PC_lon_JsB5'].values
        # Information on the beam location 
        self.C5.alpha_boresight = np.degrees(pj['angle_JnJc_B5'].values) # Angle between boresight and antenna direciton 
        self.C5.planet          = self.C5.alpha_boresight + self.beamoffset*self.C5.hpbw < self.jca   
        self.C5.indices_planet  = self.indices[np.where(self.C5.planet==True)].astype(int)    
        self.C5.sky             = self.C5.alpha_boresight > self.jca + self.beamoffset*self.C5.hpbw  
        self.C5.indices_sky     = self.indices[np.where(self.C1.sky==True)].astype(int)   
        self.C5.skyandplanet    = (self.C5.sky  == False) & (self.C5.planet == False) 
        self.C5.indices_sandp   = self.indices[np.where(self.C5.skyandplanet==True)].astype(int) 
        # Antenna temperatures 
        self.C5.T_a   = pj['R5TA'].values 
        self.C5.eangle= pj['emssn_angle_JsB5'].values 
        # Compute the uncleaned zonal average 
        self.C5.T_n, self.C5.p, self.C5.lat_g_i = self.zonalaverage(5)

        # Synchrotron
        if synchrotron: 
            self.C5.bs_lat_J2000       = np.degrees(np.arctan(pj['S3RH_z_B5'].values/np.sqrt(pj['S3RH_x_B5'].values **2+pj['S3RH_y_B5'].values **2)))
            self.C5.bs_lon_J2000       = np.degrees(np.arctan2(pj['S3RH_y_B5'].values  ,pj['S3RH_x_B5'].values)) 
            self.C5.bs_lat_mag         = np.degrees(np.arctan(pj['JMag_z_B5'].values/np.sqrt(pj['JMag_x_B5'].values **2+pj['JMag_y_B5'].values **2)))
            self.C5.bs_lon_mag         = np.degrees(np.arctan2(pj['JMag_y_B5'].values  ,pj['JMag_x_B5'].values)) 
            self.C5.sync_mag        = [[[] for i in range(360)] for i in range(180)]
            self.C5.sync_J2000      = [[[] for i in range(360)] for i in range(180)]
            # convert to 0 - 180 (south pole = 0) & 0 - 360 
            self.C5.bs_lat_J2000    += 90 
            self.C5.bs_lat_J2000[np.where((self.C5.bs_lat_J2000>179.5)==True)] = 179   
            self.C5.bs_lon_J2000    =  np.mod(self.C5.bs_lon_J2000 ,360)
            self.C5.bs_lon_J2000[np.where((self.C5.bs_lon_J2000>359.5)==True)] = 0   
            self.C5.bs_lat_mag      += 90
            self.C5.bs_lat_mag[np.where((self.C5.bs_lat_mag>179.5)==True)] = 179   
            self.C5.bs_lon_mag      = np.mod(self.C5.bs_lon_mag ,360)      
            self.C5.bs_lon_mag[np.where((self.C5.bs_lon_mag>359.5)==True)] = 0 
     


            for i in range(self.ndata): 
                if self.C5.sky[i] == True:  
                    self.C5.sync_mag[int(np.round(self.C5.bs_lat_mag[i]))][int(np.round(self.C5.bs_lon_mag[i]))].append(self.C5.T_a[i])
                    self.C5.sync_J2000[int(np.round(self.C5.bs_lat_J2000[i]))][int(np.round(self.C5.bs_lon_J2000[i]))].append(self.C5.T_a[i])
            
            self.C5.sync_m = np.zeros((180,360))
            self.C5.sync_J = np.zeros((180,360))

            for i1 in range(180):
                for i2 in range(360):  
                    if len(self.C5.sync_mag[i1][i2]) < 1:
                        self.C5.sync_m[i1,i2] = 0 
                    else: 
                         self.C5.sync_m[i1,i2] = np.mean(self.C5.sync_mag[i1][i2])

                    if len(self.C5.sync_J2000[i1][i2]) < 1:
                        self.C5.sync_J[i1,i2] = 0 
                    else: 
                         self.C5.sync_J[i1,i2] = np.mean(self.C5.sync_J2000[i1][i2])
 
        #Channel 6 
        chn = 6
        print(f'Channel {chn}')
        # Antenna information 
        self.C6.hpbw  = 10.8 # deg 
        self.C6.f     = 21.900e9  # [Hz] Center frequency  
        self.C6.ld    = 0.08
        # Boresight information
        self.C6.lat_c = pj['PC_lat_JsB6'].values
        self.C6.lat_g = pj['PG_lat_JsB6'].values
        self.C6.lon   = pj['PC_lon_JsB6'].values
        # Information on the beam location 
        self.C6.alpha_boresight = np.degrees(pj['angle_JnJc_B6'].values) # Angle between boresight and antenna direciton 
        self.C6.planet          = self.C6.alpha_boresight + self.beamoffset*self.C6.hpbw < self.jca   
        self.C6.indices_planet  = self.indices[np.where(self.C6.planet==True)].astype(int)    
        self.C6.sky             = self.C6.alpha_boresight > self.jca + self.beamoffset*self.C6.hpbw  
        self.C6.indices_sky     = self.indices[np.where(self.C1.sky==True)].astype(int)   
        self.C6.skyandplanet    = (self.C6.sky  == False) & (self.C6.planet == False) 
        self.C6.indices_sandp   = self.indices[np.where(self.C6.skyandplanet==True)].astype(int) 
        # Antenna temperatures
        self.C6.T_a   = pj['R6TA'].values 
        self.C6.eangle= pj['emssn_angle_JsB6'].values 
        # Compute the uncleaned zonal average 
        self.C6.T_n, self.C6.p, self.C6.lat_g_i = self.zonalaverage(6)
        # Synchrotron 
        if synchrotron: 
            self.C6.bs_lat_J2000       = np.degrees(np.arctan(pj['S3RH_z_B6'].values/np.sqrt(pj['S3RH_x_B6'].values **2+pj['S3RH_y_B6'].values **2)))
            self.C6.bs_lon_J2000       = np.degrees(np.arctan2(pj['S3RH_y_B6'].values  ,pj['S3RH_x_B6'].values)) 
            self.C6.bs_lat_mag         = np.degrees(np.arctan(pj['JMag_z_B6'].values/np.sqrt(pj['JMag_x_B6'].values **2+pj['JMag_y_B6'].values **2)))
            self.C6.bs_lon_mag         = np.degrees(np.arctan2(pj['JMag_y_B6'].values  ,pj['JMag_x_B6'].values)) 
            self.C6.sync_mag        = [[[] for i in range(360)] for i in range(180)]
            self.C6.sync_J2000      = [[[] for i in range(360)] for i in range(180)]
            # convert to 0 - 180 (south pole = 0) & 0 - 360 
            self.C6.bs_lat_J2000    += 90 
            self.C6.bs_lat_J2000[np.where((self.C6.bs_lat_J2000>179.5)==True)] = 179   
            self.C6.bs_lon_J2000    =  np.mod(self.C6.bs_lon_J2000 ,360)
            self.C6.bs_lon_J2000[np.where((self.C6.bs_lon_J2000>359.5)==True)] = 0   
            self.C6.bs_lat_mag      += 90
            self.C6.bs_lat_mag[np.where((self.C6.bs_lat_mag>179.5)==True)] = 179   
            self.C6.bs_lon_mag      = np.mod(self.C6.bs_lon_mag ,360)      
            self.C6.bs_lon_mag[np.where((self.C6.bs_lon_mag>359.5)==True)] = 0 

            for i in range(self.ndata): 
                if self.C6.sky[i] == True:  
                    self.C6.sync_mag[int(np.round(self.C6.bs_lat_mag[i]))][int(np.round(self.C6.bs_lon_mag[i]))].append(self.C6.T_a[i])
                    self.C6.sync_J2000[int(np.round(self.C6.bs_lat_J2000[i]))][int(np.round(self.C6.bs_lon_J2000[i]))].append(self.C6.T_a[i])
            
            self.C6.sync_m = np.zeros((180,360))
            self.C6.sync_J = np.zeros((180,360))

            for i1 in range(180):
                for i2 in range(360):  
                    if len(self.C6.sync_mag[i1][i2]) < 1:
                        self.C6.sync_m[i1,i2] = 0 
                    else: 
                         self.C6.sync_m[i1,i2] = np.mean(self.C6.sync_mag[i1][i2])

                    if len(self.C6.sync_J2000[i1][i2]) < 1:
                        self.C6.sync_J[i1,i2] = 0 
                    else: 
                         self.C6.sync_J[i1,i2] = np.mean(self.C6.sync_J2000[i1][i2])

        # Save the zonal average 
        self.savezonalaverage(pathJ + f'PJ{self.PJnumber}/' + 'ZonalAverage_v1.npz' )        

    def rot2ind(self,rotnumb):
        return np.where(np.floor(self.rotation) == np.floor(rotnumb))[0] 

    def lat2rot(self,lat,geocentric=True): 
        if geocentric:  
            indm =  np.argmin(np.abs(self.ob_lat_c - lat))
        else: 
            indm =  np.argmin(np.abs(self.ob_lat_g - lat))

        return self.rotation[indm]


    def plotTA(self,channel=6, latlim = 45, eanglelim=30, ld=[], colorize='eangle',geocentric=True): 

        # Unfiltered data 
        if geocentric: 
            lat_uf = eval('self.C{:d}.lat_g'.format(channel))
        else: 
            lat_uf = eval('self.C{:d}.lat_c'.format(channel))

        T_a_uf = eval('self.C{:d}.T_a'.format(channel))
        eangle_uf = eval('self.C{:d}.eangle'.format(channel))
        t_uf = self.t_pj
 
        idx = (np.abs(lat_uf)<latlim) & (eangle_uf < eanglelim)
        # Filter the data 
        lat   = lat_uf[idx]
        T_a     = T_a_uf[idx]
        eangle  = eangle_uf[idx]
        t       = t_uf[idx]

        #IPython.embed()

        plt.figure(figsize=(6,9)) 


        # Using contourf to provide my colorbar info, then clearing the figure
        if colorize=='eangle': 
            Cmin = int(np.nanmin(eangle)) 
            Cmax = int(np.nanmax(eangle)) 
            colorvector = eangle 
            cbarstring = r'$\mu$[$^\circ$]'

        elif colorize=='time': 
            Cmin = int(np.nanmin(t)) 
            Cmax = int(np.nanmax(t)) 
            colorvector = t 
            cbarstring = 't_PJ [s]'

        Z = [[0,0],[0,0]]
        levels = range(Cmin,Cmax,1)
        CS3 = plt.contourf(Z, levels, cmap=cmap)
        plt.clf()


        # Plot Juno DPS data 
        n = T_a.size  

        if ld: 
            n = T_a.size  

            for i in range(n):
                plt.plot(T_a[i]/np.cos(np.radians(eangle[i]))**ld,lat[i],'.',alpha = 0.5,color=cmap(((colorvector[i] - Cmin)/(Cmax - Cmin))) )
                
        else:
            for i in range(n):
                plt.plot(T_a[i],lat[i],'.',alpha = 0.5,color=cmap(((colorvector[i] - Cmin)/(Cmax - Cmin))) )


        cb = plt.colorbar(CS3) # using the colorbar info I got from contourf
        cb.set_label(cbarstring, rotation=90)

        # Plot published data 
        plt.ylim([-45,45])
        plt.xlabel('Temperature [K]')
        plt.title('Radio brightness - PJ{:d} - C{:d}'.format(self.PJnumber,channel))
        plt.ylabel('Latitude [deg]') 
        plt.show()

        return 

    def PlotSynchrotronMap(self, channel, Frame_magnetic = True,  Crange=[0,0], latsam=3, lonsam=3, outputname=False, dpi=200): 
        from matplotlib import cm,rcParams 
        cmap = cm.magma

        if Frame_magnetic: 
            M = eval(f'self.C{channel}.sync_m')
        else: 
            M = eval(f'self.C{channel}.sync_J')


        # Minimum antenna temperature 
        if Crange ==[0,0]:
            Cmax = int(np.max(M))
            Cmin = int(np.min(M))
        else: 
            Cmax = int(Crange[1])
            Cmin = int(Crange[0])

        # Synchrotron map 
        fig, faxs = plt.subplots(1, 1,figsize=(16,8))

        Z = [[0,0],[0,0]]
        levels = range(Cmin,Cmax,1)
        CS3 = faxs.contourf(Z, levels, cmap=cmap)
        plt.clf()

        print(f'Working on channel {channel}.')

        for i1 in range(0,180,latsam):
            for i2 in range(0,360,lonsam):  
                if M[i1,i2] > 0:
                    plt.scatter(i2,i1-90,color=cmap(M[i1,i2]/Cmax),alpha=0.3)

        cb = plt.colorbar(CS3) # using the colorbar info I got from contourf
        cbarstring = '[K]'
        cb.set_label(cbarstring, rotation=90)
        faxs.set_xlabel('Longitude S3RH [deg]')
        faxs.set_ylabel('Latitude S3RH [deg]')
        faxs.set_title('Synchrotron map')

        if outputname: 
            print('File written to: ' + outputname + '.png')
            plt.savefig(outputname+ '.png', format='png', transparent = False, dpi=dpi)
            plt.savefig(outputname+ '.pdf', format='png', transparent = False, dpi=dpi)
            plt.savefig(outputname+ '.png', format='png', transparent = False, dpi=dpi)

        return 




    def PlotRotationFootprints( self, mapfile, rotnum, channel,  Crange = [0,0], outputname = False, keepTA = None, dpi = 200, TAlim = None): 

        #path = '/Users/chris/GDrive-UCB/Berkeley/Research/VLA/VLA2016/Jup_x_20161211/'
        #fitsfile = 'Products/Maps/spw2~17/jup-20161211-x_lr_spw2~17_rTb.fits' 
        from pyPR.PlanetGeometry import Map

        from matplotlib import cm,rcParams 
        cmap = cm.magma
        cmap2 = cm.viridis
        rcParams.update({'font.size': 14}) 

        R = [71492., 71492., 66854.]

        V = Map('J')
        V.read_deprojected(mapfile,bandwidth=2.0,fluxcal = 1)
        
        # find the indices coresponding to the rotation 
        ind = np.where(np.floor(self.rotation) == rotnum)[0]  

        # Obtain indices corresponding to all rotations consider for this function 
        if keepTA: 
            inda = np.where(np.floor(self.rotation) == keepTA[0])[0][0]  
            indi = ind[-1] # This correspond to the last point for the given rotation 
            # This corresponds to the last point 
            indb = np.where(np.floor(self.rotation) == keepTA[1])[0][-1]  
            indTA = np.linspace(inda,indi,(indi-inda+1),dtype=np.int32)
            indTA = indTA[~np.isnan(eval('self.C{:d}.eangle[indTA]'.format(channel)))] 

        if keepTA: 
            # Find the maximum emission angles for the given range of rotations 
            eangle = eval('self.C{:d}.eangle[np.linspace(inda,indb,(indb-inda+1),dtype=np.int32)]'.format(channel))
        else: 
            # Find the maximum only within the given rotations 
            eangle = eval('self.C{:d}.eangle[ind]'.format(channel))

        # Set plotting range 
        if Crange ==[0,0]:
            Crange[0] = int(np.nanmin(eangle)- 1) 
            Crange[1] = int(np.nanmax(eangle)+ 1) 

        cbarstring = r'$\mu$[$^\circ$]'



        fig = plt.figure(figsize=(V.n_x/V.n_y*8,6))
        grid = plt.GridSpec(1,2,wspace = 0.4,hspace=-0.1,width_ratios=[4, 1])
        grid.update(hspace=0.01,wspace=0.00)
        # Map of the planet 
        ax1 = plt.subplot(grid[0, 0])
        cs = plt.contourf(V.theta ,V.phi , V.Tb_r, 50 ,cmap = 'gray')
        ax1.set_aspect('equal')
        ax1.set_ylabel('Latitude  [deg]')
        ax1.set_xlabel('Longitude  [deg]')
        ax1.set_title('PJ' + r', $\nu$ = {:2.1f} GHz, $\Delta\nu$ = {:2.1f} GHz'.format(10,2))




        # Brightness temperature plot 
        ax2 = plt.subplot(grid[0, 1])

        Z = [[0,0],[0,0]]
        levels = range(Crange[0],Crange[1],1)
        CS3 = ax2.contourf(Z, levels, cmap=cmap)
        ax2.clear()


        d2r = 180/np.pi

        for j in ind: 
            # Skip all the non relevant data 
            beam = np.radians([eval('self.C{:d}.lon[j]'.format(channel)),eval('self.C{:d}.lat_c[j]'.format(channel))]) 
            obs  = np.radians([self.ob_lon[j],self.ob_lat_c[j]]) 
            r_s = R[0]*R[2]/(np.sqrt((R[2]*np.cos(obs[1]))**2 + (R[0]*np.sin(obs[1]))**2))

            # Fraction for color and alpha 
            frac = (eval('self.C{:d}.eangle[j]'.format(channel)) - Crange[0])/(Crange[1] - Crange[0] )
            if frac>1: 
                frac=0.9999 
            if frac <0: 
                frac =0  
            # print(j,eval('self.C{:d}.eangle[j]'.format(channel)), frac)
            if np.isnan(beam).any(): 
                continue 

            [lon_fp,lat_fp],horizon = ProjectedFootprint(beam,obs,eval('self.C{:d}.hpbw'.format(channel))*2,self.range[j]*1e3,r_s=r_s*1e3)  
            ax1.scatter(lon_fp*d2r,lat_fp*d2r,color=cmap(frac), alpha=0.1) 
            ax1.scatter(lon_fp[horizon]*d2r,lat_fp[horizon]*d2r, color='red',alpha=0.15) 
            ax1.scatter(beam[0]*d2r,beam[1]*d2r,color=cmap2(frac), alpha=0.3)
            ax1.scatter(obs[0]*d2r,obs[1]*d2r,color='navy',alpha=1)
            ax1.set_ylim([-60,60])
            if not keepTA: 
                (ax2.scatter(eval('self.C{:d}.T_a[j]'.format(channel))/(np.cos(eval('self.C{:d}.eangle[j]'.format(channel))/d2r)**eval('self.C{:d}.ld'.format(channel))),eval('self.C{:d}.lat_g[j]'.format(channel))
                , alpha = 0.5*(1-frac),color=cmap(frac) ))

        # Plot the temperature data 
        
        if keepTA:
            for j in indTA: 
                print(frac)
                # Fraction for color and alpha 
                frac = (eval('self.C{:d}.eangle[j]'.format(channel)) - Crange[0])/(Crange[1] - Crange[0] )
                (ax2.scatter(eval('self.C{:d}.T_a[j]'.format(channel))/(np.cos(eval('self.C{:d}.eangle[j]'.format(channel))/d2r)**eval('self.C{:d}.ld'.format(channel))),eval('self.C{:d}.lat_g[j]'.format(channel))
                , alpha = 0.5*(1-frac), color=cmap(frac) ))

        # Adjust the map and add a color bar
        clim = [-30,30]
        locs = [-150,-120,-90,-60,-30,0,30,60,90,120,150] 
        labels=['150','120','90','60','30','0','330','300','270','240','210'] 
        ax1.set_xticks(locs, labels,)  
        #ax1.set_rasterized(True)
        cbaxes = fig.add_axes([0.15,0.1,0.25,0.03])
        cbar = plt.colorbar(cs,ticks = [clim[0],clim[0]/2, 0, clim[1]/2, clim[1]], cax = cbaxes,orientation = 'horizontal')
        cbar.set_label('[K]')

        cax = fig.add_axes([0.92, 0.275, 0.01, 0.45])
        cb = fig.colorbar(CS3,cax=cax) # using the colorbar info I got from contourf
        cb.set_label(cbarstring, rotation=90)


        # Adjust the temperature plot 
        ax2.set_ylim([-60,60])
        if TAlim: 
            ax2.set_xlim([TAlim[0],TAlim[1]])
        ax2.set_xlabel('Temperature [K]')
        ax2.set_title('C{:d}'.format(channel))
        ax2.set_yticks([])  

        plt.subplots_adjust(top=0.725, bottom=0.275)

        if outputname: 
            print('File written to: ' + outputname + '.png')
            plt.savefig(outputname+ '.png', format='png', transparent = False, dpi=dpi)


        return 


#  ffmpeg -framerate 1 -pattern_type glob -i '*_small.png' out.mp4
#    os.system('ffmpeg -f image2  -framerate {:d}  -i {:s}%0{:d}d.png {:s}.mp4'.format(framerate, path2im,width,path2gif))


    def PlotFootprint( self, mapfile, rotnum, fpnum, channel,  Crange = [0,0], outputname = False, keepTA = None, dpi = 200, TAlim = None): 

        #path = '/Users/chris/GDrive-UCB/Berkeley/Research/VLA/VLA2016/Jup_x_20161211/'
        #fitsfile = 'Products/Maps/spw2~17/jup-20161211-x_lr_spw2~17_rTb.fits' 
        from pyPR.PlanetGeometry import Map

        from matplotlib import cm,rcParams 
        cmap = cm.magma
        cmap2 = cm.viridis
        rcParams.update({'font.size': 14}) 

        R = [71492., 71492., 66854.]

        V = Map('J')
        V.read_deprojected(mapfile,bandwidth=2.0,fluxcal = 1)
        
        # find the indices coresponding to the rotation 
        ind = np.where(np.floor(self.rotation) == rotnum)[0] 
        j = ind[fpnum] 

        # Obtain indices corresponding to all rotations  
        if keepTA: 
            inda = np.where(np.floor(self.rotation) == keepTA)[0][0]  
            indb = ind[-1] 
            indTA = np.linspace(inda,indb,(indb-inda+1),dtype=np.int32)
            indTA = indTA[~np.isnan(eval('self.C{:d}.eangle[indTA]'.format(channel)))] 

        if keepTA: 
            eangle = eval('self.C{:d}.eangle[indTA]'.format(channel))
        else: 
            eangle = eval('self.C{:d}.eangle[ind]'.format(channel))

        # Set plotting range 
        if Crange ==[0,0]:
            Crange[0] = int(np.nanmin(eangle)-1)
            Crange[1] = int(np.nanmax(eangle)+1)


        cbarstring = r'$\mu$[$^\circ$]'



        fig = plt.figure(figsize=(V.n_x/V.n_y*8,6))
        grid = plt.GridSpec(1,2,wspace = 0.4,hspace=-0.1,width_ratios=[4, 1])
        grid.update(hspace=0.01,wspace=0.00)
        # Map of the planet 
        ax1 = plt.subplot(grid[0, 0])
        cs = plt.contourf(V.theta ,V.phi , V.Tb_r, 50 ,cmap = 'gray')
        ax1.set_aspect('equal')
        ax1.set_ylabel('Latitude  [deg]')
        ax1.set_xlabel('Longitude  [deg]')
        ax1.set_title('PJ' + r', $\nu$ = {:2.1f} GHz, $\Delta\nu$ = {:2.1f} GHz'.format(10,2))


        # Brightness temperature plot 
        ax2 = plt.subplot(grid[0, 1])

        Z = [[0,0],[0,0]]
        levels = range(Crange[0],Crange[1],1)
        CS3 = ax2.contourf(Z, levels, cmap=cmap)
        ax2.clear()


        d2r = 180/np.pi

        
        # Skip all the non relevant data 
        beam = np.radians([eval('self.C{:d}.lon[j]'.format(channel)),eval('self.C{:d}.lat_c[j]'.format(channel))]) 
        obs  = np.radians([self.ob_lon[j],self.ob_lat_c[j]]) 

        #print('j: {:d}, Latitude: {:2.2f}, Beam location {:2.2f},{:2.2f}'.format(j,obs[1]*d2r,*beam*d2r))
        r_s = R[0]*R[2]/(np.sqrt((R[2]*np.cos(obs[1]))**2 + (R[0]*np.sin(obs[1]))**2))

        if np.isnan(beam).any(): 
            print('No beam planet intercept')
            return 

        [lon_fp,lat_fp],horizon = ProjectedFootprint(beam,obs,eval('self.C{:d}.hpbw'.format(channel))*2,self.range[j]*1e3,r_s=r_s*1e3)  
        
        # Set Colorscheme for the foot prints 
        plotcolor = (eval('self.C{:d}.eangle[j]'.format(channel)) - Crange[0])/(Crange[1] - Crange[0] + 1)
        if  plotcolor < 0 or plotcolor > 1: 
            print('The color for the scatter points is outside of [0,1]. Check Crange')
        ax1.scatter(lon_fp*d2r,lat_fp*d2r,color=cmap(plotcolor),alpha=0.1) 
        ax1.scatter(lon_fp[horizon]*d2r,lat_fp[horizon]*d2r,color='red',alpha=0.3) 
        ax1.scatter(beam[0]*d2r,beam[1]*d2r,color=cmap2(plotcolor),alpha=0.3)
        ax1.scatter(obs[0]*d2r,obs[1]*d2r,color='navy',alpha=1)
        ax1.set_ylim([-60,60])
        

        (ax2.scatter(eval('self.C{:d}.T_a[j]'.format(channel))/(np.cos(eval('self.C{:d}.eangle[j]'.format(channel))/d2r)**eval('self.C{:d}.ld'.format(channel))),eval('self.C{:d}.lat_g[j]'.format(channel))
            , alpha = 1 - plotcolor ,color=cmap(plotcolor)) )

        # (ax2.plot(eval('self.C{:d}.T_a[j]'.format(channel))/(np.cos(eval('self.C{:d}.eangle[j]'.format(channel))/d2r)**eval('self.C{:d}.ld'.format(channel))),eval('self.C{:d}.lat_g[j]'.format(channel))
        #     ,'.', alpha = 1 - plotcolor ,color=cmap(plotcolor)) )

        # Plot the temperature data 
        
        if keepTA:
            for j in indTA: 
                (ax2.scatter(eval('self.C{:d}.T_a[j]'.format(channel))/(np.cos(eval('self.C{:d}.eangle[j]'.format(channel))/d2r)**eval('self.C{:d}.ld'.format(channel))),eval('self.C{:d}.lat_g[j]'.format(channel))
                , alpha = 1 - plotcolor,color=cmap(plotcolor) ))

        # Adjust the map and add a color bar
        clim = [-30,30]
        locs = [-150,-120,-90,-60,-30,0,30,60,90,120,150] 
        labels=['150','120','90','60','30','0','330','300','270','240','210'] 
        ax1.set_xticks(locs, labels,)  
        #ax1.set_rasterized(True)
        cbaxes = fig.add_axes([0.15,0.1,0.25,0.03])
        cbar = plt.colorbar(cs,ticks = [clim[0],clim[0]/2, 0, clim[1]/2, clim[1]], cax = cbaxes,orientation = 'horizontal')
        cbar.set_label('[K]')

        cax = fig.add_axes([0.92, 0.275, 0.01, 0.45])
        cb = fig.colorbar(CS3,cax=cax) # using the colorbar info I got from contourf
        cb.set_label(cbarstring, rotation=90)


        # Adjust the temperature plot 
        ax2.set_ylim([-60,60])
        if TAlim: 
            ax2.set_xlim([TAlim[0],TAlim[1]])
        ax2.set_xlabel('Temperature [K]')
        ax2.set_title('C{:d}'.format(channel))
        ax2.set_yticks([])  

        plt.subplots_adjust(top=0.725, bottom=0.275)

        if outputname: 
            print('File written to: ' + outputname + '.png')
            plt.savefig(outputname+ '.png', format='png', transparent = True, dpi=dpi)


        return 


    def PlotTForFullRotation(self,rotnumb,channel,magneticframe=True): 
        '''
        Make a 3D projection based on 2D slice through the magnetosphere,
        based on a spacecraft data just before and after. 
        Rotation matrix from spacecraft to Mag or IAU required for convolution 
        
        Transform axis into new frame 
        Resample matrix 
        Convolve with beam 


        ''' 

        ind = np.where(np.floor(self.rotation) == rotnumb)[0] 




        ind[np.where(self.eval(f'C{channel}').planet[ind] == True)] 

        indp = ind[np.where(self.eval(f'C{channel}').planet[ind] == True)]   
        inds = ind[np.where(self.eval(f'C{channel}').sky[ind] == True)]
        indsp = ind[np.where(self.eval(f'C{channel}').skyandplanet[ind] == True)]


        plt.figure() 
        plt.plot(self.eval(f'C{channel}').T_a[inds],self.eval(f'C{channel}').bs_lat_mag[inds]-90,'.',label='Sky')
        plt.plot(self.eval(f'C{channel}').T_a[indp],self.eval(f'C{channel}').bs_lat_mag[indp]-90,'.',label='Planet')
        plt.plot(self.eval(f'C{channel}').T_a[indsp],self.eval(f'C{channel}').bs_lat_mag[indsp]-90,'.',label='Overlap')
        plt.xlabel('Antenna temperature') 
        plt.ylabel('Latitude_mag')
        plt.title(f'Rotation: {rotnumb}')
        plt.legend()

        plt.figure() 
        plt.plot(self.eval(f'C{channel}').bs_lon_mag[inds],self.eval(f'C{channel}').bs_lat_mag[inds]-90,'.',label='Sky')
        plt.plot(self.eval(f'C{channel}').bs_lon_mag[indp],self.eval(f'C{channel}').bs_lat_mag[indp]-90,'.',label='Planet')
        plt.plot(self.eval(f'C{channel}').bs_lon_mag[indsp],self.eval(f'C{channel}').bs_lat_mag[indsp]-90,'.',label='Overlap')
        plt.xlabel('Longitude') 
        plt.ylabel('Latitude_mag')
        plt.title(f'Rotation: {rotnumb}')
        plt.legend()

        return 


    def ReduceSynchrotron(self,rotnumb): 
        '''
        Make a 3D projection based on 2D slice through the magnetosphere,
        based on a spacecraft data just before and after. 
        Rotation matrix from spacecraft to Mag or IAU required for convolution 
        
        Transform axis into new frame 
        Resample matrix 
        Convolve with beam 


        ''' 

        # Test case for PJ3 or PJ1 
        import pyPR.Junotools as jt
        import pyPR.PlanetGeometry as pg
        import numpy as np
        import matplotlib.pyplot as plt 

        PJ = jt.PJ(3)
        pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
        PJ.readdata(pathJ) 


        self = PJ

        rotnumb =int( self.lat2rot(40))


        ind = self.rot2ind(rotnumb)
        indp = ind[np.where(self.C2.planet[ind] == True)]   
        inds = ind[np.where(self.C2.sky[ind] == True)]
        indsp = ind[np.where(self.C2.skyandplanet[ind] == True)]



        # Plot this rotation 
        plt.figure() 
        plt.plot(self.C2.T_a[inds],self.C2.bs_lat_mag[inds]-90,'.',label='Sky')
        plt.plot(self.C2.T_a[indp],self.C2.bs_lat_mag[indp]-90,'.',label='Planet')
        plt.plot(self.C2.T_a[indsp],self.C2.bs_lat_mag[indsp]-90,'.',label='Overlap')
        plt.xlabel('Antenna temperature') 
        plt.ylabel('Latitude_mag')
        plt.title(f'Rotation: {int(rotnumb)}')
        plt.legend()

        # Find the corresponding angles 
        for i in range(len(indp)): 
            #i = 0
            bs_mag = self.C2.bs_lat_mag[indp][i]-90
            jca    = self.jca[indp][i]
            bs     = self.C2.alpha_boresight[indp][i]

            # Determine if we look north or south 
            if self.ob_lat_c[indp][i] < self.C2.lat_c[indp][i]: 
                lco = jca+bs
                uco = -(jca-bs) 
            elif self.ob_lat_c[indp][i] > self.C2.lat_c[indp][i]:
                lco = -(jca-bs)
                uco = jca+bs


            # Retrieve the raw data synchtrotron data 
            lat_rs = self.C2.bs_lat_mag[inds]-90 
            T_rs   =  self.C2.T_a[inds]

            # # Retrieve the raw planet data 
            # lat_rp = self.C2.bs_lat_mag[indp]-90 
            # T_rp = self.C2.T_a[indp] 


            # Read in the Juno beam  
            path_J = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/Beam/'
            beam = 'antenna1.txt' 
            path = path_J+beam
            G,t,p=jt.readJunoBeam(path,plotting=False) 


            # Cut the beam and make it a 360 cut through it
            B_cut = np.hstack([np.flipud(10**(0.1*G[:,270])),10**(0.1*G[:,90])]).T
            t_cut = np.degrees(np.hstack([np.flipud(-t),t]).T)



            # Create arraty for interpolating onto beam grid  
            lat_s = -180 - lat_rs[lat_rs<0] 
            T_s_s = T_rs[lat_rs<0]
            lat_n = 180 - lat_rs[lat_rs>0]
            T_s_n = T_rs[lat_rs>0]


            #Combine so that you avoid the discontinuity 
            lat_inp = np.hstack([lat_s,lat_n])
            T_inp   = np.hstack([T_s_s,T_s_n])
            T_scale = np.mean(T_inp)
            # T_inp = np.ones_like(lat_inp)*T_scale

            # Insert the cut offs and make the planet disappear 
            # lat_inp = np.hstack([lat_inp,[lco*1.01,lco,0,uco,uco*1.01]])
            # T_inp = np.hstack([T_inp,[T_scale,0,0,0,T_scale]])

            lat_inp = np.hstack([lat_inp,[lco,0,uco]])
            T_inp = np.hstack([T_inp,[0,0,0]])


            # Sort the data 
            ind_sorted    = np.argsort(lat_inp) 

            # Shift the beam according to the  
            t_cut_r = np.roll(t_cut,int(-bs_mag))
            
            # # Plot new beam 
            # plt.figure()
            # plt.plot((t_cut_r),np.log10(B_cut),'.')
            # plt.plot((t_cut),np.log10(B_cut),'.')

            #  Estimate the impact of the Synchrotron

            # T_s_bf = np.interp(t_cut,lat_inp[ind_sorted],T_inp[ind_sorted],left=np.mean(T_inp),right=np.mean(T_inp))
            T_s_bf = np.interp(t_cut_r,lat_inp[ind_sorted],T_inp[ind_sorted],left=0,right=0)

            fig, ax1 = plt.subplots(1, 1)
            ax1.plot((t_cut_r),np.log10(B_cut),'.',alpha=0.3,label='shifted beam')
            ax1.plot((t_cut),np.log10(B_cut),'.',alpha=0.3,label='original beam')
            ax2 = ax1.twinx()

            ax2.plot(lat_inp,T_inp,'.',color='firebrick',label='Raw data')  
            ax2.plot(t_cut_r, T_s_bf,'.',color='coral',label='Interp data')      
            
            plt.title(f'{i}')
            plt.legend()
            # Convolve the beam cut with the synchortron model 
            Tsync = np.dot(B_cut,T_s_bf)*180    
            print(i,Tsync)

        return Tsync

    def ReduceSynchrotron_old(self,rotnumb): 
        '''
        Make a 3D projection based on 2D slice through the magnetosphere,
        based on a spacecraft data just before and after. 
        Rotation matrix from spacecraft to Mag or IAU required for convolution 
        
        Transform axis into new frame 
        Resample matrix 
        Convolve with beam 


        ''' 

        # Test case for PJ3 or PJ1 
        import pyPR.Junotools as jt
        import pyPR.PlanetGeometry as pg
        import numpy as np
        import matplotlib.pyplot as plt 

        PJ = jt.PJ(3)
        pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
        PJ.readdata(pathJ) 


        self = PJ

        rotnumb =int( self.lat2rot(36))


        ind = self.rot2ind(rotnumb)
        indp = ind[np.where(self.C2.planet[ind] == True)]   
        inds = ind[np.where(self.C2.sky[ind] == True)]
        indsp = ind[np.where(self.C2.skyandplanet[ind] == True)]



        # Plot this rotation 
        plt.figure() 
        plt.plot(self.C2.T_a[inds],self.C2.bs_lat_mag[inds]-90,'.',label='Sky')
        plt.plot(self.C2.T_a[indp],self.C2.bs_lat_mag[indp]-90,'.',label='Planet')
        plt.plot(self.C2.T_a[indsp],self.C2.bs_lat_mag[indsp]-90,'.',label='Overlap')
        plt.xlabel('Antenna temperature') 
        plt.ylabel('Latitude_mag')
        plt.title(f'Rotation: {int(rotnumb)}')
        plt.legend()

        # Find the corresponding angles 
        for i in range(len(indp)): 
            #i = 0
            bs_mag = self.C2.bs_lat_mag[indp][i]-90
            jca    = self.jca[indp][i]
            bs     = self.C2.alpha_boresight[indp][i]

            # Determine if we look north or south 
            if self.ob_lat_c[indp][i] < self.C2.lat_c[indp][i]: 
                lco = jca+bs
                uco = -(jca-bs) 
            elif self.ob_lat_c[indp][i] > self.C2.lat_c[indp][i]:
                lco = -(jca-bs)
                uco = jca+bs


            # Retrieve the raw data synchtrotron data 
            lat_rs = self.C2.bs_lat_mag[inds]-90 
            T_rs   =  self.C2.T_a[inds]

            # # Retrieve the raw planet data 
            # lat_rp = self.C2.bs_lat_mag[indp]-90 
            # T_rp = self.C2.T_a[indp] 


            # Read in the Juno beam  
            path_J = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/Beam/'
            beam = 'antenna1.txt' 
            path = path_J+beam
            G,t,p=jt.readJunoBeam(path,plotting=False) 


            # Cut the beam and make it a 360 cut through it
            B_cut = np.hstack([np.flipud(10**(0.1*G[:,270])),10**(0.1*G[:,90])]).T
            t_cut = np.degrees(np.hstack([np.flipud(-t),t]).T)



            # Create arraty for interpolating onto beam grid  
            lat_s = -180 - lat_rs[lat_rs<0] 
            T_s_s = T_rs[lat_rs<0]
            lat_n = 180 - lat_rs[lat_rs>0]
            T_s_n = T_rs[lat_rs>0]


            #Combine so that you avoid the discontinuity 
            lat_inp = np.hstack([lat_s,lat_n])
            T_inp   = np.hstack([T_s_s,T_s_n])
            T_scale = np.mean(T_inp)
            # T_inp = np.ones_like(lat_inp)*T_scale

            # Insert the cut offs and make the planet disappear 
            # lat_inp = np.hstack([lat_inp,[lco*1.01,lco,0,uco,uco*1.01]])
            # T_inp = np.hstack([T_inp,[T_scale,0,0,0,T_scale]])

            lat_inp = np.hstack([lat_inp,[lco,0,uco]])
            T_inp = np.hstack([T_inp,[0,0,0]])


            # Sort the data 
            ind_sorted    = np.argsort(lat_inp) 

            # Shift the beam south 
            t_cut_r = np.roll(t_cut,int(-bs_mag))
            
            # # Plot new beam 
            # plt.figure()
            # plt.plot((t_cut_r),np.log10(B_cut),'.')
            # plt.plot((t_cut),np.log10(B_cut),'.')

            #  Estimate the impact of the Synchrotron

            # T_s_bf = np.interp(t_cut,lat_inp[ind_sorted],T_inp[ind_sorted],left=np.mean(T_inp),right=np.mean(T_inp))
            T_s_bf = np.interp(t_cut_r,lat_inp[ind_sorted],T_inp[ind_sorted],left=0,right=0)

            fig, ax1 = plt.subplots(1, 1)
            ax1.plot((t_cut_r),np.log10(B_cut),'.',alpha=0.3,label='shifted beam')
            ax1.plot((t_cut),np.log10(B_cut),'.',alpha=0.3,label='original beam')
            ax2 = ax1.twinx()

            ax2.plot(lat_inp,T_inp,'.',color='firebrick',label='Interpolated data')  
            ax2.plot(t_cut, T_s_bf,'.',color='coral',label='Raw data')      
            
            plt.title(f'{i}')
            plt.legend()
            # Convolve the beam cut with the synchortron model 
            Tsync = np.dot(B_cut,T_s_bf)*180    
            print(i,Tsync)

        return Tsync

    def limbinbeamframe(beam):
        '''
        For a given beam direction, calculate where exactly the limb of the planet is 

        '''


        # Central Angle




    def zonalaverage(self, channel, window=20, weight=2, sigmafil=2, smoothing=True, plotting=False, geocentric=True,): 
        """ Calculate the mean profile   

        Based on the emission angle, calculate the weighted mean profile by 
        using the emission angle as a weight. The closer to nadir, the more 
        valuable the information should be. 


        Parameters
        ----------
        T : Nx1 float 
            [K] Temperature, boresight, beam convolved 
        lat : Nx1 float 
            [K] lattidue, planetocentric 
        mu : Nx1 float 
            [rad] Emission angle of the measurement 

        Keyword Arguments
        ----------


        Returns
        -------
        T_z : Nx1 
            [K] Zonally averaged temperatures 

        
        Warnings
        -------
        
        
        Example
        -------

        import pyPR.Junotools as jt
        import matplotlib.pyplot as plt 
        import numpy as np 
        pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
        PJ = jt.PJ(3)
        PJ.readdata(pathJ) 
        self = PJ
    
        self.zonalaverage(1,plotting=True)
        self.zonalaverage(2,plotting=True)
        path_J = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/PJ3/'

        TB_nadir = 'PJ3_TBnadir.dat'
        PJ.p_lat_c= []
        PJ.p_C1 = []
        PJ.p_C2 = []
        PJ.p_C3 = []
        PJ.p_C4 = []
        PJ.p_C5 = []
        PJ.p_C6 = []
        with open(path_J + TB_nadir) as fp: 
            line = fp.readline()
            while line:
                info = (line.split('\t'))
                PJ.p_lat_c =np.append(PJ.p_lat_c,np.float(info[0]))
                PJ.p_C1  =np.append(PJ.p_C1 ,np.float(info[1]))
                PJ.p_C2  =np.append(PJ.p_C2 ,np.float(info[2]))
                PJ.p_C3  =np.append(PJ.p_C3 ,np.float(info[3]))
                PJ.p_C4  =np.append(PJ.p_C4 ,np.float(info[4]))
                PJ.p_C5  =np.append(PJ.p_C5 ,np.float(info[5]))
                PJ.p_C6  =np.append(PJ.p_C6 ,np.float(info[6]))
                line = fp.readline()

        fp.close()

            
        from pyPR.Junotools import * 
        import matplotlib.pyplot as plt 
        import numpy as np 
        pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
        PJ = PJ(3)
        PJ.readdata(pathJ) 
        self = PJ

        channel = 1
        window=20, 
        weight=2, 
        smoothing=True, 
        plotting=False, 
        geocentric=True, 
        p = 0.44


        References
        ------------
        
        Todo
        ----- 

        Notes
        -------
        12/16/19, CM, initial comit 
        """


        import scipy.signal as signal


        # Obtain relevant data 

        indpl   = eval(f'self.C{channel}.indices_planet') 
        if geocentric:
            lat_f   = eval(f'self.C{channel}.lat_c')[indpl]
        else: 
            lat_f = eval(f'self.C{channel}.lat_g')[indpl]
        
        T_f     = eval(f'self.C{channel}.T_a')[indpl]
        mu_f    = eval(f'self.C{channel}.eangle')[indpl]

        # For the midlatitudes remove the Main Belt looking data 
        if channel == 1 or channel ==2: 
            # Remove data from NH that are looking southwards 
            # Condition 1: beam nadir between 30 and 50 
            
            # Conidition 2: beam boresight south of ob_lat
        
            if geocentric:
                ob_lat = eval(f'self.ob_lat_c')[indpl]
            else: 
                ob_lat = eval(f'self.ob_lat_g')[indpl]

            # Find indices where condition is not true (synchrotron filter)
            
            # inp_sf = 





        # Sort the profile according to latitude and discard nan 
        inds    = np.argsort(lat_f) 
        lat_s   = lat_f[inds] 
        T_s     = T_f[inds]
        mu_s    = mu_f[inds]


        # Try fitting nadir brightness temperature and limb darkening 
        T_n, p = fit_Tn_ld(T_s,mu_s,window=window,  weights=weight,)
        T_n_raw = T_n

        # Find median and std of the data based on the 20 rotations closest to the planet 
        lat_l = self.ob_lat_c[self.rot2ind(10)[0]]
        lat_u = self.ob_lat_c[self.rot2ind(-10)[0]]

        # Find indices between those points 
        ind_mean = np.where((lat_s>lat_l) & (lat_s<lat_u))[0]

        # Find mean and std 
        p_mean = np.mean(p[ind_mean])
        p_std = np.std(p[ind_mean])

        # Fit the data using a prescribed p 
        T_n_fp =  fit_Tn(T_s,mu_s,p_mean,window=window,  weights=weight,)
        
        # Find regions where p exceeds (out of limit)
        ind_ol = np.where((p > p_mean + sigmafil*p_std) | (p < p_mean - sigmafil*p_std))

        # Merge the data 
        T_n[ind_ol] = T_n_fp[ind_ol]
        p[ind_ol]   = p_mean

        # Filter the data accordingly
        T_n_raw = T_n


        # Low bandpass filter for the data 

          
        # First, design the Buterworth filter
        N  = 2    # Filter order
        Wn = 0.04 # Cutoff frequency
        B, A = signal.butter(N, Wn, output='ba')
        T_n_lbpf = signal.filtfilt(B,A, T_n)
        p_lbpf = signal.filtfilt(B,A, p)

        if smoothing: 
            T_n = T_n_lbpf
            p = p_lbpf


        # Interpolate onto 0.5deg cut out the first and last x(window size) 
        lat_i = np.arange(-90,90.2,0.2) 
        T_n_i = np.interp(lat_i,lat_s[window:-window],T_n[window:-window],left=np.nan,right=np.nan)
        p_i = np.interp(lat_i,lat_s[window:-window],p[window:-window],left=np.nan,right=np.nan)

        if plotting: 
            path_J = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/PJ3/'

            TB_nadir = 'PJ3_TBnadir.dat'
            PJ.p_lat_c= []
            PJ.p_C1 = []
            PJ.p_C2 = []
            PJ.p_C3 = []
            PJ.p_C4 = []
            PJ.p_C5 = []
            PJ.p_C6 = []
            with open(path_J + TB_nadir) as fp: 
                line = fp.readline()
                while line:
                    info = (line.split('\t'))
                    PJ.p_lat_c =np.append(PJ.p_lat_c,np.float(info[0]))
                    PJ.p_C1  =np.append(PJ.p_C1 ,np.float(info[1]))
                    PJ.p_C2  =np.append(PJ.p_C2 ,np.float(info[2]))
                    PJ.p_C3  =np.append(PJ.p_C3 ,np.float(info[3]))
                    PJ.p_C4  =np.append(PJ.p_C4 ,np.float(info[4]))
                    PJ.p_C5  =np.append(PJ.p_C5 ,np.float(info[5]))
                    PJ.p_C6  =np.append(PJ.p_C6 ,np.float(info[6]))
                    line = fp.readline()

            fp.close()

            self.plotTA(channel=channel,ld=0.0, geocentric=False)
            plt.plot(eval(f'PJ.p_C{channel}'),PJ.p_lat_c ,linewidth = 5, label='Published',alpha=0.3)

            plt.plot(T_n_raw,lat_s,label='wi=20, we=1') 
            plt.plot(T_n_lbpf,lat_s,label='Order=2,frequency=0.04') 
            plt.plot(T_n_i,lat_i,label='Interpolated') 

            plt.legend() 

            fig, axs = plt.subplots(figsize=(4,6))
            axs.plot(p_i,lat_i) 
            if channel == 1:
                axs.set_xlim([0.4,0.6])
            elif channel == 2:
                axs.set_xlim([0.3,0.4])
            elif channel == 3:
                axs.set_xlim([0.15,0.25])
            elif channel == 4:
                axs.set_xlim([0.1,0.22])
            elif channel == 5:
                axs.set_xlim([0.05,0.15])
            elif channel == 6:
                axs.set_xlim([0.0,0.1])
            axs.set_title(f'PJ{self.PJnumber}, C{channel}')

            axs.set_ylim([-45,45])
            axs.set_ylabel('Latitude [deg]')
            axs.set_xlabel('Limb darkening [deg]')


        return  T_n_i, p_i, lat_i

    def filtersynchroton(): 
        '''

        '''




    def PlotFootprintsPerRotation(self,): 

        # For each rotation count the number of planet only footprints
        nrot = int(np.ceil(self.rotation[-1]) - np.floor(self.rotation[0]))
        fppr = np.zeros((6,nrot))
        ob_lat = np.zeros(nrot)
        for irot in range(nrot): 
            # Find indices correspoding to first rotation
            indrot    = self.rot2ind(irot + np.floor(self.rotation[0]))
            ob_lat[irot] = self.ob_lat_c[indrot][0]
            fppr[0,irot] = np.sum(self.C1.planet[indrot])
            fppr[1,irot] = np.sum(self.C2.planet[indrot])
            fppr[2,irot] = np.sum(self.C3.planet[indrot])
            fppr[3,irot] = np.sum(self.C4.planet[indrot])
            fppr[4,irot] = np.sum(self.C5.planet[indrot])
            fppr[5,irot] = np.sum(self.C6.planet[indrot])


        fig, axs = plt.subplots(1, 1)
        rot = np.arange(int(np.floor(self.rotation[0])),int(np.ceil(self.rotation[-1]))) 
        axs.plot(fppr[0,:],ob_lat,label='C1') 
        axs.plot(fppr[1,:],ob_lat,label='C2') 
        axs.plot(fppr[2,:],ob_lat,label='C3') 
        axs.plot(fppr[3,:],ob_lat,label='C4') 
        axs.plot(fppr[4,:],ob_lat,label='C5') 
        axs.plot(fppr[5,:],ob_lat,label='C6') 
        axs.legend()
        axs.set_title('Footprints per rotation')
        axs.set_xlabel('Rotation index')
        axs.set_ylabel('Footprints on the planet')

        return 


    def plot_zonalaverage(self,channel,dlat=1): 
        """ Calculate the mean profile   

        Based on the emission angle, calculate the weighted mean profile by 
        using the emission angle as a weight. The closer to nadir, the more 
        valuable the information should be. 


        Parameters
        ----------
        T : Nx1 float 
            [K] Temperature, boresight, beam convolved 
        lat : Nx1 float 
            [K] lattidue, planetocentric 
        mu : Nx1 float 
            [rad] Emission angle of the measurement 

        Keyword Arguments
        ----------


        Returns
        -------
        T_z : Nx1 
            [K] Zonally averaged temperatures 

        
        Warnings
        -------
        
        
        Example
        -------

        import pyPR.Junotools as jt
        import matplotlib.pyplot as plt 
        import numpy as np 
        pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
        PJ = jt.PJ(19)
        PJ.readdata(pathJ) 
        lat = self.C1.lat_c 
        T = self.C1.T_a
        mu = self.C1.eangle 


        # Debug mode 



        # Juno published 
        # ----------------------------------------------------------------------
        # Import the Juno Temperature - from the paper 
        path_J = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/PJ3/'

        TB_nadir = 'PJ3_TBnadir.dat'
        PJ.p_lat_c= []
        PJ.p_C1 = []
        PJ.p_C2 = []
        PJ.p_C3 = []
        PJ.p_C4 = []
        PJ.p_C5 = []
        PJ.p_C6 = []
        with open(path_J + TB_nadir) as fp: 
            line = fp.readline()
            while line:
                info = (line.split('\t'))
                PJ.p_lat_c =np.append(PJ.p_lat_c,np.float(info[0]))
                PJ.p_C1  =np.append(PJ.p_C1 ,np.float(info[1]))
                PJ.p_C2  =np.append(PJ.p_C2 ,np.float(info[2]))
                PJ.p_C3  =np.append(PJ.p_C3 ,np.float(info[3]))
                PJ.p_C4  =np.append(PJ.p_C4 ,np.float(info[4]))
                PJ.p_C5  =np.append(PJ.p_C5 ,np.float(info[5]))
                PJ.p_C6  =np.append(PJ.p_C6 ,np.float(info[6]))
                line = fp.readline()

        fp.close()

        References
        ------------
        
        Todo
        ----- 

        Notes
        -------
        12/16/19, CM, initial comit 
        """



        lat = eval(f'self.C{channel}.lat_c')  
        T = eval(f'self.C{channel}.T_a')
        mu = eval(f'self.C{channel}.eangle') 



        # Filter by emission angle 

        indpl   = eval(f'self.C{channel}.indices_planet') 
        lat_f   = lat[indpl]
        T_f     = T[indpl]
        mu_f    = mu[indpl]

        # Sort the profile according to latitude and discard nan 
        inds    = np.argsort(lat_f) 
        lat_s   = lat_f[inds] 
        T_s     = T_f[inds]
        mu_s    = mu_f[inds]


        # Try fitting nadir brightness temperature and limb darkening 
    #    T_n_41, p_41 = fit_Tn_ld(T_s,mu_s,window=400, weights=1,)
        T_n_021, p_021 = jt.fit_Tn_ld(T_s,mu_s,window=20,  weights=1,)
        T_n_022, p_022 = jt.fit_Tn_ld(T_s,mu_s,window=20,  weights=2,)
        T_n_01, p_01 = jt.fit_Tn_ld(T_s,mu_s,window=50,  weights=1,)
        T_n_02, p_02 = jt.fit_Tn_ld(T_s,mu_s,window=50,  weights=2,)
        T_n_12, p_12 = jt.fit_Tn_ld(T_s,mu_s,window=100, weights=2,)


        # # Correct for the limb darkening 
        # ld = [0.45,0.2,0.16,0.16,0.16,0.08]

        # T_ld = T_s/np.cos(np.radians(mu_s))**ld[1]

        # # Decide on the width of the averaging kernel 
        # T_a_1, w_a = jt.movingaverage(T_ld,window=100,weights=1/mu_s)
        # T_a_2, w_a = jt.movingaverage(T_ld,window=200,weights=1/mu_s)
        # T_a_4, w_a = jt.movingaverage(T_ld,window=400,weights=1/mu_s)

        # Just for developing purpose 
        # PJ3.plotTA(channel=5,ld=0,colorize='eangle',geocentric=False)
        # plt.plot(PJ3.p_C5,PJ3.p_lat_c,label='Published')
        # plt.plot(T_a_1,lat_s,label=r'$\mu$ - 0.5d')
        # plt.plot(T_a_2,lat_s,label=r'$\mu$ - 1d')
        # plt.plot(T_a_4,lat_s,label=r'$\mu$ - 2d')
        # plt.legend()
        # plt.xlim([170,210])

        PJ.plotTA(channel=channel,ld=0.0, geocentric=False)
        plt.plot(eval(f'PJ.p_C{channel}'),PJ.p_lat_c ,linewidth = 5, label='Published',alpha=0.3)
        #plt.plot(PJ3.p_C2-3,PJ3.p_lat_c,label='Published-3')
        plt.plot(T_n_021,lat_s,label='wi=20, we=1') 
        plt.plot(T_n_022,lat_s,label='wi=20, we=2') 

        plt.plot(T_n_01,lat_s,label='wi=50, we=1') 
        plt.plot(T_n_02,lat_s,label='wi=50, we=2') 
        plt.plot(T_n_12,lat_s,label='wi=100, we=2') 
        #plt.plot(T_n_42,lat_s,label='wi=4, we=2')  
        #plt.plot(T_n_43,lat_s,label='wi=4, we=3')     
       
        plt.legend()
        #plt.xlim([410,480])

        fig, axs = plt.subplots(figsize=(4,6))
        axs.plot(p_01,lat_s) 
        if channel == 1:
            axs.set_xlim([0.4,0.6])
        elif channel == 2:
            axs.set_xlim([0.3,0.4])
        elif channel == 3:
            axs.set_xlim([0.15,0.25])
        elif channel == 4:
            axs.set_xlim([0.1,0.22])
        elif channel == 5:
            axs.set_xlim([0.05,0.15])
        elif channel == 6:
            axs.set_xlim([0.0,0.1])
        axs.set_title(f'PJ{self.PJnumber}, C{channel}')

        axs.set_ylim([-45,45])
        axs.set_ylabel('Latitude [deg]')
        axs.set_xlabel('Limb darkening [deg]')

        # 

        # 

        return  
   

    def savezonalaverage(self, path2save):

        # import glob 
        # if glob.glob(path2save): 
        #     path2save =  
        # Save as npz file 
        (np.savez(path2save,lat=self.C1.lat_g_i, 
                            T1 = self.C1.T_n, 
                            T2 = self.C2.T_n, 
                            T3 = self.C3.T_n, 
                            T4 = self.C4.T_n, 
                            T5 = self.C5.T_n, 
                            T6 = self.C6.T_n, 
                            p1 = self.C1.p, 
                            p2 = self.C2.p, 
                            p3 = self.C3.p, 
                            p4 = self.C4.p, 
                            p5 = self.C5.p, 
                            p6 = self.C6.p, )) 
        
        # Save as txt file 

        return 



def latitude2rotation(self, latitude): 
    return np.argmin(self.ob_lat_c - latitude) 

def readJunoBeam(path, plotting=False): 
    """ Read the Juno beam 
    
    Each file contains contains G which consists of 181 
    columns (θ=0,...,180) and 360 rows (φ=0,...,359). 
    The pattern is given by 10^(0.1*G).

    Parameters
    ----------
    file_id : str
        [-] Name of the file to be read in  
    
    Keyword Arguments
    ----------
    Normalized : Boolean
        [-] Read in normalized beaming pattern if True (Default: False)


    Returns
    -------
    gain : 180x360 
        [Unit] Description of return value
    out2 :
        [Unit] Description of return value
    
    Warnings
    -------
    
    
    Example
    -------
    import pyPR.Junotools as jt
    path_J = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/Beam/'
    beam = 'antenna_normalized6.txt' 
    path = path_J+beam
    G,t,p=jt.readJunoBeam(path) 

    References
    ------------
    https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/JUNO/microwave.html

    Todo
    ----- 
    https://plot.ly/python/3d-surface-plots/ 
    Notes
    -------
    mm/dd/yy, Initials of Author, Short description of update
    """
    
    import mpl_toolkits.mplot3d.axes3d as axes3d


    G = np.zeros([181,360])
    with open(path) as fp: 
        line = fp.readline()
        cnter = 0 
        while line:
            info = (line.split(' '))
            G[:,cnter] = [float(i) for i in info]
            cnter += 1 
            line = fp.readline()

    fp.close()


    theta   = np.radians(np.linspace(0,180,181,endpoint=True)) 
    phi     = np.radians(np.linspace(0,359,360,endpoint=True)) 


    if plotting: 

        THETA, PHI = np.meshgrid(theta, phi)

        X,Y,Z = polar2cart(1, THETA, PHI)

        fig = plt.figure()
        # plt.contourf(X,Y,Z,100,cmap=plt.get_cmap('plasma')) 
        ax = fig.add_subplot(1,1,1, projection='3d')
        plot = ax.plot_surface(
            X, Y, G.T, rstride=1, cstride=1, cmap=plt.get_cmap('plasma'),
            linewidth=0, antialiased=False, alpha=0.05)
        # plot = ax.plot_wireframe(
        #     X, Y, Z, rstride=2, cstride=2, cmap=plt.get_cmap('plasma'),
        #     linewidth=0.01, antialiased=False, alpha=0.15)

        # plot = ax.scatter(
        #     X, Y, Z,  cmap=plt.get_cmap('plasma'),
        #     linewidth=0.01, antialiased=False, alpha=0.15)

        plt.show()

    return G, (theta), (phi)
    #def Junobeam(): 

def plotJunoBeam(path,antenna,Normalized = False ): 
    """ Read the Juno beam 
    
    Each file contains contains G which consists of 181 
    columns (θ=0,...,180) and 360 rows (φ=0,...,359). 
    The pattern is given by 10^(0.1*G).

    Parameters
    ----------
    file_id : str
        [-] Name of the file to be read in  
    
    Keyword Arguments
    ----------
    Normalized : Boolean
        [-] Read in normalized beaming pattern if True (Default: False)


    Returns
    -------
    gain : 180x360 
        [Unit] Description of return value
    out2 :
        [Unit] Description of return value
    
    Warnings
    -------
    
    
    Example
    -------
    import pyPR.Junotools as jt
    path = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/Beam/'
    for antenna in range(6): 
        jt.plotJunoBeam(path,antenna + 1) 

    References
    ------------
    https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/JUNO/microwave.html

    Todo
    ----- 
    https://plot.ly/python/3d-surface-plots/ 
    Notes
    -------
    mm/dd/yy, Initials of Author, Short description of update
    """
    import plotly.graph_objects as go
    if Normalized: 
        beam = 'antenna_normalized{:d}.txt'.format(antenna) 
    else: 
        beam = 'antenna{:d}.txt'.format(antenna) 

    G,t,p=readJunoBeam(path + beam, plotting=False) 

    # Make a mesh for printing out the antenna pattern 
    THETA, PHI = np.meshgrid(t, p)

    X,Y,Z = polar2cart(1, THETA, PHI)

    # Read data from a csv
    z_data = G.T #10**(0.1*G.T)
    fig = go.Figure(data=[go.Surface(x=X,y=Y,z=z_data)])
    fig.update_traces(contours_z=dict(show=True, usecolormap=True,
                              highlightcolor="limegreen", project_z=True))
    fig.update_layout(title='Antanna {:d}'.format(antenna), autosize=True,
                  width=750, height=750,scene_camera_eye=dict(x=1.87, y=0.88, z=-0.44),
                  margin=dict(l=80, r=80, b=65, t=90))

    fig.show()
    if Normalized:
        path2save = path + "BeamingPattern{:d}_normalized.pdf".format(antenna)
    else: 
        path2save = path + "BeamingPattern{:d}.pdf".format(antenna)

    fig.write_image(path2save)
    return 


def ProjectedFootprint(beam,obs,beamsize,dist,r_s=71492e3): 
    """ Calculate the projected footprint outline  

    This function caluclates the outline of a circular beam as seen 
    on the projected Planet. Based on the spacecraft location and the 
    boresight of the sensor, and the beam size of the spacecraft, you 
    can calculate the latitude and longitude region that corresponds to 
    the beam. Note that the input angles takes the full main lobe into 
    account by using a beam with 2 HPWB. 

    Parameters
    ----------
    beam : 2x1 float 
        [rad] beam boresight, longitude, latitude 
    obs :  2x1 float 
        [rad] spacecraft subobserver points, longitude, latitude        
    beamsize : float 
        [deg] Size of the HPBW   
    dist : float
        [m] Distance of the spacecraft from the center of the planet 

    Keyword Arguments
    ----------


    Returns
    -------
    footprint : 2x100
        [rad] Location of the footprint vertices 
    
    Warnings
    -------
    
    
    Example
    -------
    import pyPR.Junotools as jt
    #beam = np.array([-0.59377846, -0.99061398])
    #obs = np.array([-0.27264095, -0.59551803])
    beam = np.radians(np.array([-64, -66]))
    obs = np.radians(np.array([-30, -26]))
    beamsize = 20.6
    dist = 85265.62*1e3 
    [lon_fp,lat_fp],horizon = jt.ProjectedFootprint(beam,obs,beamsize,dist)  
    plt.figure()
    plt.plot(beam[0]*d2r,beam[1]*d2r,'*b',label='beam')  
    plt.plot(obs[0]*d2r,obs[1]*d2r,'*r',label='sc')
    plt.plot(lon_fp*d2r,lat_fp*d2r)
    plt.legend()


    References
    ------------
    James R Wertz. Orbit and Constellation layout and administration (OCDM)
    
    Todo
    ----- 
    Notes
    -------
    08/27/19, CM, initial comit 
    """

    B_lon, B_lat = beam 
    ob_lon,ob_lat_c = obs 

    if np.isnan(B_lon):
        return [np.nan,np.nan],np.nan

    # Obtain the relative angles between the SSP and boresight of the antenna 
    d_lat = (B_lat - ob_lat_c)
    d_lon = (B_lon - ob_lon)

    eta = np.arccos(np.cos(np.abs(d_lat))*np.cos(np.abs(d_lon)))

     # Calculate boresight location in 
    lam_bs =  (np.arccos(np.cos(np.pi/2-beam[1])*np.cos(np.pi/2-obs[1]) 
        + np.sin(np.pi/2-beam[1])*np.sin(np.pi/2-obs[1])*np.cos(d_lon)))
    
    if B_lon>ob_lon:
        phi_sign = -1 
    else:
        phi_sign = 1

    phi_bs = phi_sign*(np.arccos((np.cos(np.pi/2-beam[1]) - np.cos(np.pi/2-obs[1])*np.cos(lam_bs))/ 
            (np.sin(np.pi/2-obs[1])*np.sin(lam_bs))))



    # Calculate the the boresight angles in the spacecraft frame 
    rho = np.arcsin(r_s/dist)  # True horizon angle 
    eps = np.arccos(np.sin(eta)/np.sin(rho)) # Local elevation angle 
    lam = np.pi/2 - eta - eps # Earth angle 
    D   = r_s*np.sin(lam)/np.sin(eta) 
    
    # Convert the lam_bs into eta_bs as seen from the spacecraft 
    eta_bs = np.arctan(np.sin(rho)*np.sin(lam_bs)/(1- np.sin(rho)*np.cos(lam_bs)))


    # As seen from the SSP, calculate the angles corresponding to a circular beam shape 
    B = np.radians(beamsize) # Size of the beam 
    phi_c = np.arange(0,2*np.pi,2*np.pi/100)

    # Calculate the vertices as seen from the spacecraft 
    eta_v = np.arccos(np.cos(eta_bs)*np.cos(B) + np.sin(eta_bs)*np.sin(B)*np.cos(phi_c))
    phi_v = hemisphere(phi_c)*np.arccos((np.cos(B) - np.cos(eta_bs)*np.cos(eta_v))/(np.sin(eta_bs)*np.sin(eta_v)))
    
 
    horizon = np.nan*np.ones_like(eta_v)
    horizon[eta_v>(rho)] = True
    eta_v[eta_v>(rho)]=rho 
    
    # Adjust for the viewing direction (There maybe a sign flip here) 
    phi = (phi_bs-phi_v)

    # Convert spacecraft vertices into planet vertices  
    eps = np.arccos(np.sin(eta_v)/np.sin(rho))
    lam_v = np.pi/2 - eta_v - eps

    # [rad] colatitude of the sub solar point 
    clat_ssp = np.pi/2 - (ob_lat_c)
    clatp,dL = np.zeros_like(phi),np.zeros_like(phi) 
    for i in range(len(phi)): 
        clatp[i] = np.arccos(np.cos(lam_v[i])*np.cos(clat_ssp) + np.sin(lam_v[i])*np.sin(clat_ssp)*np.cos(phi[i])) 
        #dL[i] = np.arccos((np.cos(H[i]) - np.cos(clatp[i])*np.cos(clat_ssp))/(hemisphere(clat_ssp)*np.sin(clat_ssp)*np.sin(clatp[i]))) - np.pi/2*(hemisphere(clat_ssp) -1 ) 
        
        dL[i] = hemisphere(phi[i])*(np.arccos((np.cos(lam_v[i]) - np.cos(clatp[i])*np.cos(clat_ssp))/(hemisphere(clat_ssp)*np.sin(clat_ssp)*np.sin(clatp[i]))) )#- np.pi/2*(hemisphere(clat_ssp) -1 )) 


    lat_fp = np.pi/2 - clatp 
    lon_fp = np.mod((ob_lon) - dL + np.pi,2*np.pi) - np.pi


    return np.array([lon_fp, lat_fp]), ~np.isnan(horizon,dtype=bool)




def polar2cart(r, theta, phi):
    return [
         r * np.sin(theta) * np.cos(phi),
         r * np.sin(theta) * np.sin(phi),
         r * np.cos(theta)
    ]

def hemisphere(alpha): 
    if np.size(alpha) > 1: 
        H = -1*np.ones_like(alpha) 
        idx = (alpha >= 0) & (alpha < np.pi)   
        H[idx] = 1 
        return H 
    else: 
        if (alpha >= 0) & (alpha < np.pi):  
            return 1    
        else: 
            return -1 



def p2z(p,H=37e3): 
    '''
    Calculate the corresponding height for a given pressure above 
    the 10 bar level (-85km ) 
    '''
    return -85e3 - np.log(p/10)*H 


def fit_Tn_ld(a,mu,window=200, weights=2,): 
    '''
    Try a simple linear least squares approach 
    '''

    '''
    T = T_s[12000:13000]
    W = 1/mu_s[12000:13000]
    movingaverage(T,window=20,weights=W)
    '''
    # Make a copy of the array. The beginning and the end of the array will be copied 
    Tn = np.copy(a) 
    ld = np.zeros_like(a)

    # Loop through the array 
    for i in range(len(a)-window): 
        suba = a[i:i+window]
        subb = mu[i:i+window]
        # if weights is not None: 
        #     subw = weights[i:i+window] 
        #     weights[weights<0.03]=0
        #     sumw = np.sum(subw)
        #     subw /= sumw
        # else: 
        #     subw = np.ones(window)/window

        # Build the A matrix 
        A = np.vstack([np.ones_like(suba),np.log(np.cos(np.radians(subb)))]).T  
        # Build the weights 
        W = np.diag(1/subb**weights)
        # Build the solution matrix 
        b = np.log(suba)  

        # Left hand side 
        ATW = np.matmul(A.T,W) 
        ATWA = np.matmul(ATW,A) 
        iATWA = np.linalg.inv(ATWA)

        # Right hand sie 
        ATWb = np.matmul(ATW,b) 

        logT,p = np.matmul(iATWA,ATWb)

        Tn[i+int(window/2)] = np.exp(logT)
        ld[i+int(window/2)] = p


    return Tn,ld


def fit_Tn(a,mu,p,window=200, weights=2,): 
    '''
    Try a simple linear least squares approach 
    '''

    '''
    T = T_s[12000:13000]
    W = 1/mu_s[12000:13000]
    
    '''
    # Make a copy of the array. The beginning and the end of the array will be copied 
    Tn = np.copy(a) 
    ld = np.zeros_like(a)

    # Loop through the array 
    for i in range(len(a)-window): 
        suba = a[i:i+window]
        
        subb = mu[i:i+window]
        # if weights is not None: 
        #     subw = weights[i:i+window] 
        #     weights[weights<0.03]=0
        #     sumw = np.sum(subw)
        #     subw /= sumw
        # else: 
        #     subw = np.ones(window)/window

        # Build the A matrix 
        A = np.vstack([np.ones_like(suba)]).T  
        # Build the weights 
        W = np.diag(1/subb**weights)
        # Build the solution matrix 
        b = (suba)/np.cos(np.radians(mu[i:i+window]))**p  

        # Left hand side 
        ATW = np.matmul(A.T,W) 
        ATWA = np.matmul(ATW,A) 
        iATWA = np.linalg.inv(ATWA)

        # Right hand sie 
        ATWb = np.matmul(ATW,b) 

        T = np.matmul(iATWA,ATWb)

        Tn[i+int(window/2)] = T
        ld[i+int(window/2)] = p


    return Tn



def movingaverage(a,window=200, weights=None,): 
    '''
    T = T_s[12000:13000]
    W = 1/mu_s[12000:13000]
    movingaverage(T,window=20,weights=W)
    '''
    # Make a copy of the array. The beginning and the end of the array will be copied 
    aa = np.copy(a) 
    wa = np.ones_like(a)

    # Loop through the array 
    for i in range(len(a)-window): 
        suba = a[i:i+window]
        if weights is not None: 
            subw = weights[i:i+window] 
            weights[weights<0.03]=0
            sumw = np.sum(subw)
            subw /= sumw
        else: 
            subw = np.ones(window)/window

        aa[i+int(window/2)] = np.dot(suba,subw)
        wa[i+int(window/2)] = sumw


    return aa,wa

def geoc2geod( lat_c, r=714921e3 , f=1-66854/71492): 
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
    10/28/19, CM, Changed order of variables, and made Jupiter the go to planet 
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

    # There will be a nan for zero
    lat_d[np.isnan(lat_d)] = 0 
    
    return lat_d, h