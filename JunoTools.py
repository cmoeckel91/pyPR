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
cmap = cm.magma 
cmap2 = cm.inferno
cmap3 = cm.viridis
cmap4 = cm.plasma

# Define colors for each channel 
cc1=cmap(0/6)
cc2=cmap(1/6)
cc3=cmap(2/6)
cc4=cmap(3/6)
cc5=cmap(4/6)
cc6=cmap(5/6)

import IPython

import warnings, sys

warnings.simplefilter(action = "ignore", category = RuntimeWarning)

path_GD='/Users/chris/GDrive-UCB/'
path_radio = '/Users/chris/Documents/Research/Toolbox/radio/'

if not glob.glob(path_GD): 
    path_GD = '/Volumes/casa/chris/Google Drive/' 
    path_radio = '/Users/chris.moeckel/Documents/Research/radio/'

pathJ = path_GD + 'Berkeley/Research/Juno/'
path_J = pathJ

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



    def readdata(self,pathJ,synchrotron= True,Footprints=False):

        self.datapath=pathJ+f'PJ{self.PJnumber}/'
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
        self.pj_idx     = pj['range_JnJc'].idxmin() 
        self.jd_pj      = apt.Time(datetime.strptime( pj['t_utc_doy'].values[self.pj_idx][0],'%Y-%jT%H:%M:%S.%f')).jd     
        self.t_pj       = (self.jd - self.jd_pj) *86400
        self.rotation   = self.t_pj/30

        self.target     = 'Jupiter'
        self.range      = pj['range_JnJc'].values   
        self.r_j        = np.array([71492., 71492., 66854.])
        self.jca        = np.degrees(np.arcsin(self.r_j[0]/self.range)) # Jupiter Central Angle 
        self.beamoffset = 2 # Cut off for when the beam is off the planet 

        self.ob_lat_c   = pj['PC_lat_JsJnJc'].values # System 3 - Right handed 
        # Does not exist
        # Use PG to convert 
        self.ob_lat_g   = np.degrees(geoc2geod(np.radians(self.ob_lat_c))[0])

        self.ob_lon     = pj['PC_lon_JsJnJc'].values # System 3 - Right handed 
        self.ob_lon_lh  = np.mod(360 - self.ob_lon,360) # System 3 - Left handed 

        self.RJ         =   np.array([71492., 71492., 66854.])


        # Determine when the actual Perijove pass happens by seeing when the latitude starts wrapping around 
        signchange = np.where((np.diff(np.sign(np.gradient(self.ob_lat_c))) != 0)) 
        # Append the perijove index and sort the array 
        sc = np.sort(np.append(signchange,self.pj_idx)) 
        # Check if beginng or end of Perijove weren't read in 
        if np.max(sc) - self.pj_idx == 0: 
            print('Warning: End of perijove not read in')
            sc = np.sort(np.append(sc,self.indices[-1])) 
        elif np.min(sc) - self.pj_idx == 0:
            print('Warning: Beginning of perijove not read in')
            sc = np.sort(np.append(sc,0)) 
        elif len(sc) == 1: 
            print('Warning: Neither end nor beginning of perijove read in')
            sc = np.array([0,self.pj_idx,self.indices[-1]])


        sc_idx = np.argwhere(sc==self.pj_idx)[0][0] 

        # If pj_ids = 1, we are in the perijove pass 
        pj_ids = np.zeros_like(self.rotation)
        pj_ids[int(sc[sc_idx-1]):int(sc[sc_idx+1])]  = 1 
        # Boolean array with perojove pass information 
        self.pj_ids = np.array(pj_ids, dtype=bool)  

        print(f'Perijove pass at {self.time[self.pj_idx]} lasted from {self.time[np.where(self.pj_ids)[0][0] ]} to {self.time[np.where(self.pj_ids)[0][-1]]}')

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
        self.C1.T_n, self.C1.p, self.C1.lat_c_i, self.C1.w, self.C1.sig_T, self.C1.sig_p = self.zonalaverage(chn)
        self.C1.lat_g_i = np.degrees(geoc2geod(np.radians(self.C1.lat_c_i))[0])

        # Save the footprints 
        if Footprints:
            opname = f'NadirFP-PJ{self.PJnumber}-Ch{chn}'
            self.NadirFootprints(chn,outputname=opname) 

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
        self.C2.T_n, self.C2.p, self.C2.lat_c_i, self.C2.w, self.C2.sig_T, self.C2.sig_p = self.zonalaverage(chn)
        self.C2.lat_g_i = np.degrees(geoc2geod(np.radians(self.C2.lat_c_i))[0])


        # Save the footprints 
        if Footprints:
            opname = f'NadirFP-PJ{self.PJnumber}-Ch{chn}'
            self.NadirFootprints(chn,outputname=opname) 

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
        self.C3.T_n, self.C3.p, self.C3.lat_c_i, self.C3.w, self.C3.sig_T, self.C3.sig_p = self.zonalaverage(chn)
        self.C3.lat_g_i = np.degrees(geoc2geod(np.radians(self.C3.lat_c_i))[0])


        # Save the footprints 
        if Footprints:
            opname = f'NadirFP-PJ{self.PJnumber}-Ch{chn}'
            self.NadirFootprints(chn,outputname=opname) 

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
        self.C4.T_n, self.C4.p, self.C4.lat_c_i, self.C4.w, self.C4.sig_T, self.C4.sig_p = self.zonalaverage(chn)
        self.C4.lat_g_i = np.degrees(geoc2geod(np.radians(self.C4.lat_c_i))[0])

        # Save the footprints 
        if Footprints:
            opname = f'NadirFP-PJ{self.PJnumber}-Ch{chn}'
            self.NadirFootprints(chn,outputname=opname) 

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
        self.C5.T_n, self.C5.p, self.C5.lat_c_i, self.C5.w, self.C5.sig_T, self.C5.sig_p = self.zonalaverage(chn)
        self.C5.lat_g_i = np.degrees(geoc2geod(np.radians(self.C5.lat_c_i))[0])

        # Save the footprints 
        if Footprints:
            opname = f'NadirFP-PJ{self.PJnumber}-Ch{chn}'
            self.NadirFootprints(chn,outputname=opname) 

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
        self.C6.T_n, self.C6.p, self.C6.lat_c_i, self.C6.w, self.C6.sig_T, self.C6.sig_p = self.zonalaverage(chn)
        self.C6.lat_g_i = np.degrees(geoc2geod(np.radians(self.C6.lat_c_i))[0])

        # Save the footprints 
        if Footprints:
            opname = f'NadirFP-PJ{self.PJnumber}-Ch{chn}'
            self.NadirFootprints(chn,outputname=opname) 

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
        self.savezonalaverage(pathJ + f'PJ{self.PJnumber}/' + f'PJ{self.PJnumber}' )        



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



    def NadirFootprints(self,channel,outputname=None): 
        '''
        This function finds for each rotation where there is a beam on the planet, 
        the location of the most nadir looking beam and calculates the outline 
        of that footprint. 
        '''

        import os


        FPcenter = []
        FPoutline = []

        nrot = int(np.ceil(self.rotation[-1]) - np.floor(self.rotation[0]))
        for irot in range(int(np.floor(self.rotation[0])), int(np.ceil(self.rotation[-1]))): 
            indrot = self.rot2ind(irot)
            indpl = np.where(eval(f'self.C{channel}.planet'))    

            # Determine if there are beams on the planet 
            if np.sum(eval(f'self.C{channel}.planet[indrot]')) < 1:
                continue


            # Find the location of the nadir pointing beam 
            inadir = np.argmin(eval(f'self.C{channel}.alpha_boresight[indrot]')) + indrot[0]
            # print(channel, inadir,self.ob_lat_c[inadir], eval(f'self.C{channel}.alpha_boresight[inadir]'))

            # Skip all the non relevant data 
            beam = np.radians([eval('self.C{:d}.lon[inadir]'.format(channel)),eval('self.C{:d}.lat_c[inadir]'.format(channel))]) 
            obs  = np.radians([self.ob_lon[inadir],self.ob_lat_c[inadir]]) 

            r_s = self.RJ[0]*self.RJ[2]/(np.sqrt((self.RJ[2]*np.cos(obs[1]))**2 + (self.RJ[0]*np.sin(obs[1]))**2))

            [lon_fp,lat_fp],horizon = ProjectedFootprint(beam,obs,eval('self.C{:d}.hpbw'.format(channel)),self.range[inadir]*1e3,r_s=r_s*1e3)  
            FPcenter.append(np.degrees(beam))
            FPoutline.append(np.degrees([lon_fp,lat_fp]))
            
        if outputname: 
            fname = 'Footprints/'
            try:
                os.makedirs(self.datapath+fname)
            except FileExistsError:
                pass
       

            (np.savez(self.datapath+fname+outputname,
                        channel = channel,
                        center=FPcenter,
                        outline = FPoutline))

        return FPcenter,FPoutline

    def PlotFootprintsForRotation( self, mapfile, rotnum, channel,  Crange = [0,0], outputname = False, keepTA = None, dpi = 200, TAlim = None): 

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


    def PlotSingleFootprint( self, mapfile, rotnum, fpnum, channel,  Crange = [0,0], outputname = False, keepTA = None, dpi = 200, TAlim = None, eacontour=False): 
        '''

        import pyPR.JunoTools as jt 
        PJnumber = 3
        PJ = jt.PJ(PJnumber)
        pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
        PJ.readdata(pathJ) 
        
        
        # Mapfile 
        path = '/Users/chris/GDrive-UCB/Berkeley/Research/VLA/VLA2016/Jup_x_20161211/'
        fitsfile = 'Products/Maps/spw2~17/jup-20161211-x_lr_spw2~17_rTb.fits' 
        mapfile = path+fitsfile 

        rotnum = -1
        fprange = np.where(np.isfinite(PJ.C2.eangle[PJ.rot2ind(rotnum)]))
        print(f'Rotnumber {fprange[0][0]} - {fprange[0][-1]}')
        
        fpnum = 200
        channel = 2
    

        

        print(f'Emission angle {PJ.C2.eangle[PJ.rot2ind(rotnum)[fpnum]]}')
        PJ.PlotSingleFootprint( mapfile, rotnum, fpnum, channel,eacontour=True)


        '''


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
        ind = self.rot2ind(rotnum)
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
        dist = self.range[j]


        #print('j: {:d}, Latitude: {:2.2f}, Beam location {:2.2f},{:2.2f}'.format(j,obs[1]*d2r,*beam*d2r))
        r_s = R[0]*R[2]/(np.sqrt((R[2]*np.cos(obs[1]))**2 + (R[0]*np.sin(obs[1]))**2))

        if np.isnan(beam).any(): 
            print('No beam planet intercept')
            return 

        [lon_fp,lat_fp],horizon = ProjectedFootprint(beam,obs,eval('self.C{:d}.hpbw'.format(channel))*2,self.range[j]*1e3,r_s=r_s*1e3)  
        
        if eacontour: 
            #compute the ea outlines for the given axis 

            lon_range = np.arange(np.floor(np.nanmin(lon_fp*57.3)),np.ceil(np.nanmax(lon_fp*57.3)),0.5) 
            lat_range = np.arange(np.floor(np.nanmin(lat_fp*57.3)),np.ceil(np.nanmax(lat_fp*57.3)),0.5) 
            sc = [dist,obs[0], obs[1]]
            ea,lo,la = EmissionAngle(sc, np.radians(lon_range),np.radians(lat_range) ) 
            print(f'SC {sc}, Lo {np.floor(np.nanmin(lon_fp*57.3))}-{np.ceil(np.nanmax(lon_fp*57.3))}, La {np.floor(np.nanmin(lat_fp*57.3))} - {np.ceil(np.nanmax(lat_fp*57.3))} ')


        # Set Colorscheme for the foot prints 
        plotcolor = (eval('self.C{:d}.eangle[j]'.format(channel)) - Crange[0])/(Crange[1] - Crange[0] + 1)
        if  plotcolor < 0 or plotcolor > 1: 
            print('The color for the scatter points is outside of [0,1]. Check Crange')
        ax1.scatter(lon_fp*d2r,lat_fp*d2r,color=cmap(plotcolor),alpha=0.1) 
        ax1.scatter(lon_fp[horizon]*d2r,lat_fp[horizon]*d2r,color='red',alpha=0.3) 
        ax1.scatter(beam[0]*d2r,beam[1]*d2r,color=cmap2(plotcolor),alpha=0.3)
        ax1.scatter(obs[0]*d2r,obs[1]*d2r,color='navy',alpha=1)
        ax1.set_ylim([-60,60])
        if eacontour: 
            CS = ax1.contour(np.degrees(lo),np.degrees(la),np.degrees(ea))            
            ax1.clabel(CS, inline=1, fontsize=10,fmt='%d')


        

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


    def PlotTemperatureRotation(self,rotnumb,channel,magneticframe=True): 
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
            path_beam = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/Beam/'
            beam = 'antenna1.txt' 
            path = path_beam+beam
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
            T_scale = np.max(T_inp)
            # T_inp = np.ones_like(lat_inp)*T_scale

            # Insert the cut offs and make the planet disappear 
            lat_inp = np.hstack([lat_inp,[lco*1.01,lco,0,uco,uco*1.01]])
            T_inp = np.hstack([T_inp,[T_scale,0,0,0,T_scale]])

            # lat_inp = np.hstack([lat_inp,[lco,0,uco]])
            # T_inp = np.hstack([T_inp,[0,0,0]])


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

    def zonalaverage(self, channel, window=20, weight=2, sigmafil=2, smoothing=True, plotting=False, geocentric=True,rotlim = 10): 
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

        import pyPR.JunoTools as jt
        import matplotlib.pyplot as plt 
        import numpy as np 
        pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
        PJ = jt.PJ(3)
        PJ.readdata(pathJ) 
        self = PJ
    
        _,_,_ = self.zonalaverage(1,plotting=True)
        _,_,_ = self.zonalaverage(2,plotting=True)
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

            
        import pyPR.JunoTools as jt 
        import matplotlib.pyplot as plt 
        import numpy as np 
        pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
        PJ = jt.PJ(3)
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
        if channel == 1 or channel == 2 or channel == 3: 
            # Remove data from NH that are looking southwards 
            # Condition 1: beam nadir between 30 and 50 
            
            # Conidition 2: beam boresight south of ob_lat
        
            if geocentric:
                ob_lat = eval(f'self.ob_lat_c')[indpl]
            else: 
                ob_lat = eval(f'self.ob_lat_g')[indpl]

            # Find indices where condition is not true (synchrotron filter)
            
            ind_sf = np.where(~((np.abs(lat_f) > 30) & (np.abs(lat_f) < 60) & (np.sign(ob_lat)*(lat_f-ob_lat) < 0))) 
            ind_nsf = np.where(((np.abs(lat_f) > 30) & (np.abs(lat_f) < 60) & (np.sign(ob_lat)*(lat_f-ob_lat) < 0))) 

            lat_f = lat_f[ind_sf]
            T_f = T_f[ind_sf]
            mu_f = mu_f[ind_sf]  




        # Sort the profile according to latitude and discard nan 
        inds    = np.argsort(lat_f) 
        lat_s   = lat_f[inds] 
        T_s     = T_f[inds]
        mu_s    = mu_f[inds]


        # Try fitting nadir brightness temperature and limb darkening 
        T_n, p, sig = fit_Tn_ld(T_s,mu_s, window=window,  weights=weight,)
        T_n_raw = T_n

        # Find median and std of the data based on the 20 rotations closest to the planet 
        lat_l = self.ob_lat_c[self.rot2ind(rotlim)[0]]
        lat_u = self.ob_lat_c[self.rot2ind(-rotlim)[0]]

        # Find indices between those points 
        ind_mean = np.where((lat_s>lat_l) & (lat_s<lat_u))[0]

        # Find mean and std 
        p_mean = np.nanmean(p[ind_mean])
        p_std = np.nanstd(p[ind_mean])

        # Fit the data using a prescribed p 
        T_n_fp, sig_fp =  fit_Tn(T_s,mu_s,p_mean,window=window,  weights=weight,)
        
        # Find regions where p exceeds (out of limit)
        ind_ol = np.where((p > p_mean + sigmafil*p_std) | (p < p_mean - sigmafil*p_std))

        # Merge the data 
        T_n[ind_ol] = T_n_fp[ind_ol]
        p[ind_ol]   = p_mean
        sig[ind_ol,0] = sig_fp[ind_ol].T
        sig[ind_ol,1] = 0 # p was not fit here 

        # Apply a moving average instead of low bandpass filter 
        sig_T = movingaverage(sig[:,0], window=window)[0]
        sig_p = movingaverage(sig[:,1], window=window)[0]

        # Weigh the data according to their deviation from good data 
        weights = np.zeros_like(T_n)  
        # Calculate how far away from the mean we are 
        weights = p_std/np.abs(p - p_mean)
        # Where the mean is formed assign a weight equal to the maximum 
        weights[ind_mean] = np.min(weights)
        # Downweigh where prescribed p is used 
        weights[ind_ol] = np.max(weights)
        # Normalize weights to 1 
        weights /= np.min(weights)
        # Invert the weights 
        w = 1/weights 



        # Filter the data accordingly
        T_n_raw = T_n


        # Low bandpass filter (Buterworth filter) for the data 
        N  = 2    # Filter order
        Wn = 0.04 # Cutoff frequency
        B, A = signal.butter(N, Wn, output='ba')
        T_n_lbpf = signal.filtfilt(B,A, T_n)
        p_lbpf = signal.filtfilt(B,A, p)
        w_lbpf = signal.filtfilt(B,A, w)


        if smoothing: 
            T_n = T_n_lbpf
            p = p_lbpf
            w = w_lbpf



        # Interpolate onto 0.5deg cut out the first and last x(window size) 
        lat_i = np.arange(-90,90.2,0.2) 
        T_n_i = np.interp(lat_i,lat_s[window:-window],T_n[window:-window],left=np.nan,right=np.nan)
        p_i = np.interp(lat_i,lat_s[window:-window],p[window:-window],left=np.nan,right=np.nan)
        w_i = np.interp(lat_i,lat_s[window:-window],w[window:-window],left=np.nan,right=np.nan)
        sig_T_i = np.interp(lat_i,lat_s[window:-window],sig_T[window:-window],left=np.nan,right=np.nan)
        sig_p_i = np.interp(lat_i,lat_s[window:-window],sig_p[window:-window],left=np.nan,right=np.nan)

        # Apply weights to the measurements 


        if plotting:
            # Read in published data  
            path_J = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/PJ3/'
            TB_nadir = 'PJ3_TBnadir.dat'
            self.p_lat_c= []
            self.p_C1 = []
            self.p_C2 = []
            self.p_C3 = []
            self.p_C4 = []
            self.p_C5 = []
            self.p_C6 = []
            with open(path_J + TB_nadir) as fp: 
                line = fp.readline()
                while line:
                    info = (line.split('\t'))
                    self.p_lat_c =np.append(self.p_lat_c,np.float(info[0]))
                    self.p_C1  =np.append(self.p_C1 ,np.float(info[1]))
                    self.p_C2  =np.append(self.p_C2 ,np.float(info[2]))
                    self.p_C3  =np.append(self.p_C3 ,np.float(info[3]))
                    self.p_C4  =np.append(self.p_C4 ,np.float(info[4]))
                    self.p_C5  =np.append(self.p_C5 ,np.float(info[5]))
                    self.p_C6  =np.append(self.p_C6 ,np.float(info[6]))
                    line = fp.readline()

            fp.close()

            self.plotTA(channel=channel,ld=0.0, geocentric=False)
            plt.plot(eval(f'self.p_C{channel}'),self.p_lat_c ,linewidth = 5, label='Published',alpha=0.3)

            plt.plot(T_n_raw,lat_s,label='wi=20, we=1') 
            plt.plot(T_n_lbpf,lat_s,label='Order=2,frequency=0.04') 
            plt.plot(T_n_i,lat_i,label='Interpolated') 

            plt.legend() 

            fig, axs = plt.subplots(figsize=(4,6))
            axs.plot(p_i,lat_i) 
            if channel == 1:
                axs.set_xlim([0.4,0.6])
            elif channel == 2:
                axs.set_xlim([0.25,0.35])
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

        return  T_n_i, p_i, lat_i, w_i, sig_T_i, sig_p_i


    def PlotZonalCoverage(self,): 
        '''
        Visualizes the number of footprints as a function of latitude. 

        This function caluclates how many footprints can be found as a 
        function of latitude that are at least 2 HPBW on the planet. 
        Observation that are close to the limb are heavily contaminated by 
        Synchrotron and the Cold sky, causing the limb darkening fit to 
        deteriorate. 

        '''
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

    def PlotEangle(self,): 
        '''
        Visualizes the emission angle as a function of latitude  

        This function caluclates how many footprints can be found as a 
        function of latitude that are at least 2 HPBW on the planet. 
        Observation that are close to the limb are heavily contaminated by 
        Synchrotron and the Cold sky, causing the limb darkening fit to 
        deteriorate. 

        '''

        lat_interp = np.arange(-88,89)
        eangle_min = np.zeros((len(lat_interp),6))*np.nan 
        for cn in range(1,7):
            for i in range(len(lat_interp)): 
                ind = np.where(eval(f'(self.C{cn}.lat_c > lat_interp[i] - 1)')  & eval(f'(self.C{cn}.lat_c < lat_interp[i] + 1)') )[0]
                if ind.size!=0: 
                    eangle_min[i,cn-1] = np.nanmin(eval(f'(self.C{cn}.eangle[ind])')) 


        fig, axs = plt.subplots(1, 1)
        axs.plot(lat_interp,eangle_min[:,0], label='C1', alpha=0.3)
        axs.plot(lat_interp,eangle_min[:,1], label='C2', alpha=0.3)
        axs.plot(lat_interp,eangle_min[:,2], label='C3', alpha=0.3)
        axs.plot(lat_interp,eangle_min[:,3], label='C4', alpha=0.3)
        axs.plot(lat_interp,eangle_min[:,4], label='C5', alpha=0.3)
        axs.plot(lat_interp,eangle_min[:,5], label='C6', alpha=0.3)
        axs.legend()
        axs.set_title(f'PJ{self.PJnumber} ')
        axs.set_xlabel('Latitude')
        axs.set_ylabel('Minimum emission angle')


        return eangle_min, lat_interp



    def plot_age(self,channel,dlat=1): 
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

        # Version 2
        (np.savez(path2save + '_v2.npz',
        lat=self.C1.lat_c_i, lat_g=self.C1.lat_g_i, lat_c=self.C1.lat_c_i,   
        T = np.vstack([self.C1.T_n, self.C2.T_n, self.C3.T_n, self.C4.T_n, self.C5.T_n, self.C6.T_n]), 
        Tsig = np.vstack([self.C1.sig_T, self.C2.sig_T, self.C3.sig_T, self.C4.sig_T, self.C5.sig_T, self.C6.sig_T]), 
        p = np.vstack([self.C1.p, self.C2.p, self.C3.p, self.C4.p, self.C5.p, self.C6.p]), 
        psig = np.vstack([self.C1.sig_p, self.C2.sig_p, self.C3.sig_p, self.C4.sig_p, self.C5.sig_p, self.C6.sig_p]), 
        w = np.vstack([self.C1.w, self.C2.w, self.C3.w, self.C4.w, self.C5.w, self.C6.w]))) 
        

        # Version 1 
        (np.savez(path2save + '.npz',
        lat=self.C1.lat_c_i, lat_g=self.C1.lat_g_i, lat_c=self.C1.lat_c_i,    
        T1 = self.C1.T_n, T2 = self.C2.T_n, T3 = self.C3.T_n, T4 = self.C4.T_n, T5 = self.C5.T_n, T6 = self.C6.T_n, 
        sig_T1 = self.C1.sig_T, sig_T2 = self.C2.sig_T, sig_T3 = self.C3.sig_T, sig_T4 = self.C4.sig_T, sig_T5 = self.C5.sig_T, sig_T6 = self.C6.sig_T, 
        sig_p1 = self.C1.sig_p, sig_p2 = self.C2.sig_p, sig_p3 = self.C3.sig_p, sig_p4 = self.C4.sig_p, sig_p5 = self.C5.sig_p, sig_p6 = self.C6.sig_p, 
        p1 = self.C1.p, p2 = self.C2.p, p3 = self.C3.p, p4 = self.C4.p, p5 = self.C5.p, p6 = self.C6.p, 
        w1 = self.C1.w, w2 = self.C2.w, w3 = self.C3.w, w4 = self.C4.w, w5 = self.C5.w, w6 = self.C6.p,)) 
        



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


def BeamConvolve(): 
    '''
    Function that calculates the effect of emission angle on the beam 

    '''

    import pyPR.JunoTools as jt
    import matplotlib as mpl

    path_J = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/Beam/'
    beam = 'antenna_normalized2.txt' 
    path    = path_J+beam
    G,p,t   = jt.readJunoBeam(path,plotting=True) 

    #-- Generate Data -----------------------------------------
    # Using linspace so that the endpoint of 360 is included...
    azimuths = np.degrees(t) #theta -> outwards radial component 
    zeniths = np.degrees(p) # phi -> contours around the beam

    r, theta = np.meshgrid(zeniths, azimuths)

    #-- Plot... ------------------------------------------------
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    p1 = ax.contourf(theta, r, np.log(G.T))
    #-- obtaining the colormap limits
    vmin,vmax = p1.get_clim()
    ax.set_rlim(0,50)

    #-- Defining a normalised scale
    cNorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    #-- Creating a new axes at the right side
    ax2 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    #-- Plotting the colormap in the created axes
    cb1 = mpl.colorbar.ColorbarBase(ax2, norm=cNorm)
    fig.subplots_adjust(left=0.05,right=0.85)
    plt.show()

    # Make a map of emission angle and convert to polar 

    # Convert from Planetframe into spacecraft frame 








def ProcessZonalAverage(path, pjmin = 1, pjmax = 9, reprocess=False, LBfilter=True): 
    '''
    
    import pyPR.JunoTools as jt 
    path_GD='/Users/chris/GDrive-UCB/'
    path =  path_GD + 'Berkeley/Research/Juno/'
    jt.ProcessZonalAverage(path, pjmin=1, pjmax=9,reprocess=False)

    '''

    # ----------------------------------------------------------------------
    import glob, os 
    import scipy.signal as signal

    pjtotal = pjmax - pjmin + 1
    fname = f'PJ{pjmin}-{pjmax}'
    try: 
        import os 
        os.system(f'mkdir  {path}{fname}')
    except: 
        print(f'Folder {fname} already exists')

    if not glob.glob(path + fname + f'/PJ{pjmin}-{pjmax}.npz') or reprocess: 
        # If data need to be reproccessed and saved again 
        for pj in range(pjmin,pjmax+1): 
            PJi = PJ(pj)
            try: PJi.readdata(path,Footprints=True)
            except SystemExit: continue 

    # Pre allocate the arrays
    C1 = []; C2 = []; C3 = []; C4 = []; C5 = []; C6 = []
    sig_C1 = []; sig_C2 = []; sig_C3 = []; sig_C4 = []; sig_C5 = []; sig_C6 = []; 
    p1 = []; p2 = []; p3 = []; p4 = []; p5 = []; p6 = []  
    sig_p1 = []; sig_p2 = []; sig_p3 = []; sig_p4 = []; sig_p5 = []; sig_p6 = [] 
    w1 = []; w2 = []; w3 = []; w4 = []; w5 = []; w6 = []
    C1a = np.empty((0,901)); C2a = np.empty((0,901)); C3a = np.empty((0,901)); C4a = np.empty((0,901)); C5a = np.empty((0,901)); C6a = np.empty((0,901))
    p1a = np.empty((0,901)); p2a = np.empty((0,901)); p3a = np.empty((0,901)); p4a = np.empty((0,901)); p5a = np.empty((0,901)); p6a = np.empty((0,901))
    w1a = np.empty((0,901)); w2a = np.empty((0,901)); w3a = np.empty((0,901)); w4a = np.empty((0,901)); w5a = np.empty((0,901)); w6a = np.empty((0,901))


    data = []
    for pj in range(pjmin,pjmax+1): 
        try: 
            print(f'Appending PJ{pj}')
            T = np.load(path+f'PJ{pj}/PJ{pj}.npz')
            data.append(True)

            for cn in range(1,7): 
                # Nadir brightness temperature 
                eval(f'C{cn}.append(T[\'T{cn}\'])')
                # Limb darkening 
                eval(f'p{cn}.append(T[\'p{cn}\'])')
                # Weights of the measurement 
                eval(f'w{cn}.append(T[\'w{cn}\'])')
                # Weights of the measurement 
                eval(f'sig_C{cn}.append(T[\'sig_T{cn}\'])')
                # Weights of the measurement 
                eval(f'sig_p{cn}.append(T[\'sig_p{cn}\'])')
        except: 
            data.append(False)
            for cn in range(1,7): 
                eval(f'C{cn}.append(np.ones_like(T[\'lat\'])*np.nan)')
                eval(f'p{cn}.append(np.ones_like(T[\'lat\'])*np.nan)')
                eval(f'w{cn}.append(np.ones_like(T[\'lat\'])*np.nan)')
                eval(f'sig_C{cn}.append(np.ones_like(T[\'lat\'])*np.nan)')
                eval(f'sig_p{cn}.append(np.ones_like(T[\'lat\'])*np.nan)')

    
    lat = T['lat']
    lat_g = T['lat_g']
    lat_c = T['lat_c']

    # 
    C1 = np.array(C1); C2 = np.array(C2); C3 = np.array(C3); C4 = np.array(C4); C5 = np.array(C5); C6 = np.array(C6)
    p1 = np.array(p1); p2 = np.array(p2); p3 = np.array(p3); p4 = np.array(p4); p5 = np.array(p5); p6 = np.array(p6)   
    w1 = np.array(w1); w2 = np.array(w2); w3 = np.array(w3); w4 = np.array(w4); w5 = np.array(w5); w6 = np.array(w6) 

    for pj in range(0,pjtotal): 
        if data[pj]: 
            C1a = np.vstack([C1a,C1[pj]]); C2a = np.vstack([C2a,C2[pj]]); C3a = np.vstack([C3a,C3[pj]]); C4a = np.vstack([C4a,C4[pj]]); C5a = np.vstack([C5a,C5[pj]]); C6a = np.vstack([C6a,C6[pj]])
            p1a = np.vstack([p1a,p1[pj]]); p2a = np.vstack([p2a,p2[pj]]); p3a = np.vstack([p3a,p3[pj]]); p4a = np.vstack([p4a,p4[pj]]); p5a = np.vstack([p5a,p5[pj]]); p6a = np.vstack([p6a,p6[pj]])
            w1a = np.vstack([w1a,w1[pj]]); w2a = np.vstack([w2a,w2[pj]]); w3a = np.vstack([w3a,w3[pj]]); w4a = np.vstack([w4a,w4[pj]]); w5a = np.vstack([w5a,w5[pj]]); w6a = np.vstack([w6a,w6[pj]])

    # Calculate the statistics for the Temperature 
    C1md = np.nanmedian(C1a,axis=0); C2md = np.nanmedian(C2a,axis=0); C3md = np.nanmedian(C3a,axis=0); C4md = np.nanmedian(C4a,axis=0); C5md = np.nanmedian(C5a,axis=0); C6md = np.nanmedian(C6a,axis=0)
    C1m = np.nanmean(C1a,axis=0); C2m = np.nanmean(C2a,axis=0); C3m = np.nanmean(C3a,axis=0); C4m = np.nanmean(C4a,axis=0); C5m = np.nanmean(C5a,axis=0); C6m = np.nanmean(C6a,axis=0) 
    C1std = np.nanstd(C1a,axis=0); C2std = np.nanstd(C2a,axis=0); C3std = np.nanstd(C3a,axis=0); C4std = np.nanstd(C4a,axis=0); C5std = np.nanstd(C5a,axis=0); C6std = np.nanstd(C6a,axis=0) 

    # Remove NaN's by setting their weight to zero 
    w1a[~np.isfinite(w1a)] = 0; w2a[~np.isfinite(w2a)] = 0; w3a[~np.isfinite(w3a)] = 0; w4a[~np.isfinite(w4a)] = 0; w5a[~np.isfinite(w5a)] = 0; w6a[~np.isfinite(w6a)] = 0 
    C1wa = C1a; C2wa = C2a; C3wa = C3a; C4wa = C4a; C5wa = C5a; C6wa = C6a 
    C1wa[~np.isfinite(C1wa)] = 0; C2wa[~np.isfinite(C2wa)] = 0; C3wa[~np.isfinite(C3wa)] = 0; C4wa[~np.isfinite(C4wa)] = 0; C5wa[~np.isfinite(C5wa)] = 0; C6wa[~np.isfinite(C6wa)] = 0 
    # Create a weighted mean stack based the limb darkening 
    C1wm = np.ma.filled(np.ma.average(C1wa,axis=0,weights=w1a),np.nan); C2wm = np.ma.filled(np.ma.average(C2wa,axis=0,weights=w2a),np.nan); C3wm = np.ma.filled(np.ma.average(C3wa,axis=0,weights=w3a),np.nan); C4wm = np.ma.filled(np.ma.average(C4wa,axis=0,weights=w4a),np.nan); C5wm = np.ma.filled(np.ma.average(C5wa,axis=0,weights=w5a),np.nan); C6wm = np.ma.filled(np.ma.average(C6wa,axis=0,weights=w6a),np.nan) 
    
    # Create a weighted mean using the uncertainty in fitting the limb darkening parameter
    C1wm2 = np.ma.filled(np.ma.average(C1wa,axis=0,weights=(1/np.asarray(sig_C1))[np.where(data)]),np.nan); C2wm2 = np.ma.filled(np.ma.average(C2wa,axis=0,weights=1/np.asarray(sig_C2)[np.where(data)]),np.nan); C3wm2 = np.ma.filled(np.ma.average(C3wa,axis=0,weights=1/np.asarray(sig_C3)[np.where(data)]),np.nan); C4wm2 = np.ma.filled(np.ma.average(C4wa,axis=0,weights=1/np.asarray(sig_C4)[np.where(data)]),np.nan); C5wm2 = np.ma.filled(np.ma.average(C5wa,axis=0,weights=1/np.asarray(sig_C5)[np.where(data)]),np.nan); C6wm2 = np.ma.filled(np.ma.average(C6wa,axis=0,weights=1/np.asarray(sig_C6)[np.where(data)]),np.nan) 

    # Calculate the statistics for the limb darkening  
    p1md = np.nanmedian(p1a,axis=0); p2md = np.nanmedian(p2a,axis=0); p3md = np.nanmedian(p3a,axis=0); p4md = np.nanmedian(p4a,axis=0); p5md = np.nanmedian(p5a,axis=0); p6md = np.nanmedian(p6a,axis=0)
    p1m = np.nanmean(p1a,axis=0); p2m = np.nanmean(p2a,axis=0); p3m = np.nanmean(p3a,axis=0); p4m = np.nanmean(p4a,axis=0); p5m = np.nanmean(p5a,axis=0); p6m = np.nanmean(p6a,axis=0) 
    p1std = np.nanstd(p1a,axis=0); p2std = np.nanstd(p2a,axis=0); p3std = np.nanstd(p3a,axis=0); p4std = np.nanstd(p4a,axis=0); p5std = np.nanstd(p5a,axis=0); p6std = np.nanstd(p6a,axis=0) 

    C1wm_lbpf = np.copy(C1wm)
        # Low bandpass filter (Buterworth filter) for the Channel 1 data to remove wobbles 
    if LBfilter: 
        N  = 2    # Filter order
        Wn = 0.03 # Cutoff frequency
        B, A = signal.butter(N, Wn, output='ba')
        C1wm_lbpf[np.isfinite(C1wm)] = signal.filtfilt(B,A, C1wm[np.isfinite(C1wm)])
   
    # Overwrite! 
    C1wm = C1wm_lbpf

    # Save the variables into one big matrix  
    T = np.zeros((pjmax,6,len(lat))) 
    p = np.zeros((pjmax,6,len(lat))) 

    for pj in range(0,pjmax): 
        if pj < pjmin: 
            print(pj) 
            T[pj,:,:] = np.zeros((6,len(lat)))*np.nan
            p[pj,:,:] = np.zeros((6,len(lat)))*np.nan
        if data[pj-pjmin]: 
            for cn in range(1,7): 
                T[pj,cn-1,:] = eval(f'C{cn}[pj]')
                p[pj,cn-1,:] = eval(f'p{cn}[pj]')  
        else: 
            T[pj,:,:] = np.zeros((6,len(lat)))*np.nan
            p[pj,:,:] = np.zeros((6,len(lat)))*np.nan


    # Obtain relevant data 
    (np.savez(path + fname + f'/PJ{pjmin}-{pjmax}_v2.npz',
                data = data, lat     = lat, lat_g = lat_g, lat_c = lat_c, 
                
                # Brightness temperature 
                T = T, 
                Tmean = np.vstack([C1m, C2m, C3m, C4m, C5m, C6m]), 
                Tstd = np.vstack([C1std, C2std, C3std, C4std, C5std, C6std]), 
                Tmedian = np.vstack([C1md , C2md , C3md , C4md , C5md ,C6md ]),  
                Twm = np.vstack([C1wm , C2wm , C3wm , C4wm , C5wm ,C6wm ]),
                Twm2 = np.vstack([C1wm2 , C2wm2 , C3wm2 , C4wm2 , C5wm2 , C6wm2 ]),
                Tsig = np.vstack([sig_C1 , sig_C2 , sig_C3 , sig_C4 , sig_C5 , sig_C6 ]),
                
                # Limb darkening coefficient 
                p = p,
                pmean = np.vstack([p1m, p2m, p3m, p4m, p5m, p6m]), 
                pstd = np.vstack([p1std, p2std, p3std, p4std, p5std, p6std]), 
                pmedian = np.vstack([p1md , p2md , p3md , p4md , p5md ,p6md ]),  
                psig = np.vstack([sig_p1 , sig_p2 , sig_p3 , sig_p4 , sig_p5 , sig_p6 ]))) 



    # Save the variables 
    (np.savez(path + fname + f'/PJ{pjmin}-{pjmax}.npz',
                data = data, lat     = lat, lat_g = lat_g, lat_c = lat_c, 
                # Channel 1
                C1      = C1, C1std   = C1std, C1mean  = C1m, C1median= C1md, C1wm    = C1wm, C1wm2 = C1wm2,  C1sig = sig_C1, 
                p1std   = p1std, p1      = p1, p1mean  = p1m, p1median= p1md, p1sig = sig_p1, 
                # Channel 2
                C2      = C2, C2std   = C2std, C2mean  = C2m, C2median= C2md, C2wm    = C2wm, C2wm2 = C2wm2,  C2sig = sig_C2, 
                p2      = p2, p2std   = p2std,  p2mean  = p2m, p2median= p2md,  p2sig = sig_p2, 
                # Channel 3
                C3      = C3, C3std   = C3std, C3mean  = C3m, C3median= C3md, C3wm    = C3wm, C3wm2 = C3wm2,  C3sig = sig_C3, 
                p3      = p3, p3std   = p3std, p3mean  = p3m, p3median= p3md,  p3sig = sig_p3, 
                # Channel 4
                C4      = C4, C4std   = C4std, C4mean  = C4m, C4median= C4md, C4wm    = C4wm, C4wm2 = C4wm2,  C4sig = sig_C4, 
                p4      = p4, p4std   = p4std, p4mean  = p4m, p4median= p4md,  p4sig = sig_p4, 
                # Channel 5
                C5      = C5, C5std   = C5std, C5mean  = C5m, C5median= C5md, C5wm    = C5wm, C5wm2 = C5wm2,  C5sig = sig_C5, 
                p5      = p5, p5std   = p5std,  p5mean  = p5m, p5median= p5md,  p5sig = sig_p5, 
                # Channel 6
                C6      = C6, C6std   = C6std, C6mean  = C6m, C6median= C6md, C6wm    = C6wm, C6wm2 = C6wm2,  C6sig = sig_C5, 
                p6      = p6, p6std   = p6std, p6mean  = p6m, p6median= p6md,  p6sig = sig_p6, ))
 

    return path + fname + f'/PJ{pjmin}-{pjmax}.npz'

def PlotZonalAverage(path,pjmin=1, pjmax = 9, savefig = False): 
    '''
    pjmin = 1 
    pjmax = 9

    path_J = path_GD + 'Berkeley/Research/Juno/'
    path = path_J +  f'PJ{pjmin}-{pjmax}' + f'/PJ{pjmin}-{pjmax}_v2.npz'
    jt.PlotZonalAverage(path,pjmin = pjmin, pjmax = pjmax)
    '''

    pjtotal = pjmax - pjmin + 1 
    temp    = np.load(path,allow_pickle=True)

    lat     = temp['lat']
    # Load in the background 
    T   = temp['T']; Tm  = temp['Tmean']; Tstd= temp['Tstd']; Tmd = temp['Tmedian']; Twm = temp['Twm']; Twm2= temp['Twm2'] 
    pm  = temp['pmean']; pstd= temp['pstd']; pmd = temp['pmedian']

    labels = ['0.6GHz','1.2GHz','2.5GHz','5GHz','10GHz', '22GHz']

    # Mean nadir brightness temperature for various 
    # ----------------------------------------------------------------------
    fig, axs = plt.subplots(1, 1,figsize=(16,9))
    for cn in range(0,6): 
        axs.plot(lat,Tm[cn,:], label=labels[cn], alpha=0.9 ,linestyle='-', color=cmap(cn/6),)
        axs.plot(lat,Twm[cn,:], alpha=0.9 ,linestyle='--', color=cmap(cn/6),)
        axs.plot(lat,Twm2[cn,:], alpha=0.9 ,linestyle=':', color=cmap(cn/6),)
        axs.plot(lat,Tmd[cn,:], alpha=0.9 ,linestyle='-.', color=cmap(cn/6),)

    axs.set_xlabel('Latitude [deg]')
    axs.set_ylabel('Nadir temperature [K]')
    axs.invert_xaxis()
    axs.set_ylim([950,100])
    axs.set_xlim([-40,40])
   
    linelabels = ["mean", "ldwm", "sigwm","md"]

    lines = axs.get_lines()
    legend1 = plt.legend([lines[i] for i in [0,4,8,12,16,20]],labels, loc=7,ncol=2)
    legend2 = plt.legend([lines[i] for i in [0,1,2,3]], linelabels, loc=6)
    axs.add_artist(legend1)
    axs.add_artist(legend2)

    if savefig: 
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_meanoverview.png', format='png', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_meanoverview.pdf', format='pdf', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_meanoverview.eps', format='eps', transparent = True, dpi=500)

    
    # Plot the ld weighted mean of the first orbits 
    # ----------------------------------------------------------------------
    fig, axs = plt.subplots(1, 1,figsize=(16,9))
    for cn in range(0,6): 
        axs.plot(lat,Twm[cn,:], label=labels[cn], alpha=0.9 ,linestyle='-', color=cmap(cn/6),)
        axs.fill_between(lat,Twm[cn,:]+Tstd[cn,:],Twm[cn,:]-Tstd[cn,:], alpha=0.3, color=cmap(cn/6))

    axs.set_xlabel('Latitude [deg]')
    axs.set_ylabel('Nadir temperature [K]')
    axs.invert_xaxis()
    axs.set_ylim([950,100])
    axs.set_xlim([-40,40])
    axs.legend(loc=7,ncol=2) 

    if savefig: 
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_mean.png', format='png', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_mean.pdf', format='pdf', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_mean.eps', format='eps', transparent = True, dpi=500)



    # Plot only the mean of the first orbits 
    # ----------------------------------------------------------------------
    fig, axs = plt.subplots(1, 1,figsize=(16,9))
    for cn in range(2,6):
        axs.plot(Twm[cn,:],lat,  label=labels[cn], alpha=0.9, color=cmap(cn/6), linewidth=3)
        axs.fill_betweenx(lat,Twm[cn,:]+Tstd[cn,:],Twm[cn,:]-Tstd[cn,:], alpha=0.3, color=cmap(cn/6))

    axs.set_ylabel('Latitude [deg]')
    axs.set_xlabel('Nadir temperature [K]')
    axs.invert_xaxis()
    axs.set_xlim([350,100])
    axs.set_ylim([-40,40])

    #plt.legend(loc='center right')
    if savefig: 
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_mean.png', format='png', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_mean.pdf', format='pdf', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_mean.eps', format='eps', transparent = True, dpi=500)



    # Make a polar plot 
    fig = plt.figure(figsize=(9,16))
    ax = fig.add_subplot(111, polar=True)
    ax.set_rlim(bottom=350, top=100)
    ax.set_rorigin(700)
    ax.set_rlabel_position(45)
    plt.text(0.3, .84, 'Nadir temperature [K]', transform=ax.transAxes,rotation=45,fontsize=11)
    ax.set_thetamin(-45)
    ax.set_thetamax(45)
    for cn in range(2,6):
        ax.plot(np.radians(lat),Twm[cn,:], label=labels[cn], color=cmap(cn/6), alpha=0.75, linewidth= 2 )
        ax.fill_between(np.radians(lat),Twm[cn,:]+Tstd[cn,:],Twm[cn,:]-Tstd[cn,:], alpha=0.3, color=cmap(cn/6))

    plt.show()

    if savefig: 
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_polar_mean.png', format='png', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_polar_mean.pdf', format='pdf', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno{pjmin}-{pjmax}_zonal_polar_mean.eps', format='eps', transparent = True, dpi=500)



    # Plot mean and the variations of pjtotal orbits 
    # ----------------------------------------------------------------------
    fig, axs = plt.subplots(1, 1,figsize=(16,9)) 

    for pj in range(pjmin,pjmax): 
        try:   
            for cn in range(2,6):  
                axs.plot(C1[pj],lat,  label=f'PJ{pj+1}', alpha=0.5, color=cmap3(pj/pjtotal))
        except: pass 

    axs.plot(C1m,lat,  label='0.6GHz', alpha=0.9 ,  color=cc1, linewidth=3)
    axs.plot(C2m,lat,  label='1.2GHz', alpha=0.9 ,  color=cc2, linewidth=3)
    axs.plot(C3m,lat,  label='2.5GHz', alpha=0.9 ,  color=cc3, linewidth=3)
    axs.plot(C4m,lat,  label='5GHz', alpha=0.9 ,    color=cc4, linewidth=3)
    axs.plot(C5m,lat,  label='10GHz', alpha=0.9 ,   color=cc5, linewidth=3)
    axs.plot(C6m,lat,  label='22GHz', alpha=0.9 ,   color=cc6, linewidth=3)
    
    axs.plot(C1md, lat, alpha=0.9 ,linestyle='-.', color=cc1,)
    axs.plot(C2md, lat, alpha=0.9 ,linestyle='-.', color=cc2,)
    axs.plot(C3md, lat, alpha=0.9 ,linestyle='-.', color=cc3,)
    axs.plot(C4md, lat, alpha=0.9 ,linestyle='-.', color=cc4,)
    axs.plot(C5md, lat, alpha=0.9 ,linestyle='-.', color=cc5,)
    axs.plot(C6md, lat, alpha=0.9 ,linestyle='-.', color=cc6,)

    axs.plot(C1wm, lat, alpha=0.9 ,linestyle='--', color=cc1,)
    axs.plot(C2wm, lat, alpha=0.9 ,linestyle='--', color=cc2,)
    axs.plot(C3wm, lat, alpha=0.9 ,linestyle='--', color=cc3,)
    axs.plot(C4wm, lat, alpha=0.9 ,linestyle='--', color=cc4,)
    axs.plot(C5wm, lat, alpha=0.9 ,linestyle='--', color=cc5,)
    axs.plot(C6wm, lat, alpha=0.9 ,linestyle='--', color=cc6,)

    axs.plot(C1wm2, lat, alpha=0.9 ,linestyle=':', color=cc1,)
    axs.plot(C2wm2, lat, alpha=0.9 ,linestyle=':', color=cc2,)
    axs.plot(C3wm2, lat, alpha=0.9 ,linestyle=':', color=cc3,)
    axs.plot(C4wm2, lat, alpha=0.9 ,linestyle=':', color=cc4,)
    axs.plot(C5wm2, lat, alpha=0.9 ,linestyle=':', color=cc5,)
    axs.plot(C6wm2, lat, alpha=0.9 ,linestyle=':', color=cc6,)

    axs.set_ylabel('Latitude [deg]')
    axs.set_xlabel('Nadir temperature [K]')
    axs.invert_xaxis()
    axs.set_xlim([950,100])
    axs.set_ylim([-40,40])
    #plt.legend(loc='center right')
    if savefig: 
        plt.savefig(path2save + f'Juno_zonal_mean_{pjmax}.png', format='png', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno_zonal_mean_{pjmax}.pdf', format='pdf', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno_zonal_mean_{pjmax}.eps', format='eps', transparent = True, dpi=500)





    # Make a polar plot 
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111, polar=True)
    ax.set_rlim(bottom=350, top=100)
    ax.set_rorigin(700)
    ax.set_rlabel_position(45)
    ax.set_thetamin(-45)
    ax.set_thetamax(45)
    plt.text(0.3, .84, 'Nadir temperature [K]', transform=ax.transAxes,rotation=45,fontsize=11)

    for pj in range(pjmin,pjmax): 
        try:  
            #ax.plot(np.radians(lat), C1[pj],  alpha=0.5, color=cmap(pj/pjtotal))
            #ax.plot(np.radians(lat), C2[pj],  alpha=0.5, color=cmap(pj/pjtotal))
            ax.plot(np.radians(lat), C3[pj],  alpha=0.5, color=cmap3(pj/pjtotal))
            ax.plot(np.radians(lat), C4[pj],  alpha=0.5, color=cmap3(pj/pjtotal))
            ax.plot(np.radians(lat), C5[pj],  alpha=0.5, color=cmap3(pj/pjtotal))
            ax.plot(np.radians(lat), C6[pj],  alpha=0.5, color=cmap3(pj/pjtotal))
        except: 
            pass 
    #ax.plot(np.radians(lat), C1m, label='', alpha=0.9 , color='navy', linewidth=3)
    #ax.plot(np.radians(lat), C2m, label='', alpha=0.9 , color='navy', linewidth=3)
    ax.plot(np.radians(lat), C3m, label='2.5GHz', alpha=0.9, color=cc3, linewidth=2)
    ax.plot(np.radians(lat), C4m, label='5GHz', alpha=0.9 ,  color=cc4,  linewidth=2)
    ax.plot(np.radians(lat), C5m, label='10GHz', alpha=0.9 , color=cc5, linewidth=2)
    ax.plot(np.radians(lat), C6m, label='22GHz', alpha=0.9 , color=cc6, linewidth=2)


    #plt.legend(loc='center right')
    if savefig: 
        plt.savefig(path2save + f'Juno_zonal_polar_mean_{pjmax}.png', format='png', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno_zonal_polar_mean_{pjmax}.pdf', format='pdf', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno_zonal_polar_mean_{pjmax}.eps', format='eps', transparent = True, dpi=500)



    # Plot the limb darkening parameter 
    # Plot the mean between -30 and 30 deg 
    mean_l = 300 
    mean_h = 600
    fig, axs = plt.subplots(1, 1,figsize=(16,9))
    axs.plot(lat,p1m,  label='0.6GHz', alpha=0.9 ,       color=cc1, linewidth=3)
    axs.fill_between(lat,p1m+p1std,p1m-p1std, alpha=0.3, color=cc1)
    axs.plot(lat,p2m,  label='1.2GHz', alpha=0.9 ,       color=cc2, linewidth=3)
    axs.fill_between(lat,p2m+p2std,p2m-p2std, alpha=0.3, color=cc2)
    axs.plot(lat,p3m,  label='2.5GHz', alpha=0.9 ,       color=cc3, linewidth=3)
    axs.fill_between(lat,p3m+p3std,p3m-p3std, alpha=0.3, color=cc3)
    axs.plot(lat,p4m,  label='5GHz', alpha=0.9 ,         color=cc4, linewidth=3)
    axs.fill_between(lat,p4m+p4std,p4m-p4std, alpha=0.3, color=cc4)
    axs.plot(lat,p5m,  label='10GHz', alpha=0.9 ,        color=cc5, linewidth=3)
    axs.fill_between(lat,p5m+p5std,p5m-p5std, alpha=0.3, color=cc5)
    axs.plot(lat,p6m,  label='22GHz', alpha=0.9 ,        color=cc6, linewidth=3)
    axs.fill_between(lat,p6m+p6std,p6m-p6std, alpha=0.3, color=cc6)

    axs.text(-20,0.57,f'C1: {np.mean(p1m[300:600]):2.2f}')
    axs.text(-20,0.53,f'C2: {np.mean(p2m[300:600]):2.2f}')
    axs.text(-20,0.49,f'C3: {np.mean(p3m[300:600]):2.2f}')
    axs.text(10,0.57,f'C4: {np.mean(p4m[300:600]):2.2f}')
    axs.text(10,0.53,f'C5: {np.mean(p5m[300:600]):2.2f}')
    axs.text(10,0.49,f'C6: {np.mean(p6m[300:600]):2.2f}')

    axs.set_xlabel('Latitude [deg]')
    axs.set_ylabel('Limb darkening [-]')
    axs.invert_xaxis()
    axs.set_ylim([0,0.6])
    axs.set_xlim([-60,60])

    if savefig: 
        plt.savefig(path2save + f'Juno_ld_zonal_mean_{pjmax}.png', format='png', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno_ld_zonal_mean_{pjmax}.pdf', format='pdf', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno_ld_zonal_mean_{pjmax}.eps', format='eps', transparent = True, dpi=500)




    fig, axs = plt.subplots(1, 1,figsize=(16,9))
    for pj in range(pjmin,pjmax): 
        try:  
            axs.plot(lat,p1[pj],  label=f'PJ{pj+1}', alpha=0.5, color=cmap3(pj/pjtotal))
            axs.plot(lat,p2[pj],  alpha=0.5, color=cmap3(pj/pjtotal))
            axs.plot(lat,p3[pj],  alpha=0.5, color=cmap3(pj/pjtotal))
            axs.plot(lat,p4[pj],  alpha=0.5, color=cmap3(pj/pjtotal))
            axs.plot(lat,p5[pj],  alpha=0.5, color=cmap3(pj/pjtotal))
            axs.plot(lat,p6[pj],  alpha=0.5, color=cmap3(pj/pjtotal))
        except: 
            pass 
    axs.plot(lat,p1m,  label='0.6GHz', alpha=0.9 , color=cc1, linewidth=3)
    axs.plot(lat,p2m,  label='1.2GHz', alpha=0.9 , color=cc2, linewidth=3)
    axs.plot(lat,p3m,  label='2.5GHz', alpha=0.9 , color=cc3, linewidth=3)
    axs.plot(lat,p4m,  label='5GHz', alpha=0.9 ,   color=cc4, linewidth=3)
    axs.plot(lat,p5m,  label='10GHz', alpha=0.9 ,  color=cc5, linewidth=3)
    axs.plot(lat,p6m,  label='22GHz', alpha=0.9 ,  color=cc6, linewidth=3)

    axs.set_xlabel('Latitude [deg]')
    axs.set_ylabel('Limb darkening [-]')
    axs.invert_xaxis()
    axs.set_ylim([0,0.8])
    axs.set_xlim([-60,60])

    if savefig: 
        plt.savefig(path2save + f'Juno_ld_zonal_mean_PJ{pj}.png', format='png', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno_ld_zonal_mean_PJ{pj}.pdf', format='pdf', transparent = True, dpi=500)
        plt.savefig(path2save + f'Juno_ld_zonal_mean_PJ{pj}.eps', format='eps', transparent = True, dpi=500)

    return 

def PlotPJ(PJ,xlim=[-40,40]): 
    '''
    Plot the most important graphics for a given PJ and compare them to the orbit mean 
        
    import pyPR.JunoTools as jt 
    PJnumber = 9
    PJ = jt.PlotPJ(PJnumber)

    '''
    labels = ['0.6GHz','1.2GHz','2.5GHz','5GHz','10GHz', '22GHz']

    # Load the observed brightness temperatures  
    # ------------------------------------------------------------------
    PJbg = '1-9' 
    path2load = path_J +f'PJ{PJbg}/' 
    temp = np.load(path2load + f'PJ{PJbg}_v2.npz')
    lat     = temp['lat_c']
    Tm  = temp['Tmean']; Tstd= temp['Tstd']; Tmd = temp['Tmedian']; Twm = temp['Twm']; Twm2= temp['Twm2'] 
    pm  = temp['pmean']; pstd= temp['pstd']; pmd = temp['pmedian']

    # Read in specific PJ data 
    # ------------------------------------------------------------------
    path2load = path_J +f'PJ{PJ}/' 
    temp = np.load(path2load + f'PJ{PJ}_v2.npz')
    T = temp['T']; Tsig = temp['Tsig']; p = temp['p']; psig = temp['psig']; w = temp['w']   

    # Plot the brightness temperature 
    # ----------------------------------------------------------------------
    fig, axs = plt.subplots(1, 1,figsize=(16,9))
    for cn in range(0,6): 
        axs.plot(lat, T[cn,:],  label=labels[cn],    alpha=0.9 , color='firebrick', linewidth=3)
        axs.fill_between(lat, T[cn,:]+Tsig[cn,:],T[cn,:]-Tsig[cn,:], color='firebrick', alpha=0.3)
        axs.plot(lat, Tm[cn,:],                     alpha=0.9 , color=cmap(cn/6), linewidth=3)
        axs.plot(lat, Twm[cn,:],                     alpha=0.9 , color=cmap(cn/6), linewidth=3, linestyle = ':')
        axs.fill_between(lat, Twm[cn,:]+Tstd[cn,:], Twm[cn,:]-Tstd[cn,:], color=cmap(cn/6), alpha=0.3)

    axs.set_xlabel('Latitude [deg]')
    axs.set_ylabel('Nadir temperature [K]')
    axs.invert_yaxis()
    axs.set_xlim(xlim)




    # Plot the limb darkening  
    # ----------------------------------------------------------------------
    fig, axs = plt.subplots(1, 1,figsize=(16,9))
    for cn in range(0,6): 
        axs.plot(lat, p[cn,:],  label=labels[cn],    alpha=0.9 , color='firebrick', linewidth=3)
        axs.fill_between(lat, p[cn,:]+psig[cn,:],p[cn,:]-psig[cn,:], color='firebrick', alpha=0.3)
        axs.plot(lat, pm[cn,:],  alpha=0.9 , color=cmap(cn/6), linewidth=3)
        axs.fill_between(lat, pm[cn,:]+pstd[cn,:], pm[cn,:]-pstd[cn,:], color=cmap(cn/6), alpha=0.3)

    axs.set_xlabel('Latitude [deg]')
    axs.set_ylabel('Limb darkening coefficient [K]')
    axs.invert_yaxis()
    axs.set_xlim(xlim)

    # axs.plot(lat_j, T2m,  label='0.6GHz',   alpha=0.9 , color=cc2, linewidth=3)
    # axs.plot(lat_j, T2m,  label='1.25GHz',   alpha=0.9 , color=cc2, linewidth=3)
    # axs.plot(lat_j, T3m,  label='2.5GHz',    alpha=0.9 , color=cc3, linewidth=3)
    # axs.plot(lat_j, T4m,  label='5GHz',      alpha=0.9 , color=cc4, linewidth=3)
    # axs.plot(lat_j, T5m,  label='10GHz',     alpha=0.9 , color=cc5, linewidth=3)
    # axs.plot(lat_j, T6m,  label='22GHz',     alpha=0.9 , color=cc6, linewidth=3)

    # axs.fill_between(lat_j, T1m+T1s, T1m-T1s, color=cc1, alpha=0.3)
    # axs.fill_between(lat_j, T2m+T2s, T2m-T2s, color=cc2, alpha=0.3)
    # axs.fill_between(lat_j, T3m+T3s, T3m-T3s, color=cc3, alpha=0.3)
    # axs.fill_between(lat_j, T4m+T4s, T4m-T4s, color=cc4, alpha=0.3)
    # axs.fill_between(lat_j, T5m+T5s, T5m-T5s, color=cc5, alpha=0.3)
    # axs.fill_between(lat_j, T6m+T6s, T6m-T6s, color=cc6, alpha=0.3)

    # axs.plot(lat_j,T1w,  label='0.6GHz',alpha=0.9 , color=cc1, linestyle = ':', linewidth=3)
    # axs.plot(lat_j,T2w,  label='1.2GHz',alpha=0.9 , color=cc2, linestyle = ':', linewidth=3)
    # axs.plot(lat_j,T3w,  label='2.5GHz',alpha=0.9 , color=cc3, linestyle = ':', linewidth=3)
    # axs.plot(lat_j,T4w,  label='5GHz',  alpha=0.9 , color=cc4, linestyle = ':', linewidth=3)
    # axs.plot(lat_j,T5w,  label='10GHz', alpha=0.9 , color=cc5, linestyle = ':', linewidth=3)
    # axs.plot(lat_j,T6w,  label='22GHz', alpha=0.9 , color=cc6, linestyle = ':', linewidth=3)
 
    # axs.plot(lat_j,T1,   alpha=0.9 , color='firebrick', linewidth=3)
    # axs.plot(lat_j,T2,   alpha=0.9 , color='firebrick', linewidth=3)
    # axs.plot(lat_j,T3,   alpha=0.9 , color='firebrick', linewidth=3)
    # axs.plot(lat_j,T4,   alpha=0.9 , color='firebrick', linewidth=3)
    # axs.plot(lat_j,T5,   alpha=0.9 , color='firebrick', linewidth=3)
    # axs.plot(lat_j,T6,   alpha=0.9 , color='firebrick', linewidth=3)
    
    # axs.fill_between(lat_j, T1m+T1sig, T1-T1sig, color='firebrick', alpha=0.3)
    # axs.fill_between(lat_j, T2m+T2sig, T2-T2sig, color='firebrick', alpha=0.3)
    # axs.fill_between(lat_j, T3m+T3sig, T3-T3sig, color='firebrick', alpha=0.3)
    # axs.fill_between(lat_j, T4m+T4sig, T4-T4sig, color='firebrick', alpha=0.3)
    # axs.fill_between(lat_j, T5m+T5sig, T5-T5sig, color='firebrick', alpha=0.3)
    # axs.fill_between(lat_j, T6m+T6sig, T6-T6sig, color='firebrick', alpha=0.3)


    # axs.set_xlabel('Latitude [deg]')
    # axs.set_ylabel('Nadir temperature [K]')
    # axs.invert_yaxis()
    # axs.set_xlim(xlim)


    # # Limb darkening coefficient 
    # # ----------------------------------------------------------------------
    # fig, axs = plt.subplots(1, 1,figsize=(16,9))
    # axs.plot(lat_j,p1m,  label='0.6GHz', alpha=0.9 , color='navy', linewidth=3)
    # axs.plot(lat_j,p2m,  label='1.25GHz', alpha=0.9 , color='navy', linewidth=3)
    # axs.plot(lat_j,p3m,  label='2.5GHz', alpha=0.9 , color='navy', linewidth=3)
    # axs.plot(lat_j,p4m,  label='5GHz', alpha=0.9 , color='navy', linewidth=3)
    # axs.plot(lat_j,p5m,  label='10GHz', alpha=0.9 , color='navy', linewidth=3)
    # axs.plot(lat_j,p6m,  label='22GHz', alpha=0.9 , color='navy', linewidth=3)

    # axs.fill_between(lat_j, p1m+p1s, p1m-p1s, alpha=0.3)
    # axs.fill_between(lat_j, p2m+p2s, p2m-p2s, alpha=0.3)
    # axs.fill_between(lat_j, p3m+p3s, p3m-p3s, alpha=0.3)
    # axs.fill_between(lat_j, p4m+p4s, p4m-p4s, alpha=0.3)
    # axs.fill_between(lat_j, p5m+p5s, p5m-p5s, alpha=0.3)
    # axs.fill_between(lat_j, p6m+p6s, p6m-p6s, alpha=0.3)

    # axs.plot(lat_j,p1,   alpha=0.9 , color='firebrick', linewidth=3)
    # axs.plot(lat_j,p2,   alpha=0.9 , color='firebrick', linewidth=3)
    # axs.plot(lat_j,p3,   alpha=0.9 , color='firebrick', linewidth=3)
    # axs.plot(lat_j,p4,   alpha=0.9 , color='firebrick', linewidth=3)
    # axs.plot(lat_j,p5,   alpha=0.9 , color='firebrick', linewidth=3)
    # axs.plot(lat_j,p6,   alpha=0.9 , color='firebrick', linewidth=3)

    # axs.fill_between(lat_j, p1+p1sig, p1-p1sig, color='firebrick', alpha=0.3)
    # axs.fill_between(lat_j, p2+p2sig, p2-p2sig, color='firebrick', alpha=0.3)
    # axs.fill_between(lat_j, p3+p3sig, p3-p3sig, color='firebrick', alpha=0.3)
    # axs.fill_between(lat_j, p4+p4sig, p4-p4sig, color='firebrick', alpha=0.3)
    # axs.fill_between(lat_j, p5+p5sig, p5-p5sig, color='firebrick', alpha=0.3)
    # axs.fill_between(lat_j, p6+p6sig, p6-p6sig, color='firebrick', alpha=0.3)

    # axs.set_xlabel('Latitude [deg]')
    # axs.set_ylabel('Limb darkening coefficient [K]')
    # axs.invert_yaxis()
    # axs.set_xlim(xlim)





    return 

def PlotRays_flat(PJ, channel, lat_s = [10,15], eanglelim = [30,45], xlim = [], ylim = [] ): 
    '''
    import pyPR.JunoTools as jt 
    PJnumber = 1
    PJ = jt.PJ(PJnumber)
    pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
    PJ.readdata(pathJ) 
    lat = [0,3]
    jt.PlotRays_flat(PJ,1,lat_s = lat,xlim = [-25000, 25000],ylim=[-300,10000], eanglelim = [35,40])
    '''

    # Load Juno perijove three data 

    # PJ = jt.PJ(PJnumber)
    # pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
    # PJ.readdata(pathJ) 


    import matplotlib 


    # Depth probed by ray 

    if channel == 1: 
        zdepth = -300 
    elif channel == 2:
        zdepth = -200 
    elif channel == 3:
        zdepth = -120
    elif channel == 4:
        zdepth = -60
    elif channel == 5:
        zdepth = -30 
    elif channel == 6:
        zdepth = 10

    


    # Planet's radius 
    R = PJ.RJ

    # Range of data to plot, excluding data outside of perijove pass 
    ja = np.where(PJ.pj_ids == True)[0][0]  + np.argmin(np.abs(PJ.ob_lat_g[PJ.pj_ids] - lat_s[0])) 
    jb = np.where(PJ.pj_ids == True)[0][0]  + np.argmin(np.abs(PJ.ob_lat_g[PJ.pj_ids] - lat_s[1]))

    j_data = np.sort([ja,jb]) 

    print(j_data)



    # Indices corresponding to the chosen latitudes 
    j = np.linspace(j_data[0],j_data[1],j_data[1]-j_data[0] ,dtype=np.int32) 

    # Emission angle  
    ea = eval(f'PJ.C{channel}.eangle[j]')  

    # Filter out all rays larger than eanglelim 
    j_f = j[np.where(ea>eanglelim[0]) and np.where(ea<eanglelim[1])]       

    print(j_f) 

    # latitude of the spacecraft 
    Slat = PJ.ob_lat_c[j_f]
    # Latitude to ground track km
    gS = np.radians(PJ.ob_lat_c[j_f])*R[2]

    # Beam intercept latitude 
    Blat = eval(f'PJ.C{channel}.lat_c[j_f]') 

    # Distance from the equator 
    g = np.radians(eval(f'PJ.C{channel}.lat_c[j_f]'))*R[2]
    # Height of the spacecraft 
    h = PJ.range[j_f]-R[0]


    # Viewing angle, off spacecraft nadir 
    va = np.degrees(np.arctan(h/np.abs(g-gS)))



    # Show the location of the Zones and Belts, 
    NTrZ = np.array([17,21],   )
    NEB  = np.array([7,17],    )
    EZ   = np.array([-7,7],    )
    SEB  = np.array([-7,-17], )
    STrZ = np.array([-17,-21],)



    # Converting into local frame 
    S =np.radians([NTrZ,NEB,EZ,SEB,STrZ])*R[2] 
    # Key what clouds are present 
    SC = np.array([[1,1,1],[1,1,0],[1,1,1],[1,1,0],[1,1,1]])

    # Clouds 
    nh3i = np.array([10,44])
    nh4shs = np.array([-32,-6,]) 
    h2os =  np.array([-60,-25,]) 

    # Create list for all the error patches
    nh3     = []
    nh4sh   = []
    h2o     = []

    # Ammonia Clouds 
    for zone,clouds in zip(S,SC):
        if clouds[2]==1:
            rect = matplotlib.patches.Rectangle((zone[0], nh3i[0]), zone[1]-zone[0],nh3i[1]-nh3i[0])
            nh3.append(rect)

    for zone,clouds in zip(S,SC):
        if clouds[1]==1:
            rect = matplotlib.patches.Rectangle((zone[0], nh4shs[0]), zone[1]-zone[0],nh4shs[1]-nh4shs[0])
            nh4sh.append(rect)

    for zone,clouds in zip(S,SC):
        if clouds[0]==1:
            rect = matplotlib.patches.Rectangle((zone[0], h2os[0]), zone[1]-zone[0],h2os[1]-h2os[0])
            h2o.append(rect)

    alpha = 0.3
    # Create patch collection with specified colour/alpha
    pc1 = matplotlib.collections.PatchCollection(nh3, facecolor='green', alpha=alpha,)
    # Create patch collection with specified colour/alpha
    pc2 = matplotlib.collections.PatchCollection(nh4sh,facecolor='yellow', alpha=alpha,)
    # Create patch collection with specified colour/alpha
    pc3 = matplotlib.collections.PatchCollection(h2o, facecolor='blue', alpha=alpha,)


    # Determine xlim and ylim 
    if not xlim:
        xlim =  [np.min(gS),np.max(gS)] 
    if not ylim:
        ylim =  [zdepth,50]  
        
    



    fig,ax = plt.subplots(1,figsize=(15,6))
    # Plot first the clouds 
    # Add collection to axes
    ax.add_collection(pc1)
    ax.add_collection(pc2)
    ax.add_collection(pc3)

    plt.plot(gS,h)
    for i in range(len(va)):
        plt.plot(gS[i]+np.linspace(0,h[i]+500,10)*np.tan(np.degrees(va[i])),h[i]-np.linspace(-0,h[i]+500,10),alpha=0.2,color=cmap((90-ea[i])/90))

    # Plot the latitudes 
    for i in range(0,30): 
        plt.plot([np.radians(i)*R[2],np.radians(i)*R[2]]  ,ylim, linestyle=':',alpha=0.3,color='gray')
        plt.plot([np.radians(-i)*R[2],np.radians(-i)*R[2]],ylim, linestyle=':',alpha=0.3,color='gray')


    plt.ylabel('Altitude (km)')
    plt.xlabel('Ground track (km)')
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.show()

    return 

def Raydlat(PJ, channel, lat_range, zdepth = [], plotting=False): 
    '''
    import pyPR.JunoTools as jt 
    PJnumber = 3
    PJ = jt.PJ(PJnumber)
    pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
    PJ.readdata(pathJ) 
    lat = [-45,45]


    fig,ax = plt.subplots(1,figsize=(16,9))
    ax2 = ax.twinx()
    ax.set_ylabel('Maximum latitude coverage (deg)')
    ax2.set_ylabel('Maximum emission angle (deg)')
    ax.set_title(f'PJ{PJ.PJnumber} C1 - C6')
    ax.set_xlabel(r'Latitude$_c$ (deg)') 
    ax2.invert_yaxis()
    for channel in range(1,7):
        ea_max, dray_max, dbeam_max, lat_array = jt.Raydlat(PJ,channel,lat,plotting=False)
        ax.plot(lat_array,np.abs(dray_max),color=jt.cmap((channel-1)/6),label= f'C{channel}')
        ax2.plot(lat_array,ea_max,color=jt.cmap((channel-1)/6),linestyle='-.')
    ax.legend()
    plt.show()
    

    '''

    # Retrieval has a function to compute this 
    # rr.ProbeDepth(freq, b=[0.0,0.0], thres = [0.05 , 0.5, 1, 0.5 , 0.05 ])


    if channel == 1  and not zdepth: 
        zdepth = [-16, -80, -160, -310, -920 ]
    elif channel == 2 and not zdepth:
        zdepth = [-2, -60 ,-120 , -190, -324 ]
    elif channel == 3 and not zdepth:
        zdepth = [7, -40 , -70 , -110, -160]
    elif channel == 4 and not zdepth:
        zdepth = [11, -10 , -40 , -70, -95]
    elif channel == 5 and not zdepth:
        zdepth = [14, 7, -8 , -30, - 48]
    elif channel == 6 and not zdepth:
        zdepth = [30, 20, 12 , 8, 6]

    # Find the relevant rays to plot 
    # ------------------------------------------------------------------
    # Find all the data points corresponding to the given latitude range 

    # Find all the beam intercept points that are wihin the main PJ path (pj_ids) and that fall within lat_range
    j_alldata = np.where( PJ.pj_ids)[0][0] + np.argwhere((eval(f'PJ.C{channel}.lat_c[PJ.pj_ids]') > lat_range[0]) & (eval(f'PJ.C{channel}.lat_c[PJ.pj_ids]') < lat_range[1]+1)) 
    
    # Filter out data where the main beam is not on the planet 
    j_data = j_alldata[eval(f'PJ.C{channel}.planet[j_alldata]')]

    # Find the corresponding emission angles 
    ea = eval(f'PJ.C{channel}.eangle[j_data]')

    # Spacecraft l
    SC_lat = PJ.ob_lat_c[j_data]
    # Height of the spacecraft 
    SC_h = PJ.range[j_data]-PJ.RJ[0]

    # Beam intercept latitude at 0 km 
    B_lat = eval(f'PJ.C{channel}.lat_c[j_data]') 
    B_0 = np.zeros_like(SC_h) 

    # Interpolate to the right height 
    dldh = (SC_lat - B_lat)/SC_h
    dlat_ray = (dldh)*(np.min(zdepth) - np.max(zdepth))
    dlat_beam = (dldh)*(zdepth[int(len(zdepth) / 2)]) 

    # Bin the data in 1 degree steps 
    lat_array =np.arange(lat_range[0],lat_range[-1]+1) 
    bincount = np.digitize(B_lat, lat_array)

    ea_max   = np.zeros(len(lat_array)) 
    dray_max = np.zeros(len(lat_array))  
    dbeam_max = np.zeros(len(lat_array)) 
    for i in range(len(lat_array)):
        ind_lat = np.where(bincount == i+1)[0]
        if ind_lat.tolist():
            ea_max[i] = np.max(ea[ind_lat]) 
            dray_max[i] =  dlat_ray[ind_lat][np.argmax(np.abs(dlat_ray[ind_lat]))] 
            dbeam_max[i] =  dlat_beam[ind_lat][np.argmax(np.abs(dlat_beam[ind_lat]))] 
        else: 
            continue 

    if plotting: 

        fig,ax = plt.subplots(1,figsize=(12,9))
        ax.plot(lat_array,np.abs(dray_max))
        ax2 = ax.twinx()
        ax2.plot(lat_array,ea_max,color='firebrick',linestyle='-.')
        ax.set_ylabel('Maximum latitude coverage (deg)')
        ax2.set_ylabel('Maximum emission angle (deg)')
        ax.set_title(f'PJ{PJ.PJnumber} C{channel}')
        ax.set_xlabel(r'Latitude$_c$ (deg)') 
        plt.show()


    return ea_max, dray_max, dbeam_max, lat_array 



def PlotRays(PJ, channel, lat_range, lat_type = 'beam', plotspacecraft=False, beamblend = False, eanglelim = [0,45], latlim = [], ylim = [], zdepth = [] ): 
    """ Show the latitudinal extent for rays penetrating the atmosphere 

    This function will plot all the rayes for a given lat_range, corresponding 
    to either all latitudes of the spacecraft, or all rays that are within 
    the given latitude range. 

    Parameters
    ----------
    PJ : object 
        [-] Contains all the relevant PJ data (obtained by running PJ.readdata() 
    channel : int 
        [-] Channel number 
    lat_range :  2x1 float 
        [-] latitude range for plotting          


    channel :  2x1 float 
        [-] spacecraft subobserver points, longitude, latitude        
    beamsize : float 
        [deg] Size of the HPBW   
    dist : float
        [m] Distance of the spacecraft from the center of the planet 

    Keyword Arguments
    ----------
    lat_type : str
        [-] 'beam' or 'spacecraft' determines if the lat_range is for beam or for spacecraft location 

    plotspacecraft : boolean 
        [-] Determines if the spacecraft altitude is set as the maximum height for plotting. 
            Use when the spacecraft is close to the planet 

    eanglelim :  2x1 float 
        [deg,deg] Set the limit for the emission angle filter [upper, lower] 


    Returns
    -------

    
    Warnings
    -------
    
    
    Example
    -------
    import pyPR.JunoTools as jt 
    PJnumber = 1
    PJ = jt.PJ(PJnumber)
    pathJ = '/Users/chris/GDrive-UCB/Berkeley/Research/Juno/'
    PJ.readdata(pathJ) 
    lat = [-1,1]

    for channel in range(1,4): 
        ea_max = np.max(jt.Raydlat(PJ,channel,lat,plotting=False)[0])
        jt.PlotRays(PJ,channel,lat,eanglelim=[38,40], plotspacecraft = True, beamblend=True)
        jt.PlotRays(PJ,channel,lat,eanglelim=[38,40], plotspacecraft = True, beamblend=True, lat_type='spacecraft',latlim=[-3,3])


    References
    ------------
    James R Wertz. Orbit and Constellation layout and administration (OCDM)
    
    Todo
    ----- 
    Notes
    -------
    08/27/19, CM, initial comit 


    """
    # Retrieval has a function to compute this 
    # rr.ProbeDepth(f,b=[0.0,0.0], thres = [0.05 , 0.5, 1, 0.5 , 0.05 ])

    if channel == 1  and not zdepth: 
        zdepth = [-16, -80, -160, -310, -920 ]
    elif channel == 2 and not zdepth:
        zdepth = [-2, -60 ,-120 , -190, -324 ]
    elif channel == 3 and not zdepth:
        zdepth = [7, -40 , -70 , -110, -160]
    elif channel == 4 and not zdepth:
        zdepth = [11, -10 , -40 , -70, -95]
    elif channel == 5 and not zdepth:
        zdepth = [14, 7, -8 , -30, - 48]
    elif channel == 6 and not zdepth:
        zdepth = [30, 20, 12 , 8, 6]


    dlat = lat_range[1] - lat_range[0]
    # Find the relevant rays to plot 
    # ------------------------------------------------------------------
    # Find all the data points corresponding to the given latitude range 
    if lat_type.lower() == 'beam': 
        # Beam intercept latitude 
        j_data = np.where(PJ.pj_ids == True)[0][0] + np.argwhere((eval(f'PJ.C{channel}.lat_c[PJ.pj_ids]') > lat_range[0]- 1 ) & (eval(f'PJ.C{channel}.lat_c[PJ.pj_ids]') < lat_range[1] + 1)) 
        # Find the corresponding emission angles 
        ea = eval(f'PJ.C{channel}.eangle[j_data]')
        # Filter out all rays larger than eanglelim 
        j_f = j_data[np.where((ea>eanglelim[0]) & (ea<eanglelim[1])) ]   

    elif lat_type.lower() == 'spacecraft': 
        # Spacecraft latitude 
        ja = np.where(PJ.pj_ids == True)[0][0] +   np.argmin(np.abs(PJ.ob_lat_g[PJ.pj_ids] - lat_range[0])) 
        jb = np.where(PJ.pj_ids == True)[0][0] + np.argmin(np.abs(PJ.ob_lat_g[PJ.pj_ids] - lat_range[1]))
        j_data = np.sort([ja,jb]) 
        j = np.linspace(j_data[0],j_data[1],j_data[1]-j_data[0] ,dtype=np.int32)
        # Emission angle  
        ea = eval(f'PJ.C{channel}.eangle[j]')  
        # Filter out all rays larger than eanglelim 
        j_f = j[np.where((ea>eanglelim[0]) & (ea<eanglelim[1])) ]   


    print(len(j_f)) 
    print(eval(f'PJ.C{channel}.eangle[j_f]') )
    print(eanglelim) 

    # Read in the Spacecraft data 
    # ------------------------------------------------------------------
    # Latitude  
    SC_lat = PJ.ob_lat_c[j_f]

    # Height of the spacecraft 
    SC_h = PJ.range[j_f]-PJ.RJ[0]

    # Read in the beam data 
    # ------------------------------------------------------------------
    # Beam intercept latitude at 0 km 
    B_lat = eval(f'PJ.C{channel}.lat_c[j_f]') 
    B_0 = np.zeros_like(SC_h)

    # Interpolate to the right height 
    dldh = (SC_lat - B_lat)/SC_h


    color = (dldh)*np.abs(np.max(zdepth) - np.min(zdepth)) 
    plotcolor = np.abs(color)/np.max(np.abs(color)) # Adjustment to only plot the darker colors in the color map

    # Propagate the ray until it terminates  
    B_lat +=  + dldh*zdepth[-1] 
    B_h  = np.ones_like(SC_h)*zdepth[-1] 


    if beamblend: 
        d_fp = np.zeros((2,len(j_f)))
        cnt = 0 
        for i in j_f: 
            beam = np.radians([eval(f'PJ.C{channel}.lon[i]'), B_lat[cnt]])
            obs = np.radians([PJ.ob_lon[i], PJ.ob_lat_c[i]])
            dist = PJ.range[i]*1e3
            [lon_fp,lat_fp],horizon = ProjectedFootprint(beam,obs,eval(f'PJ.C{channel}.hpbw')*2,dist)  
            # The algorithm works by stepping through a circle of the beam. First entry, and half way point are the points closest and furthest away from the spacecraft 

            # obtain the maximum difference in degrees  
            n = lon_fp.size
            d_fp_0 = ((lon_fp[0]-beam[0])**2 + (lat_fp[0]-beam[1])**2)**0.5 
            d_fp_pi = ((lon_fp[int(n/2)+1]-beam[0])**2 + (lat_fp[int(n/2)+1]-beam[1])**2)**0.5 
            d_fp[:,cnt] = [d_fp_0*180/np.pi,d_fp_pi*180/np.pi]
            # mu = eval(f'PJ.C{channel}.eangle[i]')
            # print(f'mu {mu:2.2f}, d_fp_0 {d_fp_0*57.3:2.2f}, d_fp_pi {d_fp_pi*57.3:2.2f}, horizon {horizon[50]}')
            cnt += 1



    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111, polar=True)
    if not plotspacecraft: 
        ax.set_rlim(bottom=zdepth[-1], top=zdepth[0])
    else: 
        ax.set_rlim(bottom=zdepth[-1], top=np.max(SC_h))

    ax.set_rorigin(-PJ.r_j[0])
    ax.set_title(f'C{channel} rays')

    if not latlim: 
        ax.set_thetamin(lat_range[0])
        ax.set_thetamax(lat_range[1])
    else: 
        ax.set_thetamin(latlim[0])
        ax.set_thetamax(latlim[1])

    ax.set_yticks(zdepth[::2])
    if not latlim: 
        ax.set_xticks(np.linspace(np.radians(lat_range[0]),np.radians(lat_range[1]),(lat_range[1]-lat_range[0]+1)))
    else: 
        ax.set_xticks(np.linspace(np.radians(latlim[0]),np.radians(latlim[1]),(latlim[1]-latlim[0]+1)))

    # ax.set_yticklabels(zdepth[::2], rotation=-90)
    yticks = ax.get_yticks() 
    ax.set_yticklabels([])
    xlimloc = ax.get_xlim() 
    xtickstep = ax.get_xticks()[1] - ax.get_xticks()[0]
    for yloc, ylabel in zip(yticks,zdepth[::2]): 
        ax.text(xlimloc[0] - 0.05*dlat*(xtickstep), yloc, ylabel, rotation=-90)

    # Plot the spacecraft 
    for i in range(len(SC_lat)): 
        if color[i] > 1: 
            ax.plot(np.radians([SC_lat[i],B_lat[i]]),[SC_h[i],B_h[i]],color = 'firebrick', alpha=1)
        else: 
            ax.plot(np.radians([SC_lat[i],B_lat[i]]),[SC_h[i],B_h[i]],color = cmap(1 - plotcolor[i]), alpha=plotcolor[i])
            if beamblend: 
                ax.plot(np.radians([SC_lat[i],B_lat[i]-d_fp[0,i]]),[SC_h[i],B_h[i]],linewidth = 1, linestyle=':',color = cmap(1 - plotcolor[i]), alpha=plotcolor[i])
                ax.plot(np.radians([SC_lat[i],B_lat[i]+d_fp[1,i]]),[SC_h[i],B_h[i]],linewidth = 1, linestyle=':',color = cmap(1 - plotcolor[i]), alpha=plotcolor[i])

    # Plot the surface 
    ax.plot(np.radians(B_lat),np.zeros_like(B_lat),color = 'gray', alpha=1)


    return  

def surface_normal(lat_g, lon_w, ob_lon):
    '''Returns the normal vector to the surface of the planet.
    Take dot product with sub-obs or sub-sun vector to find cosine of emission angle 

    lat_g = -7.0
    lon_w = 2
    ob_lon = -5.8 
    ob_lat = 4.2
    surf_n = surface_normal(lat_g, lon_w, ob_lon)
    ''' 
    nx = np.cos((lat_g))*np.cos((lon_w-ob_lon))
    ny = np.cos((lat_g))*np.sin((lon_w-ob_lon))
    nz = np.sin((lat_g))
    return np.asarray([nx,ny,nz])

def emission_angle(ob_lat, surf_n):
    '''Return the cosine of the emission angle of surface wrt observer'''
    ob = np.asarray([np.cos(np.radians(ob_lat)),0,np.sin(np.radians(ob_lat))])
    return np.dot(surf_n.T, ob).T 


def local_position_ellipsoid(lon,lat,R):
    '''
    Calculate the local position on the ellipsoid for a given latitude and longitude 
    '''

    x   =   R[0]*np.cos(lon)*np.sin(np.pi/2-lat)   
    y   =   R[1]*np.sin(lon)*np.sin(np.pi/2-lat)   
    z   =   R[2]*np.cos(np.pi/2-lat)

    return np.array([x,y,z])

def surface_normal_ellipsoid(lon,lat,R): 
    '''
    Calculate surface normal for an ellipsoid of dimension R at loaction r 
    '''
    r = local_position_ellipsoid(lon,lat,R)

    return np.array([r[0,:,:]/R[0]**2,r[1,:,:]/R[1]**2,r[2,:,:]/R[2]**2])

def local_radius(lon,lat,R):
    
    r_s = local_position_ellipsoid(lon,lat,R)

    return np.sqrt(r_s[0,:,:]**2 + r_s[1,:,:]**2 + r_s[2,:,:]**2)


def EmissionAngleFF(lon,lat,obs): 
    '''
    For a given spacecraft and beam location calculate the corresponding 
    emission angles. Farfield assumption

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

    lat = np.arange(-1,5,0.5) 
    lon = np.arange(-9,-5,0.5)
    obs = [-5.8, 4.2] 
    '''
    ob_lon, ob_lat = obs

    # Construct a grid for all surface points and compute their normals 
    lo, la = np.meshgrid(lon,lat)
    surf_n = surface_normal(la, lo, ob_lon)

    # Compute the emisson angle 
    ea = emission_angle(ob_lat, surf_n)

    return ea, lo, la 


def EmissionAngle(obs, lon, lat, RJ = [71492., 71492., 66854.]): 
    '''
    For a given spacecraft and beam location calculate the corresponding 
    emission angles. 

    Parameters
    ----------
    S : 3x1 float 
        [km,rad,rad] Spacecraft position (r,lon,lat)

    beam : 2x1 float 
        [rad,rad] beam boresight, longitude, latitude 

    obs :  2x1 float 
        [rad] spacecraft subobserver points, longitude, latitude        
    beamsize : float 
        [deg] Size of the HPBW   
    dist : float
        [m] Distance of the spacecraft from the center of the planet 

    Keyword Arguments
    ----------
    
    import pyPR.JunoTools as jt
    dist = 75601.17
    obs = [dist,np.radians(-5.79556788),  np.radians(4.43404659)]
    lat = np.radians(np.arange(-1,5,0.5) ) 
    lon = np.radians(np.arange(-9,-5,0.5)) 
    
    obs = np.array([79075.245, 0.00012249427993173018, 0.009174527763385208]) 


    ea,lo,la = jt.EmissionAngle(obs,lon,lat) 


  
    '''

    from numpy.linalg import norm 

    # Compute the surface normal centered around the observer longitude 
    ob_lon, ob_lat = obs[1],obs[2]

    # Construct a grid for all surface points centered around the spacecraft location and compute their normals 
    lo, la = np.meshgrid(lon,lat)
    surf_n = surface_normal_ellipsoid(lo, la, np.array(RJ))

    # Compute the relevant vectors (Ray, Spacecraft, Local, Ray) 
    r_s = local_radius(lo,la,RJ)

    # Spacecraft in cartesian coordinates
    s = polar2cart(obs[0],obs[1],obs[2])

    # Adjust dimensions to match L
    S = np.array([np.ones_like(r_s)*s[0],np.ones_like(r_s)*s[1],np.ones_like(r_s)*s[2]])

    # Beam intercept
    L = polar2cart(r_s,lo,la)

    # Ray From spacecraft to local intercept point 
    R = S - L 

    # Emission angle (dot protuct between local normal and the Ray 
    ea = np.arccos((np.sum([R[0]*surf_n[0] ,R[1]*surf_n[1] ,R[2]*surf_n[2] ],axis=0))/(norm(surf_n,axis=0)*norm(R,axis=0)))
 

    return ea, lo, la 


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
    There is a problem where nan is returned for phi = 180deg, so the anti spacecraft point 
    
    Example
    -------
    import pyPR.JunoTools as jt
    #beam = np.array([-0.59377846, -0.99061398])
    #obs = np.array([-0.27264095, -0.59551803])
    beam = np.radians(np.array([-64, -66]))
    obs = np.radians(np.array([-30, -26]))
    beam = np.radians([-6.74549685,  1.94814349])
    obs = np.radians([-5.79556788,  4.43404659])
    beamsize = 20.6
    dist = 75601.17*1e3 
    d2r = 180/np.pi
    [lon_fp,lat_fp],horizon = jt.ProjectedFootprint(beam,obs,beamsize,dist)  
    # Compute the emission angle for the given range 

    fig, axs = plt.subplots(figsize=(8,4))
    axs.plot(beam[0]*d2r,beam[1]*d2r,'*b',label='beam')  
    axs.plot(obs[0]*d2r,obs[1]*d2r,'*r',label='sc')
    axs.plot(lon_fp*d2r,lat_fp*d2r)
    axs.set_xlim([-10,10])

    #compute the ea outlines for the given axis 
    lon_range = np.arange(np.floor(axs.get_xlim()[0]),np.ceil(axs.get_xlim()[1])+1,0.5) 
    lat_range = np.arange(np.floor(axs.get_ylim()[0]),np.ceil(axs.get_ylim()[1])+1,0.5) 

    sc = [dist*1e-3,obs[0],obs[1]]
    ea,lo,la = jt.EmissionAngle(sc, np.radians(lon_range),np.radians(lat_range))

    
    cb = plt.contourf(lo*57.3,la*57.3,ea*57.3)
    plt.colorbar()
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
    return np.array([
         r * np.sin(np.pi/2-phi) * np.cos(theta),
         r * np.sin(np.pi/2-phi) * np.sin(theta),
         r * np.cos(np.pi/2-phi)
    ])


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
    nparam = 2
    Tn = np.copy(a) 
    ld = np.zeros_like(a)
    sig = np.zeros((len(a),nparam))

    # Loop through the array 
    for i in range(len(a)-window): 
        suba = a[i:i+window]
        subb = mu[i:i+window]

        # Build the A matrix 
        A = np.vstack([np.ones_like(suba),np.log(np.cos(np.radians(subb)))]).T  
        # Build the weights 
        W = np.diag(1/subb**weights)
        # Build the solution matrix 
        b = np.log(suba)
        db = np.log(suba*0.015)

        # Left hand side 
        ATW = A.T@W
        ATWA = ATW@A 
        iATWA = np.linalg.inv(A.T@W.T@A)
        iATWTA = np.linalg.inv(A.T@W.T@A)
        # Right hand side 
        ATWb = ATW@b 
        ATWdb = ATW@db 
        WTA = W.T@A
        logT,p = iATWA@ATWb

        Tn[i+int(window/2)] = np.exp(logT)
        ld[i+int(window/2)] = p

        # Reduced chi-squared statistic 
        My = b - A@([logT,p]) # Calculate the residuals    
        sigsq = My.T@My/len(suba) # Calculate the chi-squared 
        RCS = np.linalg.inv(A.T@A)*sigsq # Calculate the covariance 
        # Uncertainties are on the diagonal 
        sig[i+int(window/2),0] = np.exp(logT) * (np.exp(np.sqrt(RCS[0,0])) - 1) 
        sig[i+int(window/2),1] = np.sqrt(RCS[1,1]) 

        # Experiment with given measurement uncertainty first 
        SIGSQ = np.diag(np.log(suba*0.015))**2
        COV = iATWA@ATW@SIGSQ@WTA@iATWTA
        # https://www.jstor.org/stable/2346491?seq=3#metadata_info_tab_contents
        sig_prop=iATWA@ATWdb

    return Tn,ld, sig




def fit_Tn(a,mu,p,window=200, weights=2,): 
    '''
    Try a simple linear least squares approach 
    '''

    '''
    T = T_s[12000:13000]
    W = 1/mu_s[12000:13000]
    
    '''
    # Make a copy of the array. The beginning and the end of the array will be copied 
    nparam = 1
    Tn = np.copy(a) 
    ld = np.zeros_like(a)
    sig = np.zeros((len(a),nparam))

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
        ATW = A.T@W
        ATWA = ATW@A 
        iATWA = np.linalg.inv(A.T@W.T@A)
        iATWTA = np.linalg.inv(A.T@W.T@A)
        # Right hand sie 
        ATWb = ATW@b 
        WTA = W.T@A
        T = iATWA@ATWb

        Tn[i+int(window/2)] = T
        ld[i+int(window/2)] = p


        # Reduced chi-squared statistic 
        My = b - A@T  # Calculate the residuals  
        sigsq = My.T@My/len(suba) # Calculate the chi-squared 
        RCS = np.linalg.inv(A.T@A)*sigsq # Calculate the covariance 
        # Uncertainties are on the diagonal 
        sig[i+int(window/2),0] = np.sqrt(RCS)


    return Tn, sig



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
            sumw = np.sum(subw)

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