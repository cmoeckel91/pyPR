"""Initialize python repository for pyPR geometry.

Notes
-----
7/8/18,CM, Initial Commit
"""

import sys, os

import numpy as np

# Debug 
import warnings
import IPython


def locate2template(logger,output,analytics=False): 
    """ Reads out the located visibilities and transform them into 
    flagging template format 

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
    
    from pyPR.CASAtools import * 
    logger = '/Users/chris/Documents/Research/Toolbox/Testfiles/casa-20180117-011733.log'
    test  = locate2template(logger, 'test')

    # testing 
    logger = '/Users/chris/Documents/Research/Toolbox/Testfiles/Locate_format' 
    filepath = logger
    fp = open(filepath)
    line = fp.readline()
    line = fp.readline()
    line = fp.readline()
    info = line.split() 
    References
    ------------
    

    Todo
    ----- 
    
    Notes
    -------
    mm/dd/yy, Initials of Author, Short description of update
    """
    

    # Loop over all the input files 

    import re 
    import sys
    from select import select

    

    flags = []

    # Loop over all logger files 
    for logfiles in range(len(logger)):

        # Loop through the logger file 
        filepath = logger
        with open(filepath) as fp:  
            line = fp.readline()
            cnt = 1
            corrupteddatacnt = 0
            while line:
                info = line.split() 
                # Find locate information 
                if info[3] == 'PlotMS::locate+': 

                    for i in range(4,len(info)): 
                        subline = info[i]
                        subinfo = subline.split('=', 1)

                        # Scan 
                        if subinfo[0] == 'Scan':
                            try: 
                                #scan = int(re.search(r'\d+',info[4],re.I)[0])
                                scan = subinfo[1]
                            except: 
                                sys.warnings('Failed to identify scan in line %d'.format(cnt))
                        
                        # Field 
                        elif subinfo[0] == 'Field':
                            try:
                                # field = int(re.search(r'\[\d\]',info[5],re.I)[0][1:-1])
                                field = subinfo[1].split('[',1)[0]
                                fieldid = re.search(r'\[(.*)\]', subinfo[1]).group(1)
                            except: 
                                sys.warnings('Failed to identify field in line %d'.format(cnt))
                        
                        # Time 
                        elif subinfo[0] == 'Time':
                            time = subinfo[1]
                        
                        # Antenna 
                        elif subinfo[0] == 'ANT1': 
                            ant1 = subinfo[1].split('@',1)[0]#[-2:]
                            ant2 = float('NaN') 

                        # Baseline 
                        elif subinfo[0] == 'BL': 
                            ant1 = subinfo[1].split('@',1)[0]#[-2:]
                            i +=2 
                            subline = info[i]
                            ant2 = subline.split('@',1)[0]#[-2:]
                        
                        # Spectral window 
                        elif subinfo[0] == 'Spw': 
                            spw_r = re.search(r'\<(.*)\>', subinfo[1])
                            if spw_r: 
                                spw_r = spw_r.group(1).split('~',1)
                                spw = arange(int(spw_r[0]),int(spw_r[1]))
                            else: 
                                spw = subinfo[1] 
                        
                        # Channel 
                        elif subinfo[0] == 'Chan': 
                            ch_r = re.search(r'\<(.*)\>', subinfo[1])
                            if ch_r: 
                                ch_r = ch_r.group(1).split('~',1)
                                if ch_r[0] == ch_r[1]:
                                    chan = ch_r[0]
                                else:
                                    chan = ch_r[0]+'~'+ch_r[1]
                            else: 
                                chan = subinfo[1] 

                        # Correlation 
                        elif subinfo[0] == 'Corr':
                            corr = subinfo[1] 

                    try:        
                        data_point = [scan, field, fieldid, time, ant1, ant2, spw, chan, corr ]
                        flags.append(data_point)
                        cnt  += 1
                    except: 
                        corrupteddatacnt += 1 
                        print(info)

                line = fp.readline()

        fp.close()
                 


    print(cnt) 
    print(corrupteddatacnt)

    # Do analytics 
    if analytics:
        from collections import Counter
        scan    = [item[0] for item in data]
        field   = [item[1] for item in data]
        fieldid = [item[2] for item in data]
        time    = [item[3] for item in data]
        ant1    = [item[4] for item in data]
        ant2    = [item[5] for item in data]
        spw     = [item[6] for item in data]
        chan    = [item[7] for item in data]
        corr    = [item[8] for item in data]

        Counter(ant1)



    # Write into the flagging template format 
    if os.path.exists(output): 
        timeout = 10
        print("The file exists already. Are you sure you want to overwrite it?: (y/N)")
        rlist, _, _ = select([sys.stdin], [], [], timeout)
        if rlist:
            s = sys.stdin.readline()
            if s.lower() == 'y\n' or s.lower() == 'yes\n': 
                print('File will be overwritten!')
                os.system('rm -rf ' + output) 
            else: 
                output = output + '_vx'
                print('File will NOT be overwritten! Output can be found under ' + output)
        else:
            output = output + '_vx'
            print('File will NOT be overwritten! Output can be found under ' + output)

    filepath = output + '.txt'
    with open(filepath,'w') as fo:  
        fo.close() 

    # Create file to be written to 

    filepath = output + '.txt'
    with open(filepath,'a') as fo:  
        for i in range(len(flags)):

            # data = [scan, field, fieldid, time, ant1, ant2, spw, chan, corr ]
            line = ('mode=\'manual\' scan=\'{:s}\' field=\'{:s}\'' 
                    ' timerange=\'{:s}\' antenna=\'DV{:s}&DV{:s}\''
                    ' spw=\'{:s}\' chan=\'{:s}\' corr=\'{:s}\'\n'.format(
                        flags[i][0],flags[i][1],flags[i][3],flags[i][4][-2:],
                        flags[i][5][-2:],flags[i][6],flags[i][7],flags[i][8] ))
            
            fo.write(line)

    fo.close()

    return flags 

                # read out the data 
            #     if info[3] == 'PlotMS::locate+': 
            #         # Scan 
            #         try: 
            #             scan = int(re.search(r'\d+',info[4],re.I)[0])
            #         except: 
            #             sys.warnings('Failed to identify scan in line %d'.format(cnt))
            #         # Field 
            #         try:
            #             field = int(re.search(r'\[\d\]',info[5],re.I)[0][1:-1])
            #         except: 
            #             sys.warnings('Failed to identify field in line %d'.format(cnt))
            #         # time 
            #         time = info[6][5:] 
            #         # Antenna 
            #         # spw 
            #         try: 
            #             spw = int(re.search(r'\d+',info[9],re.I)[0])
            #          except: 
            #             sys.warnings('Failed to identify spw in line %d'.format(cnt))
            #         # Channel 
            #         try: 
            #             Chan = int(re.search(r'\d+',info[10],re.I)[0])
            #          except: 
            #             sys.warnings('Failed to identify channel in line %d'.format(cnt))
            #         # Corr 
            #         Corr = info[12][-2:]


            # line = fp.readline()
            # cnt += 1








def ldvis(k,amp,lam,B,a): 
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
    Josh Ivy league education

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

def shortspacingobservatory(nu,uvhole,name, obs='VLA',n_ants = 30, filepath='./'): 
    """Create an oblate spheroid 



    Returns
    ----------


    Keywords
    ----------


    Example
    -------
    import CASAtools
    filepath = '/Users/chris/GDrive-UCB/Berkeley/Research/VLA/VLA2017/Jup_x_20170111/'
    name = 'vca.tp'  
    uvhole = 18000 # [lambda] 
    nu = 10e9 # [Hz] Observing frequency 
    CASAtools.shortspacingobservatory(nu,uvhole,name,filepath = filepath,n_ants = 70) 
    os.system('cat ' + filepath+name + '.cfg') 

    References
    ------------
    Triaxial ellipsoid http://mathworld.wolfram.com/Ellipsoid.html

    Notes
    -------
    10/24/2017, CM, Initial Commit
    """

    from astropy import constants as cst 
    import os

    if obs.lower().strip() == 'vla':
        d_ant = 25 # [m] antenna diameter
    elif obs.lower().strip() == 'almba':
        d_ant = 12 # [m] antenna diameter 
    else: 
        print('Observatory not yet implemented. Should be easy to add') 

    # The radius in which the antennas can be found. 
    spread = uvhole*cst.c.value/nu/2

    # Create the inner ring 
    n_ring = 8 
    ring_ang = np.arange(0,2*np.pi,2*np.pi/n_ring)
    ring_loc = np.zeros([8,2])
    for i in range(n_ring): 
        ring_loc[i,:] = np.cos(ring_ang[i]), np.sin(ring_ang[i])

    ring_loc = ring_loc*d_ant/2
    
    loc = np.multiply(np.random.rand(n_ants,2)*spread,np.random.choice([-1, 1],  size=(n_ants,2), p=[1./2, 1./2])) 


    file = filepath + name + '.cfg'
    # print(os.system('ls '+ filepath))
    with open(file,'w+') as fo: 
        fo.write('# observatory='+obs.upper() + '\n') 
        fo.write('# coordsys=LOC (local tangent plane)' + '\n') 
        fo.write('0         0           0           ' + str(d_ant) + '          A' + '01\n' )
        for i in range(n_ring):
            fo.write('{:.0f}         {:.0f}         0          {:d}           A{:02d}\n'.format(ring_loc[i,0],ring_loc[i,1],d_ant,i + 2) )

        for j in range(n_ants-n_ring-1):
            fo.write('{:.0f}         {:.0f}         0          {:d}           A{:02d}\n'.format(loc[j,0],loc[j,1],d_ant,j+n_ring + 2) )

    fo.close() 

    return 


def parrallel_Miriad_script(m_ncore, uvfits, latrange, latint, cell, 
    planet = 'jupiter', 
    filepath_data= './', 
    filepath_script = 'Parallel_facets.bsh', 
    tmp_directory = '/Volumes/scratch/tmp'): 
    """Write a bash script that allows for parrallel execution 
    of Bob's Deprojection technique in MIRIAD 

    You need to be within the folder where the uv file sits normally. 
    The script will generate temporary directories and distribute them 

    Returns
    ----------


    Keywords
    ----------


    Example
    -------
    import pyPR
    import CASAtools
    import os

    cd '/Volumes/CASA/chris/2017-Jan-11/Deprojection/spw2~33_p0'

    m_ncore = 24
    uvfits = 'jup-x.uv.comp' 
    latrange = 120  
    latint = 2.5 
    cell = 0.039 
    planet = 'jupiter'
    CASAtools.parrallel_Miriad_script(m_ncore, uvfits, latrange, latint, cell, planet = 'jupiter',filepath_script = 'Parallel_facets_v2.bsh',tmp_directory = '/tmp') 

    References
    ------------
    

    Notes
    -------

    """



    # Create new subdirectories 
    import os
    import sys
    from select import select
    import numpy as np 

    # Calculate the corresponding latitude ranges 
    dlat = latrange/np.ceil(latrange/latint)
    # Number of latitude circles 
    nlat = latrange/dlat 

    # Find optimum number of cores (Assures that they will finish around the same time) 
    cond = True
    perf_p = 0.5
    ncore = m_ncore
    while cond: 
        if np.ceil(nlat/ncore) != np.round (nlat/ncore):  
            ncore -=1 
        else: 
            perf = np.ceil(nlat/ncore  ) - (nlat/ncore)  
            perf_r = np.ceil(nlat/(ncore-1)) - nlat/(ncore-1)
            if  perf_r < perf: 
                ncore -=1 
            else: 
                cond = False 

    # Number of latitude bands per core 
    nband = np.ceil(nlat/ncore) 

    lat_lower = -latrange/2 


    overwrite = False
    temp_folder = 'temp_p'
    if os.path.exists(temp_folder+'0'): 
        timeout = 10
        print("The folders exists already. Are you sure you want to overwrite it? Abort after 10 seconds: (y/N)")
        rlist, _, _ = select([sys.stdin], [], [], timeout)
        if rlist:
            s = sys.stdin.readline()
            if s.lower() == 'y\n' or s.lower() == 'yes\n': 
                print('Folders will be overwritten!')
                os.system('rm -rf ' + temp_folder + '*') 
                overwrite = True
            else: 
                print('Assume files are there already')
        else: 
            sys.exit('Move your files manually')
    else: 
        overwrite = True


    for i in range(ncore): 
        temp_name = temp_folder+ str(i) # Mkdir where temp data are stored 

        if overwrite:
            os.system('mkdir ' + temp_name )
        else: 
            os.system('rm -rf {:s}/params.pl'.format(temp_name) )
        # os.system('cp -r' + uvfits + ' ' + temp_name +'/'  ) # Copy the file to temp folder 
        # Write params.file 


        # Calculate the latitude range 
        lat_upper = lat_lower + nband*dlat 
        if lat_upper > latrange/2: 
            lat_upper = latrange/2


        filepath = temp_name+'/params' + '.pl' 
        with open(filepath,'w') as fo:  
            fo.write('$planet = "{:s}";# Planet name\n'.format(planet))
            fo.write('$vis = "../{:s}";     # Visibility file\n'.format(uvfits))
            fo.write('$cell = {:4.4f};      # Image pixel size, in arcseconds (!).\n'.format(cell))
            fo.write('$minlat = {:2.10f};      # Min latitude (degrees) to map\n'.format(lat_lower))
            fo.write('$maxlat = {:2.10f};       # Max latitude\n'.format(lat_upper))
            fo.write('$minlon = 0;        # Min longitude (degrees)\n')
            fo.write('$maxlon = 360;      # Max longitude\n')
            fo.write('$latint = {:2.1f};      # Increment in latitude for facets, in degrees.\n'.format(latint))
            fo.write('$imsize = 150;          # Facet pixel size.\n')
            fo.write('#$fwhm_km = 3200;       # Optional extra smoothing function.\n')
            fo.write('$robust = 0.0;          # Optional imaging weighting parameter.\n')
            fo.write('$plradius = 71492.0     # Optional planet radius in kilometers.\n')
            fo.write('# $obstime = ""  # Optional observation time used for geometry.\n')
            fo.close()

        lat_lower = lat_upper+dlat


        # Make a bin file that you can execute that points to all the correct points 

    bashcmd = 'mkdir facets && echo "Script is running" && '
    with open(filepath_script,'w') as fo:
        fo.write('#!/bin/bash\n')
        for i in range(ncore):
            # Reduce memory allocation problems by running them out of phase 
            if i < ncore/2:
                # bashcmd = bashcmd + (' \n(rm -rf ~/tmp{:d} && mkdir ~/tmp{:d} && TMPDIR="~/tmp{:d}" && cd temp_p{:d} && nohup perl /usr/local/miriad/bin/darwin/facets.pl && cp facets/* ../facets/) &'.format(i))
                bashcmd = bashcmd + ('\n(rm -rf {:s}{:d} && mkdir {:s}{:d} && TMPDIR="{:s}{:d}" && cd temp_p{:d} && perl /usr/local/miriad/bin/darwin/facets.pl &> logger.txt  && cp -r facets/* ../facets/) &'.format(tmp_directory,i,tmp_directory,i,tmp_directory,i,i))
            else: 
                bashcmd = bashcmd + ('\n(rm -rf {:s}{:d} && mkdir {:s}{:d} && TMPDIR="{:s}{:d}" && cd temp_p{:d} && sleep && perl /usr/local/miriad/bin/darwin/facets.pl &> logger.txt  && cp -r facets/* ../facets/) &'.format(tmp_directory,i,tmp_directory,i,tmp_directory,i,i))

        bashcmd = bashcmd +('&\nmail -s "Facetting done" chris.moeckel@berkeley.edu <<< " "')
        fo.write(bashcmd)
    # Change permisson   
    os.system('chmod u+x ' + filepath_script)


    # Also write a master script and place in main folder 

    filepath = 'params' + '.pl' 
    with open(filepath,'w') as fo:  
        fo.write('$planet = "{:s}";# Planet name\n'.format(planet))
        fo.write('$vis = "{:s}";     # Visibility file\n'.format(uvfits))
        fo.write('$cell = {:4.4f};      # Image pixel size, in arcseconds (!).\n'.format(cell))
        fo.write('$minlat = {:2.10f};      # Min latitude (degrees) to map\n'.format(-lat_range/2))
        fo.write('$maxlat = {:2.10f};       # Max latitude\n'.format(lat_range/2))
        fo.write('$minlon = 0;        # Min longitude (degrees)\n')
        fo.write('$maxlon = 360;      # Max longitude\n')
        fo.write('$latint = {:2.1f};      # Increment in latitude for facets, in degrees.\n'.format(latint))
        fo.write('$imsize = 150;          # Facet pixel size.\n')
        fo.write('#$fwhm_km = 3200;       # Optional extra smoothing function.\n')
        fo.write('$robust = 0.0;          # Optional imaging weighting parameter.\n')
        fo.write('$plradius = 71492.0     # Optional planet radius in kilometers.\n')
        fo.write('# $obstime = ""  # Optional observation time used for geometry.\n')
        fo.close()


    return 

# itemize in=temp_p2/facets/n50p11.icln 
