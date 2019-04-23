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


def locate2template(loggers,output,analytics=False): 
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
    logger = ['/Users/chris/Documents/Research/Toolbox/Testfiles/casa-20180117-011733.log']
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

    # Make sure that you are not re-reading the files 
    if len(loggers) > 3: 
        print('Currently reading {:d} log files. If that is not correct, \
           make sure your log files are given as a list [log1, log2,]\
        '.format(len(logger)))

    flags = []

    # Loop over all logger files 
    for logfiles in loggers:

        # Loop through the logger file 
        filepath = logfiles
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
        scan    = [item[0] for item in flags]
        field   = [item[1] for item in flags]
        fieldid = [item[2] for item in flags]
        time    = [item[3] for item in flags]
        ant1    = [item[4] for item in flags]
        ant2    = [item[5] for item in flags]
        spw     = [item[6] for item in flags]
        chan    = [item[7] for item in flags]
        corr    = [item[8] for item in flags]

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




def flagstatistics(flags, spwrange = np.arange(0,99), minfreq=10, targetfieldID=1 , plotting=True ): 



    # # Sort by spw for example 
    # sbs = sorted(flags, key=itemgetter(6)) 

    # # fsbs 
    # ssbs =  [t for t in flags if int(t[6]) == 33]

    # tlist = ssbs 


    # see how often certain baselines pop up 
    from collections import Counter
    from operator import itemgetter  
    import matplotlib.pyplot as plt


    # 
    spwrange = np.append(np.array([-1]),spwrange,)
    # Loop through each spw and look for affected data per spw 
    for ispw in spwrange: 
        if ispw == -1: 
            tlist = flags 
        else: 
            tlist = [t for t in flags if int(t[6]) == ispw] 
        
        # Skip if there are no information on this specific spws 
        if tlist ==[]: 
            continue 

        scan    = [item[0] for item in tlist]
        field   = [item[1] for item in tlist]
        fieldid = [item[2] for item in tlist]
        time    = [item[3] for item in tlist]
        ant1    = [item[4] for item in tlist]
        ant2    = [item[5] for item in tlist]
        spw     = [item[6] for item in tlist]
        chan    = [item[7] for item in tlist]
        corr    = [item[8] for item in tlist]

        bl = [(item[4]+item[5]) for item in tlist]

        # Find the most common baselines 
        # https://stackoverflow.com/questions/20950650/how-to-sort-counter-by-value-python
        blc = Counter(bl) 

        # Create statstics for the analysis across all flags 

        if ispw == -1: 
            # Find all the flagged baselines 
            temp = list(blc.items())
            u_bl = [item[0] for item in temp]  
            nbl = len(u_bl)
            # Obtain the number of spws involved 
            spwa = np.array(spw)
            _, idx = np.unique(spwa, return_index=True)
            u_spw=spwa[np.sort(idx)] # Unique values in spw) 
            nspw = len(u_spw) 
            # Create an array where we can sum the individual spw 
            Blspw = np.zeros([nbl,nspw])


        # Populate the Blspw matrix 
        if ispw != -1: 
            ind_spw = list(u_spw).index(str(ispw))
            temp = list(blc.items())
            for i in range(len(temp)): 
                ind_bl = u_bl.index(temp[i][0])
                Blspw[ind_bl,ind_spw] = temp[i][1]
            

        # Show the antenna histogram 
        # labels, values = zip(*Counter(ant1+ant2).items())
        # indexes = np.arange(len(labels))
        # width = 1
        # plt.figure()
        # plt.bar(indexes, values, width)
        # plt.xticks(indexes + width * 0.5, labels)
        # plt.title('Flagged antennas for spw{:d}'.format(ispw))
        # plt.show()


        # Baseline histogram 
        labels, heights = zip(*sorted(((k, v) for k, v in blc.items()), key=itemgetter(1), reverse=True))
        # Add the baseline sign into this 
        labels = [item[0:4]+'&'+item[4:8] for item in labels]
        # Remove baselines below a certain value
        heights_f,labels_f = np.array([]),  np.array([]), 
        for i in range(len(heights)): 
            if heights[i] > minfreq: 
                heights_f = np.append(heights_f,(heights[i]))
                labels_f = np.append(labels_f,(labels[i]))

        if  plotting:  
        # lefthand edge of each bar
            if labels != '': 
                left = np.arange(len(heights_f))
                fig, ax = plt.subplots(1, 1)
                ax.bar(left, heights_f, 1)
                ax.set_xticks(left)
                ax.set_xticklabels(labels,rotation=45,  fontsize='small')
                plt.title('Flagged baseline for spw{:d}'.format(ispw))
                plt.show()              

        flagstring = ('mode=\'manual\' field=\'{:d}\' spw=\'{:d}\' '.format(targetfieldID, ispw),('antenna=\''+'{:s};'*len(labels_f)+'\'').format(*labels_f))
        print(flagstring[0], flagstring[1][0:-2]+'\'')

    if plotting: 
    # Stacked barplot 
        fig, ax = plt.subplots(1, 1)
        plt.imshow(Blspw.T)
        plt.ylabel('Spectral window')
        plt.xlabel('Baselines')
        plt.show()

        # Flags per baseline 
        fpbl = np.sum(Blspw, axis = 1) 
        sort_ind_bl = np.argsort(fpbl)[::-1][:len(fpbl)]

        
        # Sort by spw 
        fl_spw = [float(t) for t in u_spw]  
        sort_ind_spw = np.argsort(fl_spw)  
     

        # Plot only the major baselines 
        sort_ind_bl_fil = np.array([])
        for i in range(len(sort_ind_bl)): 
            if fpbl[sort_ind_bl[i]]> minfreq*10:
                sort_ind_bl_fil = np.append(sort_ind_bl_fil,int(sort_ind_bl[i]))

        sort_ind_bl_fil = [int(t) for t in sort_ind_bl_fil ]       
            
       
        # Sort to have baselines sorted 
        Blspw_sort = Blspw[sort_ind_bl_fil,:]
        # Sort by spw 
        Blspw_sort = Blspw_sort[:,sort_ind_spw]

        series_labels = u_spw[sort_ind_spw].tolist() 

        data = Blspw_sort.T.tolist() 

        category_labels = np.array(u_bl)[sort_ind_bl_fil].tolist()# baselines 

        stacked_bar(
            data, 
            series_labels, 
            category_labels=category_labels, 
            show_values=False, 
            value_format="{:.1f}",
            y_label="Flagged amplitudes ", grid=False
        )

        plt.show()



    return 

def stacked_bar(data, series_labels, category_labels=None, 
                show_values=False, value_format="{}", y_label=None, 
                grid=True, reverse=False):
    """Plots a stacked bar chart with the data and labels provided.

    Keyword arguments:
    data            -- 2-dimensional numpy array or nested list
                       containing data for each series in rows
    series_labels   -- list of series labels (these appear in
                       the legend)
    category_labels -- list of category labels (these appear
                       on the x-axis)
    show_values     -- If True then numeric value labels will 
                       be shown on each bar
    value_format    -- Format string for numeric value labels
                       (default is "{}")
    y_label         -- Label for y-axis (str)
    grid            -- If True display grid
    reverse         -- If True reverse the order that the
                       series are displayed (left-to-right
                       or right-to-left)
    
    # reference https://stackoverflow.com/questions/44309507/stacked-bar-plot-using-matplotlib 
    """

    import numpy as np
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1)

    ny = len(data[0])
    ind = list(range(ny))

    axes = []
    cum_size = np.zeros(ny)

    data = np.array(data)

    if reverse:
        data = np.flip(data, axis=1)
        category_labels = reversed(category_labels)

    for i, row_data in enumerate(data):
        axes.append(plt.bar(ind, row_data, bottom=cum_size, 
                            label=series_labels[i]))
        cum_size += row_data

    if category_labels:
        plt.xticks(ind, category_labels)
        ax.set_xticklabels(category_labels,rotation=90,  fontsize='small')

    if y_label:
        plt.ylabel(y_label)



    plt.legend()

    if grid:
        plt.grid()

    if show_values:
        for axis in axes:
            for bar in axis:
                w, h = bar.get_width(), bar.get_height()
                plt.text(bar.get_x() + w/2, bar.get_y() + h/2, 
                         value_format.format(h), ha="center", 
                         va="center")

def ldvis(rho, amp, k,  a, ): 
    '''Computes the visibility function of a cos(theta)**p limb-darkened disk   

    
    
    Parameters
    -------
    rho : [N] float
        [lambda] Baseline length in wavelength
    k : [1] float
        [-] Limb darkening parameter 
    amp : [1] float
        [Jy] zero baseline flux amplitude    
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
    V = ldvis(B/2./lam,k,amp,a)
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

    # rho = B/2./lam

    # val = np.abs((2*jv(1,2*np.pi*rho*a)/(2*np.pi*rho*a)))
    integral_val = np.zeros_like(rho)
    for i in range(len(rho)): 
        f = lambda r: r*np.sqrt(1-r**2/a**2)**k*jv(0,2*np.pi*rho[i]*r) 
        integral_val[i] = np.abs(2*np.pi*integrate.quad(f,0,a)[0] )


    # Normalize and scale to  
    integral_val = amp*integral_val/integral_val[0]


    return integral_val

def import_vis(pickle_vis, nu, uvd_limit = 1000,): 

    import pickle 
    from astropy import constants as cst


    p = open(pickle_vis,"rb")
    vis = pickle.load(p,encoding='latin1') 

    amp = vis['amplitude']
    uvd  = vis["uvdist"]

    # 
    uvd_t = uvd[uvd<uvd_limit] 
    amp_t = amp[0,0,uvd<uvd_limit]  

    # Sort the array 
    i_s = np.argsort(uvd_t)
    uvd_s = uvd_t[i_s]
    amp_s = amp_t[i_s]

    # Corresponding wavelength
    lam = cst.c.value/nu
    uvw_s = uvd_s/2/lam #[lambda] UVwave distance 
    B = np.arange(0,uvd_limit) # [m] reference array, uvdistance 
    rho = B/2/lam       #[lambda] reference array, uvwave 

    return uvw_s, amp_s, rho

def flux_fitting(uvw, amp, p0=None, bounds=(-np.inf,np.inf),  Plotting=False):

    
    from scipy.optimize import curve_fit

    rho = np.arange(0,uvw[-1])

    popt, pcov = curve_fit(ldvis, uvw, amp, p0=p0, method='dogbox', bounds=bounds)
    mvis = ldvis(rho,*popt)

    if Plotting: 

        plt.figure(figsize=(16,9))
        plt.plot(uvw_s,amp_s,'.',label='calibrated vis')
        plt.plot(rho,mvis,label='curve fit: flux=%5.3f, p=%5.3f, size=%5.3f' % (popt[0],popt[1],popt[2]*3600*180/np.pi) )
        plt.plot(rho,mvisp0,label='initial guess: flux=%5.3f, p=%5.3f, size=%5.3f' % (p0[0],p0[1],p0[2]*3600*180/np.pi) )
        plt.xlim(0,uvd_limit/lam/2)
        plt.ylim(0,1)
        plt.xlabel('UVwave [lambda]')
        plt.ylabel('Amplitude [Jy]')
        plt.legend()
        plt.show()

    return popt, pcov


def shortspacingobservatory(nu, uvhole, name, obs='VLA',n_ants = 30, filepath='./'): 
    """ Create a short spacing observatory for CASA simobserve 



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
    elif obs.lower().strip() == 'alma':
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

def miriad_processing(uv, niter, robust=0.5, ): 
    """Prepare the uvfits file for Miriad's Deprojection. 

    This script will go through the steps required for preparing the 
    CASA derived, uvsubbed data set. It create rotationally smeared 
    images in units of brightness temperature on the way. 


    This assumes that uvsubtraction has been done in CASA already
    
    Parameters
    -------
    uv : [-] str
        [-] Name of the CASA fits file WITHOUT extension 

    niter: [-] int
        [-] Number of clean iterations for the rotationally smeared images 



    Returns
    ----------
    script : [-] str
        [-] create a script called miriad.bsh that can be executed in 
            the folder where the uv file resides  

   

    Keywords
    ----------
    robust : [-] int 
        [-] robustness parameter for cleaning. 0.5 is the sweet spot 


    Example
    -------
    import CASAtools

    uvdata = 'jup_x' # Uvdata set is called jup_x.uvfits
    niter = 10 
    CASAtools.miriad_processing(uvdata, niter, )

    References
    ------------
    

    Notes
    -------
    02/11/2019, CM, Initial Commit

    """  
    
    filepath = 'miriad.bsh'
    with open(filepath_script,'w') as fo:
        fo.write('#!/bin/bash -xef\n\n')
        
        fo.write('# Read the uvdata into Miriad\n')  
        fo.write('fits op=uvin in={:s}.fits out=temp.uv\n'.format(uvfits))
        
        fo.write('# Add planet ephemeris data to the dataset\n')  
        fo.write('uvplanet vis=temp.uv out={:s}.uv pltb=0,0\n'.format(uvfits))

        fo.write('# Average the uvdataset in Stokes i\n')  
        fo.write('uvaver vis={:s}.uv out={:s}.uv.comp stokes=i\n'.format(uvfits,uvfits)) 
        
        fo.write('# Create rotationally averaged maps\n')
        fo.write('invert vis={:s} map=pl_r{:1.1f}.imap beam=pl.ibem options=mfs robust={:1.1f}\n'.format(uvfits + '.uv',robust,robust)) 
        fo.write('cgdisp in=pl_r{:1.1f}.imap  device=LS_r{:1.1f}.ps/vps\n'.format(robust,robust)) 
        
        fo.write('#Convert to brightness temperature\n')
        if niter <= 1: 
            fo.write('clean map=pl_r{:1.1f}.imap beam=pl.ibem out=p_dirty.imod niters={:d} gain=0.01\n'.format(robust,niter)) 
            fo.write('restor map=pl_r{:1.1f}.imap beam=pl.ibem model=pl_dirty.imod out=pl.icln \n'.format(robust))

        else: 
            fo.write('clean map=pl_r{:1.1f}.imap beam=pl.ibem out=pl_cli{:d}.imod niters={:d}\n'.format(robust,niter,niter))      
            fo.write('restor map=pl_r{:1.1f}.imap beam=pl.ibem model=pl_cli{:d}.imod out=pl.icln \n'.format(robust,niter))
        
        fo.write('# Convert to brightness temperature\n')
        fo.write('tok.pl in=pl.icln out=pl_Tb.imap\n') 

        fo.write('# Export rotationally smeared images\n') 
        fo.write('fits op=xyout in=pl_Tb.imap out={:s}_Tb.fits'.format(uvfits)) 
        fo.close()

    return 



def parrallel_Miriad_script(m_ncore, uvfits, latrange, latint, cell,
    spwn=1,
    fwhm = None,
    robust = 0.0,
    planet = 'jupiter', 
    filepath_data= './', 
    filepath_script = 'Parallel_facets.bsh', 
    tmp_directory = '/tmp'): 
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

    # Establish if averaging is happening or not 
    if spwn != 1: 
        uvfitsc = uvfits + '.comp' 
    else: 
        uvfitsc = uvfits  # Has been compressed outside 


    # Split the uvfits file into their individual components, seperate by comma, add upper directory to each 




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
        lat_upper = lat_lower + (nband-1)*dlat 
        if lat_upper > latrange/2: 
            lat_upper = latrange/2


        filepath = temp_name+'/params' + '.pl' 
        with open(filepath,'w') as fo:  
            fo.write('$planet = "{:s}";# Planet name\n'.format(planet))
            fo.write('$vis = "../{:s}";     # Visibility file\n'.format(uvfitsc))
            fo.write('$cell = {:4.4f};      # Image pixel size, in arcseconds (!).\n'.format(cell))
            fo.write('$minlat = {:2.10f};      # Min latitude (degrees) to map\n'.format(lat_lower))
            fo.write('$maxlat = {:2.10f};       # Max latitude\n'.format(lat_upper))
            fo.write('$minlon = 0;        # Min longitude (degrees)\n')
            fo.write('$maxlon = 360;      # Max longitude\n')
            fo.write('$latint = {:2.1f};      # Increment in latitude for facets, in degrees.\n'.format(latint))
            fo.write('$imsize = 150;          # Facet pixel size.\n')
            if fwhm is None:
                fo.write('#$fwhm_km = 3200;       # Optional extra smoothing function.\n')
            else:
                fo.write('$fwhm_km = {:2.0f};       # Optional extra smoothing function.\n'.format(fwhm))
            if robust != 0.0:
                fo.write('$robust = {:2.1f};          # Optional imaging weighting parameter.\n'.format(robust))
            else: 
                fo.write('$robust = 0.0;          # Optional imaging weighting parameter.\n')            
            fo.write('$plradius = 71492.0     # Optional planet radius in kilometers.\n')
            fo.write('# $obstime = ""  # Optional observation time used for geometry.\n')
            fo.close()

        lat_lower = lat_upper+dlat


        # Make a bin file that you can execute that points to all the correct points 



    # Write the main script that will launch all the subscripts 

    with open(filepath_script,'w') as fo:
        fo.write('#!/bin/bash -xef\n\n')
        if spwn != 1: 
        # Reduce the size of the uvfits file if it doesn't exist already  
            spwids = list(range(1,spwn+1))
            spwstr = ''
            for i in spwids: 
                spwstr += ' {:d}'.format(i)
            fo.write('if [ ! -e {:s} ]; then\n'.format(uvfits+'.comp'))   
            fo.write('  for i in {:s} \n'.format(spwstr))
            fo.write('  do\n')  
            fo.write('  rm -rf junk$i.uv\n') 
            fo.write('    uvaver vis={:s} out=junk$i.uv "select=win($i)" stokes=i\n'.format(uvfits))
            fo.write('  done\n')
            fo.write('  uvaver vis=junk*.uv out={:s}\n'.format(uvfits+'.comp'))
            fo.write('fi\n\n')
        # Initiate temporary directory and append relevant scripts 
        fo.write('export TMPDIR={:s}\n'.format(tmp_directory))
        bashcmd = 'echo "Script is running" \nrm -rf facets; mkdir facets '
        for i in range(ncore):
            # Reduce memory allocation problems by running them out of phase 
        # if i < ncore/2:
            # bashcmd = bashcmd + (' \n(rm -rf ~/tmp{:d} && mkdir ~/tmp{:d} && TMPDIR="~/tmp{:d}" && cd temp_p{:d} && nohup perl /usr/local/miriad/bin/darwin/facets.pl && cp facets/* ../facets/) &'.format(i))
            bashcmd = bashcmd + ('\n(cd temp_p{:d} && perl /usr/local/miriad/bin/darwin/facets.pl &> logger.txt  && cp -r facets/ ../facets/) &'.format(i))
            # bashcmd = bashcmd + ('\n(rm -rf /tmp{:s}{:d} ; mkdir /tmp{:s}{:d} && TMPDIR="{:s}{:d}" && cd temp_p{:d} && perl /usr/local/miriad/bin/darwin/facets.pl &> logger.txt  && cp -r facets/* ../facets/) &'.format(tmp_directory,i,tmp_directory,i,tmp_directory,i,i))

        # else: 
            # bashcmd = bashcmd + ('\n(cd temp_p{:d} && sleep && perl /usr/local/miriad/bin/darwin/facets.pl &> logger.txt  && cp -r facets/* ../facets/) &'.format(i))
            # bashcmd = bashcmd + ('\n(rm -rf /tmp{:s}{:d} ; mkdir /tmp{:s}{:d} && TMPDIR="{:s}{:d}" && cd temp_p{:d} && sleep && perl /usr/local/miriad/bin/darwin/facets.pl &> logger.txt  && cp -r facets/* ../facets/) &'.format(tmp_directory,i,tmp_directory,i,tmp_directory,i,i))

        bashcmd = bashcmd +('&\nmail -s "Facetting done" chris.moeckel@berkeley.edu <<< " " ')
        
        fo.write(bashcmd)
        fo.write('&&\nperl /usr/local/miriad/bin/darwin/stitch.pl')


    # Change permisson   
    os.system('chmod u+x ' + filepath_script)


    # Also write a master script and place in main folder 

    filepath = 'params' + '.pl' 
    with open(filepath,'w') as fo:  
        fo.write('$planet = "{:s}";# Planet name\n'.format(planet))
        fo.write('$vis = "{:s}";     # Visibility file\n'.format(uvfitsc))
        fo.write('$cell = {:4.4f};      # Image pixel size, in arcseconds (!).\n'.format(cell))
        fo.write('$minlat = {:2.10f};      # Min latitude (degrees) to map\n'.format(-latrange/2))
        fo.write('$maxlat = {:2.10f};       # Max latitude\n'.format(latrange/2))
        fo.write('$minlon = 0;        # Min longitude (degrees)\n')
        fo.write('$maxlon = 360;      # Max longitude\n')
        fo.write('$latint = {:2.1f};      # Increment in latitude for facets, in degrees.\n'.format(latint))
        fo.write('$imsize = 150;          # Facet pixel size.\n')
        if fwhm is None:
            fo.write('#$fwhm_km = 3200;       # Optional extra smoothing function.\n')
        else: 
            fo.write('$fwhm_km = {:2.0f};       # Optional extra smoothing function.\n'.format(fwhm))
        if robust != 0.0: 
            fo.write('$robust = {:2.1f};          # Optional imaging weighting parameter.\n'.format(robust))
        else: 
            fo.write('$robust = 0.0;          # Optional imaging weighting parameter.\n')
        fo.write('$plradius = 71492.0     # Optional planet radius in kilometers.\n')
        fo.write('# $obstime = ""  # Optional observation time used for geometry.\n')
        fo.close()

    print('Data/scripts can be found here')
    print('cd ' + os.getcwd())


    return 

# itemize in=temp_p2/facets/n50p11.icln 
