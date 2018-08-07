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
            if s.lower() == 'y' or s.lower() == 'yes': 
                print('File will be overwritten!')
                os.system('rm -rf' + output) 
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
