''' ltmbcfiles.py '''

import io
import numpy as np

# module that contains functions used for LTM paper

def read_bc_file(fname):

    # the file format varies, so we need to be careful
    # first type: 
    #    YYYY-MM-DD
    # second type:
    #    YYYY-MM-DD HH:MM:SS
    # third type: 
    #    ID BC_VALUE YYYY-MM-DD HH:MM:SS
    #

    # open file, replacing '-' with space 
    s = io.BytesIO(open(fname, 'rb').read().replace(b'-',b' '))

    # read with space as delimiter
    # all file types have just 1 line header
    dayTXT = np.genfromtxt(s, delimiter=' ', skip_header=1)

    # file format depends on the number of columns (after turning date into YYYY MM DD)

    if dayTXT.shape[1] == 3:
        # YYYY MM DD
        # merge the date as one number
        day = (dayTXT[:,0]*100 + dayTXT[:,1])*100 + dayTXT[:,2]
        nn = np.arange(len(dayTXT[:,0]))
        bc = nn * np.nan

    elif dayTXT.shape[1] == 4:
        # YYYY MM DD HH:MM:SS 
        # keep only the columns for y/m/d
        dayTXT = dayTXT[:,0:3]
        # merge the date as one number
        day = (dayTXT[:,0]*100 + dayTXT[:,1])*100 + dayTXT[:,2]
        nn = np.arange(len(dayTXT[:,0]))
        bc = nn * np.nan
    
    elif dayTXT.shape[1] == 6:
        # ID BC YYYY MM DD HH:MM:SS
        nn = dayTXT[:,0]
        bc = dayTXT[:,1]
        # keep only the columns for y/m/d
        dayTXT = dayTXT[:,2:5]
        # merge the date as one number
        day = (dayTXT[:,0]*100 + dayTXT[:,1])*100 + dayTXT[:,2]

    else:
        print('Unexpected number of columns in file!')
        exit(1)
        
    print('File: ', fname)
    print('   Number of days = ' + str(len(day)))
    print('   BC values = ',np.mean(bc),' +- ', np.std(bc), '  max/min=', np.max(bc), np.min(bc))

    return(day, bc)

# end of module
