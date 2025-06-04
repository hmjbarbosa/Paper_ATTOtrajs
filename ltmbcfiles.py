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


def density_map_points(bounds, dbnd, coords, flag, dt):
    
    # grid limits
    loni = bounds[0]
    lonf = bounds[1]
    lati = bounds[2]
    latf = bounds[3]
    dlon = dbnd[0]
    dlat = dbnd[1]
    
    # define grid
    xedges = np.arange(loni, lonf + dlon, dlon)
    yedges = np.arange(lati, latf + dlat, dlat)
    #xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
    #ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
    print('Longitude: ', loni, ' to ', lonf, ' in ', len(xedges), ' steps of ', dlon)
    print('Latitude: ', lati, ' to ', latf, ' in ', len(yedges), ' steps of ', dlat)
    
    # matrix to store results
    nxy = np.zeros([len(yedges)-1, len(xedges)-1, 3])
    allflags = np.unique(flag).astype(int)
    ntrajs = np.zeros(len(allflags))
    nmax = np.zeros(len(allflags))
    nsum = np.zeros(len(allflags))

    # data
    x = coords[0]
    y = coords[1]
    
    # clean=0 polluted=1 other=2
    for ff in allflags: 
        mask = (flag==ff)

        # count the number of trajectory points inside each gridbox
        nxy[:,:,ff] = np.histogram2d(y[mask], x[mask], bins=[yedges, xedges])[0]
        
        # Each traj can contribute more than once for each gridpoint,
        # depending on the output time-step of the trajectory and the size
        # of the grid. Which means there are different ways of normalizing
        # the density map.
        
        # number of trajectories
        trajini = mask & (dt == 0)
        ntrajs[ff] = np.nansum(trajini)
        
        # maximum number of counts in any gridbox
        nmax[ff] = np.nanmax(nxy[:,:,ff])
        
        # total number of counts
        nsum[ff] = np.nansum(nxy[:,:,ff])
        
        #nxy[:,:,ff] = 100*nxy[:,:,ff]/nmax[ff] 
        #nxy[:,:,ff][nxy[:,:,ff]<0.01]=np.nan

    return(nxy, xedges, yedges, ntrajs, nmax, nsum)



def density_map_trajs(bounds, dbnd, coords, flag, tid):
    
    # grid limits
    loni = bounds[0]
    lonf = bounds[1]
    lati = bounds[2]
    latf = bounds[3]
    dlon = dbnd[0]
    dlat = dbnd[1]
    
    # define grid
    xedges = np.arange(loni, lonf + dlon, dlon)
    yedges = np.arange(lati, latf + dlat, dlat)
    #xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
    #ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
    print('Longitude: ', loni, ' to ', lonf, ' in ', len(xedges), ' steps of ', dlon)
    print('Latitude: ', lati, ' to ', latf, ' in ', len(yedges), ' steps of ', dlat)
    
    # matrix to store results
    nxy = np.zeros([len(yedges)-1, len(xedges)-1, 3])
    allflags = np.unique(flag).astype(int)
    ntrajs = np.zeros(len(allflags))
    nmax = np.zeros(len(allflags))
    nsum = np.zeros(len(allflags))

    # data
    x = coords[0]
    y = coords[1]

    # clean=0 polluted=1 other=2
    for ff in [0, 1, 2]: 
        mask = (flag==ff)

        # we have to build a histogram for each trajectory separately
        # first, get the list of trajectory IDs for the current pollution flag
        idlist = np.unique(tid[mask])

        # now loop over each trajectory
        for ii in idlist:
            # build histogram with x,y data from a single trajectory
            #NOTE: tid==ii will not work if reading multiple NPY together
            temp =  np.histogram2d(y[mask & (tid==ii)], x[mask & (tid==ii)], bins=[yedges, xedges])[0]

            # get a matrix of 0 (didn't pass) and 1 (trajectory passed there)
            temp[temp>0] = 1

            # finally, accumulate this counts
            nxy[:,:,ff] =nxySingle[:,:,ff] + temp
    
    return(nxy, xedges, yedges, ntrajs, nmax, nsum)
    
# end of module
