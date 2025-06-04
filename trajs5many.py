import cartopy.crs as ccrs
import cartopy.feature as cfeature
#import glob
import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.ticker as mticker
import matplotlib.colors as colors
#from matplotlib.dates import DateFormatter
from matplotlib import cm
#import matplotlib.dates as mdates
import numpy as np
import sys
import io

plt.ion()
plt.interactive(True)

# ========================================================================================
# LOAD FILE WITH ALL TRAJECTORIES
# ========================================================================================

# skip reading if variable already in memory
if 'listall' not in locals(): 
    print('read data')

    # FIRST RUN
    #listall = np.load('run1/listall_002_12z_96h.npy')

    # SECOND RUN
    #listall = np.load('run2/listall_001_4times_168h.npy') # -2.144  -59.100    500.0 -DLon
    #listall = np.load('run2/listall_002_4times_168h.npy') # -2.144  -59.000    500.0  Normal
    #listall = np.load('run2/listall_003_4times_168h.npy') # -2.144  -59.000    250.0  -Dz
    #listall = np.load('run2/listall_004_4times_168h.npy') # -2.144  -59.000    750.0  +Dz
    #listall = np.load('run2/listall_005_4times_168h.npy') # -2.244  -59.000    500.0 -DLat
    #listall = np.load('run2/listall_006_4times_168h.npy') # -2.044  -59.000    500.0 +DLat
    #listall = np.load('run2/listall_007_4times_168h.npy') # -2.144  -58.900    500.0 +DLon

    # THIRD RUN
    #listall = np.load('run3/listall_250m_24times_168h.npy') # -2.144  -59.100  
    #listall = np.load('run3/listall_500m_24times_168h.npy') # -2.144  -59.000  
    #listall = np.load('run3/listall_750m_24times_168h.npy') # -2.144  -59.000  
    #listall = np.load('run3/listall_1000m_24times_168h.npy') # -2.144  -59.000 
    listall = np.load('run3/listall_2000m_24times_168h.npy') # -2.244  -59.000 
    #listall = np.load('run3/listall_3000m_24times_168h.npy') # -2.044  -59.000 
    #listall = np.load('run3/listall_4000m_24times_168h.npy') # -2.144  -58.900 
    
    trajlen = 168/0.5 + 1  # number of points
    trajalt = listall[0,1] # initial altitude

    # SECOND RUN - ALL AT ONCE
    #listall = np.concatenate((np.load('run2/listall_001_4times_168h.npy'), # -2.144  -59.100    500.0 -DLon
    #                          np.load('run2/listall_002_4times_168h.npy'), # -2.144  -59.000    500.0  Normal
    #                          np.load('run2/listall_003_4times_168h.npy'), # -2.144  -59.000    250.0  -Dz
    #                          np.load('run2/listall_004_4times_168h.npy'), # -2.144  -59.000    750.0  +Dz
    #                          np.load('run2/listall_005_4times_168h.npy'), # -2.244  -59.000    500.0 -DLat
    #                          np.load('run2/listall_006_4times_168h.npy'), # -2.044  -59.000    500.0 +DLat
    #                          np.load('run2/listall_007_4times_168h.npy')), # -2.144  -58.900    500.0 +DLon
    #                         axis = 0)

    # THIRD RUN - ALL AT ONCE
    #listall = np.concatenate((np.load('run3/listall_250m_24times_168h.npy'), # -2.144  -59.000    250.0 -Dz  
    #                          np.load('run3/listall_500m_24times_168h.npy'), # -2.144  -59.000    500.0 Normal
    #                          np.load('run3/listall_750m_24times_168h.npy'), # -2.144  -59.000    750.0 +Dz
    #                          np.load('run3/listall_1000m_24times_168h.npy'), # -2.144  -59.000  1000.0 ++Dz
    #                          np.load('run3/listall_2000m_24times_168h.npy'), # -2.244  -59.000  2000.0 +++Dz
    #                          np.load('run3/listall_3000m_24times_168h.npy'), # -2.044  -59.000  3000.0 ++++Dz
    #                          np.load('run3/listall_4000m_24times_168h.npy')), # -2.144  -58.900 4000.0 +++++Dz
    #                         axis = 0)
    #
    #trajlen = 168/0.5 + 1  # number of points
    #trajalt = 'all'       # initial altitude

    
# grab just the columns we need here 
t0   = listall[:,0]   # start hour of traj (0, 6, 12, 18Z)
alt  = listall[:,1]   # start altitude (250, 500, 750 m)
tid  = listall[:,2]   # trajectory number
dt   = listall[:,3]   # backward time (0h, -0.5h, -1h, ...)
tf   = listall[:,4]   # time of day
x    = listall[:,5]   
y    = listall[:,6]   # position
z    = listall[:,7]
date = (listall[:,13]*100 + listall[:,14])*100 + listall[:,15]

# ========================================================================================
# LIST WITH CLEAN X POLLUTTED DAYS
# ========================================================================================

#daypolutTXT = np.genfromtxt('polluted_day.txt', delimiter='-', skip_header=1)
#daycleanTXT = np.genfromtxt('clean_day.txt', delimiter='-', skip_header=1)

# open file, replacing '-' with space 
s = io.BytesIO(open('data/BCe_high.txt', 'rb').read().replace(b'-',b' '))
# read with space as delimiter
daypolutTXT = np.genfromtxt(s, delimiter=' ', skip_header=1)
# keep only the columns for y/m/d
highbc = daypolutTXT[:,1]
daypolutTXT = daypolutTXT[:,2:5]
# merge the date as one number
daypolut = (daypolutTXT[:,0]*100 + daypolutTXT[:,1])*100 + daypolutTXT[:,2]
print('List of polluted days = ' + str(len(daypolut)))
print('    BCe = ',np.mean(highbc),' +- ', np.std(highbc), '  max/min=', np.max(highbc), np.min(highbc))

s = io.BytesIO(open('data/BCe_low.txt', 'rb').read().replace(b'-',b' '))
daycleanTXT = np.genfromtxt(s, delimiter=' ', skip_header=1)
lowbc = daycleanTXT[:,1]
daycleanTXT = daycleanTXT[:,2:5]
dayclean = (daycleanTXT[:,0]*100 + daycleanTXT[:,1])*100 + daycleanTXT[:,2]
print('List of clean days = ' + str(len(dayclean)))
print('    BCe = ',np.mean(lowbc),' +- ', np.std(lowbc), '  max/min=', np.max(lowbc), np.min(lowbc))

# ========================================================================================
# CREATE FLAG OF CLEAN / POLUT FOR EACH TRAJECTORY POINT
# ========================================================================================
print('create clean / polut flag')
flag = 2 * np.ones(date.shape)      # flag for other days
flag[ np.in1d(date, daypolut, assume_unique=True) & (t0>-1) ] = 1 # pollution flag
flag[ np.in1d(date, dayclean, assume_unique=True) & (t0>-1) ] = 0 # clean flag

# ========================================================================================
# DENSITY MAP
# ========================================================================================

# define grid/histogram
# old domain for 4 days
#loni=-70; lonf=-20; lati=-15; latf=25; dl=1
# larger domain for 7 days
loni=-70; lonf=10; lati=-15; latf=35; dl=1

xedges = np.arange(loni, lonf + dl, dl)
yedges = np.arange(lati, latf + dl, dl)
xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])

# matrix to store results
nxy = np.zeros([len(ycenters), len(xcenters), 3])
ntrajs = np.zeros(3)
nmax = np.zeros(3)
nsum = np.zeros(3)

# clean=0 polluted=1 other=2
print('compute density')
for ff in [0, 1, 2]: 
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


# ========================================================================================
# trying to count each trajectory just once
# ========================================================================================

#nxySingle = np.zeros([len(ycenters), len(xcenters), 3])
#
## clean=0 polluted=1 other=2
#print('compute single')
#for ff in [0, 1, 2]: 
#    mask = (flag==ff)
#
#    # we have to build a histogram for each trajectory separately
#    # first, get the list of trajectory IDs for the current pollution flag
#    idlist = np.unique(tid[mask])
#
#    # now loop over each trajectory
#    for ii in idlist:
#        # build histogram with x,y data from a single trajectory
#NOTE: tid==ii will not work if reading multiple NPY together
#        temp =  np.histogram2d(y[mask & (tid==ii)], x[mask & (tid==ii)], bins=[yedges, xedges])[0]
#
#        # get a matrix of 0 (didn't pass) and 1 (trajectory passed there)
#        temp[temp>0] = 1
#
#        # finally, accumulate this counts
#        nxySingle[:,:,ff] =nxySingle[:,:,ff] + temp
    

# ========================================================================================
# Figures
# ========================================================================================

cmap = matplotlib.cm.viridis
print('plots')
norm = ['byTrajs', 'byPoints', 'byMax', 'single']
for nn in [2]: 

    tag = ['clean', 'polluted', 'other']
    for ff in [0, 1, 2]:
    
        fig = plt.figure(num=ff+10*nn, clear=True)

        # geographic coordinates for maping
        ax1 = plt.axes(projection=ccrs.PlateCarree())
        #ax1.set_extent([loni, lonf, lati, latf])


        #Cnorm = colors.BoundaryNorm([0.1, 0.2, 0.5, 1,2,5,10,20,50,100], cmap.N)
        Cnorm = colors.BoundaryNorm([1,2,5,10,20,50,100], cmap.N)
        #Cnorm = colors.LogNorm()
        #Cnorm = colors.Normalize()
        if nn==0:
            pc = ax1.pcolormesh(xedges, yedges, 100*nxy[:,:,ff]/ntrajs[ff], cmap=cmap, norm=Cnorm)
        if nn==1:
            pc = ax1.pcolormesh(xedges, yedges, 100*nxy[:,:,ff]/nsum[ff], cmap=cmap, norm=Cnorm)
        if nn==2:
            pc = ax1.pcolormesh(xedges, yedges, 100*nxy[:,:,ff]/nmax[ff], cmap=cmap, norm=Cnorm)
        if nn==3: 
            pc = ax1.pcolormesh(xedges, yedges, 100*nxySingle[:,:,ff]/ntrajs[ff], cmap=cmap, norm=Cnorm)

        fig.colorbar(pc, orientation='horizontal',pad=0.1, aspect=30,
                     label='Density of trajectories [%]', fraction=0.06, shrink=0.9)

        ax1.set_title(tag[ff] + '(N={:.0f})'.format(ntrajs[ff]) + ' ' + norm[nn] + ' h=' + str(trajalt))

        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                           linewidth=1, color='gray', alpha=0.5, linestyle=':')
        

        # map
        ax1.add_feature(cfeature.COASTLINE, linewidth=1.5, linestyle='-', edgecolor='white', alpha=1)
        ax1.add_feature(cfeature.BORDERS, linewidth=1, linestyle='-', edgecolor='white', alpha=.75)
        ax1.add_feature(cfeature.STATES, linewidth=1, linestyle='-', edgecolor='gray', alpha=.5)
        plt.scatter(-58.999861, -2.144111, s=30, c='blue', marker='o') #T0a
        #sys.exit()

        # export PNG
        fname = 'histplot2d_traj_' + tag[ff] + '_' + norm[nn] + '_' + str(trajalt) + '_168h_counts.png'
        plt.savefig(fname)
#

#



