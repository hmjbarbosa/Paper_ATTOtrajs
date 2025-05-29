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

plt.ion()
plt.interactive(True)

# load file with all trajectories
#listall = np.load('listall_002_12z_96h.npy')

# read the longer runs, with 4x a day
print('read data')
if 'listall' not in locals(): 
    #listall = np.load('listall_001_4times_168h.npy') # -2.144  -59.100    500.0 -DLon
    listall = np.load('listall_002_4times_168h.npy') # -2.144  -59.000    500.0  Normal
    #listall = np.load('listall_003_4times_168h.npy') # -2.144  -59.000    250.0  -Dz
    #listall = np.load('listall_004_4times_168h.npy') # -2.144  -59.000    750.0  +Dz
    #listall = np.load('listall_005_4times_168h.npy') # -2.244  -59.000    500.0 -DLat
    #listall = np.load('listall_006_4times_168h.npy') # -2.044  -59.000    500.0 +DLat
    #listall = np.load('listall_007_4times_168h.npy') # -2.144  -58.900    500.0 +DLon

trajlen = 168/0.5 + 1  # number of points
trajalt = listall[0,1] # initial altitude

if True:
    trajalt = 'all'
    listall = np.concatenate((np.load('listall_001_4times_168h.npy'), # -2.144  -59.100    500.0 -DLon
                              np.load('listall_002_4times_168h.npy'), # -2.144  -59.000    500.0  Normal
                              np.load('listall_003_4times_168h.npy'), # -2.144  -59.000    250.0  -Dz
                              np.load('listall_004_4times_168h.npy'), # -2.144  -59.000    750.0  +Dz
                              np.load('listall_005_4times_168h.npy'), # -2.244  -59.000    500.0 -DLat
                              np.load('listall_006_4times_168h.npy'), # -2.044  -59.000    500.0 +DLat
                              np.load('listall_007_4times_168h.npy')), # -2.144  -58.900    500.0 +DLon
                             axis = 0)

    
# grab just the columns we need here 
t0   = listall[:,0]   # start hour of traj (0, 6, 12, 18Z)
alt  = listall[:,1]   # start altitude (250, 500, 750 m)
tid  = listall[:,2]   # trajectory number
dt   = listall[:,3]   # backward time (0h, -0.5h, -1h, ...)
tf   = listall[:,4]   # time of day
x    = listall[:,5]   
y    = listall[:,6]   # position
z    = listall[:,7]
year = listall[:,13]
mon  = listall[:,14]  # date
day  = listall[:,15]
flag = listall[:,16]  # pollution flag

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


########################################################################
# trying to count each trajectory just once

nxySingle = np.zeros([len(ycenters), len(xcenters), 3])

# clean=0 polluted=1 other=2
print('compute single')
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
    

########################################################################
# Figures

cmap = cm.viridis
print('plots')
norm = ['byTrajs', 'byPoints', 'byMax', 'single']
nn = 2
tag = ['clean', 'polluted', 'other']
    
fig = plt.figure(num=123, clear=True)

# geographic coordinates for maping
ax1 = plt.axes(projection=ccrs.PlateCarree())
#ax1.set_extent([loni, lonf, lati, latf])


#Cnorm = colors.BoundaryNorm([0.1, 0.2, 0.5, 1,2,5,10,20,50,100], cmap.N)
Cnorm = colors.BoundaryNorm([1,2,5,10,20,50,100], cmap.N)
#Cnorm = colors.LogNorm()
#Cnorm = colors.Normalize()

if nn==2:
    pcCl = ax1.contourf(xcenters, ycenters, 100*nxy[:,:,0]/nmax[0], cmap=cm.Blues.reversed(), norm=Cnorm,
                        levels=[1,2,5,10,20,50,100])
    #ax1.clabel(pc, pc.levels[:-1], inline=True)

    pcPo = ax1.contour(xcenters, ycenters, 100*nxy[:,:,1]/nmax[1], cmap=cm.Reds.reversed(), norm=Cnorm,
                       levels=[1,2,5,10,20,50,100])
    ax1.clabel(pcPo, pcPo.levels[:-1], inline=True)

clean = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor=cm.Blues(0.8))
polut = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor=cm.Reds(0.8))

ax1.legend([clean, polut], ['Clean', 'Polluted'],
           loc='upper left')#, bbox_to_anchor=(0.025, -0.1), fancybox=True)
    
fig.colorbar(pcCl, orientation='horizontal',pad=0.1, aspect=30,
             label='Density of trajectories [%]', fraction=0.06, shrink=0.9)

ff=0
ax1.set_title(norm[nn] + ' h=' + str(trajalt))

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                           linewidth=1, color='gray', alpha=0.5, linestyle=':')
        

# map
ax1.add_feature(cfeature.COASTLINE, linewidth=1.5, linestyle='-', edgecolor='black', alpha=1)
ax1.add_feature(cfeature.BORDERS, linewidth=1, linestyle='-', edgecolor='black', alpha=.75)
ax1.add_feature(cfeature.STATES, linewidth=1, linestyle='-', edgecolor='gray', alpha=.5)
plt.scatter(-58.999861, -2.144111, s=30, c='blue', marker='o') #T0a
#sys.exit()

# export PNG
fname = 'histplot2d_traj_both_' + norm[nn] + '_' + str(trajalt) + '_168h_counts.png'
plt.savefig(fname)
#

#



