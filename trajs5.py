trajs5.py
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


plt.ion()
plt.interactive(True)

# load file with all trajectories
listall = np.load('listall_002_12z_96h.npy')
trajlen = 96/0.5 + 1

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
loni=-70; lonf=-20; lati=-15; latf=25; dl=1
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
        nxySingle[:,:,ff] =nxySingle[:,:,ff] + temp
    

########################################################################
# Figures

cmap = matplotlib.cm.viridis

norm = ['byTrajs', 'byPoints', 'byMax', 'single']
for nn in [0, 1, 2, 3]: 

    tag = ['clean', 'polluted', 'other']
    for ff in [0, 1, 2]:
    
        fig = plt.figure(num=ff+10*nn, clear=True)

        # geographic coordinates for maping
        ax1 = plt.axes(projection=ccrs.PlateCarree())
        ax1.set_extent([loni, lonf, lati, latf])

    
        Cnorm = colors.BoundaryNorm([0.1, 0.2, 0.5, 1,2,5,10,20,50,100], cmap.N)
        #Cnorm = colors.BoundaryNorm([1,2,5,10,20,50,100], cmap.N)
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
    
        ax1.set_title(tag[ff] + '(N={:.0f})'.format(ntrajs[ff]) + ' ' + norm[nn])
        
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                           linewidth=1, color='gray', alpha=0.5, linestyle=':')
        

        # map
        ax1.add_feature(cfeature.COASTLINE, linewidth=1.5, linestyle='-', edgecolor='black', alpha=1)
        ax1.add_feature(cfeature.BORDERS, linewidth=1, linestyle='-', edgecolor='black', alpha=.75)
        ax1.add_feature(cfeature.STATES, linewidth=1, linestyle='-', edgecolor='gray', alpha=.5)
        plt.scatter(-58.999861, -2.144111, s=30, c='blue', marker='o') #T0a

        # export PNG
        fname = 'histplot2d_traj_' + tag[ff] + '_' + norm[nn] + '_500_96h_counts.png'
        plt.savefig(fname)
#

#



