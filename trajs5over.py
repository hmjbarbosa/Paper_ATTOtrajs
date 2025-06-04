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

import ltmbcfiles 

plt.ion()
plt.interactive(True)

# ========================================================================================
# LOAD FILE WITH ALL TRAJECTORIES
# ========================================================================================

# skip reading if variable already in memory
if 'listall' not in locals(): 
    print('=== READ DATA')
    
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
    #
    #trajlen = 168/0.5 + 1  # number of points
    #trajalt = listall[0,1] # initial altitude

    # SECOND RUN - ALL AT ONCE
    listall = np.concatenate((np.load('run2/listall_001_4times_168h.npy'), # -2.144  -59.100    500.0 -DLon
                              np.load('run2/listall_002_4times_168h.npy'), # -2.144  -59.000    500.0  Normal
                              np.load('run2/listall_003_4times_168h.npy'), # -2.144  -59.000    250.0  -Dz
                              np.load('run2/listall_004_4times_168h.npy'), # -2.144  -59.000    750.0  +Dz
                              np.load('run2/listall_005_4times_168h.npy'), # -2.244  -59.000    500.0 -DLat
                              np.load('run2/listall_006_4times_168h.npy'), # -2.044  -59.000    500.0 +DLat
                              np.load('run2/listall_007_4times_168h.npy')), # -2.144  -58.900    500.0 +DLon
                             axis = 0)
    #
    trajlen = 168/0.5 + 1  # number of points
    trajalt = 'all'
    
    
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
#flag = listall[:,16]  # pollution flag

# HMJB: 29-May-2025
#
# This script was trying to read column 17, but the files only have 16 columns.
#
# I found that the script trajs.py was creating .NPY files with this additional column,
# but it was using the file clean10_quantile.txt, which had a bug (wrong dates).
# Moreover, I could not find the .NPY files with this additional column.
#
# The figure created by this script overlays the clean and polluted trajectories
# on the same figure. The figure is on the slide deck, and LTM used it in the paper
# I'm not sure how I produced this figure without the NPY files... or the wrong dates
# for the clean days. 
#
# Below, I've added code to read the correct FLAGS files, and generate the figure for
# the paper, but the results is slightly different from what is on the slide deck. 
#

# ========================================================================================
# LIST WITH CLEAN X POLLUTTED DAYS
# ========================================================================================

# #daypolutTXT = np.genfromtxt('polluted_day.txt', delimiter='-', skip_header=1)
# #daycleanTXT = np.genfromtxt('clean_day.txt', delimiter='-', skip_header=1)
# 
# # open file, replacing '-' with space 
# #s = io.BytesIO(open('data/polluted90_quantile.txt', 'rb').read().replace(b'-',b' '))
# #s = io.BytesIO(open('data/BCe_high.txt', 'rb').read().replace(b'-',b' '))
# s = io.BytesIO(open('data/polluted_day_janfev.txt', 'rb').read().replace(b'-',b' '))
# # read with space as delimiter
# daypolutTXT = np.genfromtxt(s, delimiter=' ', skip_header=1)
# # keep only the columns for y/m/d
# highbc = daypolutTXT[:,1]
# daypolutTXT = daypolutTXT[:,2:5]
# # merge the date as one number
# daypolut = (daypolutTXT[:,0]*100 + daypolutTXT[:,1])*100 + daypolutTXT[:,2]
# print('List of polluted days = ' + str(len(daypolut)))
# print('    BCe = ',np.mean(highbc),' +- ', np.std(highbc), '  max/min=', np.max(highbc), np.min(highbc))
# 
# #s = io.BytesIO(open('data/clean10_quantile.txt', 'rb').read().replace(b'-',b' '))
# #s = io.BytesIO(open('data/BCe_low.txt', 'rb').read().replace(b'-',b' '))
# s = io.BytesIO(open('data/clean_day_janfev.txt', 'rb').read().replace(b'-',b' '))
# daycleanTXT = np.genfromtxt(s, delimiter=' ', skip_header=1)
# lowbc = daycleanTXT[:,1]
# daycleanTXT = daycleanTXT[:,2:5]
# dayclean = (daycleanTXT[:,0]*100 + daycleanTXT[:,1])*100 + daycleanTXT[:,2]
# print('List of clean days = ' + str(len(dayclean)))
# print('    BCe = ',np.mean(lowbc),' +- ', np.std(lowbc), '  max/min=', np.max(lowbc), np.min(lowbc))

# 1st classification
#daypolut, highbc = ltmbcfiles.read_bc_file('data/polluted_day.txt')
#dayclean, lowbc = ltmbcfiles.read_bc_file('data/clean_day.txt')

# 2nd classification
# there was an error in 'clean10' where the dates were equal to the polluted90 !!!
#daypolut, highbc = ltmbcfiles.read_bc_file('data/polluted90_quantile.txt')
#dayclean, lowbc = ltmbcfiles.read_bc_file('data/clean10_quantile.txt')

# 3rd classification
#daypolut, highbc = ltmbcfiles.read_bc_file('data/BCe_high.txt')
#dayclean, lowbc = ltmbcfiles.read_bc_file('data/BCe_low.txt')


# 4th classification
daypolut, highbc = ltmbcfiles.read_bc_file('data/polluted_day_janfev.txt')
dayclean, lowbc = ltmbcfiles.read_bc_file('data/clean_day_janfev.txt')


# ========================================================================================
# CREATE FLAG OF CLEAN / POLUT FOR EACH TRAJECTORY POINT
# ========================================================================================
print('=== CREATE CLEAN / POLUT FLAG')
flag = 2 * np.ones(date.shape)      # flag for other days
flag[ np.in1d(date, daypolut, assume_unique=True) & (t0>-1) ] = 1 # pollution flag
flag[ np.in1d(date, dayclean, assume_unique=True) & (t0>-1) ] = 0 # clean flag
print('Clean = ', np.sum(flag==0))
print('Polut = ', np.sum(flag==1))
print('Other = ', np.sum(flag==2))

# ========================================================================================
# DENSITY MAP
# ========================================================================================

print('=== COMPUTE DENSITY')

# define grid/histogram
# old domain for 4 days
#loni=-70; lonf=-20; lati=-15; latf=25; dl=1
# larger domain for 7 days
loni=-70; lonf=10; lati=-15; latf=35; dl=1

nxy, xedges, yedges, ntrajs, nmax, nsum = ltmbcfiles.density_map_points(
    [loni, lonf, lati, latf], [dl, dl], [x, y], flag, dt)

#xedges = np.arange(loni, lonf + dl, dl)
#yedges = np.arange(lati, latf + dl, dl)
xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
#
## matrix to store results
#nxy = np.zeros([len(ycenters), len(xcenters), 3])
#ntrajs = np.zeros(3)
#nmax = np.zeros(3)
#nsum = np.zeros(3)
#
## clean=0 polluted=1 other=2
#for ff in [0, 1, 2]: 
#    mask = (flag==ff)
#
#    # count the number of trajectory points inside each gridbox
#    nxy[:,:,ff] = np.histogram2d(y[mask], x[mask], bins=[yedges, xedges])[0]
#    
#    # Each traj can contribute more than once for each gridpoint,
#    # depending on the output time-step of the trajectory and the size
#    # of the grid. Which means there are different ways of normalizing
#    # the density map.
#    
#    # number of trajectories
#    trajini = mask & (dt == 0)
#    ntrajs[ff] = np.nansum(trajini)
#
#    # maximum number of counts in any gridbox
#    nmax[ff] = np.nanmax(nxy[:,:,ff])
#
#    # total number of counts
#    nsum[ff] = np.nansum(nxy[:,:,ff])
#    
#    #nxy[:,:,ff] = 100*nxy[:,:,ff]/nmax[ff] 
#    #nxy[:,:,ff][nxy[:,:,ff]<0.01]=np.nan
#

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
print('=== PLOTS')
norm = ['byTrajs', 'byPoints', 'byMax', 'single']

def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"

for nn in [2]:

    tag = ['clean', 'polluted', 'other']
    
    fig = plt.figure(num=123, clear=True)

    # geographic coordinates for maping
    ax1 = plt.axes(projection=ccrs.PlateCarree())
    #ax1.set_extent([loni, lonf, lati, latf])

    
    #Cnorm = colors.BoundaryNorm([0.1, 0.2, 0.5, 1,2,5,10,20,50,100], cmap.N)
    #Cnorm = colors.BoundaryNorm([1,2,5,10,20,50,100], cmap.N)
    #cl = [0.1, 0.2, 0.5, 1,2,5,10,20,50,100]
    cl = [1,2,5,10,20,50,100]
    #cl = [1e-5,5,10,15,20,25]
    Cnorm = colors.BoundaryNorm(cl, cmap.N)
    #Cnorm = colors.LogNorm()
    #Cnorm = colors.Normalize()

    if nn==2:
        # revert the color map with =cm.Blues.reversed()
        
        pcCl = ax1.contourf(xcenters, ycenters, 100*nxy[:,:,0]/nmax[0], cmap=cm.Blues.reversed(), norm=Cnorm,
                            levels=cl)
        ax1.clabel(pcCl, pcCl.levels[:-1], inline=True, fmt=fmt, fontsize=12, colors='k')
        
        pcPo = ax1.contour(xcenters, ycenters, 100*nxy[:,:,1]/nmax[1], cmap=cm.Reds.reversed(), norm=Cnorm,
                           levels=cl)
        ax1.clabel(pcPo, pcPo.levels[:-1], inline=True, fmt=fmt, fontsize=12)

    # make fake rectangle. we will use the colors for the legend
    clean = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor=cm.Blues(0.8))
    polut = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor=cm.Reds(0.8))

    # otherwise the legend will use the first color in the 'shading', which is too light
    ax1.legend([clean, polut], ['Clean', 'Polluted'],
               loc='upper left')#, bbox_to_anchor=(0.025, -0.1), fancybox=True)

    # place the colocar bar
    fig.colorbar(pcCl, orientation='horizontal',pad=0.1, aspect=30,
                 label='Density of trajectories [%]', fraction=0.06, shrink=0.9)

    #ax1.set_title(norm[nn] + ' h=' + str(trajalt))

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



