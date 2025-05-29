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
listall = np.load('listall_001_4times_168h.npy') # -2.144  -59.100    500.0 -DLon
listall = np.load('listall_002_4times_168h.npy') # -2.144  -59.000    500.0  Normal
#listall = np.load('listall_003_4times_168h.npy') # -2.144  -59.000    250.0  -Dz
#listall = np.load('listall_004_4times_168h.npy') # -2.144  -59.000    750.0  +Dz
#listall = np.load('listall_005_4times_168h.npy') # -2.244  -59.000    500.0 -DLat
#listall = np.load('listall_006_4times_168h.npy') # -2.044  -59.000    500.0 +DLat
#listall = np.load('listall_007_4times_168h.npy') # -2.144  -58.900    500.0 +DLon

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
loni=-65; lonf=0; lati=-12; latf=23; dl=1

selid=0
print(year[tid==selid][0], mon[tid==selid][0], day[tid==selid][0])
print('hour[0] = ' + str(tf[tid==selid][0]))
print('x[0] = ' + str(x[tid==selid][0]))
print('y[0] = ' + str(y[tid==selid][0]))
print('z[0] = ' + str(z[tid==selid][0]))


########################################################################
# Figures
   
fig = plt.figure(num=99, clear=True, figsize=(4.8, 6.4))
gs = fig.add_gridspec(3,1)

# geographic coordinates for maping
ax1 = fig.add_subplot(gs[0:2,0], projection=ccrs.PlateCarree())

ax1.set_extent([loni, lonf, lati, latf])
ax1.plot(x[tid==selid], y[tid==selid], 'o-')
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                           linewidth=1, color='gray', alpha=0.5, linestyle=':')
        

# map
ax1.add_feature(cfeature.COASTLINE, linewidth=1.5, linestyle='-', edgecolor='black', alpha=1)
ax1.add_feature(cfeature.BORDERS, linewidth=1, linestyle='-', edgecolor='black', alpha=.75)
ax1.add_feature(cfeature.STATES, linewidth=1, linestyle='-', edgecolor='gray', alpha=.5)
plt.scatter(-58.999861, -2.144111, s=30, c='blue', marker='o') #T0a

ax2 = fig.add_subplot(gs[2,0])
ax2.plot(-dt[tid==selid], z[tid==selid], 'o-')
ax2.grid()
ax2.set_ylim([0, 1500])

# export PNG
fname = 'singletraj.png'
#plt.savefig(fname)
#

#



