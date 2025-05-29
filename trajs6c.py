import cartopy.crs as ccrs
import cartopy.feature as cfeature
#import glob
import matplotlib
#matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.dates import DateFormatter
from matplotlib import cm
import matplotlib.dates as mdates
import numpy as np
#import datetime as dt
#import netCDF4 as nc4

#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#import pysplit
#import dill
import sys

plt.ion()
plt.interactive(True)

# arquivo com todos os pontos
listall = np.load('listall.npy')

# dados do arquivo com as trajetorias
t0   = listall[:,0]
alt  = listall[:,1]
pos  = listall[:,2]
dt   = listall[:,3]
tf   = np.mod(listall[:,4], 24)
x    = listall[:,5]
y    = listall[:,6]
z    = listall[:,7]
sw   = listall[:,8]
swh  = listall[:,9]
dayh = listall[:,10]
arrTi = 14
arrTf = 18

# define a grade/histograma
loni=-62.7; lonf=-59; lati=-4.8; latf=-2.2; dl=0.2
xedges = np.arange(loni, lonf + dl, dl)
yedges = np.arange(lati, latf + dl, dl)
xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])

# start time has to be before 14-18 (maybe previous day)
t0 = tf - dt

# selecao dos dados
mask = (alt>=100) & (alt<=500) & (tf>=arrTi) & (tf<=arrTf)
nxy = np.histogram2d(y[mask], x[mask], bins=[yedges, xedges])[0]
nxy[nxy<5]=np.nan
wxy = np.histogram2d(y[mask], x[mask], bins=[yedges, xedges],weights=dayh[mask])[0]
meanval = wxy/nxy

# figura
fig = plt.figure(1)
plt.clf()
ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.set_extent([loni, lonf, lati, latf])

cmap = matplotlib.cm.viridis
norm = matplotlib.colors.BoundaryNorm(np.arange(0,10,1), cmap.N)

pc = ax1.pcolormesh(xedges, yedges, meanval, cmap=cmap, norm=norm) 
fig.colorbar(pc, orientation='horizontal',pad=0.1,aspect=30,
             label='Day-time travel [h]',fraction=0.06,shrink=0.9)

#CS = ax1.contour(xcenters, ycenters, nxy,colors='red', levels=[50,100,200,400])#np.arange(0,10,1))
#ax1.clabel(CS, CS.levels, inline=True, fmt='%1.0f', fontsize=10)

#ax1.set_title('100-500m, Anywhere '+str(arrTi)+'-'+str(arrTf)+'h UTC ('+str(arrTi-4)+'-'+str(arrTf-4)+'h LT)')

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=1, color='gray', alpha=0.5, linestyle=':')
        
# boxes of interest
minlon, maxlon, minlat, maxlat = -61.99, -60.65, -3.72, -3.06
plt.gca().add_patch(plt.Rectangle((minlon, minlat), maxlon-minlon, maxlat-minlat,
                                  facecolor="none",edgecolor='cyan'))    

# map
ax1.add_feature(cfeature.COASTLINE, linewidth=1.5, linestyle='-', edgecolor='black', alpha=1)
ax1.add_feature(cfeature.BORDERS, linewidth=1, linestyle='-', edgecolor='black', alpha=.75)
ax1.add_feature(cfeature.STATES, linewidth=1, linestyle='-', edgecolor='gray', alpha=.5)
ax1.add_feature(cfeature.RIVERS, linewidth=1, linestyle='-', edgecolor='blue', alpha=.5)
ax1.add_feature(cfeature.LAKES, linewidth=1, linestyle='-', edgecolor='green', alpha=.5)
plt.scatter(-60.598056, -3.212972, s=30, c='blue', marker='o')#T3
plt.scatter(-60.131530, -3.139200, s=30, c='blue', marker='o')#T2
plt.scatter(-59.986700, -3.097220, s=30, c='blue', marker='o')#T1


# export PNG
fname = 'figs_time/histplot2d_traj_all_100_500_24h_dayh.png'
plt.savefig(fname)
#
# verifica a distribuicao em um certo ponto da grade
#mask = mask & (x>=-61.5) & (x<=-61.)  & (y>=-3.6) & (y<=-3.3)
#mask = mask & (x>=-62.5) & (x<=-62.)  & (y>=-4.8) & (y<=-4.5)
#mask = mask & (x>=-62.5) & (x<=-62.)  & (y>=-4.5) & (y<=-4.2)
#mask = mask & (x>=-62.5) & (x<=-62.)  & (y>=-4.2) & (y<=-3.9)
#mask = mask & (x>=-62.5) & (x<=-62.)  & (y>=-3.9) & (y<=-3.6)
#mask = mask & (x>=-62.5) & (x<=-62.)  & (y>=-3.6) & (y<=-3.3)
#mask = mask & (x>=-62.0) & (x<=-61.5) & (y>=-3.6) & (y<=-3.3)
mask = mask & (x>=-61.5) & (x<=-61.0) & (y>=-3.6) & (y<=-3.3)
#mask = mask & (x>=-61.0) & (x<=-60.5) & (y>=-3.6) & (y<=-3.3)
#mask = mask & (x>=-60.5) & (x<=-60.0) & (y>=-3.6) & (y<=-3.3)
plt.plot(np.mean(x[mask]), np.mean(y[mask]), 'om')
plt.figure(2)
plt.clf()
plt.hist(dayh[mask],bins=np.arange(-0.5,24.5,1))

#
