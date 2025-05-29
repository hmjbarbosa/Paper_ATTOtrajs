import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
import numpy as np
import datetime as dt
from PIL import Image
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pysplit

plt.ion()
plt.interactive(True)

# FINAL FIGURE FOR PAPER
# MAP OF SINGLE TRAJECTORY, WITH FIRE SPOTS AND AVG BOXES
# USED TO EXPLAIN THE METHOD

def bboxtraj(atraj, margin = 0.):
    bbox = []
    bbox.append(min(atraj.data['geometry'].x.values) - margin)
    bbox.append(max(atraj.data['geometry'].x.values) + margin)
    bbox.append(min(atraj.data['geometry'].y.values) - margin)
    bbox.append(max(atraj.data['geometry'].y.values) + margin)
    return(bbox)

def bboxtrajgroup(agroup, margin = 0.):
    bbox = []
    allx = np.concatenate([ p.data['geometry'].x.values[:] for p in agroup ])
    ally = np.concatenate([ p.data['geometry'].y.values[:] for p in agroup ])
    bbox.append(min(allx) - margin)
    bbox.append(max(allx) + margin)
    bbox.append(min(ally) - margin)
    bbox.append(max(ally) + margin)
    return(bbox)

# 1 - 1500m
# 2 -  100m
# 3 -  500m
# 4 - 1000m

print('reading trajs...')

trajgroup = pysplit.make_trajectorygroup('./002/t*')
ntrajs = trajgroup.trajcount
print('ntrajs = ',ntrajs)

# find bounding box for trajetories
print('bbox...')
minlon, maxlon, minlat, maxlat = -75, -45, -25, 0 
print('trajs bbox: ',minlon,maxlon,minlat,maxlat)

# read fire data
print('reading fires...')
fdataz = np.load('fire_data.npz', allow_pickle=True)
flats = fdataz['flats']
flons = fdataz['flons']
ftime = fdataz['ftime']
frp = fdataz['frp']

# exclude some points away from the interest region
fmask = ( (flats >= minlat) & (flats <= maxlat) &
          (flons >= minlon) & (flons <= maxlon) &
          (ftime >= dt.datetime(2018,9,1)) & (ftime <= dt.datetime(2018,9,8)) ) 

flats = flats[fmask]
flons = flons[fmask]
ftime = ftime[fmask]
frp = frp[fmask]

print('fires in bbox: ', len(frp), ' / ', len(fmask))

# for each trajectory
for i in range(0,ntrajs):
    traj = trajgroup[i]
    print(traj.fullpath)

    # position
    x = traj.data['geometry'].x.values
    y = traj.data['geometry'].y.values

    # time
    t = traj.data['DateTime'].values
    ts = (t - np.datetime64('1970-01-01'))/np.timedelta64(1,'s')
    td = [ dt.datetime.fromtimestamp(tmp, tz=dt.timezone.utc) for tmp in ts ]

    # bbox
    #mminlon, mmaxlon, mminlat, mmaxlat = bboxtraj(traj, 0.5)
    mminday = dt.datetime(*td[-1].timetuple()[:3])
    mmaxday = dt.datetime(*td[0].timetuple()[:3])
    #mmaxday = dt.datetime(*(td[0] + dt.timedelta(days=1)).timetuple()[:3])

    ffmask = ( #(flats >= mminlat) & (flats <= mmaxlat) &
               #(flons >= mminlon) & (flons <= mmaxlon) &
               (ftime >= mminday) & (ftime <= mmaxday) )
    
    fflats = flats[ffmask]
    fflons = flons[ffmask]
    fftime = ftime[ffmask]
    ffrp = frp[ffmask]

    print('fires around trajectory: ', len(ffrp), ' / ', len(ffmask))
    
    plt.figure(1)
    plt.clf()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-73, -56, -20, -5])

    # fires
    #plt.scatter(fflons, fflats, s=0.5, c=mdates.date2num(fftime), marker='s',
    #            cmap=plt.get_cmap('viridis',3), edgecolor='face', alpha=1)
    plt.scatter(fflons, fflats, s=1.5, c=ffrp, marker='s', vmin=10, vmax=50,
                cmap=plt.get_cmap('magma_r',16), edgecolor='face', alpha=1)
    #cb = plt.colorbar(ticks=mdates.DayLocator(interval=1), format=mdates.DateFormatter('%d'),fraction=0.040)
    #cb.set_label('Days')
    cb = plt.colorbar(fraction=0.040)
    cb.set_label('Fire Radiative Power [W/m2]')

    # traj
    plt.plot(x,y,linewidth=2)
    #plt.title(traj.fullpath)
    # boxes of interest
    for j in range(0, len(x), 12):
        plt.gca().add_patch(plt.Rectangle((x[j] - .5, y[j] - .5), 1, 1,
                                          facecolor="none",edgecolor='red'))

    # map
    ax.add_feature(cfeature.COASTLINE, linewidth=1.5, linestyle='-', edgecolor='black', alpha=.8)
    ax.add_feature(cfeature.BORDERS, linewidth=1.5, linestyle='-', edgecolor='black', alpha=.8)
    #ax.add_feature(cfeature.STATES, linewidth=1.5, linestyle='-', edgecolor='gray', alpha=.5)
    ax.add_feature(cfeature.OCEAN, facecolor='gray', zorder=0)

    plt.scatter(-67.8076, -9.974, s=30, c='m', marker='o')
    gl = ax.gridlines(draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle=':')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator([-80, -75, -70, -65, -60, -55, -50, -45, -40])
    gl.ylocator = mticker.FixedLocator([-20, -15, -10, -5, 0])
    # export PNG
    fname = 'figs/trajfire/figure1_fires_trajs_map_with_fires.png'
    #fname = 'figs/trajfire/figure1_fires_trajs_map_with_days.png'
    plt.savefig(fname,dpi=300, bbox_inches='tight')
    im = Image.open(fname).convert('RGB').convert('P', palette=Image.ADAPTIVE)
    im.save(fname, format='PNG')

    # for each point in current trajectory
    totnum = []
    totfrp = []
    for j in range(len(x)):

        # time mask
        day = dt.datetime(*td[j].timetuple()[:3])
        fffmask = ( (abs(x[j] - fflons) < 0.5) &
                    (abs(y[j] - fflats) < 0.5) &
                    (fftime == day)  )
        
        tmp = ffrp[fffmask]
        totnum.append(len(tmp))
        totfrp.append(sum(tmp))


    fig, ax = plt.subplots(num=2, clear=True)

    tdiff = (t-t[0])/np.timedelta64(1,'h')
    
    plt.grid('on')
    plt.plot(tdiff, totfrp, '-r')
    #plt.ylabel('Fire Radiative Power [W/m2]')
    plt.xlabel('Hours')
    plt.ylim(0,11000)
    plt.xlim(-72,0)

    plt.plot(tdiff, np.array(totnum)*10, '-b')
    #plt.set_ylabel('Number of Fire Spots (#)', color='blue')
    #plt.set_ylim(0,1100)
    plt.legend(['Sum FRP', 'Sum Nfires x 10'])

    #ticks=mdates.DayLocator(interval=1), 
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    
    fname = 'figs/trajfire/figure1_fires_trajs_lineplot.png'
    #plt.gcf().canvas.draw()
    #plt.gcf().canvas.flush_events()
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    im = Image.open(fname).convert('RGB').convert('P', palette=Image.ADAPTIVE)
    im.save(fname, format='PNG')

#fim
