import numpy as np
import datetime as dt
import netCDF4 as nc4
import pysplit
import dill
import os

# Read traj files from Hysplit and save as binary file to speed up
# reading by other scripts. This script also opens the fire and
# precipitation datasets and compute their values along each
# trajectory.

# To measure the execution time: 
# %run -t -p read_trajs.py

# Function to compute the bounding box of a single trajectory
def bboxtraj(atraj, margin = 0.):
    bbox = []
    bbox.append(min(atraj.data['geometry'].x.values) - margin)
    bbox.append(max(atraj.data['geometry'].x.values) + margin)
    bbox.append(min(atraj.data['geometry'].y.values) - margin)
    bbox.append(max(atraj.data['geometry'].y.values) + margin)
    return(bbox)

# Function to compute the bounding box of a trajectory group
def bboxtrajgroup(agroup, margin = 0.):
    bbox = []
    allx = np.concatenate([ p.data['geometry'].x.values[:] for p in agroup ])
    ally = np.concatenate([ p.data['geometry'].y.values[:] for p in agroup ])
    bbox.append(min(allx) - margin)
    bbox.append(max(allx) + margin)
    bbox.append(min(ally) - margin)
    bbox.append(max(ally) + margin)
    return(bbox)

# ----------------------------------
print('reading trajs...')

# 1 - 1500m
# 2 -  100m
# 3 -  500m
# 4 - 1000m
inlist = './002/t*'
#inlist = '../milena/backward0p5/003/t*'

outfile = 'hysplit_data_500_1deg.pkl'

if (os.path.isfile(outfile)):
    print('Output file exists! ', )
    exit()

trajgroup = pysplit.make_trajectorygroup(inlist)
#trajgroup = pysplit.make_trajectorygroup('./backward0p5/003/t180828_2330*')

ntrajs = trajgroup.trajcount
print('ntrajs = ',ntrajs)

# find bounding box for all trajetories
print('bbox...')
minlon, maxlon, minlat, maxlat = bboxtrajgroup(trajgroup, 0.5)
print('trajs bbox: ',minlon,maxlon,minlat,maxlat)

# ----------------------------------
# read fire data
#print('reading fires...')
#fdataz = np.load('fire_data.npz', allow_pickle=True)
#flats = fdataz['flats']
#flons = fdataz['flons']
#ftime = fdataz['ftime']
#frp = fdataz['frp']
#
## exclude points outside trajgroup bounding box 
#fmask = ( (flats >= minlat) & (flats <= maxlat) &
#          (flons >= minlon) & (flons <= maxlon) &
#          (ftime >= dt.datetime(2018,8,24)) & (ftime <= dt.datetime(2018,11,6)) ) 
#
#flats = flats[fmask]
#flons = flons[fmask]
#ftime = ftime[fmask]
#frp = frp[fmask]
#
#print('fires in bbox: ', len(frp), ' / ', len(fmask))

# ----------------------------------
# Precipitation
#print('reading precip...')
#
#pfile = 'TRMM_L3/3B42.2018_aug_to_nov.nc4'
#pdata = nc4.Dataset(pfile,'r')
#
#plats = np.ma.filled(pdata['latitude'][:].astype(float), np.nan)
#plons = np.ma.filled(pdata['longitude'][:].astype(float), np.nan)
## this is hours since 2018-08-01 00:00:00
#ptime = np.array([dt.datetime(2018,8,1) + dt.timedelta(hours=i) for i in pdata['time'][:]])
## this is mm/hour
#prec = np.ma.filled(pdata['pcp'][:].astype(float), np.nan)
#
## exclude points outside trajgroup bounding box
#print('cropping precip from size = ', np.shape(prec))
#prec = prec[:,:,(plons>=minlon)&(plons<=maxlon)]
#prec = prec[:,(plats>=minlat)&(plats<=maxlat),:]
#prec = prec[(ptime >= dt.datetime(2018,8,24)) & (ptime <= dt.datetime(2018,11,6)),:,:]
#
#plons = plons[(plons>=minlon)&(plons<=maxlon)]
#plats = plats[(plats>=minlat)&(plats<=maxlat)]
#ptime = ptime[(ptime >= dt.datetime(2018,8,24)) & (ptime <= dt.datetime(2018,11,6))]
#
#print('to size = ', np.shape(prec))

# for each trajectory
for i in range(0, 0): 
#for i in range(0,ntrajs):
    traj = trajgroup[i]
    print(traj.fullpath)

    # position
    x = traj.data['geometry'].x.values
    y = traj.data['geometry'].y.values
    #z = [ gg.z for gg in traj.data['geometry'] ]
        
    # time stamp (seconds from 1970-1-1)
    ts = (traj.data['DateTime'].values - np.datetime64('1970-01-01'))/np.timedelta64(1,'s')
    # datetime object with corrrect UTC time zone
    td = [ dt.datetime.fromtimestamp(tmp, tz=dt.timezone.utc) for tmp in ts ]

    # bbox
    mminlon, mmaxlon, mminlat, mmaxlat = bboxtraj(traj, 0.5)
    mminday = dt.datetime(*(td[0] - dt.timedelta(days=3)).timetuple()[:3])
    mmaxday = dt.datetime(*(td[0] + dt.timedelta(days=1)).timetuple()[:3])

    ffmask = ( (flats >= mminlat) & (flats <= mmaxlat) &
               (flons >= mminlon) & (flons <= mmaxlon) &
               (ftime >= mminday) & (ftime <= mmaxday) )
    
    fflats = flats[ffmask]
    fflons = flons[ffmask]
    fftime = ftime[ffmask]
    ffrp = frp[ffmask]
    
    print('fires around trajectory: ', len(ffrp), ' / ', len(ffmask))
    
    # for each point in current trajectory
    listnum = []
    listfrp = []
    listtrmm = []
    for j in range(len(x)):

        # time/location mask for fires
        # 100 x 100km box
        # same day as traj point
        day = dt.datetime(*td[j].timetuple()[:3])
        fffmask = ( (abs(x[j] - fflons) < 0.5) &
                    (abs(y[j] - fflats) < 0.5) &
                    (fftime == day) )
        
        tmp = ffrp[fffmask]
        listnum.append(len(tmp))
        listfrp.append(sum(tmp))

        # time/location mask
        # 100 x 100km box
        # trmm has 3h resolution, so we interpolate in time
        hh = td[j].hour + td[j].minute/60.
        # full 3h before and after trajpoint
        hlow  = day  + dt.timedelta(hours=np.floor(hh/3.)*3)
        hhigh = hlow + dt.timedelta(hours=3)
        # interpolation weigths
        wlow  = ( 3 + np.floor(hh/3.)*3 - hh ) / 3.
        whigh = 1. - wlow

        # we can't apply all masks at once in python
        prec2 = prec[:,:,abs(x[j] - plons) < 0.5]
        prec3 = prec2[:,abs(y[j] - plats) < 0.5,:]        
        mprec = ( prec3[ptime==hlow,:,:]*wlow +
                  prec3[ptime==hhigh,:,:]*whigh )

        # the poing may be outside the TRMM grid
        if (mprec.size == 0):
            print('Trying to read TRMM data outside the grid: ',x[j], y[j])
            listtrmm.append(np.nan)
        else:
            listtrmm.append(np.mean(mprec))

    # add new columns to traj object
    traj.data['num_fires'] = listnum
    traj.data['frp_fires'] = listfrp
    traj.data['trmm_prec'] = listtrmm

# Save the trajgroup to a binary file using DILL

print('Writting binary file...')
with open(outfile,'wb') as f:
    dill.dump(trajgroup, f)

# To open this file, do: 
#
# with open('hysplit_data_500_0p5.pkl','rb') as f:
#    trajgroup = dill.load(f)

#
