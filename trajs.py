import glob
#import matplotlib.pyplot as plt
import numpy as np
import datetime as dt

import pysplit
import dill
import io
import myhysplit

#plt.ion()
#plt.interactive(True)
#

# LIST WITH CLEAN X POLLUTTED DAYS
#daypolutTXT = np.genfromtxt('polluted_day.txt', delimiter='-', skip_header=1)
#daycleanTXT = np.genfromtxt('clean_day.txt', delimiter='-', skip_header=1)
#filt = 'listall'

s = io.BytesIO(open('polluted90_quantile.txt', 'rb').read().replace(b'-',b' '))
daypolutTXT = np.genfromtxt(s, delimiter=' ', skip_header=1)
daypolutTXT = daypolutTXT[:,2:5]

s = io.BytesIO(open('clean10_quantile.txt', 'rb').read().replace(b'-',b' '))
daycleanTXT = np.genfromtxt(s, delimiter=' ', skip_header=1)
daycleanTXT = daycleanTXT[:,2:5]
filt = 'listquantile'

# need to add the additional times of the day, because we use that
# to compare and decide if it matches or not
#
# daypolut = np.array([[np.datetime64('{:.0f}-{:02.0f}-{:02.0f}T00:00:00'.format(x[0],x[1],x[2])),
#                       np.datetime64('{:.0f}-{:02.0f}-{:02.0f}T06:00:00'.format(x[0],x[1],x[2])),
#                       np.datetime64('{:.0f}-{:02.0f}-{:02.0f}T12:00:00'.format(x[0],x[1],x[2])),
#                       np.datetime64('{:.0f}-{:02.0f}-{:02.0f}T18:00:00'.format(x[0],x[1],x[2]))] for x in daypolutTXT[:]])
# 
# dayclean = np.array([[np.datetime64('{:.0f}-{:02.0f}-{:02.0f}T00:00:00'.format(x[0],x[1],x[2])),
#                       np.datetime64('{:.0f}-{:02.0f}-{:02.0f}T06:00:00'.format(x[0],x[1],x[2])),
#                       np.datetime64('{:.0f}-{:02.0f}-{:02.0f}T12:00:00'.format(x[0],x[1],x[2])),
#                       np.datetime64('{:.0f}-{:02.0f}-{:02.0f}T18:00:00'.format(x[0],x[1],x[2]))] for x in daycleanTXT[:]])
# 

daypolut = np.array([[np.datetime64('{:.0f}-{:02.0f}-{:02.0f}T{:02.0f}:00:00'.format(x[0],x[1],x[2],t))
                      for t in np.arange(24)] for x in daypolutTXT[:]])

dayclean = np.array([[np.datetime64('{:.0f}-{:02.0f}-{:02.0f}T{:02.0f}:00:00'.format(x[0],x[1],x[2],t))
                      for t in np.arange(24)] for x in daycleanTXT[:]])

# Loop over all the different initial positions
#for tag in ['002']:
#for tag in ['001', '002', '003', '004', '005', '006', '007']: 
for tag in ['250m', '500m', '750m', '1000m', '2000m', '3000m', '4000m']: 
    print('reading trajs for tag = ' + tag)

    # where hysplit files are located
    filelist = glob.glob('/Users/hbarbosa/DATA/ATTOtrajs/backward_500_to_4000/' + tag + '/t*')
    #filelist = glob.glob('/Users/hbarbosa/DATA/ATTOtrajs/backward/' + tag + '/t*')
    #filelist = glob.glob('/Users/hbarbosa/DATA/ATTOtrajs/002/t*')

    # name of output file
    outf='listall_' + tag + '_24times_168h.npy'
    #outf='listall_' + tag + '_4times_168h.npy'
    #outf='listall.npy'

    filelist.sort()
    trajgroup = pysplit.make_trajectorygroup(filelist)
    ntrajs = trajgroup.trajcount
    print('ntrajs = ',ntrajs)

    # for each trajectory
    listall=[]
    nn=0
    for i in np.arange(ntrajs):
        traj = trajgroup[i]
        # if the mod is not a multiple of the # trajs / per day
        # the printed values will seem to be jumping strangely.
        # But they are still in the correct time order. 
        if np.mod(i, 24) == 0: 
            print(traj.fullpath)

        # position
        x = traj.data['geometry'].x.values
        y = traj.data['geometry'].y.values
        z = np.array([f.z for f in traj.data['geometry']])
        sw = traj.data['Solar_Radiation'].values
        rr = traj.data['Rainfall'].values
        
        # time
        t = traj.data['DateTime'].values
        ts = (t - np.datetime64('1970-01-01'))/np.timedelta64(1,'s')
        td = [ dt.datetime.fromtimestamp(tmp, tz=dt.timezone.utc) for tmp in ts ]
        backhours = (ts-ts[0])/3600.

        # find if traj is clean or polluted 
        if ( t[0] in dayclean ):
            flag = 0
        else:
            if ( t[0] in daypolut ):
                flag = 1
            else:
                flag = 2

        # write out the trajectory points
        for j,backh in enumerate(backhours):
            listall.append([td[0].hour,               # start hour of traj (0, 6, 12, 18Z)
                            z[0],                     # start altitude (250, 500, 750 m)
                            i,                        # trajectory number
                            backh,                    # backward time (0h, -0.5h, -1h, ...)
                            td[0].hour + backh,       # time of day
                            x[j], y[j], z[j],         # position 
                            sw[j],                    # SW instantaneous 
                            np.sum(sw[0:j+1]),        # SW integrated from traj begin
                            np.sum(sw[0:j+1]>100.),   # number of daylight timesteps
                            rr[j],                    # Rain instantaneous 
                            np.sum(rr[0:j+1]),        # Rain integrated from traj begin
                            td[0].year,               # year
                            td[0].month, td[0].day,   # date
                            flag] )                   # pollution flag

    # save data for current initial position
    np.save(outf, np.array(listall))


#fim
