# This script reads all trajectories and write only a few relevant
# columns of information in a numpy array, which is saved to disk.

import glob
#import matplotlib.pyplot as plt
import numpy as np
import datetime as dt

import pysplit
import dill
import io
import myhysplit

# Information about the data to be processed

# First attempt, 12Z only 500m at ATTO, 96h back
#basedir = '/Users/hbarbosa/DATA/ATTOtrajs/'
#tags = ['002']
#endtag = '12z_96h'

# Every 6h for 168h back, ATTO coords @ 500m and the 6 points around it
#basedir = '/Users/hbarbosa/DATA/ATTOtrajs/backward_500m_variations/'
#tags = ['001', '002', '003', '004', '005', '006', '007']
#endtag = '4times_168h'

# Every 1h for 168h back, 7 different altitudes, ATTO coordinates
basedir = '/Users/hbarbosa/DATA/ATTOtrajs/backward_500_to_4000/'
tags = ['250m', '500m', '750m', '1000m', '2000m', '3000m', '4000m']
endtag = '24times_168h'

# Loop over all the different initial positions
for tag in tags: 
    print('reading trajs for tag = ' + tag)

    # where hysplit files are located
    filelist = glob.glob(basedir + tag + '/t*')
    print('directory = ' + basedir + tag)

    # name of output file
    outf='listall_' + tag + '_' + endtag + '.npy'
    print('output filename = ' + outf)

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
                            td[0].month, td[0].day    # date
                            ])                   

    # save data for current initial position
    np.save(outf, np.array(listall))


#fim
