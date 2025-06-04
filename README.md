
Scripts and figures for the ATTO paper where we looked at the
meteorological conditions when BC was high and low. This started
as a discussion with LT Machado, with preliminar analysis compiled
here:

https://docs.google.com/presentation/d/1RzNtMT_wlOHpdUKBcgnZYycCBjDjOqKCmfdy508jU9o/edit?usp=sharing

# Classification

Initial list of dates based on BC, 2015 to 2022, but only Jan, Feb,
Mar, Dec. Sent by LTM on Nov-27-2023: 

* data/clean_day.txt has 161 clean days
* data/polluted_day.txt has 110 polluted days

Second list of dates, only lowest 10% BC days (clean) and top 10% BC
days (polluted). Sent LTM on Mar-22-2024:

* data/clean10_quantile.txt has 90 clean days (date column is wrong!)
* data/polluted90_quantile.txt has 90 polluted days

Third list of dates. Sent LTM on Apr-18-2024:

* data/BCe_low.txt has 90 clean days (date column fixed)
* data/BCe_high.txt has 90 polluted days (no change)

# Hysplit

## Initial setup: 

Hysplit 4, GDAS 1deg winds, ATTO coordinates, 500m agl, 96h backward
from 12Z of each day.  All days in DJFM, 2015-2022 (many more than
just the "case" days).

Files:
* run1/listall_002_12z_96h.npy

## Second run: 

Like the previous, but 168h (7 days) backward from 0, 6, 12 and 18Z of
each day.

I also varied the initial position (6 additional locations): 
* 250, 500 and 750 m agl (up an down variation)
* +- 0.1deg in lat (east-west variation)
* +-0.1deg in lon (north-south variation)

Files:
* run2/listall_001_4times_168h.npy 500m, ATTO-DLon
* run2/listall_002_4times_168h.npy 500m, ATTO Normal
* run2/listall_003_4times_168h.npy 250m, ATTO-Dz
* run2/listall_004_4times_168h.npy 750m, ATTO+Dz
* run2/listall_005_4times_168h.npy 500m, ATTO-DLat
* run2/listall_006_4times_168h.npy 500m, ATTO+DLat
* run2/listall_007_4times_168h.npy 500m, ATTO+DLon

## Third run:

Like previous, but starting trajectories every 1h (6x more
trajectories).  Also, added height 1, 2, 3, and 4km above ATTO and we
didn't simulatd +-Lat or +-Lon.

Files:
* run3/listall_250m_24times_168h.npy ,  250m, ATTO-Dz
* run3/listall_500m_24times_168h.npy ,  500m, ATTO Normal
* run3/listall_750m_24times_168h.npy ,  750m, ATTO+Dz
* run3/listall_1000m_24times_168h.npy, 1000m, ATTO++Dz
* run3/listall_2000m_24times_168h.npy, 2000m, ATTO+++Dz
* run3/listall_3000m_24times_168h.npy, 3000m, ATTO++++Dz
* run3/listall_4000m_24times_168h.npy, 4000m, ATTO+++++Dz


# Calculating trajectory density maps

Each trajectory is composed of points connected by lines. We ran
Hysplit with an output interval of 30min, so every 1/2 hour we have a
point, and each 4 day back-trajectory has a total of 193 points.

One can define a "trajectory density map" in different ways. 

1) Each trajectory count 1 time to each gridbox it crossed. In this
case you either show the # of trajectories, or divide by the total #
of trajs, and show as %. 

Note: if a trajectory stays over some location for a long or short
time, it doesn't matter. It will count just once.

2) We count the number of "points" inside each gridbox. Hence the map
will account for the "fly time" of each trajectory in each
location. This map can be normalized by:

* The total # of trajectories (note: the max value on the map may be >> 100%)
* The total # of points (note: the max value will be << 100%)
* The maximum counts on any gridbox on the map (note: the max value will be exactly 100%). 

NOTE: The normalization will change the values, but the map should be the same. 

I compared counting each trajectory just once (1) x counting the flying time (2),
and there isn't a big difference.

