
Scripts and figures for the ATTO paper where we looked at the
meteorological conditions when BC was high and low. This started
as a discussion with LT Machado, with preliminar anallysis compiled
here:

https://docs.google.com/presentation/d/1RzNtMT_wlOHpdUKBcgnZYycCBjDjOqKCmfdy508jU9o/edit?usp=sharing

# Trajectories arriving at ATTO

List of dates based on BC, 2015 to 2022, but only Jan, Feb, Mar, Dec
* 161 clean days
* 110 polluted days

Initial setup: 

Hysplit 4, GDAS 1deg winds, ATTO coordinates, 500m agl, 96h backward
from 12Z of each day.  All days in DJFM, 2015-2022 (many more than
just the "case" days).

## Trying different way of calculating it

Each trajectory is composed of points connected by lines. We ran
Hysplit with an output interval of 30min, so every ½ hour we have a
point, and each 4 day back-trajectory has a total of 193 points.

One can define a "trajectory density map" in different ways. 

1) Each trajectory count 1 time to each gridbox it crossed. In this
case you either show the # of trajectories, or divide by the total #
of trajs, and show as %

Note: if a trajectory stays over some location for a long or short
time, it doesn't matter. It will count just once.

2) We count the number of "points" inside each gridbox. Hence the map
will account for the "fly time" of each trajectory in each
location. This map can be normalized by:

* The total # of trajectories (note: the max value on the map may be >> 100%)
* The total # of points (note: the max value will be << 100%)
* The maximum counts on any gridbox on the map (note: the max value will be exactly 100%). 

The normalization will change the values, but the map should be the same. 

## Trying the different normalizations

As I suspected, the different normalizations don't
really matter. They change the scale of the plot (i.e. the numbers
on the color bar), but the map is the same.  Hence, from now on I will
normalize by the max count in any gridbox.

## Single count x Time of flight

Then I compared counting each trajectory just once x counting the flying time.
Not a big difference either.

## Different scales

Then I tried different log- and linear scales.
Linear does not work, as we can only see the trajectories very close to ATTO.
Forcing a log scale from 1 to 100 better.


