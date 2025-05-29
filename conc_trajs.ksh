#!/bin/ksh

# https://hysplitbbs.arl.noaa.gov/viewtopic.php?t=1311

echo <<EOF > LABELS.CFG
'TITLE&','Trajectory Frequency&'
'MAPID&',' Values &'
'UNITS&',' %&'
'VOLUM&',' &'
EOF

for file in clean_day polluted_day ; do
    
    /bin/ls ./$file/* > ${file}_list

    for ropt in 0 1 2 3 ; do
	out=freq_${file}_R${ropt}.bin
	echo $out
	
	../run/exec/trajfreq -f${out}.bin -g1.0 -iclean_day_list -r$ropt -s0:99999

	../run/exec/concplot -i${out}.bin -o${out}.ps -j../run/graphics/arlmap -m0 -k1
    done
done

#You can manually re-run concplot and force the desired contours using the "-v" option as follows in a Windows Command Prompt window, after cd /hysplit4/working:
#..\exec\concplot -itfreq.bin -ofreqplot.ps -jC:/hysplit4/graphics/arlmap -m0 -k1 -v10+20+30+40+50+60+70+80+90+100

#To see all the possible command line inputs to concplot, run "..\exec\concplot" in the Command Prompt window.

#To get the appropriate trajectory frequency plot labeling, you must create the file LABELS.CFG before running concplot.

#
