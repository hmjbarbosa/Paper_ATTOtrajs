#!/bin/ksh

trajdir=backward

base=`pwd`

for file in clean_day polluted_day ; do 

    cd ${base}
    rm -Rf ${file}
    mkdir ${file}
    cd ${file}

    while IFS=- read -r y m d
    do
	for sites in 1 2 3 4 5 6 7 ; do
	    #echo "${d:0:2}|"
	    ln -s ${base}/${trajdir}/00${sites}/t${y:2:2}${m}${d:0:2}_*_${trajdir}.dat.* ./
	done
    done < "${base}/${file}.txt"

done
#
