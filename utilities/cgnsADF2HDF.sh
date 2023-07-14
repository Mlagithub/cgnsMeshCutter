#!/bin/bash

if [ $2 -gt 100000 ]; then
    len=6
elif [ $2 -gt 10000 ]; then 
    len=5
elif [ $2 -gt 1000 ]; then 
    len=4
elif [ $2 -gt 100 ]; then
    len=3
elif [ $2 -gt 10 ]; then 
    len=2
else
    len=1
fi 

printf -v fmtStr "%%s_%%.%dd.cgns" $len

for i in `seq 0 $(($2-1))`
do 
    printf -v name $fmtStr $1 $i
    # name="$1_$i.cgns"
    cgnsconvert $name
    cgnscheck $name
done