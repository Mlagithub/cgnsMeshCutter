#!/bin/bash

getMem(){
    ps aux | grep cutter | grep 1156Bar | grep -v mpirun | awk '{print $6/1024/1024 "\t"}' | xargs | awk '{print $1}'
}

file="mem.log"

echo "#Residual" > $file 
echo "# Step Mem(Gb)" >> $file

n=0
while true
do
    n=$((n+1))
    echo -n "$n " >> $file
    getMem >> mem.log;
    sleep 1;
done 