#!/bin/bash
set -e

module load python/3.7.3-foss-2016b

if [ $# -ne 1 ]; then
    echo " " 
    echo "Requires a single integer argument: batch number"
    echo " " 
    exit
else
    echo "Running batches 1 to $1"
fi

for ((i=1;i<=$1;i++)); 
do
    if [[ $i == 8 ]]; then 
        echo "Skipping Batch $i"
        continue
    else
        echo "Running Batch $i"
        python /hpcdata/dir/SCRIPTS/generate_seqr_ped.py -b $i
    fi
done
