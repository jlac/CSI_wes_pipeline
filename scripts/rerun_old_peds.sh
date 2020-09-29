#!/bin/bash
set -e

module load python

if [ $# -ne 1 ]; then
    echo " " 
    echo "Requires a single integer argument: batch number"
    echo " " 
    exit
else
    echo "Running batches 1 to $1"
fi

# The first set of batches were not consistent, and don't use the rawdata subdirectory
#echo "Running batches 1-4"
for ((i=1;i<=4;i++)); 
do
    echo "Running Batch $i"
    yes | csi_to_gris.py -b $i -p -d /data/NCBR/projects/csi_test_batch/old_peds/BATCH${i}
done

#echo "Running batches 5-${1}"
for ((i=5;i<=$1;i++)); 
do
    if [[ $i == 8 ]]; then 
        echo "Skipping Batch $i"
        continue
    else
        echo "Running Batch $i"
        yes | csi_to_gris.py -b $i -p -d /data/NCBR/projects/csi_test_batch/old_peds/BATCH${i}
    fi
done
