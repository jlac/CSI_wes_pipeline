#! /bin/bash

###################
#
# Launching shell script for NIAID CSI batch processing of WES data
#
###################
module load snakemake/5.4.0-Python-3.6.7

##
## Location of snakemake
##
DIR="$( cd -P "$( dirname "$0" )" >/dev/null 2>&1 && pwd )"
echo "Running script from ${DIR}"

##
## Test commandline arguments
##
if [ $# -ne 1 ]; then
    echo " " 
    echo "Requires a single commandline argument: gris, npr, or process"
    echo " " 
    exit
fi

if [ $1 != "gris" ] && [ $1 != "npr" ] && [ $1 != "process" ] && [ $1 != "himem" ] ; then
    echo " " 
    echo "Invalid commandline option: $1"
    echo "Valid commandline options include: gris, npr, or process"
    echo " " 
    exit
fi

##
## Get batch and batch number
##
batchdir=`pwd`
batch=`echo $batchdir | sed -e 's/^.*\///' `
echo "BATCH: $batch"
batchnumber=`echo $batch | sed -e 's/BATCH//' -e 's/^0//' `
echo "Processing Batch $batchnumber"

##
## Find the raw directory
##
raw_root="/hpcdata/dir/CIDR_DATA_RAW"
raw_dir="${raw_root}/${batch}/Released_Data_Batch${batchnumber}/Holland_WES_Release_${batchnumber}"

if ! test -d "rawdata"; then
    if test -d $raw_dir; then
        echo "Linking rawdata subdirectory to $raw_dir"
        ln -s $raw_dir rawdata
    else
        echo "Unable to locate raw data directory $raw_dir"
        echo "Exiting"
        exit
    fi
else
    echo "input directory rawdata already exists"
fi

##
## Make the new output directories
##
for i in BAM VCF QC/TARGET QC/UCSC BATCH_QC HLA inbreeding CNV
do
    if ! test -d $i; then
        echo "Creating output directory: $i"
        mkdir -p $i
    else
        echo "output directory $i already exists"
    fi
done
# and be sure this directory and all subdirectories are writable
chmod -fR g+rwx .


##
## Run csi_to_gris.py
##
if [ "$1" == "gris" ] 
then
    echo "Running csi_to_gris"
    csi_to_gris.py -b $batchnumber
    exit
fi

mkdir -p snakejobs
mkdir -p BATCH_QC

##
## Run snakemake
##
echo "Run snakemake"

#CLUSTER_OPTS="qsub -e snakejobs -o snakejobs -wd $batchdir"
CLUSTER_OPTS="qsub -e snakejobs -o snakejobs -pe threaded {cluster.threads} -l {cluster.partition} -l h_vmem={cluster.vmem} -l mem_free={cluster.mem} -wd $batchdir"
#CLUSTER_OPTS="qsub -e snakejobs -o snakejobs -pe threaded {cluster.threads} -l {cluster.partition} -cwd"
#CLUSTER_OPTS_HIMEM="qsub -e snakejobs_himem -o snakejobs_himem -pe threaded 8 -l himem -l h_vmem=32G -l virtual_free=32G -wd $batchdir"

if [ "$1" == "npr" ]
then
    snakemake -npr --rerun-incomplete --snakefile CSI_wes_pipeline/scripts/csi_batch_processing_full.snakemake
fi

if [ "$1" == "process" ]
then
    snakemake --stats snakemake.stats --rerun-incomplete -j 150 --cluster "$CLUSTER_OPTS" --cluster-config CSI_wes_pipeline/resources/processing_cluster.json --keep-going --snakefile CSI_wes_pipeline/scripts/csi_batch_processing_full.snakemake 2>&1|tee -a csi_batch_processing.log
fi