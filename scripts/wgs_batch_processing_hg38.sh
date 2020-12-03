#! /bin/bash

###################
#
# Launching shell script for NIAID CSI batch processing of WES data
#
###################
module load snakemake/5.8.2-Python-3.6.7

##
## Location of snakemake
##
DIR="$( cd -P "$( dirname "$0" )" >/dev/null 2>&1 && pwd )"
echo "Running script from ${DIR}"

##
## Test commandline arguments
##
if [ $# -ne 2 ]; then
    echo " " 
    echo "Requires a single commandline argument: gris, npr, or process"
    echo " " 
    exit
fi

if [ $2 != "npr" ] && [ $2 != "process" ] ; then
    echo " " 
    echo "Invalid commandline option: $2"
    echo "Valid commandline options include: gris, npr, or process"
    echo " " 
    exit
fi

##
## Get batch and batch number
##
batchdir=`pwd`
#batch=`echo $batchdir | sed -e 's/^.*\///' `
#echo "BATCH: $batch"
#batchnumber=`echo $batch | sed -e 's/BATCH//' -e 's/^0//' `
#echo "Processing Batch $batchnumber"

##
## Find the raw directory
##
#raw_root="/data/NCBR/rawdata/csi_test_batch"
raw_dir=$1

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
for i in BAM VCF BATCH_QC HLA inbreeding snakejobs
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

SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`

echo $SCRIPT
echo $SCRIPTPATH

if ! test -e "batches"; then
  perl CSI_wes_pipeline/scripts/split_gvcf_batches.pl masterkey.txt 50 hg38
fi

##
## Run snakemake
##
echo "Run snakemake"

CLUSTER_OPTS="qsub -e snakejobs -o snakejobs -pe threaded {cluster.threads} -l {cluster.partition} -l h_vmem={cluster.vmem} -l mem_free={cluster.mem} -wd $batchdir"

if [ "$2" == "npr" ]
then
    snakemake -npr --snakefile CSI_wes_pipeline/scripts/wgs_batch_processing_hg38.snakemake
fi

if [ "$2" == "process" ]
then
    snakemake --stats snakemake.stats --restart-times 1 --rerun-incomplete -j 150 --cluster "$CLUSTER_OPTS" --cluster-config CSI_wes_pipeline/resources/processing_cluster_locus.json --keep-going --snakefile CSI_wes_pipeline/scripts/wgs_batch_processing_hg38.snakemake 2>&1|tee -a csi_batch_processing.log
fi