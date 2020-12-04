#! /bin/bash

###################
#
# Launching shell script for NIAID CSI batch processing of WES data
#
###################
module load snakemake/5.7.4
module load python/3.7

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
    echo "Requires a single commandline argument: npr or process"
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

##
## Run snakemake
##
echo "Run snakemake"

CLUSTER_OPTS="sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e snakejobs/slurm-%j_{params.rname}.out -o snakejobs/slurm-%j_{params.rname}.out --chdir=$batchdir"

if [ "$2" == "npr" ]
then
    snakemake -npr --snakefile CSI_wes_pipeline/scripts/wgs_batch_processing_hg38_biowulf.snakemake
fi

if [ "$2" == "process" ]
then
    snakemake --stats snakemake.stats --rerun-incomplete -j 150 --cluster "$CLUSTER_OPTS" --cluster-config CSI_wes_pipeline/resources/biowulf_processing_cluster_hg38.json --keep-going --snakefile CSI_wes_pipeline/scripts/wgs_batch_processing_hg38_biowulf.snakemake 2>&1|tee -a csi_batch_processing.log
fi