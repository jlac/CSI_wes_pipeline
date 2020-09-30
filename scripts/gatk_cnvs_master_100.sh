#!/bin/sh

# Job Name
#$ -N gatkCNV

# Execute the script from the Current Working Directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file

# Send the output of the script to a directory called 'UGE-output' in the current working directory (cwd)
# **NOTE: Be sure to create this directory before submitting your job, UGE scheduler will NOT create this**
# **directory for you**

# Tell the job your memory requirements
#$ -l h_vmem=24G

# threads
#$ -pe threaded 2

# Send mail when the job is submitted, and when the job completes
#$ -m be

#  Specify an email address to use
#$ -M justin.lack-@nih.gov

SAMPLE=$1

source /hpcdata/dir/software/resources/cnv_normals/gatk_cnv/conda/etc/profile.d/conda.sh
conda activate base
conda activate gatk
module load java/1.8.0_92
mkdir -p CNV_100
mkdir -p CNV_100
mkdir -p CNV_100/counts
mkdir -p CNV_100/ploidy
mkdir -p CNV_100/model
/hpcdata/dir/software/gatk-4.1.1.0/gatk CollectReadCounts -I BAM/$SAMPLE.bam -L /hpcdata/dir/CIDR_DATA_RENAMED/references/exome_targets.bed --interval-merging-rule OVERLAPPING_ONLY -O CNV_100/counts/$SAMPLE.counts.hdf5
/hpcdata/dir/software/gatk-4.1.1.0/gatk DetermineGermlineContigPloidy --input CNV_100/counts/$SAMPLE.counts.hdf5 --model /hpcdata/dir/CIDR_DATA_RENAMED/CNVBAMS/Ploidy/CohortMode/CIDR_COHORT-model --output-prefix ploidy$SAMPLE --output CNV_100/ploidy/ploidy$SAMPLE
mkdir -p CNV_100/gatk/model/cnv$SAMPLE
/hpcdata/dir/software/gatk-4.1.1.0/gatk GermlineCNVCaller --run-mode CASE --contig-ploidy-calls CNV_100/ploidy/ploidy$SAMPLE/ploidy$SAMPLE-calls --model /hpcdata/dir/software/resources/new_cnv_pon/GermlineCNVmodel_100/CIDR_CNVcalls_COHORT_100-model --input CNV_100/counts/$SAMPLE.counts.hdf5 --output CNV_100/model/cnv$SAMPLE --output-prefix cnv$SAMPLE
/hpcdata/dir/software/gatk-4.1.1.0/gatk PostprocessGermlineCNVCalls --calls-shard-path CNV_100/model/cnv$SAMPLE/cnv$SAMPLE-calls --model-shard-path /hpcdata/dir/software/resources/new_cnv_pon/GermlineCNVmodel_100/CIDR_CNVcalls_COHORT_100-model --contig-ploidy-calls CNV_100/ploidy/ploidy$SAMPLE/ploidy$SAMPLE-calls --allosomal-contig X --allosomal-contig Y --output-genotyped-intervals CNV_100/genotyped_intervals_$SAMPLE.vcf --output-genotyped-segments CNV_100/genotyped_segments_$SAMPLE.vcf