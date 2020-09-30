#!/bin/sh

# Job Name
#$ -N inbreeding

# Execute the script from the Current Working Directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file

# Send the output of the script to a directory called 'UGE-output' in the current working directory (cwd)
# **NOTE: Be sure to create this directory before submitting your job, UGE scheduler will NOT create this**
# **directory for you**

# Tell the job your memory requirements
#$ -l h_vmem=8G

# threads
#$ -pe threaded 8

BATCH=$1

module load GATK/3.7-0-Java-1.8.0_92
mkdir -p inbreeding
java -Djava.io.tmpdir=inbreeding/tmp -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants -R /hpcdata/dir/CIDR_DATA_RENAMED/references/human_g1k_v37_decoy.fasta -L /hpcdata/dir/CIDR_DATA_RENAMED/references/exome_targets.bed -V VCF/$BATCH.vcf.gz -o inbreeding/$BATCH.knownSites.onTarget.vcf --concordance /hpcdata/dir/software/resources/1k_genomes_phase3_autosomes.vcf.gz -nt 8
module purge
module load uge/8.5.5
module load vcftools/0.1.15-goolf-1.7.20-Perl-5.22.2
vcftools --vcf inbreeding/$BATCH.knownSites.onTarget.vcf --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 1 --recode --recode-INFO-all --out inbreeding/$BATCH.knownSites.onTarget.qualFiltered
module load tabix
bgzip inbreeding/$BATCH.knownSites.onTarget.qualFiltered.recode.vcf
tabix -p vcf inbreeding/$BATCH.knownSites.onTarget.qualFiltered.recode.vcf.gz
rm inbreeding/$BATCH.knownSites.onTarget.vcf
vcftools --gzvcf inbreeding/$BATCH.knownSites.onTarget.qualFiltered.recode.vcf.gz --plink --out inbreeding/$BATCH
module purge
module load uge/8.5.5
module load picard/2.17.6
java -jar ${EBROOTPICARD}/picard.jar CollectVariantCallingMetrics INPUT=inbreeding/$BATCH.knownSites.onTarget.qualFiltered.recode.vcf.gz OUTPUT=inbreeding/picardMetrics DBSNP=/hpcdata/dir/software/resources/00-All.vcf.gz THREAD_COUNT=8
/hpcdata/dir/software/plink-1.07-x86_64/plink --file inbreeding/$BATCH --homozyg --out inbreeding/ROH
module purge
module load uge/8.5.5
module load R
Rscript /hpcdata/dir/software/inbreedingPlot.R