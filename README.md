# CSI_wes_pipeline
Whole exome processing workflow for CSI data generated by CIDR and HGSC

To run, submit the master shell script as a qsub on LOCUS as follows:

qsub -wd /path/to/existing/working/dir HGSC_wes_pipeline/scripts/hgsc_batch_processing_hg19_collaborator.sh /path/to/raw/data process

to perform a dry run:

./HGSC_wes_pipeline/scripts/hgsc_batch_processing_hg19_collaborator.sh /path/to/raw/data npr