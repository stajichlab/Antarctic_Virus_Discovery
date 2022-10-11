#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -e logs/03_cdhit.log
#SBATCH -o logs/03_cdhit.log
#SBATCH --mem 100G #memory per node in Gb
#SBATCH -p intel
#SBATCH -J Ant_virus_cdhit

CPU=24
IN=clustered_vOTUs/data
OUT=clustered_vOTUs

INFILE=Antartica_viruses_checkv_combined_all_renamed.fa

#module load cd-hit/4.8.1  

/rhome/cassande/bigdata/software/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $IN/$INFILE -o $OUT/combined_all_clustered.fa -c 0.95 -aS 0.85 -d 0 -M 0 -T $CPU

