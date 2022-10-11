#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/vcontact2.ICTV.ant.A.log
#SBATCH -o logs/vcontact2.ICTV.ant.A.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p batch
#SBATCH -J ant_ictv_vcontact2


DIR=clustered_vOTUs
PREFIX=combined_all_clustered
OUTFILE=combined_all_clustered_prot.fa

#vcontact2 on vOTUs

conda activate vContact2

prodigal -i $DIR/$PREFIX.fa -o $DIR/$PREFIX.coords.gbk -a $DIR/$OUTFILE -p meta

vcontact2_gene2genome -p $DIR/$OUTFILE -o $DIR/$PREFIX.gene2genome.csv -s Prodigal-FAA
