#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/6_bamm.log
#SBATCH -o logs/6_bamm.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p batch
#SBATCH -J Ant_bamm_parse



OUT=data/fastq/bbmap
DIR=vOTU_RA


module unload miniconda2
module load miniconda3

conda activate bamm

bamm parse -c $DIR/output_file_tpmean_072622.tsv -b $OUT/*.bam -m 'tpmean'
bamm parse -c $DIR/output_file_count_072622.tsv -b $OUT/*.bam -m 'counts'

#gzip $OUT/*bam
#gzip $OUT/*bai
