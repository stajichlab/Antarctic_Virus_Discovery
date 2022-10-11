#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -e logs/02_dramv.%a.log
#SBATCH -o logs/02_dramv.%a.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p batch
#SBATCH -J Ant_virus_dramv
#SBATCH -a 1-185

CPU=24
DIR=data/Assemblies
SAMPFILE=data/samples.csv
OUT=viral_analysis


#export CHECKVDB=/rhome/cassande/bigdata/software/checkv-db-v1.0

CHECKVDB=/rhome/cassande/bigdata/software/checkv-db-v1.0

mkdir $OUT

module unload miniconda2
module load miniconda3

#conda activate checkv
#only need to do this once
#DRAM-setup.py prepare_databases --skip_uniref --output_dir /rhome/cassande/bigdata/software/db-dramv
#conda deactivate

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
 N=$1
 if [ -z $N ]; then
     echo "cannot run without a number provided either cmdline or --array in sbatch"
     exit
 fi
fi

conda activate checkv

tail -n +2 $SAMPFILE | sed -n ${N}p | while read SAMPLE
do 
	#dramv can take a long time - so we will skip in first pass
	
	# step 1 annotate
	DRAM-v.py annotate -i $OUT/$SAMPLE'_vs2_pass2'/for-dramv/final-viral-combined-for-dramv.fa -v $OUT/$SAMPLE'_vs2_pass2'/for-dramv/viral-affi-contigs-for-dramv.tab -o $OUT/$SAMPLE'_dramv-annotate' --skip_trnascan --threads $CPU --min_contig_size 1000

	#step 2 summarize anntotations
	DRAM-v.py distill -i $OUT/SAMPLE'_dramv-annotate'/annotations.tsv -o $OUT/SAMPLE'_dramv-distill'
	

done
