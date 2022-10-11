#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/05D_bbmap.%a.log
#SBATCH -o logs/05D_bbmap.%a.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p batch
#SBATCH -J Ant_map_missing
#SBATCH -a 1-9


DIR=clustered_vOTUs
CPU=16
SAMPFILE=data/fastq/missing/missing.samples.txt
READDIR=data/fastq/missing
OUT=data/fastq/bbmap

#mkdir $OUT

#/rhome/cassande/bigdata/software/bbmap/bbmap.sh ref=$DIR/combined_all_clustered.fa

module load samtools
module load bedtools

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
 N=$1
 if [ -z $N ]; then
     echo "cannot run without a number provided either cmdline or --array in sbatch"
     exit
 fi
fi

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read PREFIX
do
	#mkdir $OUT/$PREFIX

	if [ ! -f $OUT/$PREFIX'_mapped_coverge.tsv' ]; then
		#paired
		#/rhome/cassande/bigdata/software/bbmap/bbmap.sh in1=$READDIR/$PREFIX'_R1.fastq.gz' in2=$READDIR/$PREFIX'_R2.fastq.gz' out=$OUT/$PREFIX.sam minid=0.90 threads=$CPU

		#interleaved
		/rhome/cassande/bigdata/software/bbmap/bbmap.sh in=$READDIR/$PREFIX.fastq.gz  out=$OUT/$PREFIX.sam minid=0.90 threads=$CPU

		# make bam files from sam files
		samtools view -F 4 -bS $OUT/$PREFIX.sam | samtools sort > $OUT/$PREFIX'_sortedIndexed.bam'

		# index the bam file
		samtools index $OUT/$PREFIX'_sortedIndexed.bam'

		#get coverage
		bedtools genomecov -ibam $OUT/$PREFIX'_sortedIndexed.bam' -bga -max 10 > $OUT/$PREFIX'_mapped_coverge.tsv'

		#gzip $OUT/$PREFIX'_sortedIndexed.bam'
		#gzip $OUT/$PREFIX'_sortedIndexed.bam.bai'
		rm $OUT/$PREFIX.sam

	else
		echo "$PREFIX has already been mapped"
	fi


done
