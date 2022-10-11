#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/5A_bbmap.log
#SBATCH -o logs/5A_bbmap.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel
#SBATCH -J Ant_map_MRA 

DIR=clustered_vOTUs
CPU=16
SAMPFILE=data/fastq/MRA/MRA_samples.csv
READDIR=data/fastq/MRA
OUT=data/fastq/bbmap

#mkdir $OUT

/rhome/cassande/bigdata/software/bbmap/bbmap.sh ref=$DIR/combined_all_clustered.fa

module load samtools
module load bedtools

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read FILE PREFIX
do
	#mkdir $OUT/$PREFIX

	if [ ! -f $OUT/$PREFIX'_mapped_coverge.tsv' ]; then
		#paired
		#/rhome/cassande/bigdata/software/bbmap/bbmap.sh in1=$READDIR/$PREFIX'_1_filtered.fastq.gz' in2=$READDIR/$PREFIX'_2_filtered.fastq.gz' out=$OUT/$PREFIX.sam minid=0.90 threads=$CPU

		#interleaved
		/rhome/cassande/bigdata/software/bbmap/bbmap.sh in=$READDIR/$FILE  out=$OUT/$PREFIX.sam minid=0.90 threads=$CPU

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
