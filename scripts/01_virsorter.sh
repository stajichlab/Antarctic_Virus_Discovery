#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -e logs/01_virsorter_checkv.%a.log
#SBATCH -o logs/01_virsorter_checkv.%a.log
#SBATCH --mem 150G #memory per node in Gb
#SBATCH -p intel,batch
#SBATCH -J Ant_virus_virsorter
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

tail -n +2 $SAMPFILE | sed -n ${N}p | while read SAMPLE
do 
	conda activate virsorter2
	#using a loose cutoff of 0.5 for maximal sensitivity - bc will trim with checkv
	virsorter run --keep-original-seq -i $DIR/$SAMPLE.fa -w $OUT/$SAMPLE'_vs2_pass1' --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --min-length 5000 --min-score 0.5 -j $CPU all
	conda deactivate
	
	conda activate checkv
	#checkv to trim
	checkv end_to_end $OUT/$SAMPLE'_vs2_pass1'/final-viral-combined.fa $OUT/$SAMPLE.checkv -t $CPU -d $CHECKVDB
	conda deactivate
	
	#re-run virsorter to get data ready for dramv
	cat $OUT/$SAMPLE.checkv/proviruses.fna $OUT/$SAMPLE.checkv/viruses.fna > $OUT/$SAMPLE.checkv/combined.fna
	
	conda activate virsorter2
	virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -i $OUT/$SAMPLE.checkv/combined.fna -w $OUT/$SAMPLE'_vs2_pass2' --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --min-length 5000 --min-score 0.5 -j $CPU all
	conda deactivate
	
	#dramv can take a long time - so we will skip in first pass
	#conda activate checkv
	
	# step 1 annotate
	#DRAM-v.py annotate -i $OUT/$SAMPLE'_vs2_pass2'/for-dramv/final-viral-combined-for-dramv.fa -v $OUT/$SAMPLE'_vs2_pass2'/for-dramv/viral-affi-contigs-for-dramv.tab -o $OUT/$SAMPLE'_dramv-annotate' --skip_trnascan --threads $CPU --min_contig_size 1000

	#step 2 summarize anntotations
	#DRAM-v.py distill -i $OUT/SAMPLE'_dramv-annotate'/annotations.tsv -o $OUT/SAMPLE'_dramv-distill'
	
	#conda deactivate

done
