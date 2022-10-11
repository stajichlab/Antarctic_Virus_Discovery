#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/vcontact2.ICTV.ant.B.log
#SBATCH -o logs/vcontact2.ICTV.ant.B.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p batch
#SBATCH -J ant_ictv_vcontact2


DIR=ant
OUT=Ant_vcontact2_ictv
PREFIX=Ant_plus_1Aug2022_ICTV
OUTFILE=Ant_plus_1Aug2022_ICTV_prot.fa


#wget http://inphared.s3.climb.ac.uk/1Aug2022_data.tsv
#wget http://inphared.s3.climb.ac.uk/1Aug2022_vConTACT2_family_annotations.tsv
#wget http://inphared.s3.climb.ac.uk/1Aug2022_vConTACT2_gene_to_genome.csv
#wget http://inphared.s3.climb.ac.uk/1Aug2022_vConTACT2_genus_annotations.tsv
#wget http://inphared.s3.climb.ac.uk/1Aug2022_vConTACT2_host_annotations.tsv
#wget http://inphared.s3.climb.ac.uk/1Aug2022_vConTACT2_lowest_taxa_annotations.tsv
#wget http://inphared.s3.climb.ac.uk/1Aug2022_vConTACT2_proteins.faa
#wget http://inphared.s3.climb.ac.uk/1Aug2022_vConTACT2_subfamily_annotations.tsv


module load centos
centos.sh

#vcontact2 on vOTUs

conda activate vContact2

#prodigal -i $DIR/$PREFIX.fa -o $DIR/$PREFIX.coords.gbk -a $DIR/$OUTFILE -p meta

#vcontact2_gene2genome -p $DIR/$OUTFILE -o $DIR/$PREFIX.gene2genome.csv -s Prodigal-FAA

#combined with ICTV files
#https://github.com/RyanCook94/inphared#supplementing-and-annotating-vcontact2-clusters
#combine the date_vConTACT2_proteins.faa with your own fasta of file of translated ORFs, and combine date_vConTACT2_gene_to_genome.csv with your own mapping file (watch out for duplicated headers in the gene_to_genome.csv file if your file already has headers). Then run vConTACT2 as normal using the --db 'None' option, as this will avoid RefSeq duplicates.


vcontact2 --raw-proteins $DIR/$OUTFILE --rel-mode Diamond --proteins-fp $DIR/$PREFIX.gene2genome.csv --db 'None' --pcs-mode MCL --vcs-mode ClusterONE  --output-dir $DIR/$OUT