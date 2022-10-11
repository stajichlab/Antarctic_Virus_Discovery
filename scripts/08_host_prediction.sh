#!/usr/bin/bash
#SBATCH -n 8 --mem 24gb -p intel --out ncbi_makedb.combo.%A.log
#SBATCH -J cyano_ncbi

source activate blast

DIR=refseq

gunzip -c $DIR/archaea_and_bacteria_complete.fna.gz | makeblastdb -in - -dbtype nucl -out bactarch -title bactarch

QUERY=data/AntVirus_combined_all_clustered.fa
CPU=8

OUT=blast_results
#mkdir $OUT

DB=bactarch
blastn -query $QUERY -db $DB -outfmt "6 qseqid sseqid pident bitscore evalue length sgi sacc sallseqid staxids sscinames stitle" -num_threads $CPU -evalue 1e-25 -out $OUT/$DB.blastn.tab


