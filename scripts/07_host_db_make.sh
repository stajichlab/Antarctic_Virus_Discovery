#!/usr/bin/bash
#SBATCH -n 8 --mem 24gb -p intel --out ncbi.%A.log
#SBATCH -J cyano_db


CPU=8
export PATH="/rhome/cassande/.local/bin:$PATH"

#get assesions 
ncbi-genome-download --dry-run --assembly-levels complete bacteria > bacteria.txt

#download
ncbi-genome-download --parallel $CPU --formats fasta  --assembly-levels complete bacteria


#get assesions 
ncbi-genome-download --dry-run --assembly-levels complete archaea > archaea.txt

#download
ncbi-genome-download --parallel $CPU --formats fasta --assembly-levels complete archaea