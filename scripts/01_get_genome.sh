#!/bin/bash
#SBATCH --job-name=get_genome
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

######################
# get reference genome
######################

# download genome assembly for peromyscus maniculatus bairdii
# NCBI accession # GCA_003704035.3

# output directory
GENOMEDIR=../genome
mkdir -p $GENOMEDIR

# download genome
wget \
-P $GENOMEDIR \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/704/035/GCA_003704035.3_HU_Pman_2.1.3/GCA_003704035.3_HU_Pman_2.1.3_genomic.fna.gz
# decompress the genome
gunzip $GENOMEDIR/GCA_003704035.3_HU_Pman_2.1.3_genomic.fna.gz


# index the genome using bwa
module load bwa/0.7.17
bwa index \
-p $GENOMEDIR/pman \
$GENOMEDIR/GCA_003704035.3_HU_Pman_2.1.3_genomic.fna

# index the genome using samtools
module load samtools/1.10
samtools faidx $GENOMEDIR/GCA_003704035.3_HU_Pman_2.1.3_genomic.fna


