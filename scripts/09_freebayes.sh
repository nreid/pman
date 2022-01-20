#!/bin/bash 
#SBATCH --job-name=freebayes_parallel
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 30
#SBATCH --mem=90G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

##################################
# run freebayes to call variants
##################################

# this script calls variants using freebayes

#load software---------------------------------------------------------------
module load htslib/1.10.2
module load freebayes/1.3.1
module load vcflib/1.0.0-rc1
module load parallel/20180122

#input, output files, directories--------------------------------------------
INDIR=../results/aligned_ddRAD/
OUTDIR=../results/freebayes
mkdir -p $OUTDIR

#reference genome
GEN=../genome/GCA_003704035.3_HU_Pman_2.1.3_genomic.fna

# make a list of bam files
	# use popmap_total.txt to construct the list
POPMAP=../meta/popmap.txt

BAMLIST=$OUTDIR/bam.list
cat $POPMAP | cut -f 1 | sed 's/$/.bam/' | sed  "s,^,$INDIR," >$BAMLIST

# these are freebayes scripts found in the same location as the executable
MAKEREGIONS=/isg/shared/apps/freebayes/1.3.1/scripts/fasta_generate_regions.py

# run freebayes-parallel--------------------------------------------------------

bash freebayes_parallel.sh \
	<(python $MAKEREGIONS ${GEN}.fai 5000000) 30 \
	-f ${GEN} \
	--bam-list $BAMLIST \
	-m 30 \
	-q 20 \
	--skip-coverage 5000 | \
bgzip -c >$OUTDIR/freebayes.vcf.gz

# index the vcf file
tabix -p vcf $OUTDIR/freebayes.vcf.gz