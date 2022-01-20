#!/bin/bash
#SBATCH --job-name=filter_vcfs
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

##################################
# filter variant sets
##################################

module load bcftools/1.9
module load htslib/1.9
module load vcftools/0.1.16
module load vcflib/1.0.0-rc1


###############################
# set input, output directories
###############################

OUTDIR=../results/stacks
mkdir -p $OUTDIR

#############################
# filter SITES by missingness
#############################

# also remove multiallelic sites and indels

vcftools --gzvcf $OUTDIR/populations.snps.dict.vcf.gz \
	--max-missing-count 10 --mac 3 --remove-indels --max-alleles 2 --min-alleles 2 \
	--recode \
	--stdout | \
	bgzip >$OUTDIR/filtered.vcf.gz

	# output missing individual report
	vcftools --gzvcf $OUTDIR/filtered.vcf.gz --out $OUTDIR/filtered --missing-indv

##############
# make indexes
##############

for file in $OUTDIR/*vcf.gz
do tabix -f -p vcf $file
done
