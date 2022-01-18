#!/bin/bash
#SBATCH --job-name=clone_filter
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-3]

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

module load stacks/2.53

#input/output directories, supplementary files
INDIR=../rawdata
OUTDIR=../results/dedup
mkdir -p $OUTDIR

POOLS=($(ls -1 $INDIR/Davis*R1*fastq.gz))

FASTQ1=$(echo ${POOLS[$SLURM_ARRAY_TASK_ID]})
FASTQ2=$(echo $FASTQ1 | sed 's/R1/R2/')

# run clone_filter

clone_filter \
-i gzfastq \
-o $OUTDIR \
-1 $FASTQ2 \
-2 $FASTQ1 \
--inline_null \
--oligo_len_1 8