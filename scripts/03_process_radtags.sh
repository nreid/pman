#!/bin/bash
#SBATCH --job-name=process_radtags
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-3]

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

module load stacks/2.53

#input/output directories, supplementary files
INDIR=../rawdata

FQ1ARRAY=($(ls -1 $INDIR/Davis*R1*fastq.gz))
FQ2ARRAY=($(ls -1 $INDIR/Davis*R2*fastq.gz))
BARCODES=($(ls -1 ../meta/*barcode*))

# make demultiplexed directory if it doesn't exist
OUTDIR=../results/demultiplexed_fastqs/pool_$SLURM_ARRAY_TASK_ID
mkdir -p $OUTDIR

FASTQ1=${FQ1ARRAY[$SLURM_ARRAY_TASK_ID]}
FASTQ2=${FQ2ARRAY[$SLURM_ARRAY_TASK_ID]}
BC=$(echo ${BARCODES[$SLURM_ARRAY_TASK_ID]})

echo demultiplexing file pair $FASTQ1 using barcode set $BC

process_radtags \
-1 $FASTQ1 \
-2 $FASTQ2 \
-b $BC \
-o $OUTDIR \
-i gzfastq \
-y gzfastq \
-e nsiI \
-c \
-q \
-s 20 \