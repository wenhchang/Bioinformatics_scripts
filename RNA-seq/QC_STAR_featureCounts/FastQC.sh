#!/bin/bash

#!/bin/sh
#SBATCH --job-name=FastQC
#SBATCH --output=FastQC.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6400
#SBATCH --time=24:00:00

#   This SLURM batch script runs FastQC on all FASTQ files in the current directory.
#   It generates quality control reports for each FASTQ file.
# Requirements:
#   - SLURM job scheduler
#   - FastQC module installed and loaded    
#   - Input FASTQ files should be in the format *.fq.gz
#   - Output reports will be saved in the 'fast_qc_output' directory

# Create output directory       
mkdir -p fast_qc_output

# Load the FastQC module

module load fastqc

fastqc -o fast_qc_output */*.fq.gz
