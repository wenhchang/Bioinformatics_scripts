#!/bin/bash
#SBATCH --job-name=star_mapping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=star_mapping_A_%j.out
#SBATCH --error=star_mapping_A_%j.err


#   This SLURM batch script performs RNA-seq data processing for multiple samples.
#   For each sample folder (e.g., A001, A002...), it executes:
#     1. Adapter trimming using Cutadapt
#     2. Alignment to the reference genome using STAR (with two-pass mode)
#     3. Gene-level quantification using featureCounts
#
# Requirements:
#   - SLURM job scheduler
#   - STAR index pre-built in the specified genomeDir path
#   - Input paired-end fastq.gz files named as <sample>_1.fq.gz and <sample>_2.fq.gz
#   - GTF annotation file available at specified path
#
# Outputs:
#   - Trimmed fastq files
#   - STAR-aligned sorted BAM files
#   - Gene counts tables from featureCounts
#

# Load modules
module load star
module load subread 
module load cutadapt

# Define output directory
output_dir="STAR_result"
trim_dir="trimmed_fastq"
mkdir -p $output_dir

# Path to the GTF file
gtf_file="Ref/gv44_basicannotation.gtf"

# generate sample list (adjust the pattern as needed)
sample_list=($(ls -d A*/ 2>/dev/null | sed 's#/##'))

# Process each sample one at a time
for sample in "${sample_list[@]}"; do
    # Define input files for STAR
    raw_fq1="${sample}/${sample}_1.fq.gz"
    raw_fq2="${sample}/${sample}_2.fq.gz"

    # Define trimmed files
    trim_fq1="${trim_dir}/${sample}_1_trimmed.fastq"
    trim_fq2="${trim_dir}/${sample}_2_trimmed.fastq"

    # Run Cutadapt
    cutadapt \
    -j 32\
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A AATGATACGGCGACCACCGAGATCTACAC \
    -o $trim_fq1 \
    -p $trim_fq2 \
    -m 20\
    $raw_fq1 $raw_fq2

    #j defines how many cpus are used; m defines cutting out reads with lengths smaller than 20


    # Run STAR mapping for each sample
    STAR --runThreadN 30 \
         --genomeDir /work/users/w/h/whchang/RNA_seq/Reference/index \
         --readFilesIn ${trim_fq1} ${trim_fq2} \
         --outFileNamePrefix ${output_dir}/mapping_${sample}_ \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         --twopassMode Basic

    echo "STAR mapping completed for ${sample}."

    # Define the STAR output BAM file
    star_output_bam="${output_dir}/mapping_${sample}_Aligned.sortedByCoord.out.bam"

    # Run featureCounts for the STAR output BAM file
    featureCounts -T 24 -a $gtf_file -o ${output_dir}/${sample}_counts.txt -p -B -C -g gene_id -s 2 $star_output_bam

    echo "featureCounts completed for ${sample}."
done

echo "All mappings and featureCounts completed."
 
