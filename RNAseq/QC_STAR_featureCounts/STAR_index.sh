#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8            
#SBATCH --mem=64G                    
#SBATCH --time=12:00:00              

#   This SLURM batch script generates a genome index for STAR RNA-seq alignment.
#   It uses a reference genome FASTA file and a GTF annotation file to create the STAR index.
#
# Requirements:
#   - SLURM job scheduler
#   - STAR module installed and loaded
#   - Reference genome FASTA file (e.g., GRCh38.p14.genome.fa)
#   - Corresponding GTF annotation file (e.g., gv44_basicannotation.gtf)
#
# Parameters:
#   --genomeDir: Output directory for the STAR genome index
#   --genomeFastaFiles: Path to the reference genome FASTA
#   --sjdbGTFfile: GTF file with gene annotations
#   --sjdbOverhang: Read length minus 1 (e.g., 149 for 150 bp reads)
#
# Output:
#   - Indexed genome files under the specified --genomeDir
#

# Load STAR module 
module load star

# Set the paths and parameters
# genome_dir="/path/to/your/genome_directory"
# genome_fasta="/path/to/your/genome.fa"
# gtf_file="/path/to/your/annotations.gtf"

# Create the genome index
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir Ref/index  \
     --genomeFastaFiles GRCh38.p14.genome.fa \
     --sjdbGTFfile gv44_basicannotation.gtf  \
     --sjdbOverhang 149  # readlength-1
