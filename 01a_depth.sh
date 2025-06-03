#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=genotype
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-183

#load modules
conda activate samtools
module load gcc/10.1.0
module load r/4.3.0
R

# Set working directory
workdir=/lustre/scratch/sboyane/camphybrid

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames.txt | tail -n1 )

# Calculate the depth of coverage per site
samtools depth -a ${workdir}/01_bam_files/${basename_array}_final.bam > ${workdir}/depth/${basename_array}_depth.txt

# Isolate the third column with the read counts only to reduce file size
cut -f3 ${workdir}/depth/${basename_array}_depth.txt > ${workdir}/depth/${basename_array}_depth_reduced.txt

# Remove the full depth file
rm ${workdir}/depth/${basename_array}_depth.txt

# Make sure to keep the ".R" script in same folder of this "".sh" script

Rscript plot_coverage.R

