#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=cal_AF
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=32G

source activate vcftools

# define main working directory
workdir=/lustre/scratch/sboyane/camphybrid/23_fixed_diff/allele_freq

# calculate per lineage allele frequency 
# keep pop files in working directory 

for pop in pure_Her pure_ME pure_MW pure_Nov; do
vcftools --gzvcf combined_PP_snps.vcf.gz --keep ${pop}.txt --freq2 --out ${pop}
done


