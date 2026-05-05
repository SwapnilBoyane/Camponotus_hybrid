#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=cal_AF
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=64G


export PATH=~/bcftools/usr/local/bin:$PATH
source activate vcftools

# define main working directory
workdir=/lustre/scratch/sboyane/camphybrid/12_digno_snps/

for pop in pure_Her pure_MS pure_MN pure_MW pure_Nov hyb_lin; do
 vcftools --vcf combined_AF_PP_143_digno_SNWHHN.vcf --keep ${pop}.txt --freq2 --out ${pop}_final
done
