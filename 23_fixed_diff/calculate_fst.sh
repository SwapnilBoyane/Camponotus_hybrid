#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G


# Prepend bcftools path to PATH
export PATH=~/bcftools/usr/local/bin:$PATH
source activate vcftools


# define main working directory
workdir=/lustre/scratch/sboyane/camphybrid/23_fixed_diff


#run per-site Fst pure Herc vs other_ME_MW_Nov

vcftools --vcf Fst_all_pureI_camphybrid_mac10.vcf --weir-fst-pop pure_Her.txt --weir-fst-pop other_ME_MW_Nov.txt --out Fst_output_herc_vs_other

#run per-site Fst pure modoc_est vs other_MW_her_Nov

vcftools --vcf Fst_all_pureI_camphybrid_mac10.vcf --weir-fst-pop pure_ME.txt --weir-fst-pop other_MW_Her_Nov.txt --out Fst_output_ME_vs_other

#run per-site Fst pure modoc_west vs other_ME_her_Nov

vcftools --vcf Fst_all_pureI_camphybrid_mac10.vcf --weir-fst-pop pure_MW.txt --weir-fst-pop other_ME_Her_Nov.txt --out Fst_output_MW_vs_other

#run per-site Fst pure nova vs other_ME_MW_her

vcftools --vcf Fst_all_pureI_camphybrid_mac10.vcf --weir-fst-pop pure_Nov.txt --weir-fst-pop other_ME_MW_Her.txt --out Fst_output_Nov_vs_other

