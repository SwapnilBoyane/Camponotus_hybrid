#!/bin/bash
#SBATCH --chdir=.
#SBATCH --job-name=genotype
#SBATCH --partition=nocona
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

cd /lustre/scratch/sboyane/camphybrid/23_fixed_diff/allele_freq

bcftools view -R diag_all.bed AF_all_pureI_and_admix_camphybrid_mac10.vcf.gz -Oz -o diagnostic_snps_camphybrid_AF.vcf.gz
