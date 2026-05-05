#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=fil_AF
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-31

export PATH=~/bcftools/usr/local/bin:$PATH
source activate vcftools

# define main working directory
workdir=/lustre/scratch/sboyane/camphybrid

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

#Filter for Fst of pure individuals 
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keeplist_parental_diagno_final_143.txt --max-missing 0.8  --min-alleles 2 --max-alleles 2 --mac 10 \
--remove-indels --recode --recode-INFO-all --out ${workdir}/12_digno_snps/${region_array}_mac10_AF_PP_ind_final_SNWHHN

#Filter for Fst of pure and admix individuals 
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keeplist_194.txt  --max-missing 0.8  --min-alleles 2 --max-alleles 2 --mac 2  \
--remove-indels --recode --recode-INFO-all --out ${workdir}/12_digno_snps/${region_array}_mac2_AF_all_ind_final

