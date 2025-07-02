#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=fil_Fst
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-31

export PATH=~/bcftools/usr/local/bin:$PATH
source activate vcftools

# define main working directory
workdir=/lustre/scratch/sboyane/camphybrid

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

#Filter for Fst of pure individuals 

vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keeplist_pure_parental.txt --max-missing 0.8  --min-alleles 2 --max-alleles 2 --mac 10 --max-maf 0.49 \
--remove-indels --recode --recode-INFO-all --out ${workdir}/23_fixed_diff/${region_array}_mac10_Fst_PP_ind

#Filter for Fst of pure and admix individuals 

vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep all_PP_admix.txt --max-missing 0.8  --min-alleles 2 --max-alleles 2 --mac 10 --max-maf 0.49 \
--remove-indels --recode --recode-INFO-all --out ${workdir}/23_fixed_diff/${region_array}_mac10_Fst_all_ind
