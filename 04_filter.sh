#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
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

# filter based on missing data for each of the subsets of data

vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keeplist.txt  --max-missing 1.0 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --mac 2 --thin 20000 --remove-indels --recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/${region_array}


/home/sboyane/anaconda3/bin/bgzip ${workdir}/05_filtered_vcf/${region_array}.recode.vcf
/home/sboyane/anaconda3/bin/tabix ${workdir}/05_filtered_vcf/${region_array}.recode.vcf.gz
