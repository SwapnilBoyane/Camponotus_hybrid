#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=genotype
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-135


# Prepend bcftools path to PATH
export PATH=~/bcftools/usr/local/bin:$PATH
source activate vcftools
module load gcc/10.1.0
module load r/4.3.0
R

# define main working directory
workdir=/lustre/scratch/sboyane/camphybrid

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames.txt | tail -n1 )

# define the reference genome
refgenome=/lustre/work/sboyane/ref/camp_sp_genome_filtered.fasta

# run bcftools to genotype
bcftools mpileup --skip-indels -C 0 -d 200 --min-MQ 10 --threads 4 -f ${refgenome} ${workdir}/01_bam_files/${basename_array}_final.bam | bcftools call -m --threads 4 -o ${workdir}/02_vcf/${basename_array}.vcf

# bgzip
~/anaconda3/bin/bgzip ${workdir}/02_vcf/${basename_array}.vcf

#tabix
~/anaconda3/bin/tabix ${workdir}/02_vcf/${basename_array}.vcf.gz

# filter individual vcf files
bcftools view -i 'MIN(DP)>5' ${workdir}/02_vcf/${basename_array}.vcf.gz > ${workdir}/03_vcf/${basename_array}.vcf

# bgzip
~/anaconda3/bin/bgzip ${workdir}/03_vcf/${basename_array}.vcf

#tabix
~/anaconda3/bin/tabix ${workdir}/03_vcf/${basename_array}.vcf.gz
########
# contam check
# extract all heterozygous sites for this individual
vcftools --gzvcf ${workdir}/03_vcf/${basename_array}.vcf.gz --min-alleles 2 --max-alleles 2 \
--maf 0.5 --recode --recode-INFO-all --out ${workdir}/03_contam/${basename_array}

# extract the depth info for all the sites retained
bcftools query -f '%DP4\n' ${workdir}/03_contam/${basename_array}.recode.vcf > ${workdir}/03_contam/${basename_array}.dp

# make a histogram of MAF sequencing depth proportion
Rscript contam_check.r ${workdir}/03_contam/${basename_array}.dp ${basename_array}
#########
# alignment stats
echo ${basename_array} > ${basename_array}.stats

# samtools depth sum of aligned sites
echo "samtools depth sum of aligned sites" >> ${basename_array}.stats
samtools depth  ${workdir}/01_bam_files/${basename_array}_final.bam  |  awk '{sum+=$3} END { print "Sum = ",sum}' >> ${basename_array}.stats

# proportion dupes
echo "proportion duplicates" >> ${basename_array}.stats
head -n8 ${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt | tail -n1 | cut -f9 >> ${basename_array}.stats

# number of genotyped sites passing minimum depth filter
echo "sites genotyped" >> ${basename_array}.stats
gzip -cd ${workdir}/03_vcf/${basename_array}.vcf.gz | grep -v "^#" | wc -l >> ${basename_array}.stats

