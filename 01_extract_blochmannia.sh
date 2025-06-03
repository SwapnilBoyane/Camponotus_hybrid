#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=8
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-173

source activate samtools

# define main working directory
workdir=/lustre/scratch/sboyane/camphybrid

# run below for samples C064 to C084
#basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/extract_bloch9.txt | tail -n1 )

#only 173 samples processed 
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames2.txt | tail -n1 )

# define the reference genome
refgenome=/lustre/work/sboyane/ref/camp_sp_genome_filtered.fasta

# define the location of the reference mitogenomes
mito=/home/jmanthey/denovo_genomes/formicinae_mitogenomes.fasta

# define the location of the reference blochmannia genomes
bloch=/lustre/work/sboyane/ref/bloch/C019_modoc.fasta

# run bbduk
#/lustre/work/jmanthey/bbmap/bbduk.sh in1=${workdir}/00_fastq/${basename_array}_R1.fastq.gz in2=${workdir}/00_fastq/${basename_array}_R2.fastq.gz out1=${workdir}/01_cleaned/new/${basename_array}_R1.fastq.gz out2=${workdir}/01_cleaned/new/${basename_array}_R2.fastq.gz minlen=50 ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 ref=/lustre/work/jmanthey/bbmap/resources/adapters.fa hdist=1 tbo tpe

# run bbsplit mitogenomes
#/lustre/work/jmanthey/bbmap/bbsplit.sh in1=${workdir}/01_cleaned/${basename_array}_R1.fastq.gz in2=${workdir}/01_cleaned/${basename_array}_R2.fastq.gz ref=${mito} basename=${workdir}/01_mtDNA/${basename_array}_%.fastq.gz outu1=${workdir}/01_mtDNA/${basename_array}_R1.fastq.gz outu2=${workdir}/01_mtDNA/${basename_array}_R2.fastq.gz

# remove unnecessary bbsplit output files
#rm ${workdir}/01_mtDNA/${basename_array}_R1.fastq.gz
#rm ${workdir}/01_mtDNA/${basename_array}_R2.fastq.gz

# run bbsplit blochmannia
/lustre/work/jmanthey/bbmap/bbsplit.sh in1=${workdir}/01_cleaned/new/${basename_array}_R1.fastq.gz in2=${workdir}/01_cleaned/new/${basename_array}_R2.fastq.gz ref=${bloch} basename=${workdir}/01_blochmannia/${basename_array}_%.fastq.gz outu1=${workdir}/01_blochmannia/${basename_array}_R1.fastq.gz outu2=${workdir}/01_blochmannia/${basename_array}_R2.fastq.gz

# remove unnecessary bbsplit output files
rm ${workdir}/01_blochmannia/${basename_array}_R1.fastq.gz
rm ${workdir}/01_blochmannia/${basename_array}_R2.fastq.gz

# Extract the reads from the interleaved file and separate them into R1 and R2 files
seqtk seq -1 ${workdir}/01_blochmannia/${basename_array}_C019_modoc.fastq.gz > ${workdir}/01_blochmannia/01_cleaned/${basename_array}_R1.fastq
seqtk seq -2 ${workdir}/01_blochmannia/${basename_array}_C019_modoc.fastq.gz > ${workdir}/01_blochmannia/01_cleaned/${basename_array}_R2.fastq

# Compress the files to .gz
gzip ${workdir}/01_blochmannia/01_cleaned/${basename_array}_R1.fastq
gzip ${workdir}/01_blochmannia/01_cleaned/${basename_array}_R2.fastq

