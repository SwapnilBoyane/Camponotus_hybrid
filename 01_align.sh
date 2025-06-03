#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam_processing
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=8
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-183

#use module to load tools or use direct file path
source activate samtools

# define main working directory
workdir=/lustre/scratch/sboyane/camphybrid

#define file for array job
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames.txt | tail -n1 )

# define the reference genome
refgenome=/lustre/work/sboyane/ref/camp_sp_genome_filtered.fasta

# define the location of the reference mitogenomes
mito=/home/jmanthey/denovo_genomes/formicinae_mitogenomes.fasta

# define the location of the reference blochmannia genomes
bloch=/lustre/work/sboyane/ref/bloch/C019_modoc.fasta

# run bbduk
/lustre/work/jmanthey/bbmap/bbduk.sh in1=${workdir}/00_fastq/${basename_array}_R1.fastq.gz in2=${workdir}/00_fastq/${basename_array}_R2.fastq.gz out1=${workdir}/01_cleaned/${basename_array}_R1.fastq.gz out2=${workdir}/01_cleaned/${basename_array}_R2.fastq.gz minlen=50 ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 ref=/lustre/work/jmanthey/bbmap/resources/adapters.fa hdist=1 tbo tpe

# run bbsplit mitogenomes
/lustre/work/jmanthey/bbmap/bbsplit.sh in1=${workdir}/01_cleaned/${basename_array}_R1.fastq.gz in2=${workdir}/01_cleaned/${basename_array}_R2.fastq.gz ref=${mito} basename=${workdir}/01_mtDNA/${basename_array}_%.fastq.gz outu1=${workdir}/01_mtDNA/${basename_array}_R1.fastq.gz outu2=${workdir}/01_mtDNA/${basename_array}_R2.fastq.gz

# remove unnecessary bbsplit output files
rm ${workdir}/01_mtDNA/${basename_array}_R1.fastq.gz
rm ${workdir}/01_mtDNA/${basename_array}_R2.fastq.gz

#Keep selected samples to process for blochnammia (this step is done because a few samples have separate head and gaster FASTA files)
basename_array2=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames2.txt | tail -n1 )

# run bbsplit blochmannia
/lustre/work/jmanthey/bbmap/bbsplit.sh in1=${workdir}/01_cleaned/${basename_array2}_R1.fastq.gz in2=${workdir}/01_cleaned/${basename_array2}_R2.fastq.gz ref=${bloch} basename=${workdir}/01_blochmannia/${basename_array2}_%.fastq.gz outu1=${workdir}/01_blochmannia/${basename_array2}_R1.fastq.gz outu2=${workdir}/01_blochmannia/${basename_array2}_R2.fastq.gz

# remove unnecessary bbsplit output files
rm ${workdir}/01_blochmannia/${basename_array2}_R1.fastq.gz
rm ${workdir}/01_blochmannia/${basename_array2}_R2.fastq.gz

# run bwa mem
/home/sboyane/anaconda3/bin/bwa mem -t 8 ${refgenome} ${workdir}/01_cleaned/${basename_array}_R1.fastq.gz ${workdir}/01_cleaned/${basename_array}_R2.fastq.gz > ${workdir}/01_bam_files/${basename_array}.sam

# convert sam to bam
/home/sboyane/anaconda3/bin/samtools view -b -S -o ${workdir}/01_bam_files/${basename_array}.bam ${workdir}/01_bam_files/${basename_array}.sam

# remove sam
rm ${workdir}/01_bam_files/${basename_array}.sam

# clean up the bam file # used Picard.jar for sorting and cleaning bam
/usr/bin/java -jar /lustre/scratch/sboyane/camphybrid/01_bam_files/picard.jar CleanSam I=${workdir}/01_bam_files/${basename_array}.bam O=${workdir}/01_bam_files/${basename_array}_cleaned.bam

# remove the raw bam
rm ${workdir}/01_bam_files/${basename_array}.bam

# sort the cleaned bam file
/usr/bin/java -jar /lustre/scratch/sboyane/camphybrid/01_bam_files/picard.jar SortSam I=${workdir}/01_bam_files/${basename_array}_cleaned.bam O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam SORT_ORDER=coordinate

# remove the cleaned bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned.bam

# add read groups to sorted and cleaned bam file
/usr/bin/java -jar /lustre/scratch/sboyane/camphybrid/01_bam_files/picard.jar AddOrReplaceReadGroups I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam RGLB=1 RGPL=illumina RGPU=unit1 RGSM=${basename_array}

# remove cleaned and sorted bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam

# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file)
/usr/bin/java -jar /lustre/scratch/sboyane/camphybrid/01_bam_files/picard.jar MarkDuplicates REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 M=${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam O=${workdir}/01_bam_files/${basename_array}_final.bam

# remove sorted, cleaned, and read grouped bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam

# index the final bam file
/home/sboyane/anaconda3/bin/samtools index ${workdir}/01_bam_files/${basename_array}_final.bam

