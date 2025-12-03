#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=mrbayes
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-750

source activate bcftools
module load gcc/14.2.0
module load R/4.4.1

# Set the number of runs that each SLURM task should do
PER_TASK=9

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

# Run the loop of runs for this task.
for (( run=$START_NUM; run<=$END_NUM; run++ )); do
	echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run

	chrom_array=$( head -n${run} stat_helper.txt | cut -f1 | tail -n1 )

	start_array=$( head -n${run} stat_helper.txt | cut -f2 | tail -n1 )

	end_array=$( head -n${run} stat_helper.txt | cut -f3 | tail -n1 )

	gunzip -cd /lustre/scratch/sboyane/camphybrid/11_mrbayes/scaffold0031.recode.vcf.gz | grep "#" > /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf

	~/anaconda3/bin/tabix /lustre/scratch/sboyane/camphybrid/11_mrbayes/${chrom_array}.recode.vcf.gz ${chrom_array}:${start_array}-${end_array} >> /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf

	Rscript _write_mrbayes.r /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf

	rm /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf

	conda deactivate

	/home/sboyane/anaconda3/envs/mrbayes/bin/mb /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.nex

	rm /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.nex
	rm /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.nex.ckp
	rm /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.nex.ckp~
	rm /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.nex.mcmc
	rm /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.nex.parts
	rm /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.nex.trprobs
	rm /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.nex.tstat
	rm /lustre/scratch/sboyane/camphybrid/11_mrbayes/windows/${chrom_array}__${start_array}__${end_array}.nex.vstat

done
