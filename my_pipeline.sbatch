#!/bin/bash
#SBATCH --account=loi
#SBATCH --job-name=IMPACT
#SBATCH --partition=janis
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nima.hemmat@petermac.org
  
module purge
module load nextflow/23.04.1
module load cutadapt/1.18
module load samtools/1.9
module load bwa
module load htslib
module load tabix
module load bcftools/1.9
module load vcftools

pipeDir=/scratch/users/nhemmat/MultiVarCallernf

# Run data
nextflow run ${pipeDir}/MultiVarCallernf.nf \
    -profile debug \
    -config ${pipeDir}/params.config \
    --batchFile /scratch/users/nhemmat/MultiVarCallernf/batch.csv \
    --tmpDir /scratch/users/nhemmat/tmp_pipeline/tmp \
    --outputDir /researchers/nima.hemmat/Brior_IMPACT_new \
    --runCNV true \
    -work-dir /scratch/users/nhemmat/tmp_pipeline \
    -resume
