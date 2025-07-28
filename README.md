# IMPACT (Integrated Mutation Profiling and Cancer Immunoanalysis Tool)
A scalable and modular Nextflow pipline for somatic variant calling, copy number alteration analysis and immune analysis of tumours and matched normal samples using whole-exome sequencing data.

## Included tools
This pipeline contains several steps and tools as below:

__Alignemnt and base recalibration__
* Cutadapt
* BWA-MEM
* GATK4 and Picard

__Variant calling__
* MuTect2
* Strelka2 and Manta
* VarScan2
* VarDict

__Total on-target and off-tagreg CNA__
* CNVkit

__Allele-specific CNAs, Ploidy and Purity__
* FACETS
* ASCAT

__MSI__
* MSIsensor Pro

__HLA typing and functionality__
* Polysolver
* LOHHLA

__T cells fraction__
* TcellExTRECT

## Features
It is recommended to run the pipeline on a High Performance Computing Cluster which has enough resources to support parallele processes. Also, this version is used slurm to run on HPC clusters and apptainer (previosuly known as Singularity) as environment of runnin tools. This version of pipeline only support GRCh38 and you need to make sure that you have below tools installed in your environment:

* Nextflow (23.04.1 or more)
* apptainer
* cutadapt (1.18 or more)
* Samtools (1.9 or more)
* BWA
* htslib
* tabix
* bcftools (1.9 or more)
* vcftools

Also, you need to provide below resources and modify their directories in `params.config` file:

* FASTA file for Ref genome (GRCh38)
* Index and dictionary file for Ref genome possibly using `samtools faidx` and `picard CreateSequenceDictionary`
* BWA index files (using `bwa index reference.fasta`)
* Baits and regiond `.bed` file (currently supporting Agilent v6)
* dbsnp `.vcf.gz` and index files
* Mills and 1000G gold standard indels for hg38 (`.vcf.gz` and index)
* Hg38 known indels `.vcf.gz` and index file
* af-only-gnomad.hg38 and 1000g_pon.hg38 (Available on Gatk Resource Bundle hg38)

## Input
This pipeline support paired-end reads for all your samples in a single `.csv` file with below format:

| patientID      | read1                       | read2                       | sampleType      |
|----------------|-----------------------------|-----------------------------|-----------------|
| Patient 1      | /to/dir/normal_R1.fastq.gz  | /to/dir/normal_R2.fastq.gz  | normal          |
| Patient 1      | /to/dir/tumor_R1.fastq.gz   | /to/dir/tumor_R2.fastq.gz   | tumor           |
| Patient 2      | /to/dir/normal_R1.fastq.gz  | /to/dir/normal_R2.fastq.gz  | normal          |
| Patient 2      | /to/dir/tumor_R1.fastq.gz   | /to/dir/tumor_R2.fastq.gz   | tumor           |
| Patient x      | /to/dir/normal_R1.fastq.gz  | /to/dir/normal_R2.fastq.gz  | normal          |
| Patient x      | /to/dir/tumor_R1.fastq.gz   | /to/dir/tumor_R2.fastq.gz   | tumor           |

# How to run
To properly run the pipeline based on your required parameters, it is recommended to review `.config` files and change the arguments of each tool if needed.

To run the pipeline you only need to run the below command in bash:

```
pipeDir=/to/dir/IMPACT_pipeline

# Run data
nextflow run ${pipeDir}/IMPACT.nf \
    -config ${pipeDir}/params.config \
    --batchFile /scratch/users/nhemmat/MultiVarCallernf/batch.csv \
    --tmpDir /scratch/users/nhemmat/tmp_pipeline/tmp \
    --outputDir /researchers/nima.hemmat/Brior_IMPACT_new \
    --runCNV true \
    -work-dir /scratch/users/nhemmat/tmp_pipeline \
    -resume
```
`--batchFile` should be set to the directotry of `.csv` input file.
`--tmpDir` is directory of intermediate files created during the analaysis and should have enough space (more than 100 GB) based on the number of input samples.
`--outputDir` is directory of final output files and should have enough space too.
`--runCNV` indicates if you need to have CNV analysis in your outputs (recommended).

## Author
Nima Hemmat - Peter MacCallum Cancer Centre, VIC, Asutralia

## License
MIT License – See LICENSE file for details.

## Contact
For questions or support, contact:
📧 nima.hemmat@petermac.org