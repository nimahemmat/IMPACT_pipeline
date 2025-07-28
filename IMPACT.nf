#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import groovy.json.JsonSlurper
import java.nio.file.Files

log.info ""
log.info " NEXTFLOW ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
log.info "-------------------------------------------------------------------------"
log.info "          Integrated Somatic Variant Calling and Immuno-Oncogenomic Analysis          "
log.info "                    of Tumors from Whole-Exome Sequencing Data                    "
log.info "         IMPACT (Integrated Mutation Profiling and Cancer Immunoanalysis Tool)         "
log.info "-------------------------------------------------------------------------"
log.info ""

/*
________________________________________________________________________________

                            C O N F I G U R A T I O N
________________________________________________________________________________
*/

params.batchFile = null // Define parameter for batchFile


// Check if we got a batch file
if (params.batchFile == null) {
    log.error "No sample sheet specified. Please use --batchFile to pass the sample sheet."
    exit 1
}

// Read and parse the batch file
def batchCSV = file(params.batchFile).splitCsv(header: true)

def validFQfields = ["patientID", "read1", "read2", "sampleType"]
def validSampleTypes = ["tumor", "normal"]

// Validate batch file
if (batchCSV.size() == 0) {
    exit 1, "Error: Batch file is empty. Please check your batchFile."
}

// Check for correct header fields
if (!batchCSV[0].keySet().sort().equals(validFQfields.sort())) {
    exit 1, "Error: Incorrect fields in batch file. Expected fields: ${validFQfields.join(', ')}."
}

// Validate rows in the batch file
def invalidRows = batchCSV.findAll { row ->
    !(row.sampleType.toLowerCase() in validSampleTypes) || !file(row.read1).exists() || !file(row.read2).exists()
}

if (invalidRows) {
    def errorMessages = invalidRows.collect { row ->
        def errors = []
        if (!(row.sampleType.toLowerCase() in validSampleTypes)) {
            errors << "Invalid sampleType '${row.sampleType}' for patientID '${row.patientID}'"
        }
        if (!file(row.read1).exists() || !file(row.read2).exists() ) {
            errors << "Missing FASTQ file for patientID '${row.patientID}'"
        }
        errors.join('; ')
    }
    exit 1, "Error: Issues detected in batch file:\n${errorMessages.join('\n')}"
}

// Check that each patientID has both tumor and normal BAM files
def patientGroups = batchCSV.groupBy { it.patientID }
def incompletePatients = patientGroups.findAll { patientID, rows ->
    // Collect sample types for the patient and convert them to lowercase
    def sampleTypes = rows*.sampleType.collect { it.toLowerCase() }
    !(sampleTypes.contains("tumor") && sampleTypes.contains("normal"))
}

if (incompletePatients) {
    def errorMessages = incompletePatients.collect { patientID, rows ->
        def sampleTypes = rows*.sampleType.join(", ")
        "PatientID '${patientID}' is missing either a tumor or normal sample (found: ${sampleTypes})"
    }
    exit 1, "Error: Some patients are missing tumor or normal FASTQ files:\n${errorMessages.join('\n')}"
}

// Log the validated batch file information
log.info "Batch file validated successfully. Number of rows: ${batchCSV.size()}"
batchCSV.each { row ->
    log.info "PatientID: ${row.patientID}, read1: ${row.read1}, read2: ${row.read2}, SampleType: ${row.sampleType}"
}


// Collect raw data
def raw_data = []

def t_map = [:]
def n_map = [:]


for (row in batchCSV) {
    def meta = [:]  // Create a map to store metadata
    def reads = []   // List to store FASTQ file paths

    // Populate metadata
    meta.patientID = row.patientID.toString()
    meta.sampleType = row.sampleType.toLowerCase()

    t_map[meta.patientID] = (t_map[meta.patientID] == null) ? 0 : t_map[meta.patientID]
    n_map[meta.patientID] = (n_map[meta.patientID] == null) ? 0 : n_map[meta.patientID]

    if (row.sampleType in validSampleTypes) {
        t_map[meta.patientID] += (row.sampleType == "tumor") ? 1 : 0
        n_map[meta.patientID] += (row.sampleType == "normal") ? 1 : 0
    } else {
        exit 1, "Error: Incorrect sampleType [got: " + row.sampleType + " - allowed: " + validSampleTypes + "]: please check your batchFile"
    }

    if (row.read1) { 
    reads.add(file(row.read1, checkIfExists: true)) 
    }
    if (row.read2) {
    reads.add(file(row.read2, checkIfExists: true))
    }

    raw_data.add([meta, reads])
}



// check if we have T/N DNA pairs for all patients
raw_data.each {
    record ->
    if (t_map[record[0].patientID] == 0) {
        exit 1, "NO tumor DNA sample specified for: " + record[0].patientID
    }
    if (n_map[record[0].patientID] == 0) {
        exit 1, "NO normal DNA sample specified for: " + record[0].patientID
    }
}

// Make Channel of raw data List
batch_raw_data_ch = Channel.fromList(raw_data)


batch_raw_data_ch.map {
        meta, fastq ->
        [ meta, fastq ]
    }
    .groupTuple(by: [0]) // Group by meta (patientID and sampleType)
    .branch {
        meta, fastq ->
            single  : fastq.unique().size() == 1
            multiple: fastq.unique().size() > 1
    }
    .set { fastq_ch } // Assign the result to fastq_ch

// Merge FASTQ files for multiple pairs
process merge_fastq {

    label 'low_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(r1_files, stageAs: "?/*"), 
            path(r2_files, stageAs: "?/*")

    output:
    tuple val(meta), path("${meta.patientID}_${meta.sampleType}_R1.merged.fastq.gz"), path("${meta.patientID}_${meta.sampleType}_R2.merged.fastq.gz")

    script:
    def prefix = "${meta.patientID}_${meta.sampleType}"
    def input_str_R1 = r1_files instanceof List ? r1_files.join(" ") : r1_files
    def input_str_R2 = r2_files instanceof List ? r2_files.join(" ") : r2_files
    """
    cat ${input_str_R1} > ${prefix}_R1.merged.fastq.gz
    cat ${input_str_R2} > ${prefix}_R2.merged.fastq.gz
    """
}

// Cut adaptors and make trimmed fastq files
process cut_adaptors {

    label 'med_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.patientID}_${meta.sampleType}_R1_trimmed.fastq.gz"), path("${meta.patientID}_${meta.sampleType}_R2_trimmed.fastq.gz")

    script:
    """
    cutadapt -j ${task.cpus} \\
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \\
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \\
        -o ${meta.patientID}_${meta.sampleType}_R1_trimmed.fastq.gz \\
        -p ${meta.patientID}_${meta.sampleType}_R2_trimmed.fastq.gz \\
        "$r1" "$r2"
    """
}

// Mak unaligned BAM
process 'MakeUbam' {

    label 'picard'
    label 'high_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path(ubam)

    script:
    ubam = meta.patientID + "_" + meta.sampleType + "_unaligned.bam"
    def read_group = meta.patientID + "_" + meta.sampleType
    def java_opt = '-Xmx' + task.memory.toGiga() + "G"
    """
    mkdir -p ${params.tmpDir}
    picard ${java_opt} FastqToSam \\
        --TMP_DIR ${params.tmpDir} \\
        --ALLOW_AND_IGNORE_EMPTY_LINES true \\
        --MAX_RECORDS_IN_RAM 16777216 \\
        --FASTQ ${r1} \\
        --FASTQ2 ${r2} \\
        --READ_GROUP_NAME ${read_group} \\
        --SAMPLE_NAME ${read_group} \\
        --LIBRARY_NAME ${read_group} \\
        --PLATFORM ILLUMINA \\
        --OUTPUT ${ubam}
    """
}

// Align to reference genome

process 'BwaAlign' {

    label 'high_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(r1), path(r2), path(refFasta), path(refIdx), path(refDict), path(bwaInd)

    output:
    tuple val(meta), path(bam)

    script:
    bam = meta.patientID + "_" + meta.sampleType + "_aligned_sorted.bam"
    def read_group = meta.patientID + "_" + meta.sampleType
    // def sam_sort_mem = Math.max((task.memory.toGiga() - 50), 4) + "G"
    def sort_threads = (task.cpus.compareTo(8) == 1) ? 8 : task.cpus
    """
    bwa mem \\
        -R "@RG\\tID:${read_group}\\tLB:${read_group}\\tSM:${read_group}\\tPL:ILLUMINA" \\
        -M \\
        -t ${task.cpus} \\
        -Y \\
        ${refFasta} ${r1} ${r2} | \\
    samtools view -@2 -Shbu - | \\
    samtools sort \\
        -l 6 \\
        -m 8G \\
        -o ${bam} \\
        -@ ${sort_threads} \\
        /dev/stdin

    """
}

// mrge aligned and unaligned BAMs

process 'MergeBams' {

    label 'picard'
    label 'med_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(bam), path(ubam), path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path("${meta.patientID}_${meta.sampleType}_aligned_uBAM_merged.bam")

    script:
    def java_opt = '-Xmx' + task.memory.toGiga() + "G"
    """
    mkdir -p ${params.tmpDir}
    picard ${java_opt} MergeBamAlignment \\
        --TMP_DIR ${params.tmpDir} \\
        --VALIDATION_STRINGENCY SILENT \\
        --EXPECTED_ORIENTATIONS FR \\
        --ATTRIBUTES_TO_RETAIN X0 \\
        --REFERENCE_SEQUENCE ${refFasta} \\
        --PAIRED_RUN true \\
        --SORT_ORDER "queryname" \\
        --IS_BISULFITE_SEQUENCE false \\
        --ALIGNED_READS_ONLY false \\
        --CLIP_ADAPTERS false \\
        --MAX_RECORDS_IN_RAM 16777216 \\
        --ADD_MATE_CIGAR true \\
        --MAX_INSERTIONS_OR_DELETIONS -1 \\
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \\
        --UNMAPPED_READ_STRATEGY COPY_TO_TAG \\
        --ALIGNER_PROPER_PAIR_FLAGS true \\
        --UNMAP_CONTAMINANT_READS true \\
        --ALIGNED_BAM ${bam} \\
        --UNMAPPED_BAM ${ubam} \\
        --OUTPUT ${meta.patientID}_${meta.sampleType}_aligned_uBAM_merged.bam
    """
}

// Mark Duplicates
process 'MarkDuplicates' {

    label 'picard'

    label 'medhigh_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(mkdp_bam)

    script:
    mkdp_bam = meta.patientID + "_" + meta.sampleType + "_mkdp_unsorted.bam"
    def java_opt = '-Xmx' + task.memory.toGiga() + "G"
    """
    mkdir -p ${params.tmpDir}
    picard ${java_opt} MarkDuplicates \\
        -I ${bam} \\
        -M /dev/null \\
        -O ${mkdp_bam} \\
        --MAX_RECORDS_IN_RAM 16777216 \\
        --TMP_DIR ${params.tmpDir} \\
        --VALIDATION_STRINGENCY LENIENT
    """

}

// Sorting mkdp bams
process 'SamtoolsSortDups' {
    
    label 'medhigh_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(mkdp_bam)

    output:
    tuple val(meta), path(mkdp_sorted_bam)

    script:
    mkdp_sorted_bam = meta.patientID + "_" + meta.sampleType + "_mkdp_sorted.bam"
    def sam_sort_mem = task.memory.toGiga() + "G"
    """
    samtools sort \\
        -l 6 \\
        -m 8G \\
        -o ${mkdp_sorted_bam} \\
        -@${task.cpus} \\
        ${mkdp_bam}
    """
}

//Add Nm Md Uq tags
process 'SetNmMdAndUqTags' {

    label 'picard'

    label 'med_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(mkdp_sorted_bam), path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path(bam_out)

    script:
    def procSampleName = meta.patientID + "_" + meta.sampleType
    bam_out = [procSampleName + "_aligned_sort_mkdp.bam", procSampleName + "_aligned_sort_mkdp.bai"]
    def java_opt = '-Xmx' + task.memory.toGiga() + "G"
    """
    mkdir -p ${params.tmpDir}
    picard ${java_opt} SetNmMdAndUqTags \\
        -I ${mkdp_sorted_bam} \\
        -O ${procSampleName}_aligned_sort_mkdp.bam \\
        -R ${refFasta} \\
        --CREATE_INDEX true \\
        --MAX_RECORDS_IN_RAM 8388608 \\
        --TMP_DIR ${params.tmpDir} \\
        --VALIDATION_STRINGENCY LENIENT
    """
}

// Make Interval list
process 'RegionsBedToIntervalList' {

    label 'gatk'

    tag 'RegionsBedToIntervalList'

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/", mode: 'copy'

    input:
    tuple path(refDict), path(regionsBed)

    output:
    path("${regionsBed.baseName}.interval_list")

    script:
    def java_opt = '-Xmx' + task.memory.toGiga() + "G"
    """
    gatk --java-options ${java_opt} BedToIntervalList \\
        -I ${regionsBed} \\
        -O ${regionsBed.baseName}.interval_list \\
        -SD ${refDict}
    """
}

process 'BaitsBedToIntervalList' {

    label 'gatk'

    tag 'BaitsBedToIntervalList'

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/", mode: 'copy'

    input:
    tuple path(refDict), path(baitsBed)

    output:
    path("${baitsBed.baseName}.interval_list")

    script:
    def java_opt = '-Xmx' + task.memory.toGiga() + "G"
    """
    gatk --java-options ${java_opt} BedToIntervalList \\
        -I ${baitsBed} \\
        -O ${baitsBed.baseName}.interval_list \\
        -SD ${refDict}
    """
}

process 'preprocessIntervalList' {

    label 'gatk'

    label 'low_memory'

    tag 'preprocessIntervalList'

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/", mode: 'copy'

    input:
    tuple path(interval_list), path(refFasta), path(refIdx), path(refDict)

    output:
    path(outFileName)

    script:
    outFileName = interval_list.baseName + "_merged_padded.interval_list"
    def java_opt = '-Xmx' + task.memory.toGiga() + "G"
    """
    gatk --java-options ${java_opt} PreprocessIntervals \\
        -R ${refFasta} \\
        -L ${interval_list} \\
        --bin-length 0 \\
        --padding 250 \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -O ${outFileName}
    """
}

// Split interval files in 20 files for scattering Mutect2 and GATK Base recalibration
process 'SplitIntervals' {
    
    label 'gatk'

    label 'low_memory'

    tag 'SplitIntervals'

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/SplitIntervals/", mode: 'copy'

    input:
    tuple path(interval_list), path(refFasta), path(refIdx), path(refDict)

    output:
    path(intervals)

    script:
    def intervalDir = interval_list.baseName
    intervals = "${intervalDir}/*"
    """
    mkdir -p ${params.tmpDir}

    gatk SplitIntervals \\
        --tmp-dir ${params.tmpDir} \\
        -R ${refFasta}  \\
        -scatter ${params.ScatterNum} \\
        --interval-merging-rule ALL \\
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \\
        -L ${interval_list} \\
        -O ${intervalDir}
    """
}

// GATK Base Recalibration pipeline

process 'scatterBaseRecalibrator' {

    label 'gatk'

    label 'med_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val (meta), path(bam_out), path(intervals), 
    path(refFasta), path(refIdx), path(refDict),
    path(dbsnp), path(dbsnpIdx), path(millsgold),
    path(millsgoldIdx), path(knownindels), path(knownindelsIdx)

    output:
    tuple val(meta), path("${procSampleName}_${intervals}_bqsr.table")

    script:
    procSampleName = meta.patientID + "_" + meta.sampleType
    def java_opt = '-Xmx' + task.memory.toGiga() + "G"
    """
    mkdir -p ${params.tmpDir}

    gatk --java-options ${java_opt} BaseRecalibrator \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam_out[0]} \\
        -R ${refFasta} \\
        -L ${intervals} \\
        -O ${procSampleName}_${intervals}_bqsr.table \\
        --known-sites ${dbsnp} \\
        --known-sites ${knownindels} \\
        --known-sites ${millsgold}
    """
}

process 'gatherscatteredBQSRtables' {

    label 'gatk'

    label 'med_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(bqsr_table)

    output:
    tuple val(meta), path("${procSampleName}_bqsr.table")

    script:
    procSampleName = meta.patientID + "_" + meta.sampleType

    """
    mkdir -p ${params.tmpDir}

    gatk GatherBQSRReports \\
        -I ${bqsr_table.join(" -I ")} \\
        -O ${procSampleName}_bqsr.table
    """
}

process 'scatterapplyBQSRS' {

    label 'gatk'
    label 'med_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val (meta), path(bam_out), path(bqsr_table), path(intervals), 
    path(refFasta), path(refIdx), path(refDict),
    path(dbsnp), path(dbsnpIdx), path(millsgold),
    path(millsgoldIdx), path(knownindels), path(knownindelsIdx)

    output:
    tuple val(meta), path(scatter_bam)

    script:
    def procSampleName = meta.patientID + "_" + meta.sampleType
    scatter_bam = [ procSampleName + "_" + intervals + "_recal4.bam",
                    procSampleName + "_" + intervals + "_recal4.bai" ]
    def java_opt = '-Xmx' + task.memory.toGiga() + "G"
    """
    mkdir -p ${params.tmpDir}

    gatk --java-options ${java_opt} ApplyBQSR \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam_out[0]} \\
        -R ${refFasta} \\
        -L ${intervals} \\
        -O ${procSampleName}_${intervals}_recal4.bam \\
        --bqsr-recal-file ${bqsr_table}
    """
}

process 'gatherRecalBamFiles' {

    label 'gatk'
    label 'high_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path(final_bam)

    script:
    def procSampleName = meta.patientID + "_" + meta.sampleType
    final_bam = "${procSampleName}_unsorted_recalibrated.bam"
    def JAVA_Xmx = "4G"
    def java_opt = '"-Xmx' + JAVA_Xmx + ' -XX:ParallelGCThreads=2"'
    """
    mkdir -p ${params.tmpDir}
    gatk --java-options ${java_opt} GatherBamFiles \\
        --TMP_DIR ${params.tmpDir} \\
        -I ${bam.join(" -I ")} \\
        -O ${final_bam}\\
        --CREATE_INDEX false \\
        --MAX_RECORDS_IN_RAM 10485760
    """
}

process 'SamtoolsSortIdxRecal' {

    label 'medhigh_memory'

    tag "$meta.patientID : $meta.sampleType"

    publishDir "${params.outputDir}/${meta.patientID}/01_processedBams/", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.patientID}_${meta.sampleType}_recal_sorted.bam"),
                     path("${meta.patientID}_${meta.sampleType}_recal_sorted.bam.bai")


    script:
    def procSampleName = meta.patientID + "_" + meta.sampleType
    """
    samtools sort \\
        -@${task.cpus} \\
        -m 8G \\
        -o ${procSampleName}_recal_sorted.bam \\
        ${bam}
    samtools index \\
        -@${task.cpus} \\
        ${procSampleName}_recal_sorted.bam
    """
}

// Alignment metrics and coverage analysis
process 'AlignmentMetrics' {

    label 'gatk'
    label 'med_memory'

    tag "$meta.patientID : $meta.sampleType"

    publishDir "${params.outputDir}/${meta.patientID}/QC", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(targets), path(baits),
                     path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path("${meta.patientID}_${meta.sampleType}.HS.metrics.txt"),
                     path("${meta.patientID}_${meta.sampleType}.perTarget.coverage.txt")

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opt = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    gatk --java-options ${java_opt} CollectHsMetrics \\
        --INPUT ${bam} \\
        --OUTPUT ${meta.patientID}_${meta.sampleType}.HS.metrics.txt \\
        -R ${refFasta} \\
        --BAIT_INTERVALS ${baits} \\
        --TARGET_INTERVALS ${targets} \\
        --PER_TARGET_COVERAGE ${meta.patientID}_${meta.sampleType}.perTarget.coverage.txt
    """
}

// Running Mutect2
process 'Mutect2' {

    label 'gatk'
    label 'med_memory'

    tag "$meta.patientID"

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai), path(intervals), 
        path(refFasta), path(refIdx), path(refDict), path(gnomAD), path(gnomADIdx),
        path(pon), path(ponIdx)

    output:
    tuple val(meta),
            path("${meta.patientID}_${intervals}.vcf.gz"),
            path("${meta.patientID}_${intervals}.vcf.gz.tbi"),
            path("${meta.patientID}_${intervals}.vcf.gz.stats"),
            path("${meta.patientID}_${intervals}-f1r2.tar.gz")

    script:
    def tumorName  = meta.patientID + "_tumor"
    def normalName = meta.patientID + "_normal"
    """
    mkdir -p ${params.tmpDir}

    gatk Mutect2 \\
        --tmp-dir ${params.tmpDir} \\
        -I ${tumorbam} -tumor ${tumorName} \\
        -I ${normalbam} -normal ${normalName} \\
        -R ${refFasta} \\
        --germline-resource ${gnomAD} \\
        --panel-of-normals ${pon} \\
        -L ${intervals} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --f1r2-tar-gz "${meta.patientID}_${intervals}-f1r2.tar.gz" \\
        -O "${meta.patientID}_${intervals}.vcf.gz"
    """
}

process 'GatherMutect2Calls' {

    label 'gatk'
    label 'med_momory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/02_Mutect2/raw_vcf", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats), path(f1r2_tar_gz)

    output:
    tuple(
        val(meta),
        path("${meta.patientID}_mutect2_raw.{vcf.gz,vcf.gz.tbi}"),
        path("${meta.patientID}_mutect2_raw.vcf.gz.stats"),
        path("${meta.patientID}_read-orientation-model.tar.gz")
    )

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opt = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    mkdir -p ${params.tmpDir}

    gatk --java-options ${java_opt} MergeVcfs \\
        --TMP_DIR ${params.tmpDir} \\
        -I ${vcf.join(" -I ")} \\
        -O ${meta.patientID}_mutect2_raw.vcf.gz

    gatk MergeMutectStats \\
        --tmp-dir ${params.tmpDir} \\
        --stats ${stats.join(" --stats ")} \\
        -O ${meta.patientID}_mutect2_raw.vcf.gz.stats

    gatk  --java-options ${java_opt} LearnReadOrientationModel \\
        --tmp-dir ${params.tmpDir} \\
        -I ${f1r2_tar_gz.join(" -I ")} \\
        -O ${meta.patientID}_read-orientation-model.tar.gz
    """
}

process 'GetPileup' {

    label 'gatk'
    label 'med_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(bam), path(bai), path(interval_list),
        path(gnomAD), path(gnomADIdx)

    output:
    tuple val(meta), path(pileup_table)

    script:
    def procSampleName = meta.patientID + "_" + meta.sampleType
    pileup_table = "${procSampleName}_pileup.table"
    """
    mkdir -p ${params.tmpDir}
    gatk GetPileupSummaries \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -O ${pileup_table} \\
        -L ${interval_list} \\
        --variant ${gnomAD}
    """
}

process 'FilterMutect2Calls' {
    
    label 'gatk'
    label 'medhigh_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/02_Mutect2/filtered_vcf", mode: 'copy'

    input:
    tuple val(meta), path(pileuptumor), path(pileupnormal),
        path(vcf), path(vcfstats), path(f1r2_tar_gz),
        path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path("${meta.patientID}_mutect2_final.{vcf.gz,vcf.gz.tbi}")

    script:
    """
    mkdir -p ${params.tmpDir}

    gatk CalculateContamination \\
        --tmp-dir ${params.tmpDir} \\
        -I ${pileuptumor} \\
        --matched-normal ${pileupnormal} \\
        -O ${meta.patientID}_cont.table && \\
    gatk FilterMutectCalls \\
        --tmp-dir ${params.tmpDir} \\
        -R ${refFasta} \\
        -V ${vcf[0]} \\
        --contamination-table ${meta.patientID}_cont.table \\
        --ob-priors ${f1r2_tar_gz} \\
        -O ${meta.patientID}_oncefiltered.vcf.gz && \\
    gatk SelectVariants \\
        --tmp-dir ${params.tmpDir} \\
        --variant ${meta.patientID}_oncefiltered.vcf.gz \\
        -R ${refFasta} \\
        --exclude-filtered true \\
        --select 'vc.getGenotype(\"${meta.patientID}_tumor\").getAD().1 >= ${params.minAD}' \\
        --output ${meta.patientID}_mutect2_final.vcf.gz
    """
}

// Strelka2 and Manta //

// Prepare .bed.gz and .bed.gz.tbi file for strelka2
process 'IntervalListtoBed' {

    label 'gatk'

    tag 'IntervalListtoBed'

    input:
    path(padded_interval_list)

    output:
    path("${padded_interval_list.baseName}.bed")

    script:
    def java_opt = '-Xmx' + task.memory.toGiga() + "G"
    """
    gatk --java-options ${java_opt} IntervalListToBed \\
        -I ${padded_interval_list} \\
        -O ${padded_interval_list.baseName}.bed
    """
}

process 'IdxZipBedFile' {

    tag 'IdxZipBedFile'

    input:
    path(bed_file)

    output:
    path("${bed_file.baseName}.{bed.gz,bed.gz.tbi}")

    script:
    """
    bgzip -c ${bed_file} > ${bed_file.baseName}.bed.gz &&
    tabix -p bed ${bed_file.baseName}.bed.gz
    """
}

// Run Manta for Indel candidating before running Strelka2
process 'Manta' {

    label 'manta'
    label 'med_memory'

    tag "$meta.patientID"
    publishDir "${params.outputDir}/${meta.patientID}/03_Strelka2/Manta", mode: 'copy'

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai),
        path(regions_bed_gz), path(regions_bed_gz_tbi),
        path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path(indel_vcf)

    script:
    indel_vcf = meta.patientID + "_candidateSmallIndels.{vcf.gz,vcf.gz.tbi}"
    """
    configManta.py \\
            --tumorBam ${tumorbam} \\
            --normalBam ${normalbam} \\
            --referenceFasta ${refFasta} \\
            --runDir manta_${meta.patientID} \\
            --callRegions ${regions_bed_gz} \\
            --exome
    manta_${meta.patientID}/runWorkflow.py -m local -j ${task.cpus}
    cp manta_${meta.patientID}/results/variants/candidateSmallIndels.vcf.gz ${meta.patientID}_candidateSmallIndels.vcf.gz
    cp manta_${meta.patientID}/results/variants/candidateSmallIndels.vcf.gz.tbi ${meta.patientID}_candidateSmallIndels.vcf.gz.tbi
    """
}

process 'Strelka2' {
    
    label 'strelka'
    label 'med_memory'

    tag "$meta.patientID"
    publishDir "${params.outputDir}/${meta.patientID}/03_Strelka2/Strelka_raw", mode: 'copy'

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai),
        path(manta_indel), path(regions_bed_gz), path(regions_bed_gz_tbi),
        path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path(somaticSNVs), path(somaticIndels), path(runstats_tsv), path(runstats_xml)

    script:
    somaticSNVs = "${meta.patientID}_somatic.snvs.{vcf.gz,vcf.gz.tbi}"
    somaticIndels = "${meta.patientID}_somatic.indels.{vcf.gz,vcf.gz.tbi}"
    runstats_tsv = "${meta.patientID}_runStats.tsv"
    runstats_xml = "${meta.patientID}_runStats.xml"
    """
    configureStrelkaSomaticWorkflow.py \\
        --tumorBam ${tumorbam} \\
        --normalBam ${normalbam} \\
        --referenceFasta ${refFasta} \\
        --indelCandidates ${manta_indel[0]} \\
        --runDir strelka_${meta.patientID} \\
        --callRegions ${regions_bed_gz} \\
        --exome
    strelka_${meta.patientID}/runWorkflow.py -m local -j ${task.cpus}
    cp strelka_${meta.patientID}/results/variants/somatic.indels.vcf.gz ${meta.patientID}_somatic.indels.vcf.gz
    cp strelka_${meta.patientID}/results/variants/somatic.indels.vcf.gz.tbi ${meta.patientID}_somatic.indels.vcf.gz.tbi
    cp strelka_${meta.patientID}/results/variants/somatic.snvs.vcf.gz ${meta.patientID}_somatic.snvs.vcf.gz
    cp strelka_${meta.patientID}/results/variants/somatic.snvs.vcf.gz.tbi ${meta.patientID}_somatic.snvs.vcf.gz.tbi
    cp strelka_${meta.patientID}/results/stats/runStats.tsv ${meta.patientID}_runStats.tsv
    cp strelka_${meta.patientID}/results/stats/runStats.xml ${meta.patientID}_runStats.xml

    """
}

process 'Strelka2CombineVCFs' {

    label 'gatk'
    label 'med_memory'

    tag "$meta.patientID"

    input:
    tuple val(meta), path(somatic_snvs), path(somatic_indels),
        path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path(combined_sorted_vcf)

    script:
    combined_sorted_vcf = "${meta.patientID}_strelka2_combined_sorted.vcf.gz"
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opt = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    mkdir -p ${params.tmpDir}

    gatk --java-options ${java_opt} MergeVcfs \\
        --TMP_DIR ${params.tmpDir} \\
        -I ${somatic_snvs[0]} \\
        -I ${somatic_indels[0]} \\
        -O ${meta.patientID}_strelka_combined.vcf.gz \\
        --SEQUENCE_DICTIONARY ${refDict} && \\
    gatk --java-options ${java_opt} SortVcf \\
        --TMP_DIR ${params.tmpDir} \\
        -I ${meta.patientID}_strelka_combined.vcf.gz \\
        -O ${combined_sorted_vcf} \\
        --SEQUENCE_DICTIONARY ${refDict}
    """
}

process 'Strelka2Reheader' {
    
    label 'low_memory'

    tag "$meta.patientID"

    input:
    tuple val(meta), path(combined_sorted_vcf)

    output:
    tuple val(meta), path(reheadered_vcf)

    script:
    reheadered_vcf = "${meta.patientID}_strelka_combined_sorted_reheadered.{vcf.gz,vcf.gz.tbi}"
    """
    printf "TUMOR ${meta.patientID}_tumor\nNORMAL ${meta.patientID}_normal\n" > vcf_rename_${meta.patientID}_tmp

    bcftools reheader \\
        -s vcf_rename_${meta.patientID}_tmp \\
        ${combined_sorted_vcf} \\
        > ${meta.patientID}_strelka_combined_sorted_reheadered.vcf.gz
    tabix -p vcf ${meta.patientID}_strelka_combined_sorted_reheadered.vcf.gz
    rm -f vcf_rename_${meta.patientID}_tmp
    """
}

process 'FinalizeStrelka2' {

    label 'gatk'
    label 'med_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/03_Strelka2/Strelka_filtered", mode: 'copy'

    input:
    tuple val(meta), path(reheadered_vcf), path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path(final_vcf)

    script:
    final_vcf = "${meta.patientID}_strelka2_somatic_filtered.{vcf.gz,vcf.gz.tbi}"
    """
    mkdir -p ${params.tmpDir}

    gatk SelectVariants \\
        --tmp-dir ${params.tmpDir} \\
        --variant ${reheadered_vcf[0]} \\
        -R ${refFasta} \\
        --exclude-filtered true \\
        --output ${meta.patientID}_strelka2_somatic_filtered.vcf.gz
    """
}

// VarScan2 //
// Bed list prepration
process 'scatteredIntervalListToBed' {
    
    label 'gatk'
    label 'low_memory'

    tag 'scatteredIntervalListToBed'

    input:
    path(intervals)

    output:
    path("${intervals.baseName}.bed")

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    gatk --java-options ${JAVA_Xmx} IntervalListToBed \\
        -I ${intervals} \\
        -O ${intervals.baseName}.bed
    """
}

process 'SamToolsMPileup' {

    label 'med_memory'

    tag "$meta.patientID"

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai),
        path(intervals), path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path("${meta.patientID}_${intervals.baseName}.mpileup")

    script:
    """
    samtools mpileup \\
        -q 1 \\
        -f ${refFasta} \\
        -l ${intervals} \\
        ${normalbam} ${tumorbam} > ${meta.patientID}_${intervals.baseName}.mpileup
    """
}

process 'Varscan2' {

    label 'varscan'
    label 'med_memory'

    tag "$meta.patientID"

    input:
    tuple val(meta), path(mpileup)

    output:
    tuple val(meta), path("${meta.patientID}_${mpileup.baseName}_varscan.snp.vcf"),
        path("${meta.patientID}_${mpileup.baseName}_varscan.indel.vcf")

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    varscan ${JAVA_Xmx} somatic \\
        ${mpileup} \\
        ${meta.patientID}_${mpileup.baseName}_varscan_tmp \\
        --output-vcf 1 \\
        --mpileup 1 \\
        --min-coverage ${params.min_cov} \\
        --min-coverage-normal ${params.min_cov_normal} \\
        --min-coverage-tumor ${params.min_cov_tumor} \\
        --p-value ${params.pval_call_hetero} \\
        --somatic-p-value ${params.somatic_pvalue} \\
        --strand-filter ${params.strand_filter} && \\
    awk '{OFS=FS="\t"} { if(\$0 !~ /^#/) { if (\$4 ~ /[ACGT]/) { print } } else { print } }' \\
        ${meta.patientID}_${mpileup.baseName}_varscan_tmp.snp.vcf \\
        > ${meta.patientID}_${mpileup.baseName}_varscan.snp.vcf && \\
    awk '{OFS=FS="\t"} { if(\$0 !~ /^#/) { if (\$4 ~ /[ACGT]+/) { print } } else { print } }' \\
        ${meta.patientID}_${mpileup.baseName}_varscan_tmp.indel.vcf \\
        > ${meta.patientID}_${mpileup.baseName}_varscan.indel.vcf
    
    rm -f ${meta.patientID}_${mpileup.baseName}_varscan_tmp.*
    """
}

process 'GatherVarscan2scatters' {
    
    label 'gatk'
    label 'med_memory'

    tag "$meta.patientID"

    input:
    tuple val(meta), path(snp_vcf), path(indel_vcf),
        path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path("${meta.patientID}_varscan.snp.vcf"),
                     path("${meta.patientID}_varscan.indel.vcf")

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opt = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    mkdir -p ${params.tmpDir}

    gatk --java-options ${java_opt} MergeVcfs \\
        --TMP_DIR ${params.tmpDir} \\
        -I ${snp_vcf.join(" -I ")} \\
        -O ${meta.patientID}_varscan.snp.vcf \\
        --SEQUENCE_DICTIONARY ${refDict}

    gatk --java-options ${java_opt} MergeVcfs \\
        --TMP_DIR ${params.tmpDir} \\
        -I ${indel_vcf.join(" -I ")} \\
        -O ${meta.patientID}_varscan.indel.vcf \\
        --SEQUENCE_DICTIONARY ${refDict}
    """
}

process 'ProcessSomaticVarscan2Calls' {
    
    label 'varscan'
    label 'med_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/04_Varscan2/raw", mode: 'copy'

    input:
    tuple val(meta), path(snp), path(indel)

    output:
    tuple val(meta), path(somatic_snp), path(somatic_hc_snp),
            path(somatic_indel), path(somatic_hc_indel)

    script:
    somatic_snp      = "${meta.patientID}_varscan.snp.Somatic.vcf"
    somatic_hc_snp   = "${meta.patientID}_varscan.snp.Somatic.hc.vcf"
    somatic_indel    = "${meta.patientID}_varscan.indel.Somatic.vcf"
    somatic_hc_indel = "${meta.patientID}_varscan.indel.Somatic.hc.vcf"
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    varscan ${JAVA_Xmx} processSomatic \\
        ${snp} \\
        --min-tumor-freq ${params.min_tumor_freq} \\
        --max-normal-freq ${params.max_normal_freq} \\
        --p-value ${params.processSomatic_pvalue} && \\
    varscan ${JAVA_Xmx} processSomatic \\
        ${indel} \\
        --min-tumor-freq ${params.min_tumor_freq} \\
        --max-normal-freq ${params.max_normal_freq} \\
        --p-value ${params.processSomatic_pvalue}
    """
}

process 'BamReadCount' {

    label 'bam_readcount'
    label 'med_memory'

    tag "$meta.patientID"

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), 
            path(somatic_hc_snp), path(somatic_hc_indel),
            path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path(readcount_snp), path(readcount_indel)

    script:
    readcount_snp = "${meta.patientID}_varscan_hc_snp_readcount.tsv"
    readcount_indel = "${meta.patientID}_varscan_hc_indel_readcount.tsv"
    """
    cat ${somatic_hc_snp} | \\
    awk '{if (!/^#/) { x = length(\$5) - 1; print \$1,\$2,(\$2+x); }}' | \\
    bam-readcount \\
        -q${params.min_map_q} \\
        -b${params.min_base_q} \\
        -w1 \\
        -l /dev/stdin \\
        -f ${refFasta} \\
        ${tumorbam} > ${readcount_snp} && \\
    cat ${somatic_hc_indel} | \\
    awk '{if (!/^#/) { x = length(\$5) - 1; print \$1,\$2,(\$2+x); }}' | \\
    bam-readcount \\
        -q${params.min_map_q} \\
        -b${params.min_base_q} \\
        -w1 \\
        -l /dev/stdin \\
        -f ${refFasta} \\
        ${tumorbam} > ${readcount_indel}
    """
}

process 'FilterVarscan2Calls' {

    label 'varscan'
    label 'med_memory'

    tag "$meta.patientID"

    input:
    tuple val(meta), path(somatic_hc_snp), path(somatic_hc_indel),
            path(readcount_snp), path(readcount_indel)

    output:
    tuple val(meta), path(filtered_hc_snp), path(filtered_hc_indel)

    script:
    filtered_hc_snp =  "${meta.patientID}_varscan.snp.Somatic.hc.filtered.vcf"
    filtered_hc_indel = "${meta.patientID}_varscan.indel.Somatic.hc.filtered.vcf"
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    varscan ${JAVA_Xmx} fpfilter \\
        ${somatic_hc_snp} \\
        ${readcount_snp} \\
        --output-file ${filtered_hc_snp} && \\
    varscan ${JAVA_Xmx} fpfilter \\
        ${somatic_hc_indel} \\
        ${readcount_indel} \\
        --output-file ${filtered_hc_indel}
    """
}

process 'BgzipFilteredVarscanCalls' {

    tag "$meta.patientID"

    input:
    tuple val(meta), path(filtered_hc_snp), path(filtered_hc_indel)

    output:
    tuple val(meta), path(filtered_hc_snp_gz), path(filtered_hc_indel_gz)

    script:
    filtered_hc_snp_gz = "${meta.patientID}_varscan.snp.Somatic.hc.filtered.{vcf.gz,vcf.gz.tbi}"
    filtered_hc_indel_gz = "${meta.patientID}_varscan.indel.Somatic.hc.filtered.{vcf.gz,vcf.gz.tbi}"
    """
    bgzip -c ${filtered_hc_snp} > ${meta.patientID}_varscan.snp.Somatic.hc.filtered.vcf.gz
    tabix -p vcf ${meta.patientID}_varscan.snp.Somatic.hc.filtered.vcf.gz
    bgzip -c ${filtered_hc_indel} > ${meta.patientID}_varscan.indel.Somatic.hc.filtered.vcf.gz
    tabix -p vcf ${meta.patientID}_varscan.indel.Somatic.hc.filtered.vcf.gz
    """    
}

process 'CombineSortVarscan2FilteredCalls' {
    
    label 'gatk'
    label 'med_memory'

    tag "$meta.patientID"

    input:
    tuple val(meta), path(filtered_hc_snp_gz), path(filtered_hc_indel_gz),
            path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path(combined_sorted_vcf)

    script:
    combined_sorted_vcf = "${meta.patientID}_varscan_combined_sorted.vcf.gz"
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opt = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    mkdir -p ${params.tmpDir}

    gatk --java-options ${java_opt} MergeVcfs \\
        --TMP_DIR ${params.tmpDir} \\
        -I ${filtered_hc_snp_gz[0]} \\
        -I ${filtered_hc_indel_gz[0]} \\
        -O ${meta.patientID}_varscan_combined.vcf.gz \\
        --SEQUENCE_DICTIONARY ${refDict} && \\
    gatk --java-options ${java_opt} SortVcf \\
        --TMP_DIR ${params.tmpDir} \\
        -I ${meta.patientID}_varscan_combined.vcf.gz \\
        -O ${combined_sorted_vcf} \\
        --SEQUENCE_DICTIONARY ${refDict}
    """
}

process 'FinalizeVarscan2Calls' {

    label 'low_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/04_Varscan2/filtered", mode: 'copy'

    input:
    tuple val(meta), path(combined_sorted_vcf)

    output:
    tuple val(meta), path(reheadered_vcf)

    script:
    reheadered_vcf = "${meta.patientID}_varscan_somatic_filtered.{vcf.gz,vcf.gz.tbi}"
    """
    printf "TUMOR ${meta.patientID}_tumor\nNORMAL ${meta.patientID}_normal\n" > vcf_rename_${meta.patientID}_tmp

    bcftools reheader \\
        -s vcf_rename_${meta.patientID}_tmp \\
        ${combined_sorted_vcf} \\
        > ${meta.patientID}_varscan_somatic_filtered.vcf.gz
    tabix -p vcf ${meta.patientID}_varscan_somatic_filtered.vcf.gz
    rm -f vcf_rename_${meta.patientID}_tmp
    """
}

//VarDict prepration and running
process 'scatterVardict' {

    label 'vardict'
    label 'medhigh_memory'

    tag "$meta.patientID"

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai),
            path(intervals), path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path(scattered_vcf)

    script:
    scattered_vcf = "${meta.patientID}_${intervals.baseName}_vardict_somatic.vcf"
    """
    vardict-java \\
        -th ${task.cpus} \\
        -G ${refFasta} \\
        -f ${params.AF_THR} \\
        -N "${meta.patientID}_tumor" \\
        -b "${tumorbam}|${normalbam}" \\
        -c 1 -S 2 -E 3 -g 4 \\
        ${intervals} | \\
        testsomatic.R | \\
        var2vcf_paired.pl \\
        -N "${meta.patientID}_tumor|${meta.patientID}_normal" \\
        -f ${params.AF_THR} \\
        > ${meta.patientID}_${intervals.baseName}_tmp.vcf && \\
    awk '{OFS=FS="\t"} { if(\$0 !~ /^#/) { if (\$4 ~ /[ACGT]/) { print } } else { print } }' \\
        ${meta.patientID}_${intervals.baseName}_tmp.vcf \\
        > ${scattered_vcf}
    """
}

process 'GatherVardictCalls' {
    
    label 'gatk'
    label 'med_memory'

    tag "$meta.patientID"

    input:
    tuple val(meta), path(vcf),
        path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path("${meta.patientID}_vardict_somatic.vcf")

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opt = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    mkdir -p ${params.tmpDir}

    gatk --java-options ${java_opt} MergeVcfs \\
        --TMP_DIR ${params.tmpDir} \\
        -I ${vcf.join(" -I ")} \\
        -O ${meta.patientID}_vardict_somatic.vcf \\
        --SEQUENCE_DICTIONARY ${refDict}
    """
}

process 'FinalizeVardictCalls' {

    label 'med_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/05_Vardict/filtered", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path(filtered_vcf)

    script:
    filtered_vcf = "${meta.patientID}_vardict_somatic_filtered.{vcf.gz,vcf.gz.tbi}"
    """
    bcftools view -i 'FILTER="PASS" && INFO/STATUS="StrongSomatic"' \\
        ${vcf} \\
        > ${meta.patientID}_vardict_tmp.vcf
    bgzip -c ${meta.patientID}_vardict_tmp.vcf > ${meta.patientID}_vardict_somatic_filtered.vcf.gz
    tabix -p vcf ${meta.patientID}_vardict_somatic_filtered.vcf.gz
    """
}

// Intersect Variants
process 'VcfIntersection' {
    
    label 'med_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/06_VariantCaller_isec", mode: 'copy'

    input:
    tuple val(meta), path(mutect2), path(strelka), path(varscan), path(vardict)

    output:
    tuple val(meta), path("${meta.patientID}_intersected_somatic.vcf.gz")

    script:
    def num_caller = "-n +" + "${params.caller_to_isec}"
    """
    vcf-isec \\
        ${num_caller} \\
        ${mutect2[0]} \\
        ${strelka[0]} \\
        ${varscan[0]} \\
        ${vardict[0]} \\
        > "${meta.patientID}_intersected_somatic.vcf.gz"
    """
}

// MSIsensor Pro
process 'MakeMSISensorRef' {

    label 'msiSensorPro'
    label 'medhigh_memory'

    tag "MakeMSISensorRef"

    input:
    tuple path(refFasta), path(refIdx), path(refDict)

    output:
    path("microsatellites.txt")

    script:
    """
    msisensor-pro \\
        scan \\
        -d ${refFasta} \\
        -o microsatellites.txt
    """
}

process 'runMSISensorPro' {

    label 'msiSensorPro'
    label 'medhigh_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/07_MSISensorPro", mode: 'copy'

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai), path(msRef)

    output:
    tuple val(meta), path("${meta.patientID}_msi_dis"),
                     path("${meta.patientID}_msi_unstable"),
                     path("${meta.patientID}_msi_all"),
                     path("${meta.patientID}_msi")

    script:
    """
    mkdir -p results
    msisensor-pro \\
        msi \\
        -d ${msRef} \\
        -n ${normalbam} \\
        -t ${tumorbam} \\
        -o "${meta.patientID}_msi" \\
        -f 0.05 \\
        -c ${params.cov_msi} \\
        -b ${task.cpus}
    """
}


// CNV ANALYSIS
params.runCNV = false

// Prepare .bed file
process 'IntervalListtoBed_CNV' {

    label 'gatk'

    tag 'IntervalListtoBed_CNV'

    input:
    path(padded_interval_list)

    output:
    path("${padded_interval_list.baseName}.bed")

    script:
    def java_opt = '-Xmx' + task.memory.toGiga() + "G"
    """
    gatk --java-options ${java_opt} IntervalListToBed \\
        -I ${padded_interval_list} \\
        -O ${padded_interval_list.baseName}.bed
    """
}

// CNV_FACETS
process 'cnv_facets' {

    label 'facets'
    label 'medhigh_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/08_CNVs/FACETS", mode: 'copy'

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai),
            path(interval), path(common_vcf), path(common_vcf_idx)

    output:
    tuple val(meta), path("${meta.patientID}_facets.cnv.png"), path("${meta.patientID}_facets.cov.pdf"),
                        path("${meta.patientID}_facets.csv.gz"), path("${meta.patientID}_facets.spider.pdf"),
                        path(cnv_out)
    
    script:
    cnv_out = "${meta.patientID}_facets.{vcf.gz,vcf.gz.tbi}"
    """
    cnv_facets.R \\
        -o ${meta.patientID}_facets \\
        -t ${tumorbam} \\
        -n ${normalbam} \\
        -vcf ${common_vcf} \\
        --snp-mapq ${params.snp_mapq} \\
        --snp-baq ${params.snp_baq} \\
        -N ${task.cpus} \\
        --depth ${params.facets_depth} \\
        --target ${interval} \\
        --cval ${params.facets_cval} \\
        --nbhd-snp ${params.nbhd_snp} \\
        --annotation ${interval} \\
        --gbuild ${params.facets_gbuild}
    """
}

process 'SplitIntervalByChr' {

    label 'med_memory'

    tag 'SplitIntervalByChr'

    input:
    path(snps)

    output:
    path("*.txt")

    script:
    """
    awk '{print > \$1".txt"}' ${snps}
    mv X.txt 23.txt
    """
}

process 'AlleleCount' {

    label 'AlleleCounter'
    label 'med_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(bam), path(bai), path(loci), path(refFasta), path(refIdx), path(refDict)

    output:
    tuple val(meta), path(allelecountFile)

    script:
    def procSampleName = loci.baseName + "_" + meta.patientID + "_" + meta.sampleType
    allelecountFile = "${procSampleName}.txt"
    """
    alleleCounter \\
        -l ${loci} \\
        -b ${bam} \\
        -o ${allelecountFile} \\
        -r ${refIdx} \\
        --min-base-qual ${params.cnv_min_base_q} \\
        --min-map-qual ${params.cnv_min_map_q} \\
        -d
    """
}

process 'GatherAlleleCounts' {

    label 'medhigh_memory'
    label 'Dplyr'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(count_files)

    output:
    tuple val(meta), path("${meta.patientID}_${meta.sampleType}_allelecount.txt")

    script:
    """
    Rscript ${baseDir}/bin/merge_allelecounter.R \\
        "${meta.patientID}_${meta.sampleType}_allelecount.txt" \\
        ${count_files.join(' ')}
    """
}

process 'FindCountsInsideInterval' {

    label 'med_memory'
    label 'ascat'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(gathered_counts), path(interval)

    output:
    tuple val(meta), path("${meta.patientID}_${meta.sampleType}_filtered_allelecount.txt")

    script:
    """
    Rscript ${baseDir}/bin/FindCountsInsideInterval.R \\
        ${gathered_counts} \\
        ${interval} \\
        "${meta.patientID}_${meta.sampleType}_filtered_allelecount.txt"
    """
}

process 'FinalizeCounts' {
    
    label 'Dplyr'
    label 'med_memory'

    tag "$meta.patientID : $meta.sampleType"

    input:
    tuple val(meta), path(filtered_counts)

    output:
    tuple val(meta), path("${meta.patientID}_${meta.sampleType}_dedpu_filtered_allelecount.txt")

    script:
    """
    Rscript ${baseDir}/bin/FinalizeCounts.R \\
        ${filtered_counts} \\
        "${meta.patientID}_${meta.sampleType}_dedpu_filtered_allelecount.txt"
    """
}

process 'ConvertAlleleCountToLogBaf' {

    label 'ascat'
    label 'medhigh_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/08_CNVs/ASCAT/LogBaf_files", mode: 'copy'

    input:
    tuple val(meta), path(tumorcount), path(normalcount)

    output:
    tuple val(meta), path("${meta.patientID}_tumor.LogR"),
                     path("${meta.patientID}_tumor.BAF"),
                     path("${meta.patientID}_normal.LogR"),
                     path("${meta.patientID}_normal.BAF")

    script:
    """
    Rscript ${baseDir}/bin/convertAlleleCounts.r \\
        "${meta.patientID}_tumor" \\
        "${tumorcount}" \\
        "${meta.patientID}_normal" \\
        "${normalcount}" \\
        "XX"
    """
}


process 'ASCAT' {

    label 'ascat'
    label 'medhigh_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/08_CNVs/ASCAT", mode: 'copy'

    input:
    tuple val(meta), path(tumorlog), path(tumorbaf), path(normallog), path(normalbaf),
          path(gccontent), path(reptiming)

    output:
    tuple val(meta), path("*.png"),
                     path("${meta.patientID}_tumor.BAF.PCFed.txt"),
                     path("${meta.patientID}_tumor.LogR.PCFed.txt"),
                     path("ascat_object.rds"),
                     path("ASCAT_objects.Rdata"),
                     path("ASCAT_QC.txt")

    script:
    """
    Rscript ${baseDir}/bin/run_ASCAT.R \\
        "${tumorlog}" \\
        "${tumorbaf}" \\
        "${normallog}" \\
        "${normalbaf}" \\
        "XX" \\
        "${gccontent}" \\
        "${reptiming}" \\
        "${params.cnv_seg_penalty}"
    """
}

process 'EstimatePloidyPurity' {

    label 'VariantAnnotation'
    label 'low_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/08_CNVs", mode: 'copy'

    input:
    tuple val(meta), path(facets_file), path(ascat_file)

    output:
    tuple val(meta), path("${meta.patientID}_Ploidy_Purity.txt"),
                     path("${meta.patientID}_copynumber.txt")

    script:
    """
    Rscript ${baseDir}/bin/ExtractPloidyPurity.R \\
        "${facets_file[0]}" \\
        "${ascat_file}" \\
        "${meta.patientID}_Ploidy_Purity.txt" \\
        "${meta.patientID}_copynumber.txt" \\
        "${meta.patientID}_tumor_recal_sorted"
    """
}

process 'RunCNVkit' {

    label 'CNVkit'
    label 'high_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/08_CNVs/CNVkit", mode: 'copy'

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai),
          val(purity), val(ploidy), path(refFasta), path(refIdx), path(baitsbed)

    output:
    tuple val(meta), path("${tumorbam.baseName}.bintest.cns"),
                     path("${tumorbam.baseName}_breaks.tsv"),
                     path("${tumorbam.baseName}.call.cns"),
                     path("${tumorbam.baseName}.cnr"),
                     path("${tumorbam.baseName}.cns"),
                     path("${tumorbam.baseName}_diagram.pdf"),
                     path("${tumorbam.baseName}_gainloss.tsv"),
                     path("${tumorbam.baseName}_scatter.png"),
                     path("${tumorbam.baseName}.segmetrics.cns"),
                     path("output_reference.cnn")
    
    script:
    """
    cnvkit.py \\
        batch \\
        ${tumorbam} \\
        -m hybrid \\
        -p ${task.cpus} \\
        --normal ${normalbam} \\
        --fasta ${refFasta} \\
        --targets ${baitsbed} \\
        --short-names \\
        --output-reference output_reference.cnn \\
        --output-dir ./
    
    cnvkit.py segmetrics \\
        -s ${tumorbam.baseName}.cn{s,r} \\
        --ci \\
        --pi
    
    cnvkit.py call \\
        ${tumorbam.baseName}.cns \\
        --filter ci \\
        -m clonal \\
        --purity ${purity} \\
        --ploidy ${ploidy} \\
        --sample-sex "female" \\
        -o ${tumorbam.baseName}.call.cns
    
    cnvkit.py \\
        scatter \\
        ${tumorbam.baseName}.cnr \\
        -s ${tumorbam.baseName}.cns \\
        -o ${tumorbam.baseName}_scatter.png

    cnvkit.py \\
        diagram \\
        ${tumorbam.baseName}.cnr \\
        -s ${tumorbam.baseName}.cns \\
        --sample-sex "female" \\
        -o ${tumorbam.baseName}_diagram.pdf

    cnvkit.py \\
        breaks \\
        ${tumorbam.baseName}.cnr ${tumorbam.baseName}.cns \\
        -o ${tumorbam.baseName}_breaks.tsv

    cnvkit.py \\
        genemetrics \\
        ${tumorbam.baseName}.cnr \\
        -s ${tumorbam.baseName}.cns \\
        --sample-sex "female" \\
        -t 0.2 -m 5 \\
        -o ${tumorbam.baseName}_gainloss.tsv
    """
}

process 'Polysolver' {

    label 'Polysolver'
    label 'medhigh_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/09_LOH_HLA/HLA_typing", mode: 'copy'

    input:
    tuple val(meta), path(tumorbam), path(tumorbai)

    output:
    tuple val(meta), path("${meta.patientID}_hla.txt")

    script:
    """
    shell_call_hla_type \\
        ${tumorbam} \\
        Unknown \\
        1 \\
        hg38 \\
        STDFQ \\
        0 \\
        .
    
    awk '{for (i=2; i<=NF; i++) print \$i}' \\
        winners.hla.txt \\
        > ${meta.patientID}_hla.txt
    """
}

process 'RunLOHHLA' {

    label 'LOHHLA'
    label 'high_memory'
    
    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/09_LOH_HLA/LOH_HLA_outputs", mode: 'copy'

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai),
          path(patientHLA), path(hlaFasta), path(hlaIdx), path(copynumber)

    output:
    tuple val(meta), path("*")

    script:
    """
    mkdir -p bam_dir
    mv ${tumorbam} ${tumorbai} ${normalbam} ${normalbai} ./bam_dir

    lohhla \\
        --patientId ${meta.patientID} \\
        --outputDir ./ \\
        --normalBAMfile ./bam_dir/${normalbam} \\
        --BAMDir ./bam_dir/ \\
        --hlaPath ${patientHLA} \\
        --HLAfastaLoc ${hlaFasta} \\
        --CopyNumLoc ${copynumber} \\
        --minCoverageFilter ${params.min_cov_filter} \\
        --gatkDir /picard-tools \\
        --novoDir /opt/conda/bin \\
        --numMisMatch ${params.num_missMatch} \\
        --mappingStep TRUE \\
        --fishingStep TRUE \\
        --plottingStep TRUE \\
        --coverageStep TRUE \\
        --cleanUp TRUE
    """
}

process 'TcellExTRECT' {

    label 'TcellExTRECT'
    label 'medhigh_memory'

    tag "$meta.patientID"

    publishDir "${params.outputDir}/${meta.patientID}/10_TcellExTRECT", mode: 'copy'

    input:
    tuple val(meta), path(tumorbam), path(tumorbai),
          path(cnvkitseg), val(purity)

    output:
    tuple val(meta), path("${meta.patientID}_adj_TcellExTRECT.txt")

    script:
    """
    Rscript ${baseDir}/bin/runTcellExTRECT.R \\
        "${tumorbam}" \\
        "${meta.patientID}_tumor" \\
        "${params.median_th}" \\
        "${purity}" \\
        "${cnvkitseg}" \\
        "${meta.patientID}_adj_TcellExTRECT.txt" \\
        "./"
    """
}

// Define the workflow
workflow {
    log.info "Starting workflow"

    // Transform fastq_ch.multiple for input to merge_fastq
    def merge_fastq_input = fastq_ch.multiple.map { meta, nested_lists ->
        // Flatten the nested list of R1 and R2 files
        def r1_files = nested_lists.collect { it[0] }.flatten()
        def r2_files = nested_lists.collect { it[1] }.flatten()
        // Return the tuple in the correct format
        tuple(meta, r1_files, r2_files)
    }

    // Directly run the merge process
    def merged_fastq_ch = merge_fastq(merge_fastq_input)

    // make single channel ready to be merged to multiple
    def single_fastq_ch = fastq_ch.single.map { meta, fastq_list ->
        def r1 = fastq_list[0][0]
        def r2 = fastq_list[0][1]
        
        tuple(meta, r1, r2)
        }
    // Concatenate `merged_fastq_ch` and `single_fastq_ch`
    def all_fastq_ch = merged_fastq_ch.mix(single_fastq_ch)

    // Run cut_adaptor process
    // def trimmed_fastq_ch = cut_adaptors(all_fastq_ch)

    // Make uBAM
    def MakeUbamOut_ch = MakeUbam(all_fastq_ch)
    
    // Run BwaAlign processs
    def bam_file_ch = BwaAlign(all_fastq_ch.map { meta, r1, r2 ->
        tuple(meta, r1, r2, file(params.RefFasta), file(params.RefIdx),
             file(params.RefDict), file(params.BwaInd))
        })
    
    // Merge BAM and uBAM
    def MergeBamsOut_ch = MergeBams(bam_file_ch
    .join(MakeUbamOut_ch, by: [0])
    .map { meta, bam, ubam -> 
           tuple(meta, bam, ubam, file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})

    // Run MarkDuplicates
    def mkdp_unsorted_bam_ch = MarkDuplicates(MergeBamsOut_ch)

    // Sort MarkDuplicates
    def mkdp_sorted_bam_ch = SamtoolsSortDups(mkdp_unsorted_bam_ch)

    // Tag sorted MKDPs
    def final_preprocess_bam_ch = SetNmMdAndUqTags(mkdp_sorted_bam_ch.map { meta, mkdp_sorted_bam ->
        tuple(meta, mkdp_sorted_bam, file(params.RefFasta), file(params.RefIdx), file(params.RefDict))
        })

    // Convert regions and baits .bed to interval list for Mutect2 and other processes    
    def RegionsBedToIntervalList_ch = RegionsBedToIntervalList(tuple(file(params.RefDict), file(params.RegionsBed)))
    def BaitsBedToIntervalList_ch = BaitsBedToIntervalList(tuple(file(params.RefDict), file(params.BaitsBed)))

    // Preprocess regions interval list
    def preprocessIntervalList_ch = preprocessIntervalList(RegionsBedToIntervalList_ch.map { interval_list ->
        tuple(interval_list, file(params.RefFasta), file(params.RefIdx), file(params.RefDict))
        })

    // Run splitting for further scattering processes
    def SplitIntervals_ch = SplitIntervals(preprocessIntervalList_ch.map { outFileName ->
        tuple(outFileName, file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})

    // Run Base recalibration
    def scatterBaseRecal_ch = scatterBaseRecalibrator(
        final_preprocess_bam_ch.combine(SplitIntervals_ch.flatten())
        .map { meta, bam_out, intervals ->
            tuple(meta, bam_out, intervals, 
                    file(params.RefFasta), file(params.RefIdx), file(params.RefDict),
                    file(params.DBSNP), file(params.DBSNPIdx),
                    file(params.MillsGold), file(params.MillsGoldIdx),
                    file(params.KnownIndels), file(params.KnownIndelsIdx))})

    def gatherBQSRtables_ch = gatherscatteredBQSRtables(
        scatterBaseRecal_ch.groupTuple(by: [0])
    )

    def scatterapplyBQSRS_ch = scatterapplyBQSRS(
        final_preprocess_bam_ch
        .join(gatherBQSRtables_ch, by: [0])
        .combine(SplitIntervals_ch.flatten())
        .map { meta, bam_out, bqsr_table, intervals ->
            tuple(meta, bam_out, bqsr_table, intervals,
                    file(params.RefFasta), file(params.RefIdx), file(params.RefDict),
                    file(params.DBSNP), file(params.DBSNPIdx),
                    file(params.MillsGold), file(params.MillsGoldIdx),
                    file(params.KnownIndels), file(params.KnownIndelsIdx))})

    def gatherRecalBamFiles_ch = gatherRecalBamFiles(scatterapplyBQSRS_ch
        .toSortedList({a, b -> a[1][0].baseName <=> b[1][0].baseName})
        .flatten()
        .collate(3)
        .groupTuple(by: [0]))

    def SamtoolsSortIdxRecal_ch = SamtoolsSortIdxRecal(gatherRecalBamFiles_ch)

    // Alignment metrics
    def AlignmentMetricsOut_ch = AlignmentMetrics(SamtoolsSortIdxRecal_ch
    .combine(RegionsBedToIntervalList_ch)
    .combine(BaitsBedToIntervalList_ch)
    .map { meta, bam, bai, target, bait ->
           tuple(meta, bam, bai, target, bait, 
                 file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})


    // Branch SamtoolsSortIdxRecal_ch to tumor and normal
    def BamsForVariantCallers_ch = SamtoolsSortIdxRecal_ch
        .branch {
        meta_ori, bam, bai ->
            def meta = meta_ori.clone()
            tumor : meta.sampleType == "tumor"
                meta.remove('sampleType')
                return [meta, bam, bai]

            normal: meta.sampleType == "normal"
                meta.remove('sampleType')
                return [meta, bam, bai]
            }

    // Make Mutect2 input and run + filtering
    def Mutect2Inputs_ch = BamsForVariantCallers_ch.tumor
    .join(BamsForVariantCallers_ch.normal, by: [0])
    .combine(SplitIntervals_ch.flatten())
    .map { meta, tumorBam, tumorBai, normalBam, normalBai, intervals ->
        tuple(meta, tumorBam, tumorBai, normalBam, normalBai, intervals,
                file(params.RefFasta), file(params.RefIdx), file(params.RefDict),
                file(params.GnomAD), file(params.GnomADIdx),
                file(params.PON), file(params.PONIdx))}
    
    def Mutect2Outputs_ch = Mutect2(Mutect2Inputs_ch)

    def Mutect2Gathered_ch = GatherMutect2Calls(Mutect2Outputs_ch
    .toSortedList {a, b -> a[1].baseName <=> b[1].baseName}
    .flatten()
    .collate(5)
    .groupTuple(by: [0]))

    def GetPileupOut_ch = GetPileup(SamtoolsSortIdxRecal_ch.combine(preprocessIntervalList_ch)
    .map { meta, bam, bai, interval_list ->
        tuple(meta, bam, bai, interval_list, file(params.GnomAD), file(params.GnomADIdx))})

    def PileupForFilter_ch = GetPileupOut_ch
        .branch {
            meta_ori, pileup_table ->
            def meta = meta_ori.clone()
            tumor : meta.sampleType == "tumor"
                meta.remove('sampleType')
                return [meta, pileup_table]

            normal: meta.sampleType == "normal"
                meta.remove('sampleType')
                return [meta, pileup_table]
        }

    def FilterMutect2Input_ch = FilterMutect2Calls(PileupForFilter_ch.tumor
    .join(PileupForFilter_ch.normal, by: [0])
    .combine(Mutect2Gathered_ch, by: [0])
    .map {meta, pileuptumor, pileupnormal, vcf, vcfstats, f1r2_tar_gz ->
        tuple(meta, pileuptumor, pileupnormal, vcf, vcfstats, f1r2_tar_gz,
            file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})

    // Prepare .bed.gz and .bed.gz.tbi file for strelka2
    def IntervalListToBed_ch = IdxZipBedFile(IntervalListtoBed(preprocessIntervalList_ch))

    // Run Manta before strelka2
    def MantaCandIndel_ch = Manta(BamsForVariantCallers_ch.tumor
    .join(BamsForVariantCallers_ch.normal, by: [0])
    .combine(IntervalListToBed_ch)
    .map {meta, tumorbam, tumorbai, normalbam, normalbai, regions_bed_gz, regions_bed_gz_tbi ->
        tuple(meta, tumorbam, tumorbai, normalbam, normalbai, regions_bed_gz, regions_bed_gz_tbi,
                file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})
    
    // Run Strelka2
    def Strelka2PrimOut_ch = Strelka2(BamsForVariantCallers_ch.tumor
    .join(BamsForVariantCallers_ch.normal, by: [0])
    .combine(MantaCandIndel_ch, by: [0])
    .combine(IntervalListToBed_ch)
    .map {meta, tumorbam, tumorbai, normalbam, normalbai, indel_vcf, regions_bed_gz, regions_bed_gz_tbi ->
        tuple(meta, tumorbam, tumorbai, normalbam, normalbai, indel_vcf,
                regions_bed_gz, regions_bed_gz_tbi,
                file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})

    def Strelka2CombinedSorted_ch = Strelka2CombineVCFs(Strelka2PrimOut_ch
    .map {meta, somatic_snvs, somatic_indels, stats_tsv, stats_xml ->
        tuple(meta, somatic_snvs, somatic_indels,
            file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})
    
    def Strelka2Reheadered_ch = Strelka2Reheader(Strelka2CombinedSorted_ch)

    def Strelka2Final_ch = FinalizeStrelka2(Strelka2Reheadered_ch
    .map { meta, reheadered_vcf ->
            tuple(meta, reheadered_vcf,
            file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})

    // Varscan2 prepration and run
    def scatteredIntervalListToBed_ch = scatteredIntervalListToBed(SplitIntervals_ch.flatten())
    // samtools mpileup
    def MpileupOut_ch = SamToolsMPileup(BamsForVariantCallers_ch.tumor
    .join(BamsForVariantCallers_ch.normal, by: [0])
    .combine(scatteredIntervalListToBed_ch.flatten())
    .map {meta, tumorbam, tumorbai, normalbam, normalbai, interval ->
        tuple(meta, tumorbam, tumorbai, normalbam, normalbai, interval,
            file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})
    // Run scatter varscan
    def ScatteredVarscanOut_ch = Varscan2(MpileupOut_ch)
    // Gather Varscan2 scatters
    def GatherVarscan2Calls_ch = GatherVarscan2scatters(ScatteredVarscanOut_ch
    .toSortedList({a, b -> a[1].baseName <=> b[1].baseName})
    .flatten()
    .collate(3)
    .groupTuple(by: [0])
    .map { meta, snp_vcf, indel_vcf ->
        tuple(meta, snp_vcf, indel_vcf,
            file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})

    def ProcessedVarscan2Calls_ch = ProcessSomaticVarscan2Calls(GatherVarscan2Calls_ch)

    // bam-readcount for snp and indel
    def BamReadCountOut_ch = BamReadCount(BamsForVariantCallers_ch.tumor
    .combine(ProcessedVarscan2Calls_ch, by: [0])
    .map {meta, tumorbam, tumorbai, snp, snp_hc, indel, indel_hc ->
            tuple(meta, tumorbam, tumorbai, snp_hc, indel_hc,
                    file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})
    
    // Filter varscan2 hc snp and indels
    def FilteredVarscan2Out_ch = FilterVarscan2Calls(ProcessedVarscan2Calls_ch
    .combine(BamReadCountOut_ch, by: [0])
    .map { meta, snp, snp_hc, indel, indel_hc, readcount_snp, readcount_indel ->
            tuple(meta, snp_hc, indel_hc, readcount_snp, readcount_indel)})

    def BgzipFilteredVarscanCalls_ch = BgzipFilteredVarscanCalls(FilteredVarscan2Out_ch)

    def VarscanCombinedSorted_ch = CombineSortVarscan2FilteredCalls(BgzipFilteredVarscanCalls_ch
    .map { meta, snp, indel ->
            tuple(meta, snp, indel,
                file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})
    
    def FinalizeVarscan2Calls_ch = FinalizeVarscan2Calls(VarscanCombinedSorted_ch)

    // Vardict run and process
    def scatterVardict_ch = scatterVardict(BamsForVariantCallers_ch.tumor
    .join(BamsForVariantCallers_ch.normal, by: [0])
    .combine(scatteredIntervalListToBed_ch.flatten())
    .map { meta, tumorbam, tumorbai, normalbam, normalbai, intervals ->
           tuple(meta, tumorbam, tumorbai, normalbam, normalbai, intervals,
            file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})

    def GatheredVardictCalls_ch = GatherVardictCalls(scatterVardict_ch
    .toSortedList({a, b -> a[1].baseName <=> b[1].baseName})
    .flatten()
    .collate(2)
    .groupTuple(by: [0])
    .map { meta, vcf -> 
        tuple(meta, vcf, file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})

    def FinalizeVardictCalls_ch = FinalizeVardictCalls(GatheredVardictCalls_ch)
    
    // VCF intersection
    def VcfIntersection_ch = VcfIntersection(FilterMutect2Input_ch
    .join(Strelka2Final_ch, by: [0])
    .join(FinalizeVarscan2Calls_ch, by: [0])
    .join(FinalizeVardictCalls_ch, by: [0]))

    // Run MSISensor pro
    // Make reference MS
    def MakeMSISensorRef_ch = MakeMSISensorRef(tuple(file(params.RefFasta), file(params.RefIdx), file(params.RefDict)))

    // Find MSI
    def runMSISensorPro_ch = runMSISensorPro(BamsForVariantCallers_ch.tumor
    .join(BamsForVariantCallers_ch.normal, by: [0])
    .combine(MakeMSISensorRef_ch))


    // run CNV
    if (params.runCNV) {

        log.info "Running CNV analysis, FACETS, ASCAT, and CNVkit"
        def IntervalCNV_ch = IntervalListtoBed_CNV(preprocessIntervalList_ch)

        // FACETS
        def FACETS_ch = cnv_facets(BamsForVariantCallers_ch.tumor
        .join(BamsForVariantCallers_ch.normal, by: [0])
        .combine(IntervalCNV_ch)
        .map { meta, tumorbam, tumorbai, normalbam, normalbai, interval -> 
                tuple(meta, tumorbam, tumorbai, normalbam, normalbai,
                        interval, file(params.CommonSNPs), file(params.CommonSNPsIdx))})
        
        // ASCAT and AlleleCounter
        // Split SNP loci based on the chromosomes
        def SplitIntervalByChr_ch = SplitIntervalByChr(file(params.SNPChr))

        // Scatter AlleleCounter for each chromosome
        def AlleleCount_ch = AlleleCount(SamtoolsSortIdxRecal_ch
        .combine(SplitIntervalByChr_ch.flatten())
        .map { meta, bam, bai, loci ->
                tuple(meta, bam, bai, loci, 
                file(params.RefFasta), file(params.RefIdx), file(params.RefDict))})

        // Gather all CHR counts to a single counts file
        def GatherAlleleCounts_ch = GatherAlleleCounts(AlleleCount_ch
        .toSortedList({a, b -> a[1].baseName <=> b[1].baseName})
        .flatten()
        .collate(2)
        .groupTuple(by: [0])
        .map { meta, counts ->
            def sorted_counts = counts.toSorted { a, b ->
                def numA = a.baseName.split('_')[0].toInteger()
                def numB = b.baseName.split('_')[0].toInteger()
                return numA <=> numB
            }
            tuple(meta, sorted_counts)})
        
        // Find counts inside the interval
        def FindCountsInsideInterval_ch = FindCountsInsideInterval(GatherAlleleCounts_ch
        .combine(IntervalCNV_ch))

        def FinalizeCounts_ch = FinalizeCounts(FindCountsInsideInterval_ch)

        // Branch counts to tumor and normal for pairing
        def AlleCountInputForConversion_ch = FinalizeCounts_ch
        .branch {
        meta_ori, count ->
            def meta = meta_ori.clone()
            tumor : meta.sampleType == "tumor"
                meta.remove('sampleType')
                return [meta, count]

            normal: meta.sampleType == "normal"
                meta.remove('sampleType')
                return [meta, count]
            }
        
        // Run Conversion to LogR and BAF file
        def ConvertAlleleCountToLogBaf_ch = ConvertAlleleCountToLogBaf(AlleCountInputForConversion_ch.tumor
        .join(AlleCountInputForConversion_ch.normal, by: [0]))

        // Run ASCAT
        def ASCAT_ch = ASCAT(ConvertAlleleCountToLogBaf_ch
        .map { meta, tumorlog, tumorbaf, normallog, normalbaf ->
               tuple(meta, tumorlog, tumorbaf, normallog, normalbaf,
                     file(params.GCContent), file(params.RepTiming))})
        
        
        // Extract ploidy and purity from FACETS and ASCAT
        def EstimatePloidyPurity_ch = EstimatePloidyPurity(FACETS_ch
        .join(ASCAT_ch, by: [0])
        .map { meta, cnv, cov, csv, spider, vcfs, pngs, baf, logr, rds, rdata, qc ->
               tuple(meta, vcfs, qc)})

        
        // CNVKIT PROCESSING
        // extract the best purity and ploidy value from ploidy purity file
        def SamplePurity_ch = EstimatePloidyPurity_ch
        .map { meta, file1, file2 -> 
            def purity = file1.text.readLines()[3].split(/\t/)[3]
            return [meta, purity as Double]
        }

        def SamplePloidy_ch = EstimatePloidyPurity_ch
        .map { meta, file1, file2 -> 
            def ploidy = file1.text.readLines()[3].split(/\t/)[2] as double
            def roundedPloidy = Math.round(ploidy)
            return [meta, roundedPloidy] 
        }

        def CNVkitOut_ch = RunCNVkit(BamsForVariantCallers_ch.tumor
        .join(BamsForVariantCallers_ch.normal, by: [0])
        .join(SamplePurity_ch, by: [0])
        .join(SamplePloidy_ch, by: [0])
        .combine(BaitsBedToIntervalList_ch)
        .map { meta, tumorbam, tumorbai, normalbam, normalbai, purity, ploidy, baits -> 
               tuple(meta, tumorbam, tumorbai, normalbam, normalbai, purity, ploidy,
                     file(params.RefFasta), file(params.RefIdx), baits)})

        log.info "Running Polysolver, LOH-HLA, and TcellExTRECT"

        // Polysoler run
        def Polysolver_ch = Polysolver(BamsForVariantCallers_ch.tumor)

        def SampleCopyNumber_ch = EstimatePloidyPurity_ch
        .map { meta, file1, file2 -> tuple(meta, file2)}

        // Run LOH-HLA
        def RunLOHHLAOut_ch = RunLOHHLA(BamsForVariantCallers_ch.tumor
        .join(BamsForVariantCallers_ch.normal, by: [0])
        .join(Polysolver_ch, by: [0])
        .join(SampleCopyNumber_ch, by: [0])
        .map { meta, tumorbam, tumorbai, normalbam, normalbai, patienthla, copynumber ->
              tuple(meta, tumorbam, tumorbai, normalbam, normalbai, patienthla,
                    file(params.HLAFasta), file(params.HLAFastaIdx), copynumber)})
                
        // Run TcellExTRECT
        def CnvKitSeg_ch = CNVkitOut_ch
        .map { meta, bin, breaks, call, cnr, cns, diagram, gainloss, scatter, segmetrics, ref ->
              tuple(meta, call)}


        def TcellExTRECT_ch = TcellExTRECT(BamsForVariantCallers_ch.tumor
        .join(CnvKitSeg_ch, by: [0])
        .join(SamplePurity_ch, by: [0]))



    } else {
        log.info "Skipping CNV process"
    }

    log.info "Workflow completed"    
}
