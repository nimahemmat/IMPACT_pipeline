# Description: 
# R script to run ASCAT using tumor and normal LogR and BAF files
# input: tumor and normal LogR and BAF file, GC content and Rep time of SNPs
# output: segmentations and ASCAT output files and plots
# By: Nima Hemmat - Peter MacCallum Cancer Centre

library(ASCAT)

args <- commandArgs(trailingOnly = TRUE)

tumor_log <- args[1]
tumor_baf <- args[2]
normal_log <- args[3]
normal_baf <- args[4]
gender <- args[5]
gc_content <- args[6]
rep_time <- args[7]
seg_penalty <- args[8]




ascat.bc <- ascat.loadData(Tumor_LogR_file = tumor_log,
                          Tumor_BAF_file = tumor_baf,
                          Germline_LogR_file = normal_log,
                          Germline_BAF_file = normal_baf,
                          gender = "XX",
                          genomeVersion = "hg38")

ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")

ascat.bc <- ascat.correctLogR(ascat.bc, GCcontentfile = gc_content, replictimingfile = rep_time)

ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")

ascat.bc <- ascat.aspcf(ascat.bc, penalty = as.numeric(seg_penalty))
ascat.plotSegmentedData(ascat.bc)

saveRDS(ascat.bc, "ascat_object.rds")

ascat.output <- ascat.runAscat(ascat.bc, gamma = 0.3, write_segments = TRUE)

QC <- ascat.metrics(ascat.bc,ascat.output)

write.table(QC, "ASCAT_QC.txt", sep = "\t")

save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')


