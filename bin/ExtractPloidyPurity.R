# Description:
# R script to collect Ploidy and purity from FACETS and ASCAT
# Compare the values and choose best ploidy and purity value

# Passing arguments
args <- commandArgs(trailingOnly=TRUE)

facets_vcf <- args[1]
ascat_qc <- args[2]
output_file <- args[3]
output_file_2 <- args[4]
patient_id <- args[5]

# Read facets vcf file

suppressMessages(library(VariantAnnotation))

facets <- readVcf(facets_vcf, genome = "hg38")
ascat <- read.table(ascat_qc, sep = "\t" , header = T)

popu <- data.frame(
  Tool = c("FACETS", "ASCAT"),
  ploidy = c(round(as.numeric(facets@metadata[["header"]]@header@listData[["ploidy"]]@listData[["Value"]]), 2),
             ascat$ploidy),
  purity = c(round(as.numeric(facets@metadata[["header"]]@header@listData[["purity"]]@listData[["Value"]]), 2),
             ascat$purity)
)

# Choose ASCAT purity and ploidy with FACETS support
best_ploidy <- NA
best_purity <- NA

if(is.na(popu[2,3]) | popu[2,3] == 1 | popu[2,2] <= 2){
  best_ploidy <- popu[1,2]
  best_purity <- popu[1,3]
} else {
  best_ploidy <- popu[2,2]
  best_purity <- popu[2,3]
}


best_popu <- data.frame(
  Tool = "BEST",
  ploidy = best_ploidy,
  purity = best_purity
)

write.table(rbind(popu, best_popu), file = output_file, sep = "\t" , quote = FALSE , col.names = TRUE)

# Make file ready for LOH-HLA

copynumber <- data.frame(
  Ploidy = best_popu$ploidy,
  tumorPurity = best_popu$purity,
  tumorPloidy = best_popu$ploidy
)

rownames(copynumber) <- patient_id

write.table(copynumber, file = output_file_2, sep = "\t", quote = FALSE, col.names = TRUE , row.names = TRUE)


