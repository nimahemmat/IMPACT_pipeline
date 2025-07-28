# Description:
# R script to filter counts which are inside the bed file interval

args <- commandArgs(trailingOnly = TRUE)

# Pass arguments
counts_file <- args[1]
bed_file <- args[2]
out_file <- args[3]

counts <- read.table(counts_file, sep = "\t" , header = T)

bed_interval <- read.table(bed_file, 
                           sep = "\t",
                           col.names = c("CHR" , "START" , "END" , "V1" , "V2" , "V3"))

# Filter counts
suppressMessages(library(GenomicRanges))

counts_gr <- GRanges(seqnames = counts$CHR, ranges = IRanges(start = counts$POS, end = counts$POS))
bed_interval_gr <- GRanges(seqnames = bed_interval$CHR, ranges = IRanges(start = bed_interval$START, end = bed_interval$END))

overlaps <- findOverlaps(counts_gr, bed_interval_gr)

filtered_counts <- counts[queryHits(overlaps), ]

final_counts <- filtered_counts[ , c("CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth")]

write.table(final_counts,
            file = out_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

