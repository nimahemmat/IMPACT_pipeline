args <- commandArgs(trailingOnly = TRUE)

counts_file <- args[1]
out_file <- args[2]

counts <- read.table(counts_file, sep = "\t" , header = T)

suppressMessages(library(dplyr))

# Aggregate counts for duplicates
aggregated_counts <- counts %>%
  group_by(CHR, POS) %>%
  summarise(
    Count_A = sum(Count_A, na.rm = TRUE),
    Count_C = sum(Count_C, na.rm = TRUE),
    Count_G = sum(Count_G, na.rm = TRUE),
    Count_T = sum(Count_T, na.rm = TRUE),
    Good_depth = sum(Good_depth, na.rm = TRUE),
    .groups = "drop"
  )

aggregated_counts <- aggregated_counts[ ,c("CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth")]

# Write the merged data to the output file
write.table(
  aggregated_counts,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)