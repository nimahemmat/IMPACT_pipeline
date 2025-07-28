args <- commandArgs(trailingOnly = TRUE)

# First argument: Output file
output_file <- args[1]

# Remaining arguments: Input files
input_files <- args[-1]

# Convert the input files to a list of files
input_list <- as.list(input_files)

# Read each input file
files_list <- lapply(input_list, function(file) {
  read.table(file, sep = "\t", stringsAsFactors = FALSE)
})


# merge into one file
merged_counts <- do.call(rbind, files_list)

# Change column names
colnames(merged_counts) <- c("CHR","POS","Count_A","Count_C","Count_G",
                             "Count_T","Good_depth")


suppressMessages(library(dplyr))

# Aggregate counts for duplicates
aggregated_counts <- merged_counts %>%
  group_by(CHR, POS) %>%
  summarise(
    Count_A = sum(Count_A, na.rm = TRUE),
    Count_C = sum(Count_C, na.rm = TRUE),
    Count_G = sum(Count_G, na.rm = TRUE),
    Count_T = sum(Count_T, na.rm = TRUE),
    Good_depth = sum(Good_depth, na.rm = TRUE),
    .groups = "drop"
  )

# Write the merged data to the output file
write.table(
  aggregated_counts,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)