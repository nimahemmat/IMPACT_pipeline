# Decription:
# R script to run TcellExTRECT tool
# The script needs bed file, tumor bam file and CNVkit call output
# by: Nima Hemmat
# Peter MacCallum Cancer Centre, Melbourne, Australia

args <- commandArgs(trailingOnly=TRUE)

bam_file <- args[1]
sample_name <- args[2]
median_thresh <- args[3]
purity <- args[4]
cnvkit_seg <- args[5]
out_file <- args[6]
out_path <- args[7]

## Find TCRA copy number from CNVkit seg file
suppressMessages(library(GenomicRanges))

cnvkit_file <- read.table(cnvkit_seg, sep = "\t" , header = T)

tcra_pos <- data.frame(
  Gene = "TRA",
  Chr = "14",
  Start = 21621904,
  End = 22552131
)

cnvkit_gr <- GRanges(seqnames = cnvkit_file$chromosome, ranges = IRanges(start = cnvkit_file$start, end = cnvkit_file$end))
tcra_gr <- GRanges(seqnames = tcra_pos$Chr, ranges = IRanges(start = tcra_pos$Start, end = tcra_pos$End))

overlaps <- findOverlaps(cnvkit_gr, tcra_gr)

seg_tcra <- cnvkit_file[queryHits(overlaps),]

tcra_cn <- ifelse(length(seg_tcra$cn) == 0, 2, seg_tcra$cn)

## Run T-cell ExTRECT
suppressMessages(library(TcellExTRECT))
data(tcra_seg_hg38)
data(TCRA_exons_hg38)

# Step 1: Generate coverage file
cov.file <- getCovFromBam(
  bamPath = bam_file,
  outPath = out_path,
  vdj.seg = tcra_seg_hg38
)

# Step 2: Load coverage data
cov_df <- loadCov(cov.file)

# Step 3: Run TcellExTRECT
TCRA.out <- runTcellExTRECT(
  vdj.region.df = cov_df,
  exons.selected = TCRA_exons_hg38, 
  vdj.seg = tcra_seg_hg38, 
  hg19_or_38 = "hg38",
  median.thresh = as.numeric(median_thresh), 
  GC_correct = T, 
  median.k = 20,
  sample_name = sample_name
)

# Step 4: Adjust the output
TCRA.out_adj <- adjustTcellExTRECT(
  TCRA.out = TCRA.out, 
  purity = as.numeric(purity), 
  TCRA.cn = tcra_cn
)

final_estimate <- as.data.frame(TCRA.out_adj)

write.table(
  final_estimate,
  file = out_file,
  quote = FALSE, 
  sep = "\t", 
  row.names = TRUE,
  col.names = TRUE
)



