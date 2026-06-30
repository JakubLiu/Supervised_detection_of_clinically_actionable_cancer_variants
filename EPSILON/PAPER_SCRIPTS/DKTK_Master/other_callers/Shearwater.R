library(deepSNV)
library(GenomicRanges)
library(VariantAnnotation)
library(Biostrings)

# ==================================== how to run ====================================================
# Rscript Shearwater.R \
#	T1-DNA1-WES1.sorted.bam \     ---> the bam files
#	normal_list.minimal.txt \     ---> the list of paths to the normal bam files (one path per line)
#	Civic.small.bed \             ---> the bed file
#	out.csv                      ---> the output file
# ===================================================================================================

options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)

tumor_bam <- as.character(args[1])
panel_of_normals_bamlist <- as.character(args[2])
bedfile <- as.character(args[3])
output_file <- as.character(args[4])

panel_of_normals <- readLines(panel_of_normals_bamlist)

bed <- read.table(
  bedfile,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(bed)[1:3] <- c("chr", "start", "end")

padding <- 1


SAMPLE <- c()
CHROM <- c()
POS <- c()
REF <- c()
ALT <- c()



for (i in seq_len(nrow(bed))) {


  chr <- bed$chr[i]
  start <- bed$start[i]
  end <- bed$end[i]

  region <- GRanges(
    seqnames = chr,
    ranges = IRanges(
      start = max(1, start + 1 - padding),
      end   = end + padding
    )
  )

  message("Processing: ", chr, ":", start, "-", end)

  

  # the minimum basecalling qual (q param) is lefr default
  data <- deepSNV:::loadAllData(
    c(tumor_bam, panel_of_normals),
    region
  )

  
  if (is.null(data) || length(data) == 0) {
    next
  }


  bf <- deepSNV:::bbb(data, model = "OR")

  

  vcf <- deepSNV:::bf2Vcf(
    bf,
    data,
    region,
    cutoff = 0.05,  # default value
    samples = c(tumor_bam, panel_of_normals),
    prior = 0.5,   # default value
    mvcf = TRUE
  )


  if (!is.null(vcf) && length(rowRanges(vcf)) > 0) {

    gr <- rowRanges(vcf)

    SAMPLE <- c(SAMPLE, tumor_bam)
    CHROM <- c(CHROM, chr)
    POS <- c(POS, start(gr))
    REF <- c(REF, as.character(ref(vcf)))
    ALT <- c(ALT, sapply(alt(vcf), function(x) paste(as.character(x), collapse = ",")))

  }

}

results <- data.frame(
    cbind(
        SAMPLE,CHROM,POS,REF,ALT
    )
)

colnames(results) <- c(
    'SAMPLE', 'CHROM', 'POS', 'REF', 'ALT'
)


write.csv(results, output_file, row.names = FALSE)

