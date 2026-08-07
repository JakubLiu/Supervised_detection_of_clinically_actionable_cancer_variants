library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_file <- as.character(args[1])
output_file <- as.character(args[2])
subsampling_strength <- as.numeric(args[3])
data <- fread('merged.processed.mpileup.txt')
nrows <- nrow(data)*subsampling_strength
idx <- sample(1:nrow(data), nrows, replace = FALSE)
data_subsampled <- data[idx,]
write.table(data_subsampled, file = output_file, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
