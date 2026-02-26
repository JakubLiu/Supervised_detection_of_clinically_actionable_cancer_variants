suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))

library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)

annotate_reads <- function(bam, chrom, pos, alt, output_file_name = 0){
    
    chrom <- as.character(chrom)
    pos <- as.numeric(pos)
    alt <- as.character(alt)

    param <- ScanBamParam(which = GRanges(chrom, IRanges(pos, pos)), what = c("qname", "seq", "qual", "mapq"))
    gal <- readGAlignments(bam, param = param)
    locus <- GRanges(chrom, IRanges(pos, pos))
    names(gal) <- mcols(gal)$qname
    mapped <- mapToAlignments(locus, gal) # fetch the reads that map to the locus
    seqs <- mcols(gal)$seq # the sequences of the reads
    alt_read_mask <- mapply(function(s, p) {as.character(subseq(s, p, p))}, seqs, start(mapped)) == alt
    alt_reads <- gal[alt_read_mask] # extract only the alternative reads
    alt_reads_seq <- as.character(mcols(alt_reads)$seq) # the sequences of the full reads
    mapq <- mcols(alt_reads)$mapq  # extract mapping qualities
    position_in_read <- start(mapped)[alt_read_mask]  # position in the read
    read_length <- width(mcols(alt_reads)$seq)  # length of the read
    strand <- as.character(strand(alt_reads))  # strand "+" = R1, "-" = R2
    basecallig_qualities <- mapply(function(q, p) {as.integer(subseq(q, p, p))}, mcols(alt_reads)$qual, position_in_read)  # as.integer automatically changes the ASCII Phred score to the numerical score
    
    output <- data.frame(cbind(alt_reads_seq, mapq, basecallig_qualities, read_length, position_in_read, strand))
    colnames(output) <- c('full_read_sequence', 'mapping_quality', 'basecalling_quality', 'read_length',
                          'position_within_the_read', 'strand')

    if(output_file_name != 0){
        output_file_name <- as.character(output_file_name)
        write.csv(output, output_file_name, row.names = FALSE)
    }

    return(output)
}

# USAGE:_________________________________________________________________________________________________________________________

#annotate_reads(
#    bam = 'T1-DNA1-WES1.mutated.sorted.bam',
#    chrom = '7',
#    pos = 55259515,
#    alt = 'G',
#    output_file_name = 'test_annotation.csv'
#)