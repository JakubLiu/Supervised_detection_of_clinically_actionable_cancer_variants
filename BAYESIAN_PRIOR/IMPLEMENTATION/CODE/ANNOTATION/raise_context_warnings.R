suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(GenomicRanges))

library(Rsamtools)
library(Biostrings)
library(GenomicRanges)

raise_context_warnings <- function(reference_genome, chr, pos, alt, padding_upstream = 10, padding_downstream = 10,
                                    GC_content_min = 0.35, # lower bound of the 'normal' GC content [2]
                                    GC_content_max = 0.60, # upper bound for the 'normal' GC content [2]
                                    min_upstream_homopolymer_length = 5,  # look at the explanation below in the code
                                    output_file_name){  

    chr <- as.character(chr)
    pos <- as.numeric(pos)
    alt <- as.character(alt)
    padding_upstream <- as.numeric(padding_upstream)
    padding_downstream <- as.numeric(padding_downstream)
    refgen <- FaFile(reference_genome) # the fasta must also have the index (.fa.fai)
    open(refgen)
    context_window <- GRanges(seqnames = chr, ranges = IRanges(start = pos - padding_upstream, end = pos + padding_downstream))
    context_window_GC <- GRanges(seqnames = chr, ranges = IRanges(start = pos - 50000, end = pos + 50000)) # extract a 100Kb fragment for computing the GC content [2]
    context <- as.character(getSeq(refgen, context_window))  # get the genomic context
    context_GC <- as.character(getSeq(refgen, context_window_GC)) # the 10Kb genomic context for the GC content calculation
    close(refgen)
    upstream_context <- substr(context, 1, padding_upstream) # the genomic context upstream of the target locus
    downstream_context <- substr(context, padding_upstream + 2, nchar(context)) # the genomic context downstream of the target locus


    # search for error linked motifs upstream of the target locus [1]
    GG_motif_upstream <- FALSE
    if(grepl(("GG"), upstream_context)){
        GG_motif_upstream <- TRUE
    }
    GGC_motif_upstream <- FALSE
    if(grepl(("GGC"), upstream_context)){
        GGC_motif_upstream <- TRUE
    }
    GGT_motif_upstream <- FALSE
    if(grepl(("GGT"), upstream_context)){
        GGT_motif_upstream <- TRUE
    }

    GC_content <- sum(strsplit(context_GC, "")[[1]] %in% c("G", "C"))/nchar(context_GC)
    abnormal_GC_content <- FALSE
    if((GC_content < GC_content_min) | (GC_content > GC_content_max)){
        abnormal_GC_content <- TRUE
    }
    
    # check if a homopolymer of the same base as the alternative, terminates exactly at the alternative [3]
    #
    # example of a positive case:
    # AAAAAAAAAAAATGCTACC
    #            ^
    #            |
    #       locus of interest
    #
    potential_homopolymer_upstream_region <- substr(upstream_context, nchar(upstream_context)-min_upstream_homopolymer_length, nchar(upstream_context))
    homopolymer <- strrep(alt, min_upstream_homopolymer_length)
    homopolymer_upstream <- FALSE
    if(potential_homopolymer_upstream_region == homopolymer){
        if(substr(downstream_context,1,1) != alt){
            homopolymer_upstream <- TRUE
        }
    }

    result <- list(
        'GG_motif_upstream_present' = GG_motif_upstream,
        'GGT_motif_upstream_present' = GGT_motif_upstream,
        'GGC_motif_upstream_present' = GGC_motif_upstream,
        'abnormal_GC_content_warning' = abnormal_GC_content,
        'preceeding_homopolymer_present' = homopolymer_upstream
    )

    write.csv(data.frame(result), output_file_name, row.names = FALSE)
}

# [1] Meacham, F., Boffelli, D., Dhahbi, J., Martin, D. I., Singer, M., & Pachter, L. (2011). Identification and correction 
#     of systematic error in high-throughput sequence data. BMC bioinformatics, 12, 451. 
#     https://doi.org/10.1186/1471-2105-12-451
#
# [2] https://en.wikipedia.org/wiki/GC-content#Among-genome_variation
#
# [3] Stoler N, Nekrutenko A. Sequencing error profiles of Illumina sequencing instruments. NAR Genom Bioinform. 
#     2021 Mar 27;3(1):lqab019. doi: 10.1093/nargab/lqab019. PMID: 33817639; PMCID: PMC8002175.

# usage:_______________________________________________________________________

#raise_context_warnings(
#                            reference_genome = 'hs37d5.fa',
#                            chr = '7',
#                            alt = 'G',
#                            pos = 55259515,
#                            padding_downstream = 10,
#                            padding_upstream = 10,
#                            output_file_name = 'tesc.csv'
#)