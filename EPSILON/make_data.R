


suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

is_homopolymer <- function(sequence, nuc, n = 3){
  pattern <- paste0(nuc, "{", n, ",}")
  return(grepl(pattern, sequence))  # return boolean
}


max_homo_run <- function(sequence, base) {
  pattern <- paste0(base, "+")
  matches <- gregexpr(pattern, sequence)[[1]]
  if (matches[1] == -1) return(0)
  max(attr(matches, "match.length"))
}


create_data <- function(bamlist,    # a list of bam files, one path per line
                        reference_genome,  # the path to the reference genome
                        chrom_list, pos_list, ref_list, alt_list,  # R vectors
                        output,  # path to output file
                        ncores,
                        tmpdir
                        ){   


    logfile <- paste0(output, '.log') # here the progress will be printed
    full_windowsize <- 100

    cl <- makeCluster(ncores)
    registerDoParallel(cl)


    bamlist <- readLines(bamlist)
    reference_genome <- FaFile(reference_genome)
    num_loci <- length(chrom_list)


    # outer loop over the loci
    all_data <- foreach(i = 1:num_loci,
                    .packages = c('Rsamtools', 'Biostrings', 'stats', 'GenomicAlignments', 'GenomicRanges', 'data.table', 'stats'),
                                        .combine = rbind) %dopar%{

                        # these functions need to be defined here, because they must be visible to each worker
                        is_homopolymer <- function(sequence, nuc, n = 3){
                                pattern <- paste0(nuc, "{", n, ",}")
                                return(grepl(pattern, sequence))  # return boolean
                        }


                        max_homo_run <- function(sequence, base) {
                                pattern <- paste0(base, "+")
                                matches <- gregexpr(pattern, sequence)[[1]]
                                if (matches[1] == -1) return(0)
                                max(attr(matches, "match.length"))
                        }

                    
                    cat(paste("Worker", Sys.getpid(), "progress: ", i/num_loci*100, "%\n"), file = logfile, append = TRUE) # write the progress to the log file
                    # fix the locus specific variables
                    current_chromosome <- as.character(chrom_list[i])
                    current_pos <- as.numeric(pos_list[i])
                    current_ref <- as.character(ref_list[i])
                    current_alt <- as.character(alt_list[i])
                    region <- GRanges(current_chromosome, IRanges(start = current_pos, end = current_pos))
                    scan_bam_param <- ScanBamParam(which = region, what = c("qname", "seq", "qual", "mapq"), tag = 'NM')
                    pileup_param <- PileupParam(distinguish_strands = TRUE, distinguish_nucleotides = TRUE)
                    partial_output <- paste0(tmpdir, '/' ,output, "_chr_",current_chromosome,'_pos_',current_pos, Sys.getpid(), ".csv")

                    # extract the different sized contexts
                    context_full <- GRanges(seqnames = current_chromosome, ranges = IRanges(start = current_pos - full_windowsize/2, end = current_pos + full_windowsize/2))
                    context_full <- toupper(as.character(getSeq(reference_genome, context_full))) # for GC, GGT, and GGC
                    homopolymer_context <- substr(context_full, 35, nchar(context_full) - 35)  # for all homopolymer features
                    upstream_base_context <- context_full[(full_windowsize/2-3):(full_windowsize/2)]  # for the three bases upstream
                    downstream_base_contect <- context_full[(full_windowsize/2):(full_windowsize/2+3)]  # for the three bases downstream


                    # inner sequential loop over the sample
                    
                    chrom <- c()  # not for the model, just for sanity checks
                    pos <- c()    # not for the model, just for sanity checks
                    ref <- c()     # not for the model, just for sanity checks
                    alt <- c()    # not for the model, just for sanity checks
                    GGC_upstream <- c()    # done
                    CGG_downstream <- c()   # GGC read from the reverse direction
                    GGT_upstream <- c()    # done
                    TGG_downstream <- c()  # GGT read from the reverse direction
                    GC_content <- c()    # done
                    homo_percentage <- c()   # done
                    homo_overlap <- c()  # done
                    homo_end <- c()  # done
                    homo_start <- c()  # done

                    # set the feature/y vectors
                    specific_alt_counts <- c()    # done
                    general_alt_counts <- c()     # done
                    ref_counts <- c()             # done
                    sampleIDs <- c()              # done
                    mismatches <-c()             # any mismatch              done
                    ranksum_mapQ_sig <- c()       # only specific alt        done
                    ranksum_baseQ_sig <- c()      # only specific alt         done
                    fisher_sig <- c()       # only specific alt (proxy for strand bias)  done
                    median_baseQ <- c()        # median of the baseQ of all reads (ref or alt in both directions)
                    median_mapQ <- c()         # median of the mapQ of all reads (ref or alt in both directions)


                    for(j in 1:length(bamlist)){

                        # setting some parameters______________________________________________________________________
                        current_bam <- bamlist[j]
                        bam_basename <- strsplit(current_bam, "/")[[1]]
                        bam_basename <- bam_basename[length(bam_basename)]
                        local_pileup <- pileup(file = current_bam, index = paste0(current_bam, '.bai'), scanBamParam = scan_bam_param, pileupParam = pileup_param)
                        coverage <- sum(local_pileup$count)
                        ref_count <- local_pileup[local_pileup$nucleotide == current_ref,]$count
                        ref_count <- sum(ref_count)
                        chrom <- c(chrom, current_chromosome)
                        pos <- c(pos, current_pos)
                        ref <- c(ref, current_ref)
                        alt <- c(alt, current_alt)


                        # sample-specific features_________________________________________________________________________
                        specific_alt_count <- sum(local_pileup[local_pileup$nucleotide == current_alt,]$count)
                        general_alt_count <- coverage - ref_count
                        specific_alt_counts <- c(specific_alt_counts, specific_alt_count)
                        general_alt_counts <- c(general_alt_counts, general_alt_count)
                        ref_counts <- c(ref_counts, ref_count)
                        sampleIDs <- c(sampleIDs, bam_basename)
                        gal <- readGAlignments(current_bam, param = scan_bam_param)
                        hits <- findOverlaps(gal, GRanges(current_chromosome, IRanges(current_pos, current_pos)))
                        gal_hits <- gal[queryHits(hits)]  # reads mapped to the locus of interest
                        names(gal_hits) <- mcols(gal_hits)$qname
                        mapped <- mapToAlignments(GRanges(current_chromosome, IRanges(current_pos, current_pos)), gal_hits) # fetch the reads that map to the locus
                        positions_in_reads <- start(mapped) # in what position relative to each read the locus of interest is
                        whole_read_qualities <- mcols(gal)$qual   # extract basecalling qualities for the whole reads


                        # subset the basecalling qualities only for the regions in each read that overlap with the target locus
                        baseQ_at_read_position <- mapply(function(qual, pos){
                                                                as.numeric(subseq(qual, pos, pos))},
                                                                whole_read_qualities, positions_in_reads)
                        
                        mismatches <- c(mismatches, mean(mcols(gal_hits)$NM>0)) # fraction of reads with mismatches
                        
                        # get base at locus for each read
                        read_seqs <- mcols(gal_hits)$seq

                        bases_at_pos <- mapply(function(seq, pos) {
                        as.character(subseq(seq, pos, pos))
                        }, read_seqs, positions_in_reads)

                        local_strand <- as.character(strand(gal_hits))

                        local_info <- data.frame(
                                base = bases_at_pos,
                                mapq = mcols(gal_hits)$mapq, # mapping qualities of the reads 'over' the target locus
                                baseq = baseQ_at_read_position,  # basecalling qualities of the reads 'over' the target locus
                                strand = local_strand  # the directions of the reads 'over' the current locus
                                )

                        
                        
                        # the mapQ/baseQ ref/alt vectors
                        ref_baseQ <- local_info[local_info$base == current_ref,]$baseq
                        alt_baseQ <- local_info[local_info$base == current_alt,]$baseq
                        ref_mapQ <- local_info[local_info$base == current_ref,]$mapq
                        alt_mapQ <- local_info[local_info$base == current_alt,]$mapq
                        all_baseQ <- local_info$baseq
                        all_mapQ <- local_info$mapq

                        median_baseQ <- c(median_baseQ, median(all_baseQ, na.rm = TRUE))
                        median_mapQ <- c(median_mapQ, median(all_mapQ, na.rm = TRUE))



                        # counts of ref/alt R1/R2 reads (for the Fisher test)
                        ref_R1 <- sum(local_info$base == current_ref & local_info$strand == '+')
                        ref_R2 <- sum(local_info$base == current_ref & local_info$strand == '-')
                        alt_R1 <- sum(local_info$base == current_alt & local_info$strand == '+')
                        alt_R2 <- sum(local_info$base == current_alt & local_info$strand == '-')


                        # mapQ rank sum test, shows the significance of the difference in the mapQ of the ref and alt reads
                        #if(length(ref_mapQ) > 0 && length(alt_mapQ) > 0){
                        #        pval_mapQ <- wilcox.test(ref_mapQ, alt_mapQ)$p.value
                        #        if(is.na(pval_mapQ) == TRUE){
                        #                ranksum_mapQ_sig <- c(ranksum_mapQ_sig, 0)
                        #        }
                        #        else if(pval_mapQ <= 0.05){
                        #                ranksum_mapQ_sig <- c(ranksum_mapQ_sig, 1)
                        #        }
                        #        else{
                        #                ranksum_mapQ_sig <- c(ranksum_mapQ_sig, 0)
                        #        }
                        #}
                        #else{
                        #        ranksum_mapQ_sig <- c(ranksum_mapQ_sig, 0)
                        #}

                        ranksum_mapQ_sig <- c(ranksum_mapQ_sig, 0)  # discard this later


                        # baseQ rank sum test, shows the significance of the difference in the baseQ of the ref and alt reads
                        #if(length(alt_baseQ) > 0 && length(ref_baseQ) > 0){

                                #ifelse(wilcox.test(ref_baseQ, alt_baseQ)$p.value <= 0.05, ranksum_baseQ_sig <- c(ranksum_baseQ_sig, 1), ranksum_baseQ_sig <- c(ranksum_baseQ_sig, 0))
                        
                        #}
                        #else{
                        #        ranksum_baseQ_sig <- c(ranksum_baseQ_sig, 0)
                        #}


                        ranksum_baseQ_sig <- c(ranksum_baseQ_sig, 0)   # discard this later

                        #if(ref_R1+ref_R2>=1){
                        #        tab <- matrix(c(ref_R1, ref_R2, alt_R1, alt_R2), nrow = 2, ncol = 2)
                        #        ifelse(fisher.test(tab)$p.value <= 0.05, fisher_sig <- c(fisher_sig, 1), fisher_sig <- c(fisher_sig, 0))
                        #}
                        #else{
                        #        fisher_sig <- c(fisher_sig, 0)
                        #}


                        fisher_sig <- c(fisher_sig, 0)  # discard this later


                        #  context specific features_______________________________________________________________________________
                        center <- ceiling(nchar(context_full) / 2)

                        # GGC upstream
                        if (grepl("GGC", substr(context_full, 1, center))) {
                        GGC_upstream <- c(GGC_upstream, 1)
                        } else {
                        GGC_upstream <- c(GGC_upstream, 0)
                        }

                        # CGG downstream
                        if (grepl("CGG", substr(context_full, center, nchar(context_full)))) {
                        CGG_downstream <- c(CGG_downstream, 1)
                        } else {
                        CGG_downstream <- c(CGG_downstream, 0)
                        }

                        # GGT upstream
                        if (grepl("GGT", substr(context_full, 1, center))) {
                        GGT_upstream <- c(GGT_upstream, 1)
                        } else {
                        GGT_upstream <- c(GGT_upstream, 0)
                        }

                        # TGG downstream
                        if (grepl("TGG", substr(context_full, center, nchar(context_full)))) {
                        TGG_downstream <- c(TGG_downstream, 1)
                        } else {
                        TGG_downstream <- c(TGG_downstream, 0)
}
                        GC_content <- c(GC_content, sum(strsplit(context_full, "")[[1]] %in% c("G", "C"))/nchar(context_full))
                    
                        # define midpoints safely
                        seq_len <- nchar(homopolymer_context)
                        mid <- ceiling(seq_len / 2)

                        # full homopolymer percentage
                        if(is_homopolymer(homopolymer_context, current_alt, 3)){
                        homo_percentage <- c(homo_percentage, max_homo_run(homopolymer_context, current_alt)/seq_len)
                        } else {
                        homo_percentage <- c(homo_percentage, 0.0)
                        }

                        # homo_overlap: 3-base center
                        start_overlap <- max(1, mid-1)
                        end_overlap   <- min(seq_len, mid+1)
                        mid_seq <- substr(homopolymer_context, start_overlap, end_overlap)
                        if(is_homopolymer(mid_seq, current_alt, 3)){
                        homo_overlap <- c(homo_overlap, 1)
                        } else {
                        homo_overlap <- c(homo_overlap, 0)
                        }

                        # homo_end: 4 bases ending at center
                        start_end <- max(1, mid-3)
                        end_end   <- mid
                        end_seq <- substr(homopolymer_context, start_end, end_end)
                        if(is_homopolymer(end_seq, current_alt, 3)){
                        homo_end <- c(homo_end, 1)
                        } else {
                        homo_end <- c(homo_end, 0)
                        }

                        # homo_start: 4 bases starting at center
                        start_start <- mid
                        end_start   <- min(seq_len, mid+3)
                        start_seq <- substr(homopolymer_context, start_start, end_start)
                        if(is_homopolymer(start_seq, current_alt, 3)){
                        homo_start <- c(homo_start, 1)
                        } else {
                        homo_start <- c(homo_start, 0)
                        }
                
                    }   # end of loop over the samples



                    data <- data.frame(cbind(chrom, pos, ref, alt, specific_alt_counts, general_alt_counts, ref_counts,
                                sampleIDs, mismatches, ranksum_mapQ_sig, ranksum_baseQ_sig, fisher_sig,
                                GGC_upstream, CGG_downstream, GGT_upstream, TGG_downstream, GC_content,
                                homo_percentage, homo_overlap, homo_end, homo_start, median_mapQ, median_baseQ))

                    colnames(data) <- c('chrom', 'pos', 'ref', 'alt', 'specific_alt_counts', 'general_alt_counts', 'ref_counts',
                                'sampleIDs', 'mismatches', 'ranksum_mapQ_sig', 'ranksum_baseQ_sig', 'fisher_sig',
                                'GGC_upstream', 'CGG_downstream', 'GGT_upstream', 'TGG_downstream', 'GC_content',
                                'homo_percentage', 'homo_overlap', 'homo_end', 'homo_start',
                                'median_mapQ', 'median_baseQ')

                    write.csv(data, partial_output, row.names = FALSE)
                    

                    return(data)


            }
    
    write.csv(all_data, file = output, row.names = FALSE)
    stopCluster(cl)
    return(all_data)
    
}


# __________________________________________________ example usage _______________________________________


#please <- create_data(bamlist = '/data/cephfs-1/home/users/jali13_c/work/EPSILON/REGRESSION/ZIBINOMIAL/SYSTEMATIC_VALIDATION/locus1_300/simVAF_0.005/bamlist.txt',
#                        reference_genome = '/data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/TUNE_PRIOR/STEP1_PREPROCESS_BAMS/hs37d5.fa',
#                        chrom_list = c('9', '14'),
#                        pos_list = c(133748283, 105246551),
#                        ref_list = c('C', 'C'),
#                        alt_list = c('T', 'T'),  # R vectors 
#                        output = 'test.csv',  # path to output file
#                        ncores = 2,
#                        tmpdir = 'tmp_test'
#                        )


#print(head(please))

