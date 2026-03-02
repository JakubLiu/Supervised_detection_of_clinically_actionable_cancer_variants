suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))

make_pileups <- function(tumor_bam_path, negative_control_bamlist, CHROM, START, END, write_files = FALSE){

    print(paste0(
        'tumor bam: ',
        tumor_bam_path
    ))

    print(paste0(
        'negative control bamlist: ',
        negative_control_bamlist
    ))


    print('setting pileup parameters...')
    CHROM <- as.character(CHROM)
    START <- as.numeric(START)
    END <- as.numeric(END)


    actionable_region <- GRanges(CHROM, IRanges(start = START, end = END))
    scan_bam_param <- ScanBamParam(which = actionable_region)
    pileup_param <- PileupParam(distinguish_strands = TRUE, distinguish_nucleotides = TRUE)
    print('setting pileup parameters done.')


    # pileup on the target tumor sample
    print('creating pileups...')
    target_tumor_pileup <- pileup(file = tumor_bam_path, index = paste0(tumor_bam_path, ".bai"), scanBamParam = scan_bam_param, pileupParam = pileup_param)
    target_tumor_pileup_R1 <- target_tumor_pileup[target_tumor_pileup$strand == '+',]  # the pileup in the target tumor sample, only the forward reads
    target_tumor_pileup_R2 <- target_tumor_pileup[target_tumor_pileup$strand == '-',]  # the pileup in the target tumor sample, only the backward reads
    print('creating pileups done.')

    # pileup in the negative control cohort
    print('creating the negative control cohort pileup...')
    negative_control_list <- readLines(negative_control_bamlist)
    negative_control_pileup_R1 <- NULL
    negative_control_pileup_R2 <- NULL
    
    for(i in 1:length(negative_control_list)){

        current_bam <- negative_control_list[i]
        current_pileup <- pileup(file = current_bam, index = paste0(current_bam, ".bai"), scanBamParam = scan_bam_param, pileupParam = pileup_param)
        current_pileup_R1 <- current_pileup[current_pileup$strand == '+', ]
        current_pileup_R2 <- current_pileup[current_pileup$strand == '-', ]
        negative_control_pileup_R1 <- rbind(negative_control_pileup_R1, current_pileup_R1)
        negative_control_pileup_R2 <- rbind(negative_control_pileup_R2, current_pileup_R2)

    }


    negative_control_pileup_R1 <- negative_control_pileup_R1[2:nrow(negative_control_pileup_R1),]
    negative_control_pileup_R2 <- negative_control_pileup_R2[2:nrow(negative_control_pileup_R2),]
    print('creating the negative control cohort pileup done.')

    if(write_files == TRUE){

        current_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

        write.csv(target_tumor_pileup_R1, paste0(current_time, '_T.R1.csv'), row.names = F)
        write.csv(target_tumor_pileup_R2, paste0(current_time, '_T.R2.csv'), row.names = F)
        write.csv(negative_control_pileup_R1, paste0(current_time, '_N.R1.csv'), row.names = F)
        write.csv(negative_control_pileup_R2, paste0(current_time, '_N.R2.csv'), row.names = F)
    }

    return(list(
        'target_tumor_pileup_R1' = target_tumor_pileup_R1,
        'target_tumor_pileup_R2' = target_tumor_pileup_R2,
        'negative_control_pileup_R1' = negative_control_pileup_R1,
        'negative_control_pileup_R2' = negative_control_pileup_R2
    ))
    print('all done.')

    
}



estimate_and_count <- function(TR1,   # pileup of the tumor sample in the forward direction
                               TR2,   # pileup of the tumor sample in the reverse direction
                               NR1,   # pileup of the negative cohort in the forward direction
                               NR2,   # pileup of the negative cohort in the reverse direction
                               alternative_allele,  # the alternative allele as is in the actionable variant
                               write_file = FALSE  # optionally write the estimates to a file
                               ){
    print('estimating parameters from pileups...')
    alternative_allele <- as.character(alternative_allele)
    # the counts _______________________________________________________________________________________________
    # the n_i 's in the paper
    covTR1 <- sum(TR1$count)  # the coverage in the forrward direction in the tumor sample
    covTR2 <- sum(TR2$count)  # the coverage in the reverse direction in the tumor sample
    covNR1 <- sum(NR1$count)  # the total coverage (sum over all samples) in the negative cohort in the forward direction
    covNR2 <- sum(NR2$count)  # the total coverage (sum over all samples) in the negative cohort in the reverse direction

    if(covTR1 == 0 || covTR2 == 0 || covNR1 == 0 || covNR2 == 0){
      stop('Error. The coverage in at least one of the directions is either the target sample or the negative control cohort is zero.')
    }

    # the X_i 's in the paper
    countTR1 <- sum(TR1$count[TR1$nucleotide == alternative_allele]) # count of the alternative allele in the tumor in the forward direction
    countTR2 <- sum(TR2$count[TR2$nucleotide == alternative_allele]) # count of the alternative allele in the tumor in the reverse direction
    countNR1 <- sum(NR1$count[NR1$nucleotide == alternative_allele])  # the total alternative allele count (sum over all samples) in the negative cohort in the forward direction
    countNR2 <- sum(NR2$count[NR2$nucleotide == alternative_allele])  # the total alternative allele count (sum over all samples) in the negative cohort in the forward direction

    # the estimates______________________________________________________________________________________________
    vi_R1 <- countNR1/covNR1
    vi_R2 <- countNR2/covNR2
    v0_R1 <- (countTR1 + countNR1)/(covTR1 + covNR1)
    v0_R2 <- (countTR2 + countNR2)/(covTR2 + covNR2)

    mu0_R1 <- countTR1/covTR1
    mu0_R2 <- countTR2/covTR2
    mu_i <- max(c(
                (countTR1 + countTR2)/(covTR1 + covTR2),
                vi_R1,
                vi_R2))



    estimates <- list(
                        'covTR1' = covTR1, 'covTR2' = covTR2, 'covNR1' = covNR1, 'covNR2' = covNR2,
                        'countTR1' = countTR1, 'countTR2' = countTR2, 'countNR1' = countNR1, 'countNR2' = countNR2,
                        'vi_R1' = vi_R1, 'vi_R2' = vi_R2,
                        'mu0_R1' = mu0_R1, 'mu0_R2' = mu0_R2
                    )

    if(write_file == TRUE){
        current_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
        writeLines(
                        paste(names(estimates), estimates, sep = ','),
                        paste0(current_time, '_estimates.csv')
                  )
    }

    return(estimates)
    print('estimating parameters from pileups done.')

}


calc_posterior <- function(countTR1, countTR2,  # the allele counts
                        covTR1, covTR2,  # the coverages
                        mu0_R1, mu0_R2,  # the tumor VAF
                        vi_R1, vi_R2,   # the backgorund error rate
                        prior,  # the user defined prior
                        pseudo = 0.00001  # a pseudo background error rate to use if the empitical one is 0.0
                        ){
      
      # if the empirical background error rate is 0, substitute it by a pseudocount
      vi_R1 <- max(vi_R1, pseudo)
      vi_R2 <- max(vi_R2, pseudo)

      # mu0 ---> the observed tumor VAF
      # vi ---> the estimated background error rate
      # BF = P(data|H1)/P(data|H0), where data are the observed allele counts 
      print('calculating bayes factor...')

      prob_H0_R1 <- dbinom(x = countTR1, size = covTR1, prob = vi_R1, log = TRUE) # probability of observing the tumor VAF given the tumor coverage and the background error rate (forward reads)
      prob_H0_R2 <- dbinom(x = countTR2, size = covTR2, prob = vi_R2, log = TRUE) # probability of observing the tumor VAF given the tumor coverage and the background error rate (reverse reads)
      prob_H1_R1 <- dbinom(x = countTR1, size = covTR1, prob = mu0_R1, log = TRUE) # probability of observing the tumor VAF given the tumor coverage and the estimated VAF in the tumor sample (forward reads)
      prob_H1_R2 <- dbinom(x = countTR2, size = covTR2, prob = mu0_R2, log = TRUE) # probability of observing the tumor VAF given the tumor coverage and the estimated VAF in the tumor sample (reverse reads)

      logBF <- (prob_H1_R1 + prob_H1_R2) - (prob_H0_R1 + prob_H0_R2)
      BF <- exp(logBF)
      print('calculating bayes factor done.')
      print('calculating posterior...')
      posterior <- BF * prior / (BF * prior + (1 - prior)) # based on the Bayes factor get the posterior
      print('calculating posterior done.')

      output <- list("prior" = prior, "posterior" = posterior, "bayes_factor" = BF)
      return(output)
}


call_variant <- function(sample, posterior_cutoff, prior, posterior, bayes_factor, chromosome, start, stop,
                        ref, alt,
                        covTR1, covTR2, countTR1, countTR2, covNR1, covNR2, countNR1, countNR2,
                        vi_R1, vi_R2, mu0_R1, mu0_R2, pseudocount, output_file){
                        
                        print('calling variants...')

                        names <- c('sample','chrom','start','stop','ref','alt',
                                   'decision', 'posterior', 'posterior_cutoff', 'prior', 'bayes_factor',
                                   'coverage_tumor_R1', 'coverage_tumor_R2',
                                   'total_coverage_normals_R1', 'total_coverage_normals_R2',
                                   'alt_count_tumor_R1', 'alt_count_tumor_R2',
                                   'total_alt_count_normals_R1', 'total_alt_count_normals_R2',
                                   'alt_rate_tumor_R1', 'alt_rate_tumor_R2',
                                   'error_rate_normals_R1', 'error_rate_normals_R2',
                                   'pseudocount')


                        if(posterior >= posterior_cutoff){Decision = 'Variant'}
                        else{Decision = 'No variant'}

                        values <- c(sample, chromosome, start, stop, ref, alt, Decision, posterior, posterior_cutoff, prior, bayes_factor,
                                    covTR1, covTR2, covNR1, covNR2, countTR1, countTR2, countNR1, countNR2,
                                    mu0_R1, mu0_R2, vi_R1, vi_R2, pseudocount)
                        
                        output <- data.frame(cbind(names, values))

                        write.csv(output, output_file, row.names = FALSE)

                        return(output)
                        print('done.')

                         }



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


raise_warnings <- function(report_table, mapQ_thresh, baseQ_thresh, output_file_name){

        raised_warnings <- c()
        n_alt <- nrow(report_table)

        strand_bias <- FALSE
        # check if all alternative reads are in one direction
        if(length(unique(report_table$strand)) == 1){
            strand_bias <- TRUE
        }

        low_mapQ_count <- 0
        low_baseQ_count <- 0
        end_of_read_count <- 0

        for(i in 1:nrow(report_table)){

            current_read <- report_table[i,]

            # catch low mapping quality reads
            if(current_read$mapping_quality < mapQ_thresh){low_mapQ_count <- low_mapQ_count + 1}

            # catch low basecalling quality reads
            if(current_read$basecalling_quality < baseQ_thresh){low_baseQ_count <- low_baseQ_count + 1}
            
            # catch bases that are in the last quarter of their respective reads
            if(current_read$position_within_the_read > as.numeric(current_read$read_length) * 0.75){
                end_of_read_count <- end_of_read_count + 1
            }
        }

        warnings <- list(
                        'strand_bias_filter_failed' = strand_bias,
                        'num_reads_failed_mapQ_filter' = low_mapQ_count,
                        'num_reads_failed_baseQ_filter' = low_baseQ_count,
                        'num_reads_failed_late_cycle_filter' = end_of_read_count
                        )

        write.csv(data.frame(warnings), output_file_name, row.names = FALSE)

}


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



#______________________________________________________________________________________________________________

args <- commandArgs(trailingOnly = TRUE)

tumor_bam_path <- as.character(args[1])
matched_normal_bam_path <- as.character(args[2])
negative_control_bamlist <- as.character(args[3])
chromosome <- as.character(args[4])
start <- as.numeric(args[5])
stop <- as.numeric(args[6])
ref <- as.character(args[7])
alt <- as.character(args[8])
prior <- as.numeric(args[9])
posterior_cutoff <- as.numeric(args[10])
pseudocount <- as.numeric(args[11])
output_file <- as.character(args[12])
mapQ_thresh <- as.numeric(args[13])
baseQ_thresh <- as.numeric(args[14])
output_annotation_file <- as.character(args[15])
reference_genome <- as.character(args[16])
padding_upstream <- as.numeric(args[17])
padding_downstream <- as.numeric(args[18])
output_context_warning_file <- as.character(args[19])
raise_normal_posterior_evidence_warning <- as.character(args[20])


# variant calling for the tumor sample______________________________________________________________________________________________________
pileups <- make_pileups(tumor_bam_path = tumor_bam_path,
                        negative_control_bamlist = negative_control_bamlist,
                        CHROM = chromosome,
                        START = start,
                        END = stop,
                        write_files = FALSE)


estimates <- estimate_and_count(TR1 = pileups$target_tumor_pileup_R1,
                                TR2 = pileups$target_tumor_pileup_R2,
                                NR1 = pileups$negative_control_pileup_R1,
                                NR2 = pileups$negative_control_pileup_R2,
                                alternative_allele = alt,
                                write_file = FALSE)

results <- calc_posterior(countTR1 = estimates$countTR1,
                        countTR2 = estimates$countTR2,
                        covTR1 = estimates$covTR1,
                        covTR2 = estimates$covTR2,
                        mu0_R1 = estimates$mu0_R1,
                        mu0_R2 = estimates$mu0_R2,
                        vi_R1 = estimates$vi_R1,
                        vi_R2 = estimates$vi_R2,
                        prior = prior,
                        pseudo = pseudocount)


output <- call_variant(sample = tumor_bam_path,
                       posterior_cutoff = posterior_cutoff,
                       prior = prior,
                       posterior = results$posterior,
                       bayes_factor = results$bayes_factor,
                       chromosome = chromosome,
                       start = start,
                       stop = stop,
                       ref = ref,
                       alt = alt,
                       covTR1 = estimates$covTR1,
                       covTR2 = estimates$covTR2,
                       covNR1 = estimates$covNR1,
                       covNR2 = estimates$covNR2,
                       countTR1 = estimates$countTR1,
                       countTR2 = estimates$countTR2,
                       countNR1 = estimates$countNR1,
                       countNR2 = estimates$countNR2,
                       vi_R1 = estimates$vi_R1,
                       vi_R2 = estimates$vi_R2,
                       mu0_R1 = estimates$mu0_R1,
                       mu0_R2 = estimates$mu0_R2,
                       pseudocount = pseudocount,
                       output_file = output_file)
# ________________________________________________________________________________________________________________________________________________________________________

# 'variant calling' in the matched normal sample__________________________________________________________________________________________________________________________
pileups <- make_pileups(tumor_bam_path = matched_normal_bam_path,
                        negative_control_bamlist = negative_control_bamlist,
                        CHROM = chromosome,
                        START = start,
                        END = stop,
                        write_files = FALSE)


estimates <- estimate_and_count(TR1 = pileups$target_tumor_pileup_R1,
                                TR2 = pileups$target_tumor_pileup_R2,
                                NR1 = pileups$negative_control_pileup_R1,
                                NR2 = pileups$negative_control_pileup_R2,
                                alternative_allele = alt,
                                write_file = FALSE)

results <- calc_posterior(countTR1 = estimates$countTR1,
                        countTR2 = estimates$countTR2,
                        covTR1 = estimates$covTR1,
                        covTR2 = estimates$covTR2,
                        mu0_R1 = estimates$mu0_R1,
                        mu0_R2 = estimates$mu0_R2,
                        vi_R1 = estimates$vi_R1,
                        vi_R2 = estimates$vi_R2,
                        prior = prior,
                        pseudo = pseudocount)
#________________________________________________________________________________________________________________________________________________________________________________

# return the matched normal posterior into the output file______________________________________________________________________________________________________________________
variant_call_report_file <- fread(output_file, sep = ',')
variant_call_report_file <- rbind(variant_call_report_file, data.frame(names = 'posterior_in_matched_normal', values = results$posterior))

if(results$posterior >= raise_normal_posterior_evidence_warning){
    variant_call_report_file <- rbind(variant_call_report_file, data.frame(names = 'normal_posterior_evidence', values = 'warning, normal posterior above threshold'))
}

write.csv(variant_call_report_file, output_file, row.names = FALSE)
#_____________________________________________________________________________________________________________________________________________________________________________


read_annotation <- annotate_reads(
                bam = tumor_bam_path,
                chrom = chromosome,
                pos = start,
                alt = alt)


read_annotation <- raise_warnings(
                                    report_table = read_annotation,
                                    mapQ_thresh = mapQ_thresh,
                                    baseQ_thresh = baseQ_thresh,
                                    output_file_name = output_annotation_file
)


raise_context_warnings(
                            reference_genome = reference_genome,
                            chr = chromosome,
                            alt = alt,
                            pos = start,
                            padding_downstream = padding_downstream,
                            padding_upstream = padding_upstream,
                            output_file_name = output_context_warning_file
)

# _____________________________________________________________________________________________________________________________________

