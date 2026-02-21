suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stats))

library(Rsamtools)
library(Biostrings)
library(stats)

make_pileups <- function(tumor_bam_path, negative_control_bamlist, CHROM, START, END, write_files = FALSE){

    CHROM <- as.character(CHROM)
    START <- as.numeric(START)
    END <- as.numeric(END)


    actionable_region <- GRanges(CHROM, IRanges(start = START, end = END))
    scan_bam_param <- ScanBamParam(which = actionable_region)
    pileup_param <- PileupParam(distinguish_strands = TRUE, distinguish_nucleotides = TRUE)


    # pileup on the target tumor sample
    target_tumor_pileup <- pileup(file = tumor_bam_path, index = paste0(tumor_bam_path, ".bai"), scanBamParam = scan_bam_param, pileupParam = pileup_param)
    target_tumor_pileup_R1 <- target_tumor_pileup[target_tumor_pileup$strand == '+',]  # the pileup in the target tumor sample, only the forward reads
    target_tumor_pileup_R2 <- target_tumor_pileup[target_tumor_pileup$strand == '-',]  # the pileup in the target tumor sample, only the backward reads

    # pileup in the negative control cohort
    negative_control_list <- readLines(negative_control_bamlist)
    negative_control_pileup_R1 <- rep(-9,6)
    negative_control_pileup_R2 <- rep(-9,6)
    
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

    if(write_files == TRUE){

        current_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

        write.csv(target_tumor_pileup_R1, paste0(current_time, '_T.R1.csv'), row.names = F)
        write.csv(target_tumor_pileup_R2, paste0(current_time, '_T.R2.csv'), row.names = F)
        write.csv(negative_control_pileup_R1, paste0(current_time, '_N.R1.csv'), row.names = F)
        write.csv(negative_control_pileup_R2, paste0(current_time, '_N.R2.csv'), row.names = F)
    }

    list(
        'target_tumor_pileup_R1' = target_tumor_pileup_R1,
        'target_tumor_pileup_R2' = target_tumor_pileup_R2,
        'negative_control_pileup_R1' = negative_control_pileup_R1,
        'negative_control_pileup_R2' = negative_control_pileup_R2
    )

    
}



estimate_and_count <- function(TR1,   # pileup of the tumor sample in the forward direction
                               TR2,   # pileup of the tumor sample in the reverse direction
                               NR1,   # pileup of the negative cohort in the forward direction
                               NR2,   # pileup of the negative cohort in the reverse direction
                               alternative_allele,  # the alternative allele as is in the actionable variant
                               write_file = FALSE  # optionally write the estimates to a file
                               ){

    alternative_allele <- as.character(alternative_allele)
    # the counts _______________________________________________________________________________________________
    # the n_i 's in the paper
    covTR1 <- sum(TR1$count)  # the coverage in the forrward direction in the tumor sample
    covTR2 <- sum(TR2$count)  # the coverage in the reverse direction in the tumor sample
    covNR1 <- sum(NR1$count)  # the total coverage (sum over all samples) in the negative cohort in the forward direction
    covNR2 <- sum(NR2$count)  # the total coverage (sum over all samples) in the negative cohort in the reverse direction

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

    estimates

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

      prob_H0_R1 <- dbinom(x = countTR1, size = covTR1, prob = vi_R1, log = TRUE) # probability of observing the tumor VAF given the tumor coverage and the background error rate (forward reads)
      prob_H0_R2 <- dbinom(x = countTR2, size = covTR2, prob = vi_R2, log = TRUE) # probability of observing the tumor VAF given the tumor coverage and the background error rate (reverse reads)
      prob_H1_R1 <- dbinom(x = countTR1, size = covTR1, prob = mu0_R1, log = TRUE) # probability of observing the tumor VAF given the tumor coverage and the estimated VAF in the tumor sample (forward reads)
      prob_H1_R2 <- dbinom(x = countTR2, size = covTR2, prob = mu0_R2, log = TRUE) # probability of observing the tumor VAF given the tumor coverage and the estimated VAF in the tumor sample (reverse reads)

      logBF <- (prob_H1_R1 + prob_H1_R2) - (prob_H0_R1 + prob_H0_R2)
      BF <- exp(logBF)
      posterior <- BF * prior / (BF * prior + (1 - prior)) # based on the Bayes factor get the posterior

      list("prior" = prior, "posterior" = posterior, "bayes_factor" = BF)
}


call_variant <- function(sample, posterior_cutoff, prior, posterior, bayes_factor, chromosome, start, stop,
                        ref, alt,
                        covTR1, covTR2, countTR1, countTR2, covNR1, covNR2, countNR1, countNR2,
                        vi_R1, vi_R2, mu0_R1, mu0_R2, pseudocount, output_file){
                        
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

                         }


#______________________________________________________________________________________________________________

args <- commandArgs(trailingOnly = TRUE)

tumor_bam_path <- as.character(args[1])
negative_control_bamlist <- as.character(args[2])
chromosome <- as.character(args[3])
start <- as.numeric(args[4])
stop <- as.numeric(args[5])
ref <- as.character(args[6])
alt <- as.character(args[7])
prior <- as.numeric(args[8])
posterior_cutoff <- as.numeric(args[9])
pseudocount <- as.numeric(args[10])
output_file <- as.character(args[11])



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


print(output)


# usage:

# Rscript Trottellumme.R 'tumor.bam' 'bamlist.txt' '7' 100 100 'A' 'G' 0.01 0.5 0.00001 'output.txt'



