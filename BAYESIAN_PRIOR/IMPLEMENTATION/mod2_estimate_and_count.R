

estimate_and_count <- function(TR1,   # pileup of the tumor sample in the forward direction
                               TR2,   # pileup of the tumor sample in the reverse direction
                               NR1,   # pileup of the negative cohort in the forward direction
                               NR2,   # pileup of the negative cohort in the reverse direction
                               alternative_allele,  # the alternative allele as is in the actionable variant
                               tumor_weight = 1.0,  # an optional weight added to the count of the alternative allele in the tumor sample
                               rho = 0.0001,  # the dispersion parameter
                               write_file = FALSE  # optionally write the estimates to a file
                               ){

    alternative_allele <- as.character(alternative_allele)
    tumor_weight <- as.numeric(tumor_weight)
    rho <- as.numeric(rho)

    # the counts _______________________________________________________________________________________________
    # the n_i 's in the paper
    covTR1 <- sum(TR1$count)  # the coverage in the forrward direction in the tumor sample
    covTR2 <- sum(TR2$count)  # the coverage in the reverse direction in the tumor sample
    covNR1 <- sum(NR1$count)  # the total coverage (sum over all samples) in the negative cohort in the forward direction
    covNR2 <- sum(NR2$count)  # the total coverage (sum over all samples) in the negative cohort in the reverse direction

    # the X_i 's in the paper
    countTR1 <- sum(TR1$count[TR1$nucleotide == alternative_allele]) * tumor_weight  # weighted count of the alternative allele in the tumor in the forward direction
    countTR2 <- sum(TR2$count[TR2$nucleotide == alternative_allele]) * tumor_weight  # weighted count of the alternative allele in the tumor in the reverse direction
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
                        'vi_R1' = vi_R1, 'vi_R2' = vi_R2, 'v0_R1' = v0_R1, 'v0_R2' = v0_R2,
                        'mu0_R1' = mu0_R1, 'mu0_R2' = mu0_R2, 'mu_i' = mu_i,
                        'rho' = rho
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

# ____________________________________________________________________________________________________________________________________________________________________________________
# usage:

#source('mod1_make_pileups.R')

#pileups <- make_pileups(
#                tumor_bam_path = 'T1-DNA1-WES1.mutated.sorted.bam',
#                negative_control_bamlist = '/all_normals.txt',
#                CHROM = '7',
#                START = 55259515,
#                END = 55259515,
#                write_files = FALSE
#                    )


#estimates <- estimate_and_count(TR1 = pileups$target_tumor_pileup_R1,
#                                TR2 = pileups$target_tumor_pileup_R2,
#                                NR1 = pileups$negative_control_pileup_R1,
#                                NR2 = pileups$negative_control_pileup_R2,
#                                alternative_allele = 'G',
#                                tumor_weight = 1.0,
#                                rho = 0.0001,
#                                write_file = FALSE)


#print(estimates)
