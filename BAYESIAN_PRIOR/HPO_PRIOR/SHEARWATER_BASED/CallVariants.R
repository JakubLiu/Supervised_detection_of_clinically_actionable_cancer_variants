# this function uses the precomputed bayes factor to call variants with a given prior and posterior cutoff
# it return a dataframe with the loci of the called variant(s)


# ========================== usage example ==========================================
#          
#            bbb <- MakeBayesFactor(...)
#
#            CallVariants(
#                bayes_factor = bbb$BayesFactor,
#                prior_value = 0.05,
#                posterior_cutoff = 0.5,
#                output_file_path = 'called_loci.txt',
#                tumor_bam = 'patient1_T1.bam',
#                normal_cohort = c('patient2_N1.bam', ..., 'patientk_N1.bam'),
#                data = bbb$data,
#                region = bbb$region
#            )
# ===================================================================================




CallVariants <- function(bayes_factor,   # an R object returned by the function MakeBayesFactor()
                        prior_value,  # [float]
                        posterior_cutoff,     # [float]
                        output_file_path,
                        tumor_bam,
                        normal_cohort,
                        data,    # R object returned by MakeBayesFactor()
                        region    # R object returned by MakeBayesFactor()
                        )
                        {
                        library(deepSNV)
                        library(GenomicRanges)
                        library(VariantAnnotation)
                        library(Biostrings)

                        vcf <- deepSNV:::bf2Vcf(
                                bayes_factor,
                                data,
                                region,
                                cutoff = posterior_cutoff,
                                sample = c(tumor, normal_cohort),
                                prior = prior_value,
                                mvcf = TRUE
                                )

                        called <- vcf[!is.na(unlist(alt(vcf)))]


                        called_loci <- data.frame(
                        chr = as.character(seqnames(rowRanges(called))),
                        pos = start(rowRanges(called))
                        )

                        return(called_loci)
                        }