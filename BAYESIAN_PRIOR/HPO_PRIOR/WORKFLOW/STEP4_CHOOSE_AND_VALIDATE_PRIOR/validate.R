source('MakeBayesFactor.R')
source('CallVariants.R')
library(deepSNV)
library(GenomicRanges)
library(VariantAnnotation)
library(Biostrings)
library(data.table)



F_beta <- function(tp,fp,fn,beta){

    top <- (1+beta^2) * tp
    down <- (1+beta^2)*tp + beta^2*fn + fp

    return(top/down)
}


# load the data
full_normals <- readLines('/normal_cohort_bams.txt') # paths to all normal bam files
positives <- readLines('mutated_bams_validate.txt')   # paths to only the mutated bam files
negatives <- readLines('nonmutated_bams_validate.txt')  # paths to only the unmutated bam files


# define the actionable variant of interest (and the minimum quality)
chrom <- 7
position <- '55259515'
position_num <- 55259515
padding <- 10
ref <- 'T'
alt <- 'G'
minq <- 10


prior <- 0.5
posterior_cutoff <- 0.5


TP = TN = FP = FN = 0 # set the counts to 0


for(mut in 1:length(positives)){ # loop over the samples with SNV's

        current_target_sample <- positives[mut] # fix the target sample

        results <- MakeBayesFactor(tumor_bam = current_target_sample,
                           normal_cohort = full_normals,
                           chromosome = chrom,
                           position = position,
                           padding = padding,
                           min_basecalling_quality = minq)

        bayes_factor <- results$BayesFactor
        data <- results$data
        region <- results$region

        called_loci <- CallVariants(
                                bayes_factor = bayes_factor,
                                prior_value = prior,    # the current prior
                                posterior_cutoff = posterior_cutoff,   # the current posterior
                                output_file_path = FALSE,
                                tumor_bam = current_target_sample,
                                normal_cohort = full_normals,
                                data = data,
                                region = region
                                )

        # a boolean if the given actionable variant is called or not
        actionable_called <- any(called_loci$chr == chrom &
                                    called_loci$pos == position_num &
                                    called_loci$ref == ref &
                                    called_loci$alt == alt)

        
        if(actionable_called == TRUE){
            TP = TP + 1
        }
        else{
            FN = FN + 1
        }

} # end of mutated SNV loop


for(nonmut in 1:length(negatives)){   # loop over the samples without SNV's

        current_target_sample <- negatives[nonmut]  # fix the target sample

        results <- MakeBayesFactor(tumor_bam = current_target_sample,
                           normal_cohort = full_normals,
                           chromosome = chrom,
                           position = position,
                           padding = padding,
                           min_basecalling_quality = minq)

        bayes_factor <- results$BayesFactor
        data <- results$data
        region <- results$region

        called_loci <- CallVariants(
                                bayes_factor = bayes_factor,
                                prior_value = prior,    # the current prior
                                posterior_cutoff = posterior_cutoff,   # the current posterior
                                output_file_path = FALSE,
                                tumor_bam = current_target_sample,
                                normal_cohort = full_normals,
                                data = data,
                                region = region
                                )

        # a boolean if the given actionable variant is called or not
        actionable_called <- any(called_loci$chr == chrom &
                                    called_loci$pos == position_num &
                                    called_loci$ref == ref &
                                    called_loci$alt == alt)

        
        if(actionable_called == TRUE){
            FP = FP + 1
        }
        else{
            TN = TN + 1
        }

    }  # end of nonmutated SNV loop


print(c(TP,TN,FP,FN))
print(paste0('F1 score: ', F_beta(tp = TP, fp = FP, fn = FN, beta = 1)))
print(paste0('F2 score: ', F_beta(tp = TP, fp = FP, fn = FN, beta = 2)))
