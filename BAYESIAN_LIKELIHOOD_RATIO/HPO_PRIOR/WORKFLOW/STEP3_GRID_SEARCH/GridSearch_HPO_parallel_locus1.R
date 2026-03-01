source('MakeBayesFactor.R')
source('CallVariants.R')
library(deepSNV)
library(GenomicRanges)
library(VariantAnnotation)
library(Biostrings)
library(data.table)
library(foreach)
library(doParallel)

num_cores <- 70
cl <- makeCluster(num_cores)
registerDoParallel(cl)

F_beta <- function(tp,fp,fn,beta){

    top <- (1+beta^2) * tp
    down <- (1+beta^2)*tp + beta^2*fn + fp

    return(top/down)
}


# load the data
full_normals <- readLines('normal_cohort_bams.txt') # paths to all normal bam files
positives <- readLines('mutated_bamlist_train.txt')   # paths to only the mutated bam files
negatives <- readLines('nonmutated_bamlist_train.txt')  # paths to only the unmutated bam files


# define the actionable variant of interest (and the minimum quality)
chrom <- 7
position <- '140453136'
position_num <- 140453136
padding <- 1
ref <- 'A'
alt <- 'T'
minq <- 10


# make the grid for the grid search
binned_prior_values <- seq(from = 0.1, to = 0.9, by = 0.05)
current_posterior_cutoff <- 0.5
grid <- expand.grid(binned_prior_values,current_posterior_cutoff)  # the first column are the priors and the 2nd column are the posterior cutoffs
output <- data.frame(matrix(-9, nrow = nrow(grid), ncol = 9))
colnames(output) <- c('prior', 'posterior_cutoff', 'F_1', 'F_1.5', 'F_2', 'F_2.5', 'TP', 'FP', 'FN')

#output_file <- file('output2_medium.csv', open = 'w')
#writeLines('PRIOR,POSTERIOR_CUTOFF,F1,F1.5,F2,F2.5,TP,FP,FN,TN', 'output2_medium.csv')



results <- foreach(i=1:nrow(grid),
                    .packages = c('deepSNV', 'GenomicRanges', 'VariantAnnotation', 'Biostrings', 'data.table'),
                    .combine = rbind) %dopar%{


    print(paste0(i/nrow(grid)*100, '%'))

    # fix the parameters
    current_prior <- grid[i,1]
    current_posterior_cutoff <- 0.5

    TP = TN = FP = FN = 0 # reset the counts to 0


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
                                prior_value = current_prior,    # the current prior
                                posterior_cutoff = current_posterior_cutoff,   # the current posterior
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
                                prior_value = current_prior,    # the current prior
                                posterior_cutoff = current_posterior_cutoff,   # the current posterior
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


    f_1 <- F_beta(tp = TP, fp = FP, fn = FN, beta = 1.0)
    f_1_5 <- F_beta(tp = TP, fp = FP, fn = FN, beta = 1.5)
    f_2 <- F_beta(tp = TP, fp = FP, fn = FN, beta = 2.0)
    f_2_5 <- F_beta(tp = TP, fp = FP, fn = FN, beta = 2.5)

    #colnames(results) <- c('prior', 'posterior_cutoff', 'F_1', 'F_1.5', 'F_2', 'F_2.5')
    output$prior[i] <- current_prior
    output$posterior_cutoff[i] <- current_posterior_cutoff
    output$F_1[i] <- f_1
    output$F_1.5[i] <- f_1_5
    output$F_2[i] <- f_2
    output$F_2.5[i] <- f_2_5
    output$TP <- TP
    output$FP <- FP
    output$FN <- FN

    c(current_prior, current_posterior_cutoff, f_1, f_1_5, f_2, f_2_5, TP, FP, FN)


    #output_line <- paste(current_prior, current_posterior_cutoff, f_1, f_1_5, f_2, f_2_5, TP, FP, FN, TN, sep = ",")
    #cat(output_line, "\n", file = 'output2_medium.csv', append = TRUE)

} # end of HP configuration loop

#close(output_file)
#print('done.')


colnames(results) <- c('PRIOR','POSTERIOR_CUTOFF','F1','F1.5','F2','F2.5','TP','FP','FN')

print(results)

write.csv(results, "results_parallel.csv", row.names = FALSE)


stopCluster(cl)

