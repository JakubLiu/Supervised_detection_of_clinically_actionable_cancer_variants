source('MakeBayesFactor.R')
source('CallVariants.R')
library(data.table)



# make the grid for the grid search
binned_prior_values <- seq(from = 0.2, to = 0.9, by = 0.1)
binned_posterior_cutoff_values <- seq(from = 0.1, to = 0.7, by = 0.1)
grid <- expand.grid(binned_prior_values,binned_posterior_cutoff_values)  # the first column are the priors and the 2nd column are the posterior cutoffs









# read in the labeled training data
false_negatives <- fread('LABELED_DATA/false_negatives.txt', header = FALSE, sep = ';')
true_negatives <- fread('LABELED_DATA/true_negatives.txt', header = FALSE, sep = ';')
colnames(false_negatives) <- c('chr', 'pos', 'ref')
colnames(true_negatives) <- colnames(false_negatives)









# read in the bam files and calculate the Bayes factor (this has to be done only once and stays constant for all hyperparameter combinations)
tumor_bam <- 'path to the tumor bam'

normal_cohort <- 'a vector of paths to the normal samples'
chrom <- 10
position <- '27966664'
padding <- 100
minq <- 10

results <- MakeBayesFactor(tumor_bam = tumor_bam,
                           normal_cohort = normal_cohort,
                           chromosome = chrom,
                           position = position,
                           padding = padding,
                           min_basecalling_quality = minq)


bayes_factor <- results$BayesFactor
data <- results$data
region <- results$region









# loop over the hyperparameter combinations
#       I know that a for loop is inefficient, but it is good for development

pseudocount <- 0.0000000001
n_combinations <- nrow(grid)
results <- matrix(-9, nrow = n_combinations, ncol = 3)

for(i in 1:n_combinations){

    current_prior <- grid[i,1]
    current_posterior_cutoff <- grid[i,2]

    called_loci <- CallVariants(
                                bayes_factor = bayes_factor,
                                prior_value = current_prior,    # the current prior
                                posterior_cutoff = current_posterior_cutoff,   # the current posterior
                                output_file_path = FALSE,
                                tumor_bam = tumor_bam,
                                normal_cohort = normal_cohort,
                                data = data,
                                region = region
                                )


    called_loci$alt <- NULL

    # number of the called variants that are labelled as false negatives, normalized by ...
    # ... the total number of called variants
    score <- nrow(merge(called_loci, false_negatives)) / nrow(called_loci) + pseudocount

    results[i,1] <- current_prior
    results[i,2] <- current_posterior_cutoff
    results[i,3] <- score



}

# find the combination of hyperparameters that corresponds to the best score
best_hyperparam_combination <- results[which.max(results[, 3]), ]

print(best_hyperparam_combination)

