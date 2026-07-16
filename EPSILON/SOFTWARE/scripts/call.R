library(data.table)
library(MASS)
library(glmmTMB)

# this is for later, when the output is a VCF, for now I will just print the variant to stdout
output_dir <- "variant_calls"

if (dir.exists(output_dir)) {
    # Remove previous fitted models and scaling parameters
    old_files <- list.files(output_dir, full.names = TRUE)
    if (length(old_files) > 0) {
        file.remove(old_files)
    }
} else {
    dir.create(output_dir, showWarnings = FALSE)
}


# ==========================================================================================================================================================
args <- commandArgs(trailingOnly = TRUE)

# Required: fitted model RDS
if (length(args) >= 1) {
    fitted_model_rds <- args[1]
} else {
    stop("Missing required argument: fitted_model_rds")
}

# Required: input data
if (length(args) >= 2) {
    input_data <- args[2]
} else {
    stop("Missing required argument: input_data")
}

# Required: output VCF
if (length(args) >= 3) {
    output_vcf <- args[3]
} else {
    stop("Missing required argument: output_vcf")
}

# Optional: variant calling mode
if (length(args) >= 4) {
    variant_calling_mode <- args[4]
} else {
    variant_calling_mode <- "binomial_test"
}

if (variant_calling_mode == "binomial_test") {

    if (length(args) >= 5) {
        multiple_testing_correction_method <- args[5]
    } else {
        multiple_testing_correction_method <- "BH"
    }

    if (length(args) >= 6) {
        significance_level_alpha <- as.numeric(args[6])
    } else {
        significance_level_alpha <- 0.05
    }

} else if (variant_calling_mode == "bayes_posterior") {

    if (length(args) >= 5) {
        prior <- as.numeric(args[5])
    } else {
        prior <- 0.3
    }

    if (length(args) >= 6) {
        posterior_cutoff <- as.numeric(args[6])
    } else {
        posterior_cutoff <- 0.5
    }

} else {
    stop("variant_calling_mode must be 'binomial_test' or 'bayes_posterior'")
}

# Optional: alt mode
if (length(args) >= 7) {
    alt_mode <- args[7]
} else {
    alt_mode <- "generic"
}

if (!(alt_mode %in% c("specific", "generic"))) {
    stop("alt_mode must be 'specific' or 'generic'")
}

# =======================================================================================================================================================

model <- readRDS(fitted_model_rds)
data <- fread(input_data, head = TRUE, sep = ',')

# optionally scale the predictors
script_path <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])
script_dir <- dirname(normalizePath(script_path))
project_dir <- normalizePath(file.path(script_dir, ".."))

scaling_file <- file.path(project_dir, "fitted_error_model", "scaling_params.rds")

if (file.exists(scaling_file)) {
    
  scaling_params <- readRDS(scaling_file)
  data$GC_content <- (data$GC_content - scaling_params$GC_content$center) / scaling_params$GC_content$scale
  data$mismatches <- (data$mismatches - scaling_params$mismatches$center) / scaling_params$mismatches$scale
  data$median_mapQ <- (data$median_mapQ - scaling_params$median_mapQ$center) / scaling_params$median_mapQ$scale
  data$median_baseQ <- (data$median_baseQ - scaling_params$median_baseQ$center) / scaling_params$median_baseQ$scale
}



if(alt_mode == 'specific'){
    observed_alt_count <- data$specific_alt_counts
    alt <- data$alt
}else{
    observed_alt_count <- data$general_alt_counts
    alt <- rep('.', length(observed_alt_count))
}


EPS_PRED <- predict(model, newdata=data, allow.new.levels=TRUE, type='response')  # the predicted error rate
observed_coverage <- data$coverage
observed_vaf <- observed_alt_count/observed_coverage


if(variant_calling_mode == 'binomial_test'){
    PVALUES <- pbinom(
    observed_alt_count - 1,
    size = observed_coverage,
    prob = EPS_PRED,
    lower.tail = FALSE)

    PVALUES_ADJUSTED <- p.adjust(PVALUES, method = multiple_testing_correction_method, n = length(PVALUES))
    significant_idx <- which(PVALUES_ADJUSTED < significance_level_alpha)
    sig_data <- data[significant_idx,]
    
    header <- c(
    "##fileformat=VCFv4.3",
    "##source=Epsilon Binomial test mode",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    '##INFO=<ID=AD,Number=1,Type=Integer,Description="Alternate Allele Count">',
    '##INFO=<ID=AF,Number=1,Type=Float,Description="Observed Variant Allele Fraction">',
    "#CHROM\tPOS\tID\tREF\tALT\tPVAL_ADJ\tPVAL_RAW\tSIG_LEVEL\tINFO"
    )
    
    writeLines(header, output_vcf)


    vcf <- data.table(
    CHROM = sig_data$chrom,
    POS = sig_data$pos,
    ID = sig_data$sampleID,
    REF = sig_data$ref,
    ALT = alt[significant_idx],
    PVAL_ADJ = PVALUES_ADJUSTED[significant_idx],
    PVAL_RAW = PVALUES[significant_idx],
    SIG_LEVEL = rep(significance_level_alpha, length(significant_idx)),
    INFO = sprintf(
        "DP=%d;AD=%d;AF=%.6f",
        sig_data$coverage,
        observed_alt_count[significant_idx],
        observed_alt_count[significant_idx] / sig_data$coverage
     )
    )

    fwrite(
        vcf,
        file = output_vcf,
        append = TRUE,
        sep = "\t",
        quote = FALSE,
        col.names = FALSE
    )

    print(paste0('Output vcf written to: ', output_vcf))

}else if (variant_calling_mode == 'bayes_posterior'){

    P_H0 <- dbinom(x=observed_alt_count, size=observed_coverage, prob=EPS_PRED, log=TRUE)
    P_H1 <- dbinom(x=observed_alt_count, size=observed_coverage, prob=observed_vaf, log=TRUE)
    BF <- exp(P_H1 - P_H0)
    POSTERIOR <- BF * prior / (BF * prior + (1 - prior))
    significant_idx <- which(POSTERIOR > posterior_cutoff)
    sig_data <- data[significant_idx,]

    header <- c(
    "##fileformat=VCFv4.3",
    "##source=Epsilon Bayesian posterior mode",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    '##INFO=<ID=AD,Number=1,Type=Integer,Description="Alternate Allele Count">',
    '##INFO=<ID=AF,Number=1,Type=Float,Description="Observed Variant Allele Fraction">',
    "#CHROM\tPOS\tID\tREF\tALT\tPRIOR\tPOSTERIOR\tPOSTERIOR_CUTOFF\tINFO"
    )

    writeLines(header, output_vcf)

    vcf <- data.table(
    CHROM = sig_data$chrom,
    POS = sig_data$pos,
    ID = sig_data$sampleID,
    REF = sig_data$ref,
    ALT = alt[significant_idx],
    PRIOR = rep(prior, length(significant_idx)),
    POSTERIOR = POSTERIOR[significant_idx],
    POSTERIOR_CUTOFF = rep(posterior_cutoff, length(significant_idx)),
    INFO = sprintf(
        "DP=%d;AD=%d;AF=%.6f",
        sig_data$coverage,
        observed_alt_count[significant_idx],
        observed_alt_count[significant_idx] / sig_data$coverage
     )
    )

    fwrite(
        vcf,
        file = output_vcf,
        append = TRUE,
        sep = "\t",
        quote = FALSE,
        col.names = FALSE
    )

    print(paste0('Output vcf written to: ', output_vcf))




}


