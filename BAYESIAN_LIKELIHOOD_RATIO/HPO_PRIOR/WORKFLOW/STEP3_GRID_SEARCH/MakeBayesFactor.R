# This function loads the data, and calculates the Bayes Factor.
# It does not call variants.
# This is irrespective of the prior and should be run once per target (tumor) sample and region.
# This function returns a BayesFactor object and other R objects that will be used by downstream functions

MakeBayesFactor <- function(tumor_bam, # a file path [string]
                            normal_cohort, # a vector of file pathts [strings]
                            chromosome,  # [string]  
                            position,   # [int]
                            padding,  # how many bases upstream and downstream should be extracted aroung the locus [int]
                            min_basecalling_quality,  # minimum basecalling quality, below which bases will be masked as 'N'  [float]
                            model_type = 'OR'  # 'OR' allows for the evidence to come only from one read direction, 'AND' needs the evidence to come from both read directions
                            )
                            
    {

    #library(deepSNV)
    #library(GenomicRanges)
    #library(VariantAnnotation)
    #library(Biostrings)


    chromosome <- as.character(chromosome)
    position <- as.numeric(position)
    padding <- as.numeric(padding)
    min_basecalling_quality <- as.numeric(min_basecalling_quality)

    region <- GRanges(chromosome, IRanges(start = position - padding,
                                end = position + padding))

    data <- deepSNV:::loadAllData(
                                c(tumor_bam, normal_cohort),
                                region,
                                q = min_basecalling_quality)

    BayesFactor <- deepSNV:::bbb(
                                data,
                                model = model_type
    )

    list(
        "BayesFactor" = BayesFactor,
        "data" = data,
        "region" = region
    )
}