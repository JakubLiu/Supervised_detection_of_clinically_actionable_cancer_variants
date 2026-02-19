# This function takes as input a genomic region, a path to the target tumor bam file and a .txt file of paths to the bam files that are
# in the negative cohort (one path per line). It outputs 4 pileups. 2 for the target tumor sample (one in both read directions)
# and 2 for the negative control cohort (one in both directions).


# * I intetionally do not do any bam filtering here, because I want it to be done before entering R with the bams.

library(Rsamtools)
library(Biostrings)

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

# _____________________________________________________________________________________________________________________________________________________________________________

# Usage:

#make_pileups(
#                tumor_bam_path = '/data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/TUNE_PRIOR/STEP3_SIMULATE_SNV/TRAIN/MUTATED/inserted_SNV/bwa.BIH_183-T1-DNA1-WES1.mutated.sorted.bam',
#                negative_control_bamlist = '/data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/TUNE_PRIOR/STEP3_SIMULATE_SNV/all_normals.txt',
#                CHROM = '7',
#                START = 55259515,
#                END = 55259515,
#                write_files = TRUE
#)
