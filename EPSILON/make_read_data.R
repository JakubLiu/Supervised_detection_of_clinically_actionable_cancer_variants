suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))


create_data <- function(bamlist,    # a list of bam files, one path per line
                        reference_genome,  # the path to the reference genome
                        chrom_list, pos_list, ref_list, alt_list,  # R vectors
                        output,  # path to output file
                        ncores,
                        tmpdir
                        ){
                            cl <- makeCluster(ncores)
                            registerDoParallel(cl)
                            bamlist <- readLines(bamlist)
                            reference_genome <- FaFile(reference_genome)
                            num_loci <- length(chrom_list)
                            chrom_list <- as.character(chrom_list)

                            all_data <- foreach(i = 1:num_loci,
                                    .packages = c('Rsamtools', 'Biostrings', 'stats', 'GenomicAlignments', 'GenomicRanges', 'data.table', 'stats'),
                                                        .combine = rbind) %dopar%{

                                                            current_chromosome <- as.character(chrom_list[i])
                                                            current_pos <- as.numeric(pos_list[i])
                                                            current_ref <- as.character(ref_list[i])
                                                            current_alt <- as.character(alt_list[i])
                                                            region <- GRanges(current_chromosome, IRanges(start = current_pos, end = current_pos))
                                                            scan_bam_param <- ScanBamParam(which = region, what = c("qname", "seq", "qual", "mapq"), tag = 'NM')
                                                            pileup_param <- PileupParam(distinguish_strands = TRUE, distinguish_nucleotides = TRUE)
                                        
                                        
                                        partial_output <- paste0(tmpdir, '/' ,output, "_chr_",current_chromosome,'_pos_',current_pos, Sys.getpid(), ".csv")

                                        read_data_list <- list()
                                        counter <- 1

                                        for(j in seq_along(bamlist)) {

                                            current_bam <- bamlist[j]
                                            bam_basename <- basename(current_bam)

                                           
                                            gal <- readGAlignments(current_bam, param = scan_bam_param)
                                            

                                            if(length(gal) == 0) next

                                            # assign names (needed!)
                                            names(gal) <- mcols(gal)$qname

                                            # map locus to reads
                                            mapped <- mapToAlignments(region, gal)

                                            # subset reads overlapping locus
                                            hits <- findOverlaps(gal, region)
                                            gal_hits <- gal[queryHits(hits)]

                                            # read length
                                            read_lengths <- width(gal_hits)

                                            if(length(gal_hits) == 0) next

                                            names(gal_hits) <- mcols(gal_hits)$qname

                                            pos_in_reads <- start(mapToAlignments(region, gal_hits))  # positions on the locus withi each read

                                            # sequences and qualities
                                            read_seqs <- mcols(gal_hits)$seq

                                            valid <- pos_in_reads >= 1 & pos_in_reads <= read_lengths  # remove those weird places where the position within the read falls outside of the read
                                            pos_in_reads <- pos_in_reads[valid]
                                            gal_hits <- gal_hits[valid]


                                            read_quals <- as.character(mcols(gal_hits)$qual)  # read qualities at every position in the read
                                            mapQ <- mcols(gal_hits)$mapq      # mapping quality of the read
                                            strands <- as.character(strand(gal_hits))   # the strand of the read


                                            # base quality at locus
                                            baseQ_at_pos <- mapply(function(qual, pos)
                                                as.integer(charToRaw(substr(qual, pos, pos))) - 33,
                                                read_quals, pos_in_reads)


                                            # distances
                                            dist_to_start <- pos_in_reads
                                            dist_to_end <- read_lengths - pos_in_reads + 1

                                            # build table
                                            read_dt <- data.table(
                                                chrom = current_chromosome,
                                                pos = current_pos,
                                                sampleID = bam_basename,
                                                baseQ = baseQ_at_pos,
                                                mapQ = mapQ,
                                                strand = strands,
                                                dist_to_read_start = dist_to_start,
                                                dist_to_read_end = dist_to_end
                                            )

                                            read_data_list[[counter]] <- read_dt
                                            counter <- counter + 1

                                            }

                                        fwrite(rbindlist(read_data_list), file = partial_output)  # write partial file
                                        return(rbindlist(read_data_list))
                                        }  # end of the loop over the loci

                                        fwrite(all_data, file = output)

stopCluster(cl)

} # end of the function definition loop