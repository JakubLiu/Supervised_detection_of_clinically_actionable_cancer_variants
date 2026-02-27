library(data.table)

raise_warnings <- function(report_table, mapQ_thresh, baseQ_thresh, output_file_name){

        raised_warnings <- c()
        n_alt <- nrow(report_table)

        strand_bias <- FALSE
        # check if all alternative reads are in one direction
        if(length(unique(report_table$strand)) == 1){
            strand_bias <- TRUE
        }

        low_mapQ_count <- 0
        low_baseQ_count <- 0
        end_of_read_count <- 0

        for(i in 1:nrow(report_table)){

            current_read <- report_table[i,]

            # catch low mapping quality reads
            if(current_read$mapping_quality < mapQ_thresh){low_mapQ_count <- low_mapQ_count + 1}

            # catch low basecalling quality reads
            if(current_read$basecalling_quality < baseQ_thresh){low_baseQ_count <- low_baseQ_count + 1}
            
            # catch bases that are in the last quarter of their respective reads
            if(current_read$position_within_the_read > current_read$read_length * 0.75){
                end_of_read_count <- end_of_read_count + 1
            }
        }

        warnings <- list(
                        'strand_bias_filter_failed' = strand_bias,
                        'num_reads_failed_mapQ_filter' = low_mapQ_count,
                        'num_reads_failed_baseQ_filter' = low_baseQ_count,
                        'num_reads_failed_late_cycle_filter' = end_of_read_count
                        )

        write.csv(data.frame(warnings), output_file_name, row.names = FALSE)


}