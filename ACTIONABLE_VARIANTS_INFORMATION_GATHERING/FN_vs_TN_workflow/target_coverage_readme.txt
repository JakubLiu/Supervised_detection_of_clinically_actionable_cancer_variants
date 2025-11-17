- as far as I know there are no tools that enable the use to subsample to a given coverage
- the tool used here (seqtk) subsamples to a given number of reads
- therefore we must know how much reads we want to keep in order to get our desired subsampled coverage
- this can be done in these steps

    1.) get the size of the reference genome
            e.g. cut -f2 $fai | paste -sd+ - | bc
                    *$fai --> the reference genome index

    2.) get the average read length
            e.g. cat $dat | awk 'NR%4==2 {sum+=length($0)} END {print sum/NR*4}'
                    *$dat --> the .fastq file (use zcat is gzipped)
                    ** another option is to look at the manual page of the sequencing platform

    3.) 
        G = size of the reference genome
        ARL = average read length
        TCOV = desired target coverage
        R = number of reads to keep

        R = (TCOV * G) / ARL