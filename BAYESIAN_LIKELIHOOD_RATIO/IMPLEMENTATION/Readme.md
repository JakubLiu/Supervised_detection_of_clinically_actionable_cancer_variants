# How to run

In order to run the model in tumor only mode please run the command below:
```
./LRB_tumor_only.sh \
                --tumor_bam <path to the tumor bam file [1]> \
                --negative_control_bamlist <path to the list of negative control bamfiles [2]> \
                --chromosome <chromosome> \
                --start <start genomic position> \
                --stop <end genomic position (actually for now start = stop)> \
                --ref_allele <the reference allele> \
                --alt_allele <the alternative allele as in the actionable variant> \
                --prior <prior (default 0.005)> \
                --posterior_cutoff <posterior cutoff (default 0.5)> \
                --pseudocount <pseudocount (defaule 0.00001)> \
                --output_call_file <path to the output variant call file> \
                --min_mapQ <minimum mapping quality threshold below which the warning will be returned> \
                --min_baseQ <minimum base calling quality threshold below which the warning will be returned> \
                --output_read_annotation_file <path to the alternative read warning file> \
                --reference_genome <path to the reference genome [3]> \
                --padding_upstream <the number of bases upstream to extract from the reference genome> \
                --padding_downstream <the number of bases downstream to extract from the reference genome> \
                --output_genomic_context_file <path to the genomic context report output file>

[1] the bam file must be sorted and indexed (.bam.bai)
[2] a .txt file with one path per line
      /path/to/sample1.bam
      /path/to/sample2.bam
      /path/to/sample3.bam
[3] the reference genome index (.fasta.fai) must be present in the same directory as the reference genome
```
Example run:
```
./LRB_tumor_only.sh \
                --tumor_bam "T1-DNA1-WES1.mutated.sorted.bam" \
                --negative_control_bamlist "negative_control_cohort.txt" \
                --chromosome "7" \
                --start "55259515" \
                --stop "55259515" \
                --ref_allele "T" \
                --alt_allele "G" \
                --prior "0.005" \
                --posterior_cutoff "0.5" \
                --pseudocount "0.00001" \
                --output_call_file "output_call_file.txt" \
                --min_mapQ "30" \
                --min_baseQ "30" \
                --output_read_annotation_file "output_read_annotation_file.txt" \
                --reference_genome "hs37d5.fa" \
                --padding_upstream "10" \
                --padding_downstream "10" \
                --output_genomic_context_file "output_genomic_context_file.txt"
```

# Output files
### variant call report file

This file return the decision on calling the variant along with the posterior, prior and Bayes factor values.
It also contains informations about the number of alternative reads and the coverage in both forward (R1) and reverse (R2)
directions in the tumor sample and in the negative control cohort.


| names | values |
|-------|--------|
| sample | T1-DNA1-WES1.mutated.sorted.bam |
| chrom | 7 |
| start | 55259515 |
| stop | 55259515 |
| ref | T |
| alt | G |
| decision | Variant |
| posterior | 0.97620021085498 |
| posterior_cutoff | 0.5 |
| prior | 0.005 |
| bayes_factor | 8162.41861540986 |
| coverage_tumor_R1 | 219 |
| coverage_tumor_R2 | 58 |
| total_coverage_normals_R1 | 11613 |
| total_coverage_normals_R2 | 6285 |
| alt_count_tumor_R1 | 2 |
| alt_count_tumor_R2 | 1 |
| total_alt_count_normals_R1 | 2 |
| total_alt_count_normals_R2 | 2 |
| alt_rate_tumor_R1 | 0.0091324200913242 |
| alt_rate_tumor_R2 | 0.0172413793103448 |
| error_rate_normals_R1 | 0.000172220787048997 |
| error_rate_normals_R2 | 0.000318217979315831 |
| pseudocount | 1e-05 |

### read report file

This file contains information about the number of reads that have failed a given filter.

| strand_bias_filter_failed | num_reads_failed_mapQ_filter | num_reads_failed_baseQ_filter | num_reads_failed_late_cycle_filter |
|---------------------------|-----------------------------|-------------------------------|-----------------------------------|
| FALSE                     | 0                           | 0                             | 2                                 |

#### filter description
**strand bias filter**<br>
Is ```TRUE``` when all the alternative reads in the tumor are in one direction, is ```FALSE``` else.<br>
<br>
**num_reads_failed_mapQ_filter**<br>
The number of alternative reads that have a mapping quality below ```min_mapQ```.<br>
<br>
**num_reads_failed_baseQ_filter**<br>
The number of alternative reads in the tumor sample that have a basecalling quality (at the given locus) below ```min_baseQ```.<br>
<br>
**num_reads_failed_late_cycle_filter**<br>
The number of alternative alleles in the tumor sample that are in the last quarter of their respective read.<br>
```if(position_in_read > read_length*0.75){flag}```

### genomic context report file

This file contains information about warnings about the genomic context around the target locus.
These filters are based on the literature [4],[5],[6] and are said to be correlated with higher background error rates in Illumina sequencing platforms.

| GG_motif_upstream_present | GGT_motif_upstream_present | GGC_motif_upstream_present | abnormal_GC_content_warning | preceeding_homopolymer_present |
|---------------------------|---------------------------|---------------------------|----------------------------|-------------------------------|
| TRUE                      | FALSE                     | TRUE                      | FALSE                      | FALSE                         |

[4] Meacham, F., Boffelli, D., Dhahbi, J., Martin, D. I., Singer, M., & Pachter, L. (2011). Identification and correction 
     of systematic error in high-throughput sequence data. BMC bioinformatics, 12, 451. 
     https://doi.org/10.1186/1471-2105-12-451<br>
     
[5] https://en.wikipedia.org/wiki/GC-content#Among-genome_variation<br>

[6] Stoler N, Nekrutenko A. Sequencing error profiles of Illumina sequencing instruments. NAR Genom Bioinform. 
     2021 Mar 27;3(1):lqab019. doi: 10.1093/nargab/lqab019. PMID: 33817639; PMCID: PMC8002175.<br>

**preceeding_homopolymer_present**<br>
checks if a homopolymer of the same base as the alternative, terminates exactly at the alternative<br>
    example of a positive case:<br>
     AAAAAAATC<br>
     123456789<br>
    If the (A) homopolymer ends at position 6 and the base at position 7 is also called as A, then it is more likely to be a sequencing error.
    
