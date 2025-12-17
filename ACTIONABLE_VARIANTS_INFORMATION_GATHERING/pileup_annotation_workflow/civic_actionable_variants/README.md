# How to run
## components
All the input files needed to run this pipeline must be supplied in the ```config.yaml``` file. They include:

- the bamfile list
  ```
  /path/to/sample1.bam
  /path/to/sample2.bam
  ```
- the file with actionable variants*
  ```
  17,7577580,7577580,T,C
  7,140481417,140481417,C,A
  17,7577093,7577093,C,A
  17,7578211,7578211,C,G
  3,178916890,178916890,C,T
  7,140453193,140453193,T,C
  17,41228590,41228590,G,A
  13,32890599,32890599,T,G
  ```
  * the columns are as follows: chromosome, start, stop, reference allele, alternative allele(s)
  * the intervals are 1-based and inclusive (so a 1 nucleotide variant can be represented as the interval start:n, stop:n)
  * the file must not have a header!
- a reference genome with all its supporting index files
- the path to the directory where the snpEff.jar file can be found **
- the path to the directory where the Snakefile and the ```scripts/``` folder can be found

** note that the correct annotation file must be present (it must match the reference genome) (```https://sourceforge.net/projects/snpeff/files/databases/```)
** moreover, it might be necessary to manually edit the snpEff.config file when downloading a new annotation file


## how to run
Assuming you are in the ```snakedir/``` directory that looks like this:
```
├── Snakefile.smk
├── config.yaml
├── scripts
   ├── annotation.sh
   ├── csv_to_bed.sh
   ├── csv_to_vcf.py
   ├── pileup.sh
   └── process_pileup.py
```
run:
```
snakemakake --snakefile Snakefile.smk --configfile config.yaml --cores <number of cores>
```

## output
After the pipeline completed correcrtly, the ```snakedir/``` should look like this:
```
├── Snakefile.smk
├── actionables.bed
├── actionables.vcf
├── annotation
│   └── actionables.annotated.vcf
├── config.yaml
├── pileup
│   └── pileup.txt
├── processed_pileup
│   └── pileup.processed.csv
├── scripts
│   ├── annotation.sh
│   ├── csv_to_bed.sh
│   ├── csv_to_vcf.py
│   ├── pileup.sh
│   └── process_pileup.py
├── snpEff_genes.txt
├── snpEff_summary.html
```
The annotated actionable variants should look like this:
```
##fileformat=VCFv4.2
##SnpEffVersion="4.3t (build 2017-11-24 10:18), by Pablo Cingolani"
##SnpEffCmd="SnpEff  GRCh37.87 actionables.vcf "
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
17	7577580	.	T	C	.	.	ANN=C|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000269305.4|protein_coding|7/11|c.701A>G|p.Tyr234Cys|891/2579|701/1182|234/393||,C|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000413465.2|protein_coding|6/7|c.701A>G|p.Tyr234Cys|701/1018|701/858|234/285||,C|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000359597.4|protein_coding|6/9|c.701A>G|p.Tyr234Cys|701/1152|701/1032|234/343||,C|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000420246.2|protein_coding|7/12|c.701A>G|p.Tyr234Cys|834/2653|701/1026|234/341||,C|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000455263.2|protein_coding|7/12|c.701A>G|p.Tyr234Cys|834/2580|701/1041|234/346||,C|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000445888.2|protein_coding|7/11|c.701A>G|p.Tyr234Cys|837/2506|701/1182|234/393||,C|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000509690.1|protein_coding|4/6|c.305A>G|p.Tyr102Cys|437/729|305/597|102/198||WARNING_TRANSCRIPT_NO_STOP_CODON,C|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000514944.1|protein_coding|6/6|c.422A>G|p.Tyr141Cys|501/546|422/467|141/154||WARNING_TRANSCRIPT_INCOMPLETE,C|upstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000576024.1|protein_coding||c.-677A>G|||||675|WARNING_TRANSCRIPT_NO_START_CODON,C|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000505014.1|retained_intron||n.*264A>G|||||264|,C|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000508793.1|protein_coding||c.*854A>G|||||854|WARNING_TRANSCRIPT_INCOMPLETE,C|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000604348.1|protein_coding||c.*900A>G|||||900|WARNING_TRANSCRIPT_NO_STOP_CODON,C|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000503591.1|protein_coding||c.*967A>G|||||967|WARNING_TRANSCRIPT_INCOMPLETE,C|non_coding_transcript_exon_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000504290.1|retained_intron|3/8|n.583A>G||||||,C|non_coding_transcript_exon_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000504937.1|retained_intron|3/7|n.583A>G||||||,C|non_coding_transcript_exon_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000510385.1|retained_intron|3/8|n.583A>G||||||,C|non_coding_transcript_exon_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000574684.1|processed_transcript|2/2|n.96A>G||||||
7	140481417	.	C	A	.	.	ANN=A|missense_variant|MODERATE|BRAF|ENSG00000157764|transcript|ENST00000288602.6|protein_coding|11/18|c.1391G>T|p.Gly464Val|1452/2480|1391/2301|464/766||,A|missense_variant|MODERATE|BRAF|ENSG00000157764|transcript|ENST00000496384.2|protein_coding|2/10|c.212G>T|p.Gly71Val|214/8294|212/1125|71/374||WARNING_TRANSCRIPT_NO_START_CODON,A|3_prime_UTR_variant|MODIFIER|BRAF|ENSG00000157764|transcript|ENST00000497784.1|nonsense_mediated_decay|12/19|c.*841G>T|||||26678|WARNING_TRANSCRIPT_NO_START_CODON
17	7577093	.	C	A	.	.	ANN=A|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000269305.4|protein_coding|8/11|c.845G>T|p.Arg282Leu|1035/2579|845/1182|282/393||,A|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000359597.4|protein_coding|7/9|c.845G>T|p.Arg282Leu|845/1152|845/1032|282/343||,A|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000420246.2|protein_coding|8/12|c.845G>T|p.Arg282Leu|978/2653|845/1026|282/341||,A|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000455263.2|protein_coding|8/12|c.845G>T|p.Arg282Leu|978/2580|845/1041|282/346||,A|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000445888.2|protein_coding|8/11|c.845G>T|p.Arg282Leu|981/2506|845/1182|282/393||,A|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000509690.1|protein_coding|5/6|c.449G>T|p.Arg150Leu|581/729|449/597|150/198||WARNING_TRANSCRIPT_NO_STOP_CODON,A|upstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000576024.1|protein_coding||c.-190G>T|||||188|WARNING_TRANSCRIPT_NO_START_CODON,A|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000514944.1|protein_coding||c.*442G>T|||||442|WARNING_TRANSCRIPT_INCOMPLETE,A|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000574684.1|processed_transcript||n.*479G>T|||||479|,A|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000505014.1|retained_intron||n.*751G>T|||||751|,A|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000508793.1|protein_coding||c.*1341G>T|||||1341|WARNING_TRANSCRIPT_INCOMPLETE,A|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000604348.1|protein_coding||c.*1387G>T|||||1387|WARNING_TRANSCRIPT_NO_STOP_CODON,A|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000503591.1|protein_coding||c.*1454G>T|||||1454|WARNING_TRANSCRIPT_INCOMPLETE,A|intron_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000413465.2|protein_coding|6/6|c.782+406G>T||||||,A|non_coding_transcript_exon_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000504290.1|retained_intron|4/8|n.727G>T||||||,A|non_coding_transcript_exon_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000504937.1|retained_intron|4/7|n.727G>T||||||,A|non_coding_transcript_exon_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000510385.1|retained_intron|4/8|n.727G>T||||||
17	7578211	.	C	G	.	.	ANN=G|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000269305.4|protein_coding|6/11|c.638G>C|p.Arg213Pro|828/2579|638/1182|213/393||,G|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000413465.2|protein_coding|5/7|c.638G>C|p.Arg213Pro|638/1018|638/858|213/285||,G|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000359597.4|protein_coding|5/9|c.638G>C|p.Arg213Pro|638/1152|638/1032|213/343||,G|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000420246.2|protein_coding|6/12|c.638G>C|p.Arg213Pro|771/2653|638/1026|213/341||,G|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000455263.2|protein_coding|6/12|c.638G>C|p.Arg213Pro|771/2580|638/1041|213/346||,G|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000445888.2|protein_coding|6/11|c.638G>C|p.Arg213Pro|774/2506|638/1182|213/393||,G|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000509690.1|protein_coding|3/6|c.242G>C|p.Arg81Pro|374/729|242/597|81/198||WARNING_TRANSCRIPT_NO_STOP_CODON,G|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000514944.1|protein_coding|5/6|c.359G>C|p.Arg120Pro|438/546|359/467|120/154||WARNING_TRANSCRIPT_INCOMPLETE,G|upstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000576024.1|protein_coding||c.-1308G>C|||||1306|WARNING_TRANSCRIPT_NO_START_CODON,G|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000508793.1|protein_coding||c.*223G>C|||||223|WARNING_TRANSCRIPT_INCOMPLETE,G|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000604348.1|protein_coding||c.*269G>C|||||269|WARNING_TRANSCRIPT_NO_STOP_CODON,G|downstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000503591.1|protein_coding||c.*336G>C|||||336|WARNING_TRANSCRIPT_INCOMPLETE,G|intron_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000574684.1|processed_transcript|1/1|n.67+160G>C||||||,G|non_coding_transcript_exon_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000504290.1|retained_intron|2/8|n.520G>C||||||,G|non_coding_transcript_exon_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000504937.1|retained_intron|2/7|n.520G>C||||||,G|non_coding_transcript_exon_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000510385.1|retained_intron|2/8|n.520G>C||||||,G|non_coding_transcript_exon_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000505014.1|retained_intron|5/5|n.894G>C||||||
3	178916890	.	C	T	.	.	ANN=T|missense_variant|MODERATE|PIK3CA|ENSG00000121879|transcript|ENST00000263967.3|protein_coding|2/21|c.277C>T|p.Arg93Trp|434/9093|277/3207|93/1068||,T|missense_variant|MODERATE|PIK3CA|ENSG00000121879|transcript|ENST00000468036.1|protein_coding|2/2|c.277C>T|p.Arg93Trp|543/620|277/354|93/117||WARNING_TRANSCRIPT_NO_STOP_CODON,T|downstream_gene_variant|MODIFIER|PIK3CA|ENSG00000121879|transcript|ENST00000477735.1|protein_coding||c.*214C>T|||||214|WARNING_TRANSCRIPT_NO_STOP_CODON
7	140453193	.	T	C	.	.	ANN=C|missense_variant&splice_region_variant|MODERATE|BRAF|ENSG00000157764|transcript|ENST00000288602.6|protein_coding|15/18|c.1742A>G|p.Asn581Ser|1803/2480|1742/2301|581/766||,C|missense_variant&splice_region_variant|MODERATE|BRAF|ENSG00000157764|transcript|ENST00000496384.2|protein_coding|6/10|c.563A>G|p.Asn188Ser|565/8294|563/1125|188/374||WARNING_TRANSCRIPT_NO_START_CODON,C|missense_variant&splice_region_variant|MODERATE|BRAF|ENSG00000157764|transcript|ENST00000479537.1|nonsense_mediated_decay|2/6|c.26A>G|p.Asn9Ser|26/743|26/309|9/102||WARNING_TRANSCRIPT_NO_START_CODON,C|splice_region_variant|LOW|BRAF|ENSG00000157764|transcript|ENST00000497784.1|nonsense_mediated_decay|16/19|c.*1192A>G||||||WARNING_TRANSCRIPT_NO_START_CODON,C|3_prime_UTR_variant|MODIFIER|BRAF|ENSG00000157764|transcript|ENST00000497784.1|nonsense_mediated_decay|16/19|c.*1192A>G|||||54902|WARNING_TRANSCRIPT_NO_START_CODON
```
The processed pileup should have the following structure:
| **SAMPLE** | **CHROMOSOME** | **START** | **STOP** | **REF** | **CIVIC_ALT** | **COVERAGE**| **N_TOTAL_ALT_READS** | **N_CIVIC_ALT_READS** | **READS** |
|------|------|------|------|------|------|------|------|------|------|
| sample1 | 1 | 100 | 100 | A | C | 100 | 2 | 1 | AAAAAATAAAAAAAACAA |
| sample2 | 1 | 100 | 100 | A | C | 83 | 0 | 0 | AAAAAAAAAAAAAAAAAA |
| sample1 | 2 | 132456 | 132456 | G | 32 | T | 3 | 3 | GGGGGGTTTGGG |


