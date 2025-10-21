# What it does
This script takes as input a file that lists paths to gzipped vcf files (one file path per line).
For more biologically stable results it is recomendet to use VCF files where multicallelic calls have been splitted, aka there is only one alternative allele
per row (see ```https://samtools.github.io/bcftools/bcftools.html#norm```)
```
path/to/sample1.vcf.gz
path/to/sample2.vcf.gz
path/to/sample3.vcf.gz
```
It outputs a csv table of dimensions _(n_total_variants x n_vcf_files)_, where each variant is identified by its chromosome, position, reference allele and alternative allele(s).
If a given variant is present in a given VCF file then its presence is denoted by _1_, else its absence is denoted by _0_.

```
variant,file_1,file_2,file_3
chr1;1058324;A;G,1,0,1
chr2;23758192;T;C,1,1,0
chr3;48273910;G;A,0,1,1
chr4;9021845;C;T,1,0,0
chr5;157382912;T;G,1,1,1
chr6;83291872;A;C,0,0,1
chr7;62938420;C;A,1,1,0
chr8;109283745;G;T,0,1,0
chr9;73910283;A;G,1,0,1
chr10;45192384;T;C,1,1,1
chr11;12384921;G;A,0,0,0
chr12;91237485;A;C,1,1,0
chr13;71928403;T;G,1,0,1
chr14;19384722;C;T,0,1,1
chr15;42837491;G;C,1,0,0
chr16;18239475;A;G,1,1,1
chr17;5738291;T;A,0,0,1
chr18;23749385;C;T,1,1,0
chr19;2918473;G;A,1,0,0
chr20;58473821;A;C,0,1,1
```
