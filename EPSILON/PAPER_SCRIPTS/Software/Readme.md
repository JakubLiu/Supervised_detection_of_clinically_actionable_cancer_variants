## Epsilon_MakeData

### prerequisites
- the bam files must be sorted and indexed (```samtools sort, samtools index```)
- the reference genome must be indexed (all necessary index files produced by ```bwa``` and ```samtools``` must be present)

### input / how to run
```
./Epsilon_MakeData.sh \
    --bamlist bamlist.txt \
    --loci_list loci_list_specific.minimal.txt\
    --reference_genome hs37d5.fa \
    --alt_mode specific \
    --nranks 2 \
    --output_prefix data
```
- ```bamlist```
    - a text file with paths to bam files, one path per line
    - it is assumed that the .bai index files are located in the same directory as the bam files
- ```loci_list```
    - a bed file that must match **exactly** one of these two structures:
      
        1.) with no specified alternative alleles
          
         | chrom | pos       | ref |
         |-------|-----------|-----|
         | 9     | 133748283 | C   |
         | 9     | 133738363 | G   |
          
        2.) with specified alternative alleles
         | chrom | pos       | ref | alt |
         |-------|-----------|-----|-----|
         | 9     | 133748283 | C   | T   |
         | 9     | 133738363 | G   | A   |
- ```reference_genome```
      - the chromosome names from the reference genome must match the names found in the ```loci_list``` files (e.g. chr1 = chr1 or 1 = 1)
   
- ```alt_mode```
     - ```specific``` must be used only in combination with the ```loci_list``` from 2.)
     - ```generic``` must be used only in combination with the ```loci_list``` from 1.)
     - In general, if we want to call actionable variants with a defined alternative allele, the ```specific``` mode should be used for increased sensitivity
   
- ```nranks```
     - the number of jobs (mpi4py)

- ```output_prefix```
     - the prefix for the output files

### output
The results will be located in the ```negative_control_data/``` directory.
