# How to run
#### 0.) create the environment based on
```save_the_environment.yaml```

#### 1.) define the list of bamfiles
- save the paths to the bamfiles in a .txt file
- one bamfile path per line

#### 2.) define the regions of interest
- they must be stored in the .bed format
- three columns: chromosome, start, stop, tab delimited
- example:
```
                    chr1    100   102
                    chr2    12    132
```

#### 3.) define the config file (in the .yaml) format. It must contain these elements
- the path to the reference genome
- the path to the list of bamfiles
- the path to the regions .bed file

#### 4.) run Snakemake
```
            snakemake --snakefile Mpileup_Snakemake.smk --configfile <your config.yaml> --profile=cubi-v1 --jobs 1
```

- if you get an error you might need to re-run the command with the additional ```--unlock``` flag and then rerun the previous command again.
- ```profile=cubi-v1``` is specific for the BIH HPC

# Where to run in (BIH HPC specific)
Since this snakemake pipeline is being run via slurm there are some steps that need to be done.

#### 1.) on a login node create a dedicated screen
```screen -S snakescreen```
#### 2.) when inside that screen move to a compute node
#### 3.) from within the compute node run the pipeline as shown above
#### 4.) optionally detach from the screen
```ctrl+a d```
#### 5.) optionally reattach to the screen
```screen -r snakescreen```

# Example setup
#### bamfile list
```
path/to/sample.bam
path/to/sample2.bam
```
#### regions bedfile
```
chr1    100   102
chr2    12    132
```
#### example output
```
CHROM,POS,REF,SAMPLE,COVERAGE,READ_BASES,N_ALT_BASES
chr1,100,A,sample1,2,AA,0
chr1,100,A,sample2,2,CA,1
chr2,96,T,sample1,5,TTTTT,0
chr2,96,T,sample2,5,TCCTT,2
```

# Other important things
- the .bed file must be zero based and the end position must be greater than the start position
- for example if we would want to encode the first position that is a 1-basepair mutation we would encode it as ```chr1  0  1```
- moreover the chromosome identifiers must match between the .bed and the .bam files (sometimes chromosome posess the "chr" prefix and sometimes they do not)
- obviously the .bed and .bam files must be based on the same reference genome
