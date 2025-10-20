This Snakemake script first creates a samtools mpileup based on a list of bam files. Then it edits this mpileup to a more readable form.
We need to run this script as a slurm job. Below it is said how to do it.

**1.) define the resources within each rule of the snakefile**

```
rule samtools_mpileup:
    input:
        ref_gen = reference_genome,
        bamlist = bamfile_list
    output:
        'samtools_mpileup/pileup_raw.txt'

    resources:
        mem='60G',         # LOOK HERE !!!
        time='20:00:00'    # LOOK HERE !!!

    shell:
        '''
        samtools mpileup -f {input.ref_gen} -b {input.bamlist} -s > {output}
        '''
```

**2.) Go onto a login node, open/attach to a screen**

```
screen -S snake_screen
```
**3.) Inside the screen run**

```
snakemake --snakefile Mpileup_Snakemake.smk --configfile config.yaml --profile=cubi-v1 --jobs 1
```
Make sure the an appropriate conda evironment is activated!
Make sure that in the new screen you are also on a compute node!
For more information go to: https://hpc-docs.cubi.bihealth.org/slurm/snakemake/

Then you can let Snakemake run in the background and detach from the screen ```ctrl+a d```.
You can reattach to the screen: ```screen -r snake_screen``` or check all available screens: ```screen -ls```.
Do remove a screen you can attach to it and when within the screen type ```exit```.
