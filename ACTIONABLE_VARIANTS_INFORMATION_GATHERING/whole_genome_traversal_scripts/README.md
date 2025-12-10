A collection of scripts that use **MPI** and **Numba** to run on the whole reference genome in a reasonable time (on an HPC with a resource manager).
An example of how to run:
```
mpirun -np <number of CPUs> python3 Calc_Shannon_entropy.optim.py <reference genome fasta> \
                                        <chromosome (with or without 'chr' prefix, depends on the reference)> \
                                        <windowsize> \
                                        <output file name>
```
```
mpirun -np 20 python3 Calc_Shannon_entropy.optim.py ref.fa \
                                        4 \
                                        100000 \
                                        output.csv
```
