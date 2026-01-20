## General steps 20.01.2026

**1.)** make a pileup of all normal and tumor samples  

**2.)** call the variants in all tumor samples using a standard caller  

**3.)** from the VCF files extract the sample name (without the T/N suffix), the locus and the ref (one extraction per VCF)

```
chr1  123  A BIH_012
chr1  1542  T BIH_012
chr1  1842  G BIH_012
```

**4.)** concatenate all the extracted files into one (paste one above the other) (shape: n_vcfs x 4)  

**5.)** postprocess the pileup  

**6.)** split the postprocessed pileup into the tumors and the normals  

**7.)** remove the T/N suffixes in both files
```
BIH_012-T1 --> BIH_012
BIH_012-N1 --> BUH_012
```

**8.)** in the tumor processed pileup keep only the sites with nonzero coverage and a given minimum number of alternative reads  

**9.)** do an inner merge with the normal processed pileup  

**10.)** do a inner merge with the output of 4.) --> these will be the candidates that are called by the standard caller  

**11.)** take the rows of 9.) that are not in 10.) --> these are the candidates that are not called by the standard caller --> the final candidates
