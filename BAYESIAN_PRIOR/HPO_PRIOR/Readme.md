# Tuning the prior
## Theory
Inspired by [1], local ow VAF (variant allele fraction) single nucleotide variants can be called by learning the local background error rate on a cohort on negative control samples.
Then the expected number of alternative reads (under the null hypothesis of no variant) is compared to the observed number of alternative reads. The authors propose a Bayesian model.
The aim of this workflow is to tune the prior of the Bayesian model in order to increase local sensitivity.

## Steps
### 1.) filtering the bam files
- removing duplicates
- removing reads based on median basecalling quality
- removing reads based on mapping quality

### 2.) simulating SNVs
Single nucleotide variants are simulated at a defined locus and a defined VAF using [2].

### 3.) define the datasets
```
data
  |
  |_train (the prior will be tuned on these samples)____________
  |                                                             |_ mutated bam files (mutations inserted with bam surgeon [2])
  |                                                             |_ non mutated bam files
  |
  |
  |_validate (the tuned prior will be validated on these samples)____________
  |                                                                          |_ mutated bam files (mutations inserted with bam surgeon [2])
  |                                                                          |_ non mutated bam files
  |_negative control cohort (the background error rate will be called on these samples)
```
### 4.) tune the prior
The prior is tuned on the training dataset. A simple (1D) grid search approach (over binned values of the prior) is used.

### 5.) validation
The tuned prior is validated on the validation dataset.

## References
[1] Gerstung, M., Papaemmanuil, E., & Campbell, P. J. (2014). Subclonal variant calling with multiple samples and prior knowledge. Bioinformatics (Oxford, England), 30(9), 1198–1204. https://doi.org/10.1093/bioinformatics/btt750 \
[2] Ewing, A., Houlahan, K., Hu, Y. et al. Combining tumor genome simulation with crowdsourcing to benchmark somatic single-nucleotide-variant detection. Nat Methods 12, 623–630 (2015). https://doi.org/10.1038/nmeth.3407
