## Theoretical background

The goal of this PhD project is to increase the senistivity of variant calling of actionable variants. Therefore loci that are not being called by existing software are of interest.
If a given locus is not being called, it's due to these two main reasons:
  - it is truly biologically absent (no variant) --> true negative (TN)
  - it is a true variant that is missed by the software --> false negative (FN) --> **we want to capture these**

To achieve this a mathematical model is needed (most probably a Deep Learning model). Therefore training data, where the targets will be the TN and FN sites, is needed.
Since there exists no ground truth for this data, we need to guess whether a given locus that is not called is a TN or an FN. This is the task of the following workflow.

We have data from patients with multiple tumor biopsies.

<p align="center">
  <img src="https://github.com/JakubLiu/Supervised_detection_of_clinically_actionable_cancer_variants/blob/main/ACTIONABLE_VARIANTS_INFORMATION_GATHERING/FN_vs_TN_workflow/venn2.png.png" width="300">
</p>

If a variant has been not consistently called across all biopsy samples from the same patient (e.g. called in 2/3 biopsies), then it is of interest to us
and we will reffer to it as a **candidate**.
We would like to know if it has not been called due to a fault of the software or due to being truly absent. One main reason why standard variant callers
miss true variants is a low read coverage depth. Therefore we create a copy of the sequencing data with a smaller read coverage and perform variant calling on this
subsampled dataset. Next, variant calling on the original sequencing data (pre downsampling) is performed. Now, for each candidate variant there are the following options:

| called in subsampled data | called in original data | explanation |
|---------------------------|--------------------------|--------------|
| YES                     | YES                  | likely a true biological variant      |
| NO                     | NO                  | likely truly biologically absent      |
| NO                     | YES                  | likely biologically preset but missed by the variant caller      |
| YES                     | NO                  | something is very wrong      |

Therefore after running the pipeline on a large cohort of patients, the goal is to classify the "likely truly biologically absent" variants as TNs and the 
"likely biologically preset but missed by the variant caller" variants as FNs. The next step would be to train a classifier on such an annotated dataset.

