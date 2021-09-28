# GBScelaneR
Error correction tool for genotype data derived from reduced representation sequencing(RRS).

## Introduction

GBScleanR is a package for quality check, filtering, and error correction
of genotype data derived from next generation sequcener (NGS) based genotyping platforms.
GBScleanR takes Variant Call Format (VCF) file as input. The main function of this package 
is "clean.geno" which estimates the true genotypes of samples from given read counts for
genotype markers using a hidden Markov model with incorporating uneven observation ratio of 
allelic reads. This implementation gives robust genotype estimation even in noisy genotype
data usually observed in Genotyping-By-Sequnencing (GBS) and similar methods, e.g. RADseq.
GBScleanR currenly only supports genotype data of biparental populations.

## Installation
You can install the package from CRAN or GitHub by running the following codes on a R console.
From CRAN:
```
install.package("GBScleanR")
```

From GitHub:
```
install.packages("devtools") # If necessary.
devtools::install_github("tomoyukif/GBScleanR", build_vignettes = TRUE)
```

For more information see vignette or run the following code on a R console.

vignette("GBScleanR-vignette", package = "GBScleanR")
