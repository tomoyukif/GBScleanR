# GBScelaneR
![GBScleanR_icon](https://github.com/tomoyukif/GBScleanR/blob/main/inst/GSBcleanR_Icon.png?raw=true)

Error correction tool for genotype data derived from reduced representation 
sequencing(RRS).

## Introduction

`GBScleanR` is a package for quality check, filtering, and error correction of 
genotype data derived from next generation sequcener (NGS) based genotyping 
platforms. GBScleanR takes Variant Call Format (VCF) file as input. The main
function of this package is `estGeno()` which estimates the true genotypes of 
samples from given read counts for genotype markers using a hidden Markov model 
with incorporating uneven observation ratio of allele reads. This implementation
gives robust genotype estimation even in noisy genotype data usually observed in
Genotyping-By-Sequnencing (GBS) and similar methods, e.g. RADseq. 
The current implementation accepts genotype data of a 
diploid population at any generation of multi-parental cross, e.g. biparental
F2 from inbred parents, biparental F2 from outbred parents, and 8-way 
recombinant inbred lines (8-way RILs) which can be refered to as MAGIC 
population.

## Installation
You can install `GBScleanR` from the Bioconductor repository with running the 
following code on an R console. Currently `GBScleanR` is listed in the 
developmental package repository of Bioconductor(https://doi.org/doi:10.18129/B9.bioc.GBScleanR). 
You need to set to use Bioc 
devel as shown below. On the other hand, `GBScleanR` on the GitHub repository
provide you the latest package with some developmental functions. If you need 
the stable release version of the package, please install the one on 
Bioconductor.

From Bioconductor:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("GBScleanR")
```

From GitHub:
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("tomoyukif/GBScleanR", build_vignettes = TRUE)
```

For more information see vignette or run the following code on a R console.
```
browseVignettes(package = "GBScleanR")
```

## Citations
"Furuta, T., Yamamoto, T., Ashikari, M., GBScleanR: robust genotyping error 
correction using a hidden Markov model with error pattern recognition. 
Genetics 2023; doi: [10.1093/genetics/iyad055](https://doi.org/10.1093/genetics/iyad055)
