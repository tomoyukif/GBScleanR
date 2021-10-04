## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.pos = 'H', fig.align = "center", warning = FALSE, message = FALSE, out.width = "70%")

## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("GBScleanR")

## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("devtools", quietly = TRUE))
#      install.packages("devtools")
#  devtools::install_github("")

## ----message=FALSE, warning=FALSE---------------------------------------------
library("GBScleanR")

## -----------------------------------------------------------------------------
vcf_fn <- system.file("extdata", "simpop.vcf", package = "GBScleanR")
gds_fn <- system.file("extdata", "simpop.gds", package = "GBScleanR")

## ----eval=FALSE---------------------------------------------------------------
#  gbsrVCF2GDS(vcf_fn = vcf_fn, # Path to the input VCF file.
#              out_fn = gds_fn) # Path to the output GDS file.

## -----------------------------------------------------------------------------
gdata <- loadGDS(gds_fn)

## -----------------------------------------------------------------------------
nsnp(gdata)
nscan(gdata)

## -----------------------------------------------------------------------------
p1 <- grep("Founder1", getScanID(gdata), value = TRUE)
p2 <- grep("Founder2", getScanID(gdata), value = TRUE)
gdata <- setParents(gdata, parents = c(p1, p2))
nsnp(gdata)

## -----------------------------------------------------------------------------
gdata <- countGenotype(gdata)

## -----------------------------------------------------------------------------
histGBSR(gdata, stats = "missing")

## -----------------------------------------------------------------------------
histGBSR(gdata, stats = "het")

## -----------------------------------------------------------------------------
histGBSR(gdata, stats = "raf")

## ----eval=FALSE---------------------------------------------------------------
#  # filter out markers with reference allele frequency
#  # less than 5% or more than 95%.
#  gdata <- setSnpFilter(gdata, maf = 0.05)

## ----eval=FALSE---------------------------------------------------------------
#  # Filter out samples with more than 90% missing genotype calls,
#  # less than 5% heterozygosity, and less than 5% minor allele frequency.
#  gdata <- setScanFilter(gdata, missing = 0.9, het = 0.05, maf = 0.05)

## -----------------------------------------------------------------------------
# Filter out genotype calls supported by reads less than 2 reads.
gdata <- setCallFilter(gdata, dp_count = c(2, Inf))

## -----------------------------------------------------------------------------
gdata <- countGenotype(gdata)

## -----------------------------------------------------------------------------
histGBSR(gdata, stats = "missing")

## -----------------------------------------------------------------------------
# Remove markers having more than 75% of missing genotype calls
gdata <- setSnpFilter(gdata, missing = 0.75) 
nsnp(gdata)
gdata <- countGenotype(gdata)

## -----------------------------------------------------------------------------
histGBSR(gdata, stats = "missing")

## -----------------------------------------------------------------------------
histGBSR(gdata, stats = "het")

## -----------------------------------------------------------------------------
histGBSR(gdata, stats = "raf")

## -----------------------------------------------------------------------------
plotGBSR(gdata, stats = "raf", coord = c(6, 2))

## -----------------------------------------------------------------------------
gdata <- setSnpFilter(gdata, maf = 0.1)
nsnp(gdata)

## -----------------------------------------------------------------------------
gdata <- countGenotype(gdata)
histGBSR(gdata, stats = "missing")

## -----------------------------------------------------------------------------
histGBSR(gdata, stats = "het")

## -----------------------------------------------------------------------------
histGBSR(gdata, stats = "raf")

## -----------------------------------------------------------------------------
# Marker density
plotGBSR(gdata, stats = "marker", coord = c(6, 2))

## -----------------------------------------------------------------------------
plotGBSR(gdata, stats = "geno", coord = c(6, 2))

## -----------------------------------------------------------------------------
subset_gdata <- subsetGDS(gdata,
                          out_fn = "sim_pop_subset.gds")

closeGDS(gdata)

## ----eval = FALSE-------------------------------------------------------------
#  gdata <- loadGDS("sim_pop_subset.gds")

## -----------------------------------------------------------------------------
closeGDS(subset_gdata)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(GBScleanR)
gdata <- loadGDS(gds_fn)

## -----------------------------------------------------------------------------
p1 <- grep("Founder1", getScanID(gdata), value = TRUE)
p2 <- grep("Founder2", getScanID(gdata), value = TRUE)
gdata <- setParents(gdata, parents = c(p1, p2))
nsnp(gdata)

## ----eval=FALSE---------------------------------------------------------------
#  gds <- initScheme(gds, crosstype = "pairing", mating = matrix(1:2, 2))
#  gds <- addScheme(gds, crosstype = "selfing")

## ----eval=FALSE---------------------------------------------------------------
#  gdata <- estGeno(gdata, iter = 4)

## ----eval=FALSE---------------------------------------------------------------
#  gdata <- estGeno(gdata, het_parent = TRUE, iter = 4)

## ----eval=FALSE---------------------------------------------------------------
#  # Following codes do the same.
#  gdata <- estGeno(gdata, iter = 1)
#  gdata <- estGeno(gdata, optim = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  est_geno <- getGenotype(gdata, node = "cor")

## ----eval=FALSE---------------------------------------------------------------
#  founder_geno <- getGenotype(gdata, node = "parents")

## ----eval=FALSE---------------------------------------------------------------
#  est_hap <- getHaplotype(gdata)

## ----eval=FALSE---------------------------------------------------------------
#  gbsrGDS2VCF(gdata, "simpop_est.vcf.gz")

## -----------------------------------------------------------------------------
closeGDS(gdata)

## -----------------------------------------------------------------------------
sessionInfo()

