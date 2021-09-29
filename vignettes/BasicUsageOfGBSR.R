## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.pos = 'H', fig.align = "center", warning = FALSE, message = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("GWASTools")
#  BiocManager::install("SNPRelate")
#  BiocManager::install("SeqArray")
#  install.packages("ggplot2")
#  install.packages("dplyr")
#  install.packages("tidyr")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("path/to/source/GBScleanR.tar.gz", repos = NULL, type = "source")

## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("devtools", quietly = TRUE))
#      install.packages("devtools")
#  devtools::install_github("")

## ----warning=FALSE, message=FALSE---------------------------------------------
library("GBScleanR")

## ----warning=FALSE, message=FALSE, eval=FALSE---------------------------------
#  gbsrVCF2GDS(vcf_fn = "./data/gbs_nbolf2.vcf.gz", # Path to the input VCF file.
#              out_fn = "./data/gbs_nbolf2.gds") # Path to the output GDS file.

## ----eval=FALSE---------------------------------------------------------------
#  exec_time <- system.time({
#    gbsrVCF2GDS(vcf_fn = "./data/gbs_nbolf2.vcf.gz", # Path to the input VCF file.
#                out_fn = "./data/gbs_nbolf2.gds") # Path to the output GDS file.
#  })
#  exec_time

## -----------------------------------------------------------------------------
gdata <- loadGDS("../inst/extdata/sim_pop.gds")

## ----eval=FALSE---------------------------------------------------------------
#  # Not run.
#  gdata <- loadGDS("./data/gbs_nbolf2.gds",
#                   non_autosomes =  list(X = 13,
#                                         Y = 14,
#                                         M = 15)) # M indicates mitochondrial chromosome.

## -----------------------------------------------------------------------------
nscan(gdata) # Number of samples

## -----------------------------------------------------------------------------
nsnp(gdata) # Number of SNPs

## -----------------------------------------------------------------------------
head(getChromosome(gdata)) # Indices of chromosome ID of all markers

## -----------------------------------------------------------------------------
head(getChromosome(gdata, name = TRUE)) # Chromosome names of all markers

## -----------------------------------------------------------------------------
getChromosome(gdata, levels = TRUE) # Unique set of chromosome names

## -----------------------------------------------------------------------------
head(getPosition(gdata)) # Position (bp) of all markers

## -----------------------------------------------------------------------------
head(getAlleleA(gdata)) # Reference allele of all markers

## -----------------------------------------------------------------------------
head(getAlleleB(gdata)) # Alternative allele of all markers

## -----------------------------------------------------------------------------
head(getSnpID(gdata)) # SNP IDs

## -----------------------------------------------------------------------------
head(getScanID(gdata)) # sample IDs

## -----------------------------------------------------------------------------
g <- getGenotype(gdata) # Genotype calls in which 0, 1, and 2 indicate the number of reference allele.

## -----------------------------------------------------------------------------
gdata <- countGenotype(gdata)
gdata <- countRead(gdata)

## -----------------------------------------------------------------------------
gdata@snpAnnot

## -----------------------------------------------------------------------------
gdata@scanAnnot

## -----------------------------------------------------------------------------
head(pData(gdata@snpAnnot), n = 3)

## -----------------------------------------------------------------------------
head(pData(gdata@scanAnnot), n = 3)

## ----fig.cap="Missing rate per marker and per sample.", out.height="35%"------
histGBSR(gdata, stats = "missing") # Histgrams of missing rate

## ----fig.cap="Heterozygosity per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "het") # Histgrams of heterozygosity

## ----fig.cap="Reference allele frequency per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "raf") # Histgrams of reference allele frequency

## ----fig.cap="Total read depth per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "dp") # Histgrams of total read depth

## ----fig.cap="Reference read depth per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "ad_ref") # Histgrams of allelic read depth

## ----fig.cap="Alternative read depth per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "ad_ref") # Histgrams of allelic read depth

## ----fig.cap="Reference read per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "rrf") # Histgrams of reference allele frequency

## -----------------------------------------------------------------------------
gdata <- calcReadStats(gdata, q = 0.5)

## ----fig.cap="Mean of reference read depth per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "mean_ref") # Histgrams of mean allelic read depth

## ----fig.cap="Mean of alternative read depth per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "mean_ref") # Histgrams of mean allelic read depth

## ----fig.cap="SD of reference read depth per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "sd_ref") # Histgrams of standard deviation of read depth

## ----fig.cap="SD of alternative read depth per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "sd_ref") # Histgrams of standard deviation of read depth

## ----fig.cap="Quantile of reference read depth per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "qtile_ref", q = 0.5) # Histgrams of quantile of read depth

## ----fig.cap="Quantile of alternative read depth per marker and per sample.", out.height="35%"----
histGBSR(gdata, stats = "qtile_ref", q = 0.5) # Histgrams of quantile of read depth

## -----------------------------------------------------------------------------
plotGBSR(gdata, stats = "missing", coord = c(6, 2)) # coord controls the number of rows and columns of facets.

## -----------------------------------------------------------------------------
plotGBSR(gdata, stats = "geno", coord = c(6, 2)) # coord controls the number of rows and columns of facets.

## -----------------------------------------------------------------------------
pairsGBSR(gdata, stats1 = "missing", stats2 = "dp")

## -----------------------------------------------------------------------------
head(getCountGenoRef(gdata, target = "snp")) # Reference genotype count per marker
head(getCountGenoRef(gdata, target = "scan")) # Reference genotype count per sample

## -----------------------------------------------------------------------------
head(getCountGenoHet(gdata, target = "snp")) # Heterozygote count per marker
head(getCountGenoHet(gdata, target = "scan")) # Heterozygote count per sample

## -----------------------------------------------------------------------------
head(getCountGenoAlt(gdata, target = "snp")) # Alternative genotype count per marker
head(getCountGenoAlt(gdata, target = "scan")) # Alternative genotype count per sample

## -----------------------------------------------------------------------------
head(getCountGenoMissing(gdata, target = "snp")) # Missing count per marker
head(getCountGenoMissing(gdata, target = "scan")) # Missing count per sample

## -----------------------------------------------------------------------------
head(getCountAlleleRef(gdata, target = "snp")) # Reference allele count per marker
head(getCountAlleleRef(gdata, target = "scan")) # Reference allele count per sample

## -----------------------------------------------------------------------------
head(getCountAlleleAlt(gdata, target = "snp")) # Alternative allele count per marker
head(getCountAlleleAlt(gdata, target = "scan")) # Alternative allele count per sample

## -----------------------------------------------------------------------------
head(getCountAlleleMissing(gdata, target = "snp")) # Missing allele count per marker
head(getCountAlleleMissing(gdata, target = "scan")) # Missing allele count per sample

## -----------------------------------------------------------------------------
head(getCountReadRef(gdata, target = "snp")) # Reference read count per marker
head(getCountReadRef(gdata, target = "scan")) # Reference read count per sample

## -----------------------------------------------------------------------------
head(getCountReadAlt(gdata, target = "snp")) # Alternative read count per marker
head(getCountReadAlt(gdata, target = "scan")) # Alternative read count per sample

## -----------------------------------------------------------------------------
head(getCountRead(gdata, target = "snp")) # Sum of reference and alternative read counts per marker
head(getCountRead(gdata, target = "scan")) # Sum of reference and alternative read counts per sample

## -----------------------------------------------------------------------------
head(getMeanReadRef(gdata, target = "snp")) # Mean of reference allele read count per marker
head(getMeanReadRef(gdata, target = "scan")) # Mean of reference allele read count per sample

## -----------------------------------------------------------------------------
head(getMeanReadAlt(gdata, target = "snp")) # Mean of Alternative allele read count per marker
head(getMeanReadAlt(gdata, target = "scan")) # Mean of Alternative allele read count per sample

## -----------------------------------------------------------------------------
head(getSDReadRef(gdata, target = "snp")) # SD of reference allele read count per marker
head(getSDReadRef(gdata, target = "scan")) # SD of reference allele read count per sample

## -----------------------------------------------------------------------------
head(getSDReadAlt(gdata, target = "snp")) # SD of Alternative allele read count per marker
head(getSDReadAlt(gdata, target = "scan")) # SD of Alternative allele read count per sample

## -----------------------------------------------------------------------------
head(getQtileReadRef(gdata, target = "snp", q = 0.5)) # Quantile of reference allele read count per marker
head(getQtileReadRef(gdata, target = "scan", q = 0.5)) # Quantile of reference allele read count per sample

## -----------------------------------------------------------------------------
head(getQtileReadAlt(gdata, target = "snp", q = 0.5)) # Quantile of Alternative allele read count per marker
head(getQtileReadAlt(gdata, target = "scan", q = 0.5)) # Quantile of Alternative allele read count per sample

## -----------------------------------------------------------------------------
head(getMAF(gdata, target = "snp")) # Minor allele frequency per marker
head(getMAF(gdata, target = "scan")) # Minor allele frequency per sample

## -----------------------------------------------------------------------------
head(getMAC(gdata, target = "snp")) # Minor allele count per marker
head(getMAC(gdata, target = "scan")) # Minor allele count per sample

## -----------------------------------------------------------------------------
head(getCountGenoRef(gdata, target = "snp", prop = TRUE))
head(getCountGenoHet(gdata, target = "snp", prop = TRUE))
head(getCountGenoAlt(gdata, target = "snp", prop = TRUE))
head(getCountGenoMissing(gdata, target = "snp", prop = TRUE))

## -----------------------------------------------------------------------------
head(getCountAlleleRef(gdata, target = "snp", prop = TRUE))
head(getCountAlleleAlt(gdata, target = "snp", prop = TRUE))
head(getCountAlleleMissing(gdata, target = "snp", prop = TRUE))

## -----------------------------------------------------------------------------
head(getCountReadRef(gdata, target = "snp", prop = TRUE))
head(getCountReadAlt(gdata, target = "snp", prop = TRUE))

## ----eval=FALSE---------------------------------------------------------------
#  # Not run
#  gdata <- setSnpFilter(
#    id,   # Specify a character vector of snpID to be removed.
#    missing = 1,   # Specify an upper limit of missing rate.
#    het = c(0, 1),   # Specify a lower and an upper limit of heterozygosity rate.
#    mac = 0,   # Specify a lower limit of minor allele count.
#    maf = 0.05,   # Specify a lower limit of minor allele frequency.
#    ad_ref = c(0, Inf),   # Specify a lower and an upper limit of reference allele count.
#    ad_alt = c(0, Inf),   # Specify a lower and an upper limit of alternative allele count.
#    dp = c(0, Inf),   # Specify a lower and an upper limit of total read count.
#    mean_ref = c(0, Inf),   # Specify a lower and an upper limit of mean reference allele count.
#    mean_alt = c(0, Inf),   # Specify a lower and an upper limit of mean alternative allele count.
#    sd_ref = Inf,   # Specify a lower and an upper limit of SD of reference allele count.
#    sd_alt = Inf    # Specify a lower and an upper limit of SD of alternative allele count.
#  )
#  
#  gdata <- setScanFilter(
#    id,   # Specify a character vector of snpID to be removed.
#    missing = 1,   # Specify an upper limit of missing rate.
#    het = c(0, 1),   # Specify a lower and an upper limit of heterozygosity rate.
#    mac = 0,   # Specify a lower limit of minor allele count.
#    maf = 0,   # Specify a lower limit of minor allele frequency.
#    ad_ref = c(0, Inf),   # Specify a lower and an upper limit of reference allele count.
#    ad_alt = c(0, Inf),   # Specify a lower and an upper limit of alternative allele count.
#    dp = c(0, Inf),   # Specify a lower and an upper limit of total read count.
#    mean_ref = c(0, Inf),   # Specify a lower and an upper limit of mean reference allele count.
#    mean_alt = c(0, Inf),   # Specify a lower and an upper limit of mean alternative allele count.
#    sd_ref = Inf,   # Specify a lower and an upper limit of SD of reference allele count.
#    sd_alt = Inf   # Specify a lower and an upper limit of SD of alternative allele count.
#  )

## ----eval=FALSE---------------------------------------------------------------
#  gdata <- setCallFilter(gdata, dp_count = c(5, Inf))

## ----eval=FALSE---------------------------------------------------------------
#  gdata <- setCallFilter(gdata, norm_dp_count = c(0, 1000))
#  gdata <- setCallFilter(gdata, norm_ref_count = c(0, 1000),
#                         norm_alt_count = c(0, 800))

## -----------------------------------------------------------------------------
gdata <- setCallFilter(gdata, dp_count = c(5, Inf))
gdata <- setSnpFilter(gdata, missing = 0.1)

## -----------------------------------------------------------------------------
thinMarker(gdata, range = 150) # Here we select only one marker from each 150 bp stretch.

## ----eval= FALSE--------------------------------------------------------------
#  gdata <- countGenotype(gdata)
#  gdata <- countRead(gdata)
#  gdata <- calcReadStats(gdata)

## -----------------------------------------------------------------------------
head(getValidSnp(gdata))
head(getValidScan(gdata))

## -----------------------------------------------------------------------------
nsnp(gdata)

## -----------------------------------------------------------------------------
nsnp(gdata, valid = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  gdata <- resetSnpFilters(gdata) # Reset the filter on markers
#  gdata <- resetScanFilters(gdata) # Reset the filter on samples
#  gdata <- resetCallFilters(gdata) # Reset the filter on calls
#  gdata <- resetFilters(gdata) # Reset all filters

## ----eval=FALSE---------------------------------------------------------------
#  subset_gdata <- subsetGDS(gdata,
#                            out_fn = "./data/gbs_nbolf2_subset.gds",
#                            snp_incl = getValidSnp(gdata),
#                            scan_incl = getValidScan(gdata))

## -----------------------------------------------------------------------------
closeGDS(gdata)

## -----------------------------------------------------------------------------
sessionInfo()

