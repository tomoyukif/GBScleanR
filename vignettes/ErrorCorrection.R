## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.pos = 'H', fig.align = "center", warning = FALSE, message = FALSE, out.width = "70%")

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
#  install.packages("cowplot")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("path/to/source/GBScleanR.tar.gz", repos = NULL, type = "source")

## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("devtools", quietly = TRUE))
#      install.packages("devtools")
#  devtools::install_github("")

## ----message=FALSE, warning=FALSE---------------------------------------------
library("GBScleanR")

## ----eval=FALSE---------------------------------------------------------------
#  gbsrVCF2GDS(vcf_fn = "./data/gbs_nbolf2.vcf.gz", # Path to the input VCF file.
#              out_fn = "./data/gbs_nbolf2.gds") # Path to the output GDS file.

## -----------------------------------------------------------------------------
gdata <- loadGDS("../inst/extdata/sim_pop.gds")

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
gdata <- calcReadStats(gdata)

## -----------------------------------------------------------------------------
pairsGBSR(gdata, stats1 = "mean_ref", stats2 = "sd_ref", target = "snp", smooth = TRUE)

## -----------------------------------------------------------------------------
pairsGBSR(gdata, stats1 = "mean_ref", stats2 = "sd_ref", target = "snp", smooth = TRUE,
      ggargs = "xlim(c(0, 400)) + ylim(c(0, 400))")

## -----------------------------------------------------------------------------
pairsGBSR(gdata, stats1 = "mean_alt", stats2 = "sd_alt", target = "snp", smooth = TRUE)

## -----------------------------------------------------------------------------
library(mgcv)

## -----------------------------------------------------------------------------
x <- getMeanReadRef(gdata, target = "snp")
y <- getSDReadRef(gdata, target = "snp")
df <- data.frame(x, y)
df <- subset(df, subset = !is.na(x) & !is.na(y))
gam_fit <- gam(formula = y ~ s(x, bs = "cs"), data = df)

## -----------------------------------------------------------------------------
library(ggplot2)
library(cowplot)

## -----------------------------------------------------------------------------
ggplot(df, aes(x = x, y = y)) + geom_point(size = 0.5, color = "darkblue") +
  geom_line(data = data.frame(x = gam_fit$model$x, y = gam_fit$fitted.values),
            mapping = aes(x = x, y = y, group = 1), color = "magenta")

## -----------------------------------------------------------------------------
p1 <- ggplot(data.frame(x = gam_fit$residuals), aes(x = x)) + geom_histogram()
p2 <- ggplot(data.frame(x = gam_fit$residuals), aes(x = x)) + geom_boxplot()
plot_grid(p1, p2, ncol = 1, rel_heights = c(3, 1), align = "v", axis = "lr")

## -----------------------------------------------------------------------------
retain_ref <- rep(FALSE, length(x))
b <- boxplot(gam_fit$residuals, plot = FALSE)
retain_ref[!is.na(x) & !is.na(y)]  <- gam_fit$residuals >= b$stats[1, 1] & gam_fit$residuals <= b$stats[5, 1]

## -----------------------------------------------------------------------------
x <- getMeanReadAlt(gdata, target = "snp")
y <- getSDReadAlt(gdata, target = "snp")
df <- data.frame(x, y)
df <- subset(df, subset = !is.na(x) & !is.na(y))
gam_fit <- gam(formula = y ~ s(x, bs = "cs"), data = df)
retain_alt <- rep(FALSE, length(x))
b <- boxplot(gam_fit$residuals, plot = FALSE)
retain_alt[!is.na(x) & !is.na(y)]  <- gam_fit$residuals >= b$stats[1, 1] & gam_fit$residuals <= b$stats[5, 1]
gdata <- setValidSnp(gdata, update = retain_ref & retain_alt)
nsnp(gdata)

## -----------------------------------------------------------------------------
pairsGBSR(gdata, stats1 = "mean_ref", stats2 = "sd_ref", target = "snp", smooth = TRUE)

## ----warning=FALSE, message=TRUE----------------------------------------------

pairsGBSR(gdata, stats1 = "mean_ref", stats2 = "sd_ref", target = "snp", smooth = TRUE,
      ggargs = "xlim(c(0, 400)) + ylim(c(0, 400))")

## -----------------------------------------------------------------------------
pairsGBSR(gdata, stats1 = "mean_alt", stats2 = "sd_alt", target = "snp", smooth = TRUE)

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
                          out_fn = "../inst/extdata/sim_pop_subset.gds")

closeGDS(gdata)

## ----eval = FALSE-------------------------------------------------------------
#  gdata <- loadGDS("../inst/extdata/sim_pop_subset.gds")

## -----------------------------------------------------------------------------
closeGDS(subset_gdata)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(GBScleanR)
gdata <- loadGDS("../inst/extdata/sim_pop.gds")
gdata

## -----------------------------------------------------------------------------
p1 <- grep("Founder1", getScanID(gdata), value = TRUE)
p2 <- grep("Founder2", getScanID(gdata), value = TRUE)
gdata <- setParents(gdata, parents = c(p1, p2))
nsnp(gdata)

## ----eval=FALSE---------------------------------------------------------------
#  gdata <- estGeno(gdata)

## ----eval=FALSE---------------------------------------------------------------
#  gdata <- countGenotype(gdata, correct = TRUE)
#  plotGBSR(gdata, stats = "geno")

## ----eval=FALSE---------------------------------------------------------------
#  gbsrGDS2VCF(gdata, "./data/gbs_nbolf2_subset_corrected.vcf.gz")

