library(GBScleanR)

vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, TRUE, FALSE)
on.exit({unlink(gds_fn)})
gds <- loadGDS(gds_fn)

test_that("Number of markers",
          {
              expect_equal(nsnp(gds), 
                           length(getVariable(gds, "snp.chromosome")))
              expect_equal(nsnp(gds, chr=1), 
                           length(getVariable(gds, "snp.chromosome") == 1))
              expect_error(nsnp(gds, chr = 100),
                           "No markers in chromosome 100")
              expect_equal(nsnp(gds, valid=FALSE), 
                           length(getVariable(gds, "snp.chromosome")))
              expect_equal(nsnp(gds, valid=TRUE), 
                           length(getVariable(gds, "snp.chromosome")))
              new_validity <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
              gds <- setValidSnp(gds, new_validity)
              expect_equal(nsnp(gds, valid=FALSE), 
                           length(getVariable(gds, "snp.chromosome")))
              expect_equal(nsnp(gds, valid=TRUE), sum(new_validity))
          })

test_that("Number of samples",
          {
              expect_equal(nscan(gds), 
                           length(getVariable(gds, "sample.id")))
              expect_equal(nscan(gds, valid=FALSE), 
                           length(getVariable(gds, "sample.id")))
              expect_equal(nscan(gds, valid=TRUE), 
                           length(getVariable(gds, "sample.id")))
              new_validity <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
              gds <- setValidScan(gds, new_validity)
              expect_equal(nscan(gds, valid=FALSE), 
                           length(getVariable(gds, "sample.id")))
              expect_equal(nscan(gds, valid=TRUE), sum(new_validity))
          })

test_that("Flipped markers", 
          {
              expect_false(hasFlipped(gds))
              expect_null(getFlipped(gds))
          })

test_that("Chromosome names", 
          {
              expect_equal(getChromosome(gds), 
                           getVariable(gds, "snp.chromosome"))
              expect_equal(getChromosome(gds, name = TRUE), 
                           getVariable(gds, "snp.chromosome.name"))
              expect_equal(getChromosome(gds, levels = TRUE), 
                           unique(getVariable(gds, "snp.chromosome")))
              
              new_validity <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
              gds <- setValidSnp(gds, new_validity)
              expect_equal(getChromosome(gds, valid = FALSE), 
                           getVariable(gds, "snp.chromosome"))
              expect_equal(getChromosome(gds, valid = TRUE), 
                           getVariable(gds, "snp.chromosome")[new_validity])
              expect_equal(getChromosome(gds, valid = FALSE, name = TRUE), 
                           getVariable(gds, "snp.chromosome.name"))
              expect_equal(getChromosome(gds, valid = TRUE, name = TRUE), 
                           getVariable(gds, "snp.chromosome.name")[new_validity])
          })


test_that("Marker positions",
          {
              expect_equal(getPosition(gds), 
                           getVariable(gds, "snp.position"))
              expect_equal(getPosition(gds, valid=FALSE), 
                           getVariable(gds, "snp.position"))
              expect_equal(getPosition(gds, valid=TRUE), 
                           getVariable(gds, "snp.position"))
              new_validity <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
              gds <- setValidSnp(gds, new_validity)
              expect_equal(getPosition(gds, valid=FALSE), 
                           getVariable(gds, "snp.position"))
              expect_equal(getPosition(gds, valid=TRUE),
                           getVariable(gds, "snp.position")[new_validity])
              chr_i <- getVariable(gds, "snp.chromosome") == 1
              expect_equal(getPosition(gds, valid=TRUE, chr = 1),
                           getVariable(gds, "snp.position")[new_validity & chr_i])
              chr_i <- getVariable(gds, "snp.chromosome") == 100
              expect_error(getPosition(gds, valid=TRUE, chr = 100),
                             "No markers in chromosome 100")
          })

test_that("Reference allele",
          {
              expect_equal(getAlleleA(gds), 
                           sub("/.*", "", getVariable(gds, "snp.allele")))
              expect_equal(getAlleleA(gds, valid=FALSE), 
                           sub("/.*", "", getVariable(gds, "snp.allele")))
              expect_equal(getAlleleA(gds, valid=TRUE), 
                           sub("/.*", "", getVariable(gds, "snp.allele")))
              new_validity <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
              gds <- setValidSnp(gds, new_validity)
              expect_equal(getAlleleA(gds, valid=FALSE), 
                           sub("/.*", "", getVariable(gds, "snp.allele")))
              expect_equal(getAlleleA(gds, valid=TRUE),
                           sub("/.*", "",
                               getVariable(gds, "snp.allele"))[new_validity])
              chr_i <- getVariable(gds, "snp.chromosome") == 1
              expect_equal(getAlleleA(gds, valid=TRUE, chr = 1),
                           sub("/.*", "", 
                               getVariable(gds, "snp.allele"))[new_validity & chr_i])
              chr_i <- getVariable(gds, "snp.chromosome") == 100
              expect_error(getAlleleA(gds, valid=TRUE, chr = 100),
                           "No markers in chromosome 100")
          })

test_that("Alternative allele",
          {
              expect_equal(getAlleleB(gds), 
                           sub(".*/", "", getVariable(gds, "snp.allele")))
              expect_equal(getAlleleB(gds, valid=FALSE), 
                           sub(".*/", "", getVariable(gds, "snp.allele")))
              expect_equal(getAlleleB(gds, valid=TRUE), 
                           sub(".*/", "", getVariable(gds, "snp.allele")))
              new_validity <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
              gds <- setValidSnp(gds, new_validity)
              expect_equal(getAlleleB(gds, valid=FALSE), 
                           sub(".*/", "", getVariable(gds, "snp.allele")))
              expect_equal(getAlleleB(gds, valid=TRUE),
                           sub(".*/", "",
                               getVariable(gds, "snp.allele"))[new_validity])
              chr_i <- getVariable(gds, "snp.chromosome") == 1
              expect_equal(getAlleleB(gds, valid=TRUE, chr = 1),
                           sub(".*/", "", 
                               getVariable(gds, "snp.allele"))[new_validity & chr_i])
              chr_i <- getVariable(gds, "snp.chromosome") == 100
              expect_error(getAlleleB(gds, valid=TRUE, chr = 100),
                           "No markers in chromosome 100")
          })

test_that("SNP marker ID",
          {
              expect_equal(getSnpID(gds), 
                           getVariable(gds, "snp.id"))
              expect_equal(getSnpID(gds, valid=FALSE), 
                           getVariable(gds, "snp.id"))
              expect_equal(getSnpID(gds, valid=TRUE), 
                           getVariable(gds, "snp.id"))
              new_validity <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
              gds <- setValidSnp(gds, new_validity)
              expect_equal(getSnpID(gds, valid=FALSE), 
                           getVariable(gds, "snp.id"))
              expect_equal(getSnpID(gds, valid=TRUE),
                           getVariable(gds, "snp.id")[new_validity])
              chr_i <- getVariable(gds, "snp.chromosome") == 1
              expect_equal(getSnpID(gds, valid=TRUE, chr = 1),
                           getVariable(gds, "snp.id")[new_validity & chr_i])
              chr_i <- getVariable(gds, "snp.chromosome") == 100
              expect_error(getSnpID(gds, valid=TRUE, chr = 100),
                           "No markers in chromosome 100")
          })

test_that("Sample ID",
          {
              expect_equal(getScanID(gds), 
                           getVariable(gds, "sample.id"))
              expect_equal(getScanID(gds, valid=FALSE), 
                           getVariable(gds, "sample.id"))
              expect_equal(getScanID(gds, valid=TRUE), 
                           getVariable(gds, "sample.id"))
              new_validity <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
              gds <- setValidScan(gds, new_validity)
              expect_equal(getScanID(gds, valid=FALSE), 
                           getVariable(gds, "sample.id"))
              expect_equal(getScanID(gds, valid=TRUE),
                           getVariable(gds, "sample.id")[new_validity])
          })

test_that("Ploidy",
          {
              expect_equal(getPloidy(gds), 
                           head(getSnpVariable(gds, "ploidy"), 1), 
                           ignore_attr = TRUE)
              expect_equal(getPloidy(gds, valid=FALSE), 
                           head(getSnpVariable(gds, "ploidy"), 1), 
                           ignore_attr = TRUE)
              expect_equal(getPloidy(gds, valid=TRUE), 
                           head(getSnpVariable(gds, "ploidy"), 1), 
                           ignore_attr = TRUE)
              new_validity <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
              gds <- setValidSnp(gds, new_validity)
              expect_equal(getPloidy(gds, valid=FALSE), 
                           head(getSnpVariable(gds, "ploidy"), 1), 
                           ignore_attr = TRUE)
              expect_equal(getPloidy(gds, valid=TRUE), 
                           head(getSnpVariable(gds, "ploidy")[new_validity], 1), 
                           ignore_attr = TRUE)
              chr_i <- getVariable(gds, "snp.chromosome") == 1
              expect_equal(getPloidy(gds, valid=TRUE), 
                           head(getSnpVariable(gds, 
                                               "ploidy")[new_validity & chr_i], 1), 
                           ignore_attr = TRUE)
              chr_i <- getVariable(gds, "snp.chromosome") == 100
              expect_error(getPloidy(gds, valid=TRUE, chr = 100),
                           "No markers in chromosome 100")
          })

library(gdsfmt)
test_that("Read data",
          {
              expect_type(getRead(gds), "list")
              expect_equal(length(getRead(gds)), 2)
              expect_equal(names(getRead(gds)), c("ref", "alt"))
              read <- read.gdsn(index.gdsn(slot(slot(gds, "data"), "handler"), 
                                                "annotation/format/AD/data"))
              ref <- read[, c(TRUE, FALSE)]
              alt <- read[, c(FALSE, TRUE)]
              get_read <- getRead(gds, valid = FALSE)
              expect_equal(get_read$ref, ref, ignore_attr = TRUE)
              expect_equal(get_read$alt, alt, ignore_attr = TRUE)
              
              get_read <- getRead(gds, valid = TRUE)
              expect_equal(get_read$ref, ref, ignore_attr = TRUE)
              expect_equal(get_read$alt, alt, ignore_attr = TRUE)
              
              valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
              gds <- setValidSnp(gds, valid_snp)
              valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
              gds <- setValidScan(gds, valid_scan)
              
              get_read <- getRead(gds, valid = FALSE)
              expect_equal(get_read$ref, ref, ignore_attr = TRUE)
              expect_equal(get_read$alt, alt, ignore_attr = TRUE)
              
              get_read <- getRead(gds, valid = TRUE)
              expect_equal(get_read$ref, ref[valid_scan, valid_snp],
                           ignore_attr = TRUE)
              expect_equal(get_read$alt, alt[valid_scan, valid_snp], 
                           ignore_attr = TRUE)
              
              get_read <- getRead(gds, valid = TRUE, chr = 1)
              expect_equal(get_read$ref, ref[valid_scan, getValidSnp(gds, 1)],
                           ignore_attr = TRUE)
              expect_equal(get_read$alt, alt[valid_scan, getValidSnp(gds, 1)], 
                           ignore_attr = TRUE)
          })


test_that("Genotype data",
          {
              expect_type(getGenotype(gds), "integer")
              geno <- read.gdsn(index.gdsn(slot(slot(gds, "data"), "handler"), 
                                           "genotype"))
              geno[geno == 3] <- NA
              expect_equal(getGenotype(gds, valid = FALSE), geno, ignore_attr = TRUE)
              expect_equal(getGenotype(gds, valid = TRUE), geno, ignore_attr = TRUE)
              
              valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
              gds <- setValidSnp(gds, valid_snp)
              valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
              gds <- setValidScan(gds, valid_scan)
              
              expect_equal(getGenotype(gds, valid = FALSE), geno, ignore_attr = TRUE)
              expect_equal(getGenotype(gds, valid = TRUE),
                           geno[getValidScan(gds), getValidSnp(gds)], 
                           ignore_attr = TRUE)
              
              expect_equal(getGenotype(gds, valid = TRUE, chr = 1),
                           geno[getValidScan(gds), getValidSnp(gds, 1)],
                           ignore_attr = TRUE)
              expect_error(getGenotype(gds, valid = TRUE, chr = 100),
                           "No markers in chromosome 100")
              })

test_that("Hapotype data",
          {
              expect_error(getHaplotype(gds), "No haplotype data")
          })
