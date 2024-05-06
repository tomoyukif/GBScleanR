library(SeqArray)

vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, verbose = FALSE)
on.exit({unlink(gds_fn)})

test_that("Number of markers",
          {
              gds <- loadGDS(gds_fn)
              expect_equal(nmar(gds),
                           length(seqGetData(gds, "chromosome")))
              expect_equal(nmar(gds, chr = 1),
                           length(seqGetData(gds, "chromosome") == 1))
              expect_equal(nmar(gds, chr = 100), 0)
              expect_equal(nmar(gds, valid=FALSE),
                           length(seqGetData(gds, "chromosome")))
              expect_equal(nmar(gds, valid=TRUE),
                           length(seqGetData(gds, "chromosome")))
              new_validity <- sample(c(TRUE, FALSE), nmar(gds, FALSE), TRUE)
              validMar(gds) <- new_validity
              expect_equal(nmar(gds, valid=FALSE),
                           length(seqGetData(gds, "chromosome")))
              expect_equal(nmar(gds, valid=TRUE), sum(new_validity))
              closeGDS(gds)
          })

test_that("Number of samples",
          {
              gds <- loadGDS(gds_fn)
              expect_equal(nsam(gds),
                           length(seqGetData(gds, "sample.id")))
              expect_equal(nsam(gds, valid=FALSE),
                           length(seqGetData(gds, "sample.id")))
              expect_equal(nsam(gds, valid=TRUE),
                           length(seqGetData(gds, "sample.id")))
              new_validity <- sample(c(TRUE, FALSE), nsam(gds, FALSE), TRUE)
              validSam(gds) <- new_validity
              expect_equal(nsam(gds, valid=FALSE),
                           length(seqGetData(gds, "sample.id")))
              expect_equal(nsam(gds, valid=TRUE), sum(new_validity))
              closeGDS(gds)
          })

test_that("Chromosome names",
          {
              gds <- loadGDS(gds_fn)
              expect_equal(getChromosome(gds),
                           seqGetData(gds, "chromosome"))
              new_validity <- sample(c(TRUE, FALSE), nmar(gds, FALSE), TRUE)
              validMar(gds) <- new_validity
              expect_equal(getChromosome(gds, valid = FALSE),
                           seqGetData(gds, "chromosome"))
              expect_equal(getChromosome(gds, valid = TRUE),
                           seqGetData(gds, "chromosome")[new_validity])
              closeGDS(gds)
          })


test_that("Marker positions",
          {
              gds <- loadGDS(gds_fn)
              expect_equal(getPosition(gds),
                           seqGetData(gds, "position"))
              expect_equal(getPosition(gds, valid=FALSE),
                           seqGetData(gds, "position"))
              expect_equal(getPosition(gds, valid=TRUE),
                           seqGetData(gds, "position"))
              new_validity <- sample(c(TRUE, FALSE), nmar(gds, FALSE), TRUE)
              validMar(gds) <- new_validity
              expect_equal(getPosition(gds, valid=FALSE),
                           seqGetData(gds, "position"))
              expect_equal(getPosition(gds, valid=TRUE),
                           seqGetData(gds, "position")[new_validity])
              chr_i <- seqGetData(gds, "chromosome") == 1
              expect_equal(getPosition(gds, valid=TRUE, chr = 1),
                           seqGetData(gds, "position")[new_validity & chr_i])
              chr_i <- seqGetData(gds, "chromosome") == 100
              expect_error(getPosition(gds, valid=TRUE, chr = 100),
                           "Nothing to return.")
              closeGDS(gds)
          })

test_that("Reference allele",
          {
              gds <- loadGDS(gds_fn)
              expect_equal(getAllele(gds),
                           seqGetData(gds, "allele"))
              expect_equal(getAllele(gds, valid=FALSE),
                           seqGetData(gds, "allele"))
              expect_equal(getAllele(gds, valid=TRUE),
                           seqGetData(gds, "allele"))
              new_validity <- sample(c(TRUE, FALSE), nmar(gds, FALSE), TRUE)
              validMar(gds) <- new_validity
              expect_equal(getAllele(gds, valid=FALSE),
                           seqGetData(gds, "allele"))
              expect_equal(getAllele(gds, valid=TRUE),
                           seqGetData(gds, "allele")[new_validity])
              chr_i <- seqGetData(gds, "chromosome") == 1
              expect_equal(getAllele(gds, valid=TRUE, chr = 1),
                           seqGetData(gds, "allele")[new_validity & chr_i])
              chr_i <- seqGetData(gds, "chromosome") == 100
              expect_error(getAllele(gds, valid=TRUE, chr = 100),
                           "Nothing to return.")
              closeGDS(gds)
          })

test_that("Marker ID",
          {
              gds <- loadGDS(gds_fn)
              expect_equal(getMarID(gds),
                           seqGetData(gds, "variant.id"))
              expect_equal(getMarID(gds, valid=FALSE),
                           seqGetData(gds, "variant.id"))
              expect_equal(getMarID(gds, valid=TRUE),
                           seqGetData(gds, "variant.id"))
              new_validity <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
              validMar(gds) <- new_validity
              expect_equal(getMarID(gds, valid=FALSE),
                           seqGetData(gds, "variant.id"))
              expect_equal(getMarID(gds, valid=TRUE),
                           seqGetData(gds, "variant.id")[new_validity])
              chr_i <- seqGetData(gds, "chromosome") == 1
              expect_equal(getMarID(gds, valid=TRUE, chr = 1),
                           seqGetData(gds, "variant.id")[new_validity & chr_i])
              chr_i <- seqGetData(gds, "chromosome") == 100
              expect_error(getMarID(gds, valid=TRUE, chr = 100),
                           "Nothing to return.")
              closeGDS(gds)
          })

test_that("Sample ID",
          {
              gds <- loadGDS(gds_fn)
              expect_equal(getSamID(gds),
                           seqGetData(gds, "sample.id"))
              expect_equal(getSamID(gds, valid=FALSE),
                           seqGetData(gds, "sample.id"))
              expect_equal(getSamID(gds, valid=TRUE),
                           seqGetData(gds, "sample.id"))
              new_validity <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
              validSam(gds) <- new_validity
              expect_equal(getSamID(gds, valid=FALSE),
                           seqGetData(gds, "sample.id"))
              expect_equal(getSamID(gds, valid=TRUE),
                           seqGetData(gds, "sample.id")[new_validity])
              closeGDS(gds)
          })

library(gdsfmt)
test_that("Read data",
          {
              gds <- loadGDS(gds_fn)
              expect_type(getRead(gds), "list")
              expect_equal(length(getRead(gds)), 2)
              expect_equal(names(getRead(gds)), c("ref", "alt"))
              read <- read.gdsn(index.gdsn(gds, "annotation/format/AD/data"))
              ref <- read[, c(TRUE, FALSE)]
              alt <- read[, c(FALSE, TRUE)]
              get_read <- getRead(gds, valid = FALSE)
              expect_equal(get_read$ref, ref, ignore_attr = TRUE)
              expect_equal(get_read$alt, alt, ignore_attr = TRUE)

              get_read <- getRead(gds, valid = TRUE)
              expect_equal(get_read$ref, ref, ignore_attr = TRUE)
              expect_equal(get_read$alt, alt, ignore_attr = TRUE)

              valid_snp <- sample(c(TRUE, FALSE), nmar(gds, FALSE), TRUE)
              validMar(gds) <- valid_snp
              valid_scan <- sample(c(TRUE, FALSE), nsam(gds, FALSE), TRUE)
              validSam(gds) <- valid_scan

              get_read <- getRead(gds, valid = FALSE)
              expect_equal(get_read$ref, ref, ignore_attr = TRUE)
              expect_equal(get_read$alt, alt, ignore_attr = TRUE)

              get_read <- getRead(gds, valid = TRUE)
              expect_equal(get_read$ref, ref[valid_scan, valid_snp],
                           ignore_attr = TRUE)
              expect_equal(get_read$alt, alt[valid_scan, valid_snp],
                           ignore_attr = TRUE)

              get_read <- getRead(gds, valid = TRUE, chr = 1)
              expect_equal(get_read$ref, ref[valid_scan, validMar(gds, chr = 1)],
                           ignore_attr = TRUE)
              expect_equal(get_read$alt, alt[valid_scan, validMar(gds, 1)],
                           ignore_attr = TRUE)
              closeGDS(gds)
          })


test_that("Genotype data",
          {
              gds <- loadGDS(gds_fn)
              geno <- read.gdsn(index.gdsn(gds, "genotype/data"))
              geno[geno == 3] <- NA
              geno <- apply(geno, 3, colSums)
              expect_equal(getGenotype(gds, valid = FALSE), geno, ignore_attr = TRUE)
              expect_equal(getGenotype(gds, valid = TRUE), geno, ignore_attr = TRUE)

              valid_snp <- sample(c(TRUE, FALSE), nmar(gds, FALSE), TRUE)
              validMar(gds) <- valid_snp
              valid_scan <- sample(c(TRUE, FALSE), nsam(gds, FALSE), TRUE)
              validSam(gds) <- valid_scan

              expect_equal(getGenotype(gds, valid = FALSE), geno, ignore_attr = TRUE)
              expect_equal(getGenotype(gds, valid = TRUE),
                           geno[validSam(gds), validMar(gds)],
                           ignore_attr = TRUE)

              expect_equal(getGenotype(gds, valid = TRUE, chr = 1),
                           geno[validSam(gds), validMar(gds, 1)],
                           ignore_attr = TRUE)
              expect_error(getGenotype(gds, valid = TRUE, chr = 100),
                           "Nothing to return.")
              closeGDS(gds)
              })

test_that("Hapotype data",
          {
              gds <- loadGDS(gds_fn)
              expect_error(getHaplotype(gds), "Nothing to return.")
              closeGDS(gds)
          })

