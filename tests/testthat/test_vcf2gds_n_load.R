library(GBScleanR)
vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
test_that("Loading the GDS file",
          {
              gds_fn <- tempfile("sample", fileext = ".gds")
              gds_fn <- gbsrVCF2GDS(vcf_fn, gds_fn, TRUE)
              gds <- loadGDS(gds_fn)
              stop()
              expect_s4_class(gds, "GbsrGenotypeData")
              expect_s4_class(getSnpAnnotation(gds), "SnpAnnotationDataFrame")
              expect_s4_class(getScanAnnotation(gds), "ScanAnnotationDataFrame")
              unlink(gds_fn)
          })

test_that("Reloading the GDS file",
          {
              gds_fn <- tempfile("sample", fileext = ".gds")
              gds_fn <- gbsrVCF2GDS(vcf_fn, gds_fn, TRUE)
              gds <- loadGDS(gds_fn)
              gds <- loadGDS(gds)
              expect_s4_class(gds, "GbsrGenotypeData")
              unlink(gds_fn)
          })
