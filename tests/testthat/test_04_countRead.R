library(GBScleanR)

vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn)
on.exit({unlink(gds_fn)})

test_that("Reference read count",
          {
              gds <- loadGDS(gds_fn)
              gds <- countRead(gds)
              read <- getRead(gds)
              ref <- read$ref
              alt <- read$alt
              expect_equal(getCountReadRef(gds, "marker"), colSums(ref),
                           ignore_attr = TRUE)
              expect_equal(getCountReadRef(gds, "sample"), rowSums(ref),
                           ignore_attr = TRUE)
              expect_equal(getCountReadRef(gds, "marker", prop = TRUE),
                           colSums(ref)/ (colSums(ref) + colSums(alt)),
                           ignore_attr = TRUE)
              expect_equal(getCountReadRef(gds, "sample", prop = TRUE),
                           rowSums(ref) / (rowSums(ref) + rowSums(alt)),
                           ignore_attr = TRUE)

              valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
              validMar(gds) <- valid_snp
              valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
              validSam(gds) <- valid_scan
              gds <- countRead(gds)
              ref <- read$ref
              alt <- read$alt

              expect_equal(getCountReadRef(gds, "marker", valid = TRUE),
                           colSums(ref[valid_scan, ])[valid_snp],
                           ignore_attr = TRUE)
              expect_equal(getCountReadRef(gds, "sample", valid = TRUE),
                           rowSums(ref[, valid_snp])[valid_scan],
                           ignore_attr = TRUE)

              expect_equal(getCountReadRef(gds, "marker", valid = TRUE, prop = TRUE),
                           (colSums(ref[valid_scan, ])/
                                (colSums(ref[valid_scan, ]) +
                                     colSums(alt[valid_scan, ])))[valid_snp],
                           ignore_attr = TRUE)
              expect_equal(getCountReadRef(gds, "sample", valid = TRUE, prop = TRUE),
                           (rowSums(ref[, valid_snp])/
                                (rowSums(ref[, valid_snp]) +
                                     rowSums(alt[, valid_snp])))[valid_scan],
                           ignore_attr = TRUE)
              closeGDS(gds)
          })


test_that("Alternative read count",
          {
              gds <- loadGDS(gds_fn)
              gds <- countRead(gds)
              read <- getRead(gds)
              ref <- read$ref
              alt <- read$alt
              expect_equal(getCountReadAlt(gds, "marker"), colSums(alt),
                           ignore_attr = TRUE)
              expect_equal(getCountReadAlt(gds, "sample"), rowSums(alt),
                           ignore_attr = TRUE)
              expect_equal(getCountReadAlt(gds, "marker", prop = TRUE),
                           colSums(alt)/ (colSums(alt) + colSums(ref)),
                           ignore_attr = TRUE)
              expect_equal(getCountReadAlt(gds, "sample", prop = TRUE),
                           rowSums(alt) / (rowSums(alt) + rowSums(ref)),
                           ignore_attr = TRUE)

              valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
              validMar(gds) <- valid_snp
              valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
              validSam(gds) <- valid_scan
              gds <- countRead(gds)
              ref <- read$ref
              alt <- read$alt

              expect_equal(getCountReadAlt(gds, "marker", valid = TRUE),
                           colSums(alt[valid_scan, ])[valid_snp],
                           ignore_attr = TRUE)
              expect_equal(getCountReadAlt(gds, "sample", valid = TRUE),
                           rowSums(alt[, valid_snp])[valid_scan],
                           ignore_attr = TRUE)

              expect_equal(getCountReadAlt(gds, "marker", valid = TRUE, prop = TRUE),
                           (colSums(alt[valid_scan, ])/
                                (colSums(ref[valid_scan, ]) +
                                     colSums(alt[valid_scan, ])))[valid_snp],
                           ignore_attr = TRUE)
              expect_equal(getCountReadAlt(gds, "sample", valid = TRUE, prop = TRUE),
                           (rowSums(alt[, valid_snp])/
                                (rowSums(ref[, valid_snp]) +
                                     rowSums(alt[, valid_snp])))[valid_scan],
                           ignore_attr = TRUE)
              closeGDS(gds)
          })

test_that("Total read count", {
    gds <- loadGDS(gds_fn)
    gds <- countRead(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    dp <- ref + alt
    expect_equal(getCountRead(gds, "marker"), colSums(dp),
                 ignore_attr = TRUE)
    expect_equal(getCountRead(gds, "sample"), rowSums(dp),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    validMar(gds) <- valid_snp
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    validSam(gds) <- valid_scan
    gds <- countRead(gds)
    ref <- read$ref
    alt <- read$alt

    expect_equal(getCountRead(gds, "marker", valid = TRUE),
                 colSums(dp[valid_scan, ])[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountRead(gds, "sample", valid = TRUE),
                 rowSums(dp[, valid_snp])[valid_scan],
                 ignore_attr = TRUE)
    closeGDS(gds)
})

