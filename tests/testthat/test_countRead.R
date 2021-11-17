library(GBScleanR)

vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, TRUE, FALSE)
on.exit({unlink(gds_fn)})
gds <- loadGDS(gds_fn)
gds <- countRead(gds)

read <- getRead(gds)
ref <- read$ref
alt <- read$alt

test_that("Reference read count",
          {
              expect_equal(getCountReadRef(gds, "snp"), colSums(ref),
                           ignore_attr = TRUE)
              expect_equal(getCountReadRef(gds, "scan"), rowSums(ref),
                           ignore_attr = TRUE)
              expect_equal(getCountReadRef(gds, "snp", prop = TRUE),
                           colSums(ref)/ (colSums(ref) + colSums(alt)),
                           ignore_attr = TRUE)
              expect_equal(getCountReadRef(gds, "scan", prop = TRUE),
                           rowSums(ref) / (rowSums(ref) + rowSums(alt)),
                           ignore_attr = TRUE)
              
              valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
              gds <- setValidSnp(gds, valid_snp)
              valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
              gds <- setValidScan(gds, valid_scan)
              gds <- countRead(gds)
              ref <- read$ref
              alt <- read$alt
              
              expect_equal(getCountReadRef(gds, "snp", valid = TRUE), 
                           colSums(ref[valid_scan, ])[valid_snp],
                           ignore_attr = TRUE)
              expect_equal(getCountReadRef(gds, "scan", valid = TRUE),
                           rowSums(ref[, valid_snp])[valid_scan],
                           ignore_attr = TRUE)
              
              expect_equal(getCountReadRef(gds, "snp", valid = TRUE, prop = TRUE), 
                           (colSums(ref[valid_scan, ])/
                                (colSums(ref[valid_scan, ]) + 
                                     colSums(alt[valid_scan, ])))[valid_snp],
                           ignore_attr = TRUE)
              expect_equal(getCountReadRef(gds, "scan", valid = TRUE, prop = TRUE),
                           (rowSums(ref[, valid_snp])/
                                (rowSums(ref[, valid_snp]) +
                                     rowSums(alt[, valid_snp])))[valid_scan],
                           ignore_attr = TRUE)
          })


test_that("Alternative read count",
          {
              expect_equal(getCountReadAlt(gds, "snp"), colSums(alt),
                           ignore_attr = TRUE)
              expect_equal(getCountReadAlt(gds, "scan"), rowSums(alt),
                           ignore_attr = TRUE)
              expect_equal(getCountReadAlt(gds, "snp", prop = TRUE),
                           colSums(alt)/ (colSums(alt) + colSums(ref)),
                           ignore_attr = TRUE)
              expect_equal(getCountReadAlt(gds, "scan", prop = TRUE),
                           rowSums(alt) / (rowSums(alt) + rowSums(ref)),
                           ignore_attr = TRUE)
              
              valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
              gds <- setValidSnp(gds, valid_snp)
              valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
              gds <- setValidScan(gds, valid_scan)
              gds <- countRead(gds)
              ref <- read$ref
              alt <- read$alt
              
              expect_equal(getCountReadAlt(gds, "snp", valid = TRUE), 
                           colSums(alt[valid_scan, ])[valid_snp],
                           ignore_attr = TRUE)
              expect_equal(getCountReadAlt(gds, "scan", valid = TRUE),
                           rowSums(alt[, valid_snp])[valid_scan],
                           ignore_attr = TRUE)
              
              expect_equal(getCountReadAlt(gds, "snp", valid = TRUE, prop = TRUE), 
                           (colSums(alt[valid_scan, ])/
                                (colSums(ref[valid_scan, ]) + 
                                     colSums(alt[valid_scan, ])))[valid_snp],
                           ignore_attr = TRUE)
              expect_equal(getCountReadAlt(gds, "scan", valid = TRUE, prop = TRUE),
                           (rowSums(alt[, valid_snp])/
                                (rowSums(ref[, valid_snp]) +
                                     rowSums(alt[, valid_snp])))[valid_scan],
                           ignore_attr = TRUE)
          })

dp <- ref + alt
test_that("Total read count", {
    expect_equal(getCountRead(gds, "snp"), colSums(dp),
                 ignore_attr = TRUE)
    expect_equal(getCountRead(gds, "scan"), rowSums(dp),
                 ignore_attr = TRUE)
    
    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- countRead(gds)
    ref <- read$ref
    alt <- read$alt
    
    expect_equal(getCountRead(gds, "snp", valid = TRUE), 
                 colSums(dp[valid_scan, ])[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountRead(gds, "scan", valid = TRUE),
                 rowSums(dp[, valid_snp])[valid_scan],
                 ignore_attr = TRUE)
})
