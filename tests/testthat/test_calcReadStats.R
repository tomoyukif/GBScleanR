library(GBScleanR)

vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, TRUE, FALSE)
on.exit({unlink(gds_fn)})
gds <- loadGDS(gds_fn)
q <- 0.5
gds <- calcReadStats(gds, q = q)
read <- getRead(gds, valid = FALSE)
ref <- read$ref
alt <- read$alt
dp <- ref + alt
tot <- rowSums(dp)
ref <- ref / tot * 10^6
alt <- alt / tot * 10^6
ref[ref == 0] <- NA
alt[alt == 0] <- NA

test_that("Mean reference read", {
    expect_equal(getMeanReadRef(gds, "snp"),
                 colMeans(ref, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getMeanReadRef(gds, "scan"),
                 rowMeans(ref, na.rm = TRUE),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- calcReadStats(gds, q = 0.5)
    read <- getRead(gds, valid = FALSE)
    ref <- read$ref
    alt <- read$alt
    dp <- ref + alt
    tot <- rowSums(dp)
    ref <- ref / tot * 10^6
    alt <- alt / tot * 10^6

    valid_ref <- ref[valid_scan, valid_snp]
    valid_ref[valid_ref == 0] <- NA
    expect_equal(getMeanReadRef(gds, "snp", valid = TRUE),
                 colMeans(valid_ref, na.rm = TRUE),
                 ignore_attr = TRUE)

    expect_equal(getMeanReadRef(gds, "scan", valid = TRUE),
                 rowMeans(valid_ref, na.rm = TRUE),
                 ignore_attr = TRUE)
})

test_that("Mean alternative read", {
    expect_equal(getMeanReadAlt(gds, "snp"),
                 colMeans(alt, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getMeanReadAlt(gds, "scan"),
                 rowMeans(alt, na.rm = TRUE),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- calcReadStats(gds, q = 0.5)
    read <- getRead(gds, valid = FALSE)
    ref <- read$ref
    alt <- read$alt
    dp <- ref + alt
    tot <- rowSums(dp)
    ref <- ref / tot * 10^6
    alt <- alt / tot * 10^6

    valid_alt <- alt[valid_scan, valid_snp]
    valid_alt[valid_alt == 0] <- NA
    expect_equal(getMeanReadAlt(gds, "snp", valid = TRUE),
                 colMeans(valid_alt, na.rm = TRUE),
                 ignore_attr = TRUE)

    expect_equal(getMeanReadAlt(gds, "scan", valid = TRUE),
                 rowMeans(valid_alt, na.rm = TRUE),
                 ignore_attr = TRUE)
})

test_that("SD of reference read", {
    expect_equal(as.integer(getSDReadRef(gds, "snp")),
                 as.integer(apply(ref, 2, sd, na.rm = TRUE)),
                 ignore_attr = TRUE)
    expect_equal(as.integer(getSDReadRef(gds, "scan")),
                 as.integer(apply(ref, 1, sd, na.rm = TRUE)),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- calcReadStats(gds, q = 0.5)
    read <- getRead(gds, valid = FALSE)
    ref <- read$ref
    alt <- read$alt
    dp <- ref + alt
    tot <- rowSums(dp)
    ref <- ref / tot * 10^6
    alt <- alt / tot * 10^6

    valid_ref <- ref[valid_scan, valid_snp]
    valid_ref[valid_ref == 0] <- NA
    expect_equal(as.integer(getSDReadRef(gds, "snp", valid = TRUE)),
                 as.integer(apply(valid_ref, 2, sd, na.rm = TRUE)),
                 ignore_attr = TRUE)

    expect_equal(as.integer(getSDReadRef(gds, "scan", valid = TRUE)),
                 as.integer(apply(valid_ref, 1, sd, na.rm = TRUE)),
                 ignore_attr = TRUE)
})

test_that("SD of alternative read", {
    expect_equal(as.integer(getSDReadAlt(gds, "snp")),
                 as.integer(apply(alt, 2, sd, na.rm = TRUE)),
                 ignore_attr = TRUE)
    expect_equal(as.integer(getSDReadAlt(gds, "scan")),
                 as.integer(apply(alt, 1, sd, na.rm = TRUE)),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- calcReadStats(gds, q = 0.5)
    read <- getRead(gds, valid = FALSE)
    ref <- read$ref
    alt <- read$alt
    dp <- ref + alt
    tot <- rowSums(dp)
    ref <- ref / tot * 10^6
    alt <- alt / tot * 10^6

    valid_alt <- alt[valid_scan, valid_snp]
    valid_alt[valid_alt == 0] <- NA
    expect_equal(as.integer(getSDReadAlt(gds, "snp", valid = TRUE)),
                 as.integer(apply(valid_alt, 2, sd, na.rm = TRUE)),
                 ignore_attr = TRUE)

    expect_equal(as.integer(getSDReadAlt(gds, "scan", valid = TRUE)),
                 as.integer(apply(valid_alt, 1, sd, na.rm = TRUE)),
                 ignore_attr = TRUE)
})

test_that("Quantile of reference read", {
    expect_equal(as.integer(getQtileReadRef(gds, "snp")),
                 as.integer(apply(ref, 2, quantile, probs = q, na.rm = TRUE)),
                 ignore_attr = TRUE)
    expect_equal(as.integer(getQtileReadRef(gds, "scan")),
                 as.integer(apply(ref, 1, quantile, probs = q, na.rm = TRUE)),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- calcReadStats(gds, q = 0.5)
    read <- getRead(gds, valid = FALSE)
    ref <- read$ref
    alt <- read$alt
    dp <- ref + alt
    tot <- rowSums(dp)
    ref <- ref / tot * 10^6
    alt <- alt / tot * 10^6

    valid_ref <- ref[valid_scan, valid_snp]
    valid_ref[valid_ref == 0] <- NA
    expect_equal(as.integer(getQtileReadRef(gds, "snp", valid = TRUE)),
                 as.integer(apply(valid_ref, 2, quantile, probs = q, na.rm = TRUE)),
                 ignore_attr = TRUE)

    expect_equal(as.integer(getQtileReadRef(gds, "scan", valid = TRUE)),
                 as.integer(apply(valid_ref, 1, quantile, probs = q, na.rm = TRUE)),
                 ignore_attr = TRUE)
})


test_that("Quantile of alternative read", {
    expect_equal(as.integer(getQtileReadAlt(gds, "snp")),
                 as.integer(apply(alt, 2, quantile, probs = q, na.rm = TRUE)),
                 ignore_attr = TRUE)
    expect_equal(as.integer(getQtileReadAlt(gds, "scan")),
                 as.integer(apply(alt, 1, quantile, probs = q, na.rm = TRUE)),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- calcReadStats(gds, q = 0.5)
    read <- getRead(gds, valid = FALSE)
    ref <- read$ref
    alt <- read$alt
    dp <- ref + alt
    tot <- rowSums(dp)
    ref <- ref / tot * 10^6
    alt <- alt / tot * 10^6

    valid_alt <- alt[valid_scan, valid_snp]
    valid_alt[valid_alt == 0] <- NA
    expect_equal(as.integer(getQtileReadAlt(gds, "snp", valid = TRUE)),
                 as.integer(apply(valid_alt, 2, quantile, probs = q, na.rm = TRUE)),
                 ignore_attr = TRUE)

    expect_equal(as.integer(getQtileReadAlt(gds, "scan", valid = TRUE)),
                 as.integer(apply(valid_alt, 1, quantile, probs = q, na.rm = TRUE)),
                 ignore_attr = TRUE)
})
