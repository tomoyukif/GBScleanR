vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, verbose = FALSE)
on.exit({unlink(gds_fn)})

test_that("Reference genotype", {
    gds <- loadGDS(gds_fn)
    gds <- countGenotype(gds)
    geno <- getGenotype(gds)
    expect_equal(getCountGenoRef(gds, "marker"), colSums(geno == 0, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "sample"), rowSums(geno == 0, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "marker", prop = TRUE),
                 colSums(geno == 0, na.rm = TRUE) / colSums(!is.na(geno)),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "sample", prop = TRUE),
                 rowSums(geno == 0, na.rm = TRUE)/ rowSums(!is.na(geno)),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    validMar(gds) <- valid_snp
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    validSam(gds) <- valid_scan
    gds <- countGenotype(gds)

    expect_equal(getCountGenoRef(gds, "marker", valid = TRUE),
                 colSums(geno[valid_scan, ] == 0, na.rm = TRUE)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "sample", valid = TRUE),
                 rowSums(geno[, valid_snp] == 0, na.rm = TRUE)[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "marker", valid = TRUE, prop = TRUE),
                 (colSums(geno[valid_scan, ] == 0, na.rm = TRUE) /
                     colSums(!is.na(geno[valid_scan, ])))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "sample", valid = TRUE, prop = TRUE),
                 (rowSums(geno[, valid_snp] == 0, na.rm = TRUE)/
                     rowSums(!is.na(geno[, valid_snp])))[valid_scan],
                 ignore_attr = TRUE)
    closeGDS(gds)
})

test_that("Heterozygous genotype", {
    gds <- loadGDS(gds_fn)
    gds <- countGenotype(gds)
    geno <- getGenotype(gds)
    expect_equal(getCountGenoHet(gds, "marker"), colSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "sample"), rowSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "marker", prop = TRUE),
                 colSums(geno == 1, na.rm = TRUE) / colSums(!is.na(geno)),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "sample", prop = TRUE),
                 rowSums(geno == 1, na.rm = TRUE)/ rowSums(!is.na(geno)),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    validMar(gds) <- valid_snp
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    validSam(gds) <- valid_scan
    gds <- countGenotype(gds)

    expect_equal(getCountGenoHet(gds, "marker", valid = TRUE),
                 colSums(geno[valid_scan, ] == 1, na.rm = TRUE)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "sample", valid = TRUE),
                 rowSums(geno[, valid_snp] == 1, na.rm = TRUE)[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "marker", valid = TRUE, prop = TRUE),
                 (colSums(geno[valid_scan, ] == 1, na.rm = TRUE) /
                      colSums(!is.na(geno[valid_scan, ])))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "sample", valid = TRUE, prop = TRUE),
                 (rowSums(geno[, valid_snp] == 1, na.rm = TRUE)/
                      rowSums(!is.na(geno[, valid_snp])))[valid_scan],
                 ignore_attr = TRUE)
    closeGDS(gds)
})


test_that("Alternative genotype", {
    gds <- loadGDS(gds_fn)
    gds <- countGenotype(gds)
    geno <- getGenotype(gds)
    expect_equal(getCountGenoAlt(gds, "marker"), colSums(geno == 2, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "sample"), rowSums(geno == 2, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "marker", prop = TRUE),
                 colSums(geno == 2, na.rm = TRUE) / colSums(!is.na(geno)),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "sample", prop = TRUE),
                 rowSums(geno == 2, na.rm = TRUE)/ rowSums(!is.na(geno)),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    validMar(gds) <- valid_snp
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    validSam(gds) <- valid_scan
    gds <- countGenotype(gds)

    expect_equal(getCountGenoAlt(gds, "marker", valid = TRUE),
                 colSums(geno[valid_scan, ] == 2, na.rm = TRUE)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "sample", valid = TRUE),
                 rowSums(geno[, valid_snp] == 2, na.rm = TRUE)[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "marker", valid = TRUE, prop = TRUE),
                 (colSums(geno[valid_scan, ] == 2, na.rm = TRUE) /
                      colSums(!is.na(geno[valid_scan, ])))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "sample", valid = TRUE, prop = TRUE),
                 (rowSums(geno[, valid_snp] == 2, na.rm = TRUE)/
                      rowSums(!is.na(geno[, valid_snp])))[valid_scan],
                 ignore_attr = TRUE)
    closeGDS(gds)
})

test_that("Missing genotype", {
    gds <- loadGDS(gds_fn)
    gds <- countGenotype(gds)
    geno <- getGenotype(gds)
    expect_equal(getCountGenoMissing(gds, "marker"), colSums(is.na(geno)),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "sample"), rowSums(is.na(geno)),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "marker", prop = TRUE),
                 colSums(is.na(geno)) / nrow(geno),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "sample", prop = TRUE),
                 rowSums(is.na(geno))/ ncol(geno),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    validMar(gds) <- valid_snp
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    validSam(gds) <- valid_scan
    gds <- countGenotype(gds)

    expect_equal(getCountGenoMissing(gds, "marker"),
                 colSums(is.na(geno[valid_scan, ]))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "sample"),
                 rowSums(is.na(geno[, valid_snp]))[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "marker", prop = TRUE),
                 colSums(is.na(geno[valid_scan, ]))[valid_snp] /
                     nrow(geno[valid_scan, ]),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "sample", prop = TRUE),
                 rowSums(is.na(geno[, valid_snp]))[valid_scan] /
                     ncol(geno[, valid_snp]),
                 ignore_attr = TRUE)
    closeGDS(gds)
})

test_that("Reference allele", {
    gds <- loadGDS(gds_fn)
    gds <- countGenotype(gds)
    geno <- getGenotype(gds)
    expect_equal(getCountAlleleRef(gds, "marker"),
                 colSums(geno == 0, na.rm = TRUE) * 2 +
                     colSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "sample"),
                 rowSums(geno == 0, na.rm = TRUE) * 2 +
                     rowSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "marker", prop = TRUE),
                 (colSums(geno == 0, na.rm = TRUE) * 2 +
                     colSums(geno == 1, na.rm = TRUE)) / colSums(!is.na(geno)) / 2,
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "sample", prop = TRUE),
                 (rowSums(geno == 0, na.rm = TRUE) * 2 +
                     rowSums(geno == 1, na.rm = TRUE)) / rowSums(!is.na(geno)) / 2,
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    validMar(gds) <- valid_snp
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    validSam(gds) <- valid_scan
    gds <- countGenotype(gds)

    expect_equal(getCountAlleleRef(gds, "marker", valid = TRUE),
                 (colSums(geno[valid_scan, ] == 0, na.rm = TRUE) * 2 +
                      colSums(geno[valid_scan, ] == 1, na.rm = TRUE))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "sample", valid = TRUE),
                 (rowSums(geno[, valid_snp] == 0, na.rm = TRUE) * 2 +
                      rowSums(geno[, valid_snp] == 1, na.rm = TRUE))[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "marker", prop = TRUE),
                 ((colSums(geno[valid_scan, ] == 0, na.rm = TRUE) * 2 +
                      colSums(geno[valid_scan, ] == 1, na.rm = TRUE)) /
                     colSums(!is.na(geno[valid_scan, ])) / 2)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "sample", prop = TRUE),
                 ((rowSums(geno[, valid_snp] == 0, na.rm = TRUE) * 2 +
                      rowSums(geno[, valid_snp] == 1, na.rm = TRUE)) /
                     rowSums(!is.na(geno[, valid_snp])) / 2)[valid_scan],
                 ignore_attr = TRUE)
    closeGDS(gds)
})


test_that("Alternative allele", {
    gds <- loadGDS(gds_fn)
    gds <- countGenotype(gds)
    geno <- getGenotype(gds)
    expect_equal(getCountAlleleAlt(gds, "marker"),
                 colSums(geno == 2, na.rm = TRUE) * 2 +
                     colSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "sample"),
                 rowSums(geno == 2, na.rm = TRUE) * 2 +
                     rowSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "marker", prop = TRUE),
                 (colSums(geno == 2, na.rm = TRUE) * 2 +
                      colSums(geno == 1, na.rm = TRUE)) / colSums(!is.na(geno)) / 2,
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "sample", prop = TRUE),
                 (rowSums(geno == 2, na.rm = TRUE) * 2 +
                      rowSums(geno == 1, na.rm = TRUE)) / rowSums(!is.na(geno)) / 2,
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    validMar(gds) <- valid_snp
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    validSam(gds) <- valid_scan
    gds <- countGenotype(gds)

    expect_equal(getCountAlleleAlt(gds, "marker", valid = TRUE),
                 (colSums(geno[valid_scan, ] == 2, na.rm = TRUE) * 2 +
                      colSums(geno[valid_scan, ] == 1, na.rm = TRUE))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "sample", valid = TRUE),
                 (rowSums(geno[, valid_snp] == 2, na.rm = TRUE) * 2 +
                      rowSums(geno[, valid_snp] == 1, na.rm = TRUE))[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "marker", prop = TRUE),
                 ((colSums(geno[valid_scan, ] == 2, na.rm = TRUE) * 2 +
                       colSums(geno[valid_scan, ] == 1, na.rm = TRUE)) /
                      colSums(!is.na(geno[valid_scan, ])) / 2)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "sample", prop = TRUE),
                 ((rowSums(geno[, valid_snp] == 2, na.rm = TRUE) * 2 +
                       rowSums(geno[, valid_snp] == 1, na.rm = TRUE)) /
                      rowSums(!is.na(geno[, valid_snp])) / 2)[valid_scan],
                 ignore_attr = TRUE)
    closeGDS(gds)
})


test_that("Missing allele", {
    gds <- loadGDS(gds_fn)
    gds <- countGenotype(gds)
    geno <- getGenotype(gds)
    expect_equal(getCountAlleleMissing(gds, "marker"),
                 colSums(is.na(geno)) * 2,
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "sample"),
                 rowSums(is.na(geno)) * 2,
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "marker", prop = TRUE),
                 (colSums(is.na(geno)) * 2) / (nrow(geno) * 2),
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "sample", prop = TRUE),
                 (rowSums(is.na(geno)) * 2) / (ncol(geno) * 2),
                 ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    validMar(gds) <- valid_snp
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    validSam(gds) <- valid_scan
    gds <- countGenotype(gds)

    expect_equal(getCountAlleleMissing(gds, "marker", valid = TRUE),
                 (colSums(is.na(geno[valid_scan, ])) * 2)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "sample", valid = TRUE),
                 (rowSums(is.na(geno[, valid_snp])) * 2)[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "marker", prop = TRUE),
                 ((colSums(is.na(geno[valid_scan, ])) * 2) /
                     (nrow(geno[valid_scan, ]) * 2))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "sample", prop = TRUE),
                 ((rowSums(is.na(geno[, valid_snp])) * 2) /
                     (ncol(geno[, valid_snp]) * 2))[valid_scan],
                 ignore_attr = TRUE)
    closeGDS(gds)
})

test_that("Minor allele frequency", {
    gds <- loadGDS(gds_fn)
    gds <- countGenotype(gds)
    geno <- getGenotype(gds)
    ref_f <- colSums(geno == 2, na.rm = TRUE) * 2 + colSums(geno == 1, na.rm = TRUE)
    alt_f <- colSums(geno == 0, na.rm = TRUE) * 2 + colSums(geno == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    maf <- ref_f
    maf[alt_minor] <- alt_f[alt_minor]
    maf <- maf / (ref_f + alt_f)
    expect_equal(getMAF(gds, "marker"), maf, ignore_attr = TRUE)

    ref_f <- rowSums(geno == 2, na.rm = TRUE) * 2 + rowSums(geno == 1, na.rm = TRUE)
    alt_f <- rowSums(geno == 0, na.rm = TRUE) * 2 + rowSums(geno == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    maf <- ref_f
    maf[alt_minor] <- alt_f[alt_minor]
    maf <- maf / (ref_f + alt_f)
    expect_equal(getMAF(gds, "sample"), maf, ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    validMar(gds) <- valid_snp
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    validSam(gds) <- valid_scan
    gds <- countGenotype(gds)

    ref_f <- colSums(geno[valid_scan, ] == 2, na.rm = TRUE) * 2 +
        colSums(geno[valid_scan, ] == 1, na.rm = TRUE)
    alt_f <- colSums(geno[valid_scan, ] == 0, na.rm = TRUE) * 2 +
        colSums(geno[valid_scan, ] == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    maf <- ref_f
    maf[alt_minor] <- alt_f[alt_minor]
    maf <- maf / (ref_f + alt_f)
    expect_equal(getMAF(gds, "marker", valid = TRUE),
                 maf[valid_snp], ignore_attr = TRUE)

    ref_f <- rowSums(geno[, valid_snp] == 2, na.rm = TRUE) * 2 +
        rowSums(geno[, valid_snp] == 1, na.rm = TRUE)
    alt_f <- rowSums(geno[, valid_snp] == 0, na.rm = TRUE) * 2 +
        rowSums(geno[, valid_snp] == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    maf <- ref_f
    maf[alt_minor] <- alt_f[alt_minor]
    maf <- maf / (ref_f + alt_f)
    expect_equal(getMAF(gds, "sample", valid = TRUE),
                 maf[valid_scan], ignore_attr = TRUE)
    closeGDS(gds)
})


test_that("Minor allele count", {
    gds <- loadGDS(gds_fn)
    gds <- countGenotype(gds)
    geno <- getGenotype(gds)
    ref_f <- colSums(geno == 2, na.rm = TRUE) * 2 + colSums(geno == 1, na.rm = TRUE)
    alt_f <- colSums(geno == 0, na.rm = TRUE) * 2 + colSums(geno == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    mac <- ref_f
    mac[alt_minor] <- alt_f[alt_minor]
    expect_equal(getMAC(gds, "marker"), mac, ignore_attr = TRUE)

    ref_f <- rowSums(geno == 2, na.rm = TRUE) * 2 + rowSums(geno == 1, na.rm = TRUE)
    alt_f <- rowSums(geno == 0, na.rm = TRUE) * 2 + rowSums(geno == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    mac <- ref_f
    mac[alt_minor] <- alt_f[alt_minor]
    expect_equal(getMAC(gds, "sample"), mac, ignore_attr = TRUE)

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    validMar(gds) <- valid_snp
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    validSam(gds) <- valid_scan
    gds <- countGenotype(gds)

    ref_f <- colSums(geno[valid_scan, ] == 2, na.rm = TRUE) * 2 +
        colSums(geno[valid_scan, ] == 1, na.rm = TRUE)
    alt_f <- colSums(geno[valid_scan, ] == 0, na.rm = TRUE) * 2 +
        colSums(geno[valid_scan, ] == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    mac <- ref_f
    mac[alt_minor] <- alt_f[alt_minor]
    expect_equal(getMAC(gds, "marker", valid = TRUE),
                 mac[valid_snp], ignore_attr = TRUE)

    ref_f <- rowSums(geno[, valid_snp] == 2, na.rm = TRUE) * 2 +
        rowSums(geno[, valid_snp] == 1, na.rm = TRUE)
    alt_f <- rowSums(geno[, valid_snp] == 0, na.rm = TRUE) * 2 +
        rowSums(geno[, valid_snp] == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    mac <- ref_f
    mac[alt_minor] <- alt_f[alt_minor]
    expect_equal(getMAC(gds, "sample", valid = TRUE),
                 mac[valid_scan], ignore_attr = TRUE)
    closeGDS(gds)
})


