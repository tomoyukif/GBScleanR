library(GBScleanR)

vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, TRUE, FALSE)
on.exit({unlink(gds_fn)})
gds <- loadGDS(gds_fn)
gds <- countGenotype(gds)

geno <- getGenotype(gds)

test_that("Reference genotype", {
    expect_equal(getCountGenoRef(gds, "snp"), colSums(geno == 2, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "scan"), rowSums(geno == 2, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "snp", prop = TRUE),
                 colSums(geno == 2, na.rm = TRUE) / colSums(!is.na(geno)),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "scan", prop = TRUE),
                 rowSums(geno == 2, na.rm = TRUE)/ rowSums(!is.na(geno)),
                 ignore_attr = TRUE)
    
    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- countGenotype(gds)
    
    expect_equal(getCountGenoRef(gds, "snp", valid = TRUE), 
                 colSums(geno[valid_scan, ] == 2, na.rm = TRUE)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "scan", valid = TRUE),
                 rowSums(geno[, valid_snp] == 2, na.rm = TRUE)[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "snp", valid = TRUE, prop = TRUE),
                 (colSums(geno[valid_scan, ] == 2, na.rm = TRUE) /
                     colSums(!is.na(geno[valid_scan, ])))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoRef(gds, "scan", valid = TRUE, prop = TRUE),
                 (rowSums(geno[, valid_snp] == 2, na.rm = TRUE)/ 
                     rowSums(!is.na(geno[, valid_snp])))[valid_scan],
                 ignore_attr = TRUE)
})

test_that("Heterozygous genotype", {
    expect_equal(getCountGenoHet(gds, "snp"), colSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "scan"), rowSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "snp", prop = TRUE),
                 colSums(geno == 1, na.rm = TRUE) / colSums(!is.na(geno)),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "scan", prop = TRUE),
                 rowSums(geno == 1, na.rm = TRUE)/ rowSums(!is.na(geno)),
                 ignore_attr = TRUE)
    
    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- countGenotype(gds)
    
    expect_equal(getCountGenoHet(gds, "snp", valid = TRUE), 
                 colSums(geno[valid_scan, ] == 1, na.rm = TRUE)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "scan", valid = TRUE),
                 rowSums(geno[, valid_snp] == 1, na.rm = TRUE)[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "snp", valid = TRUE, prop = TRUE),
                 (colSums(geno[valid_scan, ] == 1, na.rm = TRUE) /
                      colSums(!is.na(geno[valid_scan, ])))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoHet(gds, "scan", valid = TRUE, prop = TRUE),
                 (rowSums(geno[, valid_snp] == 1, na.rm = TRUE)/ 
                      rowSums(!is.na(geno[, valid_snp])))[valid_scan],
                 ignore_attr = TRUE)
})


test_that("Alternative genotype", {
    expect_equal(getCountGenoAlt(gds, "snp"), colSums(geno == 0, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "scan"), rowSums(geno == 0, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "snp", prop = TRUE),
                 colSums(geno == 0, na.rm = TRUE) / colSums(!is.na(geno)),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "scan", prop = TRUE),
                 rowSums(geno == 0, na.rm = TRUE)/ rowSums(!is.na(geno)),
                 ignore_attr = TRUE)
    
    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- countGenotype(gds)
    
    expect_equal(getCountGenoAlt(gds, "snp", valid = TRUE), 
                 colSums(geno[valid_scan, ] == 0, na.rm = TRUE)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "scan", valid = TRUE),
                 rowSums(geno[, valid_snp] == 0, na.rm = TRUE)[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "snp", valid = TRUE, prop = TRUE),
                 (colSums(geno[valid_scan, ] == 0, na.rm = TRUE) /
                      colSums(!is.na(geno[valid_scan, ])))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoAlt(gds, "scan", valid = TRUE, prop = TRUE),
                 (rowSums(geno[, valid_snp] == 0, na.rm = TRUE)/ 
                      rowSums(!is.na(geno[, valid_snp])))[valid_scan],
                 ignore_attr = TRUE)
})

test_that("Missing genotype", {
    expect_equal(getCountGenoMissing(gds, "snp"), colSums(is.na(geno)),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "scan"), rowSums(is.na(geno)),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "snp", prop = TRUE),
                 colSums(is.na(geno)) / nrow(geno),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "scan", prop = TRUE),
                 rowSums(is.na(geno))/ ncol(geno),
                 ignore_attr = TRUE)
    
    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- countGenotype(gds)
    
    expect_equal(getCountGenoMissing(gds, "snp"), 
                 colSums(is.na(geno[valid_scan, ]))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "scan"), 
                 rowSums(is.na(geno[, valid_snp]))[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "snp", prop = TRUE),
                 colSums(is.na(geno[valid_scan, ]))[valid_snp] / 
                     nrow(geno[valid_scan, ]),
                 ignore_attr = TRUE)
    expect_equal(getCountGenoMissing(gds, "scan", prop = TRUE),
                 rowSums(is.na(geno[, valid_snp]))[valid_scan] / 
                     ncol(geno[, valid_snp]),
                 ignore_attr = TRUE)
})

test_that("Reference allele", {
    expect_equal(getCountAlleleRef(gds, "snp"), 
                 colSums(geno == 2, na.rm = TRUE) * 2 +
                     colSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "scan"), 
                 rowSums(geno == 2, na.rm = TRUE) * 2 + 
                     rowSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "snp", prop = TRUE),
                 (colSums(geno == 2, na.rm = TRUE) * 2 +
                     colSums(geno == 1, na.rm = TRUE)) / colSums(!is.na(geno)) / 2,
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "scan", prop = TRUE),
                 (rowSums(geno == 2, na.rm = TRUE) * 2 + 
                     rowSums(geno == 1, na.rm = TRUE)) / rowSums(!is.na(geno)) / 2,
                 ignore_attr = TRUE)
    
    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- countGenotype(gds)
    
    expect_equal(getCountAlleleRef(gds, "snp", valid = TRUE), 
                 (colSums(geno[valid_scan, ] == 2, na.rm = TRUE) * 2 +
                      colSums(geno[valid_scan, ] == 1, na.rm = TRUE))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "scan", valid = TRUE), 
                 (rowSums(geno[, valid_snp] == 2, na.rm = TRUE) * 2 + 
                      rowSums(geno[, valid_snp] == 1, na.rm = TRUE))[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "snp", prop = TRUE),
                 ((colSums(geno[valid_scan, ] == 2, na.rm = TRUE) * 2 +
                      colSums(geno[valid_scan, ] == 1, na.rm = TRUE)) /
                     colSums(!is.na(geno[valid_scan, ])) / 2)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleRef(gds, "scan", prop = TRUE),
                 ((rowSums(geno[, valid_snp] == 2, na.rm = TRUE) * 2 + 
                      rowSums(geno[, valid_snp] == 1, na.rm = TRUE)) / 
                     rowSums(!is.na(geno[, valid_snp])) / 2)[valid_scan],
                 ignore_attr = TRUE)
})


test_that("Alternative allele", {
    expect_equal(getCountAlleleAlt(gds, "snp"), 
                 colSums(geno == 0, na.rm = TRUE) * 2 +
                     colSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "scan"), 
                 rowSums(geno == 0, na.rm = TRUE) * 2 + 
                     rowSums(geno == 1, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "snp", prop = TRUE),
                 (colSums(geno == 0, na.rm = TRUE) * 2 +
                      colSums(geno == 1, na.rm = TRUE)) / colSums(!is.na(geno)) / 2,
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "scan", prop = TRUE),
                 (rowSums(geno == 0, na.rm = TRUE) * 2 + 
                      rowSums(geno == 1, na.rm = TRUE)) / rowSums(!is.na(geno)) / 2,
                 ignore_attr = TRUE)
    
    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- countGenotype(gds)
    
    expect_equal(getCountAlleleAlt(gds, "snp", valid = TRUE), 
                 (colSums(geno[valid_scan, ] == 0, na.rm = TRUE) * 2 +
                      colSums(geno[valid_scan, ] == 1, na.rm = TRUE))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "scan", valid = TRUE), 
                 (rowSums(geno[, valid_snp] == 0, na.rm = TRUE) * 2 + 
                      rowSums(geno[, valid_snp] == 1, na.rm = TRUE))[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "snp", prop = TRUE),
                 ((colSums(geno[valid_scan, ] == 0, na.rm = TRUE) * 2 +
                       colSums(geno[valid_scan, ] == 1, na.rm = TRUE)) /
                      colSums(!is.na(geno[valid_scan, ])) / 2)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleAlt(gds, "scan", prop = TRUE),
                 ((rowSums(geno[, valid_snp] == 0, na.rm = TRUE) * 2 + 
                       rowSums(geno[, valid_snp] == 1, na.rm = TRUE)) / 
                      rowSums(!is.na(geno[, valid_snp])) / 2)[valid_scan],
                 ignore_attr = TRUE)
})


test_that("Missing allele", {
    expect_equal(getCountAlleleMissing(gds, "snp"), 
                 colSums(is.na(geno)) * 2,
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "scan"), 
                 rowSums(is.na(geno)) * 2,
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "snp", prop = TRUE),
                 (colSums(is.na(geno)) * 2) / (nrow(geno) * 2),
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "scan", prop = TRUE),
                 (rowSums(is.na(geno)) * 2) / (ncol(geno) * 2),
                 ignore_attr = TRUE)
    
    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- countGenotype(gds)
    
    expect_equal(getCountAlleleMissing(gds, "snp", valid = TRUE), 
                 (colSums(is.na(geno[valid_scan, ])) * 2)[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "scan", valid = TRUE), 
                 (rowSums(is.na(geno[, valid_snp])) * 2)[valid_scan],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "snp", prop = TRUE),
                 ((colSums(is.na(geno[valid_scan, ])) * 2) /
                     (nrow(geno[valid_scan, ]) * 2))[valid_snp],
                 ignore_attr = TRUE)
    expect_equal(getCountAlleleMissing(gds, "scan", prop = TRUE),
                 ((rowSums(is.na(geno[, valid_snp])) * 2) / 
                     (ncol(geno[, valid_snp]) * 2))[valid_scan],
                 ignore_attr = TRUE)
})

test_that("Minor allele frequency", {
    ref_f <- colSums(geno == 0, na.rm = TRUE) * 2 + colSums(geno == 1, na.rm = TRUE)
    alt_f <- colSums(geno == 2, na.rm = TRUE) * 2 + colSums(geno == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    maf <- ref_f
    maf[alt_minor] <- alt_f[alt_minor]
    maf <- maf / (ref_f + alt_f)
    expect_equal(getMAF(gds, "snp"), maf, ignore_attr = TRUE)
    
    ref_f <- rowSums(geno == 0, na.rm = TRUE) * 2 + rowSums(geno == 1, na.rm = TRUE)
    alt_f <- rowSums(geno == 2, na.rm = TRUE) * 2 + rowSums(geno == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    maf <- ref_f
    maf[alt_minor] <- alt_f[alt_minor]
    maf <- maf / (ref_f + alt_f)
    expect_equal(getMAF(gds, "scan"), maf, ignore_attr = TRUE)
    
    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- countGenotype(gds)
    
    ref_f <- colSums(geno[valid_scan, ] == 0, na.rm = TRUE) * 2 +
        colSums(geno[valid_scan, ] == 1, na.rm = TRUE)
    alt_f <- colSums(geno[valid_scan, ] == 2, na.rm = TRUE) * 2 +
        colSums(geno[valid_scan, ] == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    maf <- ref_f
    maf[alt_minor] <- alt_f[alt_minor]
    maf <- maf / (ref_f + alt_f)
    expect_equal(getMAF(gds, "snp", valid = TRUE),
                 maf[valid_snp], ignore_attr = TRUE)
    
    ref_f <- rowSums(geno[, valid_snp] == 0, na.rm = TRUE) * 2 +
        rowSums(geno[, valid_snp] == 1, na.rm = TRUE)
    alt_f <- rowSums(geno[, valid_snp] == 2, na.rm = TRUE) * 2 + 
        rowSums(geno[, valid_snp] == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    maf <- ref_f
    maf[alt_minor] <- alt_f[alt_minor]
    maf <- maf / (ref_f + alt_f)
    expect_equal(getMAF(gds, "scan", valid = TRUE), 
                 maf[valid_scan], ignore_attr = TRUE)
})


test_that("Minor allele count", {
    ref_f <- colSums(geno == 0, na.rm = TRUE) * 2 + colSums(geno == 1, na.rm = TRUE)
    alt_f <- colSums(geno == 2, na.rm = TRUE) * 2 + colSums(geno == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    mac <- ref_f
    mac[alt_minor] <- alt_f[alt_minor]
    expect_equal(getMAC(gds, "snp"), mac, ignore_attr = TRUE)
    
    ref_f <- rowSums(geno == 0, na.rm = TRUE) * 2 + rowSums(geno == 1, na.rm = TRUE)
    alt_f <- rowSums(geno == 2, na.rm = TRUE) * 2 + rowSums(geno == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    mac <- ref_f
    mac[alt_minor] <- alt_f[alt_minor]
    expect_equal(getMAC(gds, "scan"), mac, ignore_attr = TRUE)
    
    valid_snp <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, valid_snp)
    valid_scan <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, valid_scan)
    gds <- countGenotype(gds)
    
    ref_f <- colSums(geno[valid_scan, ] == 0, na.rm = TRUE) * 2 +
        colSums(geno[valid_scan, ] == 1, na.rm = TRUE)
    alt_f <- colSums(geno[valid_scan, ] == 2, na.rm = TRUE) * 2 +
        colSums(geno[valid_scan, ] == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    mac <- ref_f
    mac[alt_minor] <- alt_f[alt_minor]
    expect_equal(getMAC(gds, "snp", valid = TRUE),
                 mac[valid_snp], ignore_attr = TRUE)
    
    ref_f <- rowSums(geno[, valid_snp] == 0, na.rm = TRUE) * 2 +
        rowSums(geno[, valid_snp] == 1, na.rm = TRUE)
    alt_f <- rowSums(geno[, valid_snp] == 2, na.rm = TRUE) * 2 + 
        rowSums(geno[, valid_snp] == 1, na.rm = TRUE)
    alt_minor <- ref_f > alt_f
    mac <- ref_f
    mac[alt_minor] <- alt_f[alt_minor]
    expect_equal(getMAC(gds, "scan", valid = TRUE), 
                 mac[valid_scan], ignore_attr = TRUE)
})


