vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, verbose = FALSE)
on.exit({unlink(gds_fn)})

# Mean/SD/median are CPM of non-zero allele depths (see count_read in C++).
.cpm_sample_stat <- function(x_mat, other_mat, fun) {
    vapply(seq_len(nrow(x_mat)), function(i) {
        r <- x_mat[i, ]
        a <- other_mat[i, ]
        r[r == 0] <- NA_real_
        a[a == 0] <- NA_real_
        r2 <- stats::na.omit(r)
        a2 <- stats::na.omit(a)
        dp <- sum(r2) + sum(a2)
        fun(r2 / dp * 1e6)
    }, numeric(1))
}

.cpm_marker_stat <- function(x_mat, tot_read, fun) {
    vapply(seq_len(ncol(x_mat)), function(j) {
        v <- x_mat[, j]
        v[v == 0] <- NA_real_
        v <- v / tot_read * 1e6
        fun(stats::na.omit(v))
    }, numeric(1))
}

test_that("read summary getters error before countRead", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    expect_error(getMeanReadRef(gds), "countRead")
    expect_error(getMeanReadAlt(gds), "countRead")
    expect_error(getSDReadRef(gds), "countRead")
    expect_error(getSDReadAlt(gds), "countRead")
    expect_error(getMedianReadRef(gds), "countRead")
    expect_error(getMedianReadAlt(gds), "countRead")
    closeGDS(gds, verbose = FALSE)
})

test_that("CPM mean/SD/median of reference and alternative reads", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    gds <- countRead(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    tot <- rowSums(ref) + rowSums(alt)

    expect_equal(as.numeric(getMeanReadRef(gds, "sample")),
                 .cpm_sample_stat(ref, alt, mean),
                 ignore_attr = TRUE)
    expect_equal(as.numeric(getMeanReadAlt(gds, "sample")),
                 .cpm_sample_stat(alt, ref, mean),
                 ignore_attr = TRUE)
    expect_equal(as.numeric(getMeanReadRef(gds, "marker")),
                 .cpm_marker_stat(ref, tot, mean),
                 ignore_attr = TRUE)
    expect_equal(as.numeric(getMeanReadAlt(gds, "marker")),
                 .cpm_marker_stat(alt, tot, mean),
                 ignore_attr = TRUE)

    expect_equal(as.numeric(getSDReadRef(gds, "sample")),
                 .cpm_sample_stat(ref, alt, stats::sd),
                 ignore_attr = TRUE)
    expect_equal(as.numeric(getSDReadAlt(gds, "sample")),
                 .cpm_sample_stat(alt, ref, stats::sd),
                 ignore_attr = TRUE)
    expect_equal(as.numeric(getSDReadRef(gds, "marker")),
                 .cpm_marker_stat(ref, tot, stats::sd),
                 ignore_attr = TRUE)
    expect_equal(as.numeric(getSDReadAlt(gds, "marker")),
                 .cpm_marker_stat(alt, tot, stats::sd),
                 ignore_attr = TRUE)

    expect_equal(as.numeric(getMedianReadRef(gds, "sample")),
                 .cpm_sample_stat(ref, alt, stats::median),
                 ignore_attr = TRUE)
    expect_equal(as.numeric(getMedianReadAlt(gds, "sample")),
                 .cpm_sample_stat(alt, ref, stats::median),
                 ignore_attr = TRUE)
    expect_equal(as.numeric(getMedianReadRef(gds, "marker")),
                 .cpm_marker_stat(ref, tot, stats::median),
                 ignore_attr = TRUE)
    expect_equal(as.numeric(getMedianReadAlt(gds, "marker")),
                 .cpm_marker_stat(alt, tot, stats::median),
                 ignore_attr = TRUE)
    closeGDS(gds, verbose = FALSE)
})

test_that("read stats respect marker/sample validity", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    validMar(gds) <- valid_snp
    validSam(gds) <- valid_scan
    gds <- countRead(gds)

    read <- getRead(gds, valid = FALSE)
    ref <- read$ref[valid_scan, valid_snp, drop = FALSE]
    alt <- read$alt[valid_scan, valid_snp, drop = FALSE]
    tot <- rowSums(ref) + rowSums(alt)

    expect_equal(as.numeric(getMeanReadRef(gds, "marker", valid = TRUE)),
                 .cpm_marker_stat(ref, tot, mean),
                 ignore_attr = TRUE)
    expect_equal(as.numeric(getMeanReadRef(gds, "sample", valid = TRUE)),
                 .cpm_sample_stat(ref, alt, mean),
                 ignore_attr = TRUE)
    closeGDS(gds, verbose = FALSE)
})
