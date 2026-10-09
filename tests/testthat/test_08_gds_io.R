.make_temp_gds <- function() {
    vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
    gds_fn <- tempfile("sample", fileext = ".gds")
    gbsrVCF2GDS(vcf_fn, gds_fn, force = TRUE, verbose = FALSE)
    gds_fn
}

test_that("gbsrVCF2GDS creates a loadable GDS", {
    gds_fn <- .make_temp_gds()
    on.exit(unlink(gds_fn))
    expect_true(file.exists(gds_fn))
    gds <- loadGDS(gds_fn, verbose = FALSE)
    on.exit({
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
        unlink(gds_fn)
    })
    expect_s4_class(gds, "GbsrGenotypeData")
    expect_gt(nsam(gds), 0)
    expect_gt(nmar(gds), 0)
    expect_true(isOpenGDS(gds))
})

test_that("force overwrite of existing GDS", {
    gds_fn <- .make_temp_gds()
    on.exit(unlink(gds_fn))
    out <- gbsrVCF2GDS(
        system.file("extdata", "sample.vcf", package = "GBScleanR"),
        gds_fn, force = TRUE, verbose = FALSE
    )
    expect_equal(out, gds_fn)
    expect_true(file.exists(gds_fn))
})

test_that("loadGDS rejects invalid input", {
    expect_error(loadGDS(1L), "x must be a file path")
})

test_that("closeGDS, reopenGDS, and isOpenGDS", {
    gds_fn <- .make_temp_gds()
    gds <- loadGDS(gds_fn, verbose = FALSE)
    on.exit({
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
        unlink(gds_fn)
    })
    expect_true(isOpenGDS(gds))
    closeGDS(gds, verbose = FALSE)
    expect_false(isOpenGDS(gds))
    gds <- reopenGDS(gds)
    expect_true(isOpenGDS(gds))
    expect_equal(nsam(gds), length(SeqArray::seqGetData(gds, "sample.id")))
})

test_that("reload GDS from GbsrGenotypeData object", {
    gds_fn <- .make_temp_gds()
    gds <- loadGDS(gds_fn, verbose = FALSE)
    on.exit({
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
        unlink(gds_fn)
    })
    n_sam <- nsam(gds)
    n_mar <- nmar(gds)
    gds <- loadGDS(gds, verbose = FALSE)
    expect_equal(nsam(gds), n_sam)
    expect_equal(nmar(gds), n_mar)
})

test_that("getPloidy and setPloidy", {
    gds_fn <- .make_temp_gds()
    gds <- loadGDS(gds_fn, ploidy = 2, verbose = FALSE)
    on.exit({
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
        unlink(gds_fn)
    })
    expect_equal(getPloidy(gds), 2)
    gds <- setPloidy(gds, ploidy = 4)
    expect_equal(getPloidy(gds), 4)
})

test_that("show/print do not error", {
    gds_fn <- .make_temp_gds()
    gds <- loadGDS(gds_fn, verbose = FALSE)
    on.exit({
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
        unlink(gds_fn)
    })
    expect_no_error(capture.output(show(gds)))
    expect_no_error(capture.output(print(gds)))
})
