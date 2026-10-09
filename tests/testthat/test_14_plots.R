vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, verbose = FALSE)
on.exit({unlink(gds_fn)})

test_that("histGBSR without data returns stats catalogue", {
    info <- histGBSR()
    expect_s3_class(info, "data.frame")
    expect_true(all(c("stats", "description") %in% names(info)))
    expect_true("missing" %in% info$stats)
})

test_that("QC plot helpers return ggplot objects", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    gds <- countGenotype(gds)
    gds <- countRead(gds)

    expect_s3_class(histGBSR(gds, stats = "missing"), "ggplot")
    expect_s3_class(histGBSR(gds, stats = "het"), "ggplot")
    expect_s3_class(histGBSR(gds, stats = "raf"), "ggplot")
    expect_s3_class(histGBSR(gds, stats = "dp"), "ggplot")

    expect_s3_class(plotGBSR(gds, stats = "missing"), "ggplot")
    expect_s3_class(plotGBSR(gds, stats = "raf"), "ggplot")
    expect_s3_class(plotGBSR(gds, stats = "geno"), "ggplot")
    expect_s3_class(plotGBSR(gds, stats = "marker"), "ggplot")

    expect_s3_class(boxplotGBSR(gds, stats = "missing"), "ggplot")
    expect_s3_class(pairsGBSR(gds, stats1 = "missing", stats2 = "het"),
                    "ggplot")
    closeGDS(gds, verbose = FALSE)
})
