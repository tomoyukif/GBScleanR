library(GBScleanR)

vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, verbose = FALSE)
on.exit({unlink(gds_fn)})

test_that("Set valid markers", {
    gds <- loadGDS(gds_fn)
    all_true <- validMar(gds)
    expect_true(all(all_true))
    rundom_true <- sample(c(TRUE, FALSE), nmar(gds, FALSE), TRUE)
    validMar(gds) <- rundom_true
    expect_equal(validMar(gds), rundom_true)
    expect_equal(validMar(gds), rundom_true)
    closeGDS(gds)
})

test_that("Set valid samples", {
    gds <- loadGDS(gds_fn)
    all_true <- validSam(gds)
    expect_true(all(all_true))
    rundom_true <- sample(c(TRUE, FALSE), nsam(gds, FALSE), TRUE)
    validSam(gds) <- rundom_true
    expect_equal(validSam(gds), rundom_true)
    expect_equal(validSam(gds), rundom_true)
    closeGDS(gds)
})

test_that("Set parents", {
    gds <- loadGDS(gds_fn)
    parents <- c("Founder1", "Founder2")
    gds <- setParents(gds, parents)
    p_bool <- getParents(gds, TRUE)
    p_info <- getParents(gds)
    expect_equal(validSam(gds), !p_bool)
    expect_equal(p_info$memberID, 1:2)
    expect_equal(p_info$sampleID, parents)

    gds <- setParents(gds, parents)
    p_bool <- getParents(gds, TRUE)
    p_info <- getParents(gds)
    expect_equal(validSam(gds), !p_bool)
    expect_equal(p_info$memberID, 1:2)
    expect_equal(p_info$sampleID, parents)
    closeGDS(gds)
})

test_that("Set parents with checking monomorphic marker",{
    gds <- loadGDS(gds_fn)
    parents <- c("Founder1", "Founder2")
    gds <- setParents(gds, parents)
    read <- getRead(gds, parents = "only")
    missing <- read$ref == 0 & read$alt == 0
    het <- read$ref > 0 & read$alt > 0
    mono <- !het & !missing
    mono <- colSums(mono) == sum(gds@sample$parents != 0)
    gds <- setParents(gds, parents, mono = TRUE, bi = FALSE)
    expect_equal(validMar(gds), mono, ignore_attr = TRUE)
    gds <- setParents(gds, parents, mono = TRUE, bi = FALSE)
    expect_equal(validMar(gds), mono, ignore_attr = TRUE)
    closeGDS(gds)
})

test_that("Set parents with checking biallelic marker",{
    gds <- loadGDS(gds_fn)
    parents <- c("Founder1", "Founder2")
    gds <- setParents(gds, parents)
    geno <- getGenotype(gds, parents = "only")
    read <- getRead(gds, parents = "only")
    missing <- read$ref == 0 & read$alt == 0
    het <- read$ref > 0 & read$alt > 0
    bi <- read$ref > 0 & read$alt == 0 & !het & !missing
    bi <- colSums(bi) != sum(gds@sample$parents != 0)
    gds <- setParents(gds, parents, mono = FALSE, bi = TRUE)
    expect_equal(validMar(gds), bi, ignore_attr = TRUE)
    gds <- setParents(gds, parents, mono = FALSE, bi = TRUE)
    expect_equal(validMar(gds), bi, ignore_attr = TRUE)
    closeGDS(gds)
})

test_that("Set parents with checking marker filtering",{
    gds <- loadGDS(gds_fn)
    parents <- c("Founder1", "Founder2")
    gds <- setParents(gds, parents)
    read <- getRead(gds, parents = "only")
    missing <- read$ref == 0 & read$alt == 0
    het <- read$ref > 0 & read$alt > 0
    mono <- !het & !missing
    mono <- colSums(mono) == sum(gds@sample$parents != 0)
    bi <- read$ref > 0 & read$alt == 0 & !het & !missing
    bi <- colSums(bi) != sum(gds@sample$parents != 0)

    gds <- setParents(gds, parents, mono = TRUE, bi = TRUE)
    expect_equal(validMar(gds), bi & mono, ignore_attr = TRUE)
    gds <- setParents(gds, parents, mono = TRUE, bi = TRUE)
    expect_equal(validMar(gds), bi & mono, ignore_attr = TRUE)
    closeGDS(gds)
})

