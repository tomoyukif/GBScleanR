library(GBScleanR)

vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, TRUE, FALSE)
gds <- loadGDS(gds_fn)
on.exit({unlink(gds_fn)})

test_that("Set valid markers", {
    all_true <- getValidSnp(gds)
    expect_true(all(all_true))
    rundom_true <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, new = rundom_true)
    expect_equal(getValidSnp(gds), rundom_true)
    update_true <- sample(c(TRUE, FALSE), nsnp(gds), TRUE)
    gds <- setValidSnp(gds, update = update_true)
    rundom_true[rundom_true][!update_true] <- FALSE
    expect_equal(getValidSnp(gds), rundom_true)
})

test_that("Set valid samples", {
    all_true <- getValidScan(gds)
    expect_true(all(all_true))
    rundom_true <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, new = rundom_true)
    expect_equal(getValidScan(gds), rundom_true)
    update_true <- sample(c(TRUE, FALSE), nscan(gds), TRUE)
    gds <- setValidScan(gds, update = update_true)
    rundom_true[rundom_true][!update_true] <- FALSE
    expect_equal(getValidScan(gds), rundom_true)
})

test_that("Set parents", {
    parents <- c("Founder1", "Founder2")
    gds <- setParents(gds, parents)
    p_bool <- getParents(gds, TRUE)
    p_info <- getParents(gds)
    expect_equal(getValidScan(gds), !p_bool)
    expect_equal(p_info$memberID, 1:2)
    expect_equal(p_info$scanID, parents)
    
    gds <- setParents(gds, parents)
    p_bool <- getParents(gds, TRUE)
    p_info <- getParents(gds)
    expect_equal(getValidScan(gds), !p_bool)
    expect_equal(p_info$memberID, 1:2)
    expect_equal(p_info$scanID, parents)
})

test_that("Set parents with checking flipped marker",{
    parents <- c("Founder1", "Founder2")
    gds <- setParents(gds, parents)
    geno <- getGenotype(gds, parents = "only")
    
    gds <- setParents(gds, parents, flip = TRUE, mono = FALSE, bi = FALSE)
    expect_true(hasFlipped(gds))
    flipped <- getFlipped(gds)
    expect_equal(flipped, geno[1, ] == 0, ignore_attr = TRUE)
    gds <- setParents(gds, parents, flip = TRUE, mono = FALSE, bi = FALSE)
    flipped <- getFlipped(gds)
    expect_equal(flipped, geno[1, ] == 0, ignore_attr = TRUE)
})

test_that("Set parents with checking monomorphic marker",{
    parents <- c("Founder1", "Founder2")
    gds <- setParents(gds, parents)
    geno <- getGenotype(gds, parents = "only")
    mono <- geno[1, ] %in% c(0, 2) & geno[2, ] %in% c(0, 2)
    
    gds <- setParents(gds, parents, flip = FALSE, mono = TRUE, bi = FALSE)
    expect_equal(getValidSnp(gds), mono)
    gds <- setParents(gds, parents, flip = FALSE, mono = TRUE, bi = FALSE)
    expect_equal(getValidSnp(gds), mono)
})

test_that("Set parents with checking biallelic marker",{
    parents <- c("Founder1", "Founder2")
    gds <- setParents(gds, parents)
    geno <- getGenotype(gds, parents = "only")
    bi <- (geno[1, ] == 0 & geno[2, ] == 2) | (geno[1, ] == 2 & geno[2, ] == 0)
    
    gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = TRUE)
    expect_equal(getValidSnp(gds), bi, ignore_attr = TRUE)
    gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = TRUE)
    expect_equal(getValidSnp(gds), bi, ignore_attr = TRUE)
})

test_that("Set parents with checking marker filtering",{
    parents <- c("Founder1", "Founder2")
    gds <- setParents(gds, parents)
    geno <- getGenotype(gds, parents = "only")
    bi <- (geno[1, ] == 0 & geno[2, ] == 2) | (geno[1, ] == 2 & geno[2, ] == 0)
    mono <- geno[1, ] %in% c(0, 2) & geno[2, ] %in% c(0, 2)
    
    gds <- setParents(gds, parents, flip = TRUE, mono = TRUE, bi = TRUE)
    expect_true(hasFlipped(gds))
    flipped <- getFlipped(gds)
    expect_equal(flipped, geno[1, bi & mono] == 0, ignore_attr = TRUE)
    expect_equal(getValidSnp(gds), bi & mono, ignore_attr = TRUE)
    gds <- setParents(gds, parents, flip = TRUE, mono = TRUE, bi = TRUE)
    expect_true(hasFlipped(gds))
    flipped <- getFlipped(gds)
    expect_equal(flipped, geno[1, bi & mono] == 0, ignore_attr = TRUE)
    expect_equal(getValidSnp(gds), bi & mono, ignore_attr = TRUE)
})
