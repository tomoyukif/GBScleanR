vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, verbose = FALSE)
on.exit({unlink(gds_fn)})

test_that("makeScheme requires parents", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    expect_error(makeScheme(gds, generation = 2, crosstype = "selfing"),
                 "Parents have not been set")
    closeGDS(gds, verbose = FALSE)
})

test_that("makeScheme builds a biparental selfing pedigree", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    parents <- grep("Founder", getSamID(gds), value = TRUE)
    gds <- setParents(gds, parents)
    gds <- makeScheme(gds, generation = 2, crosstype = "selfing")

    scheme <- gds@scheme
    expect_equal(scheme@parents, 1:2)
    expect_equal(length(scheme@mating), 2)
    expect_equal(scheme@mating[[1]], cbind(c(1, 2)), ignore_attr = TRUE)
    expect_true(all(scheme@crosstype[[2]] == "selfing"))
    expect_equal(unlist(scheme@progenies[[2]]), 4)
    closeGDS(gds, verbose = FALSE)
})

test_that("initScheme and addScheme for F2 selfing", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    parents <- grep("Founder", getSamID(gds), value = TRUE)
    gds <- setParents(gds, parents)
    gds <- initScheme(gds, mating = cbind(c(1, 2)))
    expect_equal(gds@scheme@progenies[[1]], 3)

    gds <- addScheme(gds, crosstype = "selfing")
    expect_equal(length(gds@scheme@mating), 2)
    expect_true(all(gds@scheme@crosstype[[2]] == "selfing"))
    expect_equal(gds@scheme@mating[[2]], cbind(c(3, 3)), ignore_attr = TRUE)

    expect_error(addScheme(gds, crosstype = "pairing",
                           mating = cbind(c(1, 2))),
                 "Invalid member ID")
    closeGDS(gds, verbose = FALSE)
})

test_that("assignScheme validates member IDs and length", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    parents <- grep("Founder", getSamID(gds), value = TRUE)
    gds <- setParents(gds, parents)
    gds <- makeScheme(gds, generation = 2, crosstype = "selfing")

    expect_error(assignScheme(gds, id = 4), "length\\(id\\)")
    expect_error(assignScheme(gds, id = rep(1, nsam(gds))),
                 "Invalid specification")

    gds <- assignScheme(gds, id = rep(4, nsam(gds)))
    expect_equal(gds@scheme@samples, rep(4, nsam(gds)))
    expect_equal(gds@sample$pedigree[validSam(gds)], rep(4, nsam(gds)))
    closeGDS(gds, verbose = FALSE)
})

test_that("showScheme runs without error", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    parents <- grep("Founder", getSamID(gds), value = TRUE)
    gds <- setParents(gds, parents)
    gds <- makeScheme(gds, generation = 2, crosstype = "selfing")
    expect_message(showScheme(gds))
    closeGDS(gds, verbose = FALSE)
})

test_that("scheme mating validation rejects inconsistent crosses", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    parents <- grep("Founder", getSamID(gds), value = TRUE)
    gds <- setParents(gds, parents)
    gds <- initScheme(gds, mating = cbind(c(1, 2)))
    expect_error(addScheme(gds, crosstype = "selfing",
                           mating = cbind(c(3, 4))),
                 "Invalid member ID|Different member")
    closeGDS(gds, verbose = FALSE)
})
