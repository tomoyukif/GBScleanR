vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, verbose = FALSE)
on.exit({unlink(gds_fn)})

test_that("setReplicates and getReplicates", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    n <- nsam(gds, valid = FALSE)
    default_rep <- getReplicates(gds)
    expect_equal(as.numeric(default_rep), seq_len(n), ignore_attr = TRUE)

    expect_error(setReplicates(gds, replicates = 1:3),
                 "length of the vector")

    reps <- rep(seq_len(ceiling(n / 2)), each = 2)[seq_len(n)]
    gds <- setReplicates(gds, replicates = reps)
    got <- getReplicates(gds)
    expect_equal(length(unique(got)), length(unique(reps)))
    expect_equal(as.numeric(factor(got)), as.numeric(factor(reps)),
                 ignore_attr = TRUE)

    # Reset to unique IDs
    gds <- setReplicates(gds, replicates = seq_len(n))
    expect_equal(as.numeric(getReplicates(gds)), seq_len(n),
                 ignore_attr = TRUE)
    closeGDS(gds, verbose = FALSE)
})

test_that("setFixedParameter and getFixedParameter", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    expect_null(getFixedParameter(gds))

    n_valid <- nmar(gds)
    bias <- runif(n_valid, 0.2, 0.8)
    mismap <- cbind(runif(n_valid, 0, 0.05), runif(n_valid, 0, 0.05))
    gds <- setFixedParameter(gds, bias = bias, mismap = mismap)
    fp <- getFixedParameter(gds)
    expect_true(all(c("bias", "mismap_ref", "mismap_alt") %in% names(fp)))
    expect_equal(fp$bias, bias, ignore_attr = TRUE)
    expect_equal(fp$mismap_ref, mismap[, 1], ignore_attr = TRUE)
    expect_equal(fp$mismap_alt, mismap[, 2], ignore_attr = TRUE)

    gds <- setFixedParameter(gds, bias = NULL, mismap = NULL)
    expect_null(getFixedParameter(gds))
    closeGDS(gds, verbose = FALSE)
})

test_that("resetFilter clears sample, marker, and call filters", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    gds <- countRead(gds)
    gds <- countGenotype(gds)

    omit_sam <- sample(getSamID(gds), 5)
    omit_mar <- sample(getMarID(gds), 5)
    gds <- setSamFilter(gds, id = omit_sam)
    gds <- setMarFilter(gds, id = omit_mar)
    gds <- setCallFilter(gds, dp_count = c(5, Inf))
    expect_lt(nsam(gds), nsam(gds, valid = FALSE))
    expect_lt(nmar(gds), nmar(gds, valid = FALSE))
    expect_true(gdsfmt::exist.gdsn(gds, "annotation/format/FAD") ||
                    gdsfmt::exist.gdsn(gds, "annotation/format/FGT"))

    gds <- resetFilter(gds)
    expect_true(all(validSam(gds)))
    expect_true(all(validMar(gds)))
    closeGDS(gds, verbose = FALSE)
})
