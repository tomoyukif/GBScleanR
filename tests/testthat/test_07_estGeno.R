.make_temp_gds <- function() {
    vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
    gds_fn <- tempfile("sample", fileext = ".gds")
    gbsrVCF2GDS(vcf_fn, gds_fn, force = TRUE, verbose = FALSE)
    gds_fn
}

.test_estGeno_outputs <- function(gds) {
    p_hap <- getHaplotype(gds, parents = "only")
    expect_true(all(p_hap[, 1, ] == 1))
    expect_true(all(p_hap[, 2, ] == 2))

    hap <- getHaplotype(gds)
    expect_equal(dim(hap), c(2, nsam(gds), nmar(gds)))
    # Haplotype codes are founder member IDs (here 1/2); missing may be NA.
    expect_true(all(is.na(hap) | hap %in% c(0, 1, 2)))

    p_geno <- getGenotype(gds, node = "parents")
    expect_equal(nrow(p_geno), sum(getParents(gds, TRUE)) * 2)
    expect_equal(ncol(p_geno), nmar(gds, FALSE))

    cor_geno <- getGenotype(gds, "cor")
    expect_equal(dim(cor_geno), c(nsam(gds), nmar(gds)))
    expect_true(all(cor_geno %in% c(0, 1, 2, NA)))

    # Non-missing corrected calls should cover the three genotype classes
    expect_gt(sum(cor_geno == 0, na.rm = TRUE), 0)
    expect_gt(sum(cor_geno == 1, na.rm = TRUE), 0)
    expect_gt(sum(cor_geno == 2, na.rm = TRUE), 0)
}

test_that("estGeno with initScheme/addScheme", {
    gds_fn <- .make_temp_gds()
    gds <- loadGDS(gds_fn, verbose = FALSE)
    on.exit({
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
        unlink(gds_fn)
    })
    gds <- setParents(gds, grep("Founder", getSamID(gds), value = TRUE))
    gds <- initScheme(gds, cbind(c(1:2)))
    gds <- addScheme(gds, "self")
    gds <- estGeno(gds)
    .test_estGeno_outputs(gds)
})

test_that("estGeno with makeScheme", {
    gds_fn <- .make_temp_gds()
    gds <- loadGDS(gds_fn, verbose = FALSE)
    on.exit({
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
        unlink(gds_fn)
    })
    gds <- setParents(gds, grep("Founder", getSamID(gds), value = TRUE))
    gds <- makeScheme(gds, generation = 2, crosstype = "selfing")
    gds <- estGeno(gds)
    .test_estGeno_outputs(gds)

    dosage <- getGenotype(gds, node = "dosage")
    expect_equal(dim(dosage), c(nsam(gds), nmar(gds)))
    expect_true(all(dosage %in% c(0, 1, 2, NA)))
})

test_that("estGeno after call filtering uses filt node", {
    gds_fn <- .make_temp_gds()
    gds <- loadGDS(gds_fn, verbose = FALSE)
    on.exit({
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
        unlink(gds_fn)
    })
    gds <- setParents(gds, grep("Founder", getSamID(gds), value = TRUE),
                      mono = TRUE, bi = TRUE)
    gds <- setCallFilter(gds, dp_count = c(2, Inf))
    gds <- makeScheme(gds, generation = 2, crosstype = "selfing")
    gds <- estGeno(gds, node = "filt")
    cor_geno <- getGenotype(gds, "cor")
    expect_equal(dim(cor_geno), c(nsam(gds), nmar(gds)))
    expect_true(all(cor_geno %in% c(0, 1, 2, NA)))
})

test_that("corrected genotypes export to CSV and plotting works", {
    gds_fn <- .make_temp_gds()
    gds <- loadGDS(gds_fn, verbose = FALSE)
    on.exit({
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
        unlink(gds_fn)
    })
    gds <- setParents(gds, grep("Founder", getSamID(gds), value = TRUE))
    gds <- makeScheme(gds, generation = 2, crosstype = "selfing")
    gds <- estGeno(gds)

    csv_out <- tempfile("cor", fileext = ".csv")
    vcf_out <- tempfile("cor", fileext = ".vcf")
    on.exit(unlink(c(csv_out, vcf_out)), add = TRUE)

    suppressWarnings(gbsrGDS2CSV(gds, csv_out, node = "cor",
                                 incl_parents = FALSE))
    expect_true(file.exists(csv_out))
    tab <- utils::read.csv(csv_out, header = FALSE, check.names = FALSE)
    expect_equal(ncol(tab), nmar(gds) + 1)
    expect_equal(nrow(tab), nsam(gds) + 2)

    suppressWarnings(gbsrGDS2VCF(gds, vcf_out, node = "cor"))
    expect_true(file.exists(vcf_out))
    expect_gt(file.info(vcf_out)$size, 0)

    expect_s3_class(plotDosage(gds, ind = 1), "ggplot")
    expect_s3_class(plotReadRatio(gds, ind = 1), "ggplot")
})

test_that("estGeno with replicates yields identical genotypes", {
    gds_fn <- .make_temp_gds()
    gds <- loadGDS(gds_fn, verbose = FALSE)
    on.exit({
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
        unlink(gds_fn)
    })
    parents <- grep("Founder", getSamID(gds, valid = FALSE), value = TRUE)
    gds <- setParents(gds, parents)

    n <- nsam(gds, valid = FALSE)
    reps <- seq_len(n)
    offspring <- which(!getParents(gds, TRUE))
    reps[offspring[1:2]] <- reps[offspring[1]]
    gds <- setReplicates(gds, replicates = reps)
    gds <- makeScheme(gds, generation = 2, crosstype = "selfing")
    gds <- estGeno(gds)

    cor_geno <- getGenotype(gds, "cor")
    id_valid <- getSamID(gds)
    all_id <- getSamID(gds, valid = FALSE)
    idx <- match(all_id[offspring[1:2]], id_valid)
    expect_false(any(is.na(idx)))
    expect_equal(cor_geno[idx[1], ], cor_geno[idx[2], ], ignore_attr = TRUE)
})
