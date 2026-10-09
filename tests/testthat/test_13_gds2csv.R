vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, verbose = FALSE)
on.exit({unlink(gds_fn)})

.read_gbsr_csv <- function(path) {
    utils::read.csv(path, header = FALSE, check.names = FALSE,
                    na.strings = c("NA", ""))
}

test_that("gbsrGDS2CSV exports raw genotypes", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    csv_fn <- tempfile("out", fileext = ".csv")
    on.exit({
        unlink(csv_fn)
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
    })

    out <- suppressWarnings(gbsrGDS2CSV(gds, csv_fn, node = "raw"))
    expect_equal(out, csv_fn)
    expect_true(file.exists(csv_fn))

    tab <- .read_gbsr_csv(csv_fn)
    expect_equal(nrow(tab), nsam(gds) + 2) # Chr, Pos, then samples
    expect_equal(ncol(tab), nmar(gds) + 1) # row label + markers
    expect_equal(as.character(tab[1, 1]), "Chr")
    expect_equal(as.character(tab[2, 1]), "Pos")
    expect_equal(as.character(unlist(tab[1, -1])),
                 as.character(getChromosome(gds)),
                 ignore_attr = TRUE)
    expect_equal(as.numeric(unlist(tab[2, -1])),
                 as.numeric(getPosition(gds)),
                 ignore_attr = TRUE)

    geno <- getGenotype(gds)
    exported <- as.matrix(tab[-(1:2), -1])
    storage.mode(exported) <- "numeric"
    expect_equal(exported, geno, ignore_attr = TRUE)
})

test_that("gbsrGDS2CSV respects sample and marker filters", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    on.exit({
        if (isOpenGDS(gds)) closeGDS(gds, verbose = FALSE)
    })
    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), TRUE)
    valid_snp[1] <- TRUE
    valid_scan[1] <- TRUE
    validMar(gds) <- valid_snp
    validSam(gds) <- valid_scan

    csv_fn <- tempfile("filt", fileext = ".csv")
    on.exit(unlink(csv_fn), add = TRUE)
    suppressWarnings(gbsrGDS2CSV(gds, csv_fn, node = "raw"))
    tab <- .read_gbsr_csv(csv_fn)
    expect_equal(nrow(tab), nsam(gds) + 2)
    expect_equal(ncol(tab), nmar(gds) + 1)
})
