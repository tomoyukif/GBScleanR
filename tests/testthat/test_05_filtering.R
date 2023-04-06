library(GBScleanR)

vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, verbose = FALSE)
on.exit({unlink(gds_fn)})

test_that("Thin out markers", {
    gds <- loadGDS(gds_fn)
    gds <- countGenotype(gds)
    range <- 10000
    expect_error(thinMarker(gds, range = -1000), "range should be")
    expect_error(thinMarker(gds, range = "1000"),"range should be")

    gds <- thinMarker(gds, range = range)
    pos <- getPosition(gds, FALSE)
    chr <- getChromosome(gds, FALSE)
    missing <- getCountGenoMissing(gds, valid = FALSE)
    valid <- tapply(seq_along(pos), chr , function(x){
        x <- seq_along(pos)
        n_snp <- length(x)
        x_pos <- pos[x]
        x_missing <- missing[x]
        valid <- rep(TRUE, n_snp)
        i <- 1
        j <- 2
        while (TRUE){
            mar1 <- x_pos[i]
            mar2 <- x_pos[j]
            if(mar2 - mar1 <= range){
                if(x_missing[i] <= x_missing[j]){
                    valid[j] <- FALSE
                    j <- j + 1
                } else {
                    valid[i] <- FALSE
                    i <- j
                    j <- j + 1
                }
            } else {
                i <- j
                j <- j + 1
            }
            if(j > n_snp){
                break
            }
        }
        return(valid)
    })
    expect_equal(validMar(gds), unlist(valid), ignore_attr = TRUE)
    closeGDS(gds)
})

test_that("Sample filtering", {
    gds <- loadGDS(gds_fn)
    gds <- countRead(gds)
    gds <- countGenotype(gds)

    expect_error(setSamFilter(gds, id = 1:100))
    omit_id <- sample(getSamID(gds), 10)
    gds <- setSamFilter(gds, id = omit_id)
    expect_false(any(getSamID(gds) %in% omit_id))

    gds <- resetSamFilter(gds)
    missing <- getCountGenoMissing(gds, "sample", FALSE, TRUE)
    expect_error(setSamFilter(gds, missing = 10))
    expect_error(setSamFilter(gds, missing = "0.2"))
    gds <- setSamFilter(gds, missing = 0.2)
    expect_equal(validSam(gds), missing <= 0.2)

    gds <- resetSamFilter(gds)
    het <- getCountGenoHet(gds, "sample", FALSE, TRUE)
    expect_error(setSamFilter(gds, het = 10))
    expect_error(setSamFilter(gds, het = "0.2"))
    gds <- setSamFilter(gds, het = c(0.25, 0.75))
    expect_equal(validSam(gds), het >= 0.25 & het <= 0.75)

    gds <- resetSamFilter(gds)
    mac <- getMAC(gds, "sample", FALSE)
    expect_error(setSamFilter(gds, mac = -10))
    expect_error(setSamFilter(gds, mac = "100"))
    gds <- setSamFilter(gds, mac = 100)
    expect_equal(validSam(gds), mac >= 100)

    gds <- resetSamFilter(gds)
    maf <- getMAF(gds, "sample", FALSE)
    expect_error(setSamFilter(gds, maf = 10))
    expect_error(setSamFilter(gds, maf = "0.2"))
    gds <- setSamFilter(gds, maf = 0.25)
    expect_equal(validSam(gds), maf >= 0.25)

    gds <- resetSamFilter(gds)
    ad_ref <- getCountReadRef(gds, "sample", FALSE)
    expect_error(setSamFilter(gds, ad_ref = -1000))
    expect_error(setSamFilter(gds, ad_ref = "1000"))
    gds <- setSamFilter(gds, ad_ref = c(1000, 2000))
    expect_equal(validSam(gds), ad_ref >= 1000 & ad_ref <= 2000)

    gds <- resetSamFilter(gds)
    ad_alt <- getCountReadAlt(gds, "sample", FALSE)
    expect_error(setSamFilter(gds, ad_alt = -100))
    expect_error(setSamFilter(gds, ad_alt = "1000"))
    gds <- setSamFilter(gds, ad_alt = c(1000, 1500))
    expect_equal(validSam(gds), ad_alt >= 1000 & ad_alt <= 1500)

    gds <- resetSamFilter(gds)
    dp <- getCountRead(gds, "sample", FALSE)
    expect_error(setSamFilter(gds, dp = -100))
    expect_error(setSamFilter(gds, dp = "1000"))
    gds <- setSamFilter(gds, dp = c(1000, 2000))
    expect_equal(validSam(gds), dp >= 1000 & dp <= 2000)

    gds <- resetSamFilter(gds)
    mean_ref <- getMeanReadRef(gds, "sample", FALSE)
    expect_error(setSamFilter(gds, mean_ref = -100))
    expect_error(setSamFilter(gds, mean_ref = "1000"))
    gds <- setSamFilter(gds, mean_ref = c(1000, 2000))
    expect_equal(validSam(gds), mean_ref >= 1000 & mean_ref <= 2000)

    gds <- resetSamFilter(gds)
    mean_alt <- getMeanReadAlt(gds, "sample", FALSE)
    expect_error(setSamFilter(gds, mean_alt = -100))
    expect_error(setSamFilter(gds, mean_alt = "1000"))
    gds <- setSamFilter(gds, mean_alt = c(1000, 1500))
    expect_equal(validSam(gds), mean_alt >= 1000 & mean_alt <= 1500)

    gds <- resetSamFilter(gds)
    sd_ref <- getSDReadRef(gds, "sample", FALSE)
    expect_error(setSamFilter(gds, sd_ref = -100))
    expect_error(setSamFilter(gds, sd_ref = "1000"))
    gds <- setSamFilter(gds, sd_ref = 2000)
    expect_equal(validSam(gds), sd_ref <= 2000)

    gds <- resetSamFilter(gds)
    sd_alt <- getSDReadAlt(gds, "sample", FALSE)
    expect_error(setSamFilter(gds, sd_alt = -100))
    expect_error(setSamFilter(gds, sd_alt = "1000"))
    gds <- setSamFilter(gds, sd_alt = 2000)
    expect_equal(validSam(gds), sd_alt <= 2000)

    gds <- resetSamFilter(gds)
    valid <- !(getSamID(gds) %in% omit_id) & missing <= 0.2 &
        het >= 0.25 & het <= 0.75 & maf >= 0.25
    gds <- setSamFilter(gds, id = omit_id, missing = 0.2,
                         het = c(0.25, 0.75), maf = 0.25)
    expect_equal(validSam(gds), valid)
    closeGDS(gds)
})


test_that("Marker filtering", {
    gds <- loadGDS(gds_fn)
    gds <- countRead(gds)
    gds <- countGenotype(gds)

    expect_error(setMarFilter(gds, id = "100"))
    omit_id <- sample(getMarID(gds), 10)
    gds <- setMarFilter(gds, id = omit_id)
    expect_false(any(getMarID(gds) %in% omit_id))

    gds <- resetMarFilter(gds)
    missing <- getCountGenoMissing(gds, "marker", FALSE, TRUE)
    expect_error(setMarFilter(gds, missing = 10))
    expect_error(setMarFilter(gds, missing = "0.2"))
    gds <- setMarFilter(gds, missing = 0.2)
    expect_equal(validMar(gds), missing <= 0.2)

    gds <- resetMarFilter(gds)
    het <- getCountGenoHet(gds, "marker", FALSE, TRUE)
    expect_error(setMarFilter(gds, het = 10))
    expect_error(setMarFilter(gds, het = "0.2"))
    gds <- setMarFilter(gds, het = c(0.25, 0.75))
    expect_equal(validMar(gds), het >= 0.25 & het <= 0.75)

    gds <- resetMarFilter(gds)
    mac <- getMAC(gds, "marker", FALSE)
    expect_error(setMarFilter(gds, mac = -10))
    expect_error(setMarFilter(gds, mac = "100"))
    gds <- setMarFilter(gds, mac = 100)
    expect_equal(validMar(gds), mac >= 100)

    gds <- resetMarFilter(gds)
    maf <- getMAF(gds, "marker", FALSE)
    expect_error(setMarFilter(gds, maf = 10))
    expect_error(setMarFilter(gds, maf = "0.2"))
    gds <- setMarFilter(gds, maf = 0.25)
    expect_equal(validMar(gds), maf >= 0.25)

    gds <- resetMarFilter(gds)
    ad_ref <- getCountReadRef(gds, "marker", FALSE)
    expect_error(setMarFilter(gds, ad_ref = -1000))
    expect_error(setMarFilter(gds, ad_ref = "1000"))
    gds <- setMarFilter(gds, ad_ref = c(1000, 2000))
    expect_equal(validMar(gds), ad_ref >= 1000 & ad_ref <= 2000)

    gds <- resetMarFilter(gds)
    ad_alt <- getCountReadAlt(gds, "marker", FALSE)
    expect_error(setMarFilter(gds, ad_alt = -100))
    expect_error(setMarFilter(gds, ad_alt = "1000"))
    gds <- setMarFilter(gds, ad_alt = c(1000, 1500))
    expect_equal(validMar(gds), ad_alt >= 1000 & ad_alt <= 1500)

    gds <- resetMarFilter(gds)
    dp <- getCountRead(gds, "marker", FALSE)
    expect_error(setMarFilter(gds, dp = -100))
    expect_error(setMarFilter(gds, dp = "1000"))
    gds <- setMarFilter(gds, dp = c(1000, 2000))
    expect_equal(validMar(gds), dp >= 1000 & dp <= 2000)

    gds <- resetMarFilter(gds)
    mean_ref <- getMeanReadRef(gds, "marker", FALSE)
    expect_error(setMarFilter(gds, mean_ref = -100))
    expect_error(setMarFilter(gds, mean_ref = "1000"))
    gds <- setMarFilter(gds, mean_ref = c(1000, 2000))
    expect_equal(validMar(gds), mean_ref >= 1000 & mean_ref <= 2000)

    gds <- resetMarFilter(gds)
    mean_alt <- getMeanReadAlt(gds, "marker", FALSE)
    expect_error(setMarFilter(gds, mean_alt = -100))
    expect_error(setMarFilter(gds, mean_alt = "1000"))
    gds <- setMarFilter(gds, mean_alt = c(1000, 1500))
    expect_equal(validMar(gds), mean_alt >= 1000 & mean_alt <= 1500)

    gds <- resetMarFilter(gds)
    sd_ref <- getSDReadRef(gds, "marker", FALSE)
    expect_error(setMarFilter(gds, sd_ref = -100))
    expect_error(setMarFilter(gds, sd_ref = "1000"))
    gds <- setMarFilter(gds, sd_ref = 2000)
    expect_equal(validMar(gds), sd_ref <= 2000)

    gds <- resetMarFilter(gds)
    sd_alt <- getSDReadAlt(gds, "marker", FALSE)
    expect_error(setMarFilter(gds, sd_alt = -100))
    expect_error(setMarFilter(gds, sd_alt = "1000"))
    gds <- setMarFilter(gds, sd_alt = 2000)
    expect_equal(validMar(gds), sd_alt <= 2000)

    gds <- resetMarFilter(gds)
    valid <- !(getMarID(gds) %in% omit_id) & missing <= 0.2 &
        het >= 0.25 & het <= 0.75 & maf >= 0.25
    gds <- setMarFilter(gds, id = omit_id, missing = 0.2,
                         het = c(0.25, 0.75), maf = 0.25)
    expect_equal(validMar(gds), valid)
    closeGDS(gds)
})

test_that("Genotype call filtering", {
    gds <- loadGDS(gds_fn)
    gds <- countRead(gds)
    gds <- countGenotype(gds)

    gds <- setCallFilter(gds, dp_count = c(5, 10))
    geno <- getGenotype(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    dp <- ref + alt
    valid <- dp >= 5 & dp <= 10
    geno[!valid] <- NA
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno, ignore_attr = TRUE)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref, ignore_attr = TRUE)
    expect_equal(fread$alt, alt, ignore_attr = TRUE)

    gds <- resetCallFilter(gds)
    gds <- setCallFilter(gds, ref_count = c(5, 10))
    geno <- getGenotype(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    valid <- ref >= 5 & ref <= 10
    geno[!valid] <- NA
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno, ignore_attr = TRUE)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref, ignore_attr = TRUE)
    expect_equal(fread$alt, alt, ignore_attr = TRUE)

    gds <- resetCallFilter(gds)
    gds <- setCallFilter(gds, alt_count = c(5, 10))
    geno <- getGenotype(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    valid <- alt >= 5 & alt <= 10
    geno[!valid] <- NA
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno, ignore_attr = TRUE)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref, ignore_attr = TRUE)
    expect_equal(fread$alt, alt, ignore_attr = TRUE)

    gds <- resetCallFilter(gds)
    gds <- setCallFilter(gds, dp_qtile = c(0.1, 0.9))
    geno <- getGenotype(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    dp <- ref + alt
    valid <- sapply(1:nrow(dp), function(i){
        thr <- quantile(dp[i, ], c(0.1, 0.9))
        return(dp[i, ] >= thr[1] & dp[i, ] <= thr[2])
    })
    valid <- t(valid)
    geno[!valid] <- NA
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno, ignore_attr = TRUE)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref, ignore_attr = TRUE)
    expect_equal(fread$alt, alt, ignore_attr = TRUE)

    gds <- resetCallFilter(gds)
    gds <- setCallFilter(gds, ref_qtile = c(0.1, 0.9))
    geno <- getGenotype(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    valid <- sapply(1:nrow(ref), function(i){
        thr <- quantile(ref[i, ], c(0.1, 0.9))
        return(ref[i, ] >= thr[1] & ref[i, ] <= thr[2])
    })
    valid <- t(valid)
    geno[!valid] <- NA
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno, ignore_attr = TRUE)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref, ignore_attr = TRUE)
    expect_equal(fread$alt, alt, ignore_attr = TRUE)

    gds <- resetCallFilter(gds)
    gds <- setCallFilter(gds, alt_qtile = c(0.1, 0.9))
    geno <- getGenotype(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    valid <- sapply(1:nrow(alt), function(i){
        thr <- quantile(alt[i, ], c(0.1, 0.9))
        return(alt[i, ] >= thr[1] & alt[i, ] <= thr[2])
    })
    valid <- t(valid)
    geno[!valid] <- NA
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno, ignore_attr = TRUE)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref, ignore_attr = TRUE)
    expect_equal(fread$alt, alt, ignore_attr = TRUE)
    closeGDS(gds)
})


test_that("Save and load filter info", {
    gds <- loadGDS(gds_fn)
    gds <- countRead(gds)
    gds <- countGenotype(gds)

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), replace = TRUE)
    valid_scan <- sample(c(TRUE, FALSE), nsam(gds), replace = TRUE)
    validMar(gds) <- valid_snp
    validSam(gds) <- valid_scan
    closeGDS(gds, save_filter = TRUE)
    gds <- loadGDS(gds, load_filter = TRUE)
    expect_equal(validMar(gds), valid_snp, ignore_attr = TRUE)
    expect_equal(validSam(gds), valid_scan, ignore_attr = TRUE)
    closeGDS(gds)
})

