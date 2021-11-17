library(GBScleanR)

vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn, gds_fn, TRUE, FALSE)
gds <- loadGDS(gds_fn)
gds <- countGenotype(gds)
gds <- countRead(gds)
gds <- calcReadStats(gds)
on.exit({unlink(gds_fn)})

test_that("Thin out markers", {
    gds <- countGenotype(gds)
    range <- 10000
    
    expect_error(thinMarker(gds, range = -1000),
                 "range should be")
    expect_error(thinMarker(gds, range = "1000"),
                 "range should be")
    
    gds <- thinMarker(gds, range = range)
    pos <- getPosition(gds, FALSE)
    chr <- getChromosome(gds, FALSE)
    missing <- getCountGenoMissing(gds, valid = FALSE)
    valid <- tapply(seq_along(pos), chr, function(x){
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
    expect_equal(getValidSnp(gds), unlist(valid), ignore_attr = TRUE)
})

test_that("Sample filtering", {
    expect_error(setScanFilter(gds, id = 1:100))
    omit_id <- sample(getScanID(gds), 10)
    gds <- setScanFilter(gds, id = omit_id)
    expect_false(any(getScanID(gds) %in% omit_id))
    
    gds <- resetScanFilters(gds)
    missing <- getCountGenoMissing(gds, "scan", FALSE, TRUE)
    expect_error(setScanFilter(gds, missing = 10))
    expect_error(setScanFilter(gds, missing = "0.2"))
    gds <- setScanFilter(gds, missing = 0.2)
    expect_equal(getValidScan(gds), missing <= 0.2)
    
    gds <- resetScanFilters(gds)
    het <- getCountGenoHet(gds, "scan", FALSE, TRUE)
    expect_error(setScanFilter(gds, het = 10))
    expect_error(setScanFilter(gds, het = "0.2"))
    gds <- setScanFilter(gds, het = c(0.25, 0.75))
    expect_equal(getValidScan(gds), het >= 0.25 & het <= 0.75)
    
    gds <- resetScanFilters(gds)
    mac <- getMAC(gds, "scan", FALSE)
    expect_error(setScanFilter(gds, mac = -10))
    expect_error(setScanFilter(gds, mac = "100"))
    gds <- setScanFilter(gds, mac = 100)
    expect_equal(getValidScan(gds), mac >= 100)
    
    gds <- resetScanFilters(gds)
    maf <- getMAF(gds, "scan", FALSE)
    expect_error(setScanFilter(gds, maf = 10))
    expect_error(setScanFilter(gds, maf = "0.2"))
    gds <- setScanFilter(gds, maf = 0.25)
    expect_equal(getValidScan(gds), maf >= 0.25)
    
    gds <- resetScanFilters(gds)
    ad_ref <- getCountReadRef(gds, "scan", FALSE)
    expect_error(setScanFilter(gds, ad_ref = -1000))
    expect_error(setScanFilter(gds, ad_ref = "1000"))
    gds <- setScanFilter(gds, ad_ref = c(1000, 2000))
    expect_equal(getValidScan(gds), ad_ref >= 1000 & ad_ref <= 2000)
    
    gds <- resetScanFilters(gds)
    ad_alt <- getCountReadAlt(gds, "scan", FALSE)
    expect_error(setScanFilter(gds, ad_alt = -100))
    expect_error(setScanFilter(gds, ad_alt = "1000"))
    gds <- setScanFilter(gds, ad_alt = c(1000, 1500))
    expect_equal(getValidScan(gds), ad_alt >= 1000 & ad_alt <= 1500)
    
    gds <- resetScanFilters(gds)
    dp <- getCountRead(gds, "scan", FALSE)
    expect_error(setScanFilter(gds, dp = -100))
    expect_error(setScanFilter(gds, dp = "1000"))
    gds <- setScanFilter(gds, dp = c(1000, 2000))
    expect_equal(getValidScan(gds), dp >= 1000 & dp <= 2000)
    
    gds <- resetScanFilters(gds)
    mean_ref <- getMeanReadRef(gds, "scan", FALSE)
    expect_error(setScanFilter(gds, mean_ref = -100))
    expect_error(setScanFilter(gds, mean_ref = "1000"))
    gds <- setScanFilter(gds, mean_ref = c(1000, 2000))
    expect_equal(getValidScan(gds), mean_ref >= 1000 & mean_ref <= 2000)
    
    gds <- resetScanFilters(gds)
    mean_alt <- getMeanReadAlt(gds, "scan", FALSE)
    expect_error(setScanFilter(gds, mean_alt = -100))
    expect_error(setScanFilter(gds, mean_alt = "1000"))
    gds <- setScanFilter(gds, mean_alt = c(1000, 1500))
    expect_equal(getValidScan(gds), mean_alt >= 1000 & mean_alt <= 1500)
    
    gds <- resetScanFilters(gds)
    sd_ref <- getSDReadRef(gds, "scan", FALSE)
    expect_error(setScanFilter(gds, sd_ref = -100))
    expect_error(setScanFilter(gds, sd_ref = "1000"))
    gds <- setScanFilter(gds, sd_ref = 2000)
    expect_equal(getValidScan(gds), sd_ref <= 2000)
    
    gds <- resetScanFilters(gds)
    sd_alt <- getSDReadAlt(gds, "scan", FALSE)
    expect_error(setScanFilter(gds, sd_alt = -100))
    expect_error(setScanFilter(gds, sd_alt = "1000"))
    gds <- setScanFilter(gds, sd_alt = 2000)
    expect_equal(getValidScan(gds), sd_alt <= 2000)
    
    gds <- resetScanFilters(gds)
    valid <- !(getScanID(gds) %in% omit_id) & missing <= 0.2 &
        het >= 0.25 & het <= 0.75 & maf >= 0.25
    gds <- setScanFilter(gds, id = omit_id, missing = 0.2,
                         het = c(0.25, 0.75), maf = 0.25)
    expect_equal(getValidScan(gds), valid)
})


test_that("Marker filtering", {
    expect_error(setSnpFilter(gds, id = "100"))
    omit_id <- sample(getSnpID(gds), 10)
    gds <- setSnpFilter(gds, id = omit_id)
    expect_false(any(getSnpID(gds) %in% omit_id))
    
    gds <- resetSnpFilters(gds)
    missing <- getCountGenoMissing(gds, "snp", FALSE, TRUE)
    expect_error(setSnpFilter(gds, missing = 10))
    expect_error(setSnpFilter(gds, missing = "0.2"))
    gds <- setSnpFilter(gds, missing = 0.2)
    expect_equal(getValidSnp(gds), missing <= 0.2)
    
    gds <- resetSnpFilters(gds)
    het <- getCountGenoHet(gds, "snp", FALSE, TRUE)
    expect_error(setSnpFilter(gds, het = 10))
    expect_error(setSnpFilter(gds, het = "0.2"))
    gds <- setSnpFilter(gds, het = c(0.25, 0.75))
    expect_equal(getValidSnp(gds), het >= 0.25 & het <= 0.75)
    
    gds <- resetSnpFilters(gds)
    mac <- getMAC(gds, "snp", FALSE)
    expect_error(setSnpFilter(gds, mac = -10))
    expect_error(setSnpFilter(gds, mac = "100"))
    gds <- setSnpFilter(gds, mac = 100)
    expect_equal(getValidSnp(gds), mac >= 100)
    
    gds <- resetSnpFilters(gds)
    maf <- getMAF(gds, "snp", FALSE)
    expect_error(setSnpFilter(gds, maf = 10))
    expect_error(setSnpFilter(gds, maf = "0.2"))
    gds <- setSnpFilter(gds, maf = 0.25)
    expect_equal(getValidSnp(gds), maf >= 0.25)
    
    gds <- resetSnpFilters(gds)
    ad_ref <- getCountReadRef(gds, "snp", FALSE)
    expect_error(setSnpFilter(gds, ad_ref = -1000))
    expect_error(setSnpFilter(gds, ad_ref = "1000"))
    gds <- setSnpFilter(gds, ad_ref = c(1000, 2000))
    expect_equal(getValidSnp(gds), ad_ref >= 1000 & ad_ref <= 2000)
    
    gds <- resetSnpFilters(gds)
    ad_alt <- getCountReadAlt(gds, "snp", FALSE)
    expect_error(setSnpFilter(gds, ad_alt = -100))
    expect_error(setSnpFilter(gds, ad_alt = "1000"))
    gds <- setSnpFilter(gds, ad_alt = c(1000, 1500))
    expect_equal(getValidSnp(gds), ad_alt >= 1000 & ad_alt <= 1500)
    
    gds <- resetSnpFilters(gds)
    dp <- getCountRead(gds, "snp", FALSE)
    expect_error(setSnpFilter(gds, dp = -100))
    expect_error(setSnpFilter(gds, dp = "1000"))
    gds <- setSnpFilter(gds, dp = c(1000, 2000))
    expect_equal(getValidSnp(gds), dp >= 1000 & dp <= 2000)
    
    gds <- resetSnpFilters(gds)
    mean_ref <- getMeanReadRef(gds, "snp", FALSE)
    expect_error(setSnpFilter(gds, mean_ref = -100))
    expect_error(setSnpFilter(gds, mean_ref = "1000"))
    gds <- setSnpFilter(gds, mean_ref = c(1000, 2000))
    expect_equal(getValidSnp(gds), mean_ref >= 1000 & mean_ref <= 2000)
    
    gds <- resetSnpFilters(gds)
    mean_alt <- getMeanReadAlt(gds, "snp", FALSE)
    expect_error(setSnpFilter(gds, mean_alt = -100))
    expect_error(setSnpFilter(gds, mean_alt = "1000"))
    gds <- setSnpFilter(gds, mean_alt = c(1000, 1500))
    expect_equal(getValidSnp(gds), mean_alt >= 1000 & mean_alt <= 1500)
    
    gds <- resetSnpFilters(gds)
    sd_ref <- getSDReadRef(gds, "snp", FALSE)
    expect_error(setSnpFilter(gds, sd_ref = -100))
    expect_error(setSnpFilter(gds, sd_ref = "1000"))
    gds <- setSnpFilter(gds, sd_ref = 2000)
    expect_equal(getValidSnp(gds), sd_ref <= 2000)
    
    gds <- resetSnpFilters(gds)
    sd_alt <- getSDReadAlt(gds, "snp", FALSE)
    expect_error(setSnpFilter(gds, sd_alt = -100))
    expect_error(setSnpFilter(gds, sd_alt = "1000"))
    gds <- setSnpFilter(gds, sd_alt = 2000)
    expect_equal(getValidSnp(gds), sd_alt <= 2000)
    
    gds <- resetSnpFilters(gds)
    valid <- !(getSnpID(gds) %in% omit_id) & missing <= 0.2 &
        het >= 0.25 & het <= 0.75 & maf >= 0.25
    gds <- setSnpFilter(gds, id = omit_id, missing = 0.2,
                         het = c(0.25, 0.75), maf = 0.25)
    expect_equal(getValidSnp(gds), valid)
})

    

test_that("Genotype call filtering", {
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
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
    
    gds <- setRawGenotype(gds)
    gds <- setCallFilter(gds, ref_count = c(5, 10))
    geno <- getGenotype(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    valid <- ref >= 5 & ref <= 10
    geno[!valid] <- NA
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
    
    gds <- setRawGenotype(gds)
    gds <- setCallFilter(gds, alt_count = c(5, 10))
    geno <- getGenotype(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    valid <- alt >= 5 & alt <= 10
    geno[!valid] <- NA
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
    
    gds <- setRawGenotype(gds)
    gds <- setCallFilter(gds, norm_dp_count = c(1000, 5000))
    geno <- getGenotype(gds)
    read <- getRead(gds, node = "norm")
    ref <- read$ref
    alt <- read$alt
    dp <- ref + alt
    valid <- dp >= 1000 & dp <= 5000
    geno[!valid] <- NA
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
    
    gds <- setRawGenotype(gds)
    gds <- setCallFilter(gds, norm_ref_count = c(1000, 5000))
    geno <- getGenotype(gds)
    read <- getRead(gds, node = "norm")
    ref <- read$ref
    alt <- read$alt
    valid <- ref >= 1000 & ref <= 5000
    geno[!valid] <- NA
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
    
    gds <- setRawGenotype(gds)
    gds <- setCallFilter(gds, norm_alt_count = c(1000, 5000))
    geno <- getGenotype(gds)
    read <- getRead(gds, node = "norm")
    ref <- read$ref
    alt <- read$alt
    valid <- alt >= 1000 & alt <= 5000
    geno[!valid] <- NA
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
    
    gds <- setRawGenotype(gds)
    gds <- setCallFilter(gds, scan_dp_qtile = c(0.1, 0.9))
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
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
    
    gds <- setRawGenotype(gds)
    gds <- setCallFilter(gds, scan_ref_qtile = c(0.1, 0.9))
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
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
    
    gds <- setRawGenotype(gds)
    gds <- setCallFilter(gds, scan_alt_qtile = c(0.1, 0.9))
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
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
    
    gds <- setRawGenotype(gds)
    gds <- setCallFilter(gds, snp_dp_qtile = c(0.1, 0.9))
    geno <- getGenotype(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    dp <- ref + alt
    valid <- sapply(1:ncol(dp), function(i){
        thr <- quantile(dp[, i], c(0.1, 0.9))
        return(dp[, i] >= thr[1] & dp[, i] <= thr[2])
    })
    geno[!valid] <- NA
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
    
    gds <- setRawGenotype(gds)
    gds <- setCallFilter(gds, snp_ref_qtile = c(0.1, 0.9))
    geno <- getGenotype(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    valid <- sapply(1:ncol(ref), function(i){
        thr <- quantile(ref[, i], c(0.1, 0.9))
        return(ref[, i] >= thr[1] & ref[, i] <= thr[2])
    })
    geno[!valid] <- NA
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
    
    gds <- setRawGenotype(gds)
    gds <- setCallFilter(gds, snp_alt_qtile = c(0.1, 0.9))
    geno <- getGenotype(gds)
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    valid <- sapply(1:ncol(alt), function(i){
        thr <- quantile(alt[, i], c(0.1, 0.9))
        return(alt[, i] >= thr[1] & alt[, i] <= thr[2])
    })
    geno[!valid] <- NA
    read <- getRead(gds)
    ref <- read$ref
    alt <- read$alt
    ref[!valid] <- 0
    alt[!valid] <- 0
    expect_equal(getGenotype(gds, node = "filt"), geno)
    fread <- getRead(gds, node = "filt")
    expect_equal(fread$ref, ref)
    expect_equal(fread$alt, alt)
})
