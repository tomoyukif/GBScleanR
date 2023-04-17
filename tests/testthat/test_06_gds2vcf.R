library(GBScleanR)

test_that("gbsrGDS2VCF", {
    if(Sys.info()["sysname"] != "Windows"){
        vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
        gds_fn <- tempfile("sample", fileext = ".gds")
        gbsrVCF2GDS(vcf_fn, gds_fn)
        on.exit({unlink(gds_fn)})
        gds <- loadGDS(gds_fn)
        
        out_fn <- tempfile("out", fileext = ".vcf")
        
        valid_snp <- sample(c(TRUE, FALSE), nmar(gds), replace = TRUE)
        valid_scan <- sample(c(TRUE, FALSE), nsam(gds), replace = TRUE)
        validMar(gds) <- valid_snp
        validSam(gds) <- valid_scan
        out_fn <- tempfile("out", fileext = ".vcf")
        gbsrGDS2VCF(gds, out_fn)
        new_gds <- tempfile("newgds", fileext = ".gds")
        gbsrVCF2GDS(out_fn, new_gds)
        newgds <- loadGDS(new_gds)
        expect_equal(nmar(newgds), nmar(gds), ignore_attr = TRUE)
        expect_equal(nsam(newgds), nsam(gds), ignore_attr = TRUE)
        expect_equal(getRead(newgds), getRead(gds), ignore_attr = TRUE)
        expect_equal(getGenotype(newgds), getGenotype(gds), ignore_attr = TRUE)
        
        gds <- setCallFilter(gds, dp_qtile = c(0, 0.8))
        out_fn <- tempfile("out", fileext = ".vcf")
        gbsrGDS2VCF(gds, out_fn)
        new_gds <- tempfile("newgds", fileext = ".gds")
        gbsrVCF2GDS(out_fn, new_gds)
        newgds <- loadGDS(new_gds)
        expect_equal(nmar(newgds), nmar(gds), ignore_attr = TRUE)
        expect_equal(nsam(newgds), nsam(gds), ignore_attr = TRUE)
        expect_equal(getRead(newgds, node = "filt"),
                     getRead(gds, node = "filt"),
                     ignore_attr = TRUE)
        expect_equal(getGenotype(newgds, node = "filt"),
                     getGenotype(gds, node = "filt"),
                     ignore_attr = TRUE)
        closeGDS(newgds)
        closeGDS(gds)
    }
})
