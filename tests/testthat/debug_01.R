# library(GBScleanR)
# url <- "https://github.com/tomoyukif/GBSR_SupFiles/raw/main/realdata/gbs_nbolf2.vcf.gz"
# temp_dir <- tempdir()
# vcf_fn <- tempfile(pattern = "nbol_f2", tmpdir = temp_dir, fileext = ".vcf.gz")
# gds_fn <- sub(".vcf.gz", ".gds", vcf_fn)
# download.file(url = url, destfile = vcf_fn)
# on.exit({unlink(temp_dir, force = TRUE)})
# 
# test_that("debug", {
#   gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn,
#               info.import = character(length = 0L),
#               fmt.import = "AD", force = TRUE)
#   gds <- loadGDS(x = gds_fn)
#   setCallFilter(object = gds, ref_qtile = c(0, 0.9), alt_qtile = c(0, 0.9), dp_qtile = c(0, 0.95),
#                 dp_count = c(0, 5), ref_count = c(0, 8), alt_count = c(0, 5))
#   closeGDS(object = gds, save_filter = TRUE)
# })
#   
#   