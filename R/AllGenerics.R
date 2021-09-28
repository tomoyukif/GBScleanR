setGeneric("openGDS", function(object, ...)
  standardGeneric("openGDS"))
setGeneric("isOpenGDS", function(object, ...)
  standardGeneric("isOpenGDS"))
setGeneric("closeGDS", function(object, ...)
  standardGeneric("closeGDS"))
setGeneric("saveSnpAnnot", function(object, ...)
  standardGeneric("saveSnpAnnot"))
setGeneric("saveScanAnnot", function(object, ...)
  standardGeneric("saveScanAnnot"))
setGeneric("loadSnpAnnot", function(object, ...)
  standardGeneric("loadSnpAnnot"))
setGeneric("loadScanAnnot", function(object, ...)
  standardGeneric("loadScanAnnot"))
setGeneric("gbsrGDS2VCF", function(object, out_fn, node = "raw", valid = TRUE,
                                   out_fmt = NULL, out_info = NULL, ...)
  standardGeneric("gbsrGDS2VCF"))

setGeneric("nsnp", function(object, valid = TRUE, ...)
  standardGeneric("nsnp"))
setGeneric("nscan", function(object, valid = TRUE, ...)
  standardGeneric("nscan"))
setGeneric("getValidSnp", function(object, ...)
  standardGeneric("getValidSnp"))
setGeneric("getValidScan", function(object, parents = FALSE, ...)
  standardGeneric("getValidScan"))
setGeneric("setValidSnp", function(object, new, update, ...)
  standardGeneric("setValidSnp"))
setGeneric("setValidScan", function(object, new, update, ...)
  standardGeneric("setValidScan"))
setGeneric("setFlipped", function(object, flipped, ...)
  standardGeneric("setFlipped"))
setGeneric("getFlipped", function(object, valid = TRUE, ...)
  standardGeneric("getFlipped"))
setGeneric("haveFlipped", function(object, ...)
  standardGeneric("haveFlipped"))
setGeneric("flipData", function(object, ...)
  standardGeneric("flipData"))
setGeneric("getGenotype", function(object,
                                   chr = NULL,
                                   node = "raw",
                                   parents = FALSE,
                                   ...)
  standardGeneric("getGenotype"))
setGeneric("getHaplotype", function(object,
                                    chr = NULL,
                                    parents = FALSE,
                                    ...)
  standardGeneric("getHaplotype"))
setGeneric("getRead", function(object,
                               chr = NULL,
                               node = "raw",
                               parents = FALSE,
                               ...)
  standardGeneric("getRead"))
setGeneric("getChromosome", function(object,
                                     valid = TRUE,
                                     levels = FALSE,
                                     name = FALSE,
                                     ...)
  standardGeneric("getChromosome"))
setGeneric("getPosition", function(object, valid = TRUE, ...)
  standardGeneric("getPosition"))
setGeneric("getAlleleA", function(object, valid = TRUE, ...)
  standardGeneric("getAlleleA"))
setGeneric("getAlleleB", function(object, valid = TRUE, ...)
  standardGeneric("getAlleleB"))
setGeneric("getSnpID", function(object, valid = TRUE, ...)
  standardGeneric("getSnpID"))
setGeneric("getScanID", function(object, valid = TRUE, ...)
  standardGeneric("getScanID"))
setGeneric("getPloidy", function(object, valid = TRUE, ...)
  standardGeneric("getPloidy"))
setGeneric("getInfo", function(object, var, valid = TRUE, ...)
  standardGeneric("getInfo"))
setGeneric("getCountReadRef", function(object,
                                       target = "snp",
                                       valid = TRUE,
                                       prop = FALSE,
                                       ...)
  standardGeneric("getCountReadRef"))
setGeneric("getCountReadAlt", function(object,
                                       target = "snp",
                                       valid = TRUE,
                                       prop = FALSE,
                                       ...)
  standardGeneric("getCountReadAlt"))
setGeneric("getCountRead", function(object,
                                    target = "snp",
                                    valid = TRUE,
                                    ...)
  standardGeneric("getCountRead"))
setGeneric("getCountGenoRef", function(object,
                                       target = "snp",
                                       valid = TRUE,
                                       prop = FALSE,
                                       ...)
  standardGeneric("getCountGenoRef"))
setGeneric("getCountGenoHet", function(object,
                                       target = "snp",
                                       valid = TRUE,
                                       prop = FALSE,
                                       ...)
  standardGeneric("getCountGenoHet"))
setGeneric("getCountGenoAlt", function(object,
                                       target = "snp",
                                       valid = TRUE,
                                       prop = FALSE,
                                       ...)
  standardGeneric("getCountGenoAlt"))
setGeneric("getCountGenoMissing", function(object,
                                           target = "snp",
                                           valid = TRUE,
                                           prop = FALSE,
                                           ...)
  standardGeneric("getCountGenoMissing"))
setGeneric("getCountAlleleRef", function(object,
                                         target = "snp",
                                         valid = TRUE,
                                         prop = FALSE,
                                         ...)
  standardGeneric("getCountAlleleRef"))
setGeneric("getCountAlleleAlt", function(object,
                                         target = "snp",
                                         valid = TRUE,
                                         prop = FALSE,
                                         ...)
  standardGeneric("getCountAlleleAlt"))
setGeneric("getCountAlleleMissing", function(object,
                                             target = "snp",
                                             valid = TRUE,
                                             prop = FALSE,
                                             ...)
  standardGeneric("getCountAlleleMissing"))
setGeneric("getMeanReadRef", function(object,
                                      target = "snp",
                                      valid = TRUE,
                                      ...)
  standardGeneric("getMeanReadRef"))
setGeneric("getMeanReadAlt", function(object,
                                      target = "snp",
                                      valid = TRUE,
                                      ...)
  standardGeneric("getMeanReadAlt"))
setGeneric("getSDReadRef", function(object,
                                    target = "snp",
                                    valid = TRUE,
                                    ...)
  standardGeneric("getSDReadRef"))
setGeneric("getSDReadAlt", function(object,
                                    target = "snp",
                                    valid = TRUE,
                                    ...)
  standardGeneric("getSDReadAlt"))
setGeneric("getQtileReadRef", function(object,
                                       target = "snp",
                                       q = 0.5,
                                       valid = TRUE,
                                       ...)
  standardGeneric("getQtileReadRef"))
setGeneric("getQtileReadAlt", function(object,
                                       target = "snp",
                                       q = 0.5,
                                       valid = TRUE,
                                       ...)
  standardGeneric("getQtileReadAlt"))
setGeneric("getMAF", function(object,
                              target = "snp",
                              valid = TRUE,
                              ...)
  standardGeneric("getMAF"))
setGeneric("getMAC", function(object,
                              target = "snp",
                              valid = TRUE,
                              ...)
  standardGeneric("getMAC"))

setGeneric("countGenotype", function(object,
                                     target = "both",
                                     node = "",
                                     ...)
  standardGeneric("countGenotype"))
setGeneric("countRead", function(object,
                                 target = "both",
                                 node = "",
                                 ...)
  standardGeneric("countRead"))
setGeneric("calcReadStats", function(object,
                                     target = "both",
                                     q = NULL,
                                     ...)
  standardGeneric("calcReadStats"))
setGeneric("setParents", function(object,
                                  parents,
                                  flip = FALSE,
                                  mono = FALSE,
                                  bi = TRUE,
                                  ...)
  standardGeneric("setParents"))
setGeneric("getParents", function(object,
                                  ...)
  standardGeneric("getParents"))
setGeneric("swapAlleles", function(object, allele, ...)
  standardGeneric("swapAlleles"))
setGeneric("setCallFilter", function(object,
                                     dp_count = c(0, Inf),
                                     ref_count = c(0, Inf),
                                     alt_count = c(0, Inf),
                                     norm_dp_count = c(0, Inf),
                                     norm_ref_count = c(0, Inf),
                                     norm_alt_count = c(0, Inf),
                                     scan_dp_qtile = c(0, 1),
                                     scan_ref_qtile = c(0, 1),
                                     scan_alt_qtile = c(0, 1),
                                     snp_dp_qtile = c(0, 1),
                                     snp_ref_qtile = c(0, 1),
                                     snp_alt_qtile = c(0, 1),
                                     omit_geno = NULL,
                                     ...)
           standardGeneric("setCallFilter"))
setGeneric("setScanFilter", function(object,
                                     id,
                                     missing = 1,
                                     het = c(0, 1),
                                     mac = 0,
                                     maf = 0,
                                     ad_ref = c(0, Inf),
                                     ad_alt = c(0, Inf),
                                     dp = c(0, Inf),
                                     mean_ref = c(0, Inf),
                                     mean_alt = c(0, Inf),
                                     sd_ref = Inf,
                                     sd_alt = Inf,
                                     ...)
           standardGeneric("setScanFilter"))
setGeneric("setSnpFilter", function(object,
                                    id,
                                    missing = 1,
                                    het = c(0, 1),
                                    mac = 0,
                                    maf = 0,
                                    ad_ref = c(0, Inf),
                                    ad_alt = c(0, Inf),
                                    dp = c(0, Inf),
                                    mean_ref = c(0, Inf),
                                    mean_alt = c(0, Inf),
                                    sd_ref = Inf,
                                    sd_alt = Inf,
                                    ...)
           standardGeneric("setSnpFilter"))
setGeneric("setInfoFilter", function(object,
                                     mq = 0 ,
                                     fs = Inf,
                                     qd = 0,
                                     sor = Inf,
                                     mqranksum = c(-Inf, Inf),
                                     readposranksum = c(-Inf, Inf),
                                     baseqranksum = c(-Inf, Inf),
                                     ...)
  standardGeneric("setInfoFilter"))
setGeneric("resetCallFilters", function(object, ...)
  standardGeneric("resetCallFilters"))
setGeneric("resetScanFilters", function(object, ...)
  standardGeneric("resetScanFilters"))
setGeneric("resetSnpFilters", function(object, ...)
  standardGeneric("resetSnpFilters"))
setGeneric("resetFilters", function(object, ...)
  standardGeneric("resetFilters"))
setGeneric("thinMarker", function(object, range = 150, ...)
  standardGeneric("thinMarker"))
setGeneric("subsetGDS", function(object,
                                 out_fn,
                                 snp_incl,
                                 scan_incl,
                                 incl_parents = TRUE,
                                 ...)
  standardGeneric("subsetGDS"))
setGeneric("setFiltGenotype", function(object, ...)
  standardGeneric("setFiltGenotype"))
setGeneric("setRawGenotype", function(object, ...)
  standardGeneric("setRawGenotype"))
setGeneric("replaceGDSdata", function(object,
                                      target,
                                      node = "filt",
                                      ...)
  standardGeneric("replaceGDSdata"))
setGeneric("clean", function(object,
                             scheme,
                             chr,
                             recomb_rate = 0.04,
                             error_rate = 0.0025,
                             call_threshold = 0.9,
                             het_parent = FALSE,
                             optim = TRUE,
                             iter = 2,
                             n_threads = NULL,
                             ...)
  standardGeneric("clean"))

setGeneric("initScheme", function(object, crosstype, mating, ...)
  standardGeneric("initScheme"))
setGeneric("addScheme", function(object, crosstype, mating, pop_size, ...)
  standardGeneric("addScheme"))
setGeneric("showScheme", function(object, ...)
  standardGeneric("showScheme"))
