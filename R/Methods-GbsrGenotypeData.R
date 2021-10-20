setMethod("show",
          "GbsrGenotypeData",
          function(object) {
            message('Data in GDS file...')
            message('GDS file name')
            message(object@data@handler$filename)
            if (isOpenGDS(object)) {
              print(object@data)
            } else {
              message('Connection to GDS file has been closed.')
            }
            message('-----------------------------------------------------')
            message('SnpAnnotationDataSet')
            print(object@snpAnnot)
            message('-----------------------------------------------------')
            message('ScanAnnotationDataSet')
            print(object@scanAnnot)
          })

# Internal function to replace read count data to publish to a VCF file.
.insertAnnot <- function(out_gds, in_gds, out_fmt, out_info) {
  info_node <- gdsfmt::index.gdsn(in_gds, "annotation/info")
  fmt_node <- gdsfmt::index.gdsn(in_gds, "annotation/format")
  new_info <- gdsfmt::index.gdsn(out_gds, "annotation/info")
  ls_gdsn <- gdsfmt::ls.gdsn(info_node)
  for (gdsn_i in ls_gdsn) {
    if (!is.null(out_info)) {
      if (!gdsn_i %in% out_info) {
        next
      }
    }
    i_node <- gdsfmt::index.gdsn(info_node, gdsn_i)
    gdsfmt::copyto.gdsn(new_info, i_node)
  }
  
  n_snp <-
    gdsfmt::objdesp.gdsn(gdsfmt::index.gdsn(in_gds, "genotype"))$dim[2]
  
  ls_gdsn <- gdsfmt::ls.gdsn(fmt_node)
  new_fmt <- gdsfmt::index.gdsn(out_gds, "annotation/format")
  for (gdsn_i in ls_gdsn) {
    if (!is.null(out_fmt)) {
      if (!gdsn_i %in% out_fmt) {
        next
      }
    }
    old_data_gdsn <-
      gdsfmt::index.gdsn(fmt_node, paste0(gdsn_i, "/data"))
    old_data <- gdsfmt::read.gdsn(old_data_gdsn)
    old_data[is.na(old_data)] <- "."
    n_num <- ncol(old_data) / n_snp
    
    if (n_num > 1) {
      if (gdsn_i == "HAP") {
        sep_symbol <- "/"
      } else {
        sep_symbol <- ","
      }
      old_data <- apply(old_data, 1, function(x) {
        x <- as.character(x)
        x <- matrix(x, nrow = n_num)
        return(apply(x, 2, paste, collapse = sep_symbol))
      })
      old_data <- t(old_data)
    }
    
    new_tmp <-
      gdsfmt::addfolder.gdsn(new_fmt, gdsn_i, replace = TRUE)
    new_data <- gdsfmt::add.gdsn(
      new_tmp,
      "data",
      old_data,
      storage = "string",
      replace = TRUE
    )
    
    old_attr <-
      gdsfmt::get.attr.gdsn(gdsfmt::index.gdsn(fmt_node, gdsn_i))
    if (!is.null(old_attr)) {
      for (i in seq_along(old_attr))
        gdsfmt::put.attr.gdsn(new_tmp, name = names(old_attr)[i], val = old_attr[[i]])
    }
    old_attr <- gdsfmt::get.attr.gdsn(old_data_gdsn)
    if (!is.null(old_attr)) {
      for (i in seq_along(old_attr))
        gdsfmt::put.attr.gdsn(new_data, name = names(old_attr)[i], val = old_attr[[i]])
    }
  }
}

.insertHaplotype <- function(out_gds, in_gds) {
  new_fmt <- gdsfmt::index.gdsn(out_gds, "annotation/format")
  new_tmp <- gdsfmt::addfolder.gdsn(new_fmt, "HAP", replace = TRUE)
  gdsfmt::put.attr.gdsn(new_tmp, "Number", val = "1")
  gdsfmt::put.attr.gdsn(new_tmp, "Type", val = "Sting")
  gdsfmt::put.attr.gdsn(new_tmp, "Description", val = "Haplotype estimated by GBScleanR.")
  hap_gdsn <- gdsfmt::index.gdsn(in_gds, "estimated.haplotype")
  hap <- gdsfmt::read.gdsn(hap_gdsn)
  gdsfmt::add.gdsn(
    new_tmp,
    "data",
    hap,
    storage = "bit6",
    replace = TRUE
  )
  
  new_info <- gdsfmt::index.gdsn(out_gds, "annotation/info")
  
  pgt_gdsn <- gdsfmt::index.gdsn(in_gds, "parents.genotype")
  pgt <- gdsfmt::readex.gdsn(pgt_gdsn)
  pgt <- apply(pgt, 2, paste, collapse = ",")
  
  new_pgt <- gdsfmt::add.gdsn(
    new_info,
    "PGT",
    pgt,
    storage = "string",
    replace = TRUE
  )
  gdsfmt::put.attr.gdsn(new_pgt, "Number", val = "1")
  gdsfmt::put.attr.gdsn(new_pgt, "Type", val = "String")
  gdsfmt::put.attr.gdsn(new_pgt, "Description", val = "Genotype of each haplotype estimated by GBScleanR.")
}

.recalcDP <- function(object) {
  reads <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(object@data@handler,
                                                "annotation/format/AD/data"))
  dp <- reads[, c(TRUE, FALSE)] + reads[, c(FALSE, TRUE)]
  filtdp <- gdsfmt::add.gdsn(
    gdsfmt::index.gdsn(object@data@handler,
                       "annotation/format/DP"),
    name = "filt.data",
    val = dp,
    storage = "int16",
    replace = TRUE
  )
  dpdata <- gdsfmt::index.gdsn(object@data@handler,
                               "annotation/format/DP/data")
  attr <- gdsfmt::get.attr.gdsn(dpdata)
  if (!is.null(attr)) {
    for (i in seq_along(attr))
      gdsfmt::put.attr.gdsn(filtdp,
                            name = names(attr)[i],
                            val = attr[[i]])
  }
  gdsfmt::delete.gdsn(dpdata)
  gdsfmt::rename.gdsn(
    gdsfmt::index.gdsn(object@data@handler,
                       "annotation/format/DP/filt.data"),
    "data"
  )
}


.replaceGDSdata <- function(object, target, node) {
  if ("sample.id" %in% target) {
    id <- getScanID(object, valid = FALSE)
    gdsfmt::add.gdsn(
      object@data@handler,
      "sample.id",
      id,
      "string",
      replace = TRUE
    )
  }
  
  if ("genotype" %in% target) {
    if (node == "filt" &
        "filt.genotype" %in% gdsfmt::ls.gdsn(object@data@handler)) {
      gt <- gdsfmt::index.gdsn(object@data@handler, "genotype")
      gt_attr <- gdsfmt::get.attr.gdsn(gt)
      gdsfmt::delete.gdsn(gt)
      
      newgt <-
        gdsfmt::index.gdsn(object@data@handler, "filt.genotype")
      gdsfmt::rename.gdsn(newgt, "genotype")
      gdsfmt::put.attr.gdsn(newgt, names(gt_attr)[1], gt_attr[[1]])
      
    } else if (node == "cor" &
               "corrected.genotype" %in% gdsfmt::ls.gdsn(object@data@handler)) {
      gt <- gdsfmt::index.gdsn(object@data@handler, "genotype")
      gt_attr <- gdsfmt::get.attr.gdsn(gt)
      gdsfmt::delete.gdsn(gt)
      newgt <-
        gdsfmt::index.gdsn(object@data@handler, "corrected.genotype")
      gdsfmt::rename.gdsn(newgt, "genotype")
      gdsfmt::put.attr.gdsn(newgt, names(gt_attr)[1], gt_attr[[1]])
    } else {
      warning('Nothing to replace.')
    }
  }
  if ("ad" %in% target) {
    ad_gdsn <-
      gdsfmt::index.gdsn(object@data@handler, "annotation/format/AD")
    ls_gdsn <- gdsfmt::ls.gdsn(ad_gdsn)
    if ("filt.data" %in% ls_gdsn) {
      ad <- gdsfmt::index.gdsn(ad_gdsn, "data")
      gdsfmt::delete.gdsn(ad)
      ad_gdsn <-
        gdsfmt::index.gdsn(object@data@handler, "annotation/format/AD")
      fad <- gdsfmt::index.gdsn(ad_gdsn, "filt.data")
      gdsfmt::rename.gdsn(fad, "data")
      gdsfmt::readmode.gdsn(fad)
      if (length(ls_gdsn) > 2) {
        ls_gdsn <- ls_gdsn[!ls_gdsn %in% c("data", "filt.data")]
        for (i in ls_gdsn) {
          ad_gdsn <-
            gdsfmt::index.gdsn(object@data@handler, "annotation/format/AD")
          gdsfmt::delete.gdsn(gdsfmt::index.gdsn(ad_gdsn, i))
        }
      }
      
    } else {
      warning('Nothing to replace.')
    }
  }
}


# Internally used function to flip genotype data based on the flipped marker information.
# Flipped markers are markers where the alleles expected as reference allele are called as
# alternative allele.
.flipGeno <- function(object, var) {
  flipped <- getFlipped(object, valid = FALSE)
  gt_node <- gdsfmt::index.gdsn(object@data@handler, var)
  gt <- gdsfmt::read.gdsn(gt_node)
  gdsfmt::compression.gdsn(gt_node, "")
  flip_gt <- gt[, flipped]
  flip_0 <- flip_gt == 0
  flip_2 <- flip_gt == 2
  flip_gt[flip_0] <- 2
  flip_gt[flip_2] <- 0
  gt[, flipped] <- flip_gt
  gt_attr <-
    gdsfmt::get.attr.gdsn(gdsfmt::index.gdsn(object@data@handler, var))
  gt_gdsn <- gdsfmt::add.gdsn(object@data@handler,
                              var,
                              gt,
                              "bit2",
                              replace = TRUE)
  if (!is.null(gt_attr)) {
    for (i in seq_along(gt_attr)) {
      gdsfmt::put.attr.gdsn(node = gt_gdsn,
                            name = names(gt_attr)[i],
                            val = gt_attr[[i]])
    }
  }
}

# Internally used function to flip AD data based on the flipped marker information.
# Flipped markers are markers where the alleles expected as reference allele are called as
# alternative allele.
.flipAD <- function(object, var) {
  flipped <- getFlipped(object, valid = FALSE)
  ad <-
    gdsfmt::index.gdsn(object@data@handler, "annotation/format/AD")
  data <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(ad, var))
  gdsfmt::compression.gdsn(data, "")
  ref <- data[, c(TRUE, FALSE)]
  alt <- data[, c(FALSE, TRUE)]
  tmp <- ref[, flipped]
  ref[, flipped] <- alt[, flipped]
  alt[, flipped] <- tmp
  data[, c(TRUE, FALSE)] <- ref
  data[, c(FALSE, TRUE)] <- alt
  gdsfmt::add.gdsn(
    node = ad,
    name = var,
    val = data,
    storage = "int16",
    replace = TRUE
  )
}

.flipData <- function(object) {
  if (!haveFlipped(object)) {
    stop('Nothing to flip.')
  }
  
  # Genotypes
  .flipGeno(object, "genotype")
  ls_gdsn <- gdsfmt::ls.gdsn(object@data@handler)
  if ("filt.genotype" %in% ls_gdsn) {
    .flipGeno(object, "filt.genotype")
  }
  
  # Alleles
  allele <- paste(getAlleleA(object, valid = FALSE),
                  getAlleleB(object, valid = FALSE),
                  sep = "/")
  allele_node <- gdsfmt::index.gdsn(object@data@handler, "snp.allele")
  gdsfmt::compression.gdsn(allele_node, "")
  gdsfmt::add.gdsn(
    node = object@data@handler,
    name = "snp.allele",
    val = allele,
    storage = "string",
    replace = TRUE
  )
  
  # Allele reads
  .flipAD(object, "data")
  
  # Filtered allele reads
  ls_gdsn <-
    gdsfmt::ls.gdsn(gdsfmt::index.gdsn(object@data@handler,
                                       "annotation/format/AD"))
  if ("filt.data" %in% ls_gdsn) {
    .flipAD(object, "filt.data")
  }
  
  object@snpAnnot$flipped <- NULL
  return(object)
}

#' @rdname gbsrGDS2VCF
setMethod("gbsrGDS2VCF",
          "GbsrGenotypeData",
          function(object,
                   out_fn,
                   node,
                   valid,
                   out_fmt,
                   out_info) {
            
            have_hap <-
              "estimated.haplotype" %in% gdsfmt::ls.gdsn(object@data@handler)
            if(have_hap){
              hap_dim <- gdsfmt::objdesp.gdsn(gdsfmt::index.gdsn(object@data@handler, "estimated.haplotype"))$dim
              if(length(hap_dim) != 3)
              {
                have_hap <- FALSE
              }
            }
            suppressMessages(closeGDS(object))
            gds_file <- object@data@filename
            if (!grepl(".gds$", gds_file)) {
              gds_file <- paste0(gds_file, ".gds")
            }
            tmp_gds <- sub(".gds", ".tmp.gds", gds_file)
            file.copy(gds_file, tmp_gds, overwrite = TRUE)
            tmp_gds <- suppressMessages(loadGDS(tmp_gds))
            tmp_gds@snpAnnot <- object@snpAnnot
            tmp_gds@scanAnnot <- object@scanAnnot
            
            if (haveFlipped(object)) {
              tmp_gds <- .flipData(tmp_gds)
            }
            if (node %in% c("filt", "cor")) {
              suppressWarnings(.replaceGDSdata(tmp_gds, "genotype", node))
              suppressWarnings(.replaceGDSdata(tmp_gds, "ad", "filt"))
              if ("DP" %in% out_fmt) {
                .recalcDP(tmp_gds)
              }
            }
            
            if (valid == TRUE) {
              subset_gds <- sub(".gds", ".tmp.subset.gds", gds_file)
              tmp_subset_gds <-
                suppressMessages(subsetGDS(
                  object = tmp_gds,
                  out_fn = subset_gds,
                  incl_parents = TRUE
                ))
              if (is.null(tmp_subset_gds)) {
                tmp_gds_file <- tmp_gds@data@filename
                
              } else {
                tmp_subset_file <- tmp_subset_gds@data@filename
                suppressMessages(closeGDS(tmp_subset_gds))
                tmp_gds_file <- tmp_subset_gds@data@filename
              }
            } else {
              tmp_gds_file <- tmp_gds@data@filename
            }
            suppressMessages(closeGDS(tmp_gds))
            
            if (!grepl(".vcf$", out_fn)) {
              out_fn <- paste0(out_fn, ".vcf")
            }
            out_fn_tmp <- sub(".vcf", ".tmpout.gds", out_fn)
            SeqArray::seqSNP2GDS(gds.fn = tmp_gds_file,
                                 out.fn = out_fn_tmp,
                                 verbose = FALSE)
            
            out_gds <-
              gdsfmt::openfn.gds(out_fn_tmp, readonly = FALSE)
            if(have_hap) {
              in_gds <- gdsfmt::openfn.gds(gds_file, readonly = FALSE)
              .insertHaplotype(out_gds, in_gds)
              gdsfmt::closefn.gds(in_gds)
            }
            
            in_gds <- gdsfmt::openfn.gds(tmp_gds_file, readonly = FALSE)
            .insertAnnot(out_gds, in_gds, out_fmt, out_info)
            gdsfmt::closefn.gds(in_gds)
            gdsfmt::closefn.gds(out_gds)
            
            SeqArray::seqGDS2VCF(gdsfile = out_fn_tmp,
                                 vcf.fn = out_fn,
                                 verbose = FALSE)
            on.exit({
              if (valid) {
                if (!is.null(tmp_subset_gds)) {
                  suppressWarnings(file.remove(tmp_subset_file))
                }
              }
              suppressWarnings(file.remove(tmp_gds_file))
              suppressWarnings(file.remove(out_fn_tmp))
              
            })
          })

#' @rdname isOpenGDS
setMethod("isOpenGDS",
          "GbsrGenotypeData",
          function(object) {
            tryout <-
              try(gdsfmt::openfn.gds(filename = object@data@filename, readonly = FALSE),
                  silent = TRUE)
            if (inherits(tryout, "gds.class")) {
              gdsfmt::closefn.gds(tryout)
              return(FALSE)
              
            } else if (grepl("has been created or opened.", tryout[1])) {
              return(TRUE)
              
            } else {
              warning(tryout[1])
              return(NULL)
            }
          })

#' @rdname closeGDS
setMethod("closeGDS",
          "GbsrGenotypeData",
          function(object) {
            gdsfmt::closefn.gds(object@data@handler)
            message('The connection to the GDS file was closed.')
          })

#' @rdname saveSnpAnnot
#' @importFrom Biobase pData
setMethod("saveSnpAnnot",
          "GbsrGenotypeData",
          function(object) {
            .gds_decomp(object)
            new_node <- gdsfmt::add.gdsn(
              node = object@data@handler,
              name = "snpAnnot",
              val = pData(object@snpAnnot),
              replace = TRUE
            )
            .gds_decomp(object)
            .gds_comp(object)
            return(object)
          })

#' @rdname saveScanAnnot
#' @importFrom Biobase pData
setMethod("saveScanAnnot",
          "GbsrGenotypeData",
          function(object) {
            .gds_decomp(object)
            new_node <- gdsfmt::add.gdsn(
              node = object@data@handler,
              name = "scanAnnot",
              val = pData(object@scanAnnot),
              replace = TRUE
            )
            .gds_decomp(object)
            .gds_comp(object)
            return(object)
          })

#' @rdname loadSnpAnnot
#' @importFrom Biobase pData<-
setMethod("loadSnpAnnot",
          "GbsrGenotypeData",
          function(object) {
            ls_gdsn <-
              snpAnnot_node <- gdsfmt::ls.gdsn(node = object@data@handler)
            if ("snpAnnot" %in% ls_gdsn) {
              snpAnnot_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                                  path = "snpAnnot")
              ls_gdsn <- gdsfmt::ls.gdsn(node = snpAnnot_node)
              for (i in ls_gdsn) {
                gdsfmt::readmode.gdsn(node = index.gdsn(node = snpAnnot_node, path = i))
              }
              snpAnnot <- gdsfmt::read.gdsn(node = snpAnnot_node)
              pData(object@snpAnnot) <- snpAnnot
            } else {
              message('No data of snpAnnot in the GDS file.')
            }
            return(object)
          })

#' @rdname loadScanAnnot
#' @importFrom Biobase pData<-
setMethod("loadScanAnnot",
          "GbsrGenotypeData",
          function(object) {
            ls_gdsn <-
              snpAnnot_node <- gdsfmt::ls.gdsn(node = object@data@handler)
            if ("scanAnnot" %in% ls_gdsn) {
              scanAnnot_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                                   path = "scanAnnot")
              ls_gdsn <- gdsfmt::ls.gdsn(node = scanAnnot_node)
              for (i in ls_gdsn) {
                gdsfmt::readmode.gdsn(node = index.gdsn(node = scanAnnot_node, path = i))
              }
              scanAnnot <- gdsfmt::read.gdsn(node = scanAnnot_node)
              pData(object@scanAnnot) <- scanAnnot
            } else {
              message('No data of snpAnnot in the GDS file.')
            }
            return(object)
          })

#' @rdname nsnp
setMethod("nsnp",
          "GbsrGenotypeData",
          function(object, valid = TRUE) {
            if (valid) {
              out <- sum(getValidSnp(object))
            } else {
              out <- nrow(object@snpAnnot)
            }
            return(out)
          })

#' @rdname nscan
setMethod("nscan",
          "GbsrGenotypeData",
          function(object, valid = TRUE) {
            if (valid) {
              out <- sum(getValidScan(object))
            } else {
              out <- nrow(object@scanAnnot)
            }
            return(out)
          })

#' @rdname getValidSnp
setMethod("getValidSnp",
          "GbsrGenotypeData",
          function(object) {
            return(object@snpAnnot$validMarker)
          })

#' @rdname setValidSnp
setMethod("setValidSnp",
          "GbsrGenotypeData",
          function(object, new, update) {
            if (!missing(new)) {
              if (any(is.na(new))) {
                stop('NA is not allowed for a logical vector "new".')
              }
              object@snpAnnot$validMarker <- new
            } else if (!missing(update)) {
              if (any(is.na(update))) {
                stop('NA is not allowed for a logical vector "update".')
              }
              update <-
                object@snpAnnot$validMarker[getValidSnp(object)] & update
              object@snpAnnot$validMarker[getValidSnp(object)] <-
                update
            }
            return(object)
          })

#' @rdname getValidScan
setMethod("getValidScan",
          "GbsrGenotypeData",
          function(object, parents = FALSE) {
            if (parents == "only") {
              return(object@scanAnnot$parents != 0)
            }
            if (parents) {
              out <- object@scanAnnot$validScan
              out[object@scanAnnot$parents != 0] <- TRUE
              return(out)
            } else {
              return(object@scanAnnot$validScan)
            }
          })

#' @rdname setValidScan
setMethod("setValidScan",
          "GbsrGenotypeData",
          function(object, new, update) {
            if (!missing(update)) {
              if (any(is.na(update))) {
                stop('NA is not allowed for a logical vector "update".')
              }
              update <-
                object@scanAnnot$validScan[getValidScan(object)] & update
              object@scanAnnot$validScan[getValidScan(object)] <-
                update
            } else if (!missing(new)) {
              if (any(is.na(new))) {
                stop('NA is not allowed for a logical vector "new".')
              }
              object@scanAnnot$validScan <- new
            }
            return(object)
          })

# This method is internally used.
.setFlipped <- function(object, flipped) {
  if (any(is.na(flipped))) {
    stop('NA is not allowed for a logical vector "flipped".')
  }
  valid_snp <- getValidSnp(object)
  if (length(valid_snp) == length(flipped)) {
    object@snpAnnot$flipped <- flipped
  } else {
    if (nsnp(object) == length(flipped)) {
      tmp_flipped <- rep(FALSE, nsnp(object, valid = FALSE))
      tmp_flipped[valid_snp] <- flipped
      object@snpAnnot$flipped <- tmp_flipped
    } else {
      stop('The length of "flipped" does not match with the number of SNPs.')
    }
  }
  return(object)
}


#' @rdname getFlipped
setMethod("getFlipped",
          "GbsrGenotypeData",
          function(object, valid = TRUE) {
            out <- object@snpAnnot$flipped
            if (is.null(out)) {
              message('No data of flipped genotype markers.')
              return(NULL)
            }
            if (valid) {
              out <- out[getValidSnp(object)]
            }
            return(out)
          })


#' @rdname haveFlipped
setMethod("haveFlipped",
          "GbsrGenotypeData",
          function(object) {
            return(!is.null(object@snpAnnot$flipped))
          })

#' @rdname getRead
setMethod("getRead",
          "GbsrGenotypeData",
          function(object, chr, node, parents) {
            ad_gdsn <-
              gdsfmt::index.gdsn(object@data@handler, "annotation/format/AD")
            ls_gdsn <- gdsfmt::ls.gdsn(ad_gdsn)
            if (node == "filt" & "filt.data" %in% ls_gdsn) {
              path <- "filt.data"
            } else {
              path <- "data"
            }
            ad_node <- gdsfmt::index.gdsn(ad_gdsn, path)
            valid_samples <- getValidScan(object)
            if (parents == "only") {
              valid_samples <- object@scanAnnot$parents != 0
            } else if (parents) {
              valid_samples[object@scanAnnot$parents != 0] <- TRUE
            }
            
            valid_snp <- getValidSnp(object)
            if (!is.null(chr)) {
              valid_snp <- valid_snp & getChromosome(object, valid = FALSE) == chr
            }
            ad <- gdsfmt::readex.gdsn(node = ad_node,
                                      sel = list(valid_samples,
                                                 rep(valid_snp, each = 2)))
            ref <- ad[, c(TRUE, FALSE)]
            alt <- ad[, c(FALSE, TRUE)]
            
            if (haveFlipped(object)) {
              flipped <- getFlipped(object, valid = FALSE)[valid_snp]
              flipped_ref <- ref[, flipped]
              ref[, flipped] <- alt[, flipped]
              alt[, flipped] <- flipped_ref
            }
            
            rownames(ref) <-
              getScanID(object, valid = FALSE)[valid_samples]
            rownames(alt) <-
              getScanID(object, valid = FALSE)[valid_samples]
            return(list(ref = ref, alt = alt))
          })


#' @rdname getGenotype
setMethod("getGenotype",
          "GbsrGenotypeData",
          function(object, chr, node, parents) {
            ls_gdsn <- gdsfmt::ls.gdsn(object@data@handler)
            path <- NULL
            if (node == "parents") {
              if ("parents.genotype" %in% ls_gdsn) {
                path <- "parents.genotype"
              } else {
                stop("No estimated parents genotype data.")
              }
            }
            if (node == "cor") {
              if ("corrected.genotype" %in% ls_gdsn) {
                path <- "corrected.genotype"
              } else {
                stop("No corrected genotype data.")
              }
            }
            if (node == "filt") {
              if ("filt.genotype" %in% ls_gdsn) {
                path <- "filt.genotype"
              } else {
                stop("No filtered genotype data.")
              }
            }
            if (node == "raw") {
              path <- "genotype"
            }
            if (is.null(path)) {
              stop("Please specify a valid node.")
            }
            genotype_node <-
              gdsfmt::index.gdsn(node = object@data@handler,
                                 path = path)
            valid_markers <- getValidSnp(object)
            if (!is.null(chr)) {
              valid_marker[getChromosome(object) != chr] <- FALSE
            }
            
            if (node == "parents") {
              n_row <- gdsfmt::objdesp.gdsn(genotype_node)$dim[1]
              sel <- list(rep(TRUE, n_row),
                          valid_markers)
              genotype <- gdsfmt::readex.gdsn(node = genotype_node,
                                              sel = sel)
            } else {
              valid_samples <- getValidScan(object)
              if (parents == "only") {
                valid_samples <- object@scanAnnot$parents != 0
              } else if (parents) {
                valid_samples[object@scanAnnot$parents != 0] <- TRUE
              }
              sel <- list(valid_samples, valid_markers)
              genotype <- gdsfmt::readex.gdsn(node = genotype_node,
                                              sel = sel)
              rownames(genotype) <- getScanID(object,
                                              valid = FALSE)[valid_samples]
            }
            genotype[genotype == 3] <- NA
            return(genotype)
          })


#' @rdname getHaplotype
setMethod("getHaplotype",
          "GbsrGenotypeData",
          function(object, chr, parents) {
            ls_gdsn <- gdsfmt::ls.gdsn(object@data@handler)
            if ("estimated.haplotype" %in% ls_gdsn) {
              path <- "estimated.haplotype"
            } else {
              stop('No haplotype data, Run clean() to estimate haplotype.')
            }
            haplotype_node <-
              gdsfmt::index.gdsn(node = object@data@handler,
                                 path = path)
            valid_samples <- getValidScan(object)
            if (parents == "only") {
              valid_samples <- object@scanAnnot$parents != 0
            } else if (parents) {
              valid_samples[object@scanAnnot$parents != 0] <- TRUE
            }
            valid_markers <- getValidSnp(object)
            sel <- list(valid_samples, rep(valid_markers, each = 2))
            haplotype <- gdsfmt::readex.gdsn(node = haplotype_node,
                                             sel = sel)
            haplotype[haplotype == 0] <- NA
            haplotype <-
              array(t(haplotype), c(2, sum(valid_markers), sum(valid_samples)))
            if (!is.null(chr)) {
              haplotype <- haplotype[, getChromosome(object) == chr,]
            }
            return(haplotype)
          })

#' @rdname getChromosome
setMethod("getChromosome",
          "GbsrGenotypeData",
          function(object,
                   valid = TRUE,
                   levels = FALSE,
                   name = FALSE) {
            if (name) {
              out <- object@snpAnnot$chromosome.name
            } else {
              out <- object@snpAnnot$chromosome
            }
            if (valid) {
              out <- out[getValidSnp(object)]
            }
            if (levels) {
              out <- unique(out)
            }
            return(out)
          })


#' @rdname getPosition
setMethod("getPosition",
          "GbsrGenotypeData",
          function(object, valid = TRUE) {
            out <- object@snpAnnot$position
            if (valid) {
              out <- out[getValidSnp(object)]
            }
            return(out)
          })

#' @rdname getAlleleA
setMethod("getAlleleA",
          "GbsrGenotypeData",
          function(object, valid = TRUE) {
            out <- object@snpAnnot$alleleA
            if (haveFlipped(object)) {
              flipped <- getFlipped(object, valid = FALSE)
              b <- object@snpAnnot$alleleB
              out[flipped] <- b[flipped]
            }
            if (valid) {
              out <- out[getValidSnp(object)]
            }
            return(out)
          })

#' @rdname getAlleleB
setMethod("getAlleleB",
          "GbsrGenotypeData",
          function(object, valid = TRUE) {
            out <- object@snpAnnot$alleleB
            if (haveFlipped(object)) {
              flipped <- getFlipped(object, valid = FALSE)
              a <- object@snpAnnot$alleleA
              out[flipped] <- a[flipped]
            }
            if (valid) {
              out <- out[getValidSnp(object)]
            }
            return(out)
          })

#' @rdname getSnpID
setMethod("getSnpID",
          "GbsrGenotypeData",
          function(object, valid = TRUE) {
            out <- object@snpAnnot$snpID
            if (valid) {
              out <- out[getValidSnp(object)]
            }
            return(out)
          })

#' @rdname getScanID
setMethod("getScanID",
          "GbsrGenotypeData",
          function(object, valid = TRUE) {
            out <- object@scanAnnot$scanID
            if (valid) {
              out <- out[getValidScan(object)]
            }
            return(out)
          })

#' @rdname getPloidy
setMethod("getPloidy",
          "GbsrGenotypeData",
          function(object, valid = TRUE) {
            out <- object@snpAnnot$ploidy
            if (valid) {
              out <- out[getValidSnp(object)]
            }
            return(out)
          })

#' @rdname getInfo
setMethod("getInfo",
          "GbsrGenotypeData",
          function(object, var, valid = TRUE) {
            path <- paste0("annotation/info/", var)
            ls_gdsn <-
              snpAnnot_node <- gdsfmt::ls.gdsn(
                node = object@data@handler,
                recursive = TRUE,
                include.dirs = TRUE
              )
            if (path %in% ls_gdsn) {
              info_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                              path = path)
            } else {
              return(NULL)
            }
            
            if (valid) {
              out <-
                gdsfmt::readex.gdsn(node = info_node, sel = list(getValidSnp(object)))
            } else {
              out <- gdsfmt::read.gdsn(node = info_node)
            }
            return(out)
          })

#' @rdname getCountReadRef
setMethod("getCountReadRef",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE,
                   prop = FALSE) {
            if (target == "snp") {
              out <- object@snpAnnot$countReadRef
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
              if (prop) {
                out <- out / getCountRead(object, target = "snp", valid = valid)
              }
            } else {
              out <- object@scanAnnot$countReadRef
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
              if (prop) {
                out <- out / getCountRead(object, target = "scan", valid = valid)
              }
            }
            return(out)
          })

#' @rdname getCountReadAlt
setMethod("getCountReadAlt",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE,
                   prop = FALSE) {
            if (target == "snp") {
              out <- object@snpAnnot$countReadAlt
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
              if (prop) {
                out <- out / getCountRead(object, target = "snp", valid = valid)
              }
            } else {
              out <- object@scanAnnot$countReadAlt
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
              if (prop) {
                out <- out / getCountRead(object, target = "scan", valid = valid)
              }
            }
            return(out)
          })

#' @rdname getCountRead
setMethod("getCountRead",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE) {
            if (target == "snp") {
              out1 <- getCountReadRef(object, "snp", prop = FALSE, valid = valid)
              out2 <-
                getCountReadAlt(object, "snp", prop = FALSE, valid = valid)
              if (is.null(out1) | is.null(out2)) {
                return(NULL)
              }
              
            } else {
              out1 <- getCountReadRef(object, "scan", prop = FALSE, valid = valid)
              out2 <-
                getCountReadAlt(object, "scan", prop = FALSE, valid = valid)
              if (is.null(out1) | is.null(out2)) {
                return(NULL)
              }
            }
            out <- out1 + out2
            return(out)
          })

#' @rdname getCountGenoRef
setMethod("getCountGenoRef",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE,
                   prop = FALSE) {
            if (target == "snp") {
              out <- object@snpAnnot$countGenoRef
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
              if (prop) {
                out <- out / {
                  nscan(object, valid = valid) -
                    getCountGenoMissing(object, "snp", valid = valid)
                }
              }
            } else {
              out <- object@scanAnnot$countGenoRef
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
              if (prop) {
                out <- out / {
                  nsnp(object, valid = valid) -
                    getCountGenoMissing(object, "scan", valid = valid)
                }
              }
            }
            return(out)
          })

#' @rdname getCountGenoHet
setMethod("getCountGenoHet",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE,
                   prop = FALSE) {
            if (target == "snp") {
              out <- object@snpAnnot$countGenoHet
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
              if (prop) {
                out <- out / {
                  nscan(object, valid = valid) -
                    getCountGenoMissing(object, "snp", valid = valid)
                }
              }
            } else {
              out <- object@scanAnnot$countGenoHet
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
              if (prop) {
                out <- out / {
                  nsnp(object, valid = valid) -
                    getCountGenoMissing(object, "scan", valid = valid)
                }
              }
            }
            return(out)
          })

#' @rdname getCountGenoAlt
setMethod("getCountGenoAlt",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE,
                   prop = FALSE) {
            if (target == "snp") {
              out <- object@snpAnnot$countGenoAlt
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
              if (prop) {
                out <- out / {
                  nscan(object, valid = valid) -
                    getCountGenoMissing(object, "snp", valid = valid)
                }
              }
            } else {
              out <- object@scanAnnot$countGenoAlt
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
              if (prop) {
                out <- out / {
                  nsnp(object, valid = valid) -
                    getCountGenoMissing(object, "scan", valid = valid)
                }
              }
            }
            return(out)
          })

#' @rdname getCountGenoMissing
setMethod("getCountGenoMissing",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE,
                   prop = FALSE) {
            if (target == "snp") {
              out <- object@snpAnnot$countGenoMissing
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
              if (prop) {
                out <- out / nscan(object, valid = valid)
              }
            } else {
              out <- object@scanAnnot$countGenoMissing
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
              if (prop) {
                out <- out / nsnp(object, valid = valid)
              }
            }
            return(out)
          })

#' @rdname getCountAlleleRef
setMethod("getCountAlleleRef",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE,
                   prop = FALSE) {
            if (target == "snp") {
              out <- object@snpAnnot$countAlleleRef
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
              if (prop) {
                out <- out / {
                  nscan(object, valid = valid) * 2 -
                    getCountAlleleMissing(object, "snp", valid = valid)
                }
              }
            } else {
              out <- object@scanAnnot$countAlleleRef
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
              if (prop) {
                out <- out / {
                  nsnp(object, valid = valid) * 2 -
                    getCountAlleleMissing(object, "scan", valid = valid)
                }
              }
            }
            return(out)
          })

#' @rdname getCountAlleleAlt
setMethod("getCountAlleleAlt",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE,
                   prop = FALSE) {
            if (target == "snp") {
              out <- object@snpAnnot$countAlleleAlt
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
              if (prop) {
                out <- out / {
                  nscan(object, valid = valid) * 2 -
                    getCountAlleleMissing(object, "snp", valid = valid)
                }
              }
            } else {
              out <- object@scanAnnot$countAlleleAlt
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
              if (prop) {
                out <- out / {
                  nsnp(object, valid = valid) * 2 -
                    getCountAlleleMissing(object, "scan", valid = valid)
                }
              }
            }
            return(out)
          })

#' @rdname getCountAlleleMissing
setMethod("getCountAlleleMissing",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE,
                   prop = FALSE) {
            if (target == "snp") {
              out <- object@snpAnnot$countAlleleMissing
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
              if (prop) {
                out <- out / {
                  nscan(object, valid = valid) * 2
                }
              }
            } else {
              out <- object@scanAnnot$countAlleleMissing
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
              if (prop) {
                out <- out / {
                  nsnp(object, valid = valid) * 2
                }
              }
            }
            return(out)
          })

#' @rdname getMeanReadRef
setMethod("getMeanReadRef",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE) {
            if (target == "snp") {
              out <- object@snpAnnot$meanReadRef
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
            } else {
              out <- object@scanAnnot$meanReadRef
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
            }
            return(out)
          })

#' @rdname getMeanReadAlt
setMethod("getMeanReadAlt",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE) {
            if (target == "snp") {
              out <- object@snpAnnot$meanReadAlt
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
            } else {
              out <- object@scanAnnot$meanReadAlt
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
            }
            return(out)
          })

#' @rdname getSDReadRef
setMethod("getSDReadRef",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE) {
            if (target == "snp") {
              out <- object@snpAnnot$sdReadRef
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
            } else {
              out <- object@scanAnnot$sdReadRef
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
            }
            return(out)
          })

#' @rdname getSDReadAlt
setMethod("getSDReadAlt",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE) {
            if (target == "snp") {
              out <- object@snpAnnot$sdReadAlt
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidSnp(object)]
              }
            } else {
              out <- object@scanAnnot$sdReadAlt
              if (is.null(out)) {
                return(NULL)
              }
              if (valid) {
                out <- out[getValidScan(object)]
              }
            }
            return(out)
          })

#' @rdname getQtileReadRef
#' @importFrom Biobase pData
setMethod("getQtileReadRef",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   q = 0.5,
                   valid = TRUE) {
            if (target == "snp") {
              pdata <- pData(object@snpAnnot)
              if (q == "all") {
                index <- grepl("qtileReadRef", names(pdata))
              } else {
                index <- names(pdata) %in% paste0("qtileReadRef", q)
              }
              if (all(!index)) {
                return(NULL)
              }
              out <- pdata[, index]
              if (valid) {
                out <- subset(out, subset = getValidSnp(object))
              }
            } else {
              pdata <- pData(object@scanAnnot)
              if (q == "all") {
                index <- grepl("qtileReadRef", names(pdata))
              } else {
                index <- names(pdata) %in% paste0("qtileReadRef", q)
              }
              if (all(!index)) {
                return(NULL)
              }
              out <- pdata[, index]
              if (valid) {
                out <- subset(out, subset = getValidScan(object))
              }
            }
            if (length(out) == 1) {
              out <- unlist(out)
            }
            return(out)
          })

#' @rdname getQtileReadAlt
#' @importFrom Biobase pData
setMethod("getQtileReadAlt",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   q = 0.5,
                   valid = TRUE) {
            if (target == "snp") {
              pdata <- pData(object@snpAnnot)
              if (q == "all") {
                index <- grepl("qtileReadAlt", names(pdata))
              } else {
                index <- names(pdata) %in% paste0("qtileReadAlt", q)
              }
              if (all(!index)) {
                return(NULL)
              }
              out <- pdata[, index]
              if (valid) {
                out <- subset(out, subset = getValidSnp(object))
              }
            } else {
              pdata <- pData(object@scanAnnot)
              if (q == "all") {
                index <- grepl("qtileReadAlt", names(pdata))
              } else {
                index <- names(pdata) %in% paste0("qtileReadAlt", q)
              }
              if (all(!index)) {
                return(NULL)
              }
              out <- pdata[, index]
              if (valid) {
                out <- subset(out, subset = getValidScan(object))
              }
            }
            if (length(out) == 1) {
              out <- unlist(out)
            }
            return(out)
          })

#' @rdname getMAF
setMethod("getMAF",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE) {
            if (target == "snp") {
              out <- getCountAlleleRef(object, "snp", prop = TRUE, valid = valid)
              if (is.null(out)) {
                return(NULL)
              }
              out <- 0.5 - abs(out - 0.5)
            } else {
              out <- getCountAlleleRef(object, "scan", prop = TRUE, valid = valid)
              if (is.null(out)) {
                return(NULL)
              }
              out <- 0.5 - abs(out - 0.5)
            }
            return(out)
          })

#' @rdname getMAC
setMethod("getMAC",
          "GbsrGenotypeData",
          function(object,
                   target = "snp",
                   valid = TRUE) {
            if (target == "snp") {
              ac_ref <- getCountAlleleRef(object, "snp", prop = FALSE, valid = valid)
              ac_alt <-
                getCountAlleleAlt(object, "snp", prop = FALSE, valid = valid)
              out <- ac_ref
              if (is.null(out)) {
                return(NULL)
              }
              alt_minor <- ac_ref > ac_alt
              out[alt_minor] <- ac_alt[alt_minor]
            } else {
              ac_ref <- getCountAlleleRef(object, "scan", prop = FALSE, valid = valid)
              ac_alt <-
                getCountAlleleAlt(object, "scan", prop = FALSE, valid = valid)
              out <- ac_ref
              if (is.null(out)) {
                return(NULL)
              }
              alt_minor <- ac_ref > ac_alt
              out[alt_minor] <- ac_alt[alt_minor]
            }
            return(out)
          })

#' @rdname countGenotype
setMethod("countGenotype",
          "GbsrGenotypeData",
          function(object,
                   target = "both",
                   node = "") {
            ls_gdsn <- gdsfmt::ls.gdsn(object@data@handler)
            if (node == "cor" & "corrected.genotype" %in% ls_gdsn) {
              genotype_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                                  path = "corrected.genotype")
            } else if (node != "raw" &
                       object@data@genotypeVar == "filt.genotype") {
              genotype_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                                  path = "filt.genotype")
            } else {
              genotype_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                                  path = "genotype")
            }
            
            valid_markers <- getValidSnp(object)
            valid_scans <- getValidScan(object)
            if (node == "cor") {
              have_flipped <- FALSE
            } else {
              have_flipped <- haveFlipped(object)
            }
            
            if (have_flipped) {
              valid_flipped <- getFlipped(object, valid = TRUE)
            }
            sel <- list(valid_scans, valid_markers)
            
            # Counts per sample
            if (target %in% c("both", "scan")) {
              geno_table <- gdsfmt::apply.gdsn(
                node = genotype_node,
                margin = 1,
                selection = sel,
                FUN = function(x) {
                  if (have_flipped) {
                    ref <- x == 2
                    alt <- x == 0
                    x[valid_flipped &
                        ref] <- 0
                    x[valid_flipped &
                        alt] <- 2
                  }
                  x <-
                    factor(x, levels = 0:3)
                  return(as.integer(table(x)))
                },
                as.is = "list"
              )
              geno_table <- matrix(unlist(geno_table), nrow = 4)
              
              ## Summarize genotype call counts
              object@scanAnnot$countGenoRef <- NA
              object@scanAnnot$countGenoHet <- NA
              object@scanAnnot$countGenoAlt <- NA
              object@scanAnnot$countGenoMissing <- NA
              object@scanAnnot$countGenoRef[valid_scans] <-
                geno_table[3,]
              object@scanAnnot$countGenoHet[valid_scans] <-
                geno_table[2,]
              object@scanAnnot$countGenoAlt[valid_scans] <-
                geno_table[1,]
              object@scanAnnot$countGenoMissing[valid_scans] <-
                geno_table[4,]
              
              
              ## Summarize allele counts
              object@scanAnnot$countAlleleRef <- NA
              object@scanAnnot$countAlleleAlt <- NA
              object@scanAnnot$countAlleleMissing <- NA
              object@scanAnnot$countAlleleRef[valid_scans] <-
                colSums(geno_table * c(0, 1, 2, 0))
              object@scanAnnot$countAlleleAlt[valid_scans] <-
                colSums(geno_table * c(2, 1, 0, 0))
              object@scanAnnot$countAlleleMissing[valid_scans] <-
                colSums(geno_table * c(0, 0, 0, 2))
            }
            
            # Counts per marker
            if (target %in% c("both", "snp")) {
              geno_table <- gdsfmt::apply.gdsn(
                node = genotype_node,
                margin = 2,
                selection = sel,
                FUN = function(x) {
                  x <- factor(x, levels = 0:3)
                  return(as.integer(table(x)))
                },
                as.is = "list"
              )
              geno_table <- matrix(unlist(geno_table), nrow = 4)
              
              if (have_flipped) {
                flip_ref <- geno_table[1, valid_flipped]
                geno_table[1, valid_flipped] <-
                  geno_table[3, valid_flipped]
                geno_table[3, valid_flipped] <- flip_ref
              }
              
              ## Summarize genotype call counts
              object@snpAnnot$countGenoRef <- NA
              object#' @param q A numeric value [0-1] to indicate quantile to obtain.@snpAnnot$countGenoHet <- NA
              object@snpAnnot$countGenoAlt <- NA
              object@snpAnnot$countGenoMissing <- NA
              object@snpAnnot$countGenoRef[valid_markers] <-
                geno_table[3,]
              object@snpAnnot$countGenoHet[valid_markers] <-
                geno_table[2,]
              object@snpAnnot$countGenoAlt[valid_markers] <-
                geno_table[1,]
              object@snpAnnot$countGenoMissing[valid_markers] <-
                geno_table[4,]
              
              ## Summarize allele counts
              object@snpAnnot$countAlleleRef <- NA
              object@snpAnnot$countAlleleAlt <- NA
              object@snpAnnot$countAlleleMissing <- NA
              object@snpAnnot$countAlleleRef[valid_markers] <-
                colSums(geno_table * c(0, 1, 2, 0))
              object@snpAnnot$countAlleleAlt[valid_markers] <-
                colSums(geno_table * c(2, 1, 0, 0))
              object@snpAnnot$countAlleleMissing[valid_markers] <-
                colSums(geno_table * c(0, 0, 0, 2))
            }
            
            return(object)
          })

#' @rdname countRead
setMethod("countRead",
          "GbsrGenotypeData",
          function(object,
                   target = "both",
                   node = "") {
            if (node != "raw" & object@data@genotypeVar == "filt.genotype") {
              ad_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                            path = "annotation/format/AD/filt.data")
            } else {
              ad_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                            path = "annotation/format/AD/data")
            }
            valid_markers <- getValidSnp(object)
            valid_scans <- getValidScan(object)
            have_flipped <- haveFlipped(object)
            if (have_flipped) {
              valid_flipped <- getFlipped(object, valid = TRUE)
            }
            
            sel <- list(valid_scans, rep(valid_markers, each = 2))
            # Counts per sample
            if (target %in% c("both", "scan")) {
              # scan_sel <- list(rep(TRUE, nscan(object, valid=FALSE)), rep(valid_markers, each=2))
              read_count <- gdsfmt::apply.gdsn(
                node = ad_node,
                margin = 1,
                # selection=scan_sel,
                selection = sel,
                FUN = function(x) {
                  ref <- x[c(TRUE, FALSE)]
                  alt <-
                    x[c(FALSE, TRUE)]
                  if (have_flipped) {
                    flipped <- ref[valid_flipped]
                    ref[valid_flipped] <-
                      alt[valid_flipped]
                    alt[valid_flipped] <-
                      flipped
                  }
                  
                  return(c(sum(ref), sum(alt)))
                },
                as.is = "list"
              )
              read_count <- unlist(read_count)
              # read_count <- read_count[rep(valid_scans, each=2)]
              
              ## Summarize allelic read counts
              object@scanAnnot$countReadRef <- NA
              object@scanAnnot$countReadAlt <- NA
              object@scanAnnot$countReadRef[valid_scans] <-
                read_count[c(TRUE, FALSE)]
              object@scanAnnot$countReadAlt[valid_scans] <-
                read_count[c(FALSE, TRUE)]
            }
            
            
            # Counts per marker
            if (target %in% c("both", "snp")) {
              # sel_snp <- list(valid_scans, rep(valid_markers, each=2))
              read_count <- gdsfmt::apply.gdsn(
                node = ad_node,
                margin = 2,
                # selection=sel_snp,
                selection = sel,
                FUN = sum,
                as.is = "list"
              )
              read_count <- unlist(read_count)
              
              if (have_flipped) {
                valid_flipped <- rep(valid_flipped, each = 2)
                flipped <- read_count[valid_flipped]
                ref <- flipped[c(TRUE, FALSE)]
                flipped[c(TRUE, FALSE)] <- flipped[c(FALSE, TRUE)]
                flipped[c(FALSE, TRUE)] <- ref
                read_count[valid_flipped] <- flipped
              }
              
              ## Summarize allelic read counts
              object@snpAnnot$countReadRef <- NA
              object@snpAnnot$countReadAlt <- NA
              object@snpAnnot$countReadRef[valid_markers] <-
                read_count[c(TRUE, FALSE)]
              object@snpAnnot$countReadAlt[valid_markers] <-
                read_count[c(FALSE, TRUE)]
            }
            
            return(object)
          })


#' @rdname calcReadStats
#' @importFrom Biobase pData pData<-
#'
setMethod("calcReadStats",
          "GbsrGenotypeData",
          function(object, target = "both", q = NULL) {
            .gds_decomp(object)
            
            valid_markers <- getValidSnp(object)
            valid_scans <- getValidScan(object)
            have_flipped <- haveFlipped(object)
            if (have_flipped) {
              valid_flipped <- getFlipped(object, valid = TRUE)
            }
            
            sel <- list(rep(valid_markers, each = 2), valid_scans)
            
            # Calculate normarized allelic counts (counts per million).
            ad_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                          path = "annotation/format/AD")
            ad_data_node <- gdsfmt::index.gdsn(node = ad_node,
                                               path = "data")
            ls_gdsn <- gdsfmt::ls.gdsn(node = ad_data_node)
            if ("norm" %in% ls_gdsn) {
              norm_ad <- gdsfmt::index.gdsn(node = ad_data_node,
                                            path = "norm")
            } else {
              norm_ad <- gdsfmt::add.gdsn(
                node = ad_node,
                name = "norm",
                storage = "float32",
                replace = TRUE
              )
              
              gdsfmt::apply.gdsn(
                node = ad_data_node,
                margin = 1,
                FUN = function(x) {
                  denomi <- sum(x)
                  if (denomi == 0) {
                    return(x)
                  } else {
                    return(x / denomi * 10 ^ 6)
                  }
                },
                as.is = "gdsnode",
                target.node = norm_ad
              )
              gdsfmt::setdim.gdsn(node = norm_ad,
                                  valdim = gdsfmt::objdesp.gdsn(node = ad_data_node)$dim[2:1])
              gdsfmt::readmode.gdsn(node = norm_ad)
            }
            
            # read dist per sample
            if (target %in% c("both", "scan")) {
              read_dist <- gdsfmt::apply.gdsn(
                node = norm_ad,
                margin = 2,
                selection = sel,
                FUN = function(x) {
                  # x <- x[rep(valid_markers, each=2)]
                  ref <-
                    x[c(TRUE, FALSE)]
                  alt <-
                    x[c(FALSE, TRUE)]
                  if (have_flipped) {
                    flipped <- ref[valid_flipped]
                    ref[valid_flipped] <-
                      alt[valid_flipped]
                    alt[valid_flipped] <-
                      flipped
                  }
                  ref[ref == 0] <- NA
                  alt[alt == 0] <- NA
                  return(c(
                    mean(ref, na.rm = TRUE),
                    sd(ref, na.rm = TRUE),
                    quantile(ref, probs = q, na.rm =
                               TRUE),
                    mean(alt, na.rm = TRUE),
                    sd(alt, na.rm = TRUE),
                    quantile(alt, probs = q, na.rm =
                               TRUE)
                  ))
                },
                as.is = "list"
              )
              read_dist <-
                matrix(unlist(read_dist),
                       ncol = 4 + length(q) * 2,
                       byrow = TRUE)
              # read_dist <- read_dist[valid_scans, ]
              
              ## Summarize read count distribution information.
              if (is.null(q)) {
                col_names <- c("meanReadRef", "sdReadRef",
                               "meanReadAlt", "sdReadAlt")
              } else {
                col_names <-
                  c(
                    "meanReadRef",
                    "sdReadRef",
                    paste0("qtileReadRef", q),
                    "meanReadAlt",
                    "sdReadAlt",
                    paste0("qtileReadAlt", q)
                  )
              }
              pdata <- pData(object@scanAnnot)
              check <- names(pdata) %in% col_names
              if (any(check)) {
                pdata <- subset(pdata, select = !check)
              }
              df <-
                matrix(NA, nrow = nrow(pdata), ncol = ncol(read_dist))
              df[pdata$validScan,] <- read_dist
              df <- as.data.frame(df)
              names(df) <- col_names
              pData(object@scanAnnot) <- cbind(pdata, df)
            }
            
            
            # Counts per marker
            if (target %in% c("both", "snp")) {
              read_dist <- gdsfmt::apply.gdsn(
                node = norm_ad,
                margin = 1,
                selection = sel,
                FUN = function(x) {
                  # x <- x[valid_scans]
                  x[x == 0] <- NA
                  return(c(
                    mean(x, na.rm = TRUE),
                    sd(x, na.rm = TRUE),
                    quantile(x, probs = q, na.rm =
                               TRUE)
                  ))
                },
                as.is = "list"
              )
              read_dist <-
                matrix(unlist(read_dist),
                       ncol = 2 + length(q),
                       byrow = TRUE)
              
              ## Summarize allelic read counts
              if (have_flipped) {
                valid_flipped <- rep(valid_flipped, each = 2)
                flipped <- read_dist[valid_flipped,]
                ref <- flipped[c(TRUE, FALSE),]
                flipped[c(TRUE, FALSE),] <-
                  flipped[c(FALSE, TRUE),]
                flipped[c(FALSE, TRUE),] <- ref
                read_dist[valid_flipped,] <- flipped
              }
              
              pdata <- pData(object@snpAnnot)
              check <- names(pdata) %in% col_names
              if (any(check)) {
                pdata <- subset(pdata, select = !check)
              }
              read_dist <-
                cbind(read_dist[c(TRUE, FALSE),], read_dist[c(FALSE, TRUE),])
              df <-
                matrix(NA, nrow = nrow(pdata), ncol = ncol(read_dist))
              df[pdata$validMarker,] <- read_dist
              df <- as.data.frame(df)
              names(df) <- col_names
              pData(object@snpAnnot) <- cbind(pdata, df)
            }
            
            .gds_decomp(object)
            .gds_comp(object)
            return(object)
          })

#' @rdname setParents
setMethod("setParents",
          "GbsrGenotypeData",
          function(object, parents, flip, mono, bi) {
            if (length(parents) == 0 | any(is.na(parents))) {
              stop('Specify valid sample names as parents.')
            }
            if (inherits(parents, "character")) {
              id <- getScanID(object, valid = FALSE)
              p_index <- integer()
              for (i in parents) {
                p <- which(id %in% i)
                if (length(p) == 0) {
                  msg <- paste0('No sample named ', i)
                  stop(msg)
                }
                p_index <- c(p_index, p)
              }
            } else if (inherits(parents, "numeric")) {
              n_scan <- nscan(object, valid = FALSE)
              check <- n_scan < parents
              if (any(check)) {
                msg <- paste0("Total number of samples is ",
                              n_scan,
                              ". But you specified ",
                              parents[check])
                stop(msg)
              }
              p_index <- parents
            }
            
            n_parents <- length(p_index)
            object@scanAnnot$parents <- 0
            for (i in seq_len(n_parents)) {
              object@scanAnnot$parents[p_index[i]] <- i
            }
            
            object@scanAnnot$validScan[p_index] <- FALSE
            
            genotype_node <-
              gdsfmt::index.gdsn(object@data@handler, object@data@genotypeVar)
            
            geno <- integer()
            
            for (i in p_index) {
              geno <- rbind(geno,
                            gdsfmt::read.gdsn(
                              genotype_node,
                              start = c(i, 1),
                              count = c(1,-1)
                            ))
            }
            
            # Find markers which are homozygous in each parent and biallelic
            if (mono) {
              monomorphic <- apply(geno, 2, function(x) {
                return(all(x %in% c(0, 2)))
              })
            } else {
              monomorphic <- TRUE
            }
            
            if (bi) {
              biallelic <- apply(geno, 2, function(x) {
                return(length(unique(x)) != 1)
              })
            } else {
              biallelic <- TRUE
            }
            
            object <-
              setValidSnp(object, update = monomorphic & biallelic)
            
            # Find markers which p1 has an alternative allele while p2 has a reference allele.
            if (n_parents == 2 & flip) {
              message('Check flipped markers which p1 is called as alternative homozygote.')
              object@snpAnnot$flipped <- geno[1,] == 0
            }
            
            return(object)
          })

#' @rdname getParents
setMethod("getParents",
          "GbsrGenotypeData",
          function(object) {
            parents <- object@scanAnnot$parents
            if (is.null(parents)) {
              message('No information of parents.')
              return(NULL)
            }
            p_index <- which(parents != 0)
            p_id <- parents[p_index]
            p_name <- getScanID(object, valid = FALSE)[p_index]
            return(data.frame(
              scanID = p_name,
              memberID = p_id,
              indexes = p_index
            ))
          })
#'
#' setMethod("swapAlleles",
#'           "GbsrGenotypeData",
#'           function(object, allele){
#'             if(is.matrix(allele)){
#'               if(ncol(allele) != 2){
#'                 stop('The object passed to "allele" should be a data.frame with two columns for reference alleles and alternative alleles.')
#'               }
#'               if(nrow(allele) != nsnp(object)){
#'                 stop('The number of markers given as "allele" does not match with that in the data object.')
#'               }
#'               a <- getAlleleA(object)
#'               b <- getAlleleB(object)
#'               valid <- a == allele[, 1] & b == allele[, 2]
#'               flipped <- a == allele[, 2] & b == allele[, 1]
#'               invalid <- !valid & !flipped
#'
#'             } else if(is.vector(allele)){
#'               valid <- getAlleleA(object) == allele
#'               flipped <- getAlleleB(object) == allele
#'               invalid <- !valid & !flipped
#'
#'             } else {
#'               stop('Specify either of allele and genome.')
#'             }
#'             message(paste0(sum(invalid), ' SNPs were found as invalid markers which were found to have the 3rd alleles in your genotype data.'))
#'             message(paste0(sum(flipped), ' SNPs were found as flipped markers.'))
#'
#'             object <- setValidSnp(object, update = !invalid)
#'             flipped <- flipped[!invalid]
#'             # Find markers of which reference allele was called as alternative allele.
#'             object <- .setFlipped(object, flipped)
#'             return(object)
#'           })

#' @rdname thinMarker
setMethod("thinMarker",
          "GbsrGenotypeData",
          function(object, range = 150) {
            if (is.null(object@snpAnnot$countGenoMissing)) {
              stop('Run countGenotype first.')
            }
            
            missing_count <- getCountGenoMissing(object)
            
            chr <- getChromosome(object)
            pos <- getPosition(object)
            n_snp <- nsnp(object)
            valid <- rep(TRUE, n_snp)
            i <- 1
            j <- 2
            while (TRUE) {
              if (chr[i] == chr[j]) {
                mar1 <- pos[i]
                mar2 <- pos[j]
                if (mar2 - mar1 <= range) {
                  if (missing_count[i] >= missing_count[j]) {
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
              } else {
                i <- j
                j <- j + 1
              }
              if (j > n_snp) {
                break
              }
            }
            object <- setValidSnp(object, update = valid)
            return(object)
          })

#' @rdname setCallFilter
setMethod("setCallFilter",
          "GbsrGenotypeData",
          function(object,
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
                   omit_geno = NULL) {
            
            
            .gds_decomp(object)
            
            .initFilt(object)
            
            ## Quantile filtering on each genotype call
            check1 <-
              .setScanCallFilter(
                object = object,
                dp_count = dp_count,
                ref_count = ref_count,
                alt_count = alt_count,
                norm_dp_count = norm_dp_count,
                norm_ref_count = norm_ref_count,
                norm_alt_count = norm_alt_count,
                dp_qtile = scan_dp_qtile,
                ref_qtile = scan_ref_qtile,
                alt_qtile = scan_alt_qtile,
                omit_geno = omit_geno
              )
            
            check2 <- .setSnpCallFilter(
              object = object,
              dp_qtile = snp_dp_qtile,
              ref_qtile = snp_ref_qtile,
              alt_qtile = snp_alt_qtile
            )
            
            # Generate filtered AD data
            if (check1 | check2) {
              .makeCallFilteredData(object)
              object@data@genotypeVar <- "filt.genotype"
            }
            
            .gds_decomp(object)
            .gds_comp(object)
            
            return(object)
          })

.initFilt <- function(object) {
  gdsfmt::add.gdsn(object@data@handler,
                   "callfilt", 
                   val = 0,
                   valdim = c(nscan(object, valid = FALSE), nsnp(object, valid = FALSE)),
                   replace = TRUE)
}

#' @rdname setScanFilter
#' @importFrom Biobase pData pData<-
#'
setMethod("setScanFilter",
          "GbsrGenotypeData",
          function(object,
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
                   sd_alt = Inf) {
            snpID <- ploidy <- NULL
            if (missing(id)) {
              id <- rep(TRUE, nscan(object, valid = TRUE))
            } else {
              id <- getScanID(object, valid = TRUE) %in% id
            }
            
            # Filtering for samples.
            ## Missing rate
            scan_missing <-
              getCountGenoMissing(object, "scan", prop = TRUE, valid = TRUE)
            scan_missing <-
              .getSubFilter(scan_missing, missing, 1, FALSE, TRUE)
            
            ## Heterozygosity
            scan_het <-
              getCountGenoHet(object, "scan", prop = TRUE, valid = TRUE)
            scan_het <-
              .getSubFilter(scan_het, het, c(0, 1), TRUE, TRUE)
            
            ## Minor allele count
            scan_mac <- getMAC(object, "scan", valid = TRUE)
            scan_mac <- .getSubFilter(scan_mac, mac, 0, TRUE, TRUE)
            
            ## Minor allele frequency
            scan_maf <- getMAF(object, "scan", valid = TRUE)
            scan_maf <- .getSubFilter(scan_maf, maf, 0, TRUE, TRUE)
            
            ## Reference allele read count
            scan_ad_ref <-
              getCountReadRef(object, "scan", prop = FALSE, valid = TRUE)
            scan_ad_ref <-
              .getSubFilter(scan_ad_ref, ad_ref, c(0, Inf), TRUE, TRUE)
            
            
            ## Alternative allele read count
            scan_ad_alt <-
              getCountReadAlt(object, "scan", prop = FALSE, valid = TRUE)
            scan_ad_alt <-
              .getSubFilter(scan_ad_alt, ad_alt, c(0, Inf), TRUE, TRUE)
            
            ## Total read count
            scan_dp <- getCountRead(object, "scan", valid = TRUE)
            scan_dp <-
              .getSubFilter(scan_dp, dp, c(0, Inf), TRUE, TRUE)
            
            ## Mean reference allele read count
            scan_mean_ref <-
              getMeanReadRef(object, "scan", valid = TRUE)
            scan_mean_ref <-
              .getSubFilter(scan_mean_ref, mean_ref, c(0, Inf), TRUE, TRUE)
            
            ## Mean alternative allele read count
            scan_mean_alt <-
              getMeanReadAlt(object, "scan", valid = TRUE)
            scan_mean_alt <-
              .getSubFilter(scan_mean_alt, mean_alt, c(0, Inf), TRUE, TRUE)
            
            ## SD of reference allele read count
            scan_sd_ref <- getSDReadRef(object, "scan", valid = TRUE)
            scan_sd_ref <-
              .getSubFilter(scan_sd_ref, sd_ref, Inf, FALSE, TRUE)
            
            ## SD of alternative allele read count
            scan_sd_alt <- getSDReadAlt(object, "scan", valid = TRUE)
            scan_sd_alt <-
              .getSubFilter(scan_sd_alt, sd_alt, Inf, FALSE, TRUE)
            
            valid <-
              id & scan_missing & scan_het & scan_mac & scan_maf & scan_ad_ref &
              scan_ad_alt &
              scan_dp & scan_mean_ref & scan_mean_alt & scan_sd_ref & scan_sd_alt
            object <- setValidScan(object, update = valid)
            
            pData(object@snpAnnot) <-
              subset(pData(object@snpAnnot), select = snpID:ploidy)
            
            return(object)
          })

#' @rdname setSnpFilter
#' @importFrom Biobase pData pData<-
setMethod("setSnpFilter",
          "GbsrGenotypeData",
          function(object,
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
                   sd_alt = Inf) {
            scanID <- validScan <- parents <- NULL
            
            if (missing(id)) {
              id <- rep(TRUE, nsnp(object, valid = TRUE))
            } else {
              id <- !getSnpID(object, valid = TRUE) %in% id
            }
            
            # Filtering for samples.
            ## Missing rate
            snp_missing <-
              getCountGenoMissing(object, "snp", prop = TRUE, valid = TRUE)
            snp_missing <-
              .getSubFilter(snp_missing, missing, 1, FALSE, TRUE)
            
            ## Heterozygosity
            snp_het <-
              getCountGenoHet(object, "snp", prop = TRUE, valid = TRUE)
            snp_het <-
              .getSubFilter(snp_het, het, c(0, 1), TRUE, TRUE)
            
            ## Minor allele count
            snp_mac <- getMAC(object, "snp", valid = TRUE)
            snp_mac <- .getSubFilter(snp_mac, mac, 0, TRUE, TRUE)
            
            ## Minor allele frequency
            snp_maf <- getMAF(object, "snp", valid = TRUE)
            snp_maf <- .getSubFilter(snp_maf, maf, 0, TRUE, TRUE)
            
            ## Reference allele read count
            snp_ad_ref <-
              getCountReadRef(object, "snp", prop = FALSE, valid = TRUE)
            snp_ad_ref <-
              .getSubFilter(snp_ad_ref, ad_ref, c(0, Inf), TRUE, TRUE)
            
            
            ## Alternative allele read countalt_qtile = scan_alt_qtile
            snp_ad_alt <-
              getCountReadAlt(object, "snp", prop = FALSE, valid = TRUE)
            snp_ad_alt <-
              .getSubFilter(snp_ad_alt, ad_alt, c(0, Inf), TRUE, TRUE)
            
            ## Total read count
            snp_dp <- getCountRead(object, "snp", valid = TRUE)
            snp_dp <-
              .getSubFilter(snp_dp, dp, c(0, Inf), TRUE, TRUE)
            
            ## Mean reference allele read count
            snp_mean_ref <-
              getMeanReadRef(object, "snp", valid = TRUE)
            snp_mean_ref <-
              .getSubFilter(snp_mean_ref, mean_ref, c(0, Inf), TRUE, TRUE)
            
            ## Mean alternative allele read count
            snp_mean_alt <-
              getMeanReadAlt(object, "snp", valid = TRUE)
            snp_mean_alt <-
              .getSubFilter(snp_mean_alt, mean_alt, c(0, Inf), TRUE, TRUE)
            
            ## SD of reference allele read count
            snp_sd_ref <- getSDReadRef(object, "snp", valid = TRUE)
            snp_sd_ref <-
              .getSubFilter(snp_sd_ref, sd_ref, Inf, FALSE, TRUE)
            
            ## SD of alternative allele read count
            snp_sd_alt <- getSDReadAlt(object, "snp", valid = TRUE)
            snp_sd_alt <-
              .getSubFilter(snp_sd_alt, sd_alt, Inf, FALSE, TRUE)
            
            valid <-
              id & snp_missing & snp_het & snp_mac & snp_maf & snp_ad_ref &
              snp_ad_alt &
              snp_dp & snp_mean_ref & snp_mean_alt & snp_sd_ref & snp_sd_alt
            object <- setValidSnp(object, update = valid)
            pdata <- pData(object@scanAnnot)
            if ("parents" %in% names(pdata)) {
              pData(object@scanAnnot) <-
                subset(pdata, select = c(scanID, validScan, parents))
            } else {
              pData(object@scanAnnot) <-
                subset(pdata, select = c(scanID, validScan))
            }
            return(object)
          })

#' @rdname setInfoFilter
#' @importFrom Biobase pData pData<-
setMethod("setInfoFilter",
          "GbsrGenotypeData",
          function(object,
                   mq = 0 ,
                   fs = Inf,
                   qd = 0,
                   sor = Inf,
                   mqranksum = c(-Inf, Inf),
                   readposranksum = c(-Inf, Inf),
                   baseqranksum = c(-Inf, Inf)) {
            scanID <- validScan <- parents <- pdata <- NULL
            # Filtering for samples.
            ## MQ
            snp_mq <- getInfo(object, "MQ")
            snp_mq <- .getSubFilter(snp_mq, mq, 0, TRUE, TRUE)
            
            ## FS
            snp_fs <- getInfo(object, "FS")
            snp_fs <- .getSubFilter(snp_fs, fs, Inf, FALSE, TRUE)
            
            ## QD
            snp_qd <- getInfo(object, "QD")
            snp_qd <- .getSubFilter(snp_qd, qd, 0, TRUE, TRUE)
            
            ## SOR
            snp_sor <- getInfo(object, "SOR")
            snp_sor <- .getSubFilter(snp_sor, sor, Inf, FALSE, TRUE)
            
            ## MQRankSum
            snp_mqranksum <- getInfo(object, "MQRankSum")
            snp_mqranksum <-
              .getSubFilter(snp_mqranksum, mqranksum, c(-Inf, Inf), TRUE, TRUE)
            
            ## Alternative allele read count
            snp_readposranksum <- getInfo(object, "ReadPosRankSum")
            snp_readposranksum <-
              .getSubFilter(snp_readposranksum, readposranksum, c(-Inf, Inf), TRUE, TRUE)
            
            ## Total read count
            snp_baseqranksum <- getInfo(object, "BaseQRankSum")
            snp_baseqranksum <-
              .getSubFilter(snp_baseqranksum, baseqranksum, c(-Inf, Inf), TRUE, TRUE)
            
            valid <- snp_mq & snp_fs & snp_qd & snp_sor &
              snp_mqranksum & snp_readposranksum & snp_baseqranksum
            object <- setValidSnp(object, update = valid)
            return(object)
          })

#' @rdname resetScanFilters
#' @importFrom Biobase pData pData<-
setMethod("resetScanFilters",
          "GbsrGenotypeData",
          function(object) {
            scanID <- validScan <- parents <- NULL
            object <- setValidScan(object, new = TRUE)
            pdata <- pData(object@scanAnnot)
            if ("parents" %in% names(pdata)) {
              pData(object@scanAnnot) <-
                subset(pdata, select = c(scanID, validScan, parents))
            } else {
              pData(object@scanAnnot) <-
                subset(pdata, select = c(scanID, validScan))
            }
            return(object)
          })

#' @rdname resetSnpFilters
#' @importFrom Biobase pData pData<-
#'
setMethod("resetSnpFilters",
          "GbsrGenotypeData",
          function(object) {
            snpID <- ploidy <- NULL
            object <- setValidSnp(object, new = TRUE)
            pdata <- pData(object@snpAnnot)
            pdata <- subset(pdata, select = snpID:ploidy)
            pData(object@snpAnnot) <- pdata
            return(object)
          })

#' @rdname resetFilters
setMethod("resetFilters",
          "GbsrGenotypeData",
          function(object) {
            object <- setRawGenotype(object)
            object <- resetScanFilters(object)
            object <- resetSnpFilters(object)
            return(object)
          })

.getSubFilter <-
  function(variable,
           threshold,
           default,
           greater,
           equal) {
    if (is.null(variable)) {
      return(TRUE)
    }
    if (length(default) == 1) {
      if ("comp" %in% threshold[1]) {
        output <- .compareValues(variable, default, greater, FALSE)
        
      } else if (threshold != default) {
        output <- .compareValues(variable, threshold, greater, equal)
        
      } else {
        output <- TRUE
      }
    } else if (length(default) == 2) {
      if ("comp" %in% threshold[1]) {
        output <- .compareValues(variable, default[1], greater, FALSE) &
          .compareValues(variable, default[2],!greater, FALSE)
        
      } else if (any(threshold != default)) {
        output <- .compareValues(variable, threshold[1], greater, equal) &
          .compareValues(variable, threshold[2],!greater, equal)
        
      } else {
        output <- TRUE
      }
    }
    return(output)
  }

.compareValues <- function(x, y, greater, equal, na2f = TRUE) {
  if (greater & equal) {
    out <- x >= y
  } else if (greater & !equal) {
    out <- x > y
  } else if (!greater & equal) {
    out <- x <= y
  } else if (!greater & !equal) {
    out <- x < y
  }
  if (na2f) {
    out[is.na(out)] <- FALSE
  }
  return(out)
}

.setScanCallFilter <- function(object,
                               dp_count,
                               ref_count,
                               alt_count,
                               norm_dp_count,
                               norm_ref_count,
                               norm_alt_count,
                               dp_qtile,
                               ref_qtile,
                               alt_qtile,
                               omit_geno,
                               count_default = c(0, Inf),
                               qtile_default = c(0, 1)) {
  ad_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                path = "annotation/format/AD")
  ad_data_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                     path = "annotation/format/AD/data")
  callfilt <- gdsfmt::index.gdsn(object@data@handler, 
                                 "callfilt")
  
  check <- all(ref_qtile == qtile_default) &
    all(alt_qtile == qtile_default) &
    all(dp_qtile == qtile_default) &
    all(ref_count == count_default) &
    all(alt_count == count_default) &
    all(dp_count == count_default) &
    all(norm_dp_count == count_default) &
    all(norm_ref_count == count_default) &
    all(norm_alt_count == count_default) &
    is.null(omit_geno)
  
  if (!check) {
    for(i in which(getValidScan(object))){
      ref_count_i <- ref_qtile_i <-
        alt_count_i <- alt_qtile_i <-
        dp_count_i <- dp_qtile_i <-
        norm_ref_count_i <-
        norm_alt_count_i <-
        norm_dp_count_i <- omit_geno_i <- TRUE
      
      x <- gdsfmt::read.gdsn(ad_data_node, start = c(i, 1), count = c(1, -1))
      ref <- x[c(TRUE, FALSE)]
      alt <- x[c(FALSE, TRUE)]
      dp <- ref + alt
      ref[!getValidSnp(object)] <- NA
      alt[!getValidSnp(object)] <- NA
      dp[!getValidSnp(object)] <- NA
      tmp <- rep(TRUE, length(dp))
      
      if (haveFlipped(object)) {
        flipped <- getFlipped(object, valid = FALSE)
        flipped_ref <- ref[flipped]
        ref[flipped] <- alt[flipped]
        alt[flipped] <- flipped_ref
      }
      
      if (!is.null(omit_geno)) {
        if ("ref" %in% omit_geno) {
          omit_ref <- ref > 0 & alt == 0
        } else {
          omit_ref <- FALSE
        }
        if ("alt" %in% omit_geno) {
          omit_alt <- ref == 0 & alt > 0
        } else {
          omit_alt <- FALSE
        }
        if ("het" %in% omit_geno) {
          omit_het <- ref > 0 & alt > 0
        } else {
          omit_het <- FALSE
        }
        omit_geno_i <-
          !omit_ref & !omit_alt & !omit_het
      }
      
      ref_count_i <-
        .getSubFilter(ref, ref_count, count_default, TRUE, TRUE)
      alt_count_i <-
        .getSubFilter(alt, alt_count, count_default, TRUE, TRUE)
      dp_count_i <-
        .getSubFilter(dp, dp_count, count_default, TRUE, TRUE)
      
      if (!all(dp_qtile == qtile_default)) {
        threshold <- c(
          quantile(dp, probs = dp_qtile[1], na.rm = TRUE),
          quantile(dp, probs = dp_qtile[2], na.rm = TRUE)
        )
        dp_qtile_i <-
          .getSubFilter(dp, threshold, c(-1,-1), TRUE, TRUE)
      }
      
      if (!all(ref_qtile == qtile_default)) {
        threshold <- c(
          quantile(ref, probs = ref_qtile[1], na.rm = TRUE),
          quantile(ref, probs = ref_qtile[2], na.rm =
                     TRUE)
        )
        ref_qtile_i <-
          .getSubFilter(ref, threshold, c(-1,-1), TRUE, TRUE)
      }
      
      if (!all(alt_qtile == qtile_default)) {
        threshold <- c(
          quantile(alt, probs = alt_qtile[1], na.rm = TRUE),
          quantile(alt, probs = alt_qtile[2], na.rm =
                     TRUE)
        )
        alt_qtile_i <-
          .getSubFilter(alt, threshold, c(-1,-1), TRUE, TRUE)
      }
      
      if (!all(norm_ref_count == count_default)) {
        denomi <- sum(x)
        if (denomi != 0) {
          ref <- ref / denomi * 10 ^ 6
        }
        norm_ref_count_i <-
          .getSubFilter(ref, norm_ref_count, count_default, TRUE, TRUE)
      }
      
      if (!all(norm_alt_count == count_default)) {
        denomi <- sum(x)
        if (denomi != 0) {
          alt <- alt / denomi * 10 ^ 6
        }
        norm_alt_count_i <-
          .getSubFilter(alt, norm_alt_count, count_default, TRUE, TRUE)
      }
      
      if (!all(norm_dp_count == count_default)) {
        denomi <- sum(x)
        if (denomi != 0) {
          dp <- dp / denomi * 10 ^ 6
        }
        norm_dp_count_i <-
          .getSubFilter(dp, norm_dp_count, count_default, TRUE, TRUE)
      }
      
      x <- tmp & ref_count_i & ref_qtile_i &
        alt_count_i & alt_qtile_i &
        dp_count_i & dp_qtile_i &
        norm_ref_count_i &
        norm_alt_count_i &
        norm_dp_count_i & omit_geno_i
      x <- as.numeric(!x)
      
      i_callfilt <- gdsfmt::read.gdsn(callfilt, start = c(i, 1), count = c(1, -1))
      i_callfilt <- as.numeric(i_callfilt | x)
      gdsfmt::write.gdsn(callfilt,
                         val = i_callfilt,
                         start = c(i, 1), 
                         count = c(1, -1))
    }
    return(TRUE)
  } else {
    return(FALSE)
  }
}

.setSnpCallFilter <- function(object,
                              dp_qtile,
                              ref_qtile,
                              alt_qtile,
                              qtile_default = c(0, 1)) {
  ad_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                path = "annotation/format/AD")
  ad_data_node <- gdsfmt::index.gdsn(node = object@data@handler,
                                     path = "annotation/format/AD/data")
  callfilt <- gdsfmt::index.gdsn(object@data@handler, 
                                 "callfilt")
  
  n_snp <- nsnp(object, valid = FALSE)
  check <- all(dp_qtile == qtile_default) &
    all(ref_qtile == qtile_default) &
    all(alt_qtile == qtile_default)
  
  if (!check) {
    for(i in which(getValidSnp(object))){
      x <- gdsfmt::read.gdsn(ad_data_node, start = c(1, i*2-1), count = c(-1, 2))
      
      ref <- x[, c(TRUE, FALSE)]
      alt <- x[, c(FALSE, TRUE)]
      if (haveFlipped(object)) {
        flipped <- getFlipped(object, valid = FALSE)
        flipped <- rep(flipped, each = 2)
        tmp_alt <- ref[flipped]
        ref[flipped] <- alt[flipped]
        alt[flipped] <- tmp_alt
      }
      
      dp <- ref + alt
      ref[!getValidScan(object)] <- NA
      alt[!getValidScan(object)] <- NA
      dp[!getValidScan(object)] <- NA
      tmp <- rep(TRUE, length(dp))
      
      if (!all(ref_qtile == qtile_default)) {
        threshold <-
          c(quantile(ref, probs = ref_qtile[1], na.rm = TRUE),
            quantile(ref, probs = ref_qtile[2], na.rm = TRUE))
        if(all(is.na(threshold))){
          ref_qtile_i <- rep(TRUE, length(ref))
        }
        ref_qtile_i <- .getSubFilter(ref, threshold, c(-1,-1), TRUE, TRUE)
      } else {
        ref_qtile_i <- TRUE
      }
      
      if (!all(alt_qtile == qtile_default)) {
        threshold <-
          c(quantile(alt, probs = alt_qtile[1], na.rm = TRUE),
            quantile(alt, probs = alt_qtile[2], na.rm = TRUE))
        if(all(is.na(threshold))){
          alt_qtile_i <- rep(TRUE, length(alt))
        }
        alt_qtile_i <- .getSubFilter(alt, threshold, c(-1,-1), TRUE, TRUE)
      } else {
        alt_qtile_i <- TRUE
      }
      
      if (!all(dp_qtile == qtile_default)) {
        threshold <-
          c(quantile(dp, probs = dp_qtile[1], na.rm = TRUE),
            quantile(dp, probs = dp_qtile[2], na.rm = TRUE))
        if(all(is.na(threshold))){
          dp_qtile_i <- rep(TRUE, length(dp))
        }
        dp_qtile_i <- .getSubFilter(dp, threshold, c(-1,-1), TRUE, TRUE)
      } else {
        dp_qtile_i <- TRUE
      }
      x <- tmp & ref_qtile_i & alt_qtile_i & dp_qtile_i
      x <- as.numeric(!x)
      
      i_callfilt <- gdsfmt::read.gdsn(callfilt, start = c(1, i), count = c(-1, 1))
      i_callfilt <- as.numeric(i_callfilt | x)
      
      gdsfmt::write.gdsn(callfilt,
                         val = i_callfilt,
                         start = c(1, i), 
                         count = c(-1, 1))
    }
    out <- TRUE
  } else {
    out <- FALSE
  }
  return(out)
}

.makeCallFilteredData <- function(object) {
  ad_node <-
    gdsfmt::index.gdsn(node = object@data@handler, path = "annotation/format/AD")
  ad_data <-
    gdsfmt::index.gdsn(node = object@data@handler, path = "annotation/format/AD/data")
  callfilt <- gdsfmt::index.gdsn(object@data@handler, 
                                 "callfilt")
  
  ad_filt_data <- gdsfmt::add.gdsn(
    node = ad_node,
    name = "filt.data",
    storage = "int16",
    replace = TRUE
  )
  gdsfmt::assign.gdsn(ad_filt_data, src.node = ad_data)
  
  gt_data <-
    gdsfmt::index.gdsn(node = object@data@handler, path = "genotype")
  
  gt_filt_data <- gdsfmt::add.gdsn(
    node = object@data@handler,
    name = "filt.genotype",
    storage = "bit2",
    replace = TRUE
  )
  gdsfmt::assign.gdsn(gt_filt_data, src.node = gt_data)
  
  n_scan <- nscan(object, valid = FALSE)
  n_snp <- nsnp(object, valid = FALSE)
  valid_scan_indexes <- which(getValidScan(object))
  i_ad <- i_gt <- NULL
  for (i in valid_scan_indexes) {
    i_ad <- gdsfmt::read.gdsn(ad_filt_data, start = c(i, 1), count = c(1, -1))
    i_callfilt <- gdsfmt::read.gdsn(callfilt, start = c(i, 1), count = c(1, -1))
    i_ad[rep(i_callfilt, each = 2) == 1] <- 0
    gdsfmt::write.gdsn(ad_filt_data, i_ad, start = c(i, 1), count = c(1, -1))
    i_ad <- NULL
    
    i_gt <- gdsfmt::read.gdsn(gt_filt_data, start = c(i, 1), count = c(1, -1))
    i_gt[i_callfilt == 1] <- 3
    gdsfmt::write.gdsn(gt_filt_data, i_gt, start = c(i, 1), count = c(1, -1))
    i_gt <- NULL
  }
}

#' @rdname subsetGDS
setMethod("subsetGDS",
          "GbsrGenotypeData",
          function(object,
                   out_fn = "./susbet.gds",
                   snp_incl,
                   scan_incl,
                   incl_parents = TRUE) {
            n_scan <- nscan(object, valid = FALSE)
            n_snp <- nsnp(object, valid = FALSE)
            
            if (missing(snp_incl)) {
              snp_incl <- getValidSnp(object)
            } else {
              check <- n_snp == length(snp_incl)
              if (!check) {
                stop(
                  'The vector snp_incl should be the same length with the total number of SNPs in the GDS.'
                )
              }
            }
            
            if (missing(scan_incl)) {
              scan_incl <- getValidScan(object)
              if (incl_parents) {
                scan_incl[object@scanAnnot$parents != 0] <- TRUE
              }
            } else {
              check <- n_scan == length(scan_incl)
              if (!check) {
                stop(
                  'The vector scan_incl should be the same length with the total number of scans in the GDS.'
                )
              }
            }
            
            if (n_scan == sum(scan_incl) & n_snp == sum(snp_incl)) {
              message("All markers and sampels are valid.")
              message("Nothing to be subset.")
              return(NULL)
            }
            
            oldgds <- object@data@handler
            closeGDS(object)
            file.copy(object@data@filename, out_fn, overwrite = TRUE)
            newgds <- gdsfmt::openfn.gds(filename = out_fn, readonly = FALSE)
            
            .gds_decomp(newgds)
            ls_node <- gdsfmt::ls.gdsn(newgds, recursive = TRUE, include.dirs = FALSE)
            for (i_node in ls_node) {
              if(grepl("annotation/format", i_node)){
                if(!grepl("data", i_node)){
                  next
                }
              }
              newgds_i_node <- gdsfmt::index.gdsn(node = newgds, path = i_node)
              i_desc <- gdsfmt::objdesp.gdsn(newgds_i_node)
              if(any(i_desc$dim == 0)){
                next
              }
              
              if (length(i_desc$dim) == 1) {
                check <- sum(which(c(n_scan, n_snp) %in% i_desc$dim))
                if (check == 1) {
                  gdsfmt::assign.gdsn(
                    node = newgds_i_node,
                    seldim = list(which(scan_incl))
                  )
                  
                } else if (check == 2) {
                  gdsfmt::assign.gdsn(
                    node = newgds_i_node,
                    seldim = list(which(snp_incl))
                  )
                }
                
              } else if (length(i_desc$dim) == 2) {
                if (i_desc$name == "parents.genotype") {
                  times <- i_desc$dim[2] / n_snp
                  panrets_index <- seq_len(i_desc$dim[1])
                  gdsfmt::assign.gdsn(
                    node = newgds_i_node,
                    seldim = list(panrets_index,
                                  which(rep(snp_incl, each = times)))
                  )
                } else {
                  times <- i_desc$dim[2] / n_snp
                  gdsfmt::assign.gdsn(
                    node = newgds_i_node,
                    seldim = list(which(scan_incl),
                                  which(rep(snp_incl, each = times)))
                  )
                }
              }
            }
            
            .gds_comp(newgds)
            gdsfmt::closefn.gds(newgds)
            
            output <- loadGDS(gds_fn = newgds$filename)
            if (object@data@genotypeVar == "filt.genotype") {
              output <- setFiltGenotype(output)
            }
            return(output)
          })

#' @rdname setFiltGenotype
setMethod("setFiltGenotype",
          "GbsrGenotypeData",
          function(object) {
            object@data@genotypeVar <- "filt.genotype"
            return(object)
          })

#' @rdname setRawGenotype
setMethod("setRawGenotype",
          "GbsrGenotypeData",
          function(object) {
            object@data@genotypeVar <- "genotype"
            return(object)
          })


#' @rdname initScheme
setMethod("initScheme",
          "GbsrGenotypeData",
          function(object, crosstype, mating) {
            parents <- getParents(object)
            object@scheme <-
              initScheme(object@scheme, crosstype, mating, parents$memberID)
            return(object)
          })

#' @rdname addScheme
setMethod("addScheme",
          "GbsrGenotypeData",
          function(object, crosstype, mating, pop_size) {
            if (missing(pop_size)) {
              pop_size <- NA
            }
            if (missing(mating)) {
              mating <- NA
            }
            object@scheme <-
              addScheme(object@scheme, crosstype, mating, pop_size)
            return(object)
          })

#' @rdname showScheme
setMethod("showScheme",
          "GbsrGenotypeData",
          function(object) {
            parents <- getParents(object)
            showScheme(object@scheme, parents$scanID)
          })
