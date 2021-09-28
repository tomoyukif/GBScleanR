setMethod("show",
          "GbsrGenotypeData",
          function(object){
            message('Data in GDS file...')
            message('GDS file name')
            message(object@data@handler$filename)
            if(isOpenGDS(object)){
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
.insertAnnot <- function(out_gds, in_gds, out_fmt, out_info){
  info_node <- gdsfmt::index.gdsn(in_gds, "annotation/info")
  fmt_node <- gdsfmt::index.gdsn(in_gds, "annotation/format")
  new_info <- gdsfmt::index.gdsn(out_gds, "annotation/info")
  ls_gdsn <- gdsfmt::ls.gdsn(info_node)
  for(gdsn_i in ls_gdsn){
    if(!is.null(out_info)){
      if(!gdsn_i %in% out_info){
        next
      }
    }
    i_node <- gdsfmt::index.gdsn(info_node, gdsn_i)
    gdsfmt::copyto.gdsn(new_info, i_node)
  }

  n_snp <- gdsfmt::objdesp.gdsn(gdsfmt::index.gdsn(in_gds, "genotype"))$dim[2]

  ls_gdsn <- gdsfmt::ls.gdsn(fmt_node)
  new_fmt <- gdsfmt::index.gdsn(out_gds, "annotation/format")
  for(gdsn_i in ls_gdsn){
    if(!is.null(out_fmt)){
      if(!gdsn_i %in% out_fmt){
        next
      }
    }
    old_data_gdsn <- gdsfmt::index.gdsn(fmt_node, paste0(gdsn_i, "/data"))
    old_data <- gdsfmt::read.gdsn(old_data_gdsn)
    old_data[is.na(old_data)] <- "."
    n_num <- ncol(old_data) / n_snp

    if(n_num > 1){
      if(gdsn_i == "HAP"){
        sep_symbol <- "/"
      } else {
        sep_symbol <- ","
      }
      old_data <- apply(old_data, 1, function(x){
        x <- as.character(x)
        x <- matrix(x, nrow = n_num)
        return(apply(x, 2, paste, collapse = sep_symbol))
      })
      old_data <- t(old_data)
    }

    new_tmp <- gdsfmt::addfolder.gdsn(new_fmt, gdsn_i, replace = TRUE)
    new_data <- gdsfmt::add.gdsn(new_tmp,
                                 "data",
                                 old_data,
                                 storage = "string",
                                 compress = "LZMA_RA",
                                 replace = TRUE)

    old_attr <- gdsfmt::get.attr.gdsn(gdsfmt::index.gdsn(fmt_node, gdsn_i))
    if(!is.null(old_attr)){
      for(i in 1:length(old_attr))
        gdsfmt::put.attr.gdsn(new_tmp, name = names(old_attr)[i], val = old_attr[[i]])
    }
    old_attr <- gdsfmt::get.attr.gdsn(old_data_gdsn)
    if(!is.null(old_attr)){
      for(i in 1:length(old_attr))
        gdsfmt::put.attr.gdsn(new_data, name = names(old_attr)[i], val = old_attr[[i]])
    }
  }
}

.insertHaplotype <- function(out_gds, in_gds){
  new_fmt <- gdsfmt::index.gdsn(out_gds, "annotation/format")
  new_tmp <- gdsfmt::addfolder.gdsn(new_fmt, "HAP", replace = TRUE)
  gdsfmt::put.attr.gdsn(new_tmp, "Number", val = "1")
  gdsfmt::put.attr.gdsn(new_tmp, "Type", val = "Sting")
  gdsfmt::put.attr.gdsn(new_tmp, "Description", val = "Haplotype estimated by GBScleanR.")
  hap_gdsn <- gdsfmt::index.gdsn(in_gds, "estimated.haplotype")
  hap <- gdsfmt::read.gdsn(hap_gdsn)
  gdsfmt::add.gdsn(new_tmp,
                   "data",
                   hap,
                   storage = "bit6",
                   compress = "LZMA_RA",
                   replace = TRUE)

  new_info <- gdsfmt::index.gdsn(out_gds, "annotation/info")

  pgt_gdsn <- gdsfmt::index.gdsn(in_gds, "parents.genotype")
  pgt <- gdsfmt::readex.gdsn(pgt_gdsn)
  pgt <- apply(pgt, 2, paste, collapse = ",")

  new_pgt <- gdsfmt::add.gdsn(new_info,
                   "PGT",
                   pgt,
                   storage = "string",
                   compress = "LZMA_RA",
                   replace = TRUE)
  gdsfmt::put.attr.gdsn(new_pgt, "Number", val = "1")
  gdsfmt::put.attr.gdsn(new_pgt, "Type", val = "String")
  gdsfmt::put.attr.gdsn(new_pgt, "Description", val = "Genotype of each haplotype estimated by GBScleanR.")
}

.recalcDP <- function(object){
  reads <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(object@data@handler,
                                                "annotation/format/AD/data"))
  dp <- reads[, c(T, F)] + reads[, c(F, T)]
  filtdp <- gdsfmt::add.gdsn(gdsfmt::index.gdsn(object@data@handler,
                                                "annotation/format/DP"),
                             name = "filt.data",
                             val = dp,
                             storage = "vl_int",
                             compress = "LZMA_RA", replace = TRUE)
  dpdata <- gdsfmt::index.gdsn(object@data@handler,
                               "annotation/format/DP/data")
  attr <- gdsfmt::get.attr.gdsn(dpdata)
  if(!is.null(attr)){
    for(i in 1:length(attr))
      gdsfmt::put.attr.gdsn(filtdp,
                            name = names(attr)[i],
                            val = attr[[i]])
  }
  gdsfmt::delete.gdsn(dpdata)
  gdsfmt::rename.gdsn(gdsfmt::index.gdsn(object@data@handler,
                                         "annotation/format/DP/filt.data"),
                      "data")
  gdsfmt::readmode.gdsn(gdsfmt::index.gdsn(object@data@handler,
                                           "annotation/format/DP/data"))
}

#' Write a VCF file based on data in a GDS file
#'
#' Write out a VCF file with raw, filtered, or corrected genotype data
#' stored in a GDS file. The output VCF file contains only the GT filed,
#' while other annotations, AD, DP and other information will be omitted.
#'
#' @param object A GbsrGenotypeData object.
#' @param out_fn A string to specify the path to an output VCF file.
#' @param node Either one of "raw", "filt", and "cor" to output raw genotype data, filtered genotype data, or corrected genotype data, respectively.
#' @param valid A logical value to specify whether to output valid markers and samples only or all.
#' @param out_fmt A character vector to specify which variables in the annotation/format node should be output.
#' @param out_info A character vector to specify which variables in the annotation/info node should be output.
#'
#' @export
#'
#' @importFrom SeqArray seqSNP2GDS seqGDS2VCF
#'
#' @examples
#' gdata <- loadGDS("/path/to/GDS.gds")
#' gdata <- clean(gdata)
#' gbsrGDS2VCF(gdata, "/path/to/output.vcf", node = "cor")
#'
setMethod("gbsrGDS2VCF",
          "GbsrGenotypeData",
          function(object, out_fn, node, valid, out_fmt, out_info, ...){
            have_hap <- "estimated.haplotype" %in% gdsfmt::ls.gdsn(object@data@handler)
            suppressMessages(closeGDS(object))
            gds_file <- object@data@filename
            if(!grepl(".gds$", gds_file)){
              gds_file <- paste0(gds_file, ".gds")
            }
            tmp_gds <- sub(".gds", ".tmp.gds", gds_file)
            file.copy(gds_file, tmp_gds, overwrite = TRUE)
            tmp_gds <- suppressMessages(loadGDS(tmp_gds))
            tmp_gds@snpAnnot <- object@snpAnnot
            tmp_gds@scanAnnot <- object@scanAnnot

            if(haveFlipped(object)){
              tmp_gds <- flipData(tmp_gds)
            }
            if(node %in% c("filt", "cor")){
              suppressWarnings(replaceGDSdata(tmp_gds, "genotype", node))
              suppressWarnings(replaceGDSdata(tmp_gds, "ad", "filt"))
              if("DP" %in% out_fmt){
                .recalcDP(tmp_gds)
              }
            }

            if(valid == TRUE){
              subset_gds <- sub(".gds", ".tmp.subset.gds", gds_file)
              tmp_subset_gds <- suppressMessages(subsetGDS(object = tmp_gds,
                                                           out_fn = subset_gds,
                                                           incl_parents = TRUE))
              if(is.null(tmp_subset_gds)){
                gds_file <- tmp_gds@data@filename

              } else {
                suppressMessages(closeGDS(tmp_subset_gds))
                gds_file <- tmp_subset_gds@data@filename
              }
            } else {
              gds_file <- tmp_gds@data@filename
            }
            suppressMessages(closeGDS(tmp_gds))

            if(!grepl(".vcf$", out_fn)){
              out_fn <- paste0(out_fn, ".vcf")
            }
            out_fn_tmp <- sub(".vcf", ".tmp.gds", out_fn)
            SeqArray::seqSNP2GDS(gds.fn = gds_file,
                                 out.fn = out_fn_tmp,
                                 verbose = FALSE)

            out_gds <- gdsfmt::openfn.gds(out_fn_tmp, readonly = FALSE)
            if(have_hap){
              in_gds <- gdsfmt::openfn.gds(object@data@filename)
              .insertHaplotype(out_gds, in_gds)
              gdsfmt::closefn.gds(in_gds)
            }

            in_gds <- gdsfmt::openfn.gds(gds_file)
            .insertAnnot(out_gds, in_gds, out_fmt, out_info)
            gdsfmt::closefn.gds(in_gds)
            gdsfmt::closefn.gds(out_gds)

            SeqArray::seqGDS2VCF(gdsfile=out_fn_tmp,
                                 vcf.fn=out_fn,
                                 verbose = FALSE)
            on.exit({
              if(valid){
                suppressWarnings(file.remove(tmp_subset_gds@data@filename))
              }
              suppressWarnings(file.remove(tmp_gds@data@filename))
              suppressWarnings(file.remove(out_fn_tmp))
            })
            object@data@handler <- gdsfmt::openfn.gds(object@data@filename,
                                                      readonly = FALSE,
                                                      allow.duplicate = TRUE)
            return(object)
          })

#' Check if a GDS file has been opened or not.
#'
#' @param object A GbsrGenotypeData object.
#'
#' @return
#' `TRUE` if the GDS file linked to the input GbsrGenotypeData object has been opened, while `FALSE` if closed.
#'
#'@export
#'
setMethod("isOpenGDS",
          "GbsrGenotypeData",
          function(object, ...){
            tryout <- try(gdsfmt::openfn.gds(filename = object@data@filename), silent=TRUE)
            if(inherits(tryout, "gds.class")){
              gdsfmt::closefn.gds(tryout)
              return(FALSE)

            } else if(grepl("has been created or opened.", tryout[1])){
              return(TRUE)

            } else {
              warning(tryout[1])
              return(NULL)
            }
          })

#' Close the connection to the GDS file
#'
#' Close the connection to the GDS file linked to the given GbsrGenotypeData object.
#'
#' @param object A GbsrGenotypeData object.
#'
#'@export
#'
setMethod("closeGDS",
          "GbsrGenotypeData",
          function(object, ...){
            gdsfmt::closefn.gds(object@data@handler)
            message('The connection to the GDS file was closed.')
          })

#' Write out the information stored in the SnpAnnotationDataSet slot
#'
#' All the data stored in the SnpAnnotatoinDataSet slot of the GbsrGenotypeData
#' object can be saved in the GDS file linked to the given GbsrGenotypeData object.
#' You can load the saved data using [loadSnpAnnot()].
#'
#' @param object A GbsrGenotypeData object
#'
#' @export
#'
setMethod("saveSnpAnnot",
          "GbsrGenotypeData",
          function(object, ...){
            new_node <- gdsfmt::add.gdsn(node=object@data@handler,
                                         name="snpAnnot",
                                         val=pData(object@snpAnnot),
                                         compress="LZMA_RA",
                                         replace=TRUE)
            gdsfmt::readmode.gdsn(node=new_node)
          })

#' Write out the information stored in the ScanAnnotationDataSet slot
#'
#' All the data stored in the ScanAnnotationDataSet slot of the GbsrGenotypeData
#' object can be saved in the GDS file linked to the given GbsrGenotypeData object.
#' You can load the saved data using [loadSnpAnnot()].
#'
#' @param object A GbsrGenotypeData object
#'
#'@export
#'
setMethod("saveScanAnnot",
          "GbsrGenotypeData",
          function(object, ...){
            new_node <- gdsfmt::add.gdsn(node=object@data@handler,
                                         name="scanAnnot",
                                         val=pData(object@scanAnnot),
                                         compress="LZMA_RA",
                                         replace=TRUE)
            gdsfmt::readmode.gdsn(node=new_node)
          })

#' Load the stored SnpAnnotationDataSet information
#'
#' All the data stored in the SnpAnnotatoinDataSet slot of the GbsrGenotypeData
#' object can be saved in the GDS file linked to the given GbsrGenotypeData object via [saveSnpAnnot()].
#' You can load the saved data using this function.
#'
#' @param object A GbsrGenotypeData object
#'
#' @export
#'
setMethod("loadSnpAnnot",
          "GbsrGenotypeData",
          function(object, ...){
            ls_gdsn <- snpAnnot_node <- gdsfmt::ls.gdsn(node=object@data@handler)
            if("snpAnnot" %in% ls_gdsn){
              snpAnnot_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                                  path="snpAnnot")
              ls_gdsn <- gdsfmt::ls.gdsn(node=snpAnnot_node)
              for(i in ls_gdsn){
                gdsfmt::readmode.gdsn(node=index.gdsn(node=snpAnnot_node, path=i))
              }
              snpAnnot <- gdsfmt::read.gdsn(node=snpAnnot_node)
              pData(object@snpAnnot) <- snpAnnot
            } else {
              message('No data of snpAnnot in the GDS file.')
            }
            return(object)
          })

#' Load the stored ScanAnnotationDataSet information
#'
#' All the data stored in the ScanAnnotationDataSet slot of the GbsrGenotypeData
#' object can be saved in the GDS file linked to the given GbsrGenotypeData object via [saveScanAnnot()].
#' You can load the saved data using this function.
#'
#' @param object A GbsrGenotypeData object
#'
#' @export
#'
setMethod("loadScanAnnot",
          "GbsrGenotypeData",
          function(object, ...){
            ls_gdsn <- snpAnnot_node <- gdsfmt::ls.gdsn(node=object@data@handler)
            if("scanAnnot" %in% ls_gdsn){
              scanAnnot_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                                   path="scanAnnot")
              ls_gdsn <- gdsfmt::ls.gdsn(node=scanAnnot_node)
              for(i in ls_gdsn){
                gdsfmt::readmode.gdsn(node=index.gdsn(node=scanAnnot_node, path=i))
              }
              scanAnnot <- gdsfmt::read.gdsn(node=scanAnnot_node)
              pData(object@scanAnnot) <- scanAnnot
            } else {
              message('No data of snpAnnot in the GDS file.')
            }
            return(object)
          })

#' Return the number of SNPs.
#'
#' This function returns the number of SNPs recorded in the GDS file
#' connected to the given GbsrGenotypeData object.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#'
#' @details
#' If `valid = TRUE`, the number of SNPs which are labeled `TRUE` in
#' the SnpAnnotationDataSet slot will be returned. You need the number
#' of over all SNPs, set `valid = FALSE`. [getValidSnp()] tells you
#' which markers are valid.
#'
#' @seealso [getValidSnp()]
#'
#' @export
#'
setMethod("nsnp",
          "GbsrGenotypeData",
          function(object, valid=TRUE, ...){
            if(valid){
              out <- sum(getValidSnp(object))
            } else {
              out <- nrow(object@snpAnnot)
            }
            return(out)
          })

#' Return the number of scans (samples).
#'
#' This function returns the number of samples recorded in the GDS file
#' connected to the given GbsrGenotypeData object.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#'
#' @details
#' If `valid = TRUE`, the number of samples which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. You need
#' the number of over all samples, set `valid = FALSE`.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @seealso [getValidSnp()]
#'
#'@export
#'
setMethod("nscan",
          "GbsrGenotypeData",
          function(object, valid=TRUE, ...){
            if(valid){
              out <- sum(getValidScan(object))
            } else {
              out <- nrow(object@scanAnnot)
            }
            return(out)
          })

#' Return a logical vector indicating which are valid SNP markers.
#'
#' @param object A GbsrGenotypeData object.
#'
#' @seealso [setValidSnp()]
#'
#' @export
#'
setMethod("getValidSnp",
          "GbsrGenotypeData",
          function(object, ...){
            return(object@snpAnnot$validMarker)
          })

#' Manually set valid SNP markers.
#'
#' If you need manually set valid and invalid SNP markers, you can do it via this function,
#' e.g in the case you conducted a filtering on SNP markers manually by your self.
#'
#' @param object A GbsrGenotypeData object.
#' @param new A logical vector of the same length with the over all number of the SNP markers.
#' @param update A logical vector of the same length with the currently valid SNP markers.
#'
#' @details
#' To over write the current validity information, give a logical vector to `new`.
#' On the other hand, a logical vector specified to `update` will be used to
#' update validity information of the currently valid SNP markers. If you gave
#' a vector for both argument, only the vector passed to `new` will be used to
#' over write the validity information.
#'
#' @seealso [setSnpFilter()] to filter out SNP markers based on some summary statistics.
#'
#' @export
#'
setMethod("setValidSnp",
          "GbsrGenotypeData",
          function(object, new, update, ...){
            if(!missing(new)){
              if(any(is.na(new))){
                stop('NA is not allowed for a logical vector "new".')
              }
              object@snpAnnot$validMarker <- new
            } else if(!missing(update)){
              if(any(is.na(update))){
                stop('NA is not allowed for a logical vector "update".')
              }
              update <- object@snpAnnot$validMarker[getValidSnp(object)] & update
              object@snpAnnot$validMarker[getValidSnp(object)] <- update
            }
            return(object)
          })

#' Return a logical vector indicating which are valid scans (samples).
#'
#' @param object A GbsrGenotypeData object.
#'
#' @seealso [setValidScan()]
#'
#' @export
#'
setMethod("getValidScan",
          "GbsrGenotypeData",
          function(object, parents = FALSE, ...){
            if(parents == "only"){
              return(object@scanAnnot$parents != 0)
            }
            if(parents){
              out <- object@scanAnnot$validScan
              out[object@scanAnnot$parents != 0] <- TRUE
              return(out)
            } else {
              return(object@scanAnnot$validScan)
            }
          })

#' Manually set valid scans (samples).
#'
#' If you need manually set valid and invalid samples, you can do it via this function,
#' e.g in the case you conducted a filtering on samples manually by your self.
#'
#' @param object A GbsrGenotypeData object.
#' @param new A logical vector of the same length with the over all number of the samples.
#' @param update A logical vector of the same length with the currently valid samples.
#'
#' @details
#' To over write the current validity information, give a logical vector to `new`.
#' On the other hand, a logical vector specified to `update` will be used to
#' update validity information of the currently valid samples. If you gave
#' a vector for both argument, only the vector passed to `new` will be used to
#' over write the validity information.
#'
#' @seealso [setScanFilter()] to filter out samples based on some summary statistics
#'
#' @export
#'
setMethod("setValidScan",
          "GbsrGenotypeData",
          function(object, new, update, ...){
            if(!missing(update)){
              if(any(is.na(update))){
                stop('NA is not allowed for a logical vector "update".')
              }
              update <- object@scanAnnot$validScan[getValidScan(object)] & update
              object@scanAnnot$validScan[getValidScan(object)] <- update
            } else if(!missing(new)){
              if(any(is.na(new))){
                stop('NA is not allowed for a logical vector "new".')
              }
              object@scanAnnot$validScan <- new
            }
            return(object)
          })

# This method is internally used.
setMethod("setFlipped",
          "GbsrGenotypeData",
          function(object, flipped, ...){
            if(any(is.na(flipped))){
              stop('NA is not allowed for a logical vector "flipped".')
            }
            valid_snp <- getValidSnp(object)
            if(length(valid_snp) == length(flipped)){
              object@snpAnnot$flipped <- flipped
            } else {
              if(nsnp(object) == length(flipped)){
                tmp_flipped <- rep(FALSE, nsnp(object, valid=FALSE))
                tmp_flipped[valid_snp] <- flipped
                object@snpAnnot$flipped <- tmp_flipped
              } else {
                stop('The length of "flipped" does not match with the number of SNPs.')
              }
            }
            return(object)
          })

# Internally used function to flip genotype data based on the flipped marker information.
# Flipped markers are markers where the alleles expected as reference allele are called as
# alternative allele.
.flipGeno <- function(object, var){
  flipped <- getFlipped(object, valid = FALSE)
  gt <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(object@data@handler, var))
  flip_gt <- gt[, flipped]
  flip_0 <- flip_gt == 0
  flip_2 <- flip_gt == 2
  flip_gt[flip_0] <- 2
  flip_gt[flip_2] <- 0
  gt[, flipped] <- flip_gt
  gt_attr <- gdsfmt::get.attr.gdsn(gdsfmt::index.gdsn(object@data@handler, var))
  gt_gdsn <- gdsfmt::add.gdsn(object@data@handler, var, gt, "bit2",
                              compress = "LZMA_RA",
                              replace = TRUE)
  if(!is.null(gt_attr)){
    for(i in 1:length(gt_attr)){
      gdsfmt::put.attr.gdsn(node=gt_gdsn,
                            name=names(gt_attr)[i],
                            val=gt_attr[[i]])
    }
  }
  gdsfmt::readmode.gdsn(gdsfmt::index.gdsn(object@data@handler, var))
}

# Internally used function to flip AD data based on the flipped marker information.
# Flipped markers are markers where the alleles expected as reference allele are called as
# alternative allele.
.flipAD <- function(object, var){
  flipped <- getFlipped(object, valid = FALSE)
  ad <- gdsfmt::index.gdsn(object@data@handler, "annotation/format/AD")
  data <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(ad, var))
  ref <- data[, c(TRUE, FALSE)]
  alt <- data[, c(FALSE, TRUE)]
  tmp <- ref[, flipped]
  ref[, flipped] <- alt[, flipped]
  alt[, flipped] <- tmp
  data[, c(TRUE, FALSE)] <- ref
  data[, c(FALSE, TRUE)] <- alt
  gdsfmt::add.gdsn(node = ad,
                   name = var,
                   val = data,
                   storage = "vl_int",
                   compress = "LZMA_RA",
                   replace = TRUE)
  gdsfmt::readmode.gdsn(gdsfmt::index.gdsn(ad, var))
}

#' Flip genotype, allele information, and allele read counts of the flipped SNP markers.
#'
#' Genotype data, allele information, and allele read counts to be the expected reference
#' allele as actually reference allele.
#'
#' @param object A GbsrGenotypeData object.
#'
#' @details
#' Flipped markers are markers where the alleles expected as reference allele are called as
#' alternative allele. If you specify two parents in the `parents` argument of
#' [setParents()] with `flip = TRUE`, `bi = TRUE`, and `homo = TRUE`, the alleles found
#' in the parent specified as the first element to the `parents` argument are supposed as
#' reference alleles of the markers. If the "expected" reference alleles are not actually
#' called as reference alleles but alternative alleles in the given data. [setParents()] will
#' automatically labels those markers "flipped". The SnpAnnotatoinDataSet slot sores this
#' information and accessible via [getFlipped()] which gives you a logical vector
#' indicating which markers are labeled as flipped `TRUE` or not flipped `FALSE`.
#' [haveFlipped()] just tells you whether the SnpAnnotatoinDataSet slot has
#' the information of flipped markers or not.
#'
#' @seealso [setParents()], [getFlipped()], and [haveFlipped()].
#'
#' @export
#'
setMethod("flipData",
          "GbsrGenotypeData",
          function(object, ...){
            if(!haveFlipped(object)){
              stop('Nothing to flip.')
            }

            # Genotypes
            .flipGeno(object, "genotype")
            ls_gdsn <- gdsfmt::ls.gdsn(object@data@handler)
            if("filt.genotype" %in% ls_gdsn){
              .flipGeno(object, "filt.genotype")
            }

            # Alleles
            allele <- paste(getAlleleA(object, valid = FALSE),
                            getAlleleB(object, valid = FALSE), sep = "/")
            gdsfmt::add.gdsn(node = object@data@handler,
                             name = "snp.allele",
                             val = allele,
                             storage = "string",
                             compress = "LZMA_RA",
                             replace = TRUE)
            gdsfmt::readmode.gdsn(gdsfmt::index.gdsn(object@data@handler, "snp.allele"))

            # Allele reads
            .flipAD(object, "data")

            # Filtered allele reads
            ls_gdsn <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(object@data@handler,
                                                          "annotation/format/AD"))
            if("filt.data" %in% ls_gdsn){
              .flipAD(object, "filt.data")
            }

            object@snpAnnot$flipped <- NULL
            return(object)
          })


#' Get a logical vector indicating flipped SNP markers.
#'
#' @param object A GbsrGenotypeData object.
#'
#' @details
#' Flipped markers are markers where the alleles expected as reference allele are called as
#' alternative allele. If you specify two parents in the `parents` argument of
#' [setParents()] with `flip = TRUE`, `bi = TRUE`, and `homo = TRUE`, the alleles found
#' in the parent specified as the first element to the `parents` argument are supposed as
#' reference alleles of the markers. If the "expected" reference alleles are not actually
#' called as reference alleles but alternative alleles in the given data. [setParents()] will
#' automatically labels those markers "flipped". The SnpAnnotatoinDataSet slot sores this
#' information and accessible via [getFlipped()] which gives you a logical vector
#' indicating which markers are labeled as flipped `TRUE` or not flipped `FALSE`.
#' [haveFlipped()] just tells you whether the SnpAnnotatoinDataSet slot has
#' the information of flipped markers or not.
#'
#' @seealso [setParents()] and [haveFlipped()].
#'
#' @export
#'
setMethod("getFlipped",
          "GbsrGenotypeData",
          function(object, valid=TRUE, ...){
            out <- object@snpAnnot$flipped
            if(is.null(out)){
              message('No data of flipped genotype markers.')
              return(NULL)
            }
            if(valid){
              out <- out[getValidSnp(object)]
            }
            return(out)
          })


#' Get a logical value indicating flipped SNP markers whether information exists.
#'
#' @param object A GbsrGenotypeData object.
#'
#' @details
#' Flipped markers are markers where the alleles expected as reference allele are called as
#' alternative allele. If you specify two parents in the `parents` argument of
#' [setParents()] with `flip = TRUE`, `bi = TRUE`, and `homo = TRUE`, the alleles found
#' in the parent specified as the first element to the `parents` argument are supposed as
#' reference alleles of the markers. If the "expected" reference alleles are not actually
#' called as reference alleles but alternative alleles in the given data. [setParents()] will
#' automatically labels those markers "flipped". The SnpAnnotatoinDataSet slot sores this
#' information and accessible via [getFlipped()] which gives you a logical vector
#' indicating which markers are labeled as flipped `TRUE` or not flipped `FALSE`.
#' [haveFlipped()] just tells you whether the SnpAnnotatoinDataSet slot has
#' the information of flipped markers or not.
#'
#' @seealso [setParents()] and [getFlipped()].
#'
#' @export
#'
setMethod("haveFlipped",
          "GbsrGenotypeData",
          function(object, ...){
            return(!is.null(object@snpAnnot$flipped))
          })

#' Get read count data.
#'
#' Read counts for reference allele and alternative allele are retrieved
#' from the GDS file linked to the given GbsrGenotypeData object.
#'
#' @param object A gbsrGenotypeData object.
#' @param chr A integer vector of indexes indicating chromosomes to get read count data.
#' @param node Either of "raw" and "filt". See details.
#' @param parents A logical value or "only" to include data for parents or to get data only for parents.
#'
#' @details
#' Read count data can be obtained from the "annotation/format/AD/data" node or the
#' "annotation/format/AD/filt.data" node of the GDS file with `node = "raw"` or
#' `node = "filt"`, respectively. The [setCallFilter()] function generate filtered
#' read count data in the "annotation/format/AD/filt.data" node which can be accessed as
#' mentioned above.
#'
#' @seealso [setCallFilter()]
#'
#' @export
#'
setMethod("getRead",
          "GbsrGenotypeData",
          function(object, chr, node, parents, ...){
            ad_gdsn <- gdsfmt::index.gdsn(object@data@handler, "annotation/format/AD")
            ls_gdsn <- gdsfmt::ls.gdsn(ad_gdsn)
            if(node == "filt" & "filt.data" %in% ls_gdsn){
              path <- "filt.data"
            } else {
              path <- "data"
            }
            ad_node <- gdsfmt::index.gdsn(ad_gdsn, path)
            valid_samples <- getValidScan(object)
            if(parents == "only"){
              valid_samples <- object@scanAnnot$parents != 0
            } else if(parents){
              valid_samples[object@scanAnnot$parents != 0] <- TRUE
            }

            valid_snp <- getValidSnp(object)
            if(!is.null(chr)){
              valid_snp <- valid_snp & getChromosome(object, valid = FALSE) == chr
            }
            ad <- gdsfmt::readex.gdsn(node=ad_node,
                                      sel=list(valid_samples,
                                               rep(valid_snp, each = 2)))
            ref <- ad[, c(T, F)]
            alt <- ad[, c(F, T)]

            if(haveFlipped(object)){
              flipped <- getFlipped(object, valid=FALSE)[valid_snp]
              flipped_ref <- ref[, flipped]
              ref[, flipped] <- alt[, flipped]
              alt[, flipped] <- flipped_ref
            }

            rownames(ref) <- getScanID(object, valid = FALSE)[valid_samples]
            rownames(alt) <- getScanID(object, valid = FALSE)[valid_samples]
            return(list(ref = ref, alt = alt))
          })


#' Get genotype call data.
#'
#' Genotype calls are retrieved from the GDS file linked to the given
#' GbsrGenotypeData object.
#'
#' @param object A gbsrGenotypeData object.
#' @param chr A integer vector of indexes indicating chromosomes to get read count data.
#' @param node Either of "raw", "filt", and "cor. See details.
#' @param parents A logical value or "only" to include data for parents or to get data only for parents.
#'
#' @details
#' Genotype call data can be obtained from the "genotype" node, the "filt.genotype"
#' node, or the "corrected.genotype" node of the GDS file with `node = "raw"`,
#' `node = "filt"`, or `node = "raw"`, respectively. If `node = "parents`, the data in the "parents.genotype" node will be returned. The "parents.genotype" node stores phased parental genotypes estimated by the [clean()] function.
#' The [setCallFilter()] function generate filtered genotype call data in the
#' "filt.genotype" node which can be accessed as mentioned above. On the other hand, the
#' "corrected.genotype" node can be generated via the [clean()] function.
#'
#' @seealso [setCallFilter()] and [clean()]
#'
#' @export
#'
setMethod("getGenotype",
          "GbsrGenotypeData",
          function(object, chr, node, parents, ...){
            ls_gdsn <- gdsfmt::ls.gdsn(object@data@handler)
            path <- NULL
            if(node == "parents"){
             if("parents.genotype" %in% ls_gdsn){
               path <- "parents.genotype"
             } else {
               stop("No estimated parents genotype data.")
             }
            }
            if(node == "cor"){
              if("corrected.genotype" %in% ls_gdsn){
                path <- "corrected.genotype"
              } else {
                stop("No corrected genotype data.")
              }
            }
            if(node == "filt"){
              if("filt.genotype" %in% ls_gdsn){
                path <- "filt.genotype"
              } else {
                stop("No filtered genotype data.")
              }
            }
            if(node == "raw"){
              path <- "genotype"
            }
            if(is.null(path)){
              stop("Please specify a valid node.")
            }
            genotype_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                                path=path)
            valid_markers <- getValidSnp(object)
            if(!is.null(chr)){
              valid_marker[getChromosome(object) != chr] <- FALSE
            }

            if(node == "parents"){
              n_row <- gdsfmt::objdesp.gdsn(genotype_node)$dim[1]
              sel <- list(rep(TRUE, n_row),
                          valid_markers)
              genotype <- gdsfmt::readex.gdsn(node=genotype_node,
                                              sel=sel)
            } else {
              valid_samples <- getValidScan(object)
              if(parents == "only"){
                valid_samples <- object@scanAnnot$parents != 0
              } else if(parents){
                valid_samples[object@scanAnnot$parents != 0] <- TRUE
              }
              sel <- list(valid_samples, valid_markers)
              genotype <- gdsfmt::readex.gdsn(node=genotype_node,
                                              sel=sel)
              rownames(genotype) <- getScanID(object,
                                              valid = FALSE)[valid_samples]
            }
            genotype[genotype == 3] <- NA
            return(genotype)
          })


#' Get haplotype call data.
#'
#' Haplotype calls are retrieved from the GDS file linked to the given
#' GbsrGenotypeData object.
#'
#' @param object A gbsrGenotypeData object.
#' @param chr A integer vector of indexes indicating chromosomes to get read count data.
#' @param parents A logical value or "only" to include data for parents or to get data only for parents.
#'
#' @details
#' Haplotype call data can be obtained from the "estimated.haplotype" node of
#' the GDS file which can be generated via the [clean()] function. Thus, this function
#' is valid only after having executed [clean()].
#'
#' @seealso [clean()]
#'
#' @export
#'
setMethod("getHaplotype",
          "GbsrGenotypeData",
          function(object, chr, parents, ...){
            ls_gdsn <- gdsfmt::ls.gdsn(object@data@handler)
            if("estimated.haplotype" %in% ls_gdsn){
              path <- "estimated.haplotype"
            } else {
              stop('No haplotype data, Run clean() to estimate haplotype.')
            }
            haplotype_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                                 path=path)
            valid_samples <- getValidScan(object)
            if(parents == "only"){
              valid_samples <- object@scanAnnot$parents != 0
            } else if(parents){
              valid_samples[object@scanAnnot$parents != 0] <- TRUE
            }
            valid_markers <- getValidSnp(object)
            sel <- list(valid_samples, rep(valid_markers, each = 2))
            haplotype <- gdsfmt::readex.gdsn(node=haplotype_node,
                                             sel=sel)
            haplotype[haplotype == 0] <- NA
            haplotype <- array(t(haplotype), c(2, sum(valid_markers), sum(valid_samples)))
            if(!is.null(chr)){
              haplotype <- haplotype[, getChromosome(object) == chr, ]
            }
            return(haplotype)
          })

#' Obtain chromosome information of each SNP marker
#'
#' This function returns indexes or names of chromsomes of each SNP or just a
#' set of unique chromosome names.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#' @param levels A logical value. See details.
#' @param name A logical value. See details.
#'
#' @details
#' A GDS file created via GBScleanR stores chromosome names as sequential integers
#' from 1 to N, where N is the number of chromosomes. This function returns those
#' indexes as default. If you need actual names of the chromosomes, set `name = TRUE`.
#' `levels = TRUE` gives you only unique chromosome names with length N.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
#'
setMethod("getChromosome",
          "GbsrGenotypeData",
          function(object, valid=TRUE, levels=FALSE, name=FALSE, ...){
            if(name){
              out <- object@snpAnnot$chromosome.name
            } else {
              out <- object@snpAnnot$chromosome
            }
            if(valid){
              out <- out[getValidSnp(object)]
            }
            if(levels){
              out <- unique(out)
            }
            return(out)
          })


#' Obtain physical position information of each SNP marker
#'
#' This function returns physical positions of SNP markers.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getPosition",
          "GbsrGenotypeData",
          function(object, valid=TRUE, ...){
            out <- object@snpAnnot$position
            if(valid){
              out <- out[getValidSnp(object)]
            }
            return(out)
          })

#' Obtain reference allele information of each SNP marker
#'
#' This function returns reference alleles, either of A, T, G, and C, of SNP markers.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getAlleleA",
          "GbsrGenotypeData",
          function(object, valid=TRUE, ...){
            out <- object@snpAnnot$alleleA
            if(haveFlipped(object)){
              flipped <- getFlipped(object, valid = FALSE)
              b <- object@snpAnnot$alleleB
              out[flipped] <- b[flipped]
            }
            if(valid){
              out <- out[getValidSnp(object)]
            }
            return(out)
          })

#' Obtain alternative allele information of each SNP marker
#'
#' This function returns alternative alleles, either of A, T, G, and C, of SNP markers.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getAlleleB",
          "GbsrGenotypeData",
          function(object, valid=TRUE, ...){
            out <- object@snpAnnot$alleleB
            if(haveFlipped(object)){
              flipped <- getFlipped(object, valid = FALSE)
              a <- object@snpAnnot$alleleA
              out[flipped] <- a[flipped]
            }
            if(valid){
              out <- out[getValidSnp(object)]
            }
            return(out)
          })

#' Obtain SNP ID
#'
#' This function returns SNP ID of SNP markers.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getSnpID",
          "GbsrGenotypeData",
          function(object, valid=TRUE, ...){
            out <- object@snpAnnot$snpID
            if(valid){
              out <- out[getValidSnp(object)]
            }
            return(out)
          })

#' Obtain scan (sample) ID
#'
#' This function returns scan (sample) ID.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getScanID",
          "GbsrGenotypeData",
          function(object, valid=TRUE, ...){
            out <- object@scanAnnot$scanID
            if(valid){
              out <- out[getValidScan(object)]
            }
            return(out)
          })

#' Obtain ploidy information of each SNP marker
#'
#' This function returns ploidy of each SNP marker. The ploidy of all the markers in
#' a dataset is a same value and the current implementation of GBScleanR only works
#' with data having ploidy = 2 for all markers.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getPloidy",
          "GbsrGenotypeData",
          function(object, valid=TRUE, ...){
            out <- object@snpAnnot$ploidy
            if(valid){
              out <- out[getValidSnp(object)]
            }
            return(out)
          })

#' Obtain information stored in the "annotation/info" node
#'
#' The "annotation/info" node stores annotation infromation of markers obtained
#' via SNP calling tools like bcftools and GATK.
#'
#' @param object A GbsrGenotypeData object.
#' @param var A string to indicate which annotation info should be retrieved.
#' @param valid A logical value. See details.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getInfo",
          "GbsrGenotypeData",
          function(object, var, valid=TRUE, ...){
            path <- paste0("annotation/info/", var)
            ls_gdsn <- snpAnnot_node <- gdsfmt::ls.gdsn(node=object@data@handler,
                                                        recursive=TRUE,
                                                        include.dirs=TRUE)
            if(path %in% ls_gdsn){
              info_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                              path=path)
            } else {
              message(paste0('No data of ', var, '.'))
              return(NULL)
            }

            if(valid){
              out <- gdsfmt::readex.gdsn(node=info_node, sel=list(getValidSnp(object)))
            } else {
              out <- gdsfmt::read.gdsn(node=info_node)
            }
            return(out)
          })

#' Obtain total reference read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param prop A logical value whether to return values as proportions of total reference read counts in total read counts per SNP or not.
#'
#' @details
#' You need to execute [countRead()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getCountReadRef",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, prop=FALSE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$countReadRef
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
              if(prop){
                out <- out / getCountRead(object, target = "snp", valid=valid)
              }
            } else {
              out <- object@scanAnnot$countReadRef
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
              if(prop){
                out <- out / getCountRead(object, target = "scan", valid=valid)
              }
            }
            return(out)
          })

#' Obtain total alternative read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param prop A logical value whether to return values as proportions of total alternative read counts in total read counts per SNP or not.
#'
#' @details
#' You need to execute [countRead()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getCountReadAlt",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, prop=FALSE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$countReadAlt
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
              if(prop){
                out <- out / getCountRead(object, target = "snp", valid=valid)
              }
            } else {
              out <- object@scanAnnot$countReadAlt
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
              if(prop){
                out <- out / getCountRead(object, target = "scan", valid=valid)
              }
            }
            return(out)
          })

#' Obtain total read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [countRead()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getCountRead",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, ...){
            if(target == "snp"){
              out1 <- getCountReadRef(object, "snp", prop=FALSE, valid=valid)
              out2 <- getCountReadAlt(object, "snp", prop=FALSE, valid=valid)
              if(is.null(out1) | is.null(out2)) {return(NULL)}

            } else {
              out1 <- getCountReadRef(object, "scan", prop=FALSE, valid=valid)
              out2 <- getCountReadAlt(object, "scan", prop=FALSE, valid=valid)
              if(is.null(out1) | is.null(out2)) {return(NULL)}
            }
            out <- out1 + out2
            return(out)
          })

#' Obtain total reference genotype counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total reference genotype counts to total non missing genotype counts or not.
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getCountGenoRef",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, prop=FALSE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$countGenoRef
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
              if(prop){
                out <- out / {nscan(object, valid=valid) -
                    getCountGenoMissing(object, "snp", valid=valid)}
              }
            } else {
              out <- object@scanAnnot$countGenoRef
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
              if(prop){
                out <- out / {nsnp(object, valid=valid) -
                    getCountGenoMissing(object, "scan", valid=valid)}
              }
            }
            return(out)
          })

#' Obtain total heterozygote counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total heterozygote counts to total non missing genotype counts or not.
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getCountGenoHet",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, prop=FALSE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$countGenoHet
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
              if(prop){
                out <- out / {nscan(object, valid=valid) -
                    getCountGenoMissing(object, "snp", valid=valid)}
              }
            } else {
              out <- object@scanAnnot$countGenoHet
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
              if(prop){
                out <- out / {nsnp(object, valid=valid) -
                    getCountGenoMissing(object, "scan", valid=valid)}
              }
            }
            return(out)
          })

#' Obtain total alternative genotype counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total alternative genotype counts to total non missing genotype counts or not.
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getCountGenoAlt",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, prop=FALSE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$countGenoAlt
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
              if(prop){
                out <- out / {nscan(object, valid=valid) -
                    getCountGenoMissing(object, "snp", valid=valid)}
              }
            } else {
              out <- object@scanAnnot$countGenoAlt
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
              if(prop){
                out <- out / {nsnp(object, valid=valid) -
                    getCountGenoMissing(object, "scan", valid=valid)}
              }
            }
            return(out)
          })

#' Obtain total missing genotype counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total missing genotype counts to the total genotype calls or not.
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getCountGenoMissing",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, prop=FALSE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$countGenoMissing
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
              if(prop){
                out <- out / nscan(object, valid=valid)
              }
            } else {
              out <- object@scanAnnot$countGenoMissing
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
              if(prop){
                out <- out / nsnp(object, valid=valid)
              }
            }
            return(out)
          })

#' Obtain total reference allele counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total reference allele counts to total non missing allele counts or not.
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getCountAlleleRef",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, prop=FALSE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$countAlleleRef
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
              if(prop){
                out <- out / {nscan(object, valid=valid) * 2 -
                    getCountAlleleMissing(object, "snp", valid=valid)}
              }
            } else {
              out <- object@scanAnnot$countAlleleRef
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
              if(prop){
                out <- out / {nsnp(object, valid=valid) * 2 -
                    getCountAlleleMissing(object, "scan", valid=valid)}
              }
            }
            return(out)
          })

#' Obtain total alternative allele counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total alternative allele counts to total non missing allele counts or not.
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getCountAlleleAlt",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, prop=FALSE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$countAlleleAlt
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
              if(prop){
                out <- out / {nscan(object, valid=valid) * 2 -
                    getCountAlleleMissing(object, "snp", valid=valid)}
              }
            } else {
              out <- object@scanAnnot$countAlleleAlt
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
              if(prop){
                out <- out / {nsnp(object, valid=valid) * 2 -
                    getCountAlleleMissing(object, "scan", valid=valid)}
              }
            }
            return(out)
          })

#' Obtain total missing allele counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total missing allele counts to the total allele number or not.
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getCountAlleleMissing",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, prop=FALSE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$countAlleleMissing
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
              if(prop){
                out <- out / {nscan(object, valid=valid) * 2}
              }
            } else {
              out <- object@scanAnnot$countAlleleMissing
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
              if(prop){
                out <- out / {nsnp(object, valid=valid) * 2}
              }
            }
            return(out)
          })

#' Obtain mean values of total reference read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getMeanReadRef",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$meanReadRef
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
            } else {
              out <- object@scanAnnot$meanReadRef
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
            }
            return(out)
          })

#' Obtain mean values of total alternative read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getMeanReadAlt",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$meanReadAlt
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
            } else {
              out <- object@scanAnnot$meanReadAlt
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
            }
            return(out)
          })

#' Obtain standard deviations of total reference read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getSDReadRef",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$sdReadRef
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
            } else {
              out <- object@scanAnnot$sdReadRef
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
            }
            return(out)
          })

#' Obtain standard deviations of total alternative read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getSDReadAlt",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, ...){
            if(target == "snp"){
              out <- object@snpAnnot$sdReadAlt
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidSnp(object)]
              }
            } else {
              out <- object@scanAnnot$sdReadAlt
              if(is.null(out)) {return(NULL)}
              if(valid){
                out <- out[getValidScan(object)]
              }
            }
            return(out)
          })

#' Obtain quantile values of total reference read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param q A numeric value [0-1] to indicate quantile to obtain.
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getQtileReadRef",
          "GbsrGenotypeData",
          function(object, target="snp", q=0.5, valid=TRUE, ...){
            if(target == "snp"){
              pdata <- pData(object@snpAnnot)
              if(q == "all"){
                index <- grepl("qtileReadRef", names(pdata))
              } else {
                index <- names(pdata) %in% paste0("qtileReadRef", q)
              }
              if(all(!index)) {return(NULL)}
              out <- pdata[, index]
              if(valid){
                out <- subset(out, subset=getValidSnp(object))
              }
            } else {
              pdata <- pData(object@scanAnnot)
              if(q == "all"){
                index <- grepl("qtileReadRef", names(pdata))
              } else {
                index <- names(pdata) %in% paste0("qtileReadRef", q)
              }
              if(all(!index)) {return(NULL)}
              out <- pdata[, index]
              if(valid){
                out <- subset(out, subset=getValidScan(object))
              }
            }
            if(length(out) == 1){
              out <- unlist(out)
            }
            return(out)
          })

#' Obtain quantile values of total alternative read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param q A numeric value [0-1] to indicate quantile to obtain.
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getQtileReadAlt",
          "GbsrGenotypeData",
          function(object, target="snp", q=0.5, valid=TRUE, ...){
            if(target == "snp"){
              pdata <- pData(object@snpAnnot)
              if(q == "all"){
                index <- grepl("qtileReadAlt", names(pdata))
              } else {
                index <- names(pdata) %in% paste0("qtileReadAlt", q)
              }
              if(all(!index)) {return(NULL)}
              out <- pdata[, index]
              if(valid){
                out <- subset(out, subset=getValidSnp(object))
              }
            } else {
              pdata <- pData(object@scanAnnot)
              if(q == "all"){
                index <- grepl("qtileReadAlt", names(pdata))
              } else {
                index <- names(pdata) %in% paste0("qtileReadAlt", q)
              }
              if(all(!index)) {return(NULL)}
              out <- pdata[, index]
              if(valid){
                out <- subset(out, subset=getValidScan(object))
              }
            }
            if(length(out) == 1){
              out <- unlist(out)
            }
            return(out)
          })

#' Obtain minor allele frequencies per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getMAF",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, ...){
            if(target == "snp"){
              out <- getCountAlleleRef(object, "snp", prop=TRUE, valid=valid)
              if(is.null(out)) {return(NULL)}
              out <- 0.5 - abs(out - 0.5)
            } else {
              out <- getCountAlleleRef(object, "scan", prop=TRUE, valid=valid)
              if(is.null(out)) {return(NULL)}
              out <- 0.5 - abs(out - 0.5)
            }
            return(out)
          })

#' Obtain minor allele counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setMethod("getMAC",
          "GbsrGenotypeData",
          function(object, target="snp", valid=TRUE, ...){
            if(target == "snp"){
              ac_ref <- getCountAlleleRef(object, "snp", prop=FALSE, valid=valid)
              ac_alt <- getCountAlleleAlt(object, "snp", prop=FALSE, valid=valid)
              out <- ac_ref
              if(is.null(out)) {return(NULL)}
              alt_minor <- ac_ref > ac_alt
              out[alt_minor] <- ac_alt[alt_minor]
            } else {
              ac_ref <- getCountAlleleRef(object, "scan", prop=FALSE, valid=valid)
              ac_alt <- getCountAlleleAlt(object, "scan", prop=FALSE, valid=valid)
              out <- ac_ref
              if(is.null(out)) {return(NULL)}
              alt_minor <- ac_ref > ac_alt
              out[alt_minor] <- ac_alt[alt_minor]
            }
            return(out)
          })

#' Count genotype calls and alleles per sample and per marker.
#'
#' This function calculates several summary statistics of genotype calls and alleles
#' per marker and per sample. Those values will be stored in the SnpAnnotaionDataSet slot
#' and the ScanAnnotationDataSet slot and obtained via getter functions, e.g.
#' [getCountGenoRef()], [getCountAlleleRef()], and [getMAF()].
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param node Either of "raw", "filt", and "cor". See details.
#'
#' @details
#' #' Genotype call data can be obtained from the "genotype" node, the "filt.genotype"
#' node, or the "corrected.genotype" node of the GDS file with `node = "raw"`,
#' `node = "filt"`, or `node = "raw"`, respectively.
#' The [setCallFilter()] function generate filtered genotype call data in the
#' "filt.genotype" node which can be accessed as mentioned above. On the other hand, the
#' "corrected.genotype" node can be generated via the [clean()] function.
#'
#' @examples
#' gdata <- loadGDS("/path/to/GDS.gds")
#' gdata <- countGenotype(gdata)
#' sample_missing_rate <- getCountGenoMissing(gdata, target = "scan", prop = TRUE)
#' marker_minor_allele_freq <- getMAF(gdata, target = "snp")
#' hist(gdata, stats = "missing")
#'
#' @export
#'
setMethod("countGenotype",
          "GbsrGenotypeData",
          function(object, target="both", node = "",...){
            ls_gdsn <- gdsfmt::ls.gdsn(object@data@handler)
            if(node == "cor" & "corrected.genotype" %in% ls_gdsn){
              genotype_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                                  path="corrected.genotype")
            } else if(node != "raw" & object@data@genotypeVar == "filt.genotype"){
              genotype_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                                  path="filt.genotype")
            } else {
              genotype_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                                  path="genotype")
            }

            valid_markers <- getValidSnp(object)
            valid_scans <- getValidScan(object)
            if(node == "cor"){
              have_flipped <- FALSE
            } else {
              have_flipped <- haveFlipped(object)
            }

            if(have_flipped){
              valid_flipped <- getFlipped(object, valid=TRUE)
            }
            sel <- list(valid_scans, valid_markers)

            # Counts per sample
            if(target %in% c("both", "scan")){
              geno_table <- gdsfmt::apply.gdsn(node=genotype_node,
                                               margin=1,
                                               selection=sel,
                                               FUN=function(x){
                                                 if(have_flipped){
                                                   ref <- x == 2
                                                   alt <- x == 0
                                                   x[valid_flipped & ref] <- 0
                                                   x[valid_flipped & alt] <- 2
                                                 }
                                                 x <- factor(x, levels=0:3)
                                                 return(as.integer(table(x)))
                                               },
                                               as.is="list")
              geno_table <- matrix(unlist(geno_table), nrow = 4)

              ## Summarize genotype call counts
              object@scanAnnot$countGenoRef <- NA
              object@scanAnnot$countGenoHet <- NA
              object@scanAnnot$countGenoAlt <- NA
              object@scanAnnot$countGenoMissing <- NA
              object@scanAnnot$countGenoRef[valid_scans] <- geno_table[3, ]
              object@scanAnnot$countGenoHet[valid_scans] <- geno_table[2, ]
              object@scanAnnot$countGenoAlt[valid_scans] <- geno_table[1, ]
              object@scanAnnot$countGenoMissing[valid_scans] <- geno_table[4, ]


              ## Summarize allele counts
              object@scanAnnot$countAlleleRef <- NA
              object@scanAnnot$countAlleleAlt <- NA
              object@scanAnnot$countAlleleMissing <- NA
              object@scanAnnot$countAlleleRef[valid_scans] <- colSums(geno_table * c(0, 1, 2, 0))
              object@scanAnnot$countAlleleAlt[valid_scans] <- colSums(geno_table * c(2, 1, 0, 0))
              object@scanAnnot$countAlleleMissing[valid_scans] <- colSums(geno_table * c(0, 0, 0, 2))
            }

            # Counts per marker
            if(target %in% c("both", "snp")){
              geno_table <- gdsfmt::apply.gdsn(node=genotype_node,
                                               margin=2,
                                               selection=sel,
                                               FUN=function(x){
                                                 x <- factor(x, levels=0:3)
                                                 return(as.integer(table(x)))
                                               },
                                               as.is="list")
              geno_table <- matrix(unlist(geno_table), nrow = 4)

              if(have_flipped){
                flip_ref <- geno_table[1, valid_flipped]
                geno_table[1, valid_flipped] <- geno_table[3, valid_flipped]
                geno_table[3, valid_flipped] <- flip_ref
              }

              ## Summarize genotype call counts
              object@snpAnnot$countGenoRef <- NA
              object#' @param q A numeric value [0-1] to indicate quantile to obtain.@snpAnnot$countGenoHet <- NA
              object@snpAnnot$countGenoAlt <- NA
              object@snpAnnot$countGenoMissing <- NA
              object@snpAnnot$countGenoRef[valid_markers] <- geno_table[3, ]
              object@snpAnnot$countGenoHet[valid_markers] <- geno_table[2, ]
              object@snpAnnot$countGenoAlt[valid_markers] <- geno_table[1, ]
              object@snpAnnot$countGenoMissing[valid_markers] <- geno_table[4, ]

              ## Summarize allele counts
              object@snpAnnot$countAlleleRef <- NA
              object@snpAnnot$countAlleleAlt <- NA
              object@snpAnnot$countAlleleMissing <- NA
              object@snpAnnot$countAlleleRef[valid_markers] <- colSums(geno_table * c(0, 1, 2, 0))
              object@snpAnnot$countAlleleAlt[valid_markers] <- colSums(geno_table * c(2, 1, 0, 0))
              object@snpAnnot$countAlleleMissing[valid_markers] <- colSums(geno_table * c(0, 0, 0, 2))
            }

            return(object)
          }
)

#' Count reads per sample and per marker.
#'
#' This function calculates several summary statistics of read counts
#' per marker and per sample. Those values will be stored in the SnpAnnotaionDataSet slot
#' and the ScanAnnotationDataSet slot and obtained via getter functions, e.g.
#' [getCountReadRef()] and [getCountReadAlt()].
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param node Either of "raw" and "filt". See details.
#'
#' @details
#' Read count data can be obtained from the "annotation/format/AD/data" node or the
#' "annotation/format/AD/filt.data" node of the GDS file with `node = "raw"` or
#' `node = "filt"`, respectively. The [setCallFilter()] function generate filtered
#' read count data in the "annotation/format/AD/filt.data" node which can be accessed as
#' mentioned above.
#'
#' @examples
#' gdata <- loadGDS("/path/to/GDS.gds")
#' gdata <- countRead(gdata)
#' read_depth_per_marker <- getCountRead(gdata, target = "snp")
#' reference_read_freq <- getCountReadRef(gdata, target = "snp", prop = TRUE)
#' hist(gdata, stats = "ad_ref")
#'
#' @export
#'
setMethod("countRead",
          "GbsrGenotypeData",
          function(object, target="both", node = "", ...){
            if(node != "raw" & object@data@genotypeVar == "filt.genotype"){
              ad_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                            path="annotation/format/AD/filt.data")
            } else {
              ad_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                            path="annotation/format/AD/data")
            }
            valid_markers <- getValidSnp(object)
            valid_scans <- getValidScan(object)
            have_flipped <- haveFlipped(object)
            if(have_flipped){
              valid_flipped <- getFlipped(object, valid=TRUE)
            }

            sel <- list(valid_scans, rep(valid_markers, each=2))
            # Counts per sample
            if(target %in% c("both", "scan")){
              # scan_sel <- list(rep(TRUE, nscan(object, valid=FALSE)), rep(valid_markers, each=2))
              read_count <- gdsfmt::apply.gdsn(node=ad_node,
                                               margin=1,
                                               # selection=scan_sel,
                                               selection=sel,
                                               FUN=function(x){
                                                 ref <- x[c(TRUE, FALSE)]
                                                 alt <- x[c(FALSE, TRUE)]
                                                 if(have_flipped){
                                                   flipped <- ref[valid_flipped]
                                                   ref[valid_flipped] <- alt[valid_flipped]
                                                   alt[valid_flipped] <- flipped
                                                 }

                                                 return(c(sum(ref), sum(alt)))
                                               },
                                               as.is="list")
              read_count <- unlist(read_count)
              # read_count <- read_count[rep(valid_scans, each=2)]

              ## Summarize allelic read counts
              object@scanAnnot$countReadRef <- NA
              object@scanAnnot$countReadAlt <- NA
              object@scanAnnot$countReadRef[valid_scans] <- read_count[c(TRUE, FALSE)]
              object@scanAnnot$countReadAlt[valid_scans] <- read_count[c(FALSE, TRUE)]
            }


            # Counts per marker
            if(target %in% c("both", "snp")){
              # sel_snp <- list(valid_scans, rep(valid_markers, each=2))
              read_count <- gdsfmt::apply.gdsn(node=ad_node,
                                               margin=2,
                                               # selection=sel_snp,
                                               selection=sel,
                                               FUN=sum,
                                               as.is="list")
              read_count <- unlist(read_count)

              if(have_flipped){
                valid_flipped <- rep(valid_flipped, each=2)
                flipped <- read_count[valid_flipped]
                ref <- flipped[c(TRUE, FALSE)]
                flipped[c(TRUE, FALSE)] <- flipped[c(FALSE, TRUE)]
                flipped[c(FALSE, TRUE)] <- ref
                read_count[valid_flipped] <- flipped
              }

              ## Summarize allelic read counts
              object@snpAnnot$countReadRef <- NA
              object@snpAnnot$countReadAlt <- NA
              object@snpAnnot$countReadRef[valid_markers] <- read_count[c(TRUE, FALSE)]
              object@snpAnnot$countReadAlt[valid_markers] <- read_count[c(FALSE, TRUE)]
            }

            return(object)
          }
)


#' Calculate mean, standard deviation, and quantile values of reads per sample and per marker.
#'
#' This function calculates several summary statistics of read counts
#' per marker and per sample. Those values will be stored in the SnpAnnotaionDataSet slot
#' and the ScanAnnotationDataSet slot and obtained via getter functions, e.g.
#' [getMeanReadRef()] and [getQtileReadAlt()].
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param q A numeric value [0-1] to indicate quantile to obtain.
#'
#' @details
#' Read count data can be obtained from the "annotation/format/AD/data" node or the
#' "annotation/format/AD/filt.data" node of the GDS file with `node = "raw"` or
#' `node = "filt"`, respectively. The [setCallFilter()] function generate filtered
#' read count data in the "annotation/format/AD/filt.data" node which can be accessed as
#' mentioned above.
#'
#' @examples
#' gdata <- loadGDS("/path/to/GDS.gds")
#' gdata <- calcReadStats(gdata, q = 0.5)
#' mean_reference_read_depth <- getMeanReadRef(gdata, target = "snp")
#' median_reference_read_depth <- getQtileReadAlt(gdata, target = "snp", q = 0.5)
#' hist(gdata, stats = "mean_ref")
#' hist(gdata, stats = "qtile_alt", q = 0.5)
#'
#' @export
#'
setMethod("calcReadStats",
          "GbsrGenotypeData",
          function(object, target="both", q=NULL, ...){

            valid_markers <- getValidSnp(object)
            valid_scans <- getValidScan(object)
            have_flipped <- haveFlipped(object)
            if(have_flipped){
              valid_flipped <- getFlipped(object, valid=TRUE)
            }

            sel <- list(rep(valid_markers, each=2), valid_scans)

            # Calculate normarized allelic counts (counts per million).
            ad_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                          path="annotation/format/AD")
            ad_data_node <- gdsfmt::index.gdsn(node=ad_node,
                                               path="data")
            ls_gdsn <- gdsfmt::ls.gdsn(node=ad_data_node)
            if("norm" %in% ls_gdsn){
              norm_ad <- gdsfmt::index.gdsn(node=ad_data_node,
                                            path="norm")
            } else {
              norm_ad <- gdsfmt::add.gdsn(node=ad_node,
                                          name="norm",
                                          storage="float32",
                                          compress="LZMA_RA",
                                          replace = TRUE)

              gdsfmt::apply.gdsn(node=ad_data_node,
                                 margin=1,
                                 FUN=function(x){
                                   denomi <- sum(x)
                                   if(denomi == 0){
                                     return(x)
                                   } else {
                                     return(x / denomi * 10^6)
                                   }
                                 },
                                 as.is="gdsnode",
                                 target.node=norm_ad)
              gdsfmt::setdim.gdsn(node=norm_ad,
                                  valdim=gdsfmt::objdesp.gdsn(node=ad_data_node)$dim[2:1])
              gdsfmt::readmode.gdsn(node=norm_ad)
            }

            # read dist per sample
            if(target %in% c("both", "scan")){
              read_dist <- gdsfmt::apply.gdsn(node=norm_ad,
                                              margin=2,
                                              selection = sel,
                                              FUN=function(x){
                                                # x <- x[rep(valid_markers, each=2)]
                                                ref <- x[c(TRUE, FALSE)]
                                                alt <- x[c(FALSE, TRUE)]
                                                if(have_flipped){
                                                  flipped <- ref[valid_flipped]
                                                  ref[valid_flipped] <- alt[valid_flipped]
                                                  alt[valid_flipped] <- flipped
                                                }
                                                ref[ref == 0] <- NA
                                                alt[alt == 0] <- NA
                                                return(
                                                  c(mean(ref, na.rm=TRUE),
                                                    sd(ref, na.rm=TRUE),
                                                    quantile(ref, probs = q, na.rm=TRUE),
                                                    mean(alt, na.rm=TRUE),
                                                    sd(alt, na.rm=TRUE),
                                                    quantile(alt, probs = q, na.rm=TRUE))
                                                )
                                              },
                                              as.is="list")
              read_dist <- matrix(unlist(read_dist), ncol = 4 + length(q)*2, byrow=TRUE)
              # read_dist <- read_dist[valid_scans, ]

              ## Summarize read count distribution information.
              if(is.null(q)){
                col_names <- c("meanReadRef", "sdReadRef",
                               "meanReadAlt", "sdReadAlt")
              } else {
                col_names <- c("meanReadRef", "sdReadRef", paste0("qtileReadRef", q),
                               "meanReadAlt", "sdReadAlt", paste0("qtileReadAlt", q))
              }
              pdata <- pData(object@scanAnnot)
              check <- names(pdata) %in% col_names
              if(any(check)){
                pdata <- subset(pdata, select=!check)
              }
              df <- matrix(NA, nrow=nrow(pdata), ncol=ncol(read_dist))
              df[pdata$validScan, ] <- read_dist
              df <- as.data.frame(df)
              names(df) <- col_names
              pData(object@scanAnnot) <- cbind(pdata, df)
            }


            # Counts per marker
            if(target %in% c("both", "snp")){
              read_dist <- gdsfmt::apply.gdsn(node=norm_ad,
                                              margin=1,
                                              selection = sel,
                                              FUN=function(x){
                                                # x <- x[valid_scans]
                                                x[x == 0] <- NA
                                                return(
                                                  c(mean(x, na.rm=TRUE),
                                                    sd(x, na.rm=TRUE),
                                                    quantile(x, probs = q, na.rm=TRUE))
                                                )
                                              },
                                              as.is="list")
              read_dist <- matrix(unlist(read_dist), ncol = 2 + length(q), byrow=TRUE)

              ## Summarize allelic read counts
              if(have_flipped){
                valid_flipped <- rep(valid_flipped, each=2)
                flipped <- read_dist[valid_flipped, ]
                ref <- flipped[c(TRUE, FALSE), ]
                flipped[c(TRUE, FALSE), ] <- flipped[c(FALSE, TRUE), ]
                flipped[c(FALSE, TRUE), ] <- ref
                read_dist[valid_flipped, ] <- flipped
              }

              pdata <- pData(object@snpAnnot)
              check <- names(pdata) %in% col_names
              if(any(check)){
                pdata <- subset(pdata, select=!check)
              }
              read_dist <- cbind(read_dist[c(TRUE, FALSE), ], read_dist[c(FALSE, TRUE), ])
              df <- matrix(NA, nrow=nrow(pdata), ncol=ncol(read_dist))
              df[pdata$validMarker, ] <- read_dist
              df <- as.data.frame(df)
              names(df) <- col_names
              pData(object@snpAnnot) <- cbind(pdata, df)
            }
            return(object)
          }
)

#' Set labels to samples which should be recognized as parents of the population to be subjected to error correction.
#'
#' Specify two or more samples in the dataset as parents of the population. Markers will be filtered out up on your specification.
#'
#' @param object A GbsrGenotypeData object.
#' @param parents A vector of strings with at least length two. The specified strings should match with the samples ID available via [getScanID()].
#' @param flip A logical value to indicate whether markers should be checked for "flip". See details.
#' @param mono A logical value whether to filter out markers which are not monomorphic in parents.
#' @param bi A logical value whether to filter out marekrs which are not biallelic between parents.
#'
#' @details
#' The `clean` function of `GBScleanR` uses read count information of samples and
#' their parents separately to estimate most probable genotype calls of them.
#' Therefore, you must specify proper samples as parents via this function.
#' If you would like to remove SNP markers which are not biallelic and/or
#' not monomorphic in each parent, set `mono = TRUE` and `bi = TRUE`.
#' `flip = TRUE` flips alleles of markers where the alleles expected as reference
#' allele are called as alternative allele. The alleles found in the parent specified as
#' the first element to the `parents` argument are supposed as reference alleles
#' of the markers. If the "expected" reference alleles are not actually called
#' as reference alleles but alternative alleles in the given data. setParents()
#' will automatically labels those markers "flipped".
#' The SnpAnnotatoinDataSet slot sores this information and accessible
#' via [getFlipped()] which gives you a logical vector
#' indicating which markers are labeled as flipped `TRUE` or not flipped `FALSE`.
#' [haveFlipped()] just tells you whether the SnpAnnotatoinDataSet slot has
#' the information of flipped markers or not.
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' gds <- loadGDS("/path/to/GDS.gds")
#' parents <- grep("parent", getScanID(gds), value = TRUE)
#' gds <- setParents(gds, parents = parents, mono = TRUE, bi = TRUE, flip = TRUE)
#' gds <- clean(gds)
#'
#' @export
#'
setMethod("setParents",
          "GbsrGenotypeData",
          function(object, parents, flip, mono, bi, ...){
            if(length(parents) == 0 | any(is.na(parents))){
              stop('Specify valid sample names as parents.')
            }
            if(inherits(parents, "character")){
              id <- getScanID(object, valid = FALSE)
              p_index <- integer()
              for(i in parents){
                p <- which(id %in% i)
                if(length(p) == 0){
                  stop(paste0('No sample named ', i))
                }
                p_index <- c(p_index, p)
              }
            } else if(inherits(parents, "numeric")){
              n_scan <- nscan(object, valid = FALSE)
              check <- n_scan < parents
              if(any(check)){
                stop(paste0("Total number of samples is ", n_scan,
                            ". But you specified ", parents[check]))
              }
              p_index <- parents
            }

            n_parents <- length(p_index)
            object@scanAnnot$parents <- 0
            for(i in 1:n_parents){
              object@scanAnnot$parents[p_index[i]] <- i
            }

            object@scanAnnot$validScan[p_index] <- FALSE

            genotype_node <- gdsfmt::index.gdsn(object@data@handler, object@data@genotypeVar)

            geno <- integer()

            for(i in p_index){
              geno <- rbind(geno,
                            gdsfmt::read.gdsn(genotype_node, start=c(i, 1), count=c(1, -1)))
            }

            # Find markers which are homozygous in each parent and biallelic
            if(mono){
              monomorphic <- apply(geno, 2, function(x){
                return(all(x %in% c(0, 2)))
              })
            } else {
              monomorphic <- TRUE
            }

            if(bi){
              biallelic <- apply(geno, 2, function(x){
                return(length(unique(x)) != 1)
              })
            } else {
              biallelic <- TRUE
            }

            object <- setValidSnp(object, update = monomorphic & biallelic)

            # Find markers which p1 has an alternative allele while p2 has a reference allele.
            if(n_parents == 2 & flip){
              message('Check flipped markers which p1 is called as alternative homozygote.')
              object@snpAnnot$flipped <- geno[1, ] == 0
            }

            return(object)
          }
)

#' Get parental sample information
#'
#' This function returns scan IDs, member IDs and indexes of parental samples
#' set via [setParents()]. Scan IDs are IDs given by user or obtained from the
#' original VCF file. Member IDs are serial numbers assigned by [setParents()].
#'
#' @param object A GbsrGenotypeData object.
#'
#' @export
#'
#' @examples
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = c("parent1", "parent2"))
#' getParents(gds)
#'
setMethod("getParents",
          "GbsrGenotypeData",
          function(object, ...){
            parents <- object@scanAnnot$parents
            if(is.null(parents)){
              stop('No information of parents.')
            }
            p_index <- which(parents != 0)
            p_id <- parents[p_index]
            p_name <- getScanID(object, valid = FALSE)[p_index]
            return(data.frame(scanID = p_name, memberID = p_id, indexes = p_index))
          })

#' Swap the alleles recorded in a GDS file linked to the given GbsrGenotypeData object.
#'
#' The alleles of each marker are automatically obtained to match with those
#' recorded in an input VCF file when it was converted to a GDS file. This function swap
#' those alleles.
#'
#' @param object A GbsrGenotypeData object.
#' @param allele A vector or matrix of characters each of which is either of "A", "T", "G", and "C". The length and the number of rows should be same with the number of "valid" markers. See details.
#'
#' @details
#' The `allele` argument can take a vector or a two-column matrix of characters
#' indicating reference alleles or both alleles. If A vector was given, this
#' function check the current reference and alternative allele of each marker is
#' same with the specified allele for each marker in the vector.
#' If a marker showed that the current alternative allele matched with the allele
#' in the vector, the reference allele and the alternative allele will be swapped
#' each other. In the case of a matrix, the alleles specified in the first column are
#' supposed to be reference alleles while the second column is for alternative alleles.
#' This function compares both alleles between the current record and the specified in
#' the matrix. If a marker showed one of the alleles specified in the allele matrix
#' do not exist in the current allele of the marker, this marker will be labeled as "invalid"
#' marker. In the case of that both of the specified alleles exist but swapped
#' in the current record, the current alleles will be swapped to match with those specified in
#' the allele matrix.
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' In the case of that you have a reference genome data but it was not used for the SNP call, or reference alleles in the genotype data do not match with the alleles in the reference genome, e.g. TASSEL-GBS do.
#' gds <- loadGDS("/path/to/GDS.gds")
#' ref_genome <- Biostrings::readDNAStringSet("/path/to/genome.fasta")
#' chr_names <- getChromosome(object, name = TRUE)
#' snp_pos <- GenomicRanges::GRanges(seqnames = chr_names,
#'                                   ranges = IRanges::IRanges(start = getPosition(object),
#'                                                             width = 1))
#' ref_allele <- as.character(genome[snp_pos])
#' gds <- swapAlleles(gds, allele = ref_allele)
#'
#' @export
#'
setMethod("swapAlleles",
          "GbsrGenotypeData",
          function(object, allele){
            if(is.matrix(allele)){
              if(ncol(allele) != 2){
                stop('The object passed to "allele" should be a data.frame with two columns for reference alleles and alternative alleles.')
              }
              if(nrow(allele) != nsnp(object)){
                stop('The number of markers given as "allele" does not match with that in the data object.')
              }
              a <- getAlleleA(object)
              b <- getAlleleB(object)
              valid <- a == allele[, 1] & b == allele[, 2]
              flipped <- a == allele[, 2] & b == allele[, 1]
              invalid <- !valid & !flipped

            } else if(is.vector(allele)){
              valid <- getAlleleA(object) == allele
              flipped <- getAlleleB(object) == allele
              invalid <- !valid & !flipped

            } else {
              stop('Specify either of allele and genome.')
            }
            message(paste0(sum(invalid), ' SNPs were found as invalid markers which were found to have the 3rd alleles in your genotype data.'))
            message(paste0(sum(flipped), ' SNPs were found as flipped markers.'))

            object <- setValidSnp(object, update = !invalid)
            flipped <- flipped[!invalid]
            # Find markers of which reference allele was called as alternative allele.
            object <- setFlipped(object, flipped)
            return(object)
          })

#' Remove markers potentially having redundant information.
#'
#' Markers within the length of the sequenced reads (usually ~ 150 bp, up to your sequencer)
#' potentially have redundant information and those will cause unexpected errors
#' in error correction which assumes independency of markers each other.
#' This function only retains the first marker or the least missing rate marker
#' from the markers locating within the specified stretch.
#'
#' @param object A GbsrGenotypeData object.
#' @param range A integer value to indicate the stretch to search markers.
#'
#' @details
#' This function search valid markers from the first marker of each chromosome and
#' compare its physical position with a neighbor marker. If the distance between those
#' markers are equal or less then `range`, one of them which has a larger missing rate
#' will be removed (labeled as invalid marker). When the first marker was retained and
#' the second marker was removed as invalid marker, next the distance between the first marker
#' and the third marker will be checked and this cycle is repeated until reaching the
#' end of each chromosome. Run [getValidSnp()] to check the valid SNP markers.
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- thinMarker(gds, range = 150)
#'
#' @export
#'
setMethod("thinMarker",
          "GbsrGenotypeData",
          function(object, range = 150, ...){
            if(is.null(object@snpAnnot$countGenoMissing)){
              stop('Run countGenotype first.')
            }

            missing_count <- getCountGenoMissing(object)

            chr <- getChromosome(object)
            pos <- getPosition(object)
            n_snp <- nsnp(object)
            valid <- rep(TRUE, n_snp)
            i <- 1
            j <- 2
            while(TRUE){
              if(chr[i] == chr[j]){
                mar1 <- pos[i]
                mar2 <- pos[j]
                if(mar2 - mar1<= range){
                  if(missing_count[i] >= missing_count[j]){
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
              if(j > n_snp){
                break
              }
            }
            object <- setValidSnp(object, update=valid)
            return(object)
          })

#' Filter out each genotype call meeting criteria
#'
#' Perform filtering of each genotype call, neither markers nor samples. Each genotype call
#' is supported by its read counts for the reference allele and the alternative allele of
#' a marker of a sample. `setCallFilter()` set missing to the genotype calls which are
#' not reliable enough and set zero to reference and alternative read counts of
#' the genotype calls.
#'
#' @param object A GbsrGenotypeData object.
#' @param dp_count A numeric vector with length two specifying lower and upper limit of total read counts (reference reads + alternative reads).
#' @param ref_count A numeric vector with length two specifying lower and upper limit of reference read counts.
#' @param alt_count A numeric vector with length two specifying lower and upper limit of alternative read counts.
#' @param norm_dp_count A numeric vector with length two specifying lower and upper limit of normalized total read counts (normalized reference reads + normalized alternaitve reads).
#' @param norm_ref_count A numeric vector with length two specifying lower and upper limit of normalized reference read counts.
#' @param norm_alt_count A numeric vector with length two specifying lower and upper limit of normalized alternative read counts
#' @param scan_dp_qtile A numeric vector with length two specifying lower and upper limit of quantile of total read counts in each scan (sample).
#'@param scan_ref_qtile A numeric vector with length two specifying lower and upper limit of quantile of reference read counts in each scan (sample).
#' @param scan_alt_qtile A numeric vector with length two specifying lower and upper limit of quantile of alternative read counts in each scan (sample).
#' @param snp_dp_qtile A numeric vector with length two specifying lower and upper limit of quantile of total read counts in each SNP marker
#'@param snp_ref_qtile A numeric vector with length two specifying lower and upper limit of quantile of reference read counts in each SNP marker.
#' @param snp_alt_qtile A numeric vector with length two specifying lower and upper limit of quantile of alternative read counts in each SNP marker.
#' @param omit_geno A vector of string with combinations of "ref", "het", and "alt" to remove specified genotype calls.
#'
#' @details
#' `norm_dp_count`, `norm_ref_count`, and `norm_alt_count` use normalized read counts which
#' are obtained by dividing each read count by the total read count of each sample.
#' `scan_dp_qtile`, `scan_ref_qtile`, and `scan_alt_qtile` work similarly but use quantile
#' values of read counts of each sample to decide the lower and upper limit of read counts.
#' This function generate two new nodes in the GDS file linked with the given GbsrGenotypeData
#' object. The new nodes "filt.data" in the AD node and "filt.genotype" contains read count
#' data and genotype data after filtering, respectively.
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' # Filter out genotype calls supported by less than 5 reads.
#' gds <- setCallFilter(gds, dp_count = c(5, Inf))
#'
#' # Filter out genotype calls supported by reads less than the 20 percentile of read counts per marker in each sample.
#' gds <- setCallFilter(gds, scan_dp_qtile = c(0.2, 1))
#'
#' # Filter out all reference homozygote genotype calls.
#' gds <- setCallFilter(gds, omit_geno = "ref")
#'
#' @export
#'
setMethod("setCallFilter",
          "GbsrGenotypeData",
          function(object,
                   dp_count=c(0, Inf),
                   ref_count=c(0, Inf),
                   alt_count=c(0, Inf),
                   norm_dp_count=c(0, Inf),
                   norm_ref_count=c(0, Inf),
                   norm_alt_count=c(0, Inf),
                   scan_dp_qtile=c(0, 1),
                   scan_ref_qtile=c(0, 1),
                   scan_alt_qtile=c(0, 1),
                   snp_dp_qtile=c(0, 1),
                   snp_ref_qtile=c(0, 1),
                   snp_alt_qtile=c(0, 1),
                   omit_geno = NULL,
                   ...){

            .initFilt(object)

            ## Quantile filtering on each genotype call
            check1 <- .setScanCallFilter(object = object, dp_count = dp_count,
                                         ref_count = ref_count, alt_count = alt_count,
                                         norm_dp_count = norm_dp_count,
                                         norm_ref_count = norm_ref_count,
                                         norm_alt_count = norm_alt_count,
                                         dp_qtile = scan_dp_qtile,
                                         ref_qtile = scan_ref_qtile,
                                         alt_qtile = scan_alt_qtile,
                                         omit_geno = omit_geno)
            check2 <- .setSnpCallFilter(object = object,
                                        dp_qtile = snp_dp_qtile,
                                        ref_qtile = snp_ref_qtile,
                                        alt_qtile = snp_alt_qtile)

            ## Generate filtered AD data
            if(check1 | check2){
              .makeCallFilteredData(object)
              object@data@genotypeVar <- "filt.genotype"
            }

            return(object)
          })

.initFilt <- function(object){
  ls_gdsn <- gdsfmt::ls.gdsn(object@data@handler)
  if("filt.genotype" %in% ls_gdsn){
    gdsfmt::delete.gdsn(gdsfmt::index.gdsn(object@data@handler,
                                           "filt.genotype"),
                        force = TRUE)
    gdsfmt::delete.gdsn(gdsfmt::index.gdsn(object@data@handler,
                                           "annotation/format/AD/filt.data"),
                        force = TRUE)
  }
}

#' Filter out scans (samples)
#'
#' Search samples which do not meet the criteria and label them as "invalid".
#'
#' @param object A GbsrGenotypeData object.
#' @param id A vector of strings match with scan ID which can be retrieve by `getScanID()`.
#' @param missing A numeric value [0-1] to specify the maximum missing genotype call rate per sample.
#' @param het A numeric value [0-1] to specify the maximum heterozygous genotype call rate per sample.
#' @param mac A integer value to specify the minimum minor allele count per sample.
#' @param maf A numeric value to specify the minimum minor allele frequency per sample.
#' @param ad_ref A numeric vector with length two specifying lower and upper limit of reference read counts per sample.
#' @param ad_alt A numeric vector with length two specifying lower and upper limit of alternative read counts per sample.
#' @param dp A numeric vector with length two specifying lower and upper limit of total read counts per sample.
#' @param mean_ref A numeric vector with length two specifying lower and upper limit of mean of reference read counts per sample.
#' @param mean_alt A numeric vector with length two specifying lower and upper limit of mean of alternative read counts per sample.
#' @param sd_ref A numeric value specifying the upper limit of standard deviation of reference read counts per sample.
#' @param sd_alt A numeric value specifying the upper limit of standard deviation of alternative read counts per sample.
#'
#' @details
#' For `mean_ref`, `mean_alt`, `sd_ref`, and `sd_alt`, this function calculate mean and
#' standard deviation of reads obtained at SNP markers of each sample. If a mean read counts
#' of a sample was smaller than the specified lower limit or larger than the upper limit,
#' this function labels the sample as "invalid". In the case of `sd_ref` and `sd_alt`,
#' standard deviations of read counts of each sample are checked and the samples having a
#' larger standard deviation will be labeled as "invalid". To check valid and invalid
#' samples, run [getValidScan()].
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setScanFilter(gds, id = getScanID(gds)[1:10], missing = 0.2, dp = c(5, Inf))
#'
#' @export
#'
setMethod("setScanFilter",
          "GbsrGenotypeData",
          function(object,
                   id,
                   missing=1,
                   het=c(0, 1),
                   mac=0,
                   maf=0,
                   ad_ref=c(0, Inf),
                   ad_alt=c(0, Inf),
                   dp=c(0, Inf),
                   mean_ref=c(0, Inf),
                   mean_alt=c(0, Inf),
                   sd_ref=Inf,
                   sd_alt=Inf,
                   ...){

            if(missing(id)){
              id <- rep(TRUE, nscan(object, valid=TRUE))
            } else {
              id <- getScanID(object, valid=TRUE) %in% id
            }

            # Filtering for samples.
            ## Missing rate
            scan_missing <- getCountGenoMissing(object, "scan", prop=TRUE, valid=TRUE)
            scan_missing <- .getSubFilter(scan_missing, missing, 1, F, T)

            ## Heterozygosity
            scan_het <- getCountGenoHet(object, "scan", prop=TRUE, valid=TRUE)
            scan_het <- .getSubFilter(scan_het, het, c(0, 1), T, T)

            ## Minor allele count
            scan_mac <- getMAC(object, "scan", valid=TRUE)
            scan_mac <- .getSubFilter(scan_mac, mac, 0, T, T)

            ## Minor allele frequency
            scan_maf <- getMAF(object, "scan", valid=TRUE)
            scan_maf <- .getSubFilter(scan_maf, maf, 0, T, T)

            ## Reference allele read count
            scan_ad_ref <- getCountReadRef(object, "scan", prop=FALSE, valid=TRUE)
            scan_ad_ref <- .getSubFilter(scan_ad_ref, ad_ref, c(0, Inf), T, T)


            ## Alternative allele read count
            scan_ad_alt <- getCountReadAlt(object, "scan", prop=FALSE, valid=TRUE)
            scan_ad_alt <- .getSubFilter(scan_ad_alt, ad_alt, c(0, Inf), T, T)

            ## Total read count
            scan_dp <- getCountRead(object, "scan", prop=FALSE, valid=TRUE)
            scan_dp <- .getSubFilter(scan_dp, dp, c(0, Inf), T, T)

            ## Mean reference allele read count
            scan_mean_ref <- getMeanReadRef(object, "scan", valid=TRUE)
            scan_mean_ref <- .getSubFilter(scan_mean_ref, mean_ref, c(0, Inf), T, T)

            ## Mean alternative allele read count
            scan_mean_alt <- getMeanReadAlt(object, "scan", valid=TRUE)
            scan_mean_alt <- .getSubFilter(scan_mean_alt, mean_alt, c(0, Inf), T, T)

            ## SD of reference allele read count
            scan_sd_ref <- getSDReadRef(object, "scan", valid=TRUE)
            scan_sd_ref <- .getSubFilter(scan_sd_ref, sd_ref, Inf, F, T)

            ## SD of alternative allele read count
            scan_sd_alt <- getSDReadAlt(object, "scan", valid=TRUE)
            scan_sd_alt <- .getSubFilter(scan_sd_alt, sd_alt, Inf, F, T)

            valid <- id & scan_missing & scan_het & scan_mac & scan_maf & scan_ad_ref &
              scan_ad_alt & scan_dp & scan_mean_ref & scan_mean_alt & scan_sd_ref & scan_sd_alt
            object <- setValidScan(object, update=valid)

            pData(object@snpAnnot) <- subset(pData(object@snpAnnot), select=snpID:ploidy)

            return(object)
          })

#' Filter out markers
#'
#' Search markers which do not meet the criteria and label them as "invalid".
#'
#' @param object A GbsrGenotypeData object.
#' @param id A vector of strings match with scan ID which can be retrieve by `getScanID()`.
#' @param missing A numeric value [0-1] to specify the maximum missing genotype call rate per marker
#' @param het A numeric value [0-1] to specify the maximum heterozygous genotype call rate per marker
#' @param mac A integer value to specify the minimum minor allele count per marker
#' @param maf A numeric value to specify the minimum minor allele frequency per marker.
#' @param ad_ref A numeric vector with length two specifying lower and upper limit of reference read counts per marker.
#' @param ad_alt A numeric vector with length two specifying lower and upper limit of alternative read counts per marker.
#' @param dp A numeric vector with length two specifying lower and upper limit of total read counts per marker.
#' @param mean_ref A numeric vector with length two specifying lower and upper limit of mean of reference read counts per marker.
#' @param mean_alt A numeric vector with length two specifying lower and upper limit of mean of alternative read counts per marker.
#' @param sd_ref A numeric value specifying the upper limit of standard deviation of reference read counts per marker.
#' @param sd_alt A numeric value specifying the upper limit of standard deviation of alternative read counts per marker.
#'
#' @details
#' For `mean_ref`, `mean_alt`, `sd_ref`, and `sd_alt`, this function calculate mean and
#' standard deviation of reads obtained for samples at each SNP marker. If a mean read counts
#' of a marker was smaller than the specified lower limit or larger than the upper limit,
#' this function labels the marker as "invalid". In the case of `sd_ref` and `sd_alt`,
#' standard deviations of read counts of each marker are checked and the markers having a
#' larger standard deviation will be labeled as "invalid". To check valid and invalid
#' markers, run [getValidSnp()].
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setSnpFilter(gds, id = getSnpID(gds)[1:1000], missing = 0.2, dp = c(5, Inf))
#'
#' @export
#'
setMethod("setSnpFilter",
          "GbsrGenotypeData",
          function(object,
                   id,
                   missing=1,
                   het=c(0, 1),
                   mac=0,
                   maf=0,
                   ad_ref=c(0, Inf),
                   ad_alt=c(0, Inf),
                   dp=c(0, Inf),
                   mean_ref=c(0, Inf),
                   mean_alt=c(0, Inf),
                   sd_ref=Inf,
                   sd_alt=Inf,
                   ...){

            if(missing(id)){
              id <- rep(TRUE, nsnp(object, valid=TRUE))
            } else {
              id <- !getSnpID(object, valid=TRUE) %in% id
            }

            # Filtering for samples.
            ## Missing rate
            snp_missing <- getCountGenoMissing(object, "snp", prop=TRUE, valid=TRUE)
            snp_missing <- .getSubFilter(snp_missing, missing, 1, F, T)

            ## Heterozygosity
            snp_het <- getCountGenoHet(object, "snp", prop=TRUE, valid=TRUE)
            snp_het <- .getSubFilter(snp_het, het, c(0, 1), T, T)

            ## Minor allele count
            snp_mac <- getMAC(object, "snp", valid=TRUE)
            snp_mac <- .getSubFilter(snp_mac, mac, 0, T, T)

            ## Minor allele frequency
            snp_maf <- getMAF(object, "snp", valid=TRUE)
            snp_maf <- .getSubFilter(snp_maf, maf, 0, T, T)

            ## Reference allele read count
            snp_ad_ref <- getCountReadRef(object, "snp", prop=FALSE, valid=TRUE)
            snp_ad_ref <- .getSubFilter(snp_ad_ref, ad_ref, c(0, Inf), T, T)


            ## Alternative allele read countalt_qtile = scan_alt_qtile
            snp_ad_alt <- getCountReadAlt(object, "snp", prop=FALSE, valid=TRUE)
            snp_ad_alt <- .getSubFilter(snp_ad_alt, ad_alt, c(0, Inf), T, T)

            ## Total read count
            snp_dp <- getCountRead(object, "snp", prop=FALSE, valid=TRUE)
            snp_dp <- .getSubFilter(snp_dp, dp, c(0, Inf), T, T)

            ## Mean reference allele read count
            snp_mean_ref <- getMeanReadRef(object, "snp", valid=TRUE)
            snp_mean_ref <- .getSubFilter(snp_mean_ref, mean_ref, c(0, Inf), T, T)

            ## Mean alternative allele read count
            snp_mean_alt <- getMeanReadAlt(object, "snp", valid=TRUE)
            snp_mean_alt <- .getSubFilter(snp_mean_alt, mean_alt, c(0, Inf), T, T)

            ## SD of reference allele read count
            snp_sd_ref <- getSDReadRef(object, "snp", valid=TRUE)
            snp_sd_ref <- .getSubFilter(snp_sd_ref, sd_ref, Inf, F, T)

            ## SD of alternative allele read count
            snp_sd_alt <- getSDReadAlt(object, "snp", valid=TRUE)
            snp_sd_alt <- .getSubFilter(snp_sd_alt, sd_alt, Inf, F, T)

            valid <- id & snp_missing & snp_het & snp_mac & snp_maf & snp_ad_ref &
              snp_ad_alt & snp_dp & snp_mean_ref & snp_mean_alt & snp_sd_ref & snp_sd_alt
            object <- setValidSnp(object, update=valid)
            pdata <- pData(object@scanAnnot)
            if("parents" %in% names(pdata)){
              pData(object@scanAnnot) <- subset(pdata, select=c(scanID, validScan, parents))
            } else {
              pData(object@scanAnnot) <- subset(pdata, select=c(scanID, validScan))
            }
            return(object)
          }
)

#' Filter out markers based on marker quality metrics
#'
#' A VCF file usually has marker quality metrics in the INFO filed and those are stored in
#' a GDS file created via `GBScleanR`. This function filter out markers based on those marker
#' quality metrics.
#'
#' @param object A GbsrGenotypeData object.
#' @param mq A numeric value to specify minimum mapping quality (shown as MQ in the VCF format).
#' @param fs A numeric value to specify maximum Phred-scaled p-value (strand bias) (shown as FS in the VCF format).
#' @param qd A numeric value to specify minimum Variant Quality by Depth (shown as QD in the VCF format).
#' @param sor A numeric value to specify maximum Symmetric Odds Ratio (strand bias) (shown as SOR in the VCF format).
#' @param mqranksum A numeric values to specify the lower and upper limit of Alt vs. Ref read mapping qualities (shown as MQRankSum in the VCF format).
#' @param readposranksum A numeric values to specify the lower and upper limit of Alt vs. Ref read position bias (shown as ReadPosRankSum in the VCF format).
#' @param baseqranksum A numeric values to specify the lower and upper limit of Alt Vs. Ref base qualities (shown as BaseQRankSum in the VCF format).
#'
#' @details
#' Detailed explanation of each metric can be found in [GATK's web site](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setInfoFilter(gds, mq = 40, qd = 20)
#'
#' @export
#'
setMethod("setInfoFilter",
          "GbsrGenotypeData",
          function(object,
                   mq=0 ,
                   fs=Inf,
                   qd=0,
                   sor=Inf,
                   mqranksum=c(-Inf, Inf),
                   readposranksum=c(-Inf, Inf),
                   baseqranksum=c(-Inf, Inf),
                   ...){

            # Filtering for samples.
            ## MQ
            snp_mq <- getInfo(object, "MQ")
            snp_mq <- .getSubFilter(snp_mq, mq, 0, T, T)

            ## FS
            snp_fs <- getInfo(object, "FS")
            snp_fs <- .getSubFilter(snp_fs, fs, Inf, F, T)

            ## QD
            snp_qd <- getInfo(object, "QD")
            snp_qd <- .getSubFilter(snp_qd, qd, 0, T, T)

            ## SOR
            snp_sor <- getInfo(object, "SOR")
            snp_sor <- .getSubFilter(snp_sor, sor, Inf, F, T)

            ## MQRankSum
            snp_mqranksum <- getInfo(object, "MQRankSum")
            snp_mqranksum <- .getSubFilter(snp_mqranksum, mqranksum, c(-Inf, Inf), T, T)

            ## Alternative allele read count
            snp_readposranksum <- getInfo(object, "ReadPosRankSum")
            snp_readposranksum <- .getSubFilter(snp_readposranksum, readposranksum, c(-Inf, Inf), T, T)

            ## Total read count
            snp_baseqranksum <- getInfo(object, "BaseQRankSum")
            snp_baseqranksum <- .getSubFilter(snp_baseqranksum, baseqranksum, c(-Inf, Inf), T, T)

            valid <- snp_mq & snp_fs & snp_qd & snp_sor &
              snp_mqranksum & snp_readposranksum & snp_baseqranksum
            object <- setValidSnp(object, update=valid)
            if("parents" %in% names(pdata)){
              pData(object@scanAnnot) <- subset(pdata, select=c(scanID, validScan, parents))
            } else {
              pData(object@scanAnnot) <- subset(pdata, select=c(scanID, validScan))
            }
            return(object)
          }
)

#' Reset the filter made by [setCallFiler()]
#'
#' Return genotype calls and read count data to the original data which are same with
#' those data before running [setCallFilter()].
#'
#' @param object A GbsrGenotypeData object.
#'
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setMethod("resetCallFilters",
          "GbsrGenotypeData",
          function(object, ...){
            object@data@genotypeVar <- "genotype"
            return(object)
          })

#' Reset the filter made by [setScanFiler()]
#'
#' Remove "invalid" labels put on samples and make all samples valid.
#'
#' @param object A GbsrGenotypeData object.
#'
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setMethod("resetScanFilters",
          "GbsrGenotypeData",
          function(object, ...){
            object <- setValidScan(object, new = TRUE)
            pdata <- pData(object@scanAnnot)
            if("parents" %in% names(pdata)){
              pData(object@scanAnnot) <- subset(pdata, select=c(scanID, validScan, parents))
            } else {
              pData(object@scanAnnot) <- subset(pdata, select=c(scanID, validScan))
            }
            return(object)
          })

#' Reset the filter made by [setSnpFiler()]
#'
#' Remove "invalid" labels put on markers and make all markers valid.
#'
#' @param object A GbsrGenotypeData object.
#'
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setMethod("resetSnpFilters",
          "GbsrGenotypeData",
          function(object, ...){
            object <- setValidSnp(object, new = TRUE)
            pdata <- pData(object@snpAnnot)
            pdata <- subset(pdata, select=snpID:ploidy)
            pData(object@snpAnnot) <- pdata
            return(object)
          })

#' Reset all filters made by [setScanFiler()], [setSnpFiler()], and [setCallFiler()].
#'
#' Return all data intact.
#'
#' @param object A GbsrGenotypeData object.
#'
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setMethod("resetFilters",
          "GbsrGenotypeData",
          function(object, ...){
            object <- resetCallFilters(object)
            object <- resetScanFilters(object)
            object <- resetSnpFilters(object)
            return(object)
          })

.getSubFilter <- function(variable, threshold, default, greater, equal, ...){
  if(is.null(variable)){
    return(TRUE)
  }
  if(length(default) == 1){
    if("comp" %in% threshold[1]){
      output <- .compareValues(variable, default, greater, F)

    } else if(threshold != default){
      output <- .compareValues(variable, threshold, greater, equal)

    } else {
      output <- TRUE
    }
  } else if(length(default) == 2){
    if("comp" %in% threshold[1]){
      output <- .compareValues(variable, default[1], greater, F) &
        .compareValues(variable, default[2], !greater, F)

    } else if(any(threshold != default)){
      output <- .compareValues(variable, threshold[1], greater, equal) &
        .compareValues(variable, threshold[2], !greater, equal)

    } else {
      output <- TRUE
    }
  }
  return(output)
}

.compareValues <- function(x, y, greater, equal, na2f=TRUE){
  if(greater & equal){
    out <- x >= y
  } else if(greater & !equal){
    out <- x > y
  } else if(!greater & equal){
    out <- x <= y
  } else if(!greater & !equal){
    out <- x < y
  }
  if(na2f){
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
                               qtile_default = c(0, 1),
                               ...){
  ad_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                path="annotation/format/AD")
  ad_data_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                     path="annotation/format/AD/data")

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

  if(!check){

    filt_scan <- gdsfmt::add.gdsn(node=ad_node,
                                  name="filt.scan",
                                  storage="bit1",
                                  compress="LZMA_RA",
                                  replace=TRUE)

    gdsfmt::apply.gdsn(node=ad_data_node,
                       margin=1,
                       target.node=filt_scan,
                       as.is="gdsnode",
                       FUN=function(x){
                         ref <- x[c(TRUE, FALSE)]
                         alt <- x[c(FALSE, TRUE)]
                         dp <- ref + alt
                         ref[!getValidSnp(object)] <- NA
                         alt[!getValidSnp(object)] <- NA
                         dp[!getValidSnp(object)] <- NA

                         if(haveFlipped(object)){
                           flipped <- getFlipped(object, valid=FALSE)
                           flipped_ref <- ref[flipped]
                           ref[flipped] <- alt[flipped]
                           alt[flipped] <- flipped_ref
                         }

                         if(!is.null(omit_geno)){
                           if("ref" %in% omit_geno){
                             omit_ref <- ref > 0 & alt == 0
                           } else {
                             omit_ref <- FALSE
                           }
                           if("alt" %in% omit_geno){
                             omit_alt <- ref == 0 & alt > 0
                           } else {
                             omit_alt <- FALSE
                           }
                           if("het" %in% omit_geno){
                             omit_het <- ref > 0 & alt > 0
                           } else {
                             omit_het <- FALSE
                           }
                           omit_geno <- !omit_ref & !omit_alt & !omit_het
                         } else {
                           omit_geno <- TRUE
                         }

                         ref_count <- .getSubFilter(ref, ref_count, count_default, T, T)
                         alt_count <- .getSubFilter(alt, alt_count, count_default, T, T)
                         dp_count <- .getSubFilter(dp, dp_count, count_default, T, T)

                         if(!all(dp_qtile == qtile_default)){
                           threshold <- c(quantile(dp, probs=dp_qtile[1], na.rm=TRUE),
                                          quantile(dp, probs=dp_qtile[2], na.rm=TRUE))
                           dp_qtile <- .getSubFilter(dp, threshold, c(-1, -1), T, T)
                         } else {
                           dp_qtile <- TRUE
                         }

                         if(!all(ref_qtile == qtile_default)){
                           threshold <- c(quantile(ref, probs=ref_qtile[1], na.rm=TRUE),
                                          quantile(ref, probs=ref_qtile[2], na.rm=TRUE))
                           ref_qtile <- .getSubFilter(ref, threshold, c(-1, -1), T, T)
                         } else {
                           ref_qtile <- TRUE
                         }

                         if(!all(alt_qtile == qtile_default)){
                           threshold <- c(quantile(alt, probs=alt_qtile[1], na.rm=TRUE),
                                          quantile(alt, probs=alt_qtile[2], na.rm=TRUE))
                           alt_qtile <- .getSubFilter(alt, threshold, c(-1, -1), T, T)
                         } else {
                           alt_qtile <- TRUE
                         }

                         if(!all(norm_ref_count == count_default)){
                           denomi <- sum(x)
                           if(denomi != 0){
                             ref <- ref / denomi * 10^6
                           }
                           norm_ref_count <- .getSubFilter(ref, norm_ref_count, count_default, T, T)
                         } else {
                           norm_ref_count <- TRUE
                         }

                         if(!all(norm_alt_count == count_default)){
                           denomi <- sum(x)
                           if(denomi != 0){
                             alt <- alt / denomi * 10^6
                           }
                           norm_alt_count <- .getSubFilter(alt, norm_alt_count, count_default, T, T)
                         } else {
                           norm_alt_count <- TRUE
                         }

                         if(!all(norm_dp_count == count_default)){
                           denomi <- sum(x)
                           if(denomi != 0){
                             dp <- dp / denomi * 10^6
                           }
                           norm_dp_count <- .getSubFilter(dp, norm_dp_count, count_default, T, T)
                         } else {
                           norm_dp_count <- TRUE
                         }
                         return(
                           as.integer(ref_count & ref_qtile &
                                        alt_count & alt_qtile &
                                        dp_count & dp_qtile &
                                        norm_ref_count & norm_alt_count &
                                        norm_dp_count & omit_geno)
                         )
                       }
    )
    gdsfmt::setdim.gdsn(node=filt_scan,
                        valdim=c(nsnp(object, valid=FALSE),
                                 nscan(object, valid=FALSE)))
    return(TRUE)
  } else {
    return(FALSE)
  }
}

.setSnpCallFilter <- function(object,
                              dp_qtile,
                              ref_qtile,
                              alt_qtile,
                              qtile_default = c(0, 1),
                              ...){
  ad_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                path="annotation/format/AD")
  ad_data_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                     path="annotation/format/AD/data")
  n_snp <- nsnp(object, valid=FALSE)
  check <- all(ref_qtile == qtile_default) &
    all(alt_qtile == qtile_default)

  if(!check){
    ref_markers <- rep(c(TRUE, FALSE), n_snp)

    if(haveFlipped(object)){
      flipped <- getFlipped(object, valid=FALSE)
      flipped <- rep(flipped, each=2)
      flipped_markers <- ref_markers[flipped]
      ref <- flipped_markers[c(TRUE, FALSE)]
      flipped_markers[c(TRUE, FALSE)] <- flipped_markers[c(FALSE, TRUE)]
      flipped_markers[c(FALSE, TRUE)] <- ref
      ref_markers[flipped] <- flipped_markers
    }
    sel_ref <- list(rep(TRUE, nscan(object, valid=FALSE)), ref_markers)

    filt_ref <- gdsfmt::add.gdsn(node=ad_node,
                                 name="filt.ref",
                                 storage="bit1",
                                 compress="LZMA_RA",
                                 replace=TRUE)
    gdsfmt::apply.gdsn(node=ad_data_node,
                       margin=2,
                       selection=sel_ref,
                       target.node=filt_ref,
                       as.is="gdsnode",
                       FUN=function(x){
                         x[!getValidScan(object)] <- NA
                         x[x == 0] <- NA
                         ref_qtile <- c(quantile(x, probs=ref_qtile[1], na.rm=TRUE),
                                        quantile(x, probs=ref_qtile[2], na.rm=TRUE))
                         ref_qtile <- .getSubFilter(x, ref_qtile, c(-1, -1), T, T)

                         return(
                           as.integer(ref_qtile)
                         )
                       })

    gdsfmt::setdim.gdsn(node=filt_ref,
                        valdim=c(
                          nscan(object, valid=FALSE),
                          nsnp(object, valid=FALSE)))
    out <- TRUE
  } else {
    out <- FALSE
  }

  check <- all(alt_qtile == qtile_default)

  if(!check){
    alt_markers <- rep(c(FALSE, TRUE), n_snp)

    if(haveFlipped(object)){
      flipped <- getFlipped(object, valid=FALSE)
      flipped <- rep(flipped, each=2)
      flipped_markers <- alt_markers[flipped]
      alt <- flipped_markers[c(TRUE, FALSE)]
      flipped_markers[c(TRUE, FALSE)] <- flipped_markers[c(FALSE, TRUE)]
      flipped_markers[c(FALSE, TRUE)] <- alt
      alt_markers[flipped] <- flipped_markers
    }
    sel_alt <- list(rep(TRUE, nscan(object, valid=FALSE)), alt_markers)

    filt_alt <- gdsfmt::add.gdsn(node=ad_node,
                                 name="filt.alt",
                                 storage="bit1",
                                 compress="LZMA_RA",
                                 replace=TRUE)
    gdsfmt::apply.gdsn(node=ad_data_node,
                       margin=2,
                       selection=sel_alt,
                       target.node=filt_alt,
                       as.is="gdsnode",
                       FUN=function(x){
                         x[!getValidScan(object)] <- NA
                         x[x == 0] <- NA
                         alt_qtile <- c(quantile(x, probs=alt_qtile[1], na.rm=TRUE),
                                        quantile(x, probs=alt_qtile[2], na.rm=TRUE))
                         alt_qtile <- .getSubFilter(x, alt_qtile, c(-1, -1), T, T)

                         return(
                           as.integer(alt_qtile)
                         )
                       })
    gdsfmt::setdim.gdsn(node=filt_alt,
                        valdim=c(
                          nscan(object, valid=FALSE),
                          nsnp(object, valid=FALSE)))
    out <- c(out, TRUE)
  } else {
    out <- c(out, FALSE)
  }
  return(any(out))
}

.makeCallFilteredData <- function(object, ...){
  ad_node <- gdsfmt::index.gdsn(node=object@data@handler, path="annotation/format/AD")
  ad_data <- gdsfmt::index.gdsn(node=object@data@handler, path="annotation/format/AD/data")
  ad_data_desp <- gdsfmt::objdesp.gdsn(node=ad_data)
  ad_filt_data <- gdsfmt::add.gdsn(node=ad_node,
                                   name="filt.data",
                                   storage=ad_data_desp$storage,
                                   compress=ad_data_desp$compress,
                                   replace=TRUE)
  ad_filt_data_tmp <- gdsfmt::add.gdsn(node=ad_node,
                                       name="filt.data_tmp",
                                       storage=ad_data_desp$storage,
                                       compress=ad_data_desp$compress,
                                       replace=TRUE)

  gt_data <- gdsfmt::index.gdsn(node=object@data@handler, path="genotype")
  gt_data_desp <- gdsfmt::objdesp.gdsn(node=gt_data)
  gt_filt_data <- gdsfmt::add.gdsn(node=object@data@handler,
                                   name="filt.genotype",
                                   storage=gt_data_desp$storage,
                                   compress=gt_data_desp$compress,
                                   replace=TRUE)
  gt_data_desp <- gdsfmt::objdesp.gdsn(node=gt_data)
  gt_filt_data_tmp <- gdsfmt::add.gdsn(node=object@data@handler,
                                       name="filt.genotype_tmp",
                                       storage=gt_data_desp$storage,
                                       compress=gt_data_desp$compress,
                                       replace=TRUE)

  n_scan <- nscan(object, valid=FALSE)
  n_snp <- nsnp(object, valid=FALSE)
  validmarker <- getValidSnp(object)
  validmarker_ad <- rep(validmarker, each=2)
  validscan_index <- which(getValidScan(object))
  parents_index <- which(object@scanAnnot$parents != 0)
  sel_scan <- rep(FALSE, n_scan)

  for(i in 1:n_scan){
    if(i %in% validscan_index){
      tmp_scan_i <- sel_scan
      tmp_scan_i[i] <- TRUE
      scan_i_ad <- gdsfmt::readex.gdsn(node=ad_data,
                                       sel=list(tmp_scan_i, rep(TRUE, n_snp * 2)))
      scan_i_gt <- gdsfmt::readex.gdsn(node=gt_data,
                                       sel=list(tmp_scan_i, rep(TRUE, n_snp)))

      scan_i_ad[!validmarker_ad] <- 0
      ad_ls <- gdsfmt::ls.gdsn(node=ad_node)

      if("filt.ref" %in% ad_ls){
        filt_ref <- gdsfmt::index.gdsn(node=ad_node, path="filt.ref")
        gdsfmt::readmode.gdsn(node=filt_ref)
        scan_i_filt_ref <- gdsfmt::readex.gdsn(node=filt_ref,
                                               sel=list(tmp_scan_i, rep(TRUE, n_snp)))

        scan_i_gt[scan_i_filt_ref == 0] <- 3
        scan_i_filt_ref <- rep(scan_i_filt_ref, each=2)
        scan_i_ad[scan_i_filt_ref == 0] <- 0
      }

      if("filt.alt" %in% ad_ls){
        filt_alt <- gdsfmt::index.gdsn(node=ad_node, path="filt.alt")
        gdsfmt::readmode.gdsn(node=filt_alt)
        scan_i_filt_alt <- gdsfmt::readex.gdsn(node=filt_alt,
                                               sel=list(tmp_scan_i, rep(TRUE, n_snp)))
        scan_i_gt[scan_i_filt_alt == 0] <- 3
        scan_i_filt_alt <- rep(scan_i_filt_alt, each=2)
        scan_i_ad[scan_i_filt_alt == 0] <- 0
      }

      if("filt.scan" %in% ad_ls){
        filt_scan <- gdsfmt::index.gdsn(node=ad_node, path="filt.scan")
        gdsfmt::readmode.gdsn(node=filt_scan)
        scan_i_filt_scan <- gdsfmt::readex.gdsn(node=filt_scan,
                                                sel=list(rep(TRUE, n_snp), tmp_scan_i))
        scan_i_gt[scan_i_filt_scan == 0] <- 3
        scan_i_filt_scan <- rep(scan_i_filt_scan, each=2)
        scan_i_ad[scan_i_filt_scan == 0] <- 0
        ref <- scan_i_ad[c(T, F)]
        alt <- scan_i_ad[c(F, T)]
      }
    } else if(i %in% parents_index){
      tmp_scan_i <- sel_scan
      tmp_scan_i[i] <- TRUE
      scan_i_ad <- gdsfmt::readex.gdsn(node=ad_data,
                                       sel=list(tmp_scan_i, rep(TRUE, n_snp * 2)))
      scan_i_gt <- gdsfmt::readex.gdsn(node=gt_data,
                                       sel=list(tmp_scan_i, rep(TRUE, n_snp)))
      scan_i_ad[!validmarker_ad] <- 0
    } else {
      scan_i_ad <- rep(0L, n_snp*2)
      scan_i_gt <- rep(3L, n_snp)
    }

    gdsfmt::append.gdsn(node=ad_filt_data_tmp, val=scan_i_ad)
    gdsfmt::append.gdsn(node=gt_filt_data_tmp, val=scan_i_gt)
  }
  gdsfmt::setdim.gdsn(node=ad_filt_data_tmp, valdim=ad_data_desp$dim[2:1])
  gdsfmt::setdim.gdsn(node=gt_filt_data_tmp, valdim=gt_data_desp$dim[2:1])
  gdsfmt::readmode.gdsn(node=ad_filt_data_tmp)
  gdsfmt::readmode.gdsn(node=gt_filt_data_tmp)
  gdsfmt::apply.gdsn(node=ad_filt_data_tmp, margin=1,
                     FUN=c, as.is="gdsnode", target.node=ad_filt_data)
  gdsfmt::apply.gdsn(node=gt_filt_data_tmp, margin=1,
                     FUN=c, as.is="gdsnode", target.node=gt_filt_data)
  gdsfmt::setdim.gdsn(node=ad_filt_data, valdim=ad_data_desp$dim)
  gdsfmt::setdim.gdsn(node=gt_filt_data, valdim=gt_data_desp$dim)
  gdsfmt::delete.gdsn(node=ad_filt_data_tmp, force=TRUE)
  gdsfmt::delete.gdsn(node=gt_filt_data_tmp, force=TRUE)
  gdsfmt::readmode.gdsn(node=ad_filt_data)
  gdsfmt::readmode.gdsn(node=gt_filt_data)

}

#' Create a GDS file with subset data of the current GDS file
#'
#' Create a new GDS file storing the subset data from the current GDS file linked to
#' the given GbsrGenotypeData object with keeping (or removing) information based on
#' valid markers and samples information.
#'
#' @param object A GbsrGenotypeData object.
#' @param out_fn A string to specify the path to an output GDS file.
#' @param snp_incl A logical vector having the same length with the total number of markers. The values obtained via [getValidSnp()] are used.
#' @param scan_incl A logical vector having the same length with the total number of scans (samples). The values obtained via [getValidScan()] are used.
#' @param incl_parents A logical value to specify whether parental samples should be included in a subset data or not.
#'
#' @details
#' A resultant subset data in a new GDS file includes subsets of each category of data, e.g.
#' genotype, SNP ID, scan ID, read counts, and quality metrics of SNP markers.
#'
#' @return A GbsrGenotypeData object linking to the new GDS file storing subset data.
#'
#' @examples
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setSnpFilter(gds, missing = 0.2, het = c(0.1, 0.9), maf = 0.05)
#' new_gds <- subsetGDS(gds, out_fn = "/path/to/newGDS.gds")
#'
#' @export
#'
setMethod("subsetGDS",
          "GbsrGenotypeData",
          function(object, out_fn = "./susbet.gds", snp_incl, scan_incl, incl_parents = TRUE, ...){

            n_scan <- nscan(object, valid=FALSE)
            n_snp <- nsnp(object, valid=FALSE)

            if(missing(snp_incl)){
              snp_incl <- getValidSnp(object)
            } else {
              check <- n_snp == length(snp_incl)
              if(!check){
                stop('The vector snp_incl should be the same length with the total number of SNPs in the GDS.')
              }
            }

            if(missing(scan_incl)){
              scan_incl <- getValidScan(object)
              if(incl_parents){
                scan_incl[object@scanAnnot$parents != 0] <- TRUE
              }
            } else {
              check <- n_scan == length(scan_incl)
              if(!check){
                stop('The vector scan_incl should be the same length with the total number of scans in the GDS.')
              }
            }

            if(n_scan == sum(scan_incl) & n_snp == sum(snp_incl)){
              message("All markers and sampels are valid.")
              message("Nothing to be subset.")
              return(NULL)
            }

            oldgds <- object@data@handler
            newgds <- gdsfmt::createfn.gds(filename=out_fn)
            attr_old <- gdsfmt::get.attr.gdsn(oldgds$root)
            for(i in 1:length(attr_old)){
              gdsfmt::put.attr.gdsn(node=newgds$root, name=names(attr_old)[i], val=attr_old[[i]])
            }

            ls_node <- gdsfmt::ls.gdsn(oldgds, include.dirs=FALSE)
            for(i_node in ls_node){
              oldgds_i_node <- gdsfmt::index.gdsn(node=oldgds, path=i_node)
              i_desc <- gdsfmt::objdesp.gdsn(oldgds_i_node)
              newgds_i_node <- gdsfmt::add.gdsn(node=newgds,
                                                name=i_node,
                                                storage=i_desc$storage,
                                                compress=i_desc$compress,
                                                replace=TRUE)

              if(length(i_desc$dim) == 1){
                check <- sum(which(c(n_scan, n_snp) %in% i_desc$dim))
                if(check == 1){
                  gdsfmt::assign.gdsn(node=newgds_i_node,
                                      src.node=oldgds_i_node,
                                      seldim=list(scan_incl))
                } else if(check == 2){
                  gdsfmt::assign.gdsn(node=newgds_i_node,
                                      src.node=oldgds_i_node,
                                      seldim=list(snp_incl))
                }
              } else if(length(i_desc$dim) == 2){
                if(i_desc$name == "parents.genotype"){
                  panrets_index <- rep(TRUE, i_desc$dim[1])
                  gdsfmt::assign.gdsn(node=newgds_i_node,
                                      src.node=oldgds_i_node,
                                      seldim=list(panrets_index,
                                                  rep(snp_incl, each = times)))
                } else {
                  times <- i_desc$dim[2] / n_snp
                  gdsfmt::assign.gdsn(node=newgds_i_node,
                                      src.node=oldgds_i_node,
                                      seldim=list(scan_incl,
                                                  rep(snp_incl, each = times)))
                }
              } else {
                next
              }
              newgds_i_attr <- gdsfmt::get.attr.gdsn(node=oldgds_i_node)
              if(!is.null(newgds_i_attr)){
                for(attr_i in 1:length(newgds_i_attr)){
                  gdsfmt::put.attr.gdsn(node=newgds_i_node,
                                        name=names(newgds_i_attr)[attr_i],
                                        val=newgds_i_attr[[attr_i]])
                }
              }
            }

            newgds_anno <- gdsfmt::addfolder.gdsn(node=newgds, name="annotation", replace=TRUE)
            newgds_format <- gdsfmt::addfolder.gdsn(node=newgds_anno, name="format", replace=TRUE)
            newgds_info <- gdsfmt::addfolder.gdsn(node=newgds_anno, name="info", replace=TRUE)
            oldgds_info <- gdsfmt::index.gdsn(node=oldgds, path="annotation/info")
            ls_node <- gdsfmt::ls.gdsn(oldgds_info, include.dirs=FALSE)
            for(i_node in ls_node){
              oldgds_i_node <- gdsfmt::index.gdsn(node=oldgds_info, path=i_node)
              i_desc <- gdsfmt::objdesp.gdsn(oldgds_i_node)
              if(i_desc$dim[1] == 0){
                next
              }
              newgds_i_node <- gdsfmt::add.gdsn(node=newgds_info,
                                                name=i_node,
                                                storage=i_desc$storage,
                                                compress=i_desc$compress,
                                                replace=TRUE)
              gdsfmt::assign.gdsn(node=newgds_i_node,
                                  src.node=oldgds_i_node,
                                  seldim=list(snp_incl))
              newgds_i_attr <- gdsfmt::get.attr.gdsn(node=oldgds_i_node)
              if(!is.null(newgds_i_attr)){
                for(attr_i in 1:length(newgds_i_attr)){
                  gdsfmt::put.attr.gdsn(node=newgds_i_node,
                                        name=names(newgds_i_attr)[attr_i],
                                        val=newgds_i_attr[[attr_i]])
                }
              }
            }

            newgds_format <- gdsfmt::addfolder.gdsn(node=newgds_anno, name="format", replace=TRUE)
            oldgds_format <- gdsfmt::index.gdsn(node=oldgds, path="annotation/format")
            ls_node <- gdsfmt::ls.gdsn(oldgds_format)
            for(i_node in ls_node){
              oldgds_i_node <- gdsfmt::index.gdsn(node=oldgds_format, path=i_node)
              newgds_i_node <- gdsfmt::addfolder.gdsn(node=newgds_format,
                                                      name=i_node,
                                                      replace=TRUE)

              old_attr <- gdsfmt::get.attr.gdsn(oldgds_i_node)
              if(!is.null(old_attr)){
                for(i in 1:length(old_attr))
                  gdsfmt::put.attr.gdsn(newgds_i_node,
                                        name = names(old_attr)[i],
                                        val = old_attr[[i]])
              }

              ls_oldgds_i_node <- gdsfmt::ls.gdsn(oldgds_i_node)
              for(i_node_i in ls_oldgds_i_node){
                oldgds_i_node_i <- gdsfmt::index.gdsn(node=oldgds_i_node, path=i_node_i)
                i_desc <- gdsfmt::objdesp.gdsn(oldgds_i_node_i)
                newgds_i_node_i <- gdsfmt::add.gdsn(node=newgds_i_node,
                                                    name=i_node_i,
                                                    storage=i_desc$storage,
                                                    compress=i_desc$compress,
                                                    replace=TRUE)
                check <- which(i_desc$dim %in% n_scan)
                if(length(check) == 0){ next }
                if(check == 1){
                  each <- i_desc$dim[2] / n_snp
                  seldim <- list(scan_incl, rep(snp_incl, each=each))
                } else {
                  each <- i_desc$dim[1] / n_snp
                  seldim <- list(rep(snp_incl, each=each), scan_incl)
                }
                gdsfmt::assign.gdsn(node=newgds_i_node_i,
                                    src.node=oldgds_i_node_i,
                                    seldim=seldim)
                newgds_i_attr <- gdsfmt::get.attr.gdsn(node=oldgds_i_node_i)
                if(!is.null(newgds_i_attr)){
                  for(attr_i in 1:length(newgds_i_attr)){
                    gdsfmt::put.attr.gdsn(node=newgds_i_node_i,
                                          name=names(newgds_i_attr)[attr_i],
                                          val=newgds_i_attr[[attr_i]])
                  }
                }
              }
            }
            gdsfmt::closefn.gds(newgds)
            gdsfmt::cleanup.gds(newgds$filename, verbose = FALSE)
            output <- loadGDS(gds_fn=newgds$filename)
            if(object@data@genotypeVar == "filt.genotype"){
              replaceGDSdata(output, "genotype", )
            }
            ad <- "filt.data" %in% gdsfmt::ls.gdsn(gdsfmt::index.gdsn(output@data@handler,
                                                                      "annotation/format/AD"))
            if(ad){
              replaceGDSdata(output, "ad")
            }
            return(output)
          }
)

#' Set the filtered data to be used in GBScleanR's functions
#'
#' Set the "filt.genotype" node and the "filt.data" node as primary nodes for genotype
#' data and read count data. The data stored in the primary nodes are used in the
#' functions of GBScleanR.
#'
#' @param object A GbsrGenotypeData object.
#'
#' @details
#' A GbsrGenotypeData object storing information of the primary node of genotype data and
#' read count data. All of the functions implemented in `GBScleanR` check the primary nodes
#' and use data stored in those nodes. [setCallFilter()] create new nodes storing
#' filtered genotype calls and read counts in a GDS file and change the primary nodes to
#' "filt.genotype" and "filt.data" for genotype and read count data, respectively.
#' [SetRawGenotype()] set back the nodes to the original, those are "genotype" and "data" for
#' genotype and read count data, respectively. You can set the filtered data again by
#' `SetFiltGenotype()`.
#'
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setMethod("setFiltGenotype",
          "GbsrGenotypeData",
          function(object, ...){
            object@data@genotypeVar <- "filt.genotype"
            return(object)
          })

#' Set the origina; data to be used in GBScleanR's functions
#'
#' Set the "genotype" node and the "data" node as primary nodes for genotype
#' data and read count data. The data stored in the primary nodes are used in the
#' functions of GBScleanR.
#'
#' @param object A GbsrGenotypeData object.
#'
#' @details
#' A GbsrGenotypeData object storing information of the primary node of genotype data and
#' read count data. All of the functions implemented in `GBScleanR` check the primary nodes
#' and use data stored in those nodes. [setCallFilter()] create new nodes storing
#' filtered genotype calls and read counts in a GDS file and change the primary nodes to
#' "filt.genotype" and "filt.data" for genotype and read count data, respectively.
#' `SetRawGenotype()` set back the nodes to the original, those are "genotype" and "data" for
#' genotype and read count data, respectively. You can set the filtered data again by
#' [SetFiltGenotype()].
#'
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setMethod("setRawGenotype",
          "GbsrGenotypeData",
          function(object, ...){
            object@data@genotypeVar <- "genotype"
            return(object)
          })

#' Replace data stored in the GDS file
#'
#' Replace the original data with the filtered data or replace sample IDs in the GDS file linked
#' to an input GbsrGenotypeData object.
#'
#' @param object A GbsrGenotypeData object.
#' @param target A vector of combinations of "sample.id", "genotype", and "ad".
#' @param node Either one of "filt" and "cor" to replace raw genotype data with filtered genotype data or corrected genotype data, respectively.
#'
#' @details
#' If `target = "genotype"`, replace the data stored in the "genotype" node with the
#' data in the "filt.genotype" node for genotype call data. If `target = "ad`,
#' The data stored in the "data" node is also replaced with the data in the
#' "filt.data" node for read count data. `target = "sample.id"` makes this function to
#' replace data in the "sample.id" node with the sample IDs stored in the
#' ScanAnnotationDataSet slot of the input GbsrGenotypeData object.
#'
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setMethod("replaceGDSdata",
          "GbsrGenotypeData",
          function(object, target, node, ...){
            if("sample.id" %in% target){
              id <- getScanID(object, valid = FALSE)
              gdsfmt::add.gdsn(object@data@handler, "sample.id", id, "string",
                               compress = "LZMA_RA",
                               replace = TRUE)
              gdsfmt::readmode.gdsn(gdsfmt::index.gdsn(object@data@handler, "sample.id"))
            }

            if("genotype" %in% target){
              if(node == "filt" & "filt.genotype" %in% gdsfmt::ls.gdsn(object@data@handler)){
                gt <- gdsfmt::index.gdsn(object@data@handler, "genotype")
                newgt <- gdsfmt::index.gdsn(object@data@handler, "filt.genotype")
                gt_attr <- gdsfmt::get.attr.gdsn(gt)
                gdsfmt::delete.gdsn(gt)
                gdsfmt::rename.gdsn(newgt, "genotype")
                gdsfmt::put.attr.gdsn(newgt, names(gt_attr)[1], gt_attr[[1]])
                gdsfmt::readmode.gdsn(newgt)

              } else if(node == "cor" & "corrected.genotype" %in% gdsfmt::ls.gdsn(object@data@handler)){
                gt <- gdsfmt::index.gdsn(object@data@handler, "genotype")
                newgt <- gdsfmt::index.gdsn(object@data@handler, "corrected.genotype")
                gt_attr <- gdsfmt::get.attr.gdsn(gt)
                gdsfmt::delete.gdsn(gt)
                gdsfmt::rename.gdsn(newgt, "genotype")
                gdsfmt::put.attr.gdsn(newgt, names(gt_attr)[1], gt_attr[[1]])
                gdsfmt::readmode.gdsn(newgt)
              } else {
                warning('Nothing to replace.')
              }
            }
            if("ad" %in% target){
              ad_gdsn <- gdsfmt::index.gdsn(object@data@handler, "annotation/format/AD")
              ls_gdsn <- gdsfmt::ls.gdsn(ad_gdsn)
              if("filt.data" %in% ls_gdsn){
                ad <- gdsfmt::index.gdsn(ad_gdsn, "data")
                fad <- gdsfmt::index.gdsn(ad_gdsn, "filt.data")
                gdsfmt::delete.gdsn(ad)
                gdsfmt::rename.gdsn(fad, "data")
                gdsfmt::readmode.gdsn(fad)
                if(length(ls_gdsn) > 2){
                  ls_gdsn <- ls_gdsn[!ls_gdsn %in% c("data", "filt.data")]
                  for(i in ls_gdsn){
                    gdsfmt::delete.gdsn(gdsfmt::index.gdsn(ad_gdsn, i))
                  }
                }

              } else {
                warning('Nothing to replace.')
              }
            }
          })

#' Build a GbsrScheme object
#'
#' GBScleanR uses breeding scheme information to set the expected
#' number of cross overs in a chromosome which is a required parameter
#' for the genotype error correction with the Hidden Markov model
#' implemented in the "clean" function. This function build the object storing
#' type crosses performed at each generation of breeding and population sizes.
#'
#' @param object A GbsrGenotypeData object.
#' @param crosstype A string to indicate the type of cross conducted with a given generation.
#' @param mating An integer matrix to indicate mating combinations. The each element should match with IDs of parental samples which are 1 to N. see Details.
#'
#' @return A GbsrScheme object.
#'
#' @details
#' A GbsrScheme object stores information of a population size, mating combinations and
#' a type of cross applied to each generation of the breeding process
#' to generate the population which you are going to subject to the "clean" function.
#' The first generation should be parents of the population. It is supposed that
#' [setParents()] has been already executed and parents are labeled in the
#' GbsrGenotypeData object. The number of parents are automatically recognized.
#' The "crosstype" of the first generation can be "pairing" or "random" with
#' `pop_size = N`, where N is the number of parents.
#' You need to specify a matrix indicating combinations of `mating`, in which each column shows
#' a pair of parental samples. For example, if you have only two parents, the `mating` matirix
#' should be `mating = matrix(1:2, nrow = 1, ncol = 2)`. The indices used in the matrix
#' should match with the IDs labeled to parental samples by [setParents()].
#' The created GbsrScheme object is set in the `scheme` slot of the GbsrGenotypeData object.
#'
#' @export
#'
#' @seealso [addScheme()] and [showScheme()]
#'
#' @examples
#' # Biparental F2 population.
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = c("parent1", "parent2"))
#' # setParents gave member ID 1 and 2 to parent1 and parent2, respectively.
#' gds <- initScheme(gds, crosstype = "pair", mating = matrix(1:2, nrow = 1, ncol = 2))
#' # Now the progenies of the cross above have member ID 3.
#' # If `crosstype = "selfing"` or `"sibling"`, you can omit a `mating` matrix.
#' gds <- addScheme(gds, crosstype = "self")
#'
#' 8-way RILs with sibling mating.
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = paste("parent", 1:8, sep = ""))
#' # setParents set member ID 1 to 8 to parent1 to parent8, respectively.
#'
#' # If you made crosses of parent1 x parent2, parent3 x parent4, parent5 x parent6, and parent7 x parent8, run the following.
#' gds <- initScheme(gds, crosstype = "pair", mating = matrix(1:8, nrow = 4, ncol = 2))
#'
#' # Now the progenies of the crosses above have member ID 9, 10, 11, and 12 for each combination of mating.You can check IDs with showScheme().
#' showScheme(gds)
#'
#' # Then, produce 4-way crosses.
#' gds <- addScheme(gds, crosstype = "pair", mating = matrix(9:12, nrow = 2, ncol = 2))
#' # 8-way crosses.
#' gds <- addScheme(gds, crosstype = "pair", mating = matrix(13:14, nrow = 1, ncol = 2))
#' # Inbreeding by 4 times selfing.
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#'
#' # Execute error correction
#' gds <- clean(gds)
#'
setMethod("initScheme",
          "GbsrGenotypeData",
          function(object, crosstype, mating, ...){
            parents <- getParents(object)
            object@scheme <- initScheme(object@scheme, crosstype, mating, parents$memberID)
            return(object)
          })

#' Build a GbsrScheme object
#'
#' GBScleanR uses breeding scheme information to set the expected
#' number of cross overs in a chromosome which is a required parameter
#' for the genotype error correction with the Hidden Markov model
#' implemented in the "clean" function. This function build the object storing
#' type crosses performed at each generation of breeding and population sizes.
#'
#' @param object A GbsrGenotypeData object.
#' @param crosstype A string to indicate the type of cross conducted with a given generation.
#' @param mating An integer matrix to indicate mating combinations. The each element should match with member IDs of the last generation.
#' @param pop_size An integer of the number of individuals in a given generation.
#'
#' @return A GbsrScheme object.
#'
#' @details
#' A scheme object is just a data.frame indicating a population size and
#' a type of cross applied to each generation of the breeding process
#' to generate the population which you are going to subject to the "clean" function.
#' The `crosstype` can take either of "selfing", "sibling", "pairing", and "random".
#' When you set `crosstype = "random"`, you need to specify `pop_size` to indicate how many
#' individuals were crossed in the random mating.
#' You also need to specify a matrix indicating combinations of `mating`, in which
#' each column shows a pair of member IDs indicating parental samples of the cross.
#' Member IDs are serial numbers starts from 1 and automatically assigned by
#' [initScheme()] and [addScheme()]. To check the member IDs, run [showScheme()].
#' Please see the examples section for more details of specifying a `mating` matrix.
#' The created GbsrScheme object is set in the `scheme` slot of the GbsrGenotypeData object.
#'
#' @export
#'
#' @seealso [addScheme()] and [showScheme()]
#'
#' @examples
#' # Biparental F2 population.
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = c("parent1", "parent2"))
#' # setParents gave member ID 1 and 2 to parent1 and parent2, respectively.
#' gds <- initScheme(gds, crosstype = "pair", mating = matrix(1:2, nrow = 1, ncol = 2))
#' # Now the progenies of the cross above have member ID 3.
#' # If `crosstype = "selfing"` or `"sibling"`, you can omit a `mating` matrix.
#' gds <- addScheme(gds, crosstype = "self")
#'
#' 8-way RILs with sibling mating.
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = paste("parent", 1:8, sep = ""))
#' # setParents set member ID 1 to 8 to parent1 to parent8, respectively.
#'
#' # If you made crosses of parent1 x parent2, parent3 x parent4, parent5 x parent6, and parent7 x parent8, run the following.
#' gds <- initScheme(gds, crosstype = "pair", mating = matrix(1:8, nrow = 4, ncol = 2))
#'
#' # Now the progenies of the crosses above have member ID 9, 10, 11, and 12 for each combination of mating.You can check IDs with showScheme().
#' showScheme(gds)
#'
#' # Then, produce 4-way crosses.
#' gds <- addScheme(gds, crosstype = "pair", mating = matrix(9:12, nrow = 2, ncol = 2))
#' # 8-way crosses.
#' gds <- addScheme(gds, crosstype = "pair", mating = matrix(13:14, nrow = 1, ncol = 2))
#' # Inbreeding by 4 times selfing.
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#'
#' # Execute error correction
#' gds <- clean(gds)
#'
setMethod("addScheme",
          "GbsrGenotypeData",
          function(object, crosstype, mating, pop_size, ...){
            if(missing(pop_size)){
              pop_size <- NA
            }
            if(missing(mating)){
              mating <- NA
            }
            object@scheme <- addScheme(object@scheme, crosstype, mating, pop_size)
            return(object)
          })

#' Show the information stored in a GbsrScheme object
#'
#' Print the information of each generation in a GbsrScheme object in the scheme
#' slot of a GbsrGenotypeData object.
#' A GbsrScheme object stores information of a population size, mating combinations and
#' a type of cross applied to each generation of the breeding process
#' to generate the population which you are going to subject to the "clean" function.
#'
#' @param object A GbsrGenotypeData object.
#'
#' @export
#'
#' @seealso [initScheme()] and [addScheme()]
#'
#' @examples
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = paste("parent", 1:8, sep = ""))
#' # setParents set member ID 1 to 8 to parent1 to parent8, respectively.
#' # If you made crosses of parent1 x parent2, parent3 x parent4, parent5 x parent6, and parent7 x parent8, run the following.
#' gds <- initScheme(gds, crosstype = "pair", mating = matrix(1:8, nrow = 4, ncol = 2))
#'
#' Now the progenies of the crosses above have member ID 9, 10, 11, and 12 for each combination of mating. You can check IDs with showScheme().
#' showScheme(gds)
#'
setMethod("showScheme",
          "GbsrGenotypeData",
          function(object, ...){
            parents <- getParents(object)
            showScheme(object@scheme, parents$scanID)
          })
