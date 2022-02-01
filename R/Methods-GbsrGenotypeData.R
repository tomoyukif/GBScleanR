###############################################################################
## Internally used getter functions

## Get the file name of the GDS file connected to the GbsrGenotypeData object.
.getGDSFileName <- function(object){
    return(object@data@handler$filename)
}

## Get the gds.class object in the GbsrGenotypeData object.
.getGdsfmtObj <- function(object){
    if(inherits(object, "gds.class")){
        return(object)
    } else {
        return(object@data@handler)
    }
}

## Get the index of the specified path in the GDS file.
.getNodeIndex <- function(object, path){
    return(index.gdsn(.getGdsfmtObj(object), path))
}

## Get the scheme object in the GbsrGenotypeData object.
.getSchemeObj <- function(object){
    return(object@scheme)
}

## Get the genotypeVar value in the GenotypeData object in
## the GbsrGenotypeData object.
.getGenotypeVar <- function(object){
    return(object@data@genotypeVar)
}

## Get the statistical summary values in the SnpAnnotationDataFrame and
## the ScanAnnotationDataFrame.
.getStatsData <- function(object, varname, target, valid){
    if(target == "snp"){
        f_has <- get("hasSnpVariable")
        f_get <- get("getSnpVariable")
        f_valid <- get("getValidSnp")

    } else {
        f_has <- get("hasScanVariable")
        f_get <- get("getScanVariable")
        f_valid <- get("getValidScan")
    }

    if(!f_has(object, varname)){return(NULL)}
    out <- f_get(object, varname)

    if(valid){out <- out[f_valid(object)]}
    return(out)
}

## Commonly used code for getter functions.
#' @importMethodsFrom GWASTools getSnpVariable
.getterCommon <- function(object, varname, valid, chr){
    if(!is.null(chr) & valid){
        sel <- getValidSnp(object) & getChromosome(object, FALSE) %in% chr

    } else if(is.null(chr) & valid){
        sel <- getValidSnp(object)

    } else if(!is.null(chr) & !valid){
        sel <- getChromosome(object, FALSE) %in% chr

    } else {
        return(getSnpVariable(object, varname))
    }

    if(all(!sel)){
        stop("No markers in chromosome ", paste(chr, collapse = ", "),
             call. = FALSE)
    } else {
        return(getSnpVariable(object, varname, sel))
    }
}

###############################################################################
## Interbally used getter functions which returns a boolean value.
.existGdsNode <- function(object, path){
    return(exist.gdsn(.getGdsfmtObj(object), path))
}

## Check if the GbsrScheme has valid information.
.hasScheme <- function(object){
    scheme <- .getSchemeObj(object)
    return(length(scheme@crosstype) != 0)
}
###############################################################################
## Internally used setter functions.

## Set a new scheme object to the GbsrGenotypeData object.
.setSchemeObj <- function(object, scheme){
    object@scheme <- scheme
    return(object)
}

## Set a new genotypeVar value to the GenotypeData object in
## the GbsrGenotypeData object.
.setGenotypeVar <- function(object, var){
    object@data@genotypeVar <- var
    return(object)
}

## Set values to SnpAnnotationDataFrame object.
#' @importFrom GWASTools SnpAnnotationDataFrame
#' @importMethodsFrom GWASTools getSnpAnnotation hasSnpVariable
.setSnpVariable <- function(object, varname, values){
    if(varname == "df"){
        ann <- SnpAnnotationDataFrame(values)

    } else {
        ann <- getSnpAnnotation(object)
        if(!hasSnpVariable(object, varname)){
            ann[[varname]] <- values

        } else {
            if(is.null(values)){
                ann[[varname]] <- values

            } else if(length(values) == 1){
                ann[[varname]] <- values

            } else {
                stopifnot(length(ann[[varname]]) == length(values))
                ann[[varname]] <- values
            }
        }
    }
    object@snpAnnot <- ann
    return(object)
}

## Set values to ScanAnnotationDataFrame object.
#' @importFrom GWASTools ScanAnnotationDataFrame
#' @importMethodsFrom GWASTools getScanAnnotation hasScanVariable
.setScanVariable <- function(object, varname, values){
    if(varname == "df"){
        ann <- ScanAnnotationDataFrame(values)

    } else {
        ann <- getScanAnnotation(object)
        if(!hasScanVariable(object, varname)){
            ann[[varname]] <- values

        } else {
            if(is.null(values)){
                ann[[varname]] <- values

            } else if(length(values) == 1){
                ann[[varname]] <- values

            } else {
                stopifnot(length(ann[[varname]]) == length(values))
                ann[[varname]] <- values
            }
        }
    }
    object@scanAnnot <- ann
    return(object)
}

## Set logical values to indicate which marker has flipped alleles.
.setFlipped <- function(object, flipped){
    if(is.null(flipped)){
        object <- .setSnpVariable(object, "flipped", NULL)
    }

    if(any(is.na(flipped))){
        stop('NA is not allowed for a logical vector "flipped".', call. = FALSE)
    }

    if(nsnp(object, FALSE) == length(flipped)){
        object <- .setSnpVariable(object, "flipped", flipped)

    } else {
        if(nsnp(object) == length(flipped)){
            tmp <- rep(FALSE, nsnp(object, FALSE))
            tmp[getValidSnp(object)] <- flipped
            object <- .setSnpVariable(object, "flipped", tmp)

        } else {
            stop('The length of "flipped" does not match with the number of SNPs.',
                 call. = FALSE)
        }
    }
    return(object)
}

###############################################################################
## Internally used functions to control the GDS file.

## Decompress GDS file nodes.
.gds_decomp <- function(object){
    if(inherits(object, "GbsrGenotypeData")){
        root <- object@data@handler

    } else if(inherits(object, "gds.class")){
        root <- object$root
    }

    ls_gdsn <- ls.gdsn(root, recursive = TRUE, include.dirs = FALSE)
    for(i in ls_gdsn){
        i_node <- index.gdsn(root, i)
        compression.gdsn(i_node, "")
    }
}

## Compress GDS file nodes.
.gds_comp <- function(object){
    if(inherits(object, "GbsrGenotypeData")){
        root <- object@data@handler

    } else if(inherits(object, "gds.class")){
        root <- object$root
    }

    ls_gdsn <- ls.gdsn(root, FALSE, TRUE, FALSE)
    for(i in ls_gdsn){
        i_node <- index.gdsn(root, i)
        compression.gdsn(i_node, "LZMA_RA")
    }
    for(i in ls_gdsn){
        i_node <- index.gdsn(root, i)
        readmode.gdsn(i_node)
    }
}

###############################################################################
## Exported getter functions which return basic information of the data.

## Get the number of SNPs.
#' @rdname nsnp
setMethod("nsnp",
          "GbsrGenotypeData",
          function(object, valid, chr){
              out <- getValidSnp(object)
              if(!is.null(chr)){
                  out <- out[getChromosome(object, FALSE) == chr]
                  if(length(out) == 0){
                      stop("No markers in chromosome ", paste(chr, collapse = ", "),
                           call. = FALSE)
                  }
              }
              if(valid){
                  out <- sum(out)
              } else {
                  out <- length(out)
              }
              return(out)
          })

## Get the number of scans (samples).
#' @rdname nscan
setMethod("nscan",
          "GbsrGenotypeData",
          function(object, valid){
              out <- getValidScan(object)
              if(valid){
                  out <- sum(out)
              } else {
                  out <- length(out)
              }
              return(out)
          })

## Get the logical vector indicating valid SNPs.
#' @rdname getValidSnp
#' @importMethodsFrom GWASTools getSnpVariable
setMethod("getValidSnp",
          "GbsrGenotypeData",
          function(object, chr){
              out <- getSnpVariable(object, "validMarker")
              if(!is.null(chr)){
                  sel <- getChromosome(object, FALSE) %in% chr
                  if(all(!sel)){
                      stop("No markers in chromosome ", paste(chr, collapse = ", "),
                           call. = FALSE)
                  } else {
                      out <- out[sel]
                  }
              }
              return(out)
          })

## Get the logical vector indicating valid scans (samples).
#' @rdname getValidScan
#' @importMethodsFrom GWASTools getScanVariable
setMethod("getValidScan",
          "GbsrGenotypeData",
          function(object, parents){
              if(parents == "only"){
                  return(getParents(object, TRUE))
              }
              if(parents){
                  out <- getScanVariable(object, "validScan")
                  out[getParents(object, TRUE)] <- TRUE
                  return(out)
              } else {
                  return(getScanVariable(object, "validScan"))
              }
          })

## Get the logical vector indicating flipped markers.
#' @rdname getFlipped
setMethod("getFlipped",
          "GbsrGenotypeData",
          function(object, valid, chr){
              if(!hasFlipped(object)){
                  message('No data of flipped genotype markers.')
                  return(NULL)
              }
              return(.getterCommon(object, "flipped", valid, chr))
          })

## Get the chromosome names in actual strings, indices, or levels.
#' @rdname getChromosome
setMethod("getChromosome",
          "GbsrGenotypeData",
          function(object, valid, levels, name){
              if(name){
                  varname <- "chromosome.name"
              } else {
                  varname <- "chromosome"
              }
              out <- .getterCommon(object, varname, valid, NULL)

              if(levels){
                  out <- unique(out)
              }
              return(out)
          })

## Get the marker positions.
#' @rdname getPosition
setMethod("getPosition",
          "GbsrGenotypeData",
          function(object, valid, chr){
              return(.getterCommon(object, "position", valid, chr))
          })

## Get the reference alleles of markers.
#' @rdname getAlleleA
setMethod("getAlleleA",
          "GbsrGenotypeData",
          function(object, valid, chr){
              out <- .getterCommon(object, "alleleA", valid, chr)

              if(hasFlipped(object)){
                  flipped <- getFlipped(object, valid, chr)
                  b <- .getterCommon(object, "alleleB", valid, chr)
                  out[flipped] <- b[flipped]
              }
              return(out)
          })

## Get the alternative alleles of markers.
#' @rdname getAlleleB
setMethod("getAlleleB",
          "GbsrGenotypeData",
          function(object, valid, chr){
              out <- .getterCommon(object, "alleleB", valid, chr)

              if(hasFlipped(object)){
                  flipped <- getFlipped(object, valid, chr)
                  b <- .getterCommon(object, "alleleA", valid, chr)
                  out[flipped] <- b[flipped]
              }
              return(out)
          })

## Get the SNP marker IDs.
#' @rdname getSnpID
setMethod("getSnpID",
          "GbsrGenotypeData",
          function(object, valid, chr){
              return(.getterCommon(object, "snpID", valid, chr))
          })

## Get the scan (sample) IDs.
#' @rdname getScanID
#' @importMethodsFrom GWASTools getScanVariable
setMethod("getScanID",
          "GbsrGenotypeData",
          function(object, valid){
              if(valid){
                  out <- getScanVariable(object, "scanID", getValidScan(object))

              } else {
                  out <- getScanVariable(object, "scanID")
              }
              return(out)
          })

## Get the ploidies of markers.
#' @rdname getPloidy
#' @importFrom utils head
setMethod("getPloidy",
          "GbsrGenotypeData",
          function(object, valid, chr){
              out <- .getterCommon(object, "ploidy", valid, chr)
              chr_i <- getChromosome(object, valid)
              chr_l <- getChromosome(object, valid, TRUE)
              if(!is.null(chr)){
                  chr_i <- chr_i[chr_i %in% chr]
                  chr_l <- chr_l[chr_l %in% chr]
              }
              out <- tapply(out, chr_i, head, n = 1)
              names(out) <- chr_l
              return(out)
          })

## Get the values of a variable in the "annotation/info" directory.
#' @rdname getInfo
setMethod("getInfo",
          "GbsrGenotypeData",
          function(object, var, valid, chr){
              path <- paste0("annotation/info/", var)

              if(!.existGdsNode(object, path)){
                  warning("No data at ", path)
                  return(NULL)
              }
              info_node <- .getNodeIndex(object, path)
              if(valid){
                  sel <- list(getValidSnp(object) & getChromosome(object, FALSE) == chr)
                  out <- readex.gdsn(info_node, sel)

              } else {
                  out <- read.gdsn(info_node)
              }
              return(out)
          })

## Get read count data from the annotation/format/AD/data node
## in the GDS file connected to the GbsrGenotypeData object.
#' @rdname getRead
setMethod("getRead",
          "GbsrGenotypeData",
          function(object, valid, chr, node, parents){

              if(node == "norm"){
                  path <- paste0("annotation/format/AD/", "norm")
              } else if(node == "filt"){
                  path <- paste0("annotation/format/AD/", "filt.data")
              } else {
                  path <- paste0("annotation/format/AD/", "data")
              }
              if(!.existGdsNode(object, path)){
                  stop("No data for ", node, call. = FALSE)
              }
              ad_node <- .getNodeIndex(object, path)

              if(parents == "only"){
                  valid_samples <- getParents(object, TRUE)

              } else {
                  if(valid){
                      valid_samples <- getValidScan(object)

                  } else {
                      valid_samples <- rep(TRUE, nscan(object, valid))
                  }
                  if(parents){
                      valid_samples[getParents(object, TRUE)] <- TRUE
                  }
              }

              if(valid){
                  valid_snp <- getValidSnp(object)

              } else {
                  valid_snp <- rep(TRUE, nsnp(object, valid))
              }

              if(!is.null(chr)){
                  chr_sel <- getChromosome(object, FALSE) %in% chr
                  if(all(!chr_sel)){
                      stop("No markers in chromosome ",
                           paste(chr, collapse = ", "), call. = FALSE)
                  }
                  valid_snp <- valid_snp & chr_sel
              }
              sel <- list(valid_samples, rep(valid_snp, each=2))

              ad <- readex.gdsn(ad_node, sel)
              ref <- ad[, c(TRUE, FALSE)]
              alt <- ad[, c(FALSE, TRUE)]

              if(hasFlipped(object)){
                  flipped <- getFlipped(object, valid, chr)
                  flipped_ref <- ref[, flipped]
                  ref[, flipped] <- alt[, flipped]
                  alt[, flipped] <- flipped_ref
              }

              rownames(ref) <- getScanID(object, FALSE)[valid_samples]
              rownames(alt) <- getScanID(object, FALSE)[valid_samples]
              colnames(ref) <- getSnpID(object, FALSE)[valid_snp]
              colnames(alt) <- getSnpID(object, FALSE)[valid_snp]
              return(list(ref = ref, alt = alt))
          })

## Get read count data from one of the nodes `genotype`, `filt.genotype`,
## `corrected.genotype`, or `parents.genotype` in the GDS file
## connected to the GbsrGenotypeData object.
#' @rdname getGenotype
#'
setMethod("getGenotype",
          "GbsrGenotypeData",
          function(object, valid, chr, node, parents){
              node <- match.arg(node,
                                c("raw",
                                  "parents.genotype",
                                  "corrected.genotype",
                                  "filt.genotype"))
              if(node == "raw"){node <- "genotype"}
              path <- ifelse(.existGdsNode(object, node),
                             node,
                             stop("No data for ", node, call. = FALSE))

              genotype_node <- .getNodeIndex(object, path)
              if(valid){
                  valid_snp <- getValidSnp(object)

              } else {
                  valid_snp <- rep(TRUE, nsnp(object, valid))
              }

              if(!is.null(chr)){
                  chr_sel <- getChromosome(object, FALSE) %in% chr
                  if(all(!chr_sel)){
                      stop("No markers in chromosome ",
                           paste(chr, collapse = ", "), call. = FALSE)
                  }
                  valid_snp <- valid_snp & chr_sel
              }

              if(node == "parents.genotype"){
                  n_row <- objdesp.gdsn(genotype_node)$dim[1]
                  sel <- list(rep(TRUE, n_row), valid_snp)
                  genotype <- readex.gdsn(genotype_node, sel)
                  parent_i <- getParents(object, TRUE)
                  p_id <- getScanID(object, FALSE)[parent_i]
                  rownames(genotype) <- paste(rep(p_id, each=2),
                                              seq_len(2), sep="_")
                  colnames(genotype) <- getSnpID(object, FALSE)[valid_snp]

              } else {
                  if(parents == "only"){
                      valid_samples <- getParents(object, TRUE)

                  } else {
                      if(valid){
                          valid_samples <- getValidScan(object)

                      } else {
                          valid_samples <- rep(TRUE, nscan(object, valid))
                      }
                      if(parents){
                          valid_samples[getParents(object, TRUE)] <- TRUE
                      }
                  }

                  sel <- list(valid_samples, valid_snp)
                  genotype <- readex.gdsn(genotype_node, sel)
                  rownames(genotype) <- getScanID(object, FALSE)[valid_samples]
                  colnames(genotype) <- getSnpID(object, FALSE)[valid_snp]
              }

              ploidy <- max(getPloidy(object))
              genotype[genotype > ploidy] <- NA

              if(hasFlipped(object)){
                  flipped <- getFlipped(object, valid, chr)
                  tmp <- genotype[, flipped]
                  flip_ref <- tmp == 0
                  flip_alt <- tmp == 2
                  tmp[flip_ref] <- 2
                  tmp[flip_alt] <- 0
                  genotype[, flipped] <- tmp
              }

              return(genotype)
          })

#' @rdname getHaplotype
setMethod("getHaplotype",
          "GbsrGenotypeData",
          function(object, chr, parents){
              if(!.existGdsNode(object, "estimated.haplotype")){
                  stop('No haplotype data, Run clean() to estimate haplotype.',
                       call. = FALSE)
              }

              hap_node <- .getNodeIndex(object, "estimated.haplotype")
              valid_samples <- getValidScan(object)
              if(parents == "only"){
                  valid_samples <- getParents(object, TRUE)

              } else if(parents){
                  valid_samples[getParents(object, TRUE)] <- TRUE
              }

              valid_snp <- getValidSnp(object)
              if(!is.null(chr)){
                  valid_snp[getChromosome(object) != chr] <- FALSE
              }

              sel <- list(valid_samples, rep(valid_snp, each=2))
              hap <- readex.gdsn(hap_node, sel)
              hap[hap == 0] <- NA
              hap_dim <- c(2, sum(valid_snp), sum(valid_samples))
              hap_names <- list(c("1", "2"),
                                getSnpID(object, FALSE)[valid_snp],
                                getScanID(object, FALSE)[valid_samples])
              hap <- array(t(hap), hap_dim)
              return(hap)
          })

###############################################################################
## Exported getter functions which return
## statistical summary information obtained via
## countGenotype(), countRead(), and calcReadStats().

## Get the number of reference reads per marker and per sample.
#' @rdname getCountReadRef
setMethod("getCountReadRef",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countReadRef", target, valid)
              if(is.null(out)){return(NULL)}

              if(prop){
                  tmp <- .getStatsData(object, "countReadAlt", target, valid)
                  out <- out / (out + tmp)
              }
              return(out)
          })

## Get the number of alternative reads per marker and per sample.
#' @rdname getCountReadAlt
setMethod("getCountReadAlt",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countReadAlt", target, valid)
              if(is.null(out)){return(NULL)}

              if(prop){
                  tmp <- .getStatsData(object, "countReadRef", target, valid)
                  out <- out / (out + tmp)
              }
              return(out)
          })

## Get the number of total reads per marker and per sample.
#' @rdname getCountRead
setMethod("getCountRead",
          "GbsrGenotypeData",
          function(object, target, valid){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))

              ref <- .getStatsData(object, "countReadRef", target, valid)
              alt <- .getStatsData(object, "countReadAlt", target, valid)

              if(is.null(ref) | is.null(alt)){return(NULL)}
              out <- ref + alt
              return(out)
          })

## Get the number of reference homozygous genotype calls
## per marker and per sample.
#' @rdname getCountGenoRef
setMethod("getCountGenoRef",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countGenoRef", target, valid)
              if(is.null(out)){return(NULL)}

              if(prop){
                  tmp <- .getStatsData(object, "countGenoMissing", target, valid)
                  if(target == "snp"){
                      nonmissing <- nscan(object, valid) - tmp
                  } else {
                      nonmissing <- nsnp(object, valid) - tmp
                  }
                  out <- out / nonmissing
              }
              return(out)
          })

## Get the number of heterozygous genotype calls per marker and per sample.
#' @rdname getCountGenoHet
setMethod("getCountGenoHet",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countGenoHet", target, valid)
              if(is.null(out)){return(NULL)}

              if(prop){
                  tmp <- .getStatsData(object, "countGenoMissing", target, valid)
                  if(target == "snp"){
                      nonmissing <- nscan(object, valid) - tmp
                  } else {
                      nonmissing <- nsnp(object, valid) - tmp
                  }
                  out <- out / nonmissing
              }
              return(out)
          })

## Get the number of alternative homozygous genotype calls
## per marker and per sample.
#' @rdname getCountGenoAlt
setMethod("getCountGenoAlt",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countGenoAlt", target, valid)
              if(is.null(out)){return(NULL)}

              if(prop){
                  tmp <- .getStatsData(object, "countGenoMissing", target, valid)
                  if(target == "snp"){
                      nonmissing <- nscan(object, valid) - tmp
                  } else {
                      nonmissing <- nsnp(object, valid) - tmp
                  }
                  out <- out / nonmissing
              }
              return(out)
          })

## Get the number of missing genotype calls per marker and per sample.
#' @rdname getCountGenoMissing
setMethod("getCountGenoMissing",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countGenoMissing", target, valid)
              if(is.null(out)){return(NULL)}

              if(prop){
                  if(target == "snp"){
                      out <- out / nscan(object, valid)
                  } else {
                      out <- out / nsnp(object, valid)
                  }
              }
              return(out)
          })

## Get the number of reference alleles per marker and per sample.
#' @rdname getCountAlleleRef
setMethod("getCountAlleleRef",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countAlleleRef", target, valid)
              if(is.null(out)){return(NULL)}

              if(prop){
                  tmp <- .getStatsData(object, "countAlleleMissing", target, valid)
                  if(target == "snp"){
                      nonmissing <- nscan(object, valid) * 2 - tmp
                  } else {
                      nonmissing <- nsnp(object, valid) * 2 - tmp
                  }
                  out <- out / nonmissing
              }
              return(out)
          })

## Get the number of alternative alleles per marker and per sample.
#' @rdname getCountAlleleAlt
setMethod("getCountAlleleAlt",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countAlleleAlt", target, valid)
              if(is.null(out)){return(NULL)}

              if(prop){
                  tmp <- .getStatsData(object, "countAlleleMissing", target, valid)
                  if(target == "snp"){
                      nonmissing <- nscan(object, valid) * 2 - tmp
                  } else {
                      nonmissing <- nsnp(object, valid) * 2 - tmp
                  }
                  out <- out / nonmissing
              }
              return(out)
          })

## Get the number of missing alleles per marker and per sample.
#' @rdname getCountAlleleMissing
setMethod("getCountAlleleMissing",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countAlleleMissing", target, valid)
              if(is.null(out)){return(NULL)}

              if(prop){
                  if(target == "snp"){
                      out <- out / (nscan(object, valid) * 2)
                  } else {
                      out <- out / (nsnp(object, valid) * 2)
                  }
              }
              return(out)
          })

## Get the number of mean reference allele reads per marker and per sample.
#' @rdname getMeanReadRef
setMethod("getMeanReadRef",
          "GbsrGenotypeData",
          function(object, target, valid){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))

              out <- .getStatsData(object, "meanReadRef", target, valid)
              if(is.null(out)){return(NULL)}

              return(out)
          })

## Get the number of mean alternative allele reads per marker and per sample.
#' @rdname getMeanReadAlt
setMethod("getMeanReadAlt",
          "GbsrGenotypeData",
          function(object, target, valid){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))

              out <- .getStatsData(object, "meanReadAlt", target, valid)
              if(is.null(out)){return(NULL)}

              return(out)
          })

## Get SD of the number of reference allele reads per marker and per sample.
#' @rdname getSDReadRef
setMethod("getSDReadRef",
          "GbsrGenotypeData",
          function(object, target, valid){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))

              out <- .getStatsData(object, "sdReadRef", target, valid)
              if(is.null(out)){return(NULL)}

              return(out)
          })

## Get SD of the number of alternative allele reads per marker and per sample.
#' @rdname getSDReadAlt
setMethod("getSDReadAlt",
          "GbsrGenotypeData",
          function(object, target, valid){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))

              out <- .getStatsData(object, "sdReadAlt", target, valid)
              if(is.null(out)){return(NULL)}

              return(out)
          })

## Get quantile values of the number of reference allele reads
## per marker and per sample.
#' @rdname getQtileReadRef
setMethod("getQtileReadRef",
          "GbsrGenotypeData",
          function(object, target, q, valid){
              stopifnot(length(q) != 0)
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))

              out <- .getStatsData(object,
                                   paste0("qtileReadRef", q),
                                   target,
                                   valid)
              if(is.null(out)){return(NULL)}

              return(out)
          })

## Get quantile values of the number of alternative allele reads
## per marker and per sample.
#' @rdname getQtileReadAlt
setMethod("getQtileReadAlt",
          "GbsrGenotypeData",
          function(object, target, q, valid){
              stopifnot(length(q) != 0)
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))

              out <- .getStatsData(object,
                                   paste0("qtileReadAlt", q),
                                   target,
                                   valid)
              if(is.null(out)){return(NULL)}

              return(out)
          })

## Get minor allele frequencies per marker and per sample.
#' @rdname getMAF
setMethod("getMAF",
          "GbsrGenotypeData",
          function(object, target, valid){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))

              out <- getCountAlleleRef(object, target, valid, TRUE)
              out <- 0.5 - abs(out - 0.5)
              return(out)
          })

## Get minor allele counts per marker and per sample.
#' @rdname getMAC
setMethod("getMAC",
          "GbsrGenotypeData",
          function(object, target, valid){
              stopifnot(target %in% c("snp", "scan"))
              stopifnot(is.logical(valid))

              ref <- getCountAlleleRef(object, target, valid, FALSE)
              alt <- getCountAlleleAlt(object, target, valid, FALSE)
              out <- ref
              alt_minor <- ref > alt
              out[alt_minor] <- alt[alt_minor]
              return(out)
          })

## Get the information of the parents.
#' @rdname getParents
#' @importMethodsFrom GWASTools getScanVariable
setMethod("getParents",
          "GbsrGenotypeData",
          function(object, bool){
              parents <- getScanVariable(object, "parents")
              if(is.null(parents)){
                  return(NULL)
              }

              if(bool){
                  return(parents != 0)
              }

              p_index <- which(parents != 0)
              p_id <- parents[p_index]
              p_name <- getScanID(object, FALSE)[p_index]
              return(data.frame(scanID = p_name,
                                memberID = p_id,
                                indexes = p_index
              ))
          })
###############################################################################
## Exported getter functions which return a boolen value.

## Check if the connection to the GDS file is open.
#' @rdname isOpenGDS
setMethod("isOpenGDS",
          "GbsrGenotypeData",
          function(object){
              tryout <- try(diagnosis.gds(.getGdsfmtObj(object)),
                            silent = TRUE)
              if(is.list(tryout)){
                  out <- TRUE

              } else if(grepl("The GDS file is closed", tryout)){
                  out <- FALSE

              } else {
                  out <- tryout[1]
              }
              return(out)
          })

## Check if the data include flipped SNP markers.
#' @rdname hasFlipped
#' @importMethodsFrom GWASTools hasSnpVariable
setMethod("hasFlipped",
          "GbsrGenotypeData",
          function(object){
              return(hasSnpVariable(object, "flipped"))
          })

###############################################################################
## Exported setter functions.

## Set new validity information of SNPs.
#' @rdname setValidSnp
setMethod("setValidSnp",
          "GbsrGenotypeData",
          function(object, new, update){
              if(!missing(new)){
                  if(any(is.na(new))){
                      stop('NA is not allowed for a logical vector "new".',
                           call. = FALSE)
                  }
                  object <- .setSnpVariable(object, "validMarker", new)

              } else if(!missing(update)){
                  if(any(is.na(update))){
                      stop('NA is not allowed for a logical vector "update".',
                           call. = FALSE)
                  }
                  valid_marker <- getSnpVariable(object, "validMarker")
                  valid_marker[getValidSnp(object)] <- update
                  object <- .setSnpVariable(object, "validMarker", valid_marker)
              }
              return(object)
          })

## Set new validity information of scans (samples).
#' @rdname setValidScan
setMethod("setValidScan",
          "GbsrGenotypeData",
          function(object, new, update){
              if(!missing(new)){
                  if(any(is.na(new))){
                      stop('NA is not allowed for a logical vector "new".',
                           call. = FALSE)
                  }
                  object <- .setScanVariable(object, "validScan", new)

              } else if(!missing(update)){
                  if(any(is.na(update))){
                      stop('NA is not allowed for a logical vector "update".',
                           call. = FALSE)
                  }
                  valid_marker <- getScanVariable(object, "validScan")
                  valid_marker[getValidScan(object)] <- update
                  object <- .setScanVariable(object, "validScan", valid_marker)
              }
              return(object)
          })

## Set parental samples.
#' @rdname setParents
setMethod("setParents",
          "GbsrGenotypeData",
          function(object, parents, flip, mono, bi){
              if(length(parents) == 0 | any(is.na(parents))){
                  stop('Specify valid sample names as parents.', call. = FALSE)
              }
              if(inherits(parents, "numeric")){
                  parents <- getScanID(object)[parents]
              }
              if(!inherits(parents, "character")){
                  stop('Specify valid sample names as parents.', call. = FALSE)
              }

              id <- getScanID(object, FALSE)
              p_index <- match(parents, id)
              if(any(is.na(p_index))){
                  stop("No sample named: ", parents[is.na(p_index)], call. = FALSE)
              }

              n_parents <- length(p_index)
              p_vec <- integer(nscan(object, FALSE))
              for (i in seq_len(n_parents)){
                  p_vec[p_index[i]] <- i
              }

              valid_scan <- getValidScan(object, FALSE)
              if(hasScanVariable(object, "parents")){
                  old_parents <- getParents(object, TRUE)
                  valid_scan[old_parents] <- TRUE
              }

              valid_scan[p_vec != 0] <- FALSE
              object <- setValidScan(object, valid_scan)
              object <- .setScanVariable(object, "parents", p_vec)

              if(mono | bi | flip){
                  object <- .pGenoFilt(object, mono, bi, flip)
              }
              return(object)
          })

.pGenoFilt <- function(object, mono, bi, flip){
    geno <- readex.gdsn(.getNodeIndex(object, .getGenotypeVar(object)),
                        list(getParents(object, TRUE), getValidSnp(object)))

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
            return(sum(c(0:2) %in% unique(x)) == 2)
        })
    } else {
        biallelic <- TRUE
    }
    object <- setValidSnp(object, update=monomorphic & biallelic)

    ## Find markers which p1 has an alternative allele
    ## while p2 has a reference allele.
    if(sum(getParents(object, TRUE)) == 2 & flip){
        message('Check flipped markers.')
        object <- .setFlipped(object, geno[1, ] == 0)
    }
    return(object)
}

###############################################################################
## Show the GbsrGenotypeData object.
#' @importMethodsFrom GWASTools getSnpAnnotation
#' @importMethodsFrom GWASTools getScanAnnotation
setMethod("show",
          "GbsrGenotypeData",
          function(object){
              message('Data in GDS file...')
              message('GDS file name')
              message(.getGDSFileName(object))
              if(isOpenGDS(object)){
                  print(.getGdsfmtObj(object))
              } else {
                  message('Connection to GDS file has been closed.')
              }
              message('-----------------------------------------------------')
              message('SnpAnnotationDataFrame')
              print(getSnpAnnotation(object))
              message('-----------------------------------------------------')
              message('ScanAnnotationDataFrame')
              print(getScanAnnotation(object))
          })

###############################################################################
## Functions to communicate with the GDS file.

## Close the connection to the GDS file.
#' @rdname closeGDS
setMethod("closeGDS",
          "GbsrGenotypeData",
          function(object, verbose){
              closefn.gds(.getGdsfmtObj(object))
              if(verbose){
                  message('The connection to the GDS file was closed.')
              }
          })

## Close the connection to the GDS file.
#' @rdname openGDS
#' @importClassesFrom GWASTools GenotypeData GdsGenotypeReader
#' @importFrom GWASTools GdsGenotypeReader GenotypeData
#' @importMethodsFrom GWASTools getSnpAnnotation getScanAnnotation
setMethod("openGDS",
          "GbsrGenotypeData",
          function(object){
              if(isOpenGDS(object)){
                  closeGDS(object, FALSE)
              }
              genotype_var <- .getGenotypeVar(object)

              # Leave the following connection open to build the GdsGenotypeReader
              # object with `readonly=FALSE` mode.
              gds <- GdsGenotypeReader(openfn.gds(.getGDSFileName(object), FALSE),
                                       "scan,snp")
              gds <- GenotypeData(gds,
                                  getSnpAnnotation(object),
                                  getScanAnnotation(object))
              object <- new("GbsrGenotypeData", gds)
              if(genotype_var == "filt.genotype"){
                  object <- setFiltGenotype(object)
              }
              return(object)
          })

## Save the data in the SnpAnnotationDataFrame to the GDS file.
#' @rdname saveSnpAnnot
#' @importFrom Biobase pData
#' @importMethodsFrom GWASTools getSnpAnnotation
setMethod("saveSnpAnnot",
          "GbsrGenotypeData",
          function(object){
              .gds_decomp(object)
              new_node <- add.gdsn(
                  node = .getGdsfmtObj(object),
                  name = "snpAnnot",
                  val = pData(getSnpAnnotation(object)),
                  replace = TRUE
              )
              .gds_comp(object)
              return(object)
          })

## Save the data of the ScanAnnotationDataFrame object to the GDS file.
#' @rdname saveScanAnnot
#' @importFrom Biobase pData
#' @importMethodsFrom GWASTools getScanAnnotation
setMethod("saveScanAnnot",
          "GbsrGenotypeData",
          function(object){
              .gds_decomp(object)
              new_node <- add.gdsn(
                  node = .getGdsfmtObj(object),
                  name = "scanAnnot",
                  val = pData(getScanAnnotation(object)),
                  replace = TRUE
              )
              .gds_comp(object)
              return(object)
          })

## Load the stored SnpAnnotationDataFrame data from the GDS file.
#' @rdname loadSnpAnnot
setMethod("loadSnpAnnot",
          "GbsrGenotypeData",
          function(object){
              if(.existGdsNode(object, "snpAnnot")){
                  ann_node <- .getNodeIndex(object, "snpAnnot")
                  ann <- read.gdsn(ann_node)
                  object <- .setSnpVariable(object, "df", ann)
                  return(object)

              } else {
                  stop('No data in the GDS file.', call. = FALSE)
              }
          })

## Load the stored ScanAnnotationDataFrame data from the GDS file.
#' @rdname loadScanAnnot
setMethod("loadScanAnnot",
          "GbsrGenotypeData",
          function(object){
              if(.existGdsNode(object, "scanAnnot")){
                  ann_node <- .getNodeIndex(object, "scanAnnot")
                  ann <- read.gdsn(ann_node)
                  object <- .setScanVariable(object, "df", ann)
                  return(object)

              } else {
                  stop('No data in the GDS file.', call. = FALSE)
              }
          })

###############################################################################
## Functions to calculate statistical summaries of
## genotype calls and read counts.

## Calculate the numbers of each genotype and each allele per marker and
## per sample.
#' @rdname countGenotype
setMethod("countGenotype",
          "GbsrGenotypeData",
          function(object, target, node){
              node <- match.arg(node,
                                c("raw",
                                  "corrected.genotype",
                                  "filt.genotype"))
              if(node == "raw"){node <- "genotype"}
              if(.existGdsNode(object, node)){
                  path <- .getNodeIndex(object, node)
              }

              valid_markers <- getValidSnp(object)
              valid_scans <- getValidScan(object)

              if(node != "corrected.genotype" & hasFlipped(object)){
                  has_flipped <- TRUE
                  valid_flipped <- getFlipped(object, valid = TRUE)
              } else {
                  has_flipped <- FALSE
                  valid_flipped <- FALSE
              }

              sel <- list(valid_scans, valid_markers)

              # Counts per sample
              if(target %in% c("both", "scan")){
                  object <- .countGenotypeScan(object,
                                               path,
                                               sel,
                                               has_flipped,
                                               valid_flipped)
              }

              # Counts per marker
              if(target %in% c("both", "snp")){
                  object <- .countGenotypeSnp(object,
                                              path,
                                              sel,
                                              has_flipped,
                                              valid_flipped)
              }

              return(object)
          })

.countGenotypeScan <- function(object, path, sel, has_flipped, valid_flipped){
    df <- apply.gdsn(path, 1, selection=sel, as.is="list",
                     FUN=function(x){
                         if(has_flipped){
                             flipped_x <- x[valid_flipped]
                             if(length(flipped_x) > 0){
                                 ref <- flipped_x == 0
                                 alt <- flipped_x == 2
                                 flipped_x[ref] <- 2
                                 flipped_x[alt] <- 0
                                 x[valid_flipped] <- flipped_x
                             }
                         }
                         x <- factor(x, 0:3)
                         return(as.integer(table(x)))
                     })
    df <- matrix(unlist(df), 4)
    tmp <- matrix(NA, 4, nscan(object, FALSE))
    tmp[, getValidScan(object)] <- df
    df <- tmp

    object <- .setScanVariable(object, "countGenoRef", df[3,])
    object <- .setScanVariable(object, "countGenoHet", df[2,])
    object <- .setScanVariable(object, "countGenoAlt", df[1,])
    object <- .setScanVariable(object, "countGenoMissing", df[4,])
    ref_allele <- colSums(df * c(0, 1, 2, 0))
    object <- .setScanVariable(object, "countAlleleRef", ref_allele)
    alt_allele <- colSums(df * c(2, 1, 0, 0))
    object <- .setScanVariable(object, "countAlleleAlt", alt_allele)
    missing <- colSums(df * c(0, 0, 0, 2))
    object <- .setScanVariable(object, "countAlleleMissing", missing)
    return(object)
}

.countGenotypeSnp <- function(object, path, sel, has_flipped, valid_flipped){
    df <- apply.gdsn(path, 2, selection=sel, as.is="list",
                     FUN=function(x){
                         x <- factor(x, 0:3)
                         return(as.integer(table(x)))
                     })
    df <- matrix(unlist(df), 4)

    if(has_flipped){
        tmp <- df[1, valid_flipped]
        df[1, valid_flipped] <- df[3, valid_flipped]
        df[3, valid_flipped] <- tmp
    }
    tmp <- matrix(NA, 4, nsnp(object, FALSE))
    tmp[, getValidSnp(object)] <- df
    df <- tmp

    object <- .setSnpVariable(object, "countGenoRef", df[3,])
    object <- .setSnpVariable(object, "countGenoHet", df[2,])
    object <- .setSnpVariable(object, "countGenoAlt", df[1,])
    object <- .setSnpVariable(object, "countGenoMissing", df[4,])
    ref_allele <- colSums(df * c(0, 1, 2, 0))
    object <- .setSnpVariable(object, "countAlleleRef", ref_allele)
    alt_allele <- colSums(df * c(2, 1, 0, 0))
    object <- .setSnpVariable(object, "countAlleleAlt", alt_allele)
    missing <- colSums(df * c(0, 0, 0, 2))
    object <- .setSnpVariable(object, "countAlleleMissing", missing)
    return(object)
}

## Calculate the numbers of reads of each allele per marker and per sample.
#' @rdname countRead
setMethod("countRead",
          "GbsrGenotypeData",
          function(object, target, node){
              node <- match.arg(node,
                                c("raw",
                                  "data",
                                  "filt.data"))
              if(node == "raw"){node <- "data"}
              path <- paste0("annotation/format/AD/", node)
              if(.existGdsNode(object, path)){
                  path <- .getNodeIndex(object, path)
              }

              valid_markers <- getValidSnp(object)
              valid_scans <- getValidScan(object)

              has_flipped <- hasFlipped(object)
              if(has_flipped){
                  valid_flipped <- getFlipped(object, TRUE)
              }

              sel <- list(valid_scans, rep(valid_markers, each=2))
              # Counts per sample
              if(target %in% c("both", "scan")){
                  object <- .countReadScan(object,
                                           path,
                                           sel,
                                           has_flipped,
                                           valid_flipped)
              }

              # Counts per marker
              if(target %in% c("both", "snp")){
                  object <- .countReadSnp(object,
                                          path,
                                          sel,
                                          has_flipped,
                                          valid_flipped)
              }
              return(object)
          })

.countReadScan <- function(object, path, sel, has_flipped, valid_flipped){
    df <- apply.gdsn(path, 1, selection=sel, as.is="list",
                     FUN=function(x){
                         ref <- x[c(TRUE, FALSE)]
                         alt <- x[c(FALSE, TRUE)]
                         if(has_flipped){
                             tmp <- ref[valid_flipped]
                             ref[valid_flipped] <- alt[valid_flipped]
                             alt[valid_flipped] <- tmp
                         }

                         return(c(sum(ref), sum(alt)))
                     })
    df <- do.call("rbind", df)
    tmp <- rep(NA, nscan(object, FALSE))
    tmp[getValidScan(object)] <- df[, c(TRUE, FALSE)]
    object <- .setScanVariable(object, "countReadRef", tmp)

    tmp <- rep(NA, nscan(object, FALSE))
    tmp[getValidScan(object)] <- df[, c(FALSE, TRUE)]
    object <- .setScanVariable(object, "countReadAlt", tmp)
    return(object)
}

.countReadSnp <- function(object, path, sel, has_flipped, valid_flipped){
    df <- apply.gdsn(path, 2, sum, sel, "list")
    df <- unlist(df)

    if(has_flipped){
        valid_flipped <- rep(valid_flipped, each=2)
        tmp <- df[valid_flipped]
        ref <- tmp[c(TRUE, FALSE)]
        tmp[c(TRUE, FALSE)] <- tmp[c(FALSE, TRUE)]
        tmp[c(FALSE, TRUE)] <- ref
        df[valid_flipped] <- tmp
    }

    tmp <- rep(NA, nsnp(object, FALSE))
    tmp[getValidSnp(object)] <- df[c(TRUE, FALSE)]
    object <- .setSnpVariable(object, "countReadRef", tmp)

    tmp <- rep(NA, nsnp(object, FALSE))
    tmp[getValidSnp(object)] <- df[c(FALSE, TRUE)]
    object <- .setSnpVariable(object, "countReadAlt", tmp)
    return(object)
}

## Calculate mean, SD, and quantile values of normalized allele read counts
## per marker and per sample.
#' @rdname calcReadStats
setMethod("calcReadStats",
          "GbsrGenotypeData",
          function(object, target, q){
              .gds_decomp(object)
              on.exit({
                  .gds_comp(object)
              })

              valid_markers <- getValidSnp(object)
              valid_scans <- getValidScan(object)
              has_flipped <- hasFlipped(object)
              if(has_flipped){
                  valid_flipped <- getFlipped(object, TRUE)
              }

              # Calculate normarized allele read counts (reads per million).
              object <- .calcNormRead(object)
              path <- .getNodeIndex(object, "annotation/format/AD/norm")
              sel <- list(valid_scans, rep(valid_markers, each=2))
              # read dist per sample
              if(target %in% c("both", "scan")){
                  object <- .calcReadStatsScan(object,
                                               path,
                                               sel,
                                               has_flipped,
                                               valid_flipped,
                                               q)
              }

              # Counts per marker
              if(target %in% c("both", "snp")){
                  object <- .calcReadStatsSnp(object,
                                              path,
                                              sel,
                                              has_flipped,
                                              valid_flipped,
                                              q)
              }

              return(object)
          })

.calcNormRead <- function(object){
    ad_data_node <- .getNodeIndex(object, "annotation/format/AD/data")
    ad_node <- .getNodeIndex(object, "annotation/format/AD")
    tmp_node <- add.gdsn(ad_node, "tmp", storage="float32", replace=TRUE)

    apply.gdsn(ad_data_node, 1, as.is="gdsnode", target.node=tmp_node,
               FUN = function(x){
                   denomi <- sum(x)
                   if(denomi == 0){
                       return(x)
                   } else {
                       return(x / denomi * 10^6)
                   }
               })
    setdim.gdsn(tmp_node, objdesp.gdsn(ad_data_node)$dim[2:1])
    norm_node <- add.gdsn(ad_node, "norm", storage="float32", replace=TRUE)
    apply.gdsn(tmp_node, 1, as.is="gdsnode", target.node=norm_node, FUN = c)
    setdim.gdsn(norm_node, objdesp.gdsn(ad_data_node)$dim)

    return(object)
}

.calcReadStatsScan <- function(object,
                               path,
                               sel,
                               has_flipped,
                               valid_flipped,
                               q){
    df <- apply.gdsn(path, 1, selection=sel, as.is="list",
                     FUN=function(x){
                         ref <- x[c(TRUE, FALSE)]
                         alt <-  x[c(FALSE, TRUE)]
                         if(has_flipped){
                             tmp <- ref[valid_flipped]
                             ref[valid_flipped] <- alt[valid_flipped]
                             alt[valid_flipped] <- tmp
                         }
                         ref[ref == 0] <- NA
                         alt[alt == 0] <- NA
                         return(c(mean(ref, na.rm=TRUE),
                                  mean(alt, na.rm=TRUE),
                                  sd(ref, na.rm=TRUE),
                                  sd(alt, na.rm=TRUE),
                                  quantile(ref, q, TRUE),
                                  quantile(alt, q, TRUE)))}
    )
    df <- do.call("rbind", df)
    tmp <- matrix(NA, nscan(object, FALSE), ncol(df))
    tmp[getValidScan(object), ] <- df
    df <- tmp

    object <- .setScanVariable(object, "meanReadRef", df[, 1])
    object <- .setScanVariable(object, "meanReadAlt", df[, 2])
    object <- .setScanVariable(object, "sdReadRef", df[, 3])
    object <- .setScanVariable(object, "sdReadAlt", df[, 4])
    n_q <- length(q)
    if(n_q != 0){
        for(i in seq_len(n_q)){
            object <- .setScanVariable(object, paste0("qtileReadRef", q[i]),
                                       df[, 4 + i])
            object <- .setScanVariable(object, paste0("qtileReadAlt", q[i]),
                                       df[, 4 + n_q + i])
        }
    }
    return(object)
}

.calcReadStatsSnp <- function(object,
                              path,
                              sel,
                              has_flipped,
                              valid_flipped,
                              q){
    df <- apply.gdsn(path, 2, selection=sel, as.is="list",
                     FUN = function(x){
                         x[x == 0] <- NA
                         return(c(mean(x, na.rm=TRUE),
                                  sd(x, na.rm=TRUE),
                                  quantile(x, q, TRUE)))
                     }
    )
    df <- do.call("rbind", df)
    tmp <- matrix(NA, nsnp(object, FALSE) * 2, ncol(df))
    tmp[rep(getValidSnp(object), each=2), ] <- df
    df <- tmp

    ## Summarize allelic read counts
    if(has_flipped){
        valid_flipped <- rep(valid_flipped, each=2)
        tmp <- df[valid_flipped,]
        ref <- tmp[c(TRUE, FALSE),]
        tmp[c(TRUE, FALSE),] <- tmp[c(FALSE, TRUE),]
        tmp[c(FALSE, TRUE),] <- ref
        df[valid_flipped,] <- tmp
    }

    object <- .setSnpVariable(object, "meanReadRef", df[c(TRUE, FALSE), 1])
    object <- .setSnpVariable(object, "meanReadAlt", df[c(FALSE, TRUE), 1])
    object <- .setSnpVariable(object, "sdReadRef", df[c(TRUE, FALSE), 2])
    object <- .setSnpVariable(object, "sdReadAlt", df[c(FALSE, TRUE), 2])
    n_q <- length(q)
    if(n_q != 0){
        for(i in seq_len(n_q)){
            object <- .setSnpVariable(object, paste0("qtileReadRef", q[i]),
                                      df[c(TRUE, FALSE), 2 + i])
            object <- .setSnpVariable(object, paste0("qtileReadAlt", q[i]),
                                      df[c(FALSE, TRUE), 2 + i])
        }
    }
    return(object)
}

###############################################################################
## Filtering functions.

## Thinout markers.
#' @rdname thinMarker
setMethod("thinMarker",
          "GbsrGenotypeData",
          function(object, range){
              if(!hasSnpVariable(object, "countGenoMissing")){
                  stop('Run countGenotype first.', call. = FALSE)
              }

              if(!{is.numeric(range) & length(range) == 1 & range >= 0}){
                  stop("range should be a number greater or equal to zero.",
                       call. = FALSE)
              }

              missing_count <- getCountGenoMissing(object, "snp")

              chr <- getChromosome(object)
              pos <- getPosition(object)
              n_snp <- nsnp(object)
              valid <- rep(TRUE, n_snp)
              i <- 1
              j <- 2
              while (TRUE){
                  if(chr[i] == chr[j]){
                      mar1 <- pos[i]
                      mar2 <- pos[j]
                      if(mar2 - mar1 <= range){
                          if(missing_count[i] <= missing_count[j]){
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

## Filter out genotype calls meet the criteria.
#' @rdname setCallFilter
setMethod("setCallFilter",
          "GbsrGenotypeData",
          function(object,
                   dp_count,
                   ref_count,
                   alt_count,
                   norm_dp_count,
                   norm_ref_count,
                   norm_alt_count,
                   scan_dp_qtile,
                   scan_ref_qtile,
                   scan_alt_qtile,
                   snp_dp_qtile,
                   snp_ref_qtile,
                   snp_alt_qtile){

              .gds_decomp(object)
              .initFilt(object)

              ## Quantile filtering on each genotype call
              filt_list <- list(dp_count,
                                ref_count,
                                alt_count,
                                norm_dp_count,
                                norm_ref_count,
                                norm_alt_count,
                                scan_dp_qtile,
                                scan_ref_qtile,
                                scan_alt_qtile,
                                snp_dp_qtile,
                                snp_ref_qtile,
                                snp_alt_qtile)

              check <- .makeCallFilter(object, filt_list)

              # Generate filtered AD data
              if(check){
                  .makeCallFilterData(object)
                  object <- .setGenotypeVar(object, "filt.genotype")
              }
              .gds_decomp(object)
              .gds_comp(object)
              return(object)
          })

.initFilt <- function(object){
    valdim <- c(nscan(object, FALSE), nsnp(object, FALSE))
    add.gdsn(.getGdsfmtObj(object),
             "callfilt",
             0L,
             "uint8",
             valdim,
             replace=TRUE)
}

.checkCallArgValid <- function(filt_list){
    if(!{is.numeric(filt_list$dp_count) & length(filt_list$dp_count) == 2 &
            all(filt_list$dp_count >= 0)}){
        stop("dp_count should be two numbers greater than zero.",
             call. = FALSE)
    } else if(filt_list$dp_count[1] > filt_list$dp_count[2]){
        stop("dp_count[1] should be smaller than dp_count[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$ref_count) & length(filt_list$ref_count) == 2 &
            all(filt_list$ref_count >= 0)}){
        stop("ref_count should be two numbers greater than zero.",
             call. = FALSE)
    } else if(filt_list$ref_count[1] > filt_list$ref_count[2]){
        stop("ref_count[1] should be smaller than ref_count[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$alt_count) & length(filt_list$alt_count) == 2 &
            all(filt_list$alt_count >= 0)}){
        stop("alt_count should be two numbers greater than zero.",
             call. = FALSE)
    } else if(filt_list$alt_count[1] > filt_list$alt_count[2]){
        stop("alt_count[1] should be smaller than alt_count[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$norm_dp_count) &
            length(filt_list$norm_dp_count) == 2 &
            all(filt_list$norm_dp_count >= 0)}){
        stop("norm_dp_count should be two numbers greater than zero.",
             call. = FALSE)
    } else if(filt_list$norm_dp_count[1] > filt_list$norm_dp_count[2]){
        stop("norm_dp_count[1] should be smaller than norm_dp_count[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$norm_ref_count) &
            length(filt_list$norm_ref_count) == 2 &
            all(filt_list$norm_ref_count >= 0)}){
        stop("norm_ref_count should be two numbers greater than zero.",
             call. = FALSE)
    } else if(filt_list$norm_ref_count[1] > filt_list$norm_ref_count[2]){
        stop("norm_ref_count[1] should be smaller than norm_ref_count[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$norm_alt_count) &
            length(filt_list$norm_alt_count) == 2 &
            all(filt_list$norm_alt_count >= 0)}){
        stop("norm_alt_count should be two numbers greater than zero.",
             call. = FALSE)
    } else if(filt_list$norm_alt_count[1] > filt_list$norm_alt_count[2]){
        stop("norm_alt_count[1] should be smaller than norm_alt_count[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$scan_dp_qtile) &
            length(filt_list$scan_dp_qtile) == 2 &
            all(filt_list$scan_dp_qtile <= 1) & all(filt_list$scan_dp_qtile >= 0)}){
        stop("scan_dp_qtile should be two numbers in [0-1].",
             call. = FALSE)
    } else if(filt_list$scan_dp_qtile[1] > filt_list$scan_dp_qtile[2]){
        stop("scan_dp_qtile[1] should be smaller than scan_dp_qtile[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$scan_ref_qtile) &
            length(filt_list$scan_ref_qtile) == 2 &
            all(filt_list$scan_ref_qtile <= 1) &
            all(filt_list$scan_ref_qtile >= 0)}){
        stop("scan_ref_qtile should be two numbers in [0-1].",
             call. = FALSE)
    } else if(filt_list$scan_ref_qtile[1] > filt_list$scan_ref_qtile[2]){
        stop("scan_ref_qtile[1] should be smaller than scan_ref_qtile[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$scan_alt_qtile) &
            length(filt_list$scan_alt_qtile) == 2 &
            all(filt_list$scan_alt_qtile <= 1) &
            all(filt_list$scan_alt_qtile >= 0)}){
        stop("scan_alt_qtile should be two numbers in [0-1].",
             call. = FALSE)
    } else if(filt_list$scan_alt_qtile[1] > filt_list$scan_alt_qtile[2]){
        stop("scan_alt_qtile[1] should be smaller than scan_alt_qtile[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$snp_dp_qtile) &
            length(filt_list$snp_dp_qtile) == 2 &
            all(filt_list$snp_dp_qtile <= 1) & all(filt_list$snp_dp_qtile >= 0)}){
        stop("snp_dp_qtile should be two numbers in [0-1].",
             call. = FALSE)
    } else if(filt_list$snp_dp_qtile[1] > filt_list$snp_dp_qtile[2]){
        stop("snp_dp_qtile[1] should be smaller than snp_dp_qtile[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$snp_ref_qtile) &
            length(filt_list$snp_ref_qtile) == 2 &
            all(filt_list$snp_ref_qtile <= 1) & all(filt_list$snp_ref_qtile >= 0)}){
        stop("snp_ref_qtile should be two numbers in [0-1].",
             call. = FALSE)
    } else if(filt_list$snp_ref_qtile[1] > filt_list$snp_ref_qtile[2]){
        stop("snp_ref_qtile[1] should be smaller than snp_ref_qtile[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$snp_alt_qtile) &
            length(filt_list$snp_alt_qtile) == 2 &
            all(filt_list$snp_alt_qtile <= 1) & all(filt_list$snp_alt_qtile >= 0)}){
        stop("snp_alt_qtile should be two numbers in [0-1].",
             call. = FALSE)
    } else if(filt_list$snp_alt_qtile[1] > filt_list$snp_alt_qtile[2]){
        stop("snp_alt_qtile[1] should be smaller than snp_alt_qtile[2].",
             call. = FALSE)
    }
}

.makeCallFilter <- function(object,
                            filt_list){
    check1 <- any(vapply(seq_len(9), FUN.VALUE=logical(1), FUN=function(i){
        if(i <= 6){
            return(any(filt_list[[i]] != c(0, Inf)))
        } else if(i <= 9){
            return(any(filt_list[[i]] != c(0, 1)))
        }
    }))

    if(check1){
        for(i in which(getValidScan(object))){
            .callFilterScan(object, i, filt_list)
        }
        set_filt <- TRUE

    } else {
        set_filt <- FALSE
    }

    check2 <- any(vapply(seq(10, 12), FUN.VALUE=logical(1), FUN=function(i){
        return(any(filt_list[[i]] != c(0, 1)))
    }))

    if(check2){
        .callFilterSnp(object, i, filt_list)
        set_filt <- set_filt | TRUE

    } else {
        set_filt <- set_filt | FALSE
    }
    return(set_filt)
}

.callFilterScan <- function(object, i, filt_list){
    ad_data_node <- .getNodeIndex(object, "annotation/format/AD/data")
    callfilt <- .getNodeIndex(object, "callfilt")
    count_default = c(0, Inf)
    qtile_default = c(0, 1)
    x <- read.gdsn(ad_data_node, start = c(i, 1), count = c(1, -1))

    valid_df <- matrix(TRUE, nsnp(object, FALSE), 9)

    ref <- x[c(TRUE, FALSE)]
    alt <- x[c(FALSE, TRUE)]
    dp <- ref + alt
    invalid_snp <- !getValidSnp(object)
    ref[invalid_snp] <- NA
    alt[invalid_snp] <- NA
    dp[invalid_snp] <- NA

    if(hasFlipped(object)){
        flipped <- getFlipped(object, FALSE)
        tmp <- ref[flipped]
        ref[flipped] <- alt[flipped]
        alt[flipped] <- tmp
    }
    read_list <- list(dp, ref, alt)

    for(j in seq_len(3)){
        valid_df[, j] <- .calcSubFilter(read_list[[j]],
                                        filt_list[[j]],
                                        count_default,
                                        TRUE,
                                        TRUE)
    }

    for(j in 4:6){
        denomi <- sum(x)
        if(denomi == 0){
            valid_df[, j] <- FALSE
        }
        reads <- read_list[[j-3]] / denomi * 10^6
        valid_df[, j] <- .calcSubFilter(reads,
                                        filt_list[[j]],
                                        count_default,
                                        TRUE,
                                        TRUE)
    }

    for(j in 7:9){
        if(any(filt_list[[j]] != qtile_default)){
            threshold <- c(quantile(read_list[[j-6]],
                                    filt_list[[j]][1],
                                    TRUE),
                           quantile(read_list[[j-6]],
                                    filt_list[[j]][2],
                                    TRUE))
            valid_df[, j] <- .calcSubFilter(read_list[[j-6]],
                                            threshold,
                                            c(-1,-1),
                                            TRUE,
                                            TRUE)
        }
    }
    x <- as.numeric(rowSums(valid_df) < 9)

    i_callfilt <- read.gdsn(callfilt, start = c(i, 1), count = c(1, -1))
    i_callfilt <- as.numeric(i_callfilt | x)
    write.gdsn(callfilt, i_callfilt, c(i, 1),  c(1, -1))
}

.callFilterSnp <- function(object, filt_list){
    ad_data_node <- .getNodeIndex(object, "annotation/format/AD/data")
    callfilt <- .getNodeIndex(object, "callfilt")
    qtile_default = c(0, 1)
    for(i in which(getValidSnp(object))){
        x <- read.gdsn(ad_data_node, start = c(1, i*2-1), count = c(-1, 2))
        valid_df <- matrix(TRUE, nscan(object, FALSE), 3)

        ref <- x[, c(TRUE, FALSE)]
        alt <- x[, c(FALSE, TRUE)]
        dp <- ref + alt
        invalid_scan <- !getValidScan(object)
        ref[invalid_scan] <- NA
        alt[invalid_scan] <- NA
        dp[invalid_scan] <- NA

        if(hasFlipped(object)){
            flipped <- getFlipped(object, FALSE)[i]
            tmp <- ref[flipped]
            ref[flipped] <- alt[flipped]
            alt[flipped] <- tmp
        }
        read_list <- list(dp, ref, alt)

        for(j in 10:12){
            if(any(filt_list[[j]] != qtile_default)){
                threshold <- c(quantile(read_list[[j-9]],
                                        filt_list[[j]][1],
                                        TRUE),
                               quantile(read_list[[j-9]],
                                        filt_list[[j]][2],
                                        TRUE))
                valid_df[, j-9] <- .calcSubFilter(read_list[[j-9]],
                                                  threshold,
                                                  c(-1,-1),
                                                  TRUE,
                                                  TRUE)
            }
        }
        x <- as.numeric(rowSums(valid_df) < 3)

        i_callfilt <- read.gdsn(callfilt, start = c(1, i), count = c(-1, 1))
        i_callfilt <- as.numeric(i_callfilt | x)
        write.gdsn(callfilt, i_callfilt, c(1, i),  c(-1, 1))
    }
}

.makeCallFilterData <- function(object){
    ad_node <- .getNodeIndex(object, "annotation/format/AD")
    ad_data <- .getNodeIndex(object, "annotation/format/AD/data")
    callfilt <- .getNodeIndex(object, "callfilt")

    ad_filt_data <- add.gdsn(ad_node,
                             "filt.data",
                             storage="uint32",
                             replace=TRUE)
    assign.gdsn(ad_filt_data, ad_data)

    gt_data <- .getNodeIndex(object, "genotype")

    gt_filt_data <- add.gdsn(.getGdsfmtObj(object),
                             "filt.genotype",
                             storage="bit2",
                             replace=TRUE)
    assign.gdsn(gt_filt_data, gt_data)

    n_scan <- nscan(object, FALSE)
    n_snp <- nsnp(object, FALSE)
    valid_scan_indexes <- which(getValidScan(object))
    i_ad <- i_gt <- NULL
    for (i in valid_scan_indexes){
        i_ad <- read.gdsn(ad_filt_data, c(i, 1), c(1, -1))
        i_callfilt <- read.gdsn(callfilt, c(i, 1), c(1, -1))
        i_ad[rep(i_callfilt, each = 2) == 1] <- 0
        write.gdsn(ad_filt_data, i_ad, c(i, 1), c(1, -1))
        i_ad <- NULL

        i_gt <- read.gdsn(gt_filt_data, c(i, 1), c(1, -1))
        i_gt[i_callfilt == 1] <- 3
        write.gdsn(gt_filt_data, i_gt, c(i, 1), c(1, -1))
        i_gt <- NULL
    }
}

.calcSubFilter <- function(variable, threshold, default, greater, equal){
    if(is.null(variable)){
        return(TRUE)
    }
    if(length(default) == 1){
        if("comp" %in% threshold[1]){
            output <- .compareValues(variable, default, greater, FALSE)

        } else if(threshold != default){
            output <- .compareValues(variable, threshold, greater, equal)

        } else {
            output <- TRUE
        }
    } else if(length(default) == 2){
        if("comp" %in% threshold[1]){
            output <- .compareValues(variable, default[1], greater, FALSE) &
                .compareValues(variable, default[2], !greater, FALSE)

        } else if(any(threshold != default)){
            output <- .compareValues(variable, threshold[1], greater, equal) &
                .compareValues(variable, threshold[2], !greater, equal)

        } else {
            output <- TRUE
        }
    }
    return(output)
}

.compareValues <- function(x, y, greater, equal){
    if(greater & equal){
        out <- x >= y
    } else if(greater & !equal){
        out <- x > y
    } else if(!greater & equal){
        out <- x <= y
    } else if(!greater & !equal){
        out <- x < y
    }
    out[is.na(out)] <- FALSE
    return(out)
}

## Filtering on samples.
#' @rdname setScanFilter
setMethod("setScanFilter",
          "GbsrGenotypeData",
          function(object,
                   id,
                   missing,
                   het,
                   mac,
                   maf,
                   ad_ref,
                   ad_alt,
                   dp,
                   mean_ref,
                   mean_alt,
                   sd_ref,
                   sd_alt){
              filt_list <- list(id = id, missing = missing, het = het, mac = mac,
                                maf = maf, ad_ref = ad_ref, ad_alt = ad_alt,
                                dp = dp, mean_ref = mean_ref, mean_alt = mean_alt,
                                sd_ref = sd_ref, sd_alt = sd_alt)
              update <- .makeStatsFilter(object, filt_list, "scan")
              object <- setValidScan(object, update=update)
              return(object)
          })

## Filtering on markers.
#' @rdname setSnpFilter
setMethod("setSnpFilter",
          "GbsrGenotypeData",
          function(object,
                   id,
                   missing,
                   het,
                   mac,
                   maf,
                   ad_ref,
                   ad_alt,
                   dp,
                   mean_ref,
                   mean_alt,
                   sd_ref,
                   sd_alt){
              filt_list <- list(id = id, missing = missing, het = het, mac = mac,
                                maf = maf, ad_ref = ad_ref, ad_alt = ad_alt,
                                dp = dp, mean_ref = mean_ref, mean_alt = mean_alt,
                                sd_ref = sd_ref, sd_alt = sd_alt)
              update <- .makeStatsFilter(object, filt_list, "snp")
              object <- setValidSnp(object, update=update)
              return(object)
          })

.checkArgValid <- function(filt_list){
    if(!{is.numeric(filt_list$missing) & length(filt_list$missing) == 1 &
            filt_list$missing <= 1 & filt_list$missing >= 0}){
        stop("missing should be a number in [0-1].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$het) & length(filt_list$het) == 2 &
            all(filt_list$het <= 1) & all(filt_list$het >= 0)}){
        stop("het should be two numbers in [0-1].",
             call. = FALSE)
    } else if(filt_list$het[1] > filt_list$het[2]){
        stop("het[1] should be smaller than het[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$mac) &
            length(filt_list$mac) == 1 & filt_list$mac >= 0}){
        stop("mac should be a number greater or equal to zero.",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$maf) & length(filt_list$maf) == 1 &
            filt_list$maf <= 1 & filt_list$maf >= 0}){
        stop("maf should be a number in [0-1].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$ad_ref) &
            length(filt_list$ad_ref) == 2 & all(filt_list$ad_ref >= 0)}){
        stop("ad_ref should be two numbers greater or equal to zero.",
             call. = FALSE)
    } else if(filt_list$ad_ref[1] > filt_list$ad_ref[2]){
        stop("ad_ref[1] should be smaller than ad_ref[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$ad_alt) &
            length(filt_list$ad_alt) == 2 & all(filt_list$ad_alt >= 0)}){
        stop("ad_alt should be two numbers greater or equal to zero.",
             call. = FALSE)
    } else if(filt_list$ad_alt[1] > filt_list$ad_alt[2]){
        stop("ad_alt[1] should be smaller than ad_alt[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$dp) &
            length(filt_list$dp) == 2 & all(filt_list$dp >= 0)}){
        stop("dp should be two numbers greater or equal to zero.",
             call. = FALSE)
    } else if(filt_list$dp[1] > filt_list$dp[2]){
        stop("dp[1] should be smaller than dp[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$mean_ref) &
            length(filt_list$mean_ref) == 2 & all(filt_list$mean_ref >= 0)}){
        stop("mean_ref should be two numbers greater or equal to zero.",
             call. = FALSE)
    } else if(filt_list$mean_ref[1] > filt_list$mean_ref[2]){
        stop("mean_ref[1] should be smaller than mean_ref[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$mean_alt) &
            length(filt_list$mean_alt) == 2 & all(filt_list$mean_alt >= 0)}){
        stop("mean_alt should be two numbers greater or equal to zero.",
             call. = FALSE)
    } else if(filt_list$mean_alt[1] > filt_list$mean_alt[2]){
        stop("mean_alt[1] should be smaller than mean_alt[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$sd_ref) &
            length(filt_list$sd_ref) == 1 & all(filt_list$sd_ref >= 0)}){
        stop("sd_ref should be two numbers greater or equal to zero.",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$sd_alt) &
            length(filt_list$sd_alt) == 1 & all(filt_list$sd_alt >= 0)}){
        stop("sd_alt should be two numbers greater or equal to zero.",
             call. = FALSE)
    }
}

.makeStatsFilter <- function(object, filt_list, target){
    .checkArgValid(filt_list)

    if(target == "scan"){
        if(!is.character(filt_list$id)){
            stop("id should be a character vector.", call. = FALSE)
        }
        if(all(filt_list$id == "")){
            v <- rep(TRUE, nscan(object, TRUE))
        } else {
            v <- !getScanID(object, TRUE) %in% filt_list$id
        }

    } else {
        if(!is.numeric(filt_list$id)){
            stop("id should be a integer vector.", call. = FALSE)
        }
        if(all(is.na(filt_list$id))){
            v <- rep(TRUE, nsnp(object, TRUE))
        } else {
            v <- !getSnpID(object, TRUE) %in% filt_list$id
        }
    }

    # Filtering for samples.
    ## Missing rate
    v <- v & .calcSubFilter(getCountGenoMissing(object, target, TRUE, TRUE),
                            filt_list$missing, 1, FALSE, TRUE)

    ## Heterozygosity
    v <- v & .calcSubFilter(getCountGenoHet(object, target, TRUE, TRUE),
                            filt_list$het, c(0, 1), TRUE, TRUE)

    ## Minor allele count
    v <- v & .calcSubFilter(getMAC(object, target, TRUE),
                            filt_list$mac, 0, TRUE, TRUE)

    ## Minor allele frequency
    v <- v & .calcSubFilter(getMAF(object, target, TRUE),
                            filt_list$maf, 0, TRUE, TRUE)

    ## Reference allele read count
    v <- v & .calcSubFilter(getCountReadRef(object, target, TRUE, FALSE),
                            filt_list$ad_ref, c(0, Inf), TRUE, TRUE)

    ## Alternative allele read count
    v <- v & .calcSubFilter(getCountReadAlt(object, target, TRUE, FALSE),
                            filt_list$ad_alt, c(0, Inf), TRUE, TRUE)

    ## Total read count
    v <- v & .calcSubFilter(getCountRead(object, target, TRUE),
                            filt_list$dp, c(0, Inf), TRUE, TRUE)

    ## Mean reference allele read count
    v <- v & .calcSubFilter(getMeanReadRef(object, target,  TRUE),
                            filt_list$mean_ref, c(0, Inf), TRUE, TRUE)

    ## Mean alternative allele read count
    v <- v & .calcSubFilter(getMeanReadAlt(object, target, TRUE),
                            filt_list$mean_alt, c(0, Inf), TRUE, TRUE)

    ## SD of reference allele read count
    v <- v & .calcSubFilter(getSDReadRef(object, target,  TRUE),
                            filt_list$sd_ref,  Inf,  FALSE, TRUE)

    ## SD of alternative allele read count
    v <- v & .calcSubFilter(getSDReadAlt(object, target, TRUE),
                            filt_list$sd_alt,  Inf,  FALSE,  TRUE)
    return(v)
}

## Filtering on marker based on quality scores.
#' @rdname setInfoFilter
#' @importFrom Biobase pData pData<-
setMethod("setInfoFilter",
          "GbsrGenotypeData",
          function(object,
                   mq,
                   fs,
                   qd,
                   sor,
                   mqranksum,
                   readposranksum,
                   baseqranksum){
              filt_list <- list(mq = mq, fs = fs, qd = qd, sor = sor,
                                mqranksum = mqranksum,
                                readposranksum = readposranksum,
                                baseqranksum = baseqranksum)
              .checkInfoArgValid(filt_list)
              update <- .makeInfoFilter(object, filt_list)
              object <- setValidSnp(object, update=update)
              return(object)
          })

.checkInfoArgValid <- function(filt_list){
    if(!{is.numeric(filt_list$mq) & length(filt_list$mq) == 1 &
            filt_list$mq >= 0}){
        stop("mq should be a number greater than 0.",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$fs) & length(filt_list$fs) == 1 &
            all(filt_list$fs >= 0)}){
        stop("fs should be a number greater than 0.",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$qd) & length(filt_list$qd) == 1 &
            all(filt_list$qd >= 0)}){
        stop("qd should be a number greater than 0.",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$sor) & length(filt_list$sor) == 1 &
            all(filt_list$sor >= 0)}){
        stop("sor should be a number greater than 0.",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$mqranksum) & length(filt_list$mqranksum) == 2}){
        stop("mqranksum should be two numbers.",
             call. = FALSE)
    } else if(filt_list$mqranksum[1] > filt_list$mqranksum[2]){
        stop("mqranksum[1] should be smaller than mqranksum[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$readposranksum) &
            length(filt_list$readposranksum) == 2}){
        stop("readposranksum should be two numbers.",
             call. = FALSE)
    } else if(filt_list$readposranksum[1] > filt_list$readposranksum[2]){
        stop("readposranksum[1] should be smaller than readposranksum[2].",
             call. = FALSE)
    }
    if(!{is.numeric(filt_list$baseqranksum) &
            length(filt_list$baseqranksum) == 2}){
        stop("baseqranksum should be two numbers.",
             call. = FALSE)
    } else if(filt_list$baseqranksum[1] > filt_list$baseqranksum[2]){
        stop("baseqranksum[1] should be smaller than baseqranksum[2].",
             call. = FALSE)
    }
}

.makeInfoFilter <- function(object, filt_list){
    v <- rep(TRUE, nsnp(object))
    # Filtering for samples.
    ## MQ
    v <- v & .calcSubFilter(getInfo(object, "MQ"),
                            filt_list$mq, 0, TRUE, TRUE)

    ## FS
    v <- v & .calcSubFilter(getInfo(object, "FS"),
                            filt_list$fs, Inf, FALSE, TRUE)

    ## QD
    v <- v & .calcSubFilter(getInfo(object, "QD"),
                            filt_list$qd, 0, TRUE, TRUE)

    ## SOR
    v <- v & .calcSubFilter(getInfo(object, "SOR"),
                            filt_list$sor, Inf, FALSE, TRUE)

    ## MQRankSum
    v <- v & .calcSubFilter(getInfo(object, "MQRankSum"),
                            filt_list$mqranksum, c(-Inf, Inf), TRUE, TRUE)

    ## Alternative allele read count
    v <- v & .calcSubFilter(getInfo(object,
                                    "ReadPosRankSum"),
                            filt_list$readposranksum, c(-Inf, Inf), TRUE, TRUE)

    ## Total read count
    v <- v & .calcSubFilter(getInfo(object, "BaseQRankSum"),
                            filt_list$baseqranksum, c(-Inf, Inf), TRUE, TRUE)
    return(v)
}

###############################################################################
## Functions to reset filters.

## Reset filters on samples
#' @rdname resetScanFilters
setMethod("resetScanFilters",
          "GbsrGenotypeData",
          function(object){
              return(setValidScan(object, new=TRUE))
          })

## Reset filters on markers.
#' @rdname resetSnpFilters
setMethod("resetSnpFilters",
          "GbsrGenotypeData",
          function(object){
              return(setValidSnp(object, new=TRUE))
          })

## Reset all filters.
#' @rdname resetFilters
setMethod("resetFilters",
          "GbsrGenotypeData",
          function(object){
              object <- setRawGenotype(object)
              object <- resetScanFilters(object)
              object <- resetSnpFilters(object)
              return(object)
          })

## Set raw data as getnotype, which also means reseting
## filters on genotype calls and read counts made by setCallFilter().
#' @rdname setRawGenotype
setMethod("setRawGenotype",
          "GbsrGenotypeData",
          function(object){
              return(.setGenotypeVar(object, "genotype"))
          })

## Set filtered data as genotype, which also means set again the filter
## which has been made by setCallFilter().
#' @rdname setFiltGenotype
setMethod("setFiltGenotype",
          "GbsrGenotypeData",
          function(object){
              return(.setGenotypeVar(object, "filt.genotype"))
          })

###############################################################################
## Function to create a GDS file of subset data.

## Create a new GDS file with subset data.
#' @importMethodsFrom GWASTools hasSnpVariable
#' @rdname subsetGDS
setMethod("subsetGDS",
          "GbsrGenotypeData",
          function(object,
                   out_fn,
                   snp_incl,
                   scan_incl,
                   incl_parents,
                   verbose){
              if(missing(snp_incl)){
                  snp_incl <- getValidSnp(object)
              } else {
                  stopifnot(nsnp(object, FALSE) == length(snp_incl))
              }

              if(missing(scan_incl)){
                  scan_incl <- getValidScan(object)
                  if(incl_parents){
                      if(hasSnpVariable(object, "parents")){
                          scan_incl[getParents(object, TRUE)] <- TRUE
                      }
                  }
              } else {
                  stopifnot(nscan(object, FALSE) == length(scan_incl))
              }

              n_scan <- nscan(object, FALSE)
              n_snp <- nsnp(object, FALSE)
              no_scan <- n_scan == sum(scan_incl)
              no_snp <- n_snp == sum(snp_incl)
              if(no_scan & no_snp){
                  if(verbose){
                      message("All markers and sampels are valid.")
                      message("Nothing to be subset.")
                  }
              }

              closeGDS(object, FALSE)
              if(!file.copy(.getGDSFileName(object), out_fn)){
                  stop("Failed to create a new file to the following path \n",
                       out_fn, call. = FALSE)
              }
              newgds <- openfn.gds(out_fn, FALSE)
              on.exit({closefn.gds(newgds)})

              all_data_node <- ls.gdsn(newgds,
                                       recursive = TRUE,
                                       include.dirs = FALSE)
              for (i_node in all_data_node){
                  if(grepl("annotation/format", i_node)){
                      if(!grepl("data", i_node)){
                          next
                      }
                  }
                  newgds_i_node <- index.gdsn(newgds, i_node)
                  i_desc <- objdesp.gdsn(newgds_i_node)
                  if(any(i_desc$dim == 0)){
                      next
                  }

                  if(length(i_desc$dim) == 1){
                      check <- sum(which(c(n_scan, n_snp) %in% i_desc$dim))
                      if(check == 1){
                          assign.gdsn(newgds_i_node, seldim=list(which(scan_incl)))

                      } else if(check == 2){
                          assign.gdsn(newgds_i_node, seldim=list(which(snp_incl)))
                      }

                  } else if(length(i_desc$dim) == 2){
                      if(i_desc$name == "parents.genotype"){
                          times <- i_desc$dim[2] / n_snp
                          panrets_index <- seq_len(i_desc$dim[1])
                          seldim <- list(panrets_index,
                                         which(rep(snp_incl, each=times)))
                          assign.gdsn(newgds_i_node, seldim=seldim)

                      } else {
                          times <- i_desc$dim[2] / n_snp
                          seldim <- list(which(scan_incl),
                                         which(rep(snp_incl, each=times)))
                          assign.gdsn(newgds_i_node, seldim=seldim)
                      }
                  }
              }

              .gds_comp(newgds)
              closefn.gds(newgds)
              on.exit({})

              output <- loadGDS(newgds$filename, verbose=verbose)
              if(.getGenotypeVar(object) == "filt.genotype"){
                  output <- setFiltGenotype(output)
              }
              return(output)
          })

###############################################################################
## Function to output a VCF file data stored in the GDS file.

#' @rdname gbsrGDS2VCF
setMethod("gbsrGDS2VCF",
          "GbsrGenotypeData",
          function(object,
                   out_fn,
                   node,
                   incl_parents){
              gds_fn_tmp <- tempfile(fileext = "gds")
              tmp_gds <- subsetGDS(object,
                                   gds_fn_tmp,
                                   incl_parents=incl_parents,
                                   verbose=FALSE)
              on.exit({
                  closefn.gds(.getGdsfmtObj(tmp_gds))
                  unlink(gds_fn_tmp)
              }, TRUE)

              .gds_decomp(tmp_gds)
              .replaceGDSdata(tmp_gds, node)

              out_fn_tmp <- tempfile(fileext = "out.gds")
              closefn.gds(.getGdsfmtObj(tmp_gds))
              on.exit({unlink(gds_fn_tmp)})

              seqSNP2GDS(.getGDSFileName(tmp_gds), out_fn_tmp,
                         "none", FALSE, verbose=FALSE)
              on.exit({unlink(out_fn_tmp)}, TRUE)

              out_gds <- openfn.gds(out_fn_tmp, FALSE)
              on.exit({closefn.gds(out_gds)}, TRUE)
              tmp_gds <- openfn.gds(.getGDSFileName(tmp_gds), FALSE)
              on.exit({closefn.gds(tmp_gds)}, TRUE)

              .insertAnnot(tmp_gds, out_gds)
              .formatAnnot(out_gds)
              .insertHaplotype(tmp_gds, out_gds)
              closefn.gds(out_gds)
              on.exit({
                  closefn.gds(tmp_gds)
                  unlink(gds_fn_tmp)
                  unlink(out_fn_tmp)
              })

              .checkDataLen(out_gds)

              seqGDS2VCF(gdsfile = seqOpen(out_gds$filename),
                         vcf.fn = out_fn,
                         verbose = FALSE)
              return(out_fn)
          })

.checkDataLen <- function(out_gds){
    gds <- openfn.gds(out_gds$filename, readonly = FALSE)
    geno_node <- .getNodeIndex(gds, "genotype/data")
    data_dim <- objdesp.gdsn(geno_node)$dim
    for(i in .getNodeList(gds, "annotation/info")){
        i_gdsn <- .getNodeIndex(gds, i)
        if(objdesp.gdsn(i_gdsn)$dim != data_dim[3]){
            warning("The node ", i, " has invalid data length.")
            delete.gdsn(i_gdsn, TRUE)
        }
    }

    for(i in .getNodeList(gds, "annotation/format")){
        i_gdsn <- .getNodeIndex(gds, i)
        if(grepl("/tmp", i)){
            delete.gdsn(i_gdsn, TRUE)
        } else {
            if(any(objdesp.gdsn(i_gdsn)$dim != data_dim[2:3])){
                warning("The node ", i, " has invalid data length.")
                delete.gdsn(i_gdsn, TRUE)
            }
        }
    }
    closefn.gds(gds)
}

.formatAnnot <- function(out_gds){
    format_node_list <- .getNodeList(out_gds, "annotation/format")
    format_node_list <- grep("/data$", format_node_list, value=TRUE)
    for (gdsn_i in format_node_list){
        node_data <- index.gdsn(out_gds, gdsn_i)
        fld <- index.gdsn(out_gds, sub("/data$", "", gdsn_i))
        tmp <- add.gdsn(fld, "tmp", storage="string32", replace=TRUE)
        apply.gdsn(node_data, 1, as.is="gdsnode", target.node=tmp, FUN=function(x){
            ref <- x[c(TRUE, FALSE)]
            alt <- x[c(FALSE, TRUE)]
            return(paste(ref, alt, sep = ","))
        })
        data_dim <- objdesp.gdsn(index.gdsn(out_gds, "genotype/data"))$dim
        setdim.gdsn(tmp, data_dim[3:2])
        node_data <- add.gdsn(fld, "data", storage="string32", replace=TRUE)
        apply.gdsn(tmp, 1, c, as.is="gdsnode", target.node=node_data)
        setdim.gdsn(node_data, data_dim[2:3])
    }
}

.insertHaplotype <- function(tmp_gds, out_gds){
    if(.existGdsNode(tmp_gds, "estimated.haplotype")){
        hap_fld <- addfolder.gdsn(.getNodeIndex(out_gds, "annotation/format"),
                                  "HAP")
        put.attr.gdsn(hap_fld, "Number", "1")
        put.attr.gdsn(hap_fld, "Type", "Sting")
        put.attr.gdsn(hap_fld, "Description", "Haplotype estimated by GBScleanR.")
        hap_gdsn <- .getNodeIndex(tmp_gds, "estimated.haplotype")
        hap_dim <- objdesp.gdsn(hap_gdsn)$dim
        hap_dim[2] <- hap_dim[2] / 2
        tmp_data <- add.gdsn(tmp_gds, "data", storage="string", replace=TRUE)
        apply.gdsn(hap_gdsn, 1, as.is="gdsnode", target.node=tmp_data,
                   FUN=function(x){
                       x1 <- x[c(TRUE, FALSE)]
                       x2 <- x[c(FALSE, TRUE)]
                       return(paste(x1, x2, sep="|"))
                   })
        setdim.gdsn(tmp_data, hap_dim[2:1])
        hap_data <- add.gdsn(hap_fld, "data", storage="string",
                             valdim=c(hap_dim[1], 0), replace=TRUE)
        apply.gdsn(tmp_data, 1, c, as.is="gdsnode", target.node=hap_data)

        info_gdsn <- .getNodeIndex(out_gds, "annotation/info")
        pgt_node <- add.gdsn(info_gdsn, "PGT", storage = "string", replace=TRUE)
        put.attr.gdsn(pgt_node, "Number", "1")
        put.attr.gdsn(pgt_node, "Type", "String")
        put.attr.gdsn(pgt_node, "Description",
                      "Genotype of each haplotype estimated by GBScleanR.")

        pgt_gdsn <- .getNodeIndex(tmp_gds, "parents.genotype")
        apply.gdsn(pgt_gdsn, 2, as.is="gdsnode", target.node=pgt_node,
                   FUN=function(x){
                       return(paste(x, collapse=","))
                   })
    }
}

.replaceGDSdata <- function(object, node){
    id <- getScanID(object, valid = FALSE)
    add.gdsn(.getGdsfmtObj(object),
             "sample.id",
             id,
             "string",
             replace=TRUE)

    id <- getSnpID(object, valid = FALSE)
    add.gdsn(.getGdsfmtObj(object),
             "snp.id",
             id,
             "string",
             replace=TRUE)

    node <- match.arg(node, c("raw", "filt.genotype", "corrected.genotype"))
    if(.existGdsNode(object, node)){
        gt <- .getNodeIndex(object, "genotype")
        gt_attr <- get.attr.gdsn(gt)
        delete.gdsn(gt)
        newgt <- .getNodeIndex(object, node)
        rename.gdsn(newgt, "genotype")
        put.attr.gdsn(newgt, names(gt_attr)[1], gt_attr[[1]])
    }
}

################################################################################
## Functions to modify GDS file.
setMethod("addScan", "GbsrGenotypeData",
          function(object, id, genotype, reads){
              if(any(id == "")){
                  stop("Please specify non empty sample id.")
              }

              n_snp <- nsnp(object, FALSE)
              if(length(id) == 1){
                  if(!is.numeric(genotype)){
                      stop("genotype can be an integer vector or matrix.")
                  }
                  if(!is.numeric(reads)){
                      stop("reads can be an integer vector or matrix.")
                  }
                  message(length(id), " sample was given to add.")
                  if(length(genotype) != n_snp){
                      stop("The length of genotype should match",
                           " with the number of marekrs.")
                  }
                  if(length(reads) != n_snp * 2){
                      stop("The length of reads should match ",
                           "with twice the number of marekrs.")
                  }
              }

              if(length(id) > 1){
                  if(!is.numeric(genotype)){
                      stop("genotype can be an integer vector or matrix.")
                  }
                  if(!is.numeric(reads)){
                      stop("reads can be an integer vector or matrix.")
                  }
                  message(length(id), " sample was given to add.")
                  if(ncol(genotype) != n_snp){
                      stop("The number of columns of genotype should match",
                           " with the number of marekrs.")
                  }
                  if(ncol(reads) != n_snp * 2){
                      stop("The number of columns of reads should match",
                           " with twice the number of marekrs.")
                  }
                  if(length(id) != nrow(genotype) | length(id) != nrow(reads)){
                      stop("The number of specified IDs and the numbers of ",
                           "rows of genotype and reads should be same.")
                  }
              }

              .gds_decomp(object)
              append.gdsn(.getNodeIndex(object, "sample.id"), id)
              gt_node <- .getNodeIndex(object, "genotype")
              tmp <- read.gdsn(gt_node)
              tmp <- rbind(tmp, genotype)
              gt_attr <- get.attr.gdsn(gt_node)
              gt_node <- add.gdsn(.getGdsfmtObj(object),
                                  "genotype",
                                  tmp,
                                  "bit2",
                                  c(nscan(object) + length(id), nsnp(object)),
                                  "",
                                  replace = TRUE)
              for(i in seq_along(gt_attr)){
                  put.attr.gdsn(gt_node, names(gt_attr)[i], val = gt_attr[[i]])
              }
              ad_node <- .getNodeIndex(object, "annotation/format/AD/data")
              ad_attr <- get.attr.gdsn(ad_node)
              tmp <- read.gdsn(ad_node)
              tmp <- rbind(tmp, reads)
              ad_node <- add.gdsn(.getNodeIndex(object, "annotation/format/AD"),
                                  "data",
                                  tmp,
                                  "float32",
                                  c(nscan(object) + length(id),
                                    nsnp(object) * 2),
                                  "",
                                  replace = TRUE)
              for(i in seq_along(ad_attr)){
                  put.attr.gdsn(ad_node, names(ad_attr)[i], val = ad_attr[[i]])
              }
              .gds_comp(object)
              invisible(TRUE)
          })


###############################################################################
## Functions to handle the GbsrScheme object in the GbsrGenotypeData object.

## Initialize the GbsrScheme object.
#' @rdname initScheme
setMethod("initScheme",
          "GbsrGenotypeData",
          function(object, crosstype, mating){
              parents <- getParents(object)
              scheme <- initScheme(.getSchemeObj(object),
                                   crosstype, mating, parents$memberID)
              object <- .setSchemeObj(object, scheme)
              return(object)
          })

## Add information to the GbsrScheme object.
#' @rdname addScheme
setMethod("addScheme",
          "GbsrGenotypeData",
          function(object, crosstype, mating, pop_size){
              if(missing(pop_size)){
                  pop_size <- NA
              }
              if(missing(mating)){
                  mating <- NA
              }
              scheme <- addScheme(.getSchemeObj(object),
                                  crosstype, mating, pop_size)
              object <- .setSchemeObj(object, scheme)
              return(object)
          })

## Show the information stored in the GbsrScheme object.
#' @rdname showScheme
setMethod("showScheme",
          "GbsrGenotypeData",
          function(object){
              parents <- getParents(object)
              showScheme(.getSchemeObj(object), parents$scanID)
          })
