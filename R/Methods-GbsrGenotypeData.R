###############################################################################
#' @importFrom gdsfmt exist.gdsn apply.gdsn objdesp.gdsn
#' @importFrom SeqArray seqSetFilter seqGetData
.filtData <- function(object, node, filters, reduce = FALSE){
    if(!exist.gdsn(object, node)){
        stop("The GDS node ", node, " does not exists.")
    }

    if(reduce){
        if(any(sapply(filters, sum) == 0)){
            stop('Nothing to return.')
        }
        obj <- objdesp.gdsn(index.gdsn(object, node))
        dim1 <- obj$dim[1]
        out <- apply.gdsn(index.gdsn(object, node), 3, colSums,
                          list(rep(TRUE, dim1), filters$sam, filters$mar),
                          "list")
        out <- do.call("cbind", out)
        out[out == 3*dim1] <- NA

    } else {
        seqSetFilter(object, filters$mar, filters$sam,
                     action = "push+intersect", verbose = FALSE)
        out <- seqGetData(object, node)
        seqSetFilter(object, action = "pop", verbose = FALSE)
        if(length(out) == 0){
            stop('Nothing to return.')
        }
    }
    return(out)
}

#' SeqArray seqGetData
.makeFilter <- function(object, valid = FALSE, parents = FALSE, chr = NULL){
    filt_mar <- rep(TRUE, nmar(object, FALSE))
    filt_sam <- rep(TRUE, nsam(object, FALSE))
    if(valid){
        filt_mar <- validMar(object)
        filt_sam <- validSam(object, parents)

    }
    if(parents == "only"){
        filt_sam <- validSam(object, parents)
    }
    if(!is.null(chr)){
        if(length(chr) > 1){
            warning("More than one values were specified to chr,",
                    " but use the first one")
        }
        filt_mar <- filt_mar & seqGetData(object, "chromosome") == chr[1]
    }
    return(list(mar = filt_mar, sam = filt_sam))
}

.getStatsData <- function(object, var, target, valid){
    target <- match.arg(target, c("marker", "sample"))
    stopifnot(is.logical(valid))
    df <- slot(object, target)
    out <- df[[var]]
    if(valid){
        out <- out[df$valid]
    }
    return(out)
}
###############################################################################
## Internally used functions to control the GDS file.

## Decompress GDS file nodes.
#' @importFrom gdsfmt compression.gdsn ls.gdsn index.gdsn
.gds_decomp <- function(object){
    ls_gdsn <- ls.gdsn(object, TRUE, TRUE, FALSE)
    for(i in ls_gdsn){
        compression.gdsn(index.gdsn(object, i), "")
    }
}

## Compress GDS file nodes.
#' @importFrom gdsfmt readmode.gdsn ls.gdsn index.gdsn compression.gdsn
.gds_comp <- function(object){
    ls_gdsn <- ls.gdsn(object, TRUE, TRUE, FALSE)
    for(i in ls_gdsn){
        i_node <- index.gdsn(object, i)
        compression.gdsn(i_node, "LZMA_RA")
        readmode.gdsn(i_node)
    }
}

###############################################################################
## Exported getter and setter functions which return basic information of the data.

## Get the logical vector indicating valid markers.
#' @importFrom SeqArray seqGetData
#' @rdname validMar
setMethod("validMar",
          "GbsrGenotypeData",
          function(object, chr){
              out <- slot(object, "marker")[["valid"]]
              if(!is.null(chr)){
                  if(length(chr) > 1){
                      warning("More than one values were specified to chr,",
                              " but use the first one")
                  }
                  out <- out & seqGetData(object, "chromosome") == chr[1]
              }
              return(out)
          })

#' @rdname validMar
setMethod("validMar<-",
          "GbsrGenotypeData",
          function(object, value){
              if(length(value) != nmar(object, FALSE)){
                  stop("The length of the given vector should match",
                       " with the number of markers.", call. = FALSE)
              }
              slot(object, "marker")[["valid"]] <- value
              return(object)
          })

## Get the logical vector indicating valid samples.
#' @rdname validSam
setMethod("validSam",
          "GbsrGenotypeData",
          function(object, parents){
              if(parents == "only"){
                  if(is.null(slot(object, "sample")[["parents"]])){
                      stop("No parents info.", call. = FALSE)
                  } else {
                      out <- slot(object, "sample")[["parents"]] != 0
                  }

              } else {
                  out <- slot(object, "sample")[["valid"]]
                  if(parents){
                      if(is.null(slot(object, "sample")[["parents"]])){
                          stop("No parents info.", call. = FALSE)
                      } else {
                          out[slot(object, "sample")[["parents"]]] <- TRUE
                      }
                  }
              }
              return(out)
          })

#' @rdname validSam
setMethod("validSam<-",
          "GbsrGenotypeData",
          function(object, value){
              if(length(value) != nsam(object, FALSE)){
                  stop("The length of the given vector should match",
                       " with the number of samples", call. = FALSE)
              }
              slot(object, "sample")[["valid"]] <- value
              return(object)
          })

## Get the number of SNPs.
#' @rdname nmar
setMethod("nmar",
          "GbsrGenotypeData",
          function(object, valid, chr){
              out <- validMar(object, chr)
              out <- ifelse(valid, sum(out), length(out))
              return(out)
          })

## Get the number of sams (samples).
#' @rdname nsam
setMethod("nsam",
          "GbsrGenotypeData",
          function(object, valid){
              out <- validSam(object)
              out <- ifelse(valid, sum(out), length(out))
              return(out)
          })


## Get the chromosome names in actual strings, indices, or levels.
#' @rdname getChromosome
setMethod("getChromosome",
          "GbsrGenotypeData",
          function(object, valid){
              filters <- .makeFilter(object, valid = valid)
              out <- .filtData(object, "chromosome", filters)
              return(out)
          })

## Get the marker positions.
#' @rdname getPosition
setMethod("getPosition",
          "GbsrGenotypeData",
          function(object, valid, chr){
              filters <- .makeFilter(object, valid = valid, chr = chr)
              out <- .filtData(object, "position", filters)
              return(out)
          })

## Get the reference alleles of markers.
#' @rdname getAllele
setMethod("getAllele",
          "GbsrGenotypeData",
          function(object, valid, chr){
              filters <- .makeFilter(object, valid = valid, chr = chr)
              out <- .filtData(object, "allele", filters)
              return(out)
          })

## Get the marker IDs.
#' @rdname getMarID
setMethod("getMarID",
          "GbsrGenotypeData",
          function(object, valid, chr){
              filters <- .makeFilter(object, valid = valid, chr = chr)
              out <- .filtData(object, "variant.id", filters)
              return(out)
          })

## Get the sample IDs.
#' @rdname getSamID
setMethod("getSamID",
          "GbsrGenotypeData",
          function(object, valid){
              filters <- .makeFilter(object, valid = valid, parents = FALSE)
              out <- .filtData(object, "sample.id", filters)
              return(out)
          })

## Get the values of a variable in the "annotation/info" directory.
#' @rdname getInfo
setMethod("getInfo",
          "GbsrGenotypeData",
          function(object, var, valid, chr){
              node <- paste0("annotation/info/", var)
              filters <- .makeFilter(object, valid = valid, chr = chr)
              out <- .filtData(object, node, filters)
              return(out)
          })

## Get read count data from the annotation/format/AD/data node
## in the GDS file connected to the GbsrGenotypeData object.
#' @rdname getRead
setMethod("getRead",
          "GbsrGenotypeData",
          function(object, node, parents, valid, chr){
              node <- match.arg(node, c("raw", "filt.ad"))
              if(node == "raw"){ node <- "annotation/format/AD" }
              if(node == "filt.ad"){ node <- "annotation/format/FAD" }

              if(!exist.gdsn(object, node)){
                  stop("The GDS node annotation/format/FAD does not exists.\n",
                       "Run setCallFilter() if needed.")
              }

              filters <- .makeFilter(object, parents = parents,
                                     valid = valid, chr = chr)

              out <- .filtData(object, node, filters)
              ref <- out$data[, c(TRUE, FALSE)]
              alt <- out$data[, c(FALSE, TRUE)]
              rownames(ref) <- rownames(alt) <- .filtData(object, "sample.id",
                                                          filters)
              colnames(ref) <- colnames(alt) <- .filtData(object, "variant.id",
                                                          filters)
              return(list(ref = ref, alt = alt))
          })

## Get read count data from one of the nodes `genotype`,
## `corrected.genotype`, or `parents.genotype` in the GDS file
## connected to the GbsrGenotypeData object.
#' @importFrom SeqArray seqGetData
#' @rdname getGenotype
setMethod("getGenotype",
          "GbsrGenotypeData",
          function(object, node, parents, valid, chr){
              node <- match.arg(arg = node,
                                choices =  c("raw", "filt", "cor",
                                             "parents", "ratio", "dosage"))
              if(node == "raw"){node <- "genotype/data"}
              if(node == "filt"){node <- "annotation/format/FGT/data"}
              if(node == "cor"){node <- "annotation/format/CGT/data"}
              if(node == "ratio"){node <- "annotation/format/ARR"}
              if(node == "dosage"){node <- "annotation/format/EDS"}

              if(node == "parents"){
                  node <-  "annotation/info/PGT"
                  if(!exist.gdsn(object, node)){
                      stop("Nothing to return. Run estGeno() to obtain the ",
                           "corrected genotype data.")
                  }

                  filters <- .makeFilter(object, parents = parents,
                                         valid = valid, chr = chr)
                  out <- .filtData(object, node, filters)

                  p_id <- getParents(object)
                  ploidy <- attributes(slot(object, "sample"))$ploidy
                  rownames(out) <- paste(rep(p_id$memberID, each=ploidy),
                                         seq_len(ploidy), sep="_")
                  colnames(out) <- .filtData(object, "variant.id", filters)

              } else {
                  if(!exist.gdsn(object, node)){
                      if(grepl("FGT", node)){
                          stop("Nothing to return.",
                               " Run setCallFilter() if needed.")
                      } else {
                          stop("Nothing to return. Run estGeno() to obtain the",
                               " corrected genotype data.")
                      }
                  }

                  filters <- .makeFilter(object, parents = parents,
                                         valid = valid, chr = chr)

                  if(grepl("ARR|EDS", node)){
                      out <- .filtData(object, node, filters, reduce = FALSE)

                  } else {
                      out <- .filtData(object, node, filters, reduce = TRUE)
                  }

                  rownames(out) <- .filtData(object, "sample.id", filters)
                  colnames(out) <- .filtData(object, "variant.id", filters)
                  dimnames(out) <- list(dimnames(out)$allele,
                                        dimnames(out)$variant)
              }

              return(out)
          })

#' @importFrom gdsfmt readex.gdsn index.gdsn
#' @rdname getHaplotype
setMethod("getHaplotype",
          "GbsrGenotypeData",
          function(object, parents, valid, chr){
              node <-  "annotation/format/HAP/data"
              if(!exist.gdsn(object, node)){
                  stop("Nothing to return. Run estGeno() to obtain the ",
                       "estimated haplotype data.")
              }
              filters <- .makeFilter(object, parents = parents,
                                     valid = valid, chr = chr)
              ploidy <- attributes(slot(object, "sample"))[["ploidy"]]
              out <- readex.gdsn(index.gdsn(object, node),
                                 sel = list(rep(TRUE, ploidy),
                                            filters$sam, filters$mar))

              out[out == 0] <- NA
              return(out)
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
              out <- .getStatsData(object, "countReadRef", target, valid)

              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

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
              out <- .getStatsData(object, "countReadAlt", target, valid)

              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

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
              ref <- .getStatsData(object, "countReadRef", target, valid)
              alt <- .getStatsData(object, "countReadAlt", target, valid)
              if(is.null(ref)){
                  stop("Nothing to rerurn. Run countRead().")
              }
              out <- ref + alt
              return(out)
          })

## Get the number of reference homozygous genotype calls
## per marker and per sample.
#' @rdname getCountGenoRef
setMethod("getCountGenoRef",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              out <- .getStatsData(object, "countGenoRef", target, valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              stopifnot(is.logical(prop))
              if(prop){
                  tmp <- .getStatsData(object, "countGenoMissing", target, valid)
                  if(target == "marker"){
                      nonmissing <- nsam(object, valid) - tmp
                  } else {
                      nonmissing <- nmar(object, valid) - tmp
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
              out <- .getStatsData(object, "countGenoHet", target, valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              stopifnot(is.logical(prop))
              if(prop){
                  tmp <- .getStatsData(object, "countGenoMissing", target, valid)
                  if(target == "marker"){
                      nonmissing <- nsam(object, valid) - tmp
                  } else {
                      nonmissing <- nmar(object, valid) - tmp
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
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countGenoAlt", target, valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              if(prop){
                  tmp <- .getStatsData(object, "countGenoMissing", target, valid)
                  if(target == "marker"){
                      nonmissing <- nsam(object, valid) - tmp
                  } else {
                      nonmissing <- nmar(object, valid) - tmp
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
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countGenoMissing", target, valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              if(prop){
                  if(target == "marker"){
                      out <- out / nsam(object, valid)
                  } else {
                      out <- out / nmar(object, valid)
                  }
              }
              return(out)
          })

## Get the number of reference alleles per marker and per sample.
#' @rdname getCountAlleleRef
setMethod("getCountAlleleRef",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countAlleleRef", target, valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              if(prop){
                  tmp <- .getStatsData(object, "countAlleleMissing", target, valid)
                  if(target == "marker"){
                      nonmissing <- nsam(object, valid) * 2 - tmp
                  } else {
                      nonmissing <- nmar(object, valid) * 2 - tmp
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
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countAlleleAlt", target, valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              if(prop){
                  tmp <- .getStatsData(object, "countAlleleMissing", target, valid)
                  if(target == "marker"){
                      nonmissing <- nsam(object, valid) * 2 - tmp
                  } else {
                      nonmissing <- nmar(object, valid) * 2 - tmp
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
              stopifnot(is.logical(prop))

              out <- .getStatsData(object, "countAlleleMissing", target, valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              if(prop){
                  if(target == "marker"){
                      out <- out / (nsam(object, valid) * 2)
                  } else {
                      out <- out / (nmar(object, valid) * 2)
                  }
              }
              return(out)
          })

## Get the number of mean reference allele reads per marker and per sample.
#' @rdname getMeanReadRef
setMethod("getMeanReadRef",
          "GbsrGenotypeData",
          function(object, target, valid){
              out <- .getStatsData(object, "meanReadRef", target, valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

              return(out)
          })

## Get the number of mean alternative allele reads per marker and per sample.
#' @rdname getMeanReadAlt
setMethod("getMeanReadAlt",
          "GbsrGenotypeData",
          function(object, target, valid){
              out <- .getStatsData(object, "meanReadAlt", target, valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

              return(out)
          })

## Get SD of the number of reference allele reads per marker and per sample.
#' @rdname getSDReadRef
setMethod("getSDReadRef",
          "GbsrGenotypeData",
          function(object, target, valid){
              out <- .getStatsData(object, "sdReadRef", target, valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

              return(out)
          })

## Get SD of the number of alternative allele reads per marker and per sample.
#' @rdname getSDReadAlt
setMethod("getSDReadAlt",
          "GbsrGenotypeData",
          function(object, target, valid){
              out <- .getStatsData(object, "sdReadAlt", target, valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

              return(out)
          })

## Get quantile values of the number of reference allele reads
## per marker and per sample.
#' @rdname getMedianReadRef
setMethod("getMedianReadRef",
          "GbsrGenotypeData",
          function(object, target, q, valid){
              stopifnot(length(q) != 0)
              out <- .getStatsData(object,
                                   paste0("qtileReadRef", q),
                                   target,
                                   valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

              return(out)
          })

## Get quantile values of the number of alternative allele reads
## per marker and per sample.
#' @rdname getMedianReadAlt
setMethod("getMedianReadAlt",
          "GbsrGenotypeData",
          function(object, target, q, valid){
              stopifnot(length(q) != 0)
              out <- .getStatsData(object,
                                   paste0("qtileReadAlt", q),
                                   target,
                                   valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

              return(out)
          })

## Get minor allele frequencies per marker and per sample.
#' @rdname getMAF
setMethod("getMAF",
          "GbsrGenotypeData",
          function(object, target, valid){
              out <- getCountAlleleRef(object, target, valid, TRUE)
              out <- 0.5 - abs(out - 0.5)
              return(out)
          })

## Get minor allele counts per marker and per sample.
#' @rdname getMAC
setMethod("getMAC",
          "GbsrGenotypeData",
          function(object, target, valid){
              ref <- getCountAlleleRef(object, target, valid, FALSE)
              alt <- getCountAlleleAlt(object, target, valid, FALSE)
              out <- ref
              alt_minor <- ref > alt
              out[alt_minor] <- alt[alt_minor]
              return(out)
          })

## Get the information of the parents.
#' @rdname getParents
setMethod("getParents",
          "GbsrGenotypeData",
          function(object, bool){
              if(is.null(slot(object, "sample")[["parents"]])){
                  stop("No parents specified.")
              }

              parents <- slot(object, "sample")[["parents"]]
              if(bool){
                  return(parents != 0)
              }

              p_index <- which(parents != 0)
              p_id <- parents[p_index]
              p_name <- getSamID(object, FALSE)[p_index]
              return(data.frame(sampleID = p_name,
                                memberID = p_id,
                                indexes = p_index))
          })

###############################################################################
## Exported getter functions which return a boolen value.

## Check if the connection to the GDS file is open.
#' @rdname isOpenGDS
#' @importFrom gdsfmt diagnosis.gds
setMethod("isOpenGDS",
          "GbsrGenotypeData",
          function(object){
              tryout <- try(diagnosis.gds(object),
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

## Set parental samples.
#' @rdname setParents
setMethod("setParents",
          "GbsrGenotypeData",
          function(object, parents, mono, bi){
              if(length(parents) == 0 | any(is.na(parents))){
                  stop('Specify valid sample names as parents.', call. = FALSE)
              }
              if(inherits(parents, "numeric")){
                  parents <- getSamID(object)[parents]
              }
              if(!inherits(parents, "character")){
                  stop('Specify valid sample names as parents.', call. = FALSE)
              }

              id <- getSamID(object, FALSE)
              p_index <- match(parents, id)
              if(any(is.na(p_index))){
                  stop("No sample named: ", parents[is.na(p_index)], call. = FALSE)
              }

              n_parents <- length(p_index)
              p_vec <- integer(nsam(object, FALSE))
              for (i in seq_len(n_parents)){
                  p_vec[p_index[i]] <- i
              }

              valid_sam <- validSam(object, FALSE)
              if(!is.null(slot(object, "sample")[["parents"]])){
                  old_parents <- getParents(object, TRUE)
                  valid_sam[old_parents] <- TRUE
              }
              valid_sam[p_vec != 0] <- FALSE
              validSam(object) <- valid_sam
              slot(object, "sample")[["parents"]] <- p_vec

              if(mono | bi){
                  object <- .pGenoFilt(object, mono, bi)
              }
              return(object)
          })

.pGenoFilt <- function(object, mono, bi){
    read <- getRead(object, "raw", "only", FALSE, NULL)
    ref <- read$ref
    alt <- read$alt

    # Find markers which are homozygous in each parent and biallelic
    missing <- ref == 0 & alt == 0
    het <- ref > 0 & alt > 0
    if(mono){
        mono <- !het & !missing
        mono <- colSums(mono) == sum(slot(object, "sample")[["parents"]] != 0)
    } else {
        mono <- rep(TRUE, nmar(object, FALSE))
    }

    if(bi){
        bi <- ref > 0 & alt == 0 & !het & !missing
        bi <- colSums(bi) != sum(slot(object, "sample")[["parents"]] != 0)
    } else {
        bi <- rep(TRUE, nmar(object, FALSE))
    }
    validMar(object) <- mono & bi & validMar(object)
    return(object)
}

###############################################################################
## Show the GbsrGenotypeData object.
#' @importFrom methods show
setMethod("show",
          "GbsrGenotypeData",
          function(object){
              print(object$root)
              message("No of samples: ", nsam(object))
              message("No of markers: ", nmar(object))
              message("Data in the sample slot: ", paste(names(slot(object, "sample")),
                                                         collapse = ", "))
              message("Data in the marker slot: ", paste(names(slot(object, "marker")),
                                                         collapse = ", "))
              message("No of generations in the scheme object: ",
                      length(gds@scheme@mating))
          })

###############################################################################
## Functions to communicate with the GDS file.

## Close the connection to the GDS file.
#' @rdname closeGDS
#' @importFrom SeqArray seqClose seqAddValue
setMethod("closeGDS",
          "GbsrGenotypeData",
          function(object, save_filter, verbose){
              if(save_filter){
                  seqAddValue(object, "annotation/info/GFL", validMar(object),
                              replace = TRUE, compress = "LZMA_RA",
                              verbose = FALSE)
                  seqAddValue(object, "sample.annotation",
                              data.frame(GFL = validSam(object)),
                              replace = TRUE, compress = "LZMA_RA",
                              verbose = FALSE)
              }
              seqClose(object)
              if(verbose){
                  message('The connection to the GDS file was closed.')
              }
          })

## Close the connection to the GDS file.
#' @rdname reopenGDS
#' @importFrom SeqArray seqOpen
setMethod("reopenGDS",
          "GbsrGenotypeData",
          function(object){
              gds_fn <- object$filename
              if(isOpenGDS(object)){
                  closeGDS(object, FALSE)
              }
              object <- new("GbsrGenotypeData",
                            seqOpen(gds_fn, FALSE),
                            marker = slot(object, "marker"),
                            sample = slot(object, "sample"),
                            scheme = slot(object, "scheme"))
              return(object)
          })


###############################################################################
## Functions to calculate statistical summaries of
## genotype calls and read counts.

## Calculate the numbers of each genotype and each allele per marker and
## per sample.
#' @importFrom gdsfmt exist.gdsn apply.gdsn
#' @rdname countGenotype
setMethod("countGenotype",
          "GbsrGenotypeData",
          function(object, target, node){
              node <- match.arg(node,
                                c("raw",
                                  "cor",
                                  "filt"))
              if(node == "raw"){node <- "genotype/data"}
              if(node == "cor"){node <- "annotation/format/CGT/data"}
              if(node == "filt"){node <- "annotation/format/FGT/data"}


              if(!exist.gdsn(object, node)){
                  if(node == "annotation/format/FGT/data"){
                      stop("Nothing to return.",
                           " Run setCallFilter() if needed.")
                  } else {
                      stop("Nothing to return. Run estGeno() to obtain the",
                           " corrected genotype data.")
                  }
              }

              # Counts per sample
              if(target %in% c("both", "sample")){
                  object <- .countGeno(object, node, 2, "sample")
              }

              # Counts per marker
              if(target %in% c("both", "marker")){
                  object <- .countGeno(object, node, 3, "marker")
              }

              return(object)
          })

#' @importFrom gdsfmt apply.gdsn
.countGeno <- function(object, node, margin, slotid){
    ploidy <- attributes(slot(object, "sample"))[["ploidy"]]
    sel <- list(rep(TRUE, ploidy), validSam(object), validMar(object))
    out <- apply.gdsn(index.gdsn(object, node), margin, count_geno, sel, "list")
    out <- do.call("rbind", out)
    if(margin == 2){
        df <- matrix(NA, nsam(object, FALSE), 7)
        df[validSam(object), ] <- out
    } else {
        df <- matrix(NA, nmar(object, FALSE), 7)
        df[validMar(object), ] <- out
    }
    slot(object, slotid)$countGenoRef <- df[, 1]
    slot(object, slotid)$countGenoHet <- df[, 2]
    slot(object, slotid)$countGenoAlt <- df[, 3]
    slot(object, slotid)$countGenoMissing <- df[, 4]
    slot(object, slotid)$countAlleleRef <- df[, 5]
    slot(object, slotid)$countAlleleAlt <- df[, 6]
    slot(object, slotid)$countAlleleMissing <- df[, 7]
    return(object)
}

## Calculate the numbers of reads of each allele per marker and per sample.
#' SeqArray seqSetFilter
#' @rdname countRead
setMethod("countRead",
          "GbsrGenotypeData",
          function(object, target, node, q){
              node <- match.arg(node, c("raw", "filt"))
              if(node == "raw"){node <- "annotation/format/AD"}
              if(node == "filt"){node <- "annotation/format/FAD"}

              if(!exist.gdsn(object, node)){
                  stop("Nothing to return.",
                       " Run setCallFilter() if needed.")
              }

              seqSetFilter(object, validMar(object), validSam(object),
                           action = "push+intersect", verbose = FALSE)

              tot_read <- seqApply(object, node, sum, "by.sample", "integer")

              # Counts per sample
              if(target %in% c("both", "sample")){
                  object <- .countRead(object, node, "by.sample", "sample",
                                       tot_read, q)
              }

              # Counts per marker
              if(target %in% c("both", "marker")){
                  object <- .countRead(object, node, "by.variant", "marker",
                                       tot_read, q)
              }

              seqSetFilter(object, action = "pop", verbose = FALSE)
              return(object)
          })

#' @importFrom SeqArray seqApply
.countRead <- function(object, node, margin, slotid, tot_read, q){
    if(margin == "by.sample"){
        tot_read <- 0
    }
    out <- seqApply(object, node, count_read, margin, "list",
                    tot_read = tot_read)
    out <- do.call("rbind", out)
    if(margin == "by.sample"){
        df <- matrix(NA, nsam(object, FALSE), 8)
        df[validSam(object), ] <- out
    } else {
        df <- matrix(NA, nmar(object, FALSE), 8)
        df[validMar(object), ] <- out
    }
    slot(object, slotid)$countReadRef <- df[, 1]
    slot(object, slotid)$countReadAlt <- df[, 2]
    slot(object, slotid)$meanReadRef <- df[, 3]
    slot(object, slotid)$meanReadAlt <- df[, 4]
    slot(object, slotid)$sdReadRef <- df[, 5]
    slot(object, slotid)$sdReadAlt <- df[, 6]
    slot(object, slotid)$medianReadRef <- df[, 7]
    slot(object, slotid)$medianReadAlt <- df[, 8]
    return(object)
}

###############################################################################
## Filtering functions.

## Thinout markers.
#' @rdname thinMarker
setMethod("thinMarker",
          "GbsrGenotypeData",
          function(object, range){
              if(is.null(slot(object, "marker")[["countGenoMissing"]])){
                  stop('Run countGenotype() first.', call. = FALSE)
              }

              if(!{is.numeric(range) & length(range) == 1 & range >= 0}){
                  stop("range should be a number greater or equal to zero.",
                       call. = FALSE)
              }

              missing_count <- slot(object, "marker")[["countGenoMissing"]]

              chr <- as.integer(factor(getChromosome(object)))
              pos <- getPosition(object)
              out <- thinout_marker(chr, pos, missing_count, range)
              cur_valid <- validMar(object)
              cur_valid[cur_valid] <- out
              validMar(object) <- cur_valid
              return(object)
          })

## Filter out genotype calls meet the criteria.
#' @importFrom gdsfmt put.attr.gdsn
#' @importFrom SeqArray seqOptimize
#' @rdname setCallFilter
setMethod("setCallFilter",
          "GbsrGenotypeData",
          function(object,
                   dp_count,
                   ref_count,
                   alt_count,
                   dp_qtile,
                   ref_qtile,
                   alt_qtile){

              cft_node <- addfolder.gdsn(index.gdsn(object,
                                                    "annotation/format"),
                                         "CFT", replace = TRUE)
              add.gdsn(cft_node, "data", 0L, "uint8",
                       c(nsam(object, FALSE), nmar(object, FALSE)),
                       replace = TRUE)
              put.attr.gdsn(cft_node, "Number", "1")
              put.attr.gdsn(cft_node, "Type", "Integer")
              put.attr.gdsn(cft_node, "Description",
                            "Call filter information generated by GBScleanR")

              ## Quantile filtering on each genotype call
              filt_list <- list(dp_count = dp_count, ref_count = ref_count,
                                alt_count = alt_count, dp_qtile = dp_qtile,
                                ref_qtile = ref_qtile, alt_qtile = alt_qtile)
              .checkCallArgValid(filt_list)
              for(i in which(validSam(object))){
                  .callFilterScan(object, i, filt_list)
              }

              # Generate filtered AD data
              .makeCallFilterData(object)

              closeGDS(object, verbose = FALSE)
              seqOptimize(object$filename, "by.sample",
                          c("FAD", "FGT", "CFT"), verbose = FALSE)
              object <- reopenGDS(object)
              return(object)
          })

.checkCallArgValid <- function(filt_list){
    check <- FALSE
    for(i in names(filt_list)){
        if(grepl("count", i)){
            if(!{is.numeric(filt_list[[i]]) & length(filt_list[[i]]) == 2 &
                    all(filt_list[[i]] >= 0)}){
                stop(i, " should be two numbers greater than zero.",
                     call. = FALSE)
            } else if(filt_list[[i]][1] > filt_list[[i]][2]){
                stop(i, "[1] should be smaller than", i, "[2].", call. = FALSE)
            }
            check <- check | any(filt_list[[i]] != c(0, Inf))
        }
        if(grepl("qtile", i)){
            if(!{is.numeric(filt_list[[i]]) &
                    length(filt_list[[i]]) == 2 &
                    all(filt_list[[i]] <= 1) &
                    all(filt_list[[i]] >= 0)}){
                stop(i, " should be two numbers in [0-1].", call. = FALSE)
            } else if(filt_list[[i]][1] > filt_list[[i]][2]){
                stop(i, "[1] should be smaller than", i, "[2].", call. = FALSE)
            }
            check <- check | any(filt_list[[i]] != c(0, 1))
        }
    }
    if(!check){
        stop("Nothing to filter out. Check the filtering criteria.")
    }
}


#' @importFrom gdsfmt write.gdsn read.gdsn
.callFilterScan <- function(object, i, filt_list){
    ad_data_node <- index.gdsn(object, "annotation/format/AD/data")
    callfilt <- index.gdsn(object, "annotation/format/CFT/data")
    count_default = c(0, Inf)
    qtile_default = c(0, 1)
    x <- read.gdsn(ad_data_node, start = c(i, 1), count = c(1, -1))

    ref <- x[c(TRUE, FALSE)]
    alt <- x[c(FALSE, TRUE)]
    dp <- ref + alt
    invalid_snp <- !validMar(object)
    ref[invalid_snp] <- NA
    alt[invalid_snp] <- NA
    dp[invalid_snp] <- NA

    read_list <- list(dp, ref, alt)

    valid_df <- matrix(TRUE, nmar(object, FALSE), 6)
    for(j in 1:6){
        if(j <= 3){
            valid_df[, j] <- .calcSubFilter(read_list[[j]],
                                            filt_list[[j]],
                                            count_default,
                                            TRUE,
                                            TRUE)
        } else {
            if(any(filt_list[[j]] != qtile_default)){
                threshold <- quantile(read_list[[j-3]][read_list[[j-3]] != 0],
                                      filt_list[[j]][1:2],
                                      TRUE)
                valid_df[, j] <- .calcSubFilter(read_list[[j-3]],
                                                threshold,
                                                c(-1,-1),
                                                TRUE,
                                                TRUE)
            }
        }
    }
    valid_v <- as.numeric(rowSums(valid_df) < 6)
    write.gdsn(callfilt, valid_v, c(i, 1),  c(1, -1))
}

#' @importFrom gdsfmt assign.gdsn index.gdsn addfolder.gdsn add.gdsn objdesp.gdsn read.gdsn
.makeCallFilterData <- function(object){
    ad_data <- index.gdsn(object, "annotation/format/AD/data")
    ad_index <- index.gdsn(object, "annotation/format/AD/@data")
    fad <- addfolder.gdsn(index.gdsn(object, "annotation/format"),
                          "FAD", replace = TRUE)
    put.attr.gdsn(fad, "Number", "2")
    put.attr.gdsn(fad, "Type", "Integer")
    put.attr.gdsn(fad, "Description",
                  "Filtered read counts generated by GBScleanR")
    obj <- objdesp.gdsn(ad_data)
    ad_filt_data <- add.gdsn(fad, "data", valdim = obj$dim, storage = "int64", )
    ad_filt_index <- add.gdsn(fad, "@data", storage = "int32")
    assign.gdsn(ad_filt_index, ad_index)

    gt_node <- index.gdsn(object, "genotype")
    gt_data <-index.gdsn(object, "genotype/data")
    fgt <- addfolder.gdsn(index.gdsn(object, "annotation/format"),
                          "FGT", replace = TRUE)
    obj <- objdesp.gdsn(gt_data)
    gt_filt_data <- add.gdsn(fgt, "data", valdim = obj$dim, storage = "bit2")

    callfilt <- index.gdsn(object, "annotation/format/CFT/data")
    valid_scan_indexes <- which(validSam(object))
    i_ad <- i_gt <- NULL
    for (i in valid_scan_indexes){
        i_ad <- read.gdsn(ad_data, c(i, 1), c(1, -1))
        i_callfilt <- read.gdsn(callfilt, c(i, 1), c(1, -1))
        i_ad[rep(i_callfilt, each = 2) == 1] <- 0
        write.gdsn(ad_filt_data, i_ad, c(i, 1), c(1, -1))
        i_ad <- NULL

        i_gt <- read.gdsn(gt_data, c(1, i, 1), c(-1, 1, -1))
        i_gt[, 1, i_callfilt == 1] <- 3
        write.gdsn(gt_filt_data, i_gt, c(1, i, 1), c(-1, 1, -1))
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
#' @rdname setSamFilter
setMethod("setSamFilter",
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
              filt_list <- list(id = id, missing = missing,
                                het = het, mac = mac,
                                maf = maf, ad_ref = ad_ref, ad_alt = ad_alt,
                                dp = dp, mean_ref = mean_ref,
                                mean_alt = mean_alt,
                                sd_ref = sd_ref, sd_alt = sd_alt)
              out <- .makeStatsFilter(object, filt_list, "sample")
              cur_valid <- validSam(object)
              cur_valid[cur_valid] <- out
              validSam(object) <- cur_valid
              return(object)
          })

## Filtering on markers.
#' @rdname setMarFilter
setMethod("setMarFilter",
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
              out <- .makeStatsFilter(object, filt_list, "marker")
              cur_valid <- validMar(object)
              cur_valid[cur_valid] <- out
              validMar(object) <- cur_valid
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

    if(target == "sample"){
        if(!is.character(filt_list$id)){
            stop("id should be a character vector.", call. = FALSE)
        }
        if(all(filt_list$id == "")){
            v <- rep(TRUE, nsam(object, TRUE))
        } else {
            v <- !getSamID(object) %in% filt_list$id
        }

    } else {
        if(!is.numeric(filt_list$id)){
            stop("id should be a integer vector.", call. = FALSE)
        }
        if(all(is.na(filt_list$id))){
            v <- rep(TRUE, nmar(object, TRUE))
        } else {
            v <- !getMarID(object) %in% filt_list$id
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
              out <- .makeInfoFilter(object, filt_list)
              cur_valid <- validMar(object)
              cur_valid[cur_valid] <- out
              validMar(object) <- cur_valid
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
    v <- rep(TRUE, nmar(object))
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
#' @rdname resetSamFilter
setMethod("resetSamFilter",
          "GbsrGenotypeData",
          function(object){
              slot(object, "sample")[["valid"]] <- TRUE
              return(object)
          })

## Reset filters on markers.
#' @rdname resetMarFilter
setMethod("resetMarFilter",
          "GbsrGenotypeData",
          function(object){
              slot(object, "marker")[["valid"]] <- TRUE
              return(object)
          })

## Reset all filters.
#' @importFrom SeqArray seqDelete
#' @importFrom gdsfmt exist.gdsn
#' @rdname resetCallFilter
setMethod("resetCallFilter",
          "GbsrGenotypeData",
          function(object){
              if(exist.gdsn(object, "annotation/format/CFT/data")){
                  seqDelete(object, fmt.var = c("CFT", "FGT", "FAD"),
                            verbose = FALSE)
              }
              return(object)
          })

## Reset all filters.
#' @rdname resetFilter
setMethod("resetFilter",
          "GbsrGenotypeData",
          function(object){
              object <- resetCallFilter(object)
              object <- resetSamFilter(object)
              object <- resetMarFilter(object)
              return(object)
          })

###############################################################################
## Function to output a VCF file data stored in the GDS file.

#' @rdname gbsrGDS2VCF
#' @importFrom gdsfmt ls.gdsn index.gdsn
#' @importFrom SeqArray seqSetFilter
setMethod("gbsrGDS2VCF",
          "GbsrGenotypeData",
          function(object,
                   out_fn,
                   parents){
              if(is.null(slot(object, "sample")[["parents"]])){
                  parents <- FALSE
                  warnings("No parents info.")
              }
              sam_sel <- validSam(object, parents = parents)
              mar_sel <- validMar(object)

              check <- .checkNodes(object)
              if(any(check)){
                  out_gds <- .modGDS(object, check)
              } else {
                  out_gds <- object
              }

              seqSetFilter(out_gds, mar_sel, sam_sel, action = "set",
                           verbose = FALSE)
              seqGDS2VCF(out_gds, out_fn, verbose = FALSE)
              return(out_fn)
          })

.checkNodes <- function(object){
    check <- c(exist.gdsn(object, "annotation/info/PGT"),
               exist.gdsn(object, "annotation/format/HAP"),
               exist.gdsn(object, "annotation/format/CGT"),
               exist.gdsn(object, "annotation/format/FGT"))
    return(check)
}

.modGDS <- function(object, check){
    if(any(check)){
        closeGDS(object)
        tmp_gds <- tempfile("tmp",tempdir(), ".gds")
        file.copy(object$filename, tmp_gds, overwrite = TRUE)
        tmpgds <- seqOpen(tmp_gds, FALSE)
    }
    if(check[1]){
        data_gdsn <- index.gdsn(tmpgds, "annotation/info/PGT")
        obj <- objdesp.gdsn(data_gdsn)
        if(length(obj$dim) == 2){
            tmp_gdsn <- add.gdsn(index.gdsn(tmpgds, "annotation/info"), "tmp",
                                 storage = "string", replace = TRUE)
            apply.gdsn(data_gdsn, 2, as.is = "gdsnode", target.node = tmp_gdsn,
                       FUN = function(x){
                           x[x == 3] <- "."
                           return(paste(x, collapse = "|"))
                       })
            moveto.gdsn(tmp_gdsn, data_gdsn, "replace+rename")
            fld_gdsn <- index.gdsn(tmpgds, "annotation/info/PGT")
            put.attr.gdsn(fld_gdsn, "Number", "1")
            put.attr.gdsn(fld_gdsn, "Type", "String")
            put.attr.gdsn(fld_gdsn, "Description",
                          "Estimated allele of parental haplotype by GBScleanR")
        }
    }
    if(check[2]){
        data_gdsn <- index.gdsn(tmpgds, "annotation/format/HAP/data")
        obj <- objdesp.gdsn(data_gdsn)
        if(length(obj$dim) == 3){
            tmp_gdsn <- add.gdsn(index.gdsn(tmpgds, "annotation/format/HAP"), "tmp",
                                 storage = "string", replace = TRUE)
            apply.gdsn(data_gdsn, 3, as.is = "gdsnode", target.node = tmp_gdsn,
                       FUN = function(x){
                           x[x == 0] <- "."
                           apply(x, 2, paste, collapse = "|")
                       })
            setdim.gdsn(tmp_gdsn, obj$dim[-1])
            moveto.gdsn(tmp_gdsn, data_gdsn, "replace+rename")
            fld_gdsn <- index.gdsn(tmpgds, "annotation/format/HAP")
            put.attr.gdsn(fld_gdsn, "Number", "1")
            put.attr.gdsn(fld_gdsn, "Type", "String")
            put.attr.gdsn(fld_gdsn, "Description",
                          "Estimated haplotype data by GBScleanR")
        }
    }
    if(check[3]){
        data_gdsn <- index.gdsn(tmpgds, "annotation/format/CGT/data")
        obj <- objdesp.gdsn(data_gdsn)
        if(length(obj$dim) == 3){
            tmp_gdsn <- add.gdsn(index.gdsn(tmpgds, "annotation/format/CGT"), "tmp",
                                 storage = "string", replace = TRUE)
            apply.gdsn(data_gdsn, 3, as.is = "gdsnode", target.node = tmp_gdsn,
                       FUN = function(x){
                           x[x == 3] <- "."
                           apply(x, 2, paste, collapse = "|")
                       })
            setdim.gdsn(tmp_gdsn, obj$dim[-1])
            moveto.gdsn(tmp_gdsn, data_gdsn, "replace+rename")
            fld_gdsn <- index.gdsn(tmpgds, "annotation/format/CGT")
            put.attr.gdsn(fld_gdsn, "Number", "1")
            put.attr.gdsn(fld_gdsn, "Type", "String")
            put.attr.gdsn(fld_gdsn, "Description",
                          "Corrected genotype data by GBScleanR")
        }
    }
    if(check[4]){
        data_gdsn <- index.gdsn(tmpgds, "annotation/format/FGT/data")
        obj <- objdesp.gdsn(data_gdsn)
        if(length(obj$dim) == 3){
            tmp_gdsn <- add.gdsn(index.gdsn(tmpgds, "annotation/format/FGT"), "tmp",
                                 storage = "string", replace = TRUE)
            apply.gdsn(data_gdsn, 3, as.is = "gdsnode", target.node = tmp_gdsn,
                       FUN = function(x){
                           x[x == 3] <- "."
                           apply(x, 2, paste, collapse = "|")
                       })
            setdim.gdsn(tmp_gdsn, obj$dim[-1])
            moveto.gdsn(tmp_gdsn, data_gdsn, "replace+rename")
            fld_gdsn <- index.gdsn(tmpgds, "annotation/format/FGT")
            put.attr.gdsn(fld_gdsn, "Number", "1")
            put.attr.gdsn(fld_gdsn, "Type", "String")
            put.attr.gdsn(fld_gdsn, "Description",
                          "Filtered genotype data by GBScleanR")
        }
    }
    return(tmpgds)
}

###############################################################################
## Function to output a VCF file data stored in the GDS file.

#' @rdname gbsrGDS2CSV
#' @importFrom utils write.table
#'
setMethod("gbsrGDS2CSV",
          "GbsrGenotypeData",
          function(object,
                   out_fn,
                   node,
                   incl_parents,
                   bp2cm,
                   format,
                   read){
              if(is.null(bp2cm)){
                  if(format == "qtl"){
                      bp2cm <- 4e-06
                  } else {
                      bp2cm <- 1
                  }
              }
              if(is.null(getParents(object))){
                  incl_parents <- FALSE
              }
              node <- match.arg(node,
                                c("raw", "filt.genotype", "corrected.genotype"))
              if(format != "qtl" & node == "hap"){
                  geno <- getHaplotype(object, parents = incl_parents)
                  geno <- apply(geno, c(2, 3), paste, collapse = "|")
              } else {
                  geno <- getGenotype(object, node = node,
                                      parents = incl_parents)
              }
              chr <- getChromosome(object)
              pos <- getPosition(object)

              if(format == "qtl"){
                  geno[geno == 2] <- "A"
                  geno[geno == 1] <- "H"
                  geno[geno == 0] <- "B"
                  geno <- rbind(paste(paste0("S", sprintf("%02d", chr)),
                                      pos, sep = "_"), chr, pos * bp2cm, geno)
                  geno <- cbind(rownames(geno), geno)
                  geno[1:3, 1] <- c("id", "", "")
                  write.table(geno, out_fn, quote = FALSE, row.names = FALSE,
                              col.names = FALSE, sep = ",")

              } else {
                  if(read){
                      if(grepl("filt", node)){
                          node <- "filt"
                      } else {
                          node <- "raw"
                      }
                      read <- getRead(object, node = node,
                                      parents = incl_parents)
                      dim_geno <- dim(geno)
                      geno <- vapply(seq_along(geno),
                                     FUN.VALUE = character(1),
                                     function(i){
                                         paste(geno[i],
                                               paste(read$ref[i],
                                                     read$alt[i],
                                                     sep = ","),
                                               sep = ":")
                                     })
                      geno <- matrix(geno, dim_geno[1], dim_geno[2])
                  }
                  geno <- rbind(chr, pos * bp2cm, geno)
                  rownames(geno) <- c("Chr", "Pos", getSamID(object))
                  write.table(geno, out_fn, quote = TRUE, row.names = TRUE,
                              col.names = FALSE, sep = ",")
              }

              return(out_fn)
          })

###############################################################################
## Functions to handle the GbsrScheme object in the GbsrGenotypeData object.

## Initialize the GbsrScheme object.
#' @rdname initScheme
setMethod("initScheme",
          "GbsrGenotypeData",
          function(object, crosstype, mating){
              parents <- getParents(object)
              scheme <- initScheme(slot(object, "scheme"),
                                   crosstype, mating, parents$memberID)
              slot(object, "scheme") <- scheme
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
              scheme <- addScheme(slot(object, "scheme"),
                                  crosstype, mating, pop_size)
              slot(object, "scheme") <- scheme
              return(object)
          })

## Show the information stored in the GbsrScheme object.
#' @rdname showScheme
setMethod("showScheme",
          "GbsrGenotypeData",
          function(object){
              parents <- getParents(object)
              showScheme(slot(object, "scheme"), parents$sampleID)
          })
