###############################################################################
#' @importFrom gdsfmt exist.gdsn apply.gdsn objdesp.gdsn readex.gdsn
.filtData <- function(object, node, filters, reduce = FALSE){
    # Check if the specified node exists
    if(!exist.gdsn(node = object, path = node)){
        stop("The GDS node ", node, " does not exists.", call. = FALSE)
    }
    if(any(vapply(X = filters, FUN = sum, FUN.VALUE = numeric(1)) == 0)){
        stop('Nothing to return.')
    }

    # Prepare filtered output
    target_node <- index.gdsn(node = object, path = node)
    obj <- objdesp.gdsn(target_node)
    dim <- obj$dim
    if(length(dim) == 1){
        if(grepl("sample\\.id", node)){
            sel <- list(filters$sam)

        } else {
            sel <- list(filters$mar)
        }

    } else if(length(dim) == 2){
        sel <- list(filters$sam, filters$mar)

    } else {
        sel <- list(rep(TRUE, dim[1]), filters$sam, filters$mar)
    }
    out <- readex.gdsn(node = target_node, sel = sel)

    if(reduce){
        if(length(dim(out)) == 2){
            stop("'reduce' is applicable when the data has 3 dimentions.")
        }
        # Reduce the dimensions of the output array
        out[out == 3] <- NA
        out <- apply(out, 3, colSums)
    }

    return(out)
}

#' @importFrom SeqArray seqGetData
.makeFilter <- function(object, valid = FALSE, parents = FALSE, chr = NULL){
    # Initialize filters
    filt_mar <- rep(x = TRUE, times = nmar(object = object, valid = FALSE))
    filt_sam <- rep(x = TRUE, times = nsam(object = object, valid = FALSE))

    # Set validity filter
    if(valid){
        filt_mar <- filt_mar & validMar(object = object)
        filt_sam <- filt_sam & validSam(object = object, parents = parents)
    }

    # Set parent filter
    if(parents == "only"){
        filt_sam <- filt_sam & validSam(object = object, parents = parents)
    }

    # Set chromosome filter
    if(!is.null(chr)){
        if(length(chr) > 1){
            warning("More than one values were specified to chr,",
                    " but use the first one")
        }
        chr_filter <- seqGetData(gdsfile = object, var.name = "chromosome") == chr[1]
        filt_mar <- filt_mar & chr_filter
    }
    return(list(mar = filt_mar, sam = filt_sam))
}

.getStatsData <- function(object, var, target, valid){
    # Retrieve sample-wise or marker-wise stats data
    target <- match.arg(arg = target, choices = c("marker", "sample"))
    stopifnot(is.logical(valid))
    df <- slot(object = object, name = target)
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
    ls_gdsn <- ls.gdsn(node = object, include.hidden = TRUE,
                       recursive = TRUE, include.dirs = FALSE)
    for(i in ls_gdsn){
        compression.gdsn(node = index.gdsn(node = object, path = i), compress = "")
    }
}

## Compress GDS file nodes.
#' @importFrom gdsfmt readmode.gdsn ls.gdsn index.gdsn compression.gdsn
.gds_comp <- function(object){
    ls_gdsn <- ls.gdsn(node = object, include.hidden = TRUE,
                       recursive = TRUE, include.dirs = FALSE)
    for(i in ls_gdsn){
        i_node <- index.gdsn(node = object, path = i)
        compression.gdsn(node = i_node, compress = "ZIP_RA")
        readmode.gdsn(node = i_node)
    }
}

## Create a GDS node
#' @importFrom gdsfmt addfolder.gdsn
.create_gdsn <- function(root_node,
                         target_node,
                         new_node,
                         val,
                         storage,
                         valdim,
                         replace,
                         attr){

    check <- exist.gdsn(node = root_node, path = target_node)
    if(!check){ stop(target_node, "does not exist!") }
    target_node_index <- index.gdsn(node = root_node,
                                    path = target_node)
    new_node_index <- addfolder.gdsn(node = target_node_index,
                                     name = new_node,
                                     replace = replace)
    add.gdsn(node = new_node_index, name = "data",
             val = val, storage = storage, valdim = valdim,
             replace = replace)
    .putAttrGdsn(node = new_node_index, attr = attr)
}

.putAttrGdsn <- function(node, attr){
    for(i in seq_along(attr)){
        put.attr.gdsn(node = node,
                      name = names(attr)[i], val = attr[[i]])
    }
}

# Compress GDS nodes
.compressNodes <- function(object, node, decomp = FALSE){
    for(i in seq_along(node)){
        node_index <- index.gdsn(node = object$root, path = node[i])
        if(decomp){
            compress <- ""
        } else {
            compress <- "ZIP_RA"
        }
        compression.gdsn(node = node_index, compress = compress)
        readmode.gdsn(node = node_index)
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
              out <- slot(object = object, name = "marker")[["valid"]]
              if(!is.null(chr)){
                  if(length(chr) > 1){
                      warning("More than one values were specified to chr,",
                              " but use the first one")
                  }
                  chr_data <- seqGetData(gdsfile = object, var.name = "chromosome")
                  out <- out[chr_data == chr[1]]
              }
              return(out)
          })

#' @rdname validMar
#' @importFrom methods slot<-
setMethod("validMar<-",
          "GbsrGenotypeData",
          function(object, value){
              if(length(value) != nmar(object, FALSE)){
                  stop("The length of the given vector should match",
                       " with the number of markers.", call. = FALSE)
              }
              slot(object = object, name = "marker")[["valid"]] <- value
              return(object)
          })

## Get the logical vector indicating valid samples.
#' @rdname validSam
setMethod("validSam",
          "GbsrGenotypeData",
          function(object, parents){
              if(parents == "only"){
                  if(is.null(slot(object = object, name = "sample")[["parents"]])){
                      stop("No parents info.", call. = FALSE)
                  } else {
                      out <- slot(object = object, name = "sample")[["parents"]] != 0
                  }

              } else {
                  out <- slot(object = object, name = "sample")[["valid"]]
                  if(parents){
                      parents <- slot(object = object, name = "sample")[["parents"]]
                      out[parents != 0] <- TRUE
                  }
              }
              return(out)
          })

#' @rdname validSam
#' @importFrom methods slot<-
setMethod("validSam<-",
          "GbsrGenotypeData",
          function(object, value){
              if(length(value) != nsam(object, FALSE)){
                  stop("The length of the given vector should match",
                       " with the number of samples", call. = FALSE)
              }
              slot(object = object, name = "sample")[["valid"]] <- value
              return(object)
          })

## Get the number of markers.
#' @rdname nmar
setMethod("nmar",
          "GbsrGenotypeData",
          function(object, valid, chr){
              out <- validMar(object = object, chr = chr)
              out <- ifelse(test = valid, yes = sum(out), no = length(out))
              return(out)
          })

## Get the number of samples.
#' @rdname nsam
setMethod("nsam",
          "GbsrGenotypeData",
          function(object, valid, parents){
              out <- validSam(object = object, parents = parents)
              out <- ifelse(test = valid, yes = sum(out), no = length(out))
              return(out)
          })


## Get the chromosome names in actual strings, indices, or levels.
#' @rdname getChromosome
setMethod("getChromosome",
          "GbsrGenotypeData",
          function(object, valid){
              filters <- .makeFilter(object = object, valid = valid)
              out <- .filtData(object = object,
                               node = "chromosome",
                               filters = filters)
              return(out)
          })

## Get the marker positions.
#' @rdname getPosition
setMethod("getPosition",
          "GbsrGenotypeData",
          function(object, valid, chr){
              filters <- .makeFilter(object = object, valid = valid, chr = chr)
              out <- .filtData(object = object,
                               node = "position",
                               filters =  filters)
              return(out)
          })

## Get the reference alleles of markers.
#' @rdname getAllele
setMethod("getAllele",
          "GbsrGenotypeData",
          function(object, valid, chr){
              filters <- .makeFilter(object = object, valid = valid, chr = chr)
              out <- .filtData(object = object,
                               node = "allele",
                               filters = filters)
              return(out)
          })

## Get the marker IDs.
#' @rdname getMarID
setMethod("getMarID",
          "GbsrGenotypeData",
          function(object, valid, chr){
              filters <- .makeFilter(object = object, valid = valid, chr = chr)
              out <- .filtData(object = object,
                               node = "variant.id",
                               filters = filters)
              return(out)
          })

## Get the sample IDs.
#' @rdname getSamID
setMethod("getSamID",
          "GbsrGenotypeData",
          function(object, valid, parents){
              filters <- .makeFilter(object = object,
                                     valid = valid,
                                     parents = parents)
              out <- .filtData(object = object,
                               node = "sample.id",
                               filters = filters)
              return(out)
          })

## Get the values of a variable in the "annotation/info" directory.
#' @rdname getInfo
setMethod("getInfo",
          "GbsrGenotypeData",
          function(object, var, valid, chr){
              node <- paste0("annotation/info/", var)
              filters <- .makeFilter(object = object, valid = valid, chr = chr)
              out <- .filtData(object = object, node = node, filters = filters)
              return(out)
          })

## Get read count data from the annotation/format/AD/data node
## in the GDS file connected to the GbsrGenotypeData object.
#' @rdname getRead
#' @importFrom gdsfmt index.gdsn read.gdsn
setMethod("getRead",
          "GbsrGenotypeData",
          function(object, node, parents, valid, chr){
              node <- match.arg(arg = node, choices = c("raw", "filt"))
              if(node == "raw"){ node <- "annotation/format/AD/data" }
              if(node == "filt"){ node <- "annotation/format/FAD/data" }

              if(!exist.gdsn(node = object, path = node)){
                  stop("The GDS node annotation/format/FAD does not exists.\n",
                       "Run setCallFilter() if needed.")
              }

              filters <- .makeFilter(object = object, parents = parents,
                                     valid = valid, chr = chr)

              sample_id <- .filtData(object = object,
                                     node = "sample.id",
                                     filters = filters)
              variant_id <- .filtData(object = object,
                                      node = "variant.id",
                                      filters = filters)
              at_data_node <- sub("data", "@data", node)
              at_data <- read.gdsn(index.gdsn(node = object,
                                              path = at_data_node))
              at_data <- vapply(X = seq_along(at_data),
                                FUN.VALUE = list(1),
                                FUN = function(i){
                                    return(list(rep(i, at_data[i])))
                                })
              at_data <- unlist(at_data)
              filters$mar <- filters$mar[at_data]
              out <- .filtData(object = object, node = node, filters = filters)
              ref <- subset(out, select = c(TRUE, FALSE))
              alt <- subset(out, select = c(FALSE, TRUE))
              rownames(ref) <- rownames(alt) <- sample_id
              colnames(ref) <- colnames(alt) <- variant_id

              ref[is.na(ref)] <- 0
              alt[is.na(alt)] <- 0
              return(list(ref = ref, alt = alt))
          })

## Get read count data from one of the nodes `genotype`,
## `corrected.genotype`, or `parents.genotype` in the GDS file
## connected to the GbsrGenotypeData object.
#' @importFrom SeqArray seqGetData
#' @rdname getGenotype
setMethod("getGenotype",
          "GbsrGenotypeData",
          function(object, node, parents, valid, chr, phased){
              node <- match.arg(arg = node,
                                choices =  c("raw", "filt", "cor",
                                             "parents", "dosage"))

              # If node == "parents", output estimated parental genotypes
              # stored at the PGT node
              if(node == "parents"){
                  out <- .getParentGenotype(object = object,
                                            valid = valid,
                                            chr = chr)
                  return(out)

              }

              # Set FALSE to reduce if phased == TRUE
              if(phased){
                  reduce <- FALSE
                  if(node == "raw"){
                      message("GBScleanR gives no guarantee on",
                              " phasing information for raw genotypes.")
                  }

              } else {
                  if(node == "dosage"){
                      reduce <- FALSE

                  } else {
                      reduce <- TRUE
                  }
              }

              # Set the node from which the data is retrieved
              switch(node,
                     raw = node <- "genotype/data",
                     filt = node <- "annotation/format/FGT/data",
                     cor = node <- "annotation/format/CGT/data",
                     dosage = node <- "annotation/format/EDS/data")

              # Check whether the specified node exists
              if(!exist.gdsn(node = object, path = node)){
                  if(grepl(pattern = "FGT", x = node, fixed = TRUE)){
                      stop("Nothing to return.",
                           "\nRun setCallFilter() if needed.")

                  } else {
                      stop("Nothing to return.\nRun estGeno() to obtain the",
                           " phased corrected genotype data.")
                  }
              }

              # Prepare the output data
              filters <- .makeFilter(object = object, parents = parents,
                                     valid = valid, chr = chr)
              out <- .filtData(object = object, node = node,
                               filters = filters, reduce = reduce)

              # Set numerically expressed missing values to NA
              if(grepl(pattern = "EDS", x = node)){
                  out[out == 63] <- NA

              } else {
                  out[out == 3] <- NA
              }

              # Set row and col names of the output
              sample_id <- .filtData(object = object,
                                     node = "sample.id",
                                     filters = filters)
              variant_id <- .filtData(object = object,
                                      node = "variant.id",
                                      filters = filters)

              if(length(dim(out)) == 2){
                  rownames(out) <- sample_id
                  colnames(out) <- variant_id

              } else {
                  dimnames(out) <- list(phase = seq_len(dim(out)[1]),
                                        sample = sample_id,
                                        marker = variant_id)
              }
              return(out)
          })

.getParentGenotype <- function(object, valid, chr){
    node <-  "annotation/info/PGT"

    # Check whether the node exists
    if(!exist.gdsn(node = object, path = node)){
        stop("Nothing to return. Run estGeno() to obtain the ",
             "corrected genotype data.")
    }

    # Get parental sample information
    p_id <- getParents(object = object, verbose = FALSE)
    if(is.null(p_id)){
        member_id <- c("dummy1", "dummy2")

    } else {
        member_id <- p_id$memberID
    }

    # Get the ploidy information
    sample_ann <- slot(object = object, name = "sample")
    ploidy <- attributes(sample_ann)$ploidy

    # Make filters
    if(valid){
        filt_mar <- validMar(object = object)

    } else {
        filt_mar <- rep(x = TRUE,
                        times = nmar(object = object, valid = valid))
    }

    if(!is.null(chr)){
        if(length(chr) > 1){
            warning("More than one values were specified to chr,",
                    " but use the first one")
        }
        chr_data <- seqGetData(gdsfile = object, var.name = "chromosome")
        filt_mar <- filt_mar & chr_data == chr[1]
    }

    # Prepare the output
    pgt_index <- index.gdsn(node = object, path = node)
    sel <- list(rep(x = TRUE, times = length(unique(member_id)) * ploidy),
                filt_mar)
    out <- readex.gdsn(node = pgt_index, sel = sel)
    p_info <- getParents(object = object)
    p_indices <- rbind(p_info$memberID * 2 - 1, p_info$memberID * 2)
    out <- out[as.vector(p_indices), ]

    # Set row and col names
    rownames(out) <- paste(rep(p_info$sampleID, each = ploidy),
                           seq_len(ploidy), sep = "_")
    filters <- .makeFilter(object = object, parents = TRUE,
                           valid = valid, chr = chr)
    colnames(out) <- .filtData(object = object,
                               node = "variant.id",
                               filters = filters)
    return(out)
}

#' @importFrom gdsfmt readex.gdsn index.gdsn
#' @rdname getHaplotype
setMethod("getHaplotype",
          "GbsrGenotypeData",
          function(object, parents, valid, chr){
              node <-  "annotation/format/HAP/data"

              # Check whether the node exists
              if(!exist.gdsn(node = object, path = node)){
                  stop("Nothing to return.\nRun estGeno() to obtain the ",
                       "estimated haplotype data.")
              }

              # Make filters
              filters <- .makeFilter(object = object, parents = parents,
                                     valid = valid, chr = chr)

              # Get the ploidy information
              sample_ann <- slot(object = object, name = "sample")
              ploidy <- attributes(sample_ann)[["ploidy"]]

              # Prepare the output
              out <- readex.gdsn(index.gdsn(node = object, path = node),
                                 sel = list(rep(x = TRUE, time = ploidy),
                                            filters$sam, filters$mar))

              # Set the numerically expressed missing values to NA
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
              out <- .getStatsData(object = object, var = "countReadRef",
                                   target = target, valid = valid)

              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

              if(prop){
                  tmp <- .getStatsData(object = object, var = "countReadAlt",
                                       target = target, valid = valid)
                  out <- out / (out + tmp)
              }
              return(out)
          })

## Get the number of alternative reads per marker and per sample.
#' @rdname getCountReadAlt
setMethod("getCountReadAlt",
          "GbsrGenotypeData",
          function(object, target, valid, prop){
              out <- .getStatsData(object = object, var = "countReadAlt",
                                   target = target, valid = valid)

              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

              if(prop){
                  tmp <- .getStatsData(object = object, var = "countReadRef",
                                       target = target, valid = valid)
                  out <- out / (out + tmp)
              }
              return(out)
          })

## Get the number of total reads per marker and per sample.
#' @rdname getCountRead
setMethod("getCountRead",
          "GbsrGenotypeData",
          function(object, target, valid){
              ref <- .getStatsData(object = object, var = "countReadRef",
                                   target = target, valid = valid)
              alt <- .getStatsData(object = object, var = "countReadAlt",
                                   target = target, valid = valid)
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
              out <- .getStatsData(object = object, var = "countGenoRef",
                                   target = target, valid = valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              stopifnot(is.logical(prop))
              if(prop){
                  tmp <- .getStatsData(object = object, var = "countGenoMissing",
                                       target = target, valid = valid)
                  if(target == "marker"){
                      nonmissing <- nsam(object = object, valid = valid) - tmp

                  } else {
                      nonmissing <- nmar(object = object, valid = valid) - tmp
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
              out <- .getStatsData(object = object, var = "countGenoHet",
                                   target = target, valid = valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              stopifnot(is.logical(prop))
              if(prop){
                  tmp <- .getStatsData(object = object, var = "countGenoMissing",
                                       target = target, valid = valid)
                  if(target == "marker"){
                      nonmissing <- nsam(object = object, valid = valid) - tmp

                  } else {
                      nonmissing <- nmar(object = object, valid = valid) - tmp
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

              out <- .getStatsData(object = object, var = "countGenoAlt",
                                   target = target, valid = valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              if(prop){
                  tmp <- .getStatsData(object = object, var = "countGenoMissing",
                                       target = target, valid = valid)
                  if(target == "marker"){
                      nonmissing <- nsam(object = object, valid = valid) - tmp

                  } else {
                      nonmissing <- nmar(object = object, valid = valid) - tmp
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

              out <- .getStatsData(object = object, var = "countGenoMissing",
                                   target = target, valid = valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              if(prop){
                  if(target == "marker"){
                      out <- out / nsam(object = object, valid = valid)
                  } else {
                      out <- out / nmar(object = object, valid = valid)
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

              out <- .getStatsData(object = object, var = "countAlleleRef",
                                   target = target, valid = valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              if(prop){
                  tmp <- .getStatsData(object = object, var = "countAlleleMissing",
                                       target = target, valid = valid)
                  if(target == "marker"){
                      nonmissing <- nsam(object = object, valid = valid) * 2 - tmp
                  } else {
                      nonmissing <- nmar(object = object, valid = valid) * 2 - tmp
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

              out <- .getStatsData(object = object, var = "countAlleleAlt",
                                   target = target, valid = valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              if(prop){
                  tmp <- .getStatsData(object = object, var = "countAlleleMissing",
                                       target = target, valid = valid)
                  if(target == "marker"){
                      nonmissing <- nsam(object = object, valid = valid) * 2 - tmp
                  } else {
                      nonmissing <- nmar(object = object, valid = valid) * 2 - tmp
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

              out <- .getStatsData(object = object, var = "countAlleleMissing",
                                   target = target, valid = valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countGenotype().")
              }

              if(prop){
                  if(target == "marker"){
                      out <- out / (nsam(object = object, valid = valid) * 2)
                  } else {
                      out <- out / (nmar(object = object, valid = valid) * 2)
                  }
              }
              return(out)
          })

## Get the number of mean reference allele reads per marker and per sample.
#' @rdname getMeanReadRef
setMethod("getMeanReadRef",
          "GbsrGenotypeData",
          function(object, target, valid){
              out <- .getStatsData(object = object, var = "meanReadRef",
                                   target = target, valid = valid)
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
              out <- .getStatsData(object = object, var = "meanReadAlt",
                                   target = target, valid = valid)
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
              out <- .getStatsData(object = object, var = "sdReadRef",
                                   target = target, valid = valid)
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
              out <- .getStatsData(object = object, var = "sdReadAlt",
                                   target = target, valid = valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

              return(out)
          })


## Get the number of mean reference allele reads per marker and per sample.
#' @rdname getMedianReadRef
setMethod("getMedianReadRef",
          "GbsrGenotypeData",
          function(object, target, valid){
              out <- .getStatsData(object = object, var = "medianReadRef",
                                   target = target, valid = valid)
              if(is.null(out)){
                  stop("Nothing to rerurn. Run countRead().")
              }

              return(out)
          })

## Get the number of mean alternative allele reads per marker and per sample.
#' @rdname getMedianReadAlt
setMethod("getMedianReadAlt",
          "GbsrGenotypeData",
          function(object, target, valid){
              out <- .getStatsData(object = object, var = "medianReadAlt",
                                   target = target, valid = valid)
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
              out <- getCountAlleleRef(object = object, target = target,
                                       valid = valid, prop = TRUE)
              out <- 0.5 - abs(out - 0.5)
              return(out)
          })

## Get minor allele counts per marker and per sample.
#' @rdname getMAC
setMethod("getMAC",
          "GbsrGenotypeData",
          function(object, target, valid){
              ref <- getCountAlleleRef(object = object, target = target,
                                       valid = valid, prop = FALSE)
              alt <- getCountAlleleAlt(object = object, target = target,
                                       valid = valid, prop = FALSE)
              out <- ref
              alt_minor <- ref > alt
              out[alt_minor] <- alt[alt_minor]
              return(out)
          })

## Get the information of the parents.
#' @rdname getParents
setMethod("getParents",
          "GbsrGenotypeData",
          function(object, bool, verbose = TRUE){
              if(is.null(slot(object, "sample")[["parents"]])){
                  if(verbose){
                      warning("No parents specified.", call. = FALSE)
                  }
                  return(NULL)
              }

              parents <- slot(object, "sample")[["parents"]]
              if(bool){
                  return(parents != 0)
              }

              p_index <- which(parents != 0)
              p_id <- parents[p_index]
              p_name <- getSamID(object = object, valid = FALSE)[p_index]
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
              tryout <- try(diagnosis.gds(gds = object),
                            silent = TRUE)
              if(is.list(tryout)){
                  out <- TRUE

              } else if(grepl(pattern = "The GDS file is closed", x = tryout)){
                  out <- FALSE

              } else {
                  out <- tryout[1]
              }
              return(out)
          })

###############################################################################
# Setter for sample annotation information

## Set parental samples.
#' @rdname setParents
#' @importFrom methods slot<-
setMethod("setParents",
          "GbsrGenotypeData",
          function(object, parents, nonmiss, mono, bi){
              if(length(parents) == 0 | any(is.na(parents))){
                  stop('Specify valid sample names as parents.', call. = FALSE)
              }

              if(inherits(x = parents, what = "numeric")){
                  parents <- getSamID(object = object)[parents]
              }

              if(!inherits(x = parents, what = "character")){
                  stop('Specify valid sample names as parents.', call. = FALSE)
              }

              if(!is.null(slot(object = object, name = "sample")[["parents"]])){
                  message("Overwrite the previous parents information.")
                  slot(object = object, name = "sample")[["parents"]] <- NULL
              }

              # Get indices for specified parental samples
              id <- getSamID(object = object, valid = FALSE)
              p_index <- match(x = parents, table = id)
              if(any(is.na(p_index))){
                  missing_id <- parents[is.na(p_index)]
                  stop("No sample named: ", missing_id, call. = FALSE)
              }

              # Check if specified parental samples have been set as replicates
              replicates <- getReplicates(object = object, parents = TRUE)
              p_replicates <- replicates[p_index]
              if(any(duplicated(p_replicates))){
                  stop("Some specified samples have been set as replicates.",
                       "\nSpecify one of the sample IDs for the replicates as",
                       " a parent.",
                       "\nSee the vignette for more information.")
              }

              # Set replicates to have the same index
              n_parents <- length(p_index)
              p_vec <- integer(nsam(object = object, valid = FALSE))
              if(is.null(slot(object = object, name = "sample")[["replicates"]])){
                  p_vec[p_index] <- seq_len(n_parents)

              } else {
                  replicates <- slot(object = object, name = "sample")[["replicates"]]
                  p_rep_id <- replicates[p_index]
                  rep_hit <- match(x = replicates, table = p_rep_id)
                  p_vec[!is.na(rep_hit)] <- na.omit(rep_hit)
              }

              # Set parental samples to be invalid to let GBScleanR ignore them
              # from stats calculations
              valid_sam <- validSam(object = object, parents = FALSE)
              valid_sam[p_vec != 0] <- FALSE
              validSam(object = object) <- valid_sam
              slot(object = object, name = "sample")[["parents"]] <- p_vec
              .create_gdsn(root_node = object$root,
                           target_node = "",
                           new_node = "parents",
                           val = p_vec,
                           storage = "int",
                           valdim = nsam(object = object, valid = FALSE),
                           replace = TRUE, attr = NULL)

              # Parental genotype based marker filtering
              if(mono | bi){
                  object <- .pGenoFilt(object = object,
                                       nonmiss = nonmiss, mono = mono, bi = bi)
              }

              # Update the scheme object
              scheme <- slot(object = object, name = "scheme")
              if(length(slot(object = scheme, name = "parents")) != 0){
                  parents <- getParents(object = object)
                  slot(scheme, "parents") <- parents$memberID
                  slot(object, "scheme") <- scheme
              }

              return(object)
          })

.pGenoFilt <- function(object, nonmiss, mono, bi){
    p_id <- getParents(object = object)

    # Sum up reads if replicates for parental samples exist
    if(length(unique(p_id$memberID)) != length(p_id$memberID)){
        ad <- getRead(object = object, node = "raw", parents = "only",
                      valid = FALSE, chr = NULL)
        rep_id <- getReplicates(object = object, parents = "only")
        ad <- .pileupAD(ad = ad, rep_id = rep_id)
        gt <- .recalcGT(ad = ad)

    } else {
        gt <- getGenotype(object = object, node = "raw", parents = "only",
                          valid = FALSE, chr = NULL)

    }

    if(nonmiss){
        nonmiss <- colSums(is.na(gt)) == 0
    } else {
        nonmiss <- rep(x = TRUE, times = nmar(object = object, valid = FALSE))
    }

    if(mono){
        mono <- colSums(gt == 1) == 0
        mono[is.na(mono)] <- FALSE
    } else {
        mono <- rep(x = TRUE, times = nmar(object = object, valid = FALSE))
    }

    if(bi){
        bi <- colSums(gt == 0) != nrow(gt) & colSums(gt == 2) != nrow(gt)
        bi[is.na(bi)] <- FALSE
    } else {
        bi <- rep(TRUE, nmar(object, FALSE))
    }

    validMar(object = object) <- nonmiss & mono & bi & validMar(object = object)
    return(object)
}

.pileupAD <- function(ad, rep_id){
    ad_ref <- tapply(X = seq_along(rep_id), INDEX = rep_id, FUN = function(i){
        if(length(i) == 1){
            return(ad$ref[i, ])

        } else {
            return(colSums(ad$ref[i, ]))
        }
    })
    ad_ref <- do.call(what = "rbind", args = ad_ref)
    ad_alt <- tapply(X = seq_along(rep_id), INDEX = rep_id, FUN = function(i){
        if(length(i) == 1){
            return(ad$alt[i, ])

        } else {
            return(colSums(ad$alt[i, ]))
        }
    })
    ad_alt <- do.call(what = "rbind", args = ad_alt)
    return(list(ref = ad_ref, alt = ad_alt))
}

.recalcGT <- function(ad = ad, seq_error = 0.0025){
    gt <- matrix(data = NA, nrow = nrow(ad$ref), ncol = ncol(ad$ref))
    pl_ref <- log10(1 - seq_error) * ad$ref + log10(seq_error) * ad$alt
    pl_het <- log10(0.5) * (ad$ref + ad$alt)
    pl_alt <- log10(seq_error) * ad$ref + log10(1 - seq_error) * ad$alt

    gt[pl_ref > pl_het & pl_ref > pl_alt] <- 0
    gt[pl_het > pl_ref & pl_het > pl_alt] <- 1
    gt[pl_alt > pl_het & pl_alt > pl_ref] <- 2
    return(gt)
}

## Set replicates
#' @rdname setReplicates
setMethod("setReplicates",
          "GbsrGenotypeData",
          function(object, replicates){

              n_samples <- nsam(object = object, valid = FALSE)
              if(sum(n_samples) != length(replicates)){
                  stop("\nThe length of the vector specified to the replicates",
                       " argument does not match the number of total samples ",
                       "obtained by nsam(object = object, valid = FALSE).",
                       "\nReplicate IDs should be set to all samples.",
                       "\nSee the vignette for more information.",
                       call. = FALSE)
              }
              replicates <- as.numeric(factor(replicates))
              sample <- slot(object = object, name = "sample")
              sample$replicates <- replicates
              slot(object = object, name = "sample") <- sample

              p_info <- getParents(object = object)
              n_p <- length(unique(p_info$memberID))
              if(!is.null(p_info)){
                  object <- setParents(object = object,
                                       parents = p_info$sampleID)
                  p_info <- getParents(object = object)
                  n_p_updated <- length(unique(p_info$memberID))

                  if(n_p != n_p_updated){
                      warning("The number of parents was changed",
                              " after setting replicates.",
                              "\nCheck the parent information with 'getParents()'.")
                  }
              }

              return(object)
          })

## Get replicates
#' @rdname getReplicates
setMethod("getReplicates",
          "GbsrGenotypeData",
          function(object, parents){
              if(is.null(slot(object = object, name = "sample")[["replicates"]])){
                  out <- getSamID(object = object, valid = FALSE)

              } else {
                  out <- slot(object = object, name = "sample")[["replicates"]]
              }
              out <- out[validSam(object = object, parents = parents)]
              return(out)
          })

###############################################################################
## Show the GbsrGenotypeData object.
#' @importFrom methods show
setMethod("show",
          "GbsrGenotypeData",
          function(object){
              print(object$root)
              message("No of samples: ", nsam(object = object))
              message("No of markers: ", nmar(object = object))
              sample_ann <- slot(object = object, name = "sample")
              message("Data in the sample slot: ", paste(names(sample_ann),
                                                         collapse = ", "))
              marker_ann <- slot(object = object, name = "marker")
              message("Data in the marker slot: ", paste(names(marker_ann),
                                                         collapse = ", "))
              message("No of generations in the scheme object: ",
                      length(object@scheme@mating))
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
                  seqAddValue(gdsfile = object, varnm = "annotation/info/GFL",
                              val = validMar(object = object),
                              replace = TRUE, compress = "ZIP_RA",
                              verbose = FALSE)
                  seqAddValue(gdsfile = object, varnm = "sample.annotation",
                              val = data.frame(GFL = validSam(object = object)),
                              replace = TRUE, compress = "ZIP_RA",
                              verbose = FALSE)
                  gfl_gdsn <- index.gdsn(node = object,
                                         path = "annotation/info/GFL")
                  attr <- list(Number = "1",
                               Type = "Integer",
                               Description = "Filtering information generated by GBScleanR")
                  .putAttrGdsn(node = gfl_gdsn, attr = attr)
              }
              seqClose(object = object)
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
              if(isOpenGDS(object = object)){
                  closeGDS(object = object, save_filter = FALSE)
              }
              object <- new("GbsrGenotypeData",
                            seqOpen(gds.fn = gds_fn, readonly = FALSE),
                            marker = slot(object = object, name = "marker"),
                            sample = slot(object = object, name = "sample"),
                            scheme = slot(object = object, name = "scheme"))
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
              node <- match.arg(arg = node, choices = c("raw", "cor", "filt"))
              switch(node,
                     raw = node <- "genotype/data",
                     cor = node <- "annotation/format/CGT/data",
                     filt = node <- "annotation/format/FGT/data"
              )

              if(!exist.gdsn(node = object, path = node)){
                  if(grepl(pattern = "FGT", x = node)){
                      stop("Nothing to return.",
                           " Run setCallFilter() if needed.")

                  } else {
                      stop("Nothing to return. Run estGeno() to obtain the",
                           " corrected genotype data.")
                  }
              }

              # Counts per sample
              if(target %in% c("both", "sample")){
                  object <- .countGeno(object = object, node = node,
                                       margin = 2, slotid = "sample")
              }

              # Counts per marker
              if(target %in% c("both", "marker")){
                  object <- .countGeno(object = object, node = node,
                                       margin = 3, slotid = "marker")
              }

              # Clean up RAM
              gc()
              return(object)
          })

#' @importFrom gdsfmt apply.gdsn
#' @importFrom methods slot<-
.countGeno <- function(object, node, margin, slotid){

    # Make filters
    sel <- list(rep(x = TRUE, times = 2),
                validSam(object = object),
                validMar(object = object))

    # Apply the Cpp-coded count_geno function
    out <- apply.gdsn(node = index.gdsn(node = object, path = node),
                      margin = margin, FUN = count_geno,
                      selection = sel, as.is = "list")
    out <- do.call(what = "rbind", args = out)

    # Organize output
    if(margin == 2){
        df <- matrix(data = NA,
                     nrow = nsam(object = object, valid = FALSE), ncol = 7)
        df[validSam(object = object), ] <- out

    } else {
        df <- matrix(data = NA, nmar(object = object, valid = FALSE), ncol = 7)
        df[validMar(object = object), ] <- out
    }

    # Put the output into the slot
    slot(object = object, name = slotid)$countGenoRef <- df[, 1]
    slot(object = object, name = slotid)$countGenoHet <- df[, 2]
    slot(object = object, name = slotid)$countGenoAlt <- df[, 3]
    slot(object = object, name = slotid)$countGenoMissing <- df[, 4]
    slot(object = object, name = slotid)$countAlleleRef <- df[, 5]
    slot(object = object, name = slotid)$countAlleleAlt <- df[, 6]
    slot(object = object, name = slotid)$countAlleleMissing <- df[, 7]
    return(object)
}

## Calculate the numbers of reads of each allele per marker and per sample.
#' SeqArray seqSetFilter
#' @rdname countRead
setMethod("countRead",
          "GbsrGenotypeData",
          function(object, target, node){
              node <- match.arg(arg = node, choices = c("raw", "filt"))
              switch(node,
                     raw = node <- "annotation/format/AD",
                     filt = node <- "annotation/format/FAD"
              )

              if(!exist.gdsn(node = object, path = node)){
                  stop("Nothing to return.",
                       " Run setCallFilter() if needed.")
              }

              # Set filtering to specify the markers to be counted
              seqSetFilter(object = object,
                           variant.sel = validMar(object = object),
                           sample.sel = validSam(object),
                           action = "push+intersect", verbose = FALSE)

              # Get max read counts in the data
              tot_read <- seqApply(gdsfile = object, var.name = node,
                                   FUN = sum, margin = "by.sample",
                                   as.is = "integer")

              # Counts per sample
              if(target %in% c("both", "sample")){
                  object <- .countRead(object = object, node = node,
                                       margin = "by.sample", slotid = "sample",
                                       tot_read = tot_read)
              }

              # Counts per marker
              if(target %in% c("both", "marker")){
                  object <- .countRead(object = object, node = node,
                                       margin = "by.variant", slotid = "marker",
                                       tot_read = tot_read)
              }

              # Reset filters
              seqSetFilter(object = object, action = "pop", verbose = FALSE)

              # Clean up RAM
              gc()
              return(object)
          })

#' @importFrom SeqArray seqApply
.countRead <- function(object, node, margin, slotid, tot_read){
    if(margin == "by.sample"){
        tot_read <- 0
    }

    # Apply the Cpp-coded count_read function
    out <- seqApply(gdsfile = object, var.name = node, FUN = count_read,
                    margin = margin, as.is = "list", tot_read = tot_read)
    out <- do.call(what = "rbind", args = out)

    # Organize the output
    if(margin == "by.sample"){
        df <- matrix(data = NA,
                     nrow = nsam(object = object, valid = FALSE), ncol = 8)
        df[validSam(object = object), ] <- out

    } else {
        df <- matrix(data = NA,
                     nrow = nmar(object = object, valid = FALSE), ncol = 8)
        df[validMar(object = object), ] <- out
    }

    # Put the output into the slot
    slot(object = object, name = slotid)$countReadRef <- df[, 1]
    slot(object = object, name = slotid)$countReadAlt <- df[, 2]
    slot(object = object, name = slotid)$meanReadRef <- df[, 3]
    slot(object = object, name = slotid)$meanReadAlt <- df[, 4]
    slot(object = object, name = slotid)$sdReadRef <- df[, 5]
    slot(object = object, name = slotid)$sdReadAlt <- df[, 6]
    slot(object = object, name = slotid)$medianReadRef <- df[, 7]
    slot(object = object, name = slotid)$medianReadAlt <- df[, 8]
    return(object)
}

###############################################################################
## Filtering functions.

## Thinout markers.
#' @rdname thinMarker
setMethod("thinMarker",
          "GbsrGenotypeData",
          function(object, range){

              # Check if the genotype missing rate is available
              marker_ann <- slot(object = object, name = "marker")
              if(is.null(marker_ann[["countGenoMissing"]])){
                  stop('Run countGenotype() first.', call. = FALSE)
              }

              # Check input parameter validity
              if(!{is.numeric(range) & length(range) == 1 & range >= 0}){
                  stop("range should be a number greater or equal to zero.",
                       call. = FALSE)
              }

              # Get missing count par marker
              missing_count <- getCountGenoMissing(object = object)

              # Thinout markers
              chr <- as.integer(factor(getChromosome(object = object)))
              pos <- getPosition(object = object)
              out <- thinout_marker(chr = chr, pos = pos,
                                    missing_count = missing_count, range = range)

              # Update marer validity information
              cur_valid <- validMar(object = object)
              cur_valid[cur_valid] <- out
              validMar(object = object) <- cur_valid
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

              # Create the nodes to store the filtered data
              .create_gdsn(root_node = object$root,
                           target_node = "annotation/format",
                           new_node = "FAD",
                           val = 0L,
                           storage = "int64",
                           valdim = c(nsam(object = object, valid = FALSE),
                                      nmar(object = object, valid = FALSE) * 2),
                           replace = TRUE,
                           attr = list(Number = "2",
                                       Type = "Integer",
                                       Description = "Call-filtered read counts generated by GBScleanR"))
              add.gdsn(node = index.gdsn(node = object$root,
                                         path = "annotation/format/FAD"),
                       name = "@data",
                       val = rep(x = 2,
                                 times = nmar(object = object, valid = FALSE)),
                       storage = "dInt32", compress = "ZIP_RA",
                       closezip = TRUE, visible = FALSE)

              .create_gdsn(root_node = object$root,
                           target_node = "annotation/format",
                           new_node = "FGT",
                           val = 0L,
                           storage = "bit2",
                           valdim = c(2,
                                      nsam(object = object, valid = FALSE),
                                      nmar(object = object, valid = FALSE)),
                           replace = TRUE,
                           attr = list(Number = "1",
                                       Type = "Integer",
                                       Description = "Call-filtered genotype generated by GBScleanR"))

              # Validate the given filtering settings
              filt_list <- list(dp_count = dp_count, ref_count = ref_count,
                                alt_count = alt_count, dp_qtile = dp_qtile,
                                ref_qtile = ref_qtile, alt_qtile = alt_qtile)
              .checkCallArgValid(filt_list = filt_list)

              # Generate filtered AD data
              .applySampleCallFilter(object = object, filt_list = filt_list)

              # Optimize the nodes
              closeGDS(object = object, verbose = FALSE)
              seqOptimize(gdsfn = object$filename, target = "by.sample",
                          format.var = c("FAD", "FGT"),
                          verbose = FALSE)
              object <- reopenGDS(object = object)
              .compressNodes(object = object,
                             node = c("annotation/format/FAD/data",
                                      "annotation/format/FAD/~data",
                                      "annotation/format/FGT/data",
                                      "annotation/format/FGT/~data"))
              return(object)
          })

.checkCallArgValid <- function(filt_list){
    check <- FALSE
    for(i in names(filt_list)){
        if(grepl(pattern = "count", x = i)){
            if(!{is.numeric(filt_list[[i]]) & length(filt_list[[i]]) == 2 &
                    all(filt_list[[i]] >= 0)}){
                stop(i, " should be two numbers greater than zero.",
                     call. = FALSE)

            } else if(filt_list[[i]][1] > filt_list[[i]][2]){
                stop(i, "[1] should be smaller than", i, "[2].", call. = FALSE)
            }
            check <- check | any(filt_list[[i]] != c(0, Inf))
        }

        if(grepl(pattern = "qtile", x = i)){
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

#' @importFrom stats quantile
.applySampleCallFilter <- function(object, filt_list){
    read <- getRead(object = object, node = "raw")
    read$dp <- read$ref + read$alt

    # Record valid calls
    valid_sam <- validSam(object = object)
    valid_mar <- validMar(object = object)
    filter_matrix <- matrix(data = TRUE,
                            nrow = length(valid_sam),
                            ncol = length(valid_mar))

    # Set filtering parameters
    val <- c("dp", "ref", "alt", "dp", "ref", "alt")
    filt <- c("dp_count", "ref_count", "alt_count",
              "dp_qtile", "ref_qtile", "alt_qtile")
    def <- list(c(0, Inf), c(0, Inf), c(0, Inf), c(0, 1), c(0, 1), c(0, 1))

    # Apply filtering
    for(i in seq_along(filt)){
        if(any(def[[i]] != filt_list[[filt[i]]])){
            if(grepl("qtile", filt[i])){
                threshold <- t(vapply(X = seq_len(nrow(read[[val[i]]])),
                                      FUN.VALUE = numeric(length = 2L),
                                      FUN = function(j){
                                          read_j <- read[[val[i]]][j, ]
                                          out <- quantile(x = read_j,
                                                          probs = filt_list[[filt[i]]])
                                          return(out)
                                      }))
            } else {
                threshold <- filt_list[[filt[i]]]
            }
            filter_out <- .calcSubFilter(variable = read[[val[i]]],
                                         threshold = threshold,
                                         default = def[[i]],
                                         greater = TRUE, equal = TRUE)
            filter_matrix[valid_sam, valid_mar] <- filter_matrix[valid_sam, valid_mar] & filter_out
        }
    }

    # Create filtered data
    rm(read)
    gc()
    .makeCallFilterData(object = object, filter_matrix = filter_matrix)
}

#' @importFrom gdsfmt write.gdsn index.gdsn read.gdsn
.makeCallFilterData <- function(object, filter_matrix){
    read <- getRead(object = object, node = "raw", parents = TRUE, valid = FALSE)

    # Create filtered allele read depth data
    fad_data <- index.gdsn(node = object$root,
                           path = "annotation/format/FAD/data")
    val <- matrix(data = 0L, nrow = nrow(read$ref), ncol = ncol(read$ref) * 2)
    read$ref[!filter_matrix] <- 0
    read$alt[!filter_matrix] <- 0
    val[, c(TRUE, FALSE)] <- read$ref
    val[, c(FALSE, TRUE)] <- read$alt
    write.gdsn(node = fad_data, val = val)
    rm(read, val)

    # Create filtered genotype call data
    gt_data <- index.gdsn(node = object$root, path = "genotype/data")
    gt <- read.gdsn(node = gt_data)
    gt[1,,][!filter_matrix] <- 3
    gt[2,,][!filter_matrix] <- 3
    fgt_data <- index.gdsn(node = object$root,
                           path = "annotation/format/FGT/data")
    write.gdsn(node = fgt_data, val = gt)
}

.calcSubFilter <- function(variable, threshold, default, greater, equal){
    if(is.null(variable)){
        return(TRUE)
    }

    if(is.matrix(threshold)){
        # When the input to be evaluated is a matrix,
        # vectorized evaluation will applied
        output <- .compareValues(x = variable, y = threshold[, 1],
                                 greater = greater, equal = equal) &
            .compareValues(x = variable, y = threshold[, 2],
                           greater = !greater, equal = equal)

    } else {
        if(length(default) == 1){
            # If the input to be evaluated has the length of 1,
            # handle it as following.
            if("comp" %in% threshold[1]){
                output <- .compareValues(x = variable, y = default,
                                         greater = greater, equal = FALSE)

            } else if(threshold != default){
                output <- .compareValues(x = variable, y = threshold,
                                         greater = greater, equal = equal)

            } else {
                output <- TRUE
            }

        } else if(length(default) == 2){
            if("comp" %in% threshold[1]){
                output <- .compareValues(x = variable, y = default[1],
                                         greater = greater, equal = FALSE) &
                    .compareValues(x = variable, y = default[2],
                                   greater = !greater, equal = FALSE)

            } else if(any(threshold != default)){
                output <- .compareValues(x = variable, y = threshold[1],
                                         greater = greater, equal = equal) &
                    .compareValues(x = variable, y = threshold[2],
                                   greater = !greater, equal = equal)

            } else {
                output <- TRUE
            }
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
    out[is.na(out)] <- TRUE
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
              filt_list <- list(id = id,
                                missing = missing,
                                het = het,
                                mac = mac,
                                maf = maf,
                                ad_ref = ad_ref,
                                ad_alt = ad_alt,
                                dp = dp,
                                mean_ref = mean_ref,
                                mean_alt = mean_alt,
                                sd_ref = sd_ref,
                                sd_alt = sd_alt)

              # Sample-wise stats based filtering
              out <- .makeStatsFilter(object = object,
                                      filt_list = filt_list,
                                      target = "sample")

              # Update sample filtering information
              cur_valid <- validSam(object = object)
              cur_valid[cur_valid] <- out
              validSam(object = object) <- cur_valid
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
              filt_list <- list(id = id,
                                missing = missing,
                                het = het,
                                mac = mac,
                                maf = maf,
                                ad_ref = ad_ref,
                                ad_alt = ad_alt,
                                dp = dp,
                                mean_ref = mean_ref,
                                mean_alt = mean_alt,
                                sd_ref = sd_ref,
                                sd_alt = sd_alt)

              # Marker-wise stats based filtering
              out <- .makeStatsFilter(object = object,
                                      filt_list = filt_list,
                                      target = "marker")

              # Update marker filtering infromation
              cur_valid <- validMar(object = object)
              cur_valid[cur_valid] <- out
              validMar(object = object) <- cur_valid
              return(object)
          })

.checkArgValid <- function(filt_list, target){
    if(target == "sample"){
        if(!is.character(filt_list$id)){
            stop("id should be a character vector.", call. = FALSE)
        }
    } else {
        if(!is.numeric(filt_list$id)){
            stop("id should be a integer vector.", call. = FALSE)
        }
    }

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

.checkDefault <- function(filt_list){
    return(list(id = !all(filt_list$id %in% c(NA_character_, NA_character_)),
                missing = filt_list$missing != 1,
                het = any(filt_list$het != c(0, 1)),
                mac = filt_list$mac != 0,
                maf = filt_list$maf != 0,
                ad_ref = any(filt_list$ad_ref != c(0, Inf)),
                ad_alt = any(filt_list$ad_alt != c(0, Inf)),
                dp = any(filt_list$dp != c(0, Inf)),
                mean_ref = any(filt_list$mean_ref != c(0, Inf)),
                mean_alt = any(filt_list$mean_alt != c(0, Inf)),
                sd_ref = filt_list$sd_ref != Inf,
                sd_alt = filt_list$sd_alt != Inf))
}

.makeStatsFilter <- function(object, filt_list, target){
    .checkArgValid(filt_list = filt_list, target = target)
    check_list <- .checkDefault(filt_list = filt_list)

    # Initialize the output filter object
    if(target == "sample"){
        if(check_list$id){
            v <- !getSamID(object = object) %in% filt_list$id

        } else {
            v <- rep(TRUE, nsam(object, TRUE))
        }

    } else {
        if(check_list$id){
            v <- !getMarID(object = object) %in% filt_list$id

        } else {
            v <- rep(x = TRUE, times = nmar(object = object, valid = TRUE))
        }
    }

    # Filtering for samples.
    ## Missing rate
    if(check_list$missing){
        v <- v & .calcSubFilter(variable = getCountGenoMissing(object = object,
                                                               target = target,
                                                               valid = TRUE,
                                                               prop = TRUE),
                                threshold = filt_list$missing,
                                default = 1,
                                greater = FALSE,
                                equal = TRUE)
    }

    ## Heterozygosity
    if(check_list$het){
        v <- v & .calcSubFilter(variable = getCountGenoHet(object = object,
                                                           target = target,
                                                           valid = TRUE,
                                                           prop = TRUE),
                                threshold = filt_list$het,
                                default = c(0, 1),
                                greater = TRUE, equal = TRUE)
    }

    ## Minor allele count
    if(check_list$mac){
        v <- v & .calcSubFilter(variable = getMAC(object = object,
                                                  target = target,
                                                  valid = TRUE),
                                threshold = filt_list$mac,
                                default = 0,
                                greater = TRUE,
                                equal = TRUE)
    }

    ## Minor allele frequency
    if(check_list$maf){
        v <- v & .calcSubFilter(variable = getMAF(object = object,
                                                  target = target,
                                                  valid = TRUE),
                                threshold = filt_list$maf,
                                default = 0,
                                greater = TRUE,
                                equal = TRUE)
    }

    ## Reference allele read count
    if(check_list$ad_ref){
        v <- v & .calcSubFilter(variable = getCountReadRef(object = object,
                                                           target = target,
                                                           valid = TRUE,
                                                           prop = FALSE),
                                threshold = filt_list$ad_ref,
                                default = c(0, Inf),
                                greater = TRUE,
                                equal = TRUE)
    }

    ## Alternative allele read count
    if(check_list$ad_alt){
        v <- v & .calcSubFilter(variable = getCountReadAlt(object = object,
                                                           target = target,
                                                           valid = TRUE,
                                                           prop = FALSE),
                                threshold = filt_list$ad_alt,
                                default = c(0, Inf),
                                greater = TRUE,
                                equal = TRUE)
    }

    ## Total read count
    if(check_list$dp){
        v <- v & .calcSubFilter(variable = getCountRead(object = object,
                                                        target = target,
                                                        valid = TRUE),
                                threshold = filt_list$dp,
                                default = c(0, Inf),
                                greater = TRUE,
                                equal = TRUE)
    }

    ## Mean reference allele read count
    if(check_list$mean_ref){
        v <- v & .calcSubFilter(variable = getMeanReadRef(object = object,
                                                          target = target,
                                                          valid = TRUE),
                                threshold = filt_list$mean_ref,
                                default = c(0, Inf),
                                greater = TRUE,
                                equal = TRUE)
    }

    ## Mean alternative allele read count
    if(check_list$mean_alt){
        v <- v & .calcSubFilter(variable = getMeanReadAlt(object = object,
                                                          target = target,
                                                          valid = TRUE),
                                threshold = filt_list$mean_alt,
                                default = c(0, Inf),
                                greater = TRUE,
                                equal = TRUE)
    }

    ## SD of reference allele read count
    if(check_list$sd_ref){
        v <- v & .calcSubFilter(variable = getSDReadRef(object = object,
                                                        target = target,
                                                        valid = TRUE),
                                threshold = filt_list$sd_ref,
                                default = Inf,
                                greater = FALSE,
                                equal = TRUE)
    }

    ## SD of alternative allele read count
    if(check_list$sd_alt){
        v <- v & .calcSubFilter(variable = getSDReadAlt(object = object,
                                                        target = target,
                                                        valid = TRUE),
                                threshold = filt_list$sd_alt,
                                default = Inf,
                                greater = FALSE,
                                equal = TRUE)
    }

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
              filt_list <- list(mq = mq,
                                fs = fs,
                                qd = qd,
                                sor = sor,
                                mqranksum = mqranksum,
                                readposranksum = readposranksum,
                                baseqranksum = baseqranksum)
              .checkInfoArgValid(filt_list = filt_list)

              # Apply filtering on the values recorded in the INFO field
              out <- .makeInfoFilter(object = object, filt_list = filt_list)

              # Update marker filter information
              cur_valid <- validMar(object = object)
              cur_valid[cur_valid] <- out
              validMar(object = object) <- cur_valid
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
    v <- rep(x = TRUE, times = nmar(object = object))
    # Filtering for samples.
    ## MQ
    info_gdsn <- ls.gdsn(node = index.gdsn(node = object,
                                           path = "annotation/info"))
    if(filt_list$mq != 0){
        if("MQ" %in% info_gdsn){
            v <- v & .calcSubFilter(variable = getInfo(object = object, var = "MQ"),
                                    threshold = filt_list$mq,
                                    default = 0,
                                    greater = TRUE,
                                    equal = TRUE)
        } else {
            warning("Filtering on MQ was specified but no MQ in the GDS.")
        }
    }

    ## FS
    if(filt_list$fs != Inf){
        if("FS" %in% info_gdsn){
            v <- v & .calcSubFilter(variable = getInfo(object = object, var = "FS"),
                                    threshold = filt_list$fs,
                                    default = Inf,
                                    greater = FALSE,
                                    equal = TRUE)
        } else {
            warning("Filtering on FS was specified but no FS in the GDS.")
        }
    }

    ## QD
    if(filt_list$qd != 0){
        if("QD" %in% info_gdsn){
            v <- v & .calcSubFilter(variable = getInfo(object = object, var = "QD"),
                                    threshold = filt_list$qd,
                                    default = 0,
                                    greater = TRUE,
                                    equal = TRUE)
        } else {
            warning("Filtering on QD was specified but no QD in the GDS.")
        }
    }

    ## SOR
    if(filt_list$sor != Inf){
        if("SOR" %in% info_gdsn){
            v <- v & .calcSubFilter(variable = getInfo(object = object, var = "SOR"),
                                    threshold = filt_list$sor,
                                    default = Inf,
                                    greater = FALSE,
                                    equal = TRUE)
        } else {
            warning("Filtering on SOR was specified but no SOR in the GDS.")
        }
    }

    ## MQRankSum
    if(any(filt_list$mqranksum != c(-Inf, Inf))){
        if("MQRankSum" %in% info_gdsn){
            v <- v & .calcSubFilter(variable = getInfo(object = object,
                                                       var = "MQRankSum"),
                                    threshold = filt_list$mqranksum,
                                    default = c(-Inf, Inf),
                                    greater = TRUE,
                                    equal = TRUE)
        } else {
            warning("Filtering on MQRankSum was specified but no MQRankSum in the GDS.")
        }
    }

    ## Alternative allele read count
    if(any(filt_list$readposranksum != c(-Inf, Inf))){
        if("ReadPosRankSum" %in% info_gdsn){
            v <- v & .calcSubFilter(variable = getInfo(object = object,
                                                       var = "ReadPosRankSum"),
                                    threshold = filt_list$readposranksum,
                                    default = c(-Inf, Inf),
                                    greater = TRUE,
                                    equal = TRUE)
        } else {
            warning("Filtering on ReadPosRankSum was specified but no ReadPosRankSum in the GDS.")
        }
    }

    ## Total read count
    if(any(filt_list$baseqranksum != c(-Inf, Inf))){
        if("BaseQRankSum" %in% info_gdsn){
            v <- v & .calcSubFilter(variable = getInfo(object = object,
                                                       var = "BaseQRankSum"),
                                    threshold = filt_list$baseqranksum,
                                    default = c(-Inf, Inf),
                                    greater = TRUE,
                                    equal = TRUE)
        } else {
            warning("Filtering on BaseQRankSum was specified but no BaseQRankSum in the GDS.")
        }
    }
    return(v)
}

###############################################################################
## Set dominant markers
#' @rdname setFixedBias
#' @importFrom methods slot<-
setMethod("setFixedBias",
          "GbsrGenotypeData",
          function(object, bias){
              valid_marker <- validMar(object = object)
              out <- rep(NA, length(valid_marker))
              out[valid_marker] <- bias
              slot(object = object, name = "marker")[["dominant"]] <- out
              return(object)
          })

#' @rdname getFixedBias
setMethod("getFixedBias",
          "GbsrGenotypeData",
          function(object, valid, chr){
              out <- slot(object = object, name = "marker")[["dominant"]]
              if(is.null(out)){
                  out <- rep(NA, nmar(object = object, valid = valid, chr = chr))
                  return(out)
              }
              if(valid){
                  out <- out[validMar(object = object)]
              }
              if(!is.null(chr)){
                  out <- out[getChromosome(object = object, valid = valid) %in% chr]
              }
              return(out)
          })

###############################################################################
## Functions to reset filters.

## Reset filters on samples
#' @rdname resetSamFilter
#' @importFrom methods slot<-
setMethod("resetSamFilter",
          "GbsrGenotypeData",
          function(object){
              slot(object = object, name = "sample")[["valid"]] <- TRUE
              slot(object = object, name = "sample")[["parents"]] <- NULL
              return(object)
          })

## Reset filters on markers.
#' @rdname resetMarFilter
#' @importFrom methods slot<-
setMethod("resetMarFilter",
          "GbsrGenotypeData",
          function(object){
              slot(object = object, name = "marker")[["valid"]] <- TRUE
              return(object)
          })

## Reset all filters.
#' @importFrom SeqArray seqDelete
#' @importFrom gdsfmt exist.gdsn
#' @rdname resetCallFilter
setMethod("resetCallFilter",
          "GbsrGenotypeData",
          function(object){
              if(exist.gdsn(node = object, path = "annotation/format/CGT/data")){
                  seqDelete(gdsfile = object, fmt.var = c("FGT", "FAD"),
                            verbose = FALSE)
              }
              return(object)
          })

## Reset all filters.
#' @rdname resetFilter
setMethod("resetFilter",
          "GbsrGenotypeData",
          function(object){
              object <- resetCallFilter(object = object)
              object <- resetSamFilter(object = object)
              object <- resetMarFilter(object = object)
              return(object)
          })

###############################################################################
## Function to output a VCF file data stored in the GDS file.

#' @rdname gbsrGDS2VCF
#' @importFrom gdsfmt ls.gdsn index.gdsn createfn.gds copyto.gdsn closefn.gds
#' @importFrom SeqArray seqSetFilter
setMethod("gbsrGDS2VCF",
          "GbsrGenotypeData",
          function(object,
                   out_fn,
                   node,
                   info.export,
                   fmt.export,
                   parents){

              # Check if the node exists
              node <- match.arg(arg = node, choices = c("raw", "cor"))
              if(node == "cor"){
                  if(!exist.gdsn(node = object, path = "annotation/format/CGT")){
                      stop('`node = "cor"` was specified, but no corrected',
                           ' genotype data.',
                           'Run estGeno() to get corrected genotype data.')
                  }
              }

              # Check if parents have been specified
              if(is.null(slot(object = object, name = "sample")[["parents"]])){
                  parents <- FALSE
                  message("No parents info.")
              }

              # Prepare the temporary GDS file to reorganize data to be output
              check <- .checkNodes(object = object)
              tmp_gds <- tempfile(pattern = "tmp", fileext = ".gds")
              tmpgds <- createfn.gds(filename = tmp_gds, allow.duplicate = TRUE)

              # Check data type because VL_Int causes a fatal error.
              gdsn_ls <- ls.gdsn(node = object$root, include.hidden = TRUE,
                                 recursive = TRUE, include.dirs = FALSE)
              rm_gdsn <- NULL
              for(i in gdsn_ls){
                  input_node <- index.gdsn(node = object$root,
                                           path = i)
                  objdesp <- objdesp.gdsn(node = input_node)
                  if(objdesp$trait == "VL_Int"){
                      rm_gdsn <- c(rm_gdsn, i)
                  }
              }

              if(!is.null(rm_gdsn)){
                  rm_gdsn <- unique(dirname(rm_gdsn))
                  message("VL_Int type data causes a fatal error while",
                          "writing out a VCF file from a GDS in the ",
                          "current implementation.",
                          "\nThus, the following data will be removed ",
                          "from the GDS files.\n",
                          paste(rm_gdsn, collapse = "\n"))
                  for(i in rm_gdsn){
                      input_node <- index.gdsn(node = object$root,
                                               path = i)
                      objdesp <- objdesp.gdsn(node = input_node)
                      delete.gdsn(node = input_node, force = TRUE)
                  }
              }

              gdsn_ls <- ls.gdsn(node = object$root, include.hidden = TRUE)
              for(i in gdsn_ls){
                  copyto.gdsn(node = tmpgds,
                              source = index.gdsn(object$root, i))
              }

              closefn.gds(tmpgds)
              tmpgds <- seqOpen(gds.fn = tmp_gds, readonly = FALSE)

              # Replace genotype call data if node == "cor" was specified
              if(node == "cor"){
                  .replaceGT(object = tmpgds)
                  check[3] <- FALSE
                  fmt.export <- fmt.export[fmt.export != "CGT"]
              }

              # Reformat data to be output
              .modGDS(object = tmpgds, check = check)

              # Make filters
              sam_sel <- validSam(object = object, parents = parents)
              mar_sel <- validMar(object = object)

              # Apply filters
              seqSetFilter(object = tmpgds,
                           variant.sel = mar_sel,
                           sample.sel = sam_sel,
                           action = "push+intersect",
                           verbose = FALSE)

              # Generate the output VCF file
              seqGDS2VCF(gdsfile = tmpgds,
                         vcf.fn = out_fn,
                         info.var = info.export,
                         fmt.var = fmt.export,
                         verbose = FALSE)
              seqClose(object = tmpgds)
              return(out_fn)
          })

.checkNodes <- function(object){
    check <- c(exist.gdsn(node = object, path = "annotation/info/PGT"),
               exist.gdsn(node = object, path = "annotation/format/HAP"),
               exist.gdsn(node = object, path = "annotation/format/CGT"),
               exist.gdsn(node = object, path = "annotation/format/FGT"),
               exist.gdsn(node = object, path = "annotation/info/ADB"),
               exist.gdsn(node = object, path = "annotation/info/MR"))
    return(check)
}

#' @importFrom gdsfmt moveto.gdsn setdim.gdsn append.gdsn
.modGDS <- function(object, check){

    # Reformat PGT data
    if(check[1]){
        data_gdsn <- index.gdsn(node = object, path = "annotation/info/PGT")
        obj <- objdesp.gdsn(node = data_gdsn)

        if(length(obj$dim) == 2){
            tmp_gdsn <- add.gdsn(node = index.gdsn(node = object,
                                                   path = "annotation/info"),
                                 name = "tmp",
                                 storage = "string", replace = TRUE)
            apply.gdsn(node = data_gdsn, margin = 2,
                       as.is = "gdsnode", target.node = tmp_gdsn,
                       FUN = function(x){
                           x[x == 3] <- "."
                           return(paste(x, collapse = "|"))
                       })
            moveto.gdsn(node = tmp_gdsn,
                        loc.node = data_gdsn,
                        relpos = "replace+rename")
            fld_gdsn <- index.gdsn(node = object, path = "annotation/info/PGT")
            attr <- list(Number = "1",
                         Type = "String",
                         Description = "Estimated allele of parental haplotype by GBScleanR")
            .putAttrGdsn(node = fld_gdsn, attr = attr)
        }
    }

    # Reformat HAP data
    if(check[2]){
        data_gdsn <- index.gdsn(node = object, path = "annotation/format/HAP/data")
        obj <- objdesp.gdsn(node = data_gdsn)

        if(length(obj$dim) == 3){
            tmp_gdsn <- add.gdsn(node = index.gdsn(node = object,
                                                   path = "annotation/format/HAP"),
                                 name = "tmp",
                                 storage = "string", replace = TRUE)
            apply.gdsn(node = data_gdsn, margin = 3,
                       as.is = "gdsnode", target.node = tmp_gdsn,
                       FUN = function(x){
                           x[x == 0] <- "."
                           apply(X = x, MARGIN = 2, FUN = paste, collapse = "|")
                       })
            setdim.gdsn(node = tmp_gdsn, valdim = obj$dim[-1])
            moveto.gdsn(node = tmp_gdsn,
                        loc.node = data_gdsn,
                        relpos = "replace+rename")
            fld_gdsn <- index.gdsn(node = object, path = "annotation/format/HAP")
            attr <- list(Number = "1",
                         Type = "String",
                         Description = "Estimated haplotype data by GBScleanR")
            .putAttrGdsn(node = fld_gdsn, attr = attr)
        }
    }

    # Reformat CGT data
    if(check[3]){
        data_gdsn <- index.gdsn(node = object, path = "annotation/format/CGT/data")
        obj <- objdesp.gdsn(node = data_gdsn)
        if(length(obj$dim) == 3){
            tmp_gdsn <- add.gdsn(node = index.gdsn(node = object,
                                                   path = "annotation/format/CGT"),
                                 name = "tmp",
                                 storage = "string", replace = TRUE)
            apply.gdsn(node = data_gdsn, margin = 3,
                       as.is = "gdsnode", target.node = tmp_gdsn,
                       FUN = function(x){
                           x[x == 3] <- "."
                           apply(X = x, MARGIN = 2, FUN = paste, collapse = "|")
                       })
            setdim.gdsn(node = tmp_gdsn, valdim = obj$dim[-1])
            moveto.gdsn(node = tmp_gdsn,
                        loc.node = data_gdsn,
                        relpos = "replace+rename")
            fld_gdsn <- index.gdsn(node = object, path = "annotation/format/CGT")
            attr <- list(Number = "1",
                         Type = "String",
                         Description = "Corrected genotype data by GBScleanR")
            .putAttrGdsn(node = fld_gdsn, attr = attr)
        }
    }

    # Reformat FGT data
    if(check[4]){
        data_gdsn <- index.gdsn(node = object, path = "annotation/format/FGT/data")
        obj <- objdesp.gdsn(node = data_gdsn)
        if(length(obj$dim) == 3){
            tmp_gdsn <- add.gdsn(node = index.gdsn(node = object,
                                                   path = "annotation/format/FGT"),
                                 name = "tmp",
                                 storage = "string", replace = TRUE)
            apply.gdsn(node = data_gdsn, margin = 3,
                       as.is = "gdsnode", target.node = tmp_gdsn,
                       FUN = function(x){
                           x[x == 3] <- "."
                           apply(X = x, MARGIN = 2, FUN = paste, collapse = "/")
                       })
            setdim.gdsn(node = tmp_gdsn, valdim = obj$dim[-1])
            moveto.gdsn(node = tmp_gdsn,
                        loc.node = data_gdsn,
                        relpos = "replace+rename")
            fld_gdsn <- index.gdsn(node = object, path = "annotation/format/FGT")
            attr <- list(Number = "1",
                         Type = "String",
                         Description = "Filtered genotype data by GBScleanR")
            .putAttrGdsn(node = fld_gdsn, attr = attr)
        }
    }

    # Put the header information for ADB
    if(check[5]){
        data_gdsn <- index.gdsn(object, "annotation/info/ADB")
        attr <- list(Number = "1",
                     Type = "Float",
                     Description = "Estimated allele read bias by GBScleanR")
        .putAttrGdsn(node = fld_gdsn, attr = attr)
    }

    # Put the header information for MR
    if(check[6]){
        data_gdsn <- index.gdsn(object, "annotation/info/MR")
        attr <- list(Number = "2",
                     Type = "Float",
                     Description = "Estimated mismapping rate by GBScleanR")
        .putAttrGdsn(node = fld_gdsn, attr = attr)
    }
}

#' @importFrom gdsfmt copyto.gdsn
.replaceGT <- function(object){
    delete.gdsn(node = index.gdsn(node = object, path = "genotype/data"))
    copyto.gdsn(node = index.gdsn(node = object, path = "genotype"),
                source = index.gdsn(node = object,
                                    path = "annotation/format/CGT/data"),
                name = "data")
    delete.gdsn(node = index.gdsn(node = object, path = "genotype/~data"))
    copyto.gdsn(node = index.gdsn(node = object, path = "genotype"),
                source = index.gdsn(node = object,
                                    path = "annotation/format/CGT/~data"),
                name = "~data")
    delete.gdsn(node = index.gdsn(node = object,
                                  path = "annotation/format/CGT/data"))
    delete.gdsn(node = index.gdsn(node = object,
                                  path = "annotation/format/CGT/~data"))
    delete.gdsn(node = index.gdsn(node = object,
                                  path = "annotation/format/CGT"))
}


###############################################################################
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

              # Convert physical marker distances to genetic distances (cM)
              if(is.null(bp2cm)){
                  if(format == "qtl"){
                      bp2cm <- 4e-06

                  } else {
                      bp2cm <- 1
                  }
              }

              # Check if the parents have been specified
              parents <- getParents(object = object)
              if(is.null(parents)){
                  incl_parents <- FALSE
              }
              n_parents <- length(unique(parents$memberID))
              if(n_parents != 2){
                  format <- ""
              }

              node <- match.arg(arg = node,
                                choices = c("raw", "filt", "cor", "hap", "dosage"))

              # Prepare output genotype/haplotype data
              if(node == "hap"){
                  geno <- getHaplotype(object = object, parents = incl_parents)
                  if(format == "qtl"){
                      ploidy <- slot(object = object, name = "sample")$ploidy
                      if(attributes(ploidy == 2)){
                          geno <- apply(X = geno - 1, MARGIN = c(2, 3), FUN = sum)

                      } else {
                          stop("Haplotype data for non-diploid population can",
                               " not be exported in the R/QTL format")
                      }
                  } else {
                      geno <- apply(X = geno, MARGIN = c(2, 3),
                                    FUN = paste, collapse = "|")
                  }
              } else {
                  geno <- getGenotype(object = object, node = node,
                                      parents = incl_parents)
              }

              # Prepare marker information
              chr <- getChromosome(object = object)
              pos <- getPosition(object = object)
              id <- getSamID(object = object, valid = FALSE)
              if(incl_parents){
                  id <- id[validSam(object = object, parents = TRUE)]
              } else {
                  id <- id[validSam(object = object, parents = FALSE)]
              }

              # Convert genotype data to the r/qtl format
              if(format == "qtl"){
                  geno[geno == 0] <- "A"
                  geno[geno == 1] <- "H"
                  geno[geno == 2] <- "B"
                  geno <- rbind(paste(chr, pos, sep = "_"),
                                chr, pos * bp2cm, geno)
                  geno <- cbind(c("id", "", "", id), geno)
                  write.table(x = geno, file = out_fn, quote = FALSE,
                              row.names = FALSE, col.names = FALSE, sep = ",")

              } else {
                  # Prepare read count data
                  if(read){
                      if(grepl(pattern = "filt", x = node)){
                          node <- "filt"
                      } else {
                          node <- "raw"
                      }
                      read <- getRead(object = object, node = node,
                                      parents = incl_parents)
                      dim_geno <- dim(geno)

                      # Combine genotype and read information
                      geno <- vapply(X = seq_along(geno),
                                     FUN.VALUE = character(1),
                                     FUN = function(i){
                                         paste(geno[i],
                                               paste(read$ref[i],
                                                     read$alt[i],
                                                     sep = ","),
                                               sep = ":")
                                     })
                      geno <- matrix(data = geno,
                                     nrow = dim_geno[1],
                                     ncol = dim_geno[2])
                  }
                  geno <- rbind(chr, pos * bp2cm, geno)
                  rownames(geno) <- c("Chr", "Pos", id)
                  write.table(x = geno,file = out_fn, quote = TRUE,
                              row.names = TRUE, col.names = FALSE, sep = ",")
              }

              return(out_fn)
          })

###############################################################################
## Functions to handle the GbsrScheme object in the GbsrGenotypeData object.

## Initialize the GbsrScheme object.
#' @rdname initScheme
#' @importFrom methods slot<-
setMethod("initScheme",
          "GbsrGenotypeData",
          function(object, mating){
              parents <- getParents(object = object, verbose = FALSE)

              if(is.null(parents)){
                  message("Set no parents in the given data.")
                  message("estGeno() will work in the parentless mode unless you specify parents by setParents().")
                  message("See the help of estGeno() for the details of the parentless mode.")
                  p_id <- c(1, 2)
                  mating <- cbind(c(1, 2))

              } else {
                  p_id <- parents$memberID
              }
              scheme <- initScheme(object = slot(object = object, name = "scheme"),
                                   mating = mating,
                                   parents = p_id)

              if(is.null(parents)){
                  slot(object = scheme, name = "parents") <- c(NA_real_, NA_real_)
              }
              slot(object = object, name = "scheme") <- scheme
              return(object)
          })

## Add information to the GbsrScheme object.
#' @rdname addScheme
#' @importFrom methods slot<-
setMethod("addScheme",
          "GbsrGenotypeData",
          function(object, crosstype, mating){
              if(missing(mating)){
                  mating <- cbind(c(NA, NA))
              }
              scheme <- addScheme(object = slot(object = object, name = "scheme"),
                                  crosstype = crosstype,
                                  mating = mating)
              slot(object = object, name = "scheme") <- scheme
              return(object)
          })

## Add information to the GbsrScheme object.
#' @rdname makeScheme
#' @importFrom methods slot<-
setMethod("makeScheme",
          "GbsrGenotypeData",
          function(object, generation, crosstype){
              parents <- getParents(object = object)
              if(is.null(parents)){
                  stop("Parents have not been set.\n Run setParents().")
              }
              if(!length(unique(parents$memberID)) %in% 2^(1:10)){
                  stop("The number of parents should be the Nth power of 2.")
              }
              mating <- matrix(data = seq_along(unique(parents$memberID)),
                               nrow = 2,
                               byrow = TRUE)
              object <- initScheme(object = object, mating = mating)
              if(ncol(mating) > 1){
                  while(TRUE){
                      next_id <- max(mating) + 1
                      n_pair <- ncol(mating) / 2
                      mating <- matrix(data = seq(from = next_id,
                                                  length.out = n_pair),
                                       nrow = 2,
                                       byrow = TRUE)
                      object <- addScheme(object = object,
                                          crosstype = "pairing",
                                          mating = mating)

                      if(ncol(mating) == 1){
                          break
                      }
                  }
              }
              if(generation > 1){
                  for(i in seq_len(generation - 1)){
                      object <- addScheme(object = object,
                                          crosstype = crosstype)
                  }
              }
              return(object)
          })

#'
#' @rdname assignScheme
#' @importFrom methods slot<-
setMethod("assignScheme",
          "GbsrGenotypeData",
          function(object, id) {
              if(nsam(object) != length(id)){
                  stop("`length(id)` should match with the number of valid ",
                       "samples obtained by nsam().",
                       call. = FALSE)
              }
              scheme <- assignScheme(object = slot(object = object,
                                                   name = "scheme"),
                                     id = id)
              slot(object = object, name = "scheme") <- scheme
              sample <- slot(object = object, name = "sample")
              sample$pedigree <- NA
              sample$pedigree[validSam(object = object)] <- id
              slot(object = object, name = "sample") <- sample
              return(object)
          })

## Show the information stored in the GbsrScheme object.
#' @rdname showScheme
setMethod("showScheme",
          "GbsrGenotypeData",
          function(object){
              parents <- getParents(object = object, verbose = FALSE)
              if(is.null(parents)){
                  parents_name <- c(NA, NA)
              } else {
                  parents_name <- parents$sampleID
              }
              sample <- slot(object = object, name = "sample")
              pedigree <- cbind(getSamID(object = object),
                                sample$pedigree[validSam(object = object)])
              showScheme(object = slot(object = object, name = "scheme"),
                         parents_name = parents_name,
                         pedigree = pedigree)
          })
