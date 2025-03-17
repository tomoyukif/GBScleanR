#' @importFrom SeqArray seqOptimize
#' @importFrom gdsfmt put.attr.gdsn index.gdsn
#' @rdname estGeno
setMethod("estGeno",
          "GbsrGenotypeData",
          function(object,
                   node,
                   recomb_rate,
                   error_rate,
                   call_threshold,
                   het_parent,
                   optim,
                   iter,
                   n_threads,
                   dummy_reads) {

              # Check the validity of the registered scheme
              object <- .checkScheme(object = object)

              # Check if the parents have been specified
              parents <- slot(object = slot(object = object, name = "scheme"),
                              name = "parents")
              n_parents <- length(unique(parents))
              parentless <- all(is.na(parents))

              # Set the number of threads
              .setThreads(n_threads = n_threads)

              message("Start cleaning...")

              # Validate parents information
              if(parentless){
                  n_parents <- 2

              } else {
                  parents <- getParents(object = object)
                  n_parents <- length(unique(parents$memberID))
              }
              n_samples <- length(unique(getReplicates(object = object)))
              n_alleles <- 2
              n_ploidy <- attributes(slot(object = object, name = "sample"))$ploidy

              # initialize the output nodes
              .initGDS(object = object,
                       het_parent = het_parent,
                       n_parents = n_parents)

              # Generate genotype pattern list
              message("Preparing genotype and haplotype pattern table...")
              pat <- .makePattern(n_parents = n_parents, n_ploidy = n_ploidy,
                                  n_alleles = n_alleles, n_samples = n_samples,
                                  het_parent = het_parent,
                                  scheme = slot(object = object, name = "scheme"))

              .showPattern(pat = pat)

              # Get chromosome information to loop over chromosomes
              chr <- getChromosome(object = object, valid = FALSE)
              chr_levels <- unique(chr)
              no_eds <- FALSE
              for(chr_i in chr_levels) {
                  if(nmar(object = object, valid = TRUE, chr = chr_i) < 2){
                      message("Skip the genotype estimation for ", chr_i,
                              " that has less than two valid markers.")
                      no_valid_marker <- TRUE

                  } else {
                      no_valid_marker <- FALSE

                      message("\nNow cleaning chr ", chr_i, "...")
                      clean_out <- .cleanEachChr(object = object,
                                                 chr_i = chr_i,
                                                 node = node,
                                                 error_rate = error_rate,
                                                 recomb_rate = recomb_rate,
                                                 call_threshold = call_threshold,
                                                 het_parent = het_parent,
                                                 optim = optim,
                                                 iter = iter,
                                                 parentless = parentless,
                                                 dummy_reads = dummy_reads,
                                                 pat = pat)
                      # ### Debug
                      # return(clean_out)
                      # ###
                  }


                  # Make filters to store the output into the nodes
                  if(parentless){
                      sel <- list(mar = validMar(object = object, chr = chr_i),
                                  sam = validSam(object = object))
                      clean_out$best_hap <- clean_out$best_hap[, -seq_len(2), ]
                      clean_out$best_geno <- clean_out$best_geno[, -seq_len(2), ]

                  } else {
                      sel <- list(mar = validMar(object = object, chr = chr_i),
                                  sam = validSam(object = object, parents = TRUE))
                  }

                  # Output values to the nodes
                  .saveOutput(object = object,
                              clean_out = clean_out,
                              sel = sel,
                              no_valid_marker = no_valid_marker,
                              n_parents = n_parents)

                  # Output the dosage data if available
                  if(n_parents == 2 & !het_parent){
                      .saveEDS(object = object,
                               clean_out = clean_out,
                               sel = sel,
                               no_valid_marker = no_valid_marker)
                  }
              }

              # Optimize the nodes
              if(n_parents == 2 & !het_parent){
                  opt_node <- c("HAP", "CGT", "EDS")
                  comp_node <- c("format/HAP/data", "format/CGT/data",
                                 "format/EDS/data",
                                 "info/PGT", "info/ADB", "info/MR")

              } else {
                  opt_node <- c("HAP", "CGT")
                  comp_node <- c("format/HAP/data", "format/CGT/data",
                                 "info/PGT", "info/ADB", "info/MR")
              }

              closeGDS(object = object, verbose = FALSE)
              seqOptimize(gdsfn = object$filename, target = "by.sample",
                          format.var = opt_node, verbose = FALSE)
              object <- reopenGDS(object = object)
              .compressNodes(object = object,
                             node = paste0("annotation/", comp_node))

              # Clean up RAM
              gc()
              return(object)
          })

################################################################################
# Validate the scheme object
.checkScheme <- function(object){
    scheme <- slot(object = object, name = "scheme")
    if(length(slot(object = scheme, name = "progenies")) == 0){
        stop("No scheme information.",
             "\n",
             "Build scheme information with ",
             "initScheme() and addScheme().",
             call. = FALSE)
    }

    if(length(slot(object = scheme, name = "samples")) == 0){
        id <- unlist(tail(slot(object = scheme, name = "progenies"), 1))
        if(length(id) == 1){
            message("Member IDs were not assigned to samples.",
                    "\n",
                    "Assign ", id, " to all samples as member ID.")
            object <- assignScheme(object = object,
                                   id = rep(x = id, times = nsam(object)))

        } else {
            stop("Member IDs were not assigned to samples.",
                 "\n",
                 "Assign member IDs to samples using assignScheme().",
                 call. = FALSE)
        }
    }
    return(object)
}

################################################################################
################################################################################
# Set the number of threads
.setThreads <- function(n_threads){
    max_threads <- defaultNumThreads()
    if (is.null(n_threads)) { n_threads <- max_threads / 2 }

    if (max_threads <= n_threads) {
        warning("You are going to use all threads",
                "of your computer for the calculation.")
        answer <- ""
        while (answer != "y") {
            answer <- readline("Are you sure to use all threads?(y/n)")
            if (answer == "n") { stop("Stopped by user.") }
        }
        n_threads <- max_threads
    }

    message("Set the number of threads: ", n_threads)
    setThreadOptions(numThreads = n_threads)
}

################################################################################
################################################################################
# Initialize output nodes in the GDS file
#' @importFrom gdsfmt addfolder.gdsn
.initGDS <- function(object, het_parent, n_parents) {
    hap <- addfolder.gdsn(node = index.gdsn(node = object,
                                            path = "annotation/format"),
                          name = "HAP",
                          replace = TRUE)
    add.gdsn(node = hap, name = "data",
             storage = "bit6", compress = "", replace = TRUE)

    if(n_parents == 2 & !het_parent){
        eds <- addfolder.gdsn(node = index.gdsn(node = object,
                                                path = "annotation/format"),
                              name = "EDS",
                              replace = TRUE)
        add.gdsn(node = eds, name = "data",
                 storage = "bit6", compress = "", replace = TRUE)
    }

    cgt <- addfolder.gdsn(node = index.gdsn(node = object,
                                            path = "annotation/format"),
                          name = "CGT",
                          replace = TRUE)
    add.gdsn(node = cgt, name = "data",
             storage = "bit2", compress = "", replace = TRUE)

    pgt <- add.gdsn(node = index.gdsn(node = object,
                                      path = "annotation/info"), name = "PGT",
                    storage = "bit2", compress = "", replace = TRUE)

    adb <- add.gdsn(node = index.gdsn(node = object,
                                      path = "annotation/info"), name = "ADB",
                    storage = "single", compress = "", replace = TRUE)

    mr <- add.gdsn(node = index.gdsn(node = object,
                                     path = "annotation/info"), name = "MR",
                   storage = "single", compress = "", replace = TRUE)
}

################################################################################
################################################################################
# Save the ouput to the GDS file
.saveOutput <- function(object, clean_out, sel, no_valid_marker, n_parents){
    .saveHap(object = object,
             clean_out = clean_out,
             sel = sel,
             no_valid_marker = no_valid_marker)
    .saveGeno(object = object,
              clean_out = clean_out,
              sel = sel,
              no_valid_marker = no_valid_marker)
    .savePGeno(object = object,
               clean_out = clean_out,
               sel = sel,
               no_valid_marker = no_valid_marker,
               n_parents = n_parents)
    .saveADB(object = object,
             clean_out = clean_out,
             sel = sel,
             no_valid_marker = no_valid_marker)
    .saveMR(object = object,
            clean_out = clean_out,
            sel = sel,
            no_valid_marker = no_valid_marker)
}

#' @importFrom gdsfmt append.gdsn
.saveHap <- function(object, clean_out, sel, no_valid_marker) {
    n_ploidy <- dim(clean_out$best_hap)[1]
    output <- array(data = 0,
                    dim = c(n_ploidy, length(sel$sam), length(sel$mar)))

    if(!no_valid_marker){
        rep_id <- getReplicates(object = object, parents = TRUE)
        id_hit <- match(x = rep_id, table = clean_out$mapping_id)
        output[, sel$sam, sel$mar] <- clean_out$best_hap[, id_hit,]
    }

    hap_gdsn <- index.gdsn(node = object, path = "annotation/format/HAP/data")
    gdsn_dim <- objdesp.gdsn(node = hap_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(node = index.gdsn(node = object, path = "annotation/format/HAP"),
                 name = "data", val = output,
                 storage = "bit6", compress = "", replace = TRUE)
    } else {
        append.gdsn(node = hap_gdsn, val = output)
    }
}

.saveEDS <- function(object, clean_out, sel, no_valid_marker) {
    n_ploidy <- dim(clean_out$best_hap)[1]
    output <- array(data = NA, dim = c(n_ploidy, length(sel$sam), length(sel$mar)))

    if(!no_valid_marker){
        rep_id <- getReplicates(object = object, parents = TRUE)
        id_hit <- match(x = rep_id, table = clean_out$mapping_id)
        output[, sel$sam, sel$mar] <- clean_out$best_hap[, id_hit,]
    }

    output[output == 0] <- NA
    output <- apply(X = output - 1, MARGIN = c(2, 3), FUN = sum)
    output[is.na(output)] <- 63

    eds_gdsn <- index.gdsn(node = object, path = "annotation/format/EDS/data")
    gdsn_dim <- objdesp.gdsn(node = eds_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(node = index.gdsn(node = object, path = "annotation/format/EDS"),
                 name = "data", val = output,
                 storage = "bit6", compress = "", replace = TRUE)
    } else {
        append.gdsn(node = eds_gdsn, val = output)
    }
}

#' @importFrom gdsfmt append.gdsn
.saveGeno <- function(object, clean_out, sel, no_valid_marker) {
    n_ploidy <- dim(clean_out$best_geno)[1]
    output <- array(data = 3, dim = c(n_ploidy, length(sel$sam), length(sel$mar)))

    if(!no_valid_marker){
        rep_id <- getReplicates(object = object, parents = TRUE)
        id_hit <- match(x = rep_id, table = clean_out$mapping_id)
        output[, sel$sam, sel$mar] <- clean_out$best_geno[, id_hit,]
    }

    out_gdsn <-index.gdsn(node = object, path = "annotation/format/CGT/data")
    gdsn_dim <- objdesp.gdsn(node = out_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(node = index.gdsn(node = object, path = "annotation/format/CGT"),
                 name = "data", val = output,
                 storage = "bit2", compress = "", replace = TRUE)
    } else {
        append.gdsn(node = out_gdsn, val = output)
    }
}

#' @importFrom gdsfmt append.gdsn
.savePGeno <- function(object, clean_out, sel, no_valid_marker, n_parents) {
    n_row <- dim(clean_out$p_geno)[1]
    output <- matrix(data = 3, nrow = n_row, length(sel$mar))

    if(!no_valid_marker){
        output[, sel$mar] <- clean_out$p_geno
    }

    out_gdsn <- index.gdsn(node = object, path = "annotation/info/PGT")
    gdsn_dim <- objdesp.gdsn(node = out_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(node = index.gdsn(node = object, path = "annotation/info"),
                 name = "PGT", val = output,
                 storage = "bit2", compress = "", replace = TRUE)

    } else {
        append.gdsn(node = out_gdsn, val = output)
    }
}

#' @importFrom gdsfmt append.gdsn
.saveADB <- function(object, clean_out, sel, no_valid_marker) {
    output <- rep(-1, length(sel$mar))

    if(!no_valid_marker){
        output[sel$mar] <- t(clean_out$bias)
    }

    out_gdsn <- index.gdsn(node = object, path = "annotation/info/ADB")
    gdsn_dim <- objdesp.gdsn(node = out_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(node = index.gdsn(node = object, path = "annotation/info"),
                 name = "ADB", val = output,
                 storage = "single", compress = "", replace = TRUE)

    } else {
        append.gdsn(node = out_gdsn, val = output)
    }
}

#' @importFrom gdsfmt append.gdsn
.saveMR <- function(object, clean_out, sel, no_valid_marker) {
    output <- matrix(-1, nrow = 2, ncol = length(sel$mar))

    if(!no_valid_marker){
        output[, sel$mar] <- t(clean_out$mismap)
    }

    out_gdsn <- index.gdsn(object, "annotation/info/MR")
    gdsn_dim <- objdesp.gdsn(out_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(index.gdsn(object, "annotation/info"), "MR", output,
                 "single", compress = "", replace = TRUE)
    } else {
        append.gdsn(out_gdsn, output)
    }
}

################################################################################
################################################################################
# Prepare read count data
.loadReadCounts <- function(object, chr_i, node, parentless, dummy_reads) {
    ad_node <- "raw"
    if (exist.gdsn(node = object, path = "annotation/format/FAD")) {
        if(node == "filt"){
            ad_node <- "filt"
        }
    }

    reads <- getRead(object = object, node = ad_node,
                     parents = FALSE, valid = TRUE, chr = chr_i)
    rep_id <- getReplicates(object = object, parents = FALSE)
    reads <- .pileupAD(ad = reads, rep_id = rep_id)

    if(parentless){
        ref_read <- matrix(data = c(dummy_reads, 0),
                           nrow = 2,
                           ncol = ncol(reads$ref))
        alt_read <- matrix(data = c(0, dummy_reads),
                           nrow = 2,
                           ncol = ncol(reads$ref))
        p_reads <- list(ref = ref_read, alt = alt_read)

    } else {
        p_reads <- getRead(object = object, node = ad_node,
                           parents = "only", valid = TRUE, chr = chr_i)
        rep_id <- getReplicates(object = object, parents = "only")
        p_reads <- .pileupAD(ad = p_reads, rep_id = rep_id)
        p_reads <- .orderParents(object = object, p_reads = p_reads)
    }
    # p_reads <- .bumpOverRepReads(reads = p_reads)
    # reads <- .bumpOverRepReads(reads = reads)
    return(list(p_ref = p_reads$ref,
                p_alt = p_reads$alt,
                ref = reads$ref,
                alt = reads$alt))
}
#
# .bumpOverRepReads <- function(reads){
#     r_homo <- reads$ref > 0 & reads$alt == 0
#     a_homo <- reads$ref == 0 & reads$alt > 0
#     het <- reads$ref > 0 & reads$alt > 0
#     r_homo_over <- reads$ref[r_homo] > 50
#     reads$ref[r_homo][r_homo_over] <- 50
#     a_homo_over <- reads$alt[a_homo] > 50
#     reads$alt[a_homo][a_homo_over] <- 50
#     het_dp <- reads$ref[het] + reads$alt[het]
#     het_dp_over <- het_dp > 50
#     reads$ref[het][het_dp_over] <- round(reads$ref[het][het_dp_over] / het_dp[het_dp_over] * 50)
#     reads$alt[het][het_dp_over] <- round(reads$alt[het][het_dp_over] / het_dp[het_dp_over] * 50)
#     return(reads)
# }

#' @importFrom stats na.omit
.orderParents <- function(object, p_reads){
    p_id <- getParents(object = object)
    rep_id <- getReplicates(object = object, parents = "only")
    p_indices <- validSam(object = object, parents = "only")
    sam_id <- getSamID(object = object, valid = FALSE)[p_indices]
    id_hit <- match(x = sam_id, table = p_id$sampleID)
    p_id$rep_id[na.omit(id_hit)] <- rep_id[!is.na(id_hit)]
    member_id <- p_id$memberID[match(x = rownames(p_reads$ref),
                                     table = p_id$rep_id)]
    p_reads$ref <- p_reads$ref[order(member_id), ]
    p_reads$alt <- p_reads$alt[order(member_id), ]
    return(p_reads)
}

################################################################################
################################################################################
# Make genotype and haplotype pattern list
.makeGenoPat <- function(n_ploidy, alleles){
    out <- NULL
    for (i in seq_len(n_ploidy)) {
        out <- c(out, list(alleles))
    }
    out <- expand.grid(out)
    out <- t(apply(X = out, MARGIN = 1, FUN = sort))
    out <- rowSums(out[!duplicated(out),])
    return(out)
}

.makeGenoParents <- function(n_parents, n_ploidy, alleles, het_parent){
    out <- NULL
    for (i in seq_len(n_parents * n_ploidy)) {
        out <- c(out, list(alleles))
    }
    out <- expand.grid(out, KEEP.OUT.ATTRS = FALSE)
    valid_pat <- apply(X = out, MARGIN = 1,
                       FUN = function(x) length(unique(x)) > 1)
    out <- as.matrix(out[valid_pat,])
    if (!het_parent) {
        valid <- tapply(X = seq_len(ncol(out)),
                        INDEX = rep(seq_len(n_parents), each = n_ploidy),
                        FUN = function(i){
                            mono <- apply(out[, i], 1, function(x){
                                return(length(unique(x)) == 1)
                            })
                            return(mono)
                        })
        valid <- do.call("rbind", valid)
        valid <- apply(valid, 2, all)
        out <- out[valid, ]
    }
    attributes(out) <- list(dim = dim(out))
    return(out)
}

.progenyPattern <- function(zygotes, index, mt, n_ploidy){
    gamete_ploidy <- n_ploidy / 2
    rev_index <- rev(seq_len(gamete_ploidy))

    out <- vapply(X = seq_len(ncol(mt[[index]])),
                  FUN = .mateGametes,
                  FUN.VALUE = list(0),
                  zygotes = zygotes,
                  mate = mt[[index]],
                  n_ploidy = n_ploidy,
                  gamete_ploidy = gamete_ploidy,
                  rev_index = rev_index)
    out <- c(zygotes, out)

    return(out)
}

#' @importFrom utils combn
.mateGametes <- function(index, zygotes, mate, n_ploidy, gamete_ploidy, rev_index){
    i_mate <- mate[, index]
    zygote1 <- zygotes[[i_mate[1]]]
    zygote2 <- zygotes[[i_mate[2]]]
    comb1 <- combn(seq_len(n_ploidy), gamete_ploidy)
    comb2 <- combn(seq_len(n_ploidy), gamete_ploidy)
    if(length(rev_index) > 1){
        comb1 <- cbind(comb1, comb1[rev_index, ])
        comb2 <- cbind(comb2, comb2[rev_index, ])
    }
    gamete1 <- matrix(apply(zygote1, 1, "[", comb1), nrow = gamete_ploidy)
    gamete1 <- unique(t(gamete1))
    gamete2 <- matrix(apply(zygote2, 1, "[", comb2), nrow = gamete_ploidy)
    gamete2 <- unique(t(gamete2))
    comb <- expand.grid(seq_len(nrow(gamete1)), seq_len(nrow(gamete2)))
    if(nrow(comb) == 1){
        comb <- rbind(comb, comb)
    }
    zygotes <- unique(rbind(cbind(gamete1[comb[, 1], ],
                                  gamete2[comb[, 2], ]),
                            cbind(gamete2[comb[, 1], ],
                                  gamete1[comb[, 2], ])))
    return(list(zygotes))
}

.getValidPat <- function(scheme, het_parent, n_parents, n_ploidy) {
    target_pedigree <- sort(unique(slot(object = scheme, name = "samples")))
    mt <- slot(object = scheme, name = "mating")
    xtype <- slot(object = scheme, name = "crosstype")
    pg <- slot(object = scheme, name = "progenies")
    even_ploidy <- n_ploidy %% 2 != 0
    if(even_ploidy){
        n_ploidy_minus <- n_ploidy - 1
        possible_zygotes_ids <- lapply(seq_len(n_parents), function(i){
            n <- seq(n_ploidy_minus * (i - 1) + 1, length.out = n_ploidy_minus)
            if(!het_parent){
                n <- rep(n[1], n_ploidy_minus)
            }
            return(t(n))
        })
        n_ploidy <- n_ploidy + 1
    }

    # Simulate possible zygotes and gametes
    # Founder zygotes
    zygotes <- lapply(seq_len(n_parents), function(i){
        n <- seq(n_ploidy * (i - 1) + 1, length.out = n_ploidy)
        if(!het_parent){
            n <- rep(n[1], n_ploidy)
        }
        return(t(n))
    })

    # Progeny zygotes
    for(i in seq_along(mt)) {
        zygotes <- .progenyPattern(zygotes = zygotes,
                                   index = i,
                                   mt = mt,
                                   n_ploidy = n_ploidy)
    }

    out <- zygotes[as.numeric(target_pedigree)]
    if(even_ploidy){
        id_table <- lapply(seq_along(possible_zygotes_ids), function(i){
            return(cbind(as.vector(zygotes[[i]]), as.vector(possible_zygotes_ids[[i]])))
        })
        id_table <- do.call("rbind", id_table)
        out <- lapply(out, function(x){
            hit <- match(x, id_table[, 1])
            x <- matrix(id_table[hit, 2], nrow = nrow(x))
            x <- unique(x[, -1])
        })
    }
    names(out) <- target_pedigree

    return(out)
}

.getPossibleHap <- function(hap_progeny, geno_parents, geno_pat, n_ploidy){
    hap_vec <- as.vector(t(hap_progeny))

    derived_geno <- apply(X = geno_parents, MARGIN = 1,
                          FUN = function(x){
                              x <- x[hap_vec]
                              x <- matrix(x, nrow = n_ploidy)
                              x <- colSums(x)
                              return(x)
                          })
    derived_geno <- as.vector(derived_geno)
    out <- rep(derived_geno, n_ploidy + 1)
    out <- matrix(out, nrow = n_ploidy + 1, byrow = TRUE)
    out <- out == geno_pat
    out <- apply(X = out, MARGIN = 2, FUN = which)
    return(out)
}

.getPossibleGeno <- function(geno_parents, geno_pat, n_ploidy){
    derived_geno <- apply(X = geno_parents, MARGIN = 1, FUN = function(x) {
        x <- matrix(x, nrow = n_ploidy)
        x <- colSums(x)
        return(x)
    })
    derived_geno <- as.vector(derived_geno)
    out <- rep(derived_geno, n_ploidy + 1)
    out <- matrix(out, nrow = n_ploidy + 1, byrow = TRUE)
    out <- out == geno_pat[seq_len(n_ploidy + 1)]
    out <- apply(X = out, MARGIN = 2, FUN = which)
    return(out)
}

.makePattern <- function(n_parents,
                         n_ploidy,
                         n_alleles,
                         n_samples,
                         het_parent,
                         scheme) {
    alleles <- seq(0, length.out = n_alleles)

    geno_pat <- .makeGenoPat(n_ploidy = n_ploidy, alleles = alleles)

    even_ploidy <- n_ploidy %% 2 != 0
    if(even_ploidy){
        n_ploidy <- n_ploidy - 1
    }
    geno_parents <- .makeGenoParents(n_parents = n_parents,
                                     n_ploidy = n_ploidy,
                                     alleles = alleles,
                                     het_parent = het_parent)

    possiblegeno <- .getPossibleGeno(geno_parents = geno_parents,
                                     geno_pat = geno_pat,
                                     n_ploidy = n_ploidy)

    if(even_ploidy){
        n_ploidy <- n_ploidy + 1
    }
    hap_progeny <- .getValidPat(scheme = scheme,
                                het_parent = het_parent,
                                n_parents = n_parents,
                                n_ploidy = n_ploidy)

    possiblehap <- lapply(X = hap_progeny,
                          FUN = .getPossibleHap,
                          geno_parents = geno_parents,
                          geno_pat = geno_pat,
                          n_ploidy = n_ploidy)

    n_p_pat <- nrow(geno_parents)
    n_hap_pat <- vapply(X = seq_along(hap_progeny),
                        FUN.VALUE = numeric(1),
                        FUN = function(i){
                            return(nrow(hap_progeny[[i]]))
                        })
    return(list(alleles = alleles,
                geno_pat = geno_pat,
                geno_parents = geno_parents,
                hap_progeny = hap_progeny,
                possiblehap = possiblehap,
                possiblegeno = possiblegeno,
                n_p_pat = n_p_pat,
                n_hap_pat = n_hap_pat,
                n_samples = n_samples))
}

.showPattern <- function(pat){
    message("Possible allele dosages: ", paste(pat$geno_pat, collapse = " "))
    message("Number of possible founder genotypes: ", pat$n_p_pat)
    message("Member ID(s) to be processed: ", paste(names(pat$hap_progeny), collapse = " "))
    message("Number of offspring haplotypes: ", paste(pat$n_hap_pat, collapse = " "))
}

################################################################################
################################################################################
# Prepare the transition provability matrix
.getJnum <- function(scheme, n_origin, het_parent, i) {
    xtype <- vapply(X = slot(object = scheme, name = "crosstype"),
                    FUN.VALUE = character(length = 1),
                    FUN = function(x){
                        return(x[1])
                    })
    n_x <- length(xtype)

    if(n_x > 1){
        for(j in seq(n_x, 2)){
            if(xtype[j] == "pairing"){
                j_mating <- slot(object = scheme, name = "mating")[[j]]
                i_mating <- j_mating[, i]
                prev_mating <- slot(object = scheme, name = "mating")[[j - 1]]
                prev_max <- max(prev_mating)
                j_matings <- prev_mating[, i_mating - prev_max]

                if(any(j_matings[, 1] %in% j_matings[, 2])){
                    xtype[j] <- "sibling"
                }
            }
        }
    }

    i_pairing <- grep(pattern = "pairing", x = xtype)
    if (length(i_pairing) == 0) {
        n_pairing  <- 0

    } else {
        n_pairing <- max(i_pairing)
    }

    if (n_pairing == 0) {
        jnum <- .initJnum(n_origin = n_origin)

    } else {
        next_crosstype <- xtype[n_pairing + 1]
        if(is.na(next_crosstype)){
            next_crosstype <- "selfing"
        }
        jnum <- .initJnum(n_origin = n_origin,
                          n_pairing = n_pairing + het_parent,
                          next_crosstype = next_crosstype)
    }

    s_gen <- n_pairing + 1
    n_gen <- length(xtype)
    if(s_gen <= n_gen){
        for (i in s_gen:n_gen) {
            jnum <- .calcNextJnum(crosstype = xtype[i], prob_df = jnum)
        }
    }
    return(jnum)
}

.initJnum <- function(n_origin, n_pairing, next_crosstype) {
    jnum <- data.frame(a12 = 0,
                       b12 = 1,
                       a123 = 0,
                       b123 = 0,
                       r = 0,
                       j1232 = 0,
                       k1232 = 0,
                       j1122 = 0,
                       k1122 = 0,
                       j1222 = 0)
    if (n_origin > 2) {
        jnum$b123 <- 1
    }
    if (!missing(n_pairing)) {
        jnum$a12 <- 1
        if (next_crosstype == "selfing") {
            jnum$j1232 <- jnum$r <- n_pairing - 1

        } else {
            r <- ifelse(test = n_pairing == 1, yes = 0, no = n_pairing - 2)
            jnum$j1232 <- jnum$k1232 <- jnum$r <- r

            if (n_origin > 2) {
                jnum$a123 <- 1
            } else {
                jnum$b12 <- 0.5
            }
            if (next_crosstype == "sibling") {
                while (n_pairing > 1) {
                    jnum <- .calcNextJnum(crosstype = "sibling", prob_df = jnum)
                    n_pairing <- n_pairing - 1
                }
            }
        }
    }
    return(jnum)
}

.calcNextJnum <- function(crosstype, prob_df) {
    coal <- switch(crosstype,
                   "pairing" = c(0, 0),
                   "selfing" = c(1, 0),
                   "sibling" = c(0.5, 0))
    s <- coal[1]
    q <- coal[2]
    if (crosstype == "selfing") {
        next_a12 <- 0.5 * prob_df$a12
        next_b12 <- 0
        next_a123 <- 0
        next_b123 <- 0
        next_r <- prob_df$r + prob_df$a12
        next_j1232 <- 0.5 * prob_df$j1232
        next_k1232 <- 0
        next_j1122 <- 0.5 * prob_df$j1122 + 0.5 * prob_df$r
        next_k1122 <- 0
        next_j1222 <- 0.5 * (next_r - next_j1122 - next_j1232)

    } else {
        next_a12 <- prob_df$b12
        next_b12 <- 0.5 * s * prob_df$a12 + (1 - s) * prob_df$b12
        next_a123 <- s * prob_df$a123 + (1 - 2 * s) * prob_df$b123
        next_b123 <-
            1.5 * q * prob_df$a123 + (1 - s - 2 * q) * prob_df$b123
        next_r <- prob_df$r + prob_df$a12
        next_j1232 <- prob_df$k1232 + prob_df$a123
        next_k1232 <-
            0.5 * s * prob_df$j1232 + (1 - s) * (prob_df$k1232 + prob_df$a123)
        next_j1122 <- prob_df$k1122
        next_k1122 <-
            0.5 * s * (prob_df$r + prob_df$j1122) + (1 - s) * prob_df$k1122
        next_j1222 <- 0.5 * (next_r - next_j1122 - next_j1232)
    }
    return(data.frame(a12 = next_a12,
                      b12 = next_b12,
                      a123 = next_a123,
                      b123 = next_b123,
                      r = next_r,
                      j1232 = next_j1232,
                      k1232 = next_k1232,
                      j1122 = next_j1122,
                      k1122 = next_k1122,
                      j1222 = next_j1222))
}

.getXoFreq <- function(jnum, n_origin) {
    r00 <- jnum$j1232 / jnum$a12 / (n_origin - 2)
    r01 <- jnum$j1222 / jnum$a12
    r10 <- jnum$j1222 / (1 - jnum$a12) / (n_origin - 1)
    r11 <- jnum$j1122 / (1 - jnum$a12) / (n_origin - 1)
    return(data.frame(r00, r01, r10, r11))
}

.transitionProb <- function(pat, pos, recomb_rate,
                            scheme, het_parent, n_ploidy){
    if(n_ploidy %% 2 != 0){
        even_ploidy <- TRUE

    } else {
        even_ploidy <- FALSE
    }

    out <- lapply(X = seq_along(pat$hap_progeny),
                  FUN = function(i){
                      hap_progeny_i <- pat$hap_progeny[[i]]
                      n_origin <- length(unique(as.vector(hap_progeny_i)))
                      jnum <- .getJnum(scheme = scheme,
                                       n_origin = n_origin,
                                       het_parent = het_parent,
                                       i = i)
                      jrate <- .getXoFreq(jnum = jnum, n_origin = n_origin)
                      jrate[is.na(jrate)] <- 0
                      q_mat <- apply(X = hap_progeny_i,
                                     MARGIN = 1,
                                     FUN = function(x) {
                                         apply(X = hap_progeny_i,
                                               MARGIN = 1,
                                               FUN = function(y) {

                                                   # If no change in adjacent genotypes, set NA.
                                                   # These NAs will be replaced with negated rowSums() values of each row of q_mat.
                                                   no_change <- x == y
                                                   if(sum(no_change) == n_ploidy){
                                                       out <- NA

                                                   } else {
                                                       # If more than two patterns of crossing overs were observed, set zero.
                                                       # We do not assume multiple crossing overs at one joint.
                                                       joint_pat <- paste(x, y, sep = "_")
                                                       joint_pat <- joint_pat[!no_change]
                                                       if(length(unique(joint_pat)) > 1){
                                                           out <- 0

                                                       } else {
                                                           # If only one pattern of crossing over in multiple chromosomes were found,
                                                           # set the r11 value to the power of the inverse of the number of crossing over chromosomes.
                                                           if(length(joint_pat) > 1){
                                                               out <- jrate$r11^(1/length(joint_pat))

                                                           } else {

                                                               # Check the number of the major haplotype before and after the joint
                                                               check1 <- max(table(x))
                                                               check2 <- max(table(y))
                                                               if(even_ploidy){
                                                                   if(check1 == check2){
                                                                       if(check1 != n_ploidy){
                                                                           check2 <- check2 - 1
                                                                       }
                                                                   }
                                                               }

                                                               # If the major joint decreased, set r01 that is the rate to get an additional non-ibd state.
                                                               if(check1 > check2){
                                                                   out <- jrate$r01

                                                               # If the major joint increased, set r10 that is the rate to get an additional non-ibd state.
                                                               } else if(check1 < check2){
                                                                   out <- jrate$r10

                                                               } else {
                                                                   out <- jrate$r00
                                                               }
                                                           }
                                                       }
                                                   }
                                                   return(out)
                                               })
                                     })
                      diag(q_mat) <- -rowSums(q_mat, na.rm = TRUE)
                      mar_dist <- diff(pos)
                      if(any(mar_dist < 0)){
                          stop("Negative value(s) was obtained in the calculation of marker distances \n",
                               "based on physical positions of the given markers.",
                               "\nThe order of the markers might not be valid.",
                               call. = FALSE)
                      }
                      rf <- mar_dist * 1e-6 * recomb_rate
                      prob <- vapply(X = rf,
                                     FUN = function(x) {
                                         out <- expm(q_mat * x, "Higham08")
                                         out <- out[q_mat != 0]
                                         return(out)
                                     },
                                     FUN.VALUE = numeric(length(q_mat[q_mat != 0])))
                      non_zero <- q_mat != 0
                      return(list(prob = prob,
                                  non_zero = non_zero,
                                  n_hap_pat = pat$n_hap_pat))
                  })
    return(out)
}

.getInitProb <- function(prob, n_samples) {
    full_mat <- matrix(data = 0, nrow = prob$n_hap_pat, ncol = prob$n_hap_pat)
    full_mat[prob$non_zero] <- prob$prob[, 1]
    ev1 <- eigen(t(full_mat))$vectors[, 1]
    init <- ev1 / sum(ev1)
    return(init)
}

################################################################################
################################################################################
# Make parameter lists
.getPedigree <- function(object){
    vaild_sam <- validSam(object = object)
    out <- slot(object = object, name = "sample")$pedigree[vaild_sam]
    return(out)
}

.getParams <- function(object, chr_i, node, error_rate, recomb_rate,
                       call_threshold, het_parent,
                       parentless, dummy_reads, pat) {
    reads <- .loadReadCounts(object = object, chr_i = chr_i, node = node,
                             parentless = parentless,
                             dummy_reads = dummy_reads)

    if(parentless){
        n_parents <- 2

    } else {
        parents <- getParents(object = object)
        n_parents <- length(unique(parents$memberID))
    }
    n_samples <- length(unique(getReplicates(object = object)))
    n_alleles <- 2
    n_ploidy <- attributes(slot(object = object, name = "sample"))$ploidy

    pos <- getPosition(object = object, valid = TRUE, chr = chr_i)
    n_mar <- nmar(object = object, valid = TRUE, chr = chr_i)
    trans_prob <- .transitionProb(pat = pat, pos = pos,
                                  recomb_rate = recomb_rate,
                                  scheme = slot(object = object, name = "scheme"),
                                  het_parent = het_parent,
                                  n_ploidy = n_ploidy)
    init_prob <- lapply(X = trans_prob,
                        FUN = .getInitProb,
                        n_samples = n_samples)

    fixed_param <- getFixedParameter(object = object)
    if(!is.null(fixed_param$bias)){
        bias <- fixed_param$bias

    } else {
        bias <- rep(x = 0.5, times = n_mar)
    }

    if(!is.null(fixed_param$mismap_ref)){
        mismap <- cbind(fixed_param$mismap_ref, fixed_param$mismap_alt)

    } else {
        mismap <-  matrix(data = 0.005, nrow = n_mar, ncol = 2)
    }

    if(!is.null(fixed_param$parent1)){
        p_geno_fix <- apply(fixed_param[, grepl("parent[0-9]", names(fixed_param))],
                            1,
                            function(x){
                                return(which(colSums(t(pat$geno_parents) == x) == length(x)))
                            })
    } else {
        p_geno_fix <- -1
    }


    return(list(n_parents = n_parents,
                n_samples = n_samples,
                parents_ids = rownames(reads$p_ref),
                sample_ids = rownames(reads$ref),
                n_alleles = n_alleles,
                n_ploidy = n_ploidy,
                n_mar = n_mar,
                n_pedigree = length(trans_prob),
                pedigree = .getPedigree(object = object),
                het_parent = het_parent,
                error_rate = c(1 - error_rate, error_rate),
                recomb_rate = recomb_rate,
                call_threshold = call_threshold,
                reads = reads,
                pat = pat,
                bias = bias,
                mismap = mismap,
                fixed_param = fixed_param,
                trans_prob = trans_prob,
                init_prob = init_prob,
                count = 0,
                p_geno_fix = p_geno_fix,
                flip = FALSE))
}

################################################################################
################################################################################
# Execute the Viterbi algorithm
.getBestSeq <- function(param_list, outprob) {
    trans_prob <- lapply(param_list$trans_prob, function(x){
        return(as.vector(x$prob))
    })
    trans_prob <- log10(unlist(trans_prob))
    init_prob <- log10(unlist(param_list$init_prob))
    nonzero_prob <- lapply(param_list$trans_prob, function(x){
        return(x$non_zero)
    })
    nonzero_prob <- unlist(nonzero_prob)
    n_nonzero_prob <- sapply(param_list$trans_prob, function(x){
        return(sum(x$non_zero))
    })
    possiblehap <- unlist(param_list$pat$possiblehap)
    pedigree <- match(x = param_list$pedigree,
                      table = names(param_list$pat$hap_progeny))

    out_list <- run_viterbi(p_ref = param_list$reads$p_ref,
                            p_alt = param_list$reads$p_alt,
                            ref = param_list$reads$ref,
                            alt = param_list$reads$alt,
                            eseq_in = param_list$error_rate,
                            bias = param_list$bias,
                            mismap = param_list$mismap,
                            trans_prob = trans_prob,
                            init_prob = init_prob,
                            nonzero_prob = nonzero_prob,
                            n_pgeno = param_list$pat$n_p_pat,
                            n_hap = param_list$pat$n_hap_pat,
                            n_offspring = param_list$n_samples,
                            n_founder = param_list$n_parents,
                            n_marker = param_list$n_mar,
                            n_nonzero_prob = n_nonzero_prob,
                            het = param_list$het_parent,
                            pedigree = pedigree - 1,
                            possiblehap = possiblehap - 1,
                            possiblegeno = param_list$pat$possiblegeno - 1,
                            p_geno_fix = param_list$p_geno_fix - 1,
                            ploidy = param_list$n_ploidy)

    if (outprob) {
        prob <- run_fb(ref = param_list$reads$ref,
                       alt = param_list$reads$alt,
                       eseq_in = param_list$error_rate,
                       bias = param_list$bias,
                       mismap = param_list$mismap,
                       possiblehap = possiblehap - 1,
                       trans_prob = trans_prob,
                       init_prob = init_prob,
                       nonzero_prob = nonzero_prob,
                       n_pgeno = param_list$pat$n_p_pat,
                       n_hap = param_list$pat$n_hap_pat,
                       n_offspring = param_list$n_samples,
                       n_marker = param_list$n_mar,
                       n_nonzero_prob = n_nonzero_prob,
                       pedigree = pedigree - 1,
                       p_geno = out_list$p_geno,
                       ploidy = param_list$n_ploidy)
        prob <- array(data = apply(X = prob,
                                   MARGIN = 3,
                                   FUN = function(x) return(t(x))),
                      dim = c(param_list$n_ploidy + 1,
                              param_list$n_samples,
                              param_list$n_mar))
        out_list$prob <- prob

    } else {
        out_list$prob <- NULL
    }

    out_list$best_seq <- out_list$best_seq + 1
    out_list$p_geno <- out_list$p_geno + 1
    return(out_list)
}

################################################################################
################################################################################
# Convert data format the Viterbi algorithm's output to the GBScleanR style
.pat2hap <- function(best_hap, param_list) {
    n_parents_chr <- param_list$n_parents * param_list$n_ploidy

    if(param_list$het_parent){
        parent_hap <- rep(x = seq_len(n_parents_chr), each = param_list$n_mar)

    } else {
        parent_hap <- rep(x = seq(1, by = param_list$n_ploidy,
                                  length.out = param_list$n_parents),
                          each = param_list$n_ploidy * param_list$n_mar)
    }
    parent_hap <- matrix(data = parent_hap, nrow = n_parents_chr, byrow = TRUE)
    sample_pat <- best_hap$best_seq
    pat_i <- match(x = param_list$pedigree,
                   table = names(param_list$pat$hap_progeny))
    sample_hap <- vapply(X = seq_len(param_list$n_samples),
                         FUN.VALUE = numeric(param_list$n_ploidy * param_list$n_mar),
                         FUN = function(i) {
                             hap_i <- param_list$pat$hap_progeny[[pat_i[i]]]
                             out <- as.vector(hap_i[sample_pat[, i], ])
                             return(out)
                         })
    sample_hap <- matrix(data = sample_hap,
                         nrow = param_list$n_samples * param_list$n_ploidy,
                         ncol = param_list$n_mar,
                         byrow = TRUE)
    out_hap <- rbind(parent_hap, sample_hap)
    out_hap <- array(data = out_hap,
                     dim = c(param_list$n_ploidy,
                             param_list$n_parents + param_list$n_samples,
                             param_list$n_mar))
    return(out_hap)
}

.parentPat2Geno <- function(best_hap, param_list) {
    p_pat <- best_hap$p_geno
    p_geno <- t(param_list$pat$geno_parents[p_pat, ])
    return(p_geno)
}

.hap2geno <- function(hap, p_geno, param_list) {
    n <- param_list$n_ploidy * (param_list$n_samples + param_list$n_parents)
    out_geno <- vapply(X = seq_len(param_list$n_mar),
                       FUN.VALUE = numeric(n),
                       FUN = function(x) {
                           return(p_geno[hap[,, x], x])
                       })
    out_geno <- array(data = out_geno,
                      dim = c(param_list$n_ploidy,
                              param_list$n_samples + param_list$n_parents,
                              param_list$n_mar))

    return(out_geno)
}

.halfJoint <- function(best_pat_f, best_pat_r, param_list) {
    n_mar <- param_list$n_mar
    n_sample <- param_list$n_samples
    half <- round(n_mar / 2)
    first <- (n_mar:1)[seq_len(half)]
    latter <- (half + 1):n_mar

    best_seq <- rbind(best_pat_r$best_seq[first,],
                      best_pat_f$best_seq[latter,])
    p_geno <- c(best_pat_r$p_geno[first], best_pat_f$p_geno[latter])

    prob <- c(best_pat_r$prob[, , first],
              best_pat_f$prob[, , latter])
    prob <- array(data = prob, dim = dim(best_pat_f$prob))
    prob <- apply(X = prob, MARGIN = 2, FUN = function(x) return(x))
    prob <- array(data = prob, dim = c(param_list$n_ploidy + 1, n_mar, n_sample))
    out_list <- list(best_seq = best_seq,
                     p_geno = p_geno,
                     prob = prob)
    return(out_list)
}

.solveBorderConflict <- function(best_hap, param_list){
    half <- round(param_list$n_mar / 2)
    flip_i <- vapply(X = seq_len(dim(best_hap)[2]),
                     FUN = .checkBorder,
                     best_hap = best_hap,
                     param_list = param_list,
                     half = half,
                     FUN.VALUE = logical(length = 1L))
    half_m <- -seq(1, half)
    reorder <- rev(seq_len(param_list$n_ploidy))
    best_hap[, flip_i, half_m] <- best_hap[reorder, flip_i, half_m]
    return(best_hap)
}

.checkBorder <- function(i, best_hap, param_list, half){
    if(i <= param_list$n_parents){
        return(FALSE)
    }
    i_redigree <- param_list$pedigree[i - param_list$n_parents]
    i_pedigree <- names(param_list$pat$hap_progeny) %in% i_redigree
    i_pedigree <- which(i_pedigree)
    i_hap_progeny <- t(param_list$pat$hap_progeny[[i_pedigree]])
    hap1 <- which(colSums(i_hap_progeny == best_hap[, i, half]) == param_list$n_ploidy)
    hap2 <- which(colSums(i_hap_progeny == best_hap[, i, half + 1]) == param_list$n_ploidy)
    re <- rev(seq_along(best_hap[, i, half + 1]))
    hap3 <- which(colSums(i_hap_progeny == best_hap[re, i, half + 1]) == param_list$n_ploidy)
    n_hap_pat <- param_list$trans_prob[[i_pedigree]]$n_hap_pat
    prob <- matrix(data = 0, nrow = n_hap_pat, ncol = n_hap_pat)
    non_zero <- param_list$trans_prob[[i_pedigree]]$non_zero
    prob_half <- param_list$trans_prob[[i_pedigree]]$prob[, half]
    prob[non_zero] <- prob_half
    prob1 <- prob[hap1, hap2]
    prob2 <- prob[hap1, hap3]
    if(prob1 < prob2){
        out <- TRUE

    } else {
        out <- FALSE
    }
    return(out)
}

.minimizeRecombination <- function(best_hap){
    for(i in seq_len(dim(best_hap)[2])){
        recomb <- apply(X = best_hap[, i,], MARGIN = 1, FUN = diff) != 0
        n_recomb <- rowSums(recomb)
        recomb_pos <- which(n_recomb >= 1)
        for(j in seq_along(recomb_pos)){
            if(j == 1){
                prev_xo <- recomb[recomb_pos[j], ]

            } else {
                is_homo_j <- all(diff(best_hap[, i, recomb_pos[j]]) == 0)
                if(is_homo_j){
                    j_xo <- recomb[recomb_pos[j], ]
                    if(all(j_xo == prev_xo)){
                        re_order <- c(seq(2, length(j_xo)), 1)
                        prev_xo <- j_xo[re_order]
                        flip_i <- -seq(1, recomb_pos[j])
                        best_hap[, i, flip_i] <- best_hap[re_order, i, flip_i]
                    }
                }
            }
        }
    }
    return(best_hap)
}

.summarizeEst <- function(best_hap, best_geno, pat_prob, param_list) {
    n_mar <- param_list$n_mar
    n_sample <- param_list$n_samples
    n_levels <- param_list$n_ploidy + 1
    i_sample <- -seq_len(param_list$n_parents)
    sample_geno <- apply(X = best_geno[, i_sample, ], MARGIN = 2, FUN = colSums)
    sample_geno <- as.vector(sample_geno + 1)
    sample_geno <- sample_geno + seq(0, by = n_levels, length.out = n_mar * n_sample)
    log10_th <- log10(param_list$call_threshold)
    geno_prob <- matrix(data = pat_prob[sample_geno],
                        nrow = n_mar,
                        ncol = n_sample)
    geno_prob <- t(geno_prob < log10_th)
    for(i in seq_len(param_list$n_ploidy)){
        best_geno[i, i_sample, ][geno_prob] <- 3
    }
    if (!param_list$het_parent) {
        if(param_list$n_ploidy %% 2 == 0){
            new_num <- ceiling(best_hap[best_hap != 1] / param_list$n_ploidy)

        } else {
            new_num <- ceiling(best_hap[best_hap != 1] / (param_list$n_ploidy - 1))
        }
        best_hap[best_hap != 1] <- new_num
    }

    for(i in seq_len(param_list$n_ploidy)){
        best_hap[i, i_sample, ][geno_prob] <- 0
    }
    out_list <- list(best_hap = best_hap, best_geno = best_geno)
    return(out_list)
}

################################################################################
################################################################################
# Marker specific error rate estimation
.getBias <- function(best_seq, type, ref, alt, n_ploidy) {
    if (type == 1) {
        est_het <- !best_seq %in% c(0, n_ploidy)
        ref[!est_het] <- NA
        ref <- colSums(ref, na.rm = TRUE)
        best_seq[!est_het] <- NA
        n_ref <- colSums(n_ploidy - best_seq, na.rm = TRUE)
        ref_prop <- ref / n_ref
        alt[!est_het] <- NA
        alt <- colSums(alt, na.rm = TRUE)
        n_alt <- colSums(best_seq, na.rm = TRUE)
        alt_prop <- alt / n_alt
        bias <- ref_prop / (ref_prop + alt_prop)

    } else {
        est_ref <- best_seq == 0
        ref[!est_ref] <- NA
        ref <- colSums(ref, na.rm = TRUE)
        n_ref <- colSums(est_ref, na.rm = TRUE) * n_ploidy
        ref_prop <- ref / n_ref

        est_alt <- best_seq == n_ploidy
        alt[!est_alt] <- NA
        alt <- colSums(alt, na.rm = TRUE)
        n_alt <- colSums(est_alt, na.rm = TRUE) * n_ploidy
        alt_prop <- alt / n_alt
        bias <- ref_prop / (ref_prop + alt_prop)
    }
    return(list(bias = bias,
                ref = ref,
                alt = alt,
                n_ref = n_ref,
                n_alt = n_alt))
}

.sumUpBias <- function(bias1, bias2) {
    ref_prop <- (bias1$ref + bias2$ref) / (bias1$n_ref + bias2$n_ref)
    alt_prop <- (bias1$alt + bias2$alt) / (bias1$n_alt + bias2$n_alt)
    bias <- ref_prop / (ref_prop + alt_prop)
    return(bias)
}

.calcReadBias <- function(best_seq, param_list) {
    ref <- rbind(param_list$reads$p_ref, param_list$reads$ref)
    alt <- rbind(param_list$reads$p_alt, param_list$reads$alt)

    bias1 <- .getBias(best_seq = best_seq, type = 1,
                      ref = ref, alt = alt, n_ploidy = param_list$n_ploidy)
    bias2 <- .getBias(best_seq = best_seq, type = 2,
                      ref = ref, alt = alt, n_ploidy = param_list$n_ploidy)
    bias_cor <- cor(x = bias1$bias, y = bias2$bias, use = "pair")
    if (!is.na(bias_cor) & bias_cor > 0.7) {
        bias <- .sumUpBias(bias1 = bias1, bias2 = bias2)

    } else {
        bias <- bias1$bias
    }
    bias[is.na(bias)] <- 0.5
    return(bias)
}

.calcMissmap <- function(best_seq, param_list) {
    geno_call <- get_genocall(ref = param_list$reads$ref,
                              alt = param_list$reads$alt,
                              eseq_in = param_list$error_rate,
                              bias = param_list$bias,
                              mismap = param_list$mismap,
                              n_o = param_list$n_samples,
                              n_m = param_list$n_mar,
                              ploidy = param_list$n_ploidy)
    missing <- param_list$reads$ref == 0 & param_list$reads$alt == 0
    geno_call[missing] <- FALSE

    i_samples <- -seq_len(param_list$n_parents)
    est <- best_seq[i_samples, ] == 0
    n_ref <- colSums(est, na.rm = TRUE)
    alt <- param_list$reads$alt > 0
    alt[!est] <- FALSE
    alt[!geno_call] <- FALSE
    alt_mis <- colSums(alt, na.rm = TRUE) / n_ref

    est <- best_seq[i_samples, ] == param_list$n_ploidy
    n_alt <- colSums(est, na.rm = TRUE)
    ref <- param_list$reads$ref > 0
    ref[!est] <- FALSE
    ref[!geno_call] <- FALSE
    ref_mis <- colSums(ref, na.rm = TRUE) / n_alt

    return(cbind(alt_mis, ref_mis))
}

.calcErrors <- function(best_seq, param_list) {
    best_seq <- apply(X = best_seq, MARGIN = 3, FUN = colSums)
    bias <- .calcReadBias(best_seq = best_seq, param_list = param_list)
    mismap <- .calcMissmap(best_seq = best_seq, param_list = param_list)
    return(list(bias = bias, mismap = mismap))
}

.bindErrors <- function(error_f, error_r) {
    nmar <- length(error_f$bias)
    bias <- colMeans(rbind(error_f$bias, error_r$bias), na.rm = TRUE)
    bias[is.na(bias)] <- 0.5

    mismap <- cbind(colMeans(rbind(error_f$mismap[, 1],
                                   error_r$mismap[, 1]), na.rm = TRUE),
                    colMeans(rbind(error_f$mismap[, 2],
                                   error_r$mismap[, 2]), na.rm = TRUE))
    mismap[is.na(mismap)] <- 0
    return(list(bias = bias, mismap = mismap))
}

################################################################################
################################################################################
# Flip parameter lists for the reverse-direction round
.flipParam <- function(param_list) {
    n_mar <- param_list$n_mar
    param_list$trans_prob <- lapply(X = param_list$trans_prob,
                                    FUN = function(x){
                                        x$prob <- x$prob[, (n_mar - 1):1]
                                        return(x)
                                    })
    param_list$reads$p_ref <- param_list$reads$p_ref[, n_mar:1]
    param_list$reads$p_alt <- param_list$reads$p_alt[, n_mar:1]
    param_list$reads$ref <- param_list$reads$ref[, n_mar:1]
    param_list$reads$alt <- param_list$reads$alt[, n_mar:1]
    param_list$bias <- param_list$bias[n_mar:1]
    param_list$mismap <- param_list$mismap[n_mar:1,]
    param_list$p_geno_fix <- param_list$p_geno_fix[length(param_list$p_geno_fix):1]
    return(param_list)
}


################################################################################
################################################################################
# Run the iterative estimation
.runCycle <- function(param_list, outprob, outgeno) {
    param_list$count <- param_list$count + 1
    cycle <- paste0("Cycle ", param_list$count, ": ")
    message("\r", cycle)
    message("\r", "Forward round of genotype estimation ...")

    if(!is.null(param_list$fixed_param$parent1)){
        message("\r",
                "Use fixed parental genotypes ...")
    }

    if (param_list$flip) {
        best_pat_r <- .getBestSeq(param_list = .flipParam(param_list),
                                  outprob = outprob)
        best_hap_r <- .pat2hap(best_hap = best_pat_r, param_list = param_list)
        p_geno_r <- .parentPat2Geno(best_hap = best_pat_r, param_list = param_list)
        best_geno_r <- .hap2geno(hap = best_hap_r,
                                 p_geno = p_geno_r,
                                 param_list = param_list)
        last_half <- (round(param_list$n_mar / 2) + 1):param_list$n_mar
        if(is.null(param_list$fixed_param$parent1)){
            param_list$p_geno_fix <- best_pat_r$p_geno[last_half]
        }

    } else {
        best_pat_f <- .getBestSeq(param_list = param_list, outprob = outprob)
        best_hap_f <- .pat2hap(best_hap = best_pat_f, param_list = param_list)
        p_geno_f <- .parentPat2Geno(best_hap = best_pat_f, param_list = param_list)
        best_geno_f <- .hap2geno(hap = best_hap_f,
                                 p_geno = p_geno_f,
                                 param_list = param_list)
        last_half <- (round(param_list$n_mar / 2) + 1):param_list$n_mar
        if(is.null(param_list$fixed_param$parent1)){
            param_list$p_geno_fix <- rev(best_pat_f$p_geno[last_half])
        }
    }

    message("\r", "Backward round of genotype estimation  ...")
    if (param_list$flip) {
        best_pat_f <- .getBestSeq(param_list = param_list, outprob = outprob)
        best_hap_f <- .pat2hap(best_hap = best_pat_f, param_list = param_list)
        p_geno_f <- .parentPat2Geno(best_hap = best_pat_f, param_list = param_list)
        best_geno_f <- .hap2geno(hap = best_hap_f,
                                 p_geno = p_geno_f,
                                 param_list = param_list)
        if(is.null(param_list$fixed_param$parent1)){
            param_list$p_geno_fix <- -1
        }

    } else {
        best_pat_r <- .getBestSeq(param_list = .flipParam(param_list = param_list),
                                  outprob = outprob)
        best_hap_r <- .pat2hap(best_hap = best_pat_r, param_list = param_list)
        p_geno_r <- .parentPat2Geno(best_hap = best_pat_r, param_list = param_list)
        best_geno_r <- .hap2geno(hap = best_hap_r,
                                 p_geno = p_geno_r,
                                 param_list = param_list)
        if(is.null(param_list$fixed_param$parent1)){
            param_list$p_geno_fix <- -1
        }
    }

    # ### Debug
    # return(list(best_geno_f, best_geno_r, param_list))
    # ###

    if (!outprob) {
        if(!is.null(param_list$fixed_param$bias) & !is.null(param_list$fixed_param$mismap_ref)){
            message("\r",
                    "Skip paramter optimization and use fixed parameters ...")

        } else {
            message("\r",
                    "Paramter optimization ...")
            error_f <- .calcErrors(best_seq = best_geno_f, param_list = param_list)
            error_r <- .calcErrors(best_seq = best_geno_r[,,param_list$n_mar:1],
                                   param_list = param_list)
            error <- .bindErrors(error_f = error_f, error_r = error_r)

            if(!is.null(param_list$fixed_param$bias)){
                message("\r",
                        "Use fixed bias values ...")

            } else {
                check <- error$bias > param_list$error_rate[1]
                error$bias[check] <- param_list$error_rate[1]
                check <- error$bias < param_list$error_rate[2]
                error$bias[check] <- param_list$error_rate[2]
                param_list$bias <- error$bias
            }

            if(!is.null(param_list$fixed_param$mismap_ref)){
                message("\r", "Use fixed mismapping rates ...")

            } else {
                param_list$mismap <- error$mismap
            }
        }
    }

    if (outgeno) {
        message("\r", "Summarizing output ...")
        best_pat <- .halfJoint(best_pat_f = best_pat_f,
                               best_pat_r = best_pat_r,
                               param_list = param_list)
        p_geno <- .parentPat2Geno(best_hap = best_pat, param_list = param_list)
        best_hap <- .pat2hap(best_hap = best_pat, param_list = param_list)
        best_hap <- .solveBorderConflict(best_hap = best_hap,
                                         param_list = param_list)
        best_hap <- .minimizeRecombination(best_hap = best_hap)
        best_geno <- .hap2geno(hap = best_hap,
                               p_geno = p_geno,
                               param_list = param_list)
        out_list <- .summarizeEst(best_hap = best_hap,
                                  best_geno = best_geno,
                                  pat_prob = best_pat$prob,
                                  param_list = param_list)
        out_list$pat_prob <- best_pat$prob
        out_list$p_geno <- p_geno
        out_list$bias <- param_list$bias
        out_list$mismap <- param_list$mismap
        out_list$mapping_id <- c(param_list$parents_ids, param_list$sample_ids)

        message("\r", "Done!")
        return(out_list)
    } else {
        return(param_list)
    }
}

################################################################################
################################################################################
# Check data to avoid no reads for founders at the first marker
.checkPread <- function(param_list) {
    param_list$flip <- FALSE
    p_read_s <- sum(param_list$reads$p_ref[, 1] == 0 &
                        param_list$reads$p_alt[, 1] == 0)
    p_read_e <- sum(param_list$reads$p_ref[, param_list$n_mar] == 0 &
                        param_list$reads$p_alt[, param_list$n_mar] == 0)
    if (p_read_s < p_read_e) {
        param_list$flip <- TRUE

    } else if (p_read_s == p_read_e) {
        p_read_s <- sum(param_list$reads$p_ref[, 1] +
                            param_list$reads$p_alt[, 1])
        p_read_e <- sum(param_list$reads$p_ref[, param_list$n_mar] +
                            param_list$reads$p_alt[, param_list$n_mar])
        if (p_read_s < p_read_e) {
            param_list$flip <- TRUE

        } else if (p_read_s == p_read_e) {
            p_read_s <- sum(param_list$reads$p_ref[, 1] +
                                param_list$reads$p_alt[, 1])
            p_read_e <- sum(param_list$reads$p_ref[, param_list$n_mar] +
                                param_list$reads$p_alt[, param_list$n_mar])
            if (p_read_s < p_read_e) {
                param_list$flip <- TRUE
            }
        }
    }
    return(param_list)
}

################################################################################
################################################################################
# Run genotype estimation per chromosome
.cleanEachChr <- function(object,
                          chr_i,
                          node,
                          error_rate,
                          recomb_rate,
                          call_threshold,
                          het_parent,
                          optim,
                          iter,
                          parentless,
                          dummy_reads,
                          pat) {
    param_list <- .getParams(object = object,
                             chr_i = chr_i,
                             node = node,
                             error_rate = error_rate,
                             recomb_rate = recomb_rate,
                             call_threshold = call_threshold,
                             het_parent = het_parent,
                             parentless = parentless,
                             dummy_reads = dummy_reads,
                             pat = pat)
    param_list <- .checkPread(param_list = param_list)

    # ### Debug
    # out <- .runCycle(param_list, FALSE, FALSE)
    # return(out)
    # ###

    if(iter == 1){ optim <- FALSE }

    if(optim){
        param_list <- .runCycle(param_list = param_list,
                                outprob = FALSE,
                                outgeno = FALSE)

        for (i in 2:iter) {
            if (i == iter) {
                out_list <- .runCycle(param_list = param_list,
                                      outprob = TRUE,
                                      outgeno = TRUE)

            } else {
                param_list <- .runCycle(param_list = param_list,
                                        outprob = FALSE,
                                        outgeno = FALSE)
            }
        }
    } else {
        out_list <- .runCycle(param_list = param_list,
                              outprob = TRUE,
                              outgeno = TRUE)
    }

    return(out_list)
}
