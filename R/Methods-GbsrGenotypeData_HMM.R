#' @importFrom SeqArray seqOptimize
#' @importFrom gdsfmt put.attr.gdsn index.gdsn
#' @rdname estGeno
setMethod("estGeno",
          "GbsrGenotypeData",
          function(object,
                   recomb_rate,
                   error_rate,
                   call_threshold,
                   het_parent,
                   optim,
                   iter,
                   n_threads,
                   dummy_reads,
                   fix_bias,
                   fix_mismap) {

              # Check the validity of the registered scheme
              object <- .checkScheme(object = object)

              # Check if the parents have been specified
              parents <- slot(object = slot(object = object, name = "scheme"),
                              name = "parents")
              n_parents <- length(unique(parents))
              parentless <- all(is.na(parents))
              if(parentless){
                  message("Run in the parentless mode...")
              }

              # Set the number of threads
              .setThreads(n_threads = n_threads)

              message("Start cleaning...")

              # initialize the output nodes
              .initGDS(object = object,
                       het_parent = het_parent,
                       n_parents = n_parents)

              # Get chromosome information to loop over chromosomes
              chr <- getChromosome(object = object)
              chr_levels <- unique(chr)
              no_eds <- FALSE
              for(chr_i in chr_levels) {
                  message("\nNow cleaning chr ", chr_i, "...")
                  clean_out <- .cleanEachChr(object = object,
                                             chr_i = chr_i,
                                             error_rate = error_rate,
                                             recomb_rate = recomb_rate,
                                             call_threshold = call_threshold,
                                             het_parent = het_parent,
                                             optim = optim,
                                             iter = iter,
                                             fix_bias = fix_bias,
                                             fix_mismap = fix_mismap,
                                             parentless = parentless,
                                             dummy_reads = dummy_reads)
                  # ### Debug
                  # return(clean_out)
                  # ###

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
                  .saveHap(object = object, clean_out = clean_out, sel = sel)
                  .saveGeno(object = object, clean_out = clean_out, sel = sel)
                  .savePGeno(object = object, clean_out = clean_out, sel = sel)
                  .saveADB(object = object, clean_out = clean_out, sel = sel)
                  .saveMR(object = object, clean_out = clean_out, sel = sel)

                  # Output the dosage data if available
                  if(n_parents == 2 & !het_parent){
                      .saveEDS(object = object, clean_out = clean_out, sel = sel)
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
#' @importFrom gdsfmt append.gdsn
.saveHap <- function(object, clean_out, sel) {
    output <- array(data = 0, dim = c(2, length(sel$sam), length(sel$mar)))
    rep_id <- getReplicates(object = object, parents = TRUE)
    id_hit <- match(x = rep_id, table = clean_out$mapping_id)
    output[, sel$sam, sel$mar] <- clean_out$best_hap[, id_hit,]

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

.saveEDS <- function(object, clean_out, sel) {
    output <- array(data = NA, dim = c(2, length(sel$sam), length(sel$mar)))
    rep_id <- getReplicates(object = object, parents = TRUE)
    id_hit <- match(x = rep_id, table = clean_out$mapping_id)
    output[, sel$sam, sel$mar] <- clean_out$best_hap[, id_hit,]

    output[output == 0] <- NA
    output <- apply(X = output - 1, MARGIN = c(2, 3), FUN = sum)
    output[is.na(output)] <- 8

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
.saveGeno <- function(object, clean_out, sel) {
    output <- array(data = 3, dim = c(2, length(sel$sam), length(sel$mar)))
    rep_id <- getReplicates(object = object, parents = TRUE)
    id_hit <- match(x = rep_id, table = clean_out$mapping_id)
    output[, sel$sam, sel$mar] <- clean_out$best_geno[, id_hit,]

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
.savePGeno <- function(object, clean_out, sel) {
    output <- matrix(data = 3, nrow = nrow(clean_out$p_geno), length(sel$mar))
    output[, sel$mar] <- clean_out$p_geno
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
.saveADB <- function(object, clean_out, sel) {
    output <- rep(-1, length(sel$mar))
    output[sel$mar] <- t(clean_out$bias)
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
.saveMR <- function(object, clean_out, sel) {
    output <- matrix(-1, nrow = ncol(clean_out$mismap), ncol = length(sel$mar))
    output[, sel$mar] <- t(clean_out$mismap)
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
.loadReadCounts <- function(object, chr_i, parentless, dummy_reads) {
    if (exist.gdsn(node = object, path = "annotation/format/FAD")) {
        ad_node <- "filt"

    } else {
        ad_node <- "raw"
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
    }

    p_reads <- .orderParents(object = object, p_reads = p_reads)

    return(list(p_ref = p_reads$ref,
                p_alt = p_reads$alt,
                ref = reads$ref,
                alt = reads$alt))
}

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

.makeGenoParents <- function(n_parents, alleles, het_parent){
    out <- NULL
    for (i in seq_len(n_parents * 2)) {
        out <- c(out, list(alleles))
    }
    out <- expand.grid(out, KEEP.OUT.ATTRS = FALSE)
    valid_pat <- apply(X = out, MARGIN = 1,
                       FUN = function(x) length(unique(x)) > 1)
    out <- as.matrix(out[valid_pat,])
    if (!het_parent) {
        valid <- out[, c(TRUE, FALSE)] == out[, c(FALSE, TRUE)]
        valid <- apply(X = valid, MARGIN = 1, FUN = all)
        out <- out[valid, ]
    }
    attributes(out) <- list(dim = dim(out))
    return(out)
}

.initialPattern <- function(mt, xtype, het_parent, n_origin){
    if (het_parent) {
        out <- apply(X = matrix(data = seq_len(n_origin), 2),
                     MARGIN = 2,
                     FUN = paste, collapse="/")
        out <- lapply(X = seq_along(xtype[[1]]),
                      FUN = function(i){
                          gamet <- out[match(x = mt[[1]][, i],
                                             table = seq_len(0.5 * n_origin))]
                          return(paste(gamet, collapse = "|"))
                      })

    } else {
        gamet <- mt[[1]]
        gamet[gamet != 1] <- gamet[gamet != 1] * 2 - 1
        out <- lapply(X = seq_along(xtype[[1]]),
                      FUN = function(i){
                          return(paste(gamet[, i], collapse = "|"))
                      })
    }
    return(out)
}

.progenyPattern <- function(gamet, mt, xtype, pg, het_parent){
    if(all(xtype == "pairing")){
        out <- vapply(X = seq_along(xtype),
                      FUN.VALUE = list(1),
                      FUN = function(i){
                          i_gamet <- gamet[match(mt[, i], pg)]
                          return(list(paste(i_gamet, collapse = "|")))
                      })

    } else if(all(xtype != "pairing")){
        out <- vapply(X = gamet,
                      FUN.VALUE = list(1),
                      FUN = function(x){
                          return(list(paste(x, x, sep = "|")))
                      })
    }
    return(out)
}

.getValidPat <- function(scheme, het_parent, n_origin) {
    Var1 <- Var2 <- NULL
    homo <- FALSE

    target_pedigree <- sort(unique(slot(object = scheme, name = "samples")))
    mt <- slot(object = scheme, name = "mating")
    xtype <- slot(object = scheme, name = "crosstype")
    pg <- slot(object = scheme, name = "progenies")

    pat <- NULL

    # Pairing of founders
    pat <- c(pat, list(.initialPattern(mt = mt,
                                       xtype = xtype,
                                       het_parent = het_parent,
                                       n_origin = n_origin)))

    # Progeny
    if(length(xtype) >= 2){
        for(i in seq_along(xtype)[-1]) {
            gamet <- lapply(X = pat[[i - 1]],
                            FUN = sub, pattern = "\\|", replacement = "/")
            pat <- c(pat,
                     list(.progenyPattern(gamet = gamet,
                                          mt = mt[[i]],
                                          xtype = xtype[[i]],
                                          pg = pg[[i - 1]],
                                          het_parent = het_parent)))
        }
    }
    target_index <- which(unlist(pg) %in% target_pedigree)
    target_pat <- unlist(pat)[target_index]
    gamet <- lapply(X = target_pat,
                    FUN = function(x) unlist(strsplit(x, split = "\\|")))
    out <- lapply(X = gamet,
                  FUN = function(x){
                      gamet1 <- as.numeric(unlist(strsplit(x[1], "/")))
                      gamet2 <- as.numeric(unlist(strsplit(x[2], "/")))
                      tmp <- expand.grid(gamet1, gamet2, stringsAsFactors = FALSE)
                      tmp <- as.matrix(tmp)
                      tmp <- rbind(tmp, tmp[, 2:1])
                      tmp <- unique(tmp)
                      tmp <- tmp[order(tmp[, 1], tmp[, 2]),]
                      return(tmp)
                  })

    names(out) <- target_pedigree
    return(out)
}

.getPossibleHap <- function(hap_progeny, geno_parents, geno_pat){
    hap_vec <- as.vector(t(hap_progeny))

    derived_geno <- apply(X = geno_parents, MARGIN = 1,
                          FUN = function(x){
                              x <- x[hap_vec]
                              x <- x[c(TRUE, FALSE)] + x[c(FALSE, TRUE)]
                              return(x)
                          })
    derived_geno <- as.vector(derived_geno)
    out <- rbind(derived_geno, derived_geno, derived_geno) == geno_pat
    out <- apply(X = out, MARGIN = 2, FUN = which)
    return(out)
}

.getPossibleGeno <- function(geno_parents, geno_pat){
    out <- apply(X = geno_parents, MARGIN = 1, FUN = function(x) {
        x <- x[c(TRUE, FALSE)] + x[c(FALSE, TRUE)]
        return(x)
    })
    out <- as.vector(out)
    out <- rbind(out, out, out) == geno_pat
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
    geno_parents <- .makeGenoParents(n_parents = n_parents,
                                     alleles = alleles,
                                     het_parent = het_parent)
    n_origin <- n_parents * 2^het_parent
    hap_progeny <- .getValidPat(scheme = scheme,
                                het_parent = het_parent,
                                n_origin = n_origin)

    possiblehap <- lapply(X = hap_progeny,
                          FUN = .getPossibleHap,
                          geno_parents = geno_parents,
                          geno_pat = geno_pat)

    possiblegeno <- .getPossibleGeno(geno_parents = geno_parents,
                                     geno_pat = geno_pat)

    n_p_pat <- nrow(geno_parents)
    n_hap_pat <- vapply(X = hap_progeny, FUN = nrow, FUN.VALUE = numeric(1))
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
            r <- switch(n_pairing, "1" = 0, n_pairing - 2)
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
                            scheme, het_parent){
    out <- lapply(X = seq_along(pat$hap_progeny),
                  FUN = function(i){
                      all_hap <- do.call(what = "rbind", args = pat$hap_progeny)
                      hap_progeny_i <- pat$hap_progeny[[i]]
                      n_origin <- length(unique(as.vector(hap_progeny_i)))
                      jnum <- .getJnum(scheme = scheme,
                                       n_origin = n_origin,
                                       het_parent = het_parent,
                                       i = i)
                      jrate <- .getXoFreq(jnum = jnum, n_origin = n_origin)
                      invalid_joint <- apply(X = hap_progeny_i,
                                             MARGIN = 1,
                                             FUN = function(x) {
                                                 apply(X = hap_progeny_i,
                                                       MARGIN = 1,
                                                       FUN = function(y) {
                                                           check1 <- x[1] == x[2] & y[1] == y[2]
                                                           check2 <- any(x == y)
                                                           return(!check1 & !check2)
                                                       })
                                             })

                      ibd <- apply(X = hap_progeny_i,
                                   MARGIN = 1,
                                   FUN = function(x) abs(length(unique(x)) - 2))
                      ibd <- vapply(X = ibd,
                                    FUN = function(x) paste0(ibd, x),
                                    FUN.VALUE = character(length(ibd)))
                      ibd[ibd == "00"] <- jrate$r00
                      ibd[ibd == "01"] <- jrate$r01
                      ibd[ibd == "10"] <- jrate$r10
                      ibd[ibd == "11"] <- jrate$r11
                      ibd[invalid_joint] <- 0
                      q_mat <- as.numeric(ibd)

                      mar_dist <- diff(pos)
                      rf <- mar_dist * 1e-6 * recomb_rate
                      q_mat <- matrix(q_mat,
                                      nrow(hap_progeny_i),
                                      nrow(hap_progeny_i))
                      diag(q_mat) <- NA
                      diag(q_mat) <- -rowSums(q_mat, na.rm = TRUE)
                      prob <- vapply(X = rf,
                                     FUN = function(x) expm.Higham08(q_mat * x),
                                     FUN.VALUE = numeric(length(q_mat)))
                      prob_dim <- c(pat$n_hap_pat[i], pat$n_hap_pat[i], length(mar_dist))
                      prob <- array(data = prob, dim = prob_dim)
                      return(prob)
                  })
    return(out)
}

.getInitProb <- function(prob, n_pat, n_samples) {
    ev1 <- eigen(t(prob[, , 1]))$vectors[, 1]
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

.getParams <- function(object, chr_i, error_rate, recomb_rate,
                       call_threshold, het_parent, fix_bias, fix_mismap,
                       parentless, dummy_reads) {
    reads <- .loadReadCounts(object = object, chr_i = chr_i,
                             parentless = parentless, dummy_reads = dummy_reads)
    if(parentless){
        n_parents <- 2

    } else {
        parents <- getParents(object = object)
        n_parents <- length(unique(parents$memberID))
    }
    n_samples <- length(unique(getReplicates(object = object)))
    n_alleles <- 2
    n_ploidy <- attributes(slot(object = object, name = "sample"))$ploidy
    pat <- .makePattern(n_parents = n_parents, n_ploidy = n_ploidy,
                        n_alleles = n_alleles, n_samples = n_samples,
                        het_parent = het_parent,
                        scheme = slot(object = object, name = "scheme"))

    pos <- getPosition(object = object, valid = TRUE, chr = chr_i)
    n_mar <- nmar(object = object, valid = TRUE, chr = chr_i)
    trans_prob <- .transitionProb(pat = pat, pos = pos,
                                  recomb_rate = recomb_rate,
                                  scheme = slot(object = object, name = "scheme"),
                                  het_parent = het_parent)
    init_prob <- lapply(X = trans_prob, FUN = .getInitProb,
                        n_pat = pat$n_p_pat, n_samples = n_samples)
    if(is.null(fix_mismap)){
        mismap <-  matrix(data = 0.005, nrow = n_mar, ncol = 2)
        fix_mismap <- FALSE

    } else {
        mismap <- matrix(data = fix_mismap, nrow = n_mar, ncol = 2)
        fix_mismap <- TRUE
    }
    if(is.null(fix_bias)){
        bias <- rep(x = 0.5, times = n_mar)
        fix_bias <- FALSE
    } else {
        bias <- rep(x = fix_bias, times = n_mar)
        fix_bias <- TRUE
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
                fix_mismap = fix_mismap,
                fix_bias = fix_bias,
                trans_prob = lapply(X = trans_prob, FUN = log10),
                init_prob = lapply(X = init_prob, FUN = log10),
                count = 0,
                p_geno_fix = -1,
                flip = FALSE))
}

################################################################################
################################################################################
# Execute the Viterbi algorithm
.getBestSeq <- function(param_list, outprob) {
    trans_prob <- unlist(param_list$trans_prob)
    init_prob <- unlist(param_list$init_prob)
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
                            n_pgeno = param_list$pat$n_p_pat,
                            n_hap = param_list$pat$n_hap_pat,
                            n_offspring = param_list$n_samples,
                            n_founder = param_list$n_parents,
                            n_marker = param_list$n_mar,
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
                       n_pgeno = param_list$pat$n_p_pat,
                       n_hap = param_list$pat$n_hap_pat,
                       n_offspring = param_list$n_samples,
                       n_marker = param_list$n_mar,
                       pedigree = pedigree - 1,
                       p_geno = out_list$p_geno,
                       ploidy = param_list$n_ploidy)
        prob <- array(data = apply(X = prob,
                                   MARGIN = 3,
                                   FUN = function(x) return(t(x))),
                      dim = c(3, param_list$n_samples, param_list$n_mar))
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
        parent_hap <- rep(x = seq(1, by = 2, length.out = param_list$n_parents),
                          each = param_list$n_ploidy * param_list$n_mar)
    }
    parent_hap <- matrix(data = parent_hap, nrow = n_parents_chr, byrow = TRUE)
    sample_pat <- best_hap$best_seq
    pat_i <- match(x = param_list$pedigree,
                   table = names(param_list$pat$hap_progeny))
    sample_hap <- vapply(X = seq_len(param_list$n_samples),
                         FUN.VALUE = numeric(2 * param_list$n_mar),
                         FUN = function(i) {
                             hap_i <- param_list$pat$hap_progeny[[pat_i[i]]]
                             out <- c(hap_i[sample_pat[, i], 1],
                                      hap_i[sample_pat[, i], 2])
                             return(out)
                         })
    sample_hap <- matrix(data = sample_hap,
                         nrow = param_list$n_samples * 2,
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
    prob <- array(data = prob, dim = c(3, n_mar, n_sample))
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
    best_hap[, flip_i, half_m] <- best_hap[2:1, flip_i, half_m]
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
    hap1 <- which(colSums(i_hap_progeny == best_hap[, i, half]) == 2)
    hap2 <- which(colSums(i_hap_progeny == best_hap[, i, half + 1]) == 2)
    re <- rev(seq_along(best_hap[, i, half + 1]))
    hap3 <- which(colSums(i_hap_progeny == best_hap[re, i, half + 1]) == 2)
    prob1 <- param_list$trans_prob[[i_pedigree]][hap1, hap2, half]
    prob2 <- param_list$trans_prob[[i_pedigree]][hap1, hap3, half]
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
        recomb_pos <- which(n_recomb == 1)
        for(j in seq_along(recomb_pos)){
            if(j == 1){
                prev_xo <- recomb[recomb_pos[j], ]

            } else {
                is_homo_j <- diff(best_hap[, i, recomb_pos[j]]) == 0
                if(is_homo_j){
                    j_xo <- recomb[recomb_pos[j], ]
                    if(all(j_xo == prev_xo)){
                        prev_xo <- j_xo[2:1]
                        flip_i <- -seq(1, recomb_pos[j])
                        best_hap[, i, flip_i] <- best_hap[2:1, i, flip_i]
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
    i_sample <- -seq_len(param_list$n_parents)
    sample_geno <- apply(X = best_geno[, i_sample, ], MARGIN = 2, FUN = colSums)
    sample_geno <- as.vector(sample_geno + 1)
    sample_geno <- sample_geno + seq(0, by=3, length.out = n_mar * n_sample)
    log10_th <- log10(param_list$call_threshold)
    geno_prob <- matrix(data = pat_prob[sample_geno],
                        nrow = n_mar,
                        ncol = n_sample)
    geno_prob <- t(geno_prob < log10_th)
    best_geno[1, i_sample, ][geno_prob] <- 3
    best_geno[2, i_sample, ][geno_prob] <- 3
    if (!param_list$het_parent) {
        best_hap[best_hap != 1] <- (best_hap[best_hap != 1] + 1) / 2
    }
    best_hap[1, i_sample, ][geno_prob] <- 0
    best_hap[2, i_sample, ][geno_prob] <- 0
    out_list <- list(best_hap = best_hap, best_geno = best_geno)
    return(out_list)
}

################################################################################
################################################################################
# Marker specific error rate estimation
.getBias <- function(best_seq, type, ref, alt) {
    if (type == 1) {
        est_het <- best_seq == 1
        ref[!est_het] <- NA
        ref <- colSums(ref, na.rm = TRUE)
        n_ref <- colSums(est_het, na.rm = TRUE)
        ref_prop <- ref / n_ref
        alt[!est_het] <- NA
        alt <- colSums(alt, na.rm = TRUE)
        n_alt <- colSums(est_het, na.rm = TRUE)
        alt_prop <- alt / n_alt
        bias <- ref_prop / (ref_prop + alt_prop)

    } else {
        est_ref <- best_seq == 0
        ref[!est_ref] <- NA
        ref <- colSums(ref, na.rm = TRUE)
        n_ref <- colSums(est_ref, na.rm = TRUE)
        ref_prop <- ref / n_ref

        est_alt <- best_seq == 2
        alt[!est_alt] <- NA
        alt <- colSums(alt, na.rm = TRUE)
        n_alt <- colSums(est_alt, na.rm = TRUE)
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
    ref_prop <- (bias1$ref + bias2$ref) / (bias1$n_ref + bias2$n_ref * 2)
    alt_prop <- (bias1$alt + bias2$alt) / (bias1$n_alt + bias2$n_alt * 2)
    bias <- ref_prop / (ref_prop + alt_prop)
    return(bias)
}

.calcReadBias <- function(best_seq, param_list) {
    ref <- rbind(param_list$reads$p_ref, param_list$reads$ref)
    alt <- rbind(param_list$reads$p_alt, param_list$reads$alt)

    bias1 <- .getBias(best_seq = best_seq, type = 1, ref = ref, alt = alt)
    bias2 <- .getBias(best_seq = best_seq, type = 2, ref = ref, alt = alt)
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
                              n_m = param_list$n_mar)
    missing <- param_list$reads$ref == 0 & param_list$reads$alt == 0
    geno_call[missing] <- FALSE

    i_samples <- -seq_len(param_list$n_parents)
    est <- best_seq[i_samples, ] == 0
    n_ref <- colSums(est, na.rm = TRUE)
    alt <- param_list$reads$alt > 0
    alt[!est] <- FALSE
    alt[!geno_call] <- FALSE
    alt_mis <- colSums(alt, na.rm = TRUE) / n_ref

    est <- best_seq[i_samples, ] == 2
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
                                        return(x[, , (n_mar - 1):1])
                                    })
    param_list$reads$p_ref <- param_list$reads$p_ref[, n_mar:1]
    param_list$reads$p_alt <- param_list$reads$p_alt[, n_mar:1]
    param_list$reads$ref <- param_list$reads$ref[, n_mar:1]
    param_list$reads$alt <- param_list$reads$alt[, n_mar:1]
    param_list$bias <- param_list$bias[n_mar:1]
    param_list$mismap <- param_list$mismap[n_mar:1,]
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

    if (param_list$flip) {
        best_pat_r <- .getBestSeq(param_list = .flipParam(param_list),
                                  outprob = outprob)
        best_hap_r <- .pat2hap(best_hap = best_pat_r, param_list = param_list)
        p_geno_r <- .parentPat2Geno(best_hap = best_pat_r, param_list = param_list)
        best_geno_r <- .hap2geno(hap = best_hap_r,
                                 p_geno = p_geno_r,
                                 param_list = param_list)
        last_half <- (round(param_list$n_mar / 2) + 1):param_list$n_mar
        param_list$p_geno_fix <- best_pat_r$p_geno[last_half]

    } else {
        best_pat_f <- .getBestSeq(param_list = param_list, outprob = outprob)
        best_hap_f <- .pat2hap(best_hap = best_pat_f, param_list = param_list)
        p_geno_f <- .parentPat2Geno(best_hap = best_pat_f, param_list = param_list)
        best_geno_f <- .hap2geno(hap = best_hap_f,
                                 p_geno = p_geno_f,
                                 param_list = param_list)
        last_half <- (round(param_list$n_mar / 2) + 1):param_list$n_mar
        param_list$p_geno_fix <- best_pat_f$p_geno[last_half]
    }

    message("\r", "Backward round of genotype estimation  ...")
    if (param_list$flip) {
        best_pat_f <- .getBestSeq(param_list = param_list, outprob = outprob)
        best_hap_f <- .pat2hap(best_hap = best_pat_f, param_list = param_list)
        p_geno_f <- .parentPat2Geno(best_hap = best_pat_f, param_list = param_list)
        best_geno_f <- .hap2geno(hap = best_hap_f,
                                 p_geno = p_geno_f,
                                 param_list = param_list)
        param_list$p_geno_fix <- -1

    } else {
        best_pat_r <- .getBestSeq(param_list = .flipParam(param_list = param_list),
                                  outprob = outprob)
        best_hap_r <- .pat2hap(best_hap = best_pat_r, param_list = param_list)
        p_geno_r <- .parentPat2Geno(best_hap = best_pat_r, param_list = param_list)
        best_geno_r <- .hap2geno(hap = best_hap_r,
                                 p_geno = p_geno_r,
                                 param_list = param_list)
        param_list$p_geno_fix <- -1
    }

    # ### Debug
    # return(list(best_geno_f, best_geno_r, param_list))
    # ###

    if (!outprob) {
        message("\r",
                "Paramter optimization ...")
        error_f <- .calcErrors(best_seq = best_geno_f, param_list = param_list)
        error_r <- .calcErrors(best_seq = best_geno_r[,,param_list$n_mar:1],
                               param_list = param_list)
        error <- .bindErrors(error_f = error_f, error_r = error_r)
        check <- error$bias > param_list$error_rate[1]
        error$bias[check] <- param_list$error_rate[1]
        check <- error$bias < param_list$error_rate[2]
        error$bias[check] <- param_list$error_rate[2]

        if(!param_list$fix_bias){
            param_list$bias <- error$bias
        }

        if(!param_list$fix_mismap){
            param_list$mismap <- error$mismap
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
        }
    }
    return(param_list)
}

################################################################################
################################################################################
# Run genotype estimation per chromosome
.cleanEachChr <- function(object,
                          chr_i,
                          error_rate,
                          recomb_rate,
                          call_threshold,
                          het_parent,
                          optim,
                          iter,
                          fix_bias,
                          fix_mismap,
                          parentless,
                          dummy_reads) {
    param_list <- .getParams(object = object,
                             chr_i = chr_i,
                             error_rate = error_rate,
                             recomb_rate = recomb_rate,
                             call_threshold = call_threshold,
                             het_parent = het_parent,
                             fix_bias = fix_bias,
                             fix_mismap = fix_mismap,
                             parentless = parentless,
                             dummy_reads = dummy_reads)
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
