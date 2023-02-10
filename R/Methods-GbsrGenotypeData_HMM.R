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
                   fix_mismap) {

              if(length(slot(slot(object, "scheme"), "crosstype")) == 0){
                  stop("No scheme information.",
                       "\n",
                       "Build scheme information with ",
                       "initScheme() and addScheme().")
              }

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

              message("Start cleaning...")

              .initGDS(object, het_parent)
              chr <- getChromosome(object)
              chr_levels <- unique(chr)
              for (chr_i in chr_levels) {
                  message("\nNow cleaning chr ", chr_i, "...")
                  best_seq <- .cleanEachChr(object, chr_i, error_rate,
                                            recomb_rate, call_threshold,
                                            het_parent, optim, iter, fix_mismap)

                  sel <- list(mar = validMar(object, chr_i),
                              sam = validSam(object, parents = TRUE))
                  .saveHap(object, best_seq$best_hap, sel)
                  .saveGeno(object, best_seq$best_geno, sel)
                  .savePGeno(object, best_seq$p_geno, sel)
                  .saveADB(object, t(best_seq$bias), sel)
                  .saveMR(object, t(best_seq$mismap), sel)
              }
              closeGDS(object, verbose = FALSE)
              seqOptimize(object$filename, "by.sample",
                          c("HAP", "CGT"), verbose = FALSE)
              object <- reopenGDS(object)
              return(object)
          })

.initGDS <- function(object, het_parent) {
    hap <- addfolder.gdsn(index.gdsn(object, "annotation/format"), "HAP",
                          replace = TRUE)
    add.gdsn(hap, "data", storage = "bit6", compress = "", replace = TRUE)

    cgt <- addfolder.gdsn(index.gdsn(object, "annotation/format"), "CGT",
                          replace = TRUE)
    add.gdsn(cgt, "data", storage = "bit2", compress = "", replace = TRUE)

    pgt <- add.gdsn(index.gdsn(object, "annotation/info"), "PGT",
                    storage = "bit2", compress = "", replace = TRUE)

    adb <- add.gdsn(index.gdsn(object, "annotation/info"), "ADB",
                    storage = "single", compress = "", replace = TRUE)

    mr <- add.gdsn(index.gdsn(object, "annotation/info"), "MR",
                   storage = "single", compress = "", replace = TRUE)
}

#' @importFrom gdsfmt append.gdsn
.saveHap <- function(object, best_hap, sel) {
    output <- array(0, c(2, length(sel$sam), length(sel$mar)))
    i_sample <- c(which(slot(object, "sample")$parents != 0), which(validSam(object)))
    output[, sel$sam, sel$mar][, i_sample,] <- best_hap

    hap_gdsn <- index.gdsn(object, "annotation/format/HAP/data")
    gdsn_dim <- objdesp.gdsn(hap_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(index.gdsn(object, "annotation/format/HAP"), "data", output,
                 "bit6", compress = "", replace = TRUE)
    } else {
        append.gdsn(hap_gdsn, output)
    }
}

#' @importFrom gdsfmt append.gdsn
.saveGeno <- function(object, best_geno, sel) {
    output <- array(3, c(2, length(sel$sam), length(sel$mar)))
    i_sample <- c(which(slot(object, "sample")$parents != 0), which(validSam(object)))
    output[, sel$sam, sel$mar][, i_sample,] <- best_geno

    out_gdsn <-index.gdsn(object, "annotation/format/CGT/data")
    gdsn_dim <- objdesp.gdsn(out_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(index.gdsn(object, "annotation/format/CGT"), "data", output,
                 "bit2", compress = "", replace = TRUE)
    } else {
        append.gdsn(out_gdsn, output)
    }
}

#' @importFrom gdsfmt append.gdsn
.savePGeno <- function(object, p_geno, sel) {
    output <- matrix(3, nrow(p_geno), length(sel$mar))
    output[, sel$mar] <- p_geno
    out_gdsn <- index.gdsn(object, "annotation/info/PGT")
    gdsn_dim <- objdesp.gdsn(out_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(index.gdsn(object, "annotation/info"), "PGT", output,
                 "bit2", compress = "", replace = TRUE)
    } else {
        append.gdsn(out_gdsn, output)
    }
}

#' @importFrom gdsfmt append.gdsn
.saveADB <- function(object, bias, sel) {
    output <- rep(-1, length(sel$mar))
    output[sel$mar] <- bias
    out_gdsn <- index.gdsn(object, "annotation/info/ADB")
    gdsn_dim <- objdesp.gdsn(out_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(index.gdsn(object, "annotation/info"), "ADB", output,
                 "single", compress = "", replace = TRUE)
    } else {
        append.gdsn(out_gdsn, output)
    }
}

#' @importFrom gdsfmt append.gdsn
.saveMR <- function(object, mismap, sel) {
    output <- matrix(-1, nrow(mismap), length(sel$mar))
    output[, sel$mar] <- mismap
    out_gdsn <- index.gdsn(object, "annotation/info/MR")
    gdsn_dim <- objdesp.gdsn(out_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(index.gdsn(object, "annotation/info"), "MR", output,
                 "single", compress = "", replace = TRUE)
    } else {
        append.gdsn(out_gdsn, output)
    }
}

.loadReadCounts <- function(object, chr_i) {
    if (exist.gdsn(object, "annotation/format/FAD")) {
        ad_node <- "filt"
    } else {
        ad_node <- "raw"
    }
    reads <- getRead(object, ad_node, FALSE, TRUE, chr_i)
    p_reads <- getRead(object, ad_node, "only", TRUE, chr_i)

    return(list(p_ref = p_reads$ref,
                p_alt = p_reads$alt,
                ref = reads$ref,
                alt = reads$alt))
}

.getJnum <- function(scheme, het_parent) {
    xtype <- slot(scheme, "crosstype")
    i_pairing <- grep("pairing", xtype)
    if (length(i_pairing) == 0) {
        n_pairing  <- 0
    } else {
        n_pairing <- max(i_pairing)
    }

    ps <- slot(scheme, "pop_size")
    n_founder <- ps[1] * 2^het_parent
    if (n_pairing == 0) {
        jnum <- .initJnum(n_founder)
    } else {
        next_crosstype <- xtype[n_pairing + 1]
        if(is.na(next_crosstype)){
            next_crosstype <- "selfing"
        }
        jnum <- .initJnum(n_founder, n_pairing + het_parent, next_crosstype)
    }

    s_gen <- n_pairing + 1
    n_gen <- length(xtype)
    if(s_gen <= n_gen){
        for (i in s_gen:n_gen) {
            jnum <- .calcNextJnum(xtype[i], jnum, ps[i])
        }
    }
    return(jnum)
}

.initJnum <- function(n_founder, n_pairing, next_crosstype) {
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
    if (n_founder > 2) {
        jnum$b123 <- 1
    }
    if (!missing(n_pairing)) {
        jnum$a12 <- 1
        if (next_crosstype == "selfing") {
            jnum$j1232 <- jnum$r <- n_pairing - 1

        } else {
            r <- switch(n_pairing, "1" = 0, n_pairing - 2)
            jnum$j1232 <- jnum$k1232 <- jnum$r <- r

            if (n_founder > 2) {
                jnum$a123 <- 1
            } else {
                jnum$b12 <- 0.5
            }
            if (next_crosstype == "sibling") {
                while (n_pairing > 1) {
                    jnum <- .calcNextJnum("sibling", jnum, 2)
                    n_pairing <- n_pairing - 1
                }
            }
        }
    }
    return(jnum)
}

.getCoalProb <- function(crosstype, pop_size) {
    return(switch (crosstype,
                   "pairing" = c(0, 0),
                   "selfing" = c(1, 0),
                   "sibling" = c(0.5, 0),
                   "random" = rep(0.5 / pop_size, 2)))
}

.calcNextJnum <- function(crosstype, prob_df, pop_size) {
    coal <- .getCoalProb(crosstype, pop_size)
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

.getValidPat <- function(scheme, het_parent, n_origin) {
    Var1 <- Var2 <- NULL
    homo <- FALSE

    mt <- slot(scheme, "mating")
    xtype <- slot(scheme, "crosstype")
    pg <- slot(scheme, "progenies")

    if (het_parent) {
        pat <- apply(matrix(seq_len(n_origin), 2), 2, paste, collapse="/")
        chr <- pat[match(mt[[1]], seq_len(0.5 * n_origin))]
        if (xtype[1] == "pairing") {
            pat <- apply(matrix(chr, dim(mt[[1]])), 2, paste, collapse = "|")

        } else if (xtype[1] == "random") {
            mating <- t(expand.grid(chr, chr))
            mating <- mating[, apply(mating, 2,
                                     function(x) length(unique(x)) != 1)]
            pat <- apply(mating, 2, paste, collapse="|")
        }

    } else {
        chr <- mt[[1]]
        chr[chr != 1] <- chr[chr != 1] * 2 - 1
        if (xtype[1] == "pairing") {
            pat <- apply(chr, 2, paste, collapse="|")

        } else if (xtype[1] == "random") {
            mating <- t(expand.grid(chr, chr))
            mating <- mating[, apply(mating, 2,
                                     function(x) length(unique(x)) != 1)]
            pat <- apply(mating, 2, paste, collapse="|")
        }
    }

    if(length(xtype) >= 2){
        for (i in 2:length(xtype)) {
            if (xtype[i] == "pairing") {
                chr <- sub("\\|", "/", pat)
                mating <- matrix(chr[match(mt[[i]], pg[[i - 1]])], dim(mt[[i]]))
                pat <- apply(mating, 2, paste, collapse="|")

            } else if (xtype[i] == "sibling") {
                chr <- strsplit(pat, "\\|")[[1]]
                pat <- apply(t(expand.grid(chr, chr)), 2, paste, collapse="|")
                homo <- TRUE

            } else if (xtype[i] == "random") {
                chr <- sub("\\|", "/", pat)
                mating <- t(expand.grid(chr, chr))
                mating <- mating[, apply(mating, 2,
                                         function(x) length(unique(x)) != 1)]
                pat <- apply(mating, 2, paste, collapse="|")
                homo <- TRUE
            } else if (xtype[i] == "selfing") {
                homo <- TRUE
            }
        }
    }
    chr <- strsplit(pat, "\\|")
    pat <- NULL
    pat <- do.call("rbind", lapply(chr, function(x){
        chr1 <- strsplit(x[1], "/")[[1]]
        chr2 <- strsplit(x[2], "/")[[1]]
        return(expand.grid(chr1, chr2, stringsAsFactors = FALSE))
    }))
    pat <- subset(pat, subset=Var1 != Var2)
    pat <- as.matrix(pat)
    pat <- rbind(pat, cbind(pat[, 2], pat[, 1]))

    if (homo) {
        pat <- rbind(pat, t(matrix(rep(unique(as.vector(pat)), each=2), 2)))
    }
    pat_dim <- dim(pat)
    pat <- matrix(as.numeric(pat), pat_dim[1], pat_dim[2])
    pat <- pat[order(pat[, 1], pat[, 2]),]
    pat <- subset(pat, subset=!duplicated(pat))
    return(pat)
}

.transitionProb <- function(pat, pos, recomb_rate,
                            scheme, n_origin, het_parent){
    jrate <- .getXoFreq(.getJnum(scheme, het_parent), n_origin)
    invalid_joint <- apply(pat$hap_progeny, 1, function(x) {
        apply(pat$hap_progeny, 1, function(y) {
            check1 <- x[1] == x[2] & y[1] == y[2]
            check2 <- any(x == y)
            return(!check1 & !check2)
        })
    })

    ibd <- apply(pat$hap_progeny, 1, function(x) abs(length(unique(x)) - 2))
    ibd <- vapply(ibd, function(x) paste0(ibd, x), character(length(ibd)))
    ibd[ibd == "00"] <- jrate$r00
    ibd[ibd == "01"] <- jrate$r01
    ibd[ibd == "10"] <- jrate$r10
    ibd[ibd == "11"] <- jrate$r11
    ibd[invalid_joint] <- 0
    q_mat <- as.numeric(ibd)

    mar_dist <- diff(pos)
    rf <- mar_dist * 1e-6 * recomb_rate
    q_mat <- matrix(q_mat, nrow(pat$hap_progeny), nrow(pat$hap_progeny))
    diag(q_mat) <- NA
    diag(q_mat) <- -rowSums(q_mat, na.rm = TRUE)
    prob <- vapply(rf, function(x) expm.Higham08(q_mat * x),
                   numeric(length(q_mat)))
    prob_dim <- c(pat$n_hap_pat, pat$n_hap_pat, length(mar_dist))

    prob <- array(prob, prob_dim)
    return(prob)
}

.getInitProb <- function(prob, n_pat, n_samples) {
    ev1 <- eigen(t(prob))$vectors[, 1]
    init <- ev1 / sum(ev1)
    return(init)
}

.makePattern <- function(n_parents,
                         n_ploidy,
                         n_alleles,
                         n_samples,
                         het_parent,
                         scheme,
                         n_origin) {
    alleles <- seq(0, length.out = n_alleles)

    geno_pat <- NULL
    for (i in seq_len(n_ploidy)) {
        geno_pat <- c(geno_pat, list(alleles))
    }
    geno_pat <- expand.grid(geno_pat)
    geno_pat <- t(apply(geno_pat, 1, sort))
    geno_pat <- rowSums(geno_pat[!duplicated(geno_pat),])

    geno_parents <- NULL
    for (i in seq_len(n_parents * 2)) {
        geno_parents <- c(geno_parents, list(alleles))
    }
    geno_parents <- expand.grid(geno_parents, KEEP.OUT.ATTRS = FALSE)
    valid_pat <- apply(geno_parents, 1, function(x) length(unique(x)) > 1)
    geno_parents <- as.matrix(geno_parents[valid_pat,])
    if (!het_parent) {
        valid <- geno_parents[, c(TRUE, FALSE)] == geno_parents[, c(FALSE, TRUE)]
        valid <- apply(valid, 1, all)
        geno_parents <- geno_parents[valid,]
    }
    attributes(geno_parents) <- list(dim = dim(geno_parents))

    hap_progeny <- .getValidPat(scheme, het_parent, n_origin)

    hap_vec <- as.vector(t(hap_progeny))

    derived_geno <- apply(geno_parents, 1, function(x) {
        x <- x[hap_vec]
        x <- x[c(TRUE, FALSE)] + x[c(FALSE, TRUE)]
        return(x)
    })
    derived_geno <- as.vector(derived_geno)
    possiblehap <- rbind(derived_geno, derived_geno, derived_geno) == geno_pat
    possiblehap <- apply(possiblehap, 2, which)

    possiblegeno <- apply(geno_parents, 1, function(x) {
        x <- x[c(TRUE, FALSE)] + x[c(FALSE, TRUE)]
        return(x)
    })
    possiblegeno <- as.vector(possiblegeno)
    possiblegeno <- rbind(possiblegeno, possiblegeno, possiblegeno) == geno_pat
    possiblegeno <- apply(possiblegeno, 2, which)

    n_p_pat <- nrow(geno_parents)
    n_hap_pat <- nrow(hap_progeny)
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

.getParams <- function(object, chr_i, error_rate, recomb_rate,
                       call_threshold, het_parent, fix_mismap) {
    reads <- .loadReadCounts(object, chr_i)
    parents <- getParents(object)
    n_parents <- nrow(parents)
    n_samples <- nsam(object)
    n_alleles <- 2
    n_ploidy <- attributes(slot(object, "sample"))$ploidy
    n_origin <- n_parents * 2^het_parent
    pat <- .makePattern(n_parents, n_ploidy, n_alleles, n_samples,
                        het_parent, slot(object, "scheme"), n_origin)

    pos <- getPosition(object, TRUE, chr_i)
    n_mar <- nmar(object, TRUE, chr_i)
    trans_prob <- .transitionProb(pat, pos, recomb_rate, slot(object, "scheme"),
                                  n_origin, het_parent)
    init_prob <- .getInitProb(trans_prob[, , 1], pat$n_p_pat, n_samples)
    if(is.null(fix_mismap)){
        mismap = matrix(0.005, n_mar, 2)
        fix_mismap <- FALSE
    } else {
        mismap = matrix(fix_mismap, n_mar, 2)
        fix_mismap <- TRUE
    }

    return(list(n_parents = n_parents,
                n_samples = n_samples,
                n_alleles = n_alleles,
                n_ploidy = n_ploidy,
                n_mar = n_mar,
                samples_index = validSam(object),
                parents_index = parents$indexes,
                het_parent = het_parent,
                error_rate = c(1 - error_rate, error_rate),
                recomb_rate = recomb_rate,
                call_threshold = call_threshold,
                reads = reads,
                pat = pat,
                bias = rep(0.5, n_mar),
                mismap = mismap,
                fix_mismap = fix_mismap,
                trans_prob = log10(trans_prob),
                init_prob = log10(init_prob),
                count = 0,
                p_geno_fix = -1,
                flip = FALSE))
}

.getBestSeq <- function(param_list, outprob) {
    param_list$trans_prob <- matrix(param_list$trans_prob,
                                    dim(param_list$trans_prob)[1])

    out_list <- run_viterbi(p_ref = param_list$reads$p_ref,
                            p_alt = param_list$reads$p_alt,
                            ref = param_list$reads$ref,
                            alt = param_list$reads$alt,
                            eseq_in = param_list$error_rate,
                            bias = param_list$bias,
                            mismap = param_list$mismap,
                            trans_prob = param_list$trans_prob,
                            init_prob = param_list$init_prob,
                            n_p = param_list$pat$n_p_pat,
                            n_h = param_list$pat$n_hap_pat,
                            n_f = param_list$n_parents,
                            n_o = param_list$n_samples,
                            n_m = param_list$n_mar,
                            het = param_list$het_parent,
                            possiblehap = param_list$pat$possiblehap - 1,
                            possiblegeno = param_list$pat$possiblegeno - 1,
                            p_geno_fix = param_list$p_geno_fix - 1,
                            ploidy = param_list$n_ploidy)

    if (outprob) {
        prob <- run_fb(ref = param_list$reads$ref,
                       alt = param_list$reads$alt,
                       eseq_in = param_list$error_rate,
                       bias = param_list$bias,
                       mismap = param_list$mismap,
                       trans_prob = param_list$trans_prob,
                       init_prob = param_list$init_prob,
                       p_geno = out_list$p_geno,
                       n_h = param_list$pat$n_hap_pat,
                       n_o = param_list$n_samples,
                       n_m = param_list$n_mar,
                       possiblehap = param_list$pat$possiblehap - 1,
                       ploidy = param_list$n_ploidy)
        prob <- array(apply(prob, 3, function(x) return(t(x))),
                      dim = c(3, param_list$n_samples, param_list$n_mar))
        out_list$prob <- prob
    } else {
        out_list$prob <- NULL
    }
    out_list$best_seq <- out_list$best_seq + 1
    out_list$p_geno <- out_list$p_geno + 1
    return(out_list)
}

.parentPat2Geno <- function(best_hap, param_list) {
    p_pat <- best_hap$p_geno
    p_geno <- t(param_list$pat$geno_parents[p_pat, ])
    return(p_geno)
}

.hap2geno <- function(hap, p_geno, param_list) {
    out_geno <- sapply(seq_len(param_list$n_mar), function(x) {
        return(p_geno[hap[,,x], x])
    })
    out_geno <- array(out_geno, c(param_list$n_ploidy,
                                  param_list$n_samples + param_list$n_parents,
                                  param_list$n_mar))

    return(out_geno)
}

.pat2hap <- function(best_hap, param_list) {
    n_parents_chr <- (param_list$n_parents * param_list$n_ploidy)
    if (param_list$het_parent) {
        parent_hap <- rep(seq_len(n_parents_chr), each=param_list$n_mar)
    } else {
        parent_hap <- rep(seq(1, by=2, length.out=param_list$n_parents),
                          each=param_list$n_ploidy * param_list$n_mar)
    }
    parent_hap <- matrix(parent_hap, n_parents_chr, byrow=TRUE)
    sample_pat <- best_hap$best_seq
    sample_hap <- apply(sample_pat, 1, function(x) {
        return(rbind(param_list$pat$hap_progeny[x, 1],
                     param_list$pat$hap_progeny[x, 2]))
    })
    out_hap <- rbind(parent_hap, sample_hap)
    out_hap <- array(out_hap, c(param_list$n_ploidy,
                                param_list$n_parents + param_list$n_samples,
                                param_list$n_mar))
    return(out_hap)
}

.halfJoint <- function(best_pat_f, best_pat_r, param_list) {
    n_mar <- param_list$n_mar
    n_sample <- param_list$n_samples
    half <- round(n_mar / 2)
    first <- (n_mar:1)[seq_len(half)]
    latter <- (half + 1):n_mar

    best_seq <- rbind(best_pat_r$best_seq[first,],
                      best_pat_f$best_seq[latter,])
    border_geno <- cbind(best_seq[half, ], best_seq[half + 1, ])
    diff_geno <- border_geno[, 1] != border_geno[, 2]
    diff_seq <- subset(best_seq, select = diff_geno)
    if(any(diff_geno)){
        border_geno1 <- border_geno[diff_geno, 1]
        border_geno2 <- border_geno[diff_geno, 2]
        diff_border_geno1 <- param_list$pat$hap_progeny[border_geno1, ]
        diff_border_geno2 <- param_list$pat$hap_progeny[border_geno2, ]
        if(is.matrix(diff_border_geno1)){
            check_flip <- diff_border_geno1 == diff_border_geno2[, 2:1]
            fliped_geno <- which(rowSums(check_flip) == 2)
        } else {
            check_flip <- diff_border_geno1 == diff_border_geno2[2:1]
            fliped_geno <- which(sum(check_flip) == 2)
        }

        if(length(fliped_geno) > 0){
            fliped_seq <- sapply(seq_along(fliped_geno), function(i){
                g1 <- diff_seq[, fliped_geno[i]][half]
                g2 <- diff_seq[, fliped_geno[i]][half + 1]
                g2_pos <- which(diff_seq[, fliped_geno[i]][latter] == g2)
                check <- which(diff(g2_pos) != 1)
                if(length(check) == 0){
                    diff_seq[, fliped_geno[i]][latter][g2_pos] <- g1
                } else {
                    diff_seq[, fliped_geno[i]][latter][seq_len(min(check))] <- g1
                }
                return(diff_seq[, fliped_geno[i]])
            })

            if(is.matrix(diff_border_geno1)){
                best_seq[, diff_geno][, fliped_geno] <- fliped_seq
            } else {
                best_seq[, diff_geno] <- fliped_seq
            }
        }
    }

    p_geno <- c(best_pat_r$p_geno[first], best_pat_f$p_geno[latter])

    prob <- c(best_pat_r$prob[, , first],
              best_pat_f$prob[, , latter])
    prob <- array(prob, dim(best_pat_f$prob))
    prob <- apply(prob, 2, function(x) return(x))
    prob <- array(prob, c(3, n_mar, n_sample))
    out_list <- list(best_seq = best_seq,
                     p_geno = p_geno,
                     prob = prob)
    return(out_list)
}

.summarizeEst <- function(best_hap, best_geno, pat_prob, param_list) {
    n_mar <- param_list$n_mar
    n_sample <- param_list$n_samples
    i_sample <- -seq_len(param_list$n_parents)
    sample_geno <- apply(best_geno[, i_sample, ], 2, colSums)
    sample_geno <- as.vector(sample_geno + 1)
    sample_geno <- sample_geno + seq(0, by=3, length.out=n_mar * n_sample)
    log10_th <- log10(param_list$call_threshold)
    geno_prob <- t(matrix(pat_prob[sample_geno], n_mar, n_sample) < log10_th)
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

    bias1 <- .getBias(best_seq, 1, ref, alt)
    bias2 <- .getBias(best_seq, 2, ref, alt)
    bias_cor <- cor(bias1$bias, bias2$bias, "pair")
    if (!is.na(bias_cor) & bias_cor > 0.7) {
        bias <- .sumUpBias(bias1, bias2)
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
    best_seq <- apply(best_seq, 3, colSums)
    bias <- .calcReadBias(best_seq, param_list)
    mismap <- .calcMissmap(best_seq, param_list)
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

.flipParam <- function(param_list) {
    n_mar <- param_list$n_mar
    param_list$trans_prob <- param_list$trans_prob[, , (n_mar - 1):1]
    param_list$reads$p_ref <- param_list$reads$p_ref[, n_mar:1]
    param_list$reads$p_alt <- param_list$reads$p_alt[, n_mar:1]
    param_list$reads$ref <- param_list$reads$ref[, n_mar:1]
    param_list$reads$alt <- param_list$reads$alt[, n_mar:1]
    param_list$bias <- param_list$bias[n_mar:1]
    param_list$mismap <- param_list$mismap[n_mar:1,]
    return(param_list)
}

.runCycle <- function(param_list, outprob, outgeno) {
    param_list$count <- param_list$count + 1
    cycle <- paste0("Cycle ", param_list$count, ": ")
    message("\r", cycle, "Forward path...")

    if (param_list$flip) {
        best_pat_r <- .getBestSeq(.flipParam(param_list), outprob)
        best_hap_r <- .pat2hap(best_pat_r, param_list)
        p_geno_r <- .parentPat2Geno(best_pat_r, param_list)
        best_geno_r <- .hap2geno(best_hap_r, p_geno_r, param_list)
        param_list$p_geno_fix <- best_pat_r$p_geno[param_list$n_mar]

    } else {
        best_pat_f <- .getBestSeq(param_list, outprob)
        best_hap_f <- .pat2hap(best_pat_f, param_list)
        p_geno_f <- .parentPat2Geno(best_pat_f, param_list)
        best_geno_f <- .hap2geno(best_hap_f, p_geno_f, param_list)
        param_list$p_geno_fix <- best_pat_f$p_geno[param_list$n_mar]
    }

    message("\r", cycle, "Backward path...")
    if (param_list$flip) {
        best_pat_f <- .getBestSeq(param_list, outprob)
        best_hap_f <- .pat2hap(best_pat_f, param_list)
        p_geno_f <- .parentPat2Geno(best_pat_f, param_list)
        best_geno_f <- .hap2geno(best_hap_f, p_geno_f, param_list)
        param_list$p_geno_fix <- -1

    } else {
        best_pat_r <- .getBestSeq(.flipParam(param_list), outprob)
        best_hap_r <- .pat2hap(best_pat_r, param_list)
        p_geno_r <- .parentPat2Geno(best_pat_r, param_list)
        best_geno_r <- .hap2geno(best_hap_r, p_geno_r, param_list)
        param_list$p_geno_fix <- -1
    }

    if (!outprob) {
        message("\r", cycle,
                "Estimating allele read bias and mismapping pattern...")
        error_f <- .calcErrors(best_geno_f, param_list)
        error_r <- .calcErrors(best_geno_r[,,param_list$n_mar:1], param_list)
        error <- .bindErrors(error_f, error_r)
        check <- error$bias > param_list$error_rate[1]
        error$bias[check] <- param_list$error_rate[1]
        check <- error$bias < param_list$error_rate[2]
        error$bias[check] <- param_list$error_rate[2]
        param_list$bias <- error$bias

        if(!param_list$fix_mismap){
            param_list$mismap <- error$mismap
        }
    }

    if (outgeno) {
        message("\r", "Summarizing output...")
        best_pat <- .halfJoint(best_pat_f, best_pat_r, param_list)
        best_hap <- .pat2hap(best_pat, param_list)
        p_geno <- .parentPat2Geno(best_pat, param_list)
        best_geno <- .hap2geno(best_hap, p_geno, param_list)
        out_list <- .summarizeEst(best_hap, best_geno,
                                  best_pat$prob, param_list)
        out_list$p_geno <- p_geno
        out_list$bias <- param_list$bias
        out_list$mismap <-param_list$mismap

        message("\r", "Done!")
        return(out_list)
    } else {
        return(param_list)
    }
}

.checkPread <- function(param_list) {
    param_list$flip = FALSE
    p_read_s <- sum(param_list$reads$p_ref[, 1] == 0 &
                        param_list$reads$p_alt[, 1] == 0)
    p_read_e <- sum(param_list$reads$p_ref[, param_list$n_mar] == 0 &
                        param_list$reads$p_alt[, param_list$n_mar] == 0)
    if (p_read_s < p_read_e) {
        param_list$flip = TRUE
    } else if (p_read_s == p_read_e) {
        p_read_s <- sum(param_list$reads$p_ref[, 1] +
                            param_list$reads$p_alt[, 1])
        p_read_e <- sum(param_list$reads$p_ref[, param_list$n_mar] +
                            param_list$reads$p_alt[, param_list$n_mar])
        if (p_read_s < p_read_e) {
            param_list$flip = TRUE
        }
    }
    return(param_list)
}

.cleanEachChr <- function(object,
                          chr_i,
                          error_rate,
                          recomb_rate,
                          call_threshold,
                          het_parent,
                          optim,
                          iter,
                          fix_mismap) {
    param_list <- .getParams(object, chr_i, error_rate, recomb_rate,
                             call_threshold, het_parent, fix_mismap)
    param_list <- .checkPread(param_list)

    if (iter == 1) { optim <- FALSE }

    if (optim) {
        param_list <- .runCycle(param_list, FALSE, FALSE)

        for (i in 2:iter) {
            if (i == iter) {
                out_list <- .runCycle(param_list, TRUE, TRUE)

            } else {
                param_list <- .runCycle(param_list, FALSE, FALSE)
            }
        }
    } else {
        out_list <- .runCycle(param_list, TRUE, TRUE)
    }
    return(out_list)
}
