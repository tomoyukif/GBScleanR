#' @rdname estGeno
setMethod("estGeno",
          "GbsrGenotypeData",
          function(object,
                   chr,
                   recomb_rate,
                   error_rate,
                   call_threshold,
                   het_parent,
                   optim,
                   iter,
                   n_threads,
                   sep_vit){
            max_threads <- RcppParallel::defaultNumThreads()
            if(is.null(n_threads)){
              n_threads <- max_threads / 2
            }
            if(max_threads <= n_threads){
              warning("You are going to use all threads of your computer for the calculation.")
              answer <- ""
              while(answer != "y"){
                answer <- readline(prompt = "Are you sure to use all threads?(y/n)")
                if(answer == "n"){
                  stop("Stopped by user.")
                }
              }
              n_threads <- max_threads
            }
            
            message("Set the number of threads: ", n_threads)
            RcppParallel::setThreadOptions(numThreads = n_threads)
            
            message("Start cleaning...")
            
            .initGDS(object)
            chr_all <- getChromosome(object, levels = TRUE)
            if(missing(chr)){
              chr <- chr_all
            }
            for(chr_i in chr_all){
              index <- getChromosome(object, valid = FALSE) %in% chr_i
              valid_index <- index & getValidSnp(object)
              
              if(!chr_i %in% chr){
                .saveHap(object, NA, sum(index), FALSE)
                .saveGeno(object, NA, sum(index), FALSE)
                .savePGeno(object, NA, sum(index), FALSE)
                message(paste0("Chr ", chr_i, " was skipped."))
                next
              }
              
              message(paste0("\nNow cleaning chr ", chr_i, "..."))
              
              best_seq <- .cleanEachChr(object = object,
                                        index = valid_index,
                                        error_rate = error_rate,
                                        recomb_rate = recomb_rate,
                                        call_threshold = call_threshold,
                                        het_parent = het_parent,
                                        optim = optim,
                                        iter = iter,
                                        sep_vit = sep_vit)
              .saveHap(object, best_seq$best_hap, sum(index), valid_index[index])
              .saveGeno(object, best_seq$best_geno, sum(index), valid_index[index])
              .savePGeno(object, best_seq$p_geno, sum(index), valid_index[index])
            }
            .finalizeGDS(object)
            return(object)
          })

.initGDS <- function(object){
  gdsfmt::add.gdsn(node=object@data@handler,
                   name="estimated.haplotype",
                   storage="bit6",
                   compress="LZMA_RA",
                   replace=TRUE)
  gdsfmt::add.gdsn(node=object@data@handler,
                   name="corrected.genotype",
                   storage="bit2",
                   compress="LZMA_RA",
                   replace=TRUE)
  gdsfmt::add.gdsn(node=object@data@handler,
                   name="parents.genotype",
                   storage="bit2",
                   compress="LZMA_RA",
                   replace=TRUE)
}

.finalizeGDS <- function(object){
  gdsfmt::readmode.gdsn(node=gdsfmt::index.gdsn(node=object@data@handler,
                                                path = "estimated.haplotype"))
  gdsfmt::readmode.gdsn(node=gdsfmt::index.gdsn(node=object@data@handler,
                                                path = "corrected.genotype"))
  gdsfmt::readmode.gdsn(node=gdsfmt::index.gdsn(node=object@data@handler,
                                                path = "parents.genotype"))
}

.loadReadCounts <- function(object, index){
  if(object@data@genotypeVar == "filt.genotype"){
    ad_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                  path="annotation/format/AD/filt.data")
  } else {
    ad_node <- gdsfmt::index.gdsn(node=object@data@handler,
                                  path="annotation/format/AD/data")
  }
  valid_markers <- index
  valid_markers_2 <- rep(valid_markers, each=2)
  valid_scan <- getValidScan(object)
  parents <- object@scanAnnot$parents
  p_index <- which(parents != 0)
  valid_scan[p_index] <- TRUE
  parents_index <- p_index[order(parents[p_index])]
  
  parents <- parents[valid_scan]
  samples_index <- parents == 0
  
  ad <- gdsfmt::readex.gdsn(ad_node, list(valid_scan, valid_markers_2))
  ref <- ad[, c(TRUE, FALSE)]
  alt <- ad[, c(FALSE, TRUE)]
  if(haveFlipped(object)){
    flipped <-getFlipped(object, valid = FALSE)[valid_markers]
    tmp <- ref[, flipped]
    ref[, flipped] <- alt[, flipped]
    alt[, flipped] <- tmp
  }
  
  p_ref <- ref[parents != 0, ]
  p_alt <- alt[parents != 0, ]
  ref <- ref[parents == 0, ]
  alt <- alt[parents == 0, ]
  p_ref <- p_ref[order(parents[parents != 0]), ]
  p_alt <- p_alt[order(parents[parents != 0]), ]
  return(list(p_ref = p_ref,
              p_alt = p_alt,
              ref = ref,
              alt = alt,
              samples_index = samples_index,
              parents_index = parents_index))
}

.getJnum <- function(scheme, het_parent){
  i_pairing <- grep("pairing", scheme@crosstype)
  if(length(i_pairing) == 0){
    n_pairing  <- 0
  } else {
    n_pairing <- max(i_pairing)
  }
  
  n_founder <- scheme@pop_size[1] * 2^het_parent
  if(n_pairing == 0){
    jnum <- .initJnum(n_founder)
  } else {
    next_crosstype <- scheme@crosstype[n_pairing + 1]
    jnum <- .initJnum(n_founder, n_pairing + het_parent, next_crosstype)
  }
  
  s_gen <- n_pairing + 1
  n_gen <- length(scheme@crosstype)
  for(i in s_gen:n_gen){
    jnum <- .calcNextJnum(scheme@crosstype[i], jnum, scheme@pop_size[i])
  }
  return(jnum)
}

.initJnum <- function(n_founder, n_pairing, next_crosstype){
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
  if(n_founder > 2){
    jnum$b123 <- 1
  }
  if(!missing(n_pairing)){
    jnum$a12 <- 1
    if(next_crosstype == "selfing"){
      jnum$j1232 <- jnum$r <- n_pairing - 1
      
    } else {
      jnum$j1232 <- jnum$r <- jnum$k1232 <- switch(n_pairing, "1" = 0, n_pairing - 2)
      if(n_founder > 2){
        jnum$a123 <- 1
      } else {
        jnum$b12 <- 0.5
      }
      if(next_crosstype == "sibling"){
        while(n_pairing > 1){
          jnum <- .calcNextJnum("sibling", jnum, 2)
          n_pairing <- n_pairing - 1
        }
      }
    }
  }
  return(jnum)
}

.getCoalProb <- function(crosstype, pop_size){
  return(
    switch (crosstype,
            "pairing" = c(0, 0),
            "selfing" = c(1, 0),
            "sibling" = c(0.5, 0),
            "random" = rep(0.5 / pop_size, 2),
    ))
}

.calcNextJnum <- function(crosstype, prob_df, pop_size){
  coal <- .getCoalProb(crosstype, pop_size)
  s <- coal[1]
  q <- coal[2]
  if(crosstype == "selfing"){
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
    next_b123 <- 1.5 * q * prob_df$a123 + (1 - s - 2 * q) * prob_df$b123
    next_r <- prob_df$r + prob_df$a12
    next_j1232 <- prob_df$k1232 + prob_df$a123
    next_k1232 <- 0.5 * s * prob_df$j1232 + (1 - s) * (prob_df$k1232 + prob_df$a123)
    next_j1122 <- prob_df$k1122
    next_k1122 <- 0.5 * s * (prob_df$r + prob_df$j1122) + (1 - s) * prob_df$k1122
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

.getXoFreq <- function(jnum, n_origin){
  r00 <- jnum$j1232 / jnum$a12 / (n_origin - 2)
  r01 <- jnum$j1222 / jnum$a12
  r10 <- jnum$j1222 / (1 - jnum$a12) / (n_origin - 1)
  r11 <- jnum$j1122 / (1 - jnum$a12) / (n_origin - 1)
  return(data.frame(r00, r01, r10, r11))
}

.getValidPat <- function(scheme, het_parent, n_origin){
  Var1 <- Var2 <- NULL
  homo <- FALSE
  
  if(het_parent){
    pat <- apply(matrix(1:n_origin, nrow = 2), 2, paste, collapse = "/")
    chr <- pat[match(scheme@mating[[1]], 1:(0.5*n_origin))]
    if(scheme@crosstype[1] == "pairing"){
      pat <- apply(matrix(chr, dim(scheme@mating[[1]])), 2, paste, collapse = "|")
      
    } else if(scheme@crosstype[1] == "random"){
      mating <- t(expand.grid(chr, chr))
      mating <- mating[,apply(mating, 2, function(x) length(unique(x)) != 1)]
      pat <- apply(mating, 2, paste, collapse = "|")
    }
    
  } else {
    chr <- scheme@mating[[1]]
    chr[chr != 1] <- chr[chr != 1] * 2 - 1
    if(scheme@crosstype[1] == "pairing"){
      pat <- apply(chr, 2, paste, collapse = "|")
      
    } else if(scheme@crosstype[1] == "random"){
      mating <- t(expand.grid(chr, chr))
      mating <- mating[,apply(mating, 2, function(x) length(unique(x)) != 1)]
      pat <- apply(mating, 2, paste, collapse = "|")
    }
  }
  
  for(i in 2:length(scheme@crosstype)){
    if(scheme@crosstype[i] == "pairing"){
      chr <- sub("\\|", "/", pat)
      mating <- matrix(chr[match(scheme@mating[[i]], scheme@progenies[[i - 1]])],
                       dim(scheme@mating[[i]]))
      pat <- apply(mating, 2, paste, collapse = "|")
      
    } else if(scheme@crosstype[i] == "sibling"){
      chr <- strsplit(pat, "\\|")[[1]]
      pat <- apply(t(expand.grid(chr, chr)), 2, paste, collapse = "|")
      homo <- TRUE
      
    } else if(scheme@crosstype[i] == "random"){
      chr <- sub("\\|", "/", pat)
      mating <- t(expand.grid(chr, chr))
      mating <- mating[,apply(mating, 2, function(x) length(unique(x)) != 1)]
      pat <- apply(mating, 2, paste, collapse = "|")
      homo <- TRUE
    } else if(scheme@crosstype[i] == "selfing"){
      homo <- TRUE
    }
  }
  chr <- strsplit(pat, "\\|")
  pat <- NULL
  for(i in 1:length(chr)){
    chr1 <- strsplit(chr[[i]][1], "/")[[1]]
    chr2 <- strsplit(chr[[i]][2], "/")[[1]]
    pat <- rbind(pat, expand.grid(chr1, chr2, stringsAsFactors = F))
  }
  pat <- subset(pat, subset = Var1 != Var2)
  pat <- as.matrix(pat)
  pat <- rbind(pat, cbind(pat[, 2], pat[, 1]))
  
  if(homo){
    pat <- rbind(pat, t(matrix(rep(unique(as.vector(pat)), each = 2), 2)))
  }
  pat_dim <- dim(pat)
  pat <- matrix(as.numeric(pat), pat_dim[1], pat_dim[2])
  pat <- pat[order(pat[,1], pat[,2]), ]
  pat <- subset(pat, subset = !duplicated(pat))
  return(pat)
}

.transitionProb <- function(pat, pos, recomb_rate, scheme, n_origin, het_parent){
  jrate <- .getXoFreq(.getJnum(scheme, het_parent), n_origin)
  invalid_joint <- apply(pat$hap_progeny, 1, function(x){
    apply(pat$hap_progeny, 1, function(y){
      check1 <- x[1] == x[2] & y[1] == y[2]
      check2 <- any(x == y)
      return(!check1 & !check2)
    })
  })
  
  ibd <- apply(pat$hap_progeny, 1, function(x)abs(length(unique(x)) - 2))
  ibd <- sapply(ibd, function(x) paste0(ibd, x))
  ibd[ibd == "00"] <- jrate$r00
  ibd[ibd == "01"] <- jrate$r01
  ibd[ibd == "10"] <- jrate$r10
  ibd[ibd == "11"] <- jrate$r11
  ibd[invalid_joint] <- 0
  q_mat <- as.numeric(ibd)
  
  snp_dist <- diff(pos)
  rf <- snp_dist * 1e-6 * recomb_rate
  q_mat <- matrix(q_mat, nrow(pat$hap_progeny), nrow(pat$hap_progeny))
  diag(q_mat) <- NA
  diag(q_mat) <- -rowSums(q_mat, na.rm = T)
  prob <- sapply(rf, function(x){
    x <- expm::expm.Higham08(q_mat*x)
    return(x)
  })
  prob_dim <- c(pat$n_hap_pat, pat$n_hap_pat, length(snp_dist))
  
  prob <- array(prob, prob_dim)
  return(prob)
}

.traceback <- function(x){
  n_t <- ncol(x)
  x <- x[, n_t:1]
  best_i <- x[1, 1]
  for(i in 2:n_t){
    best_i <- c(best_i, x[best_i[i - 1], i])
  }
  best_i <- best_i[n_t:1]
  return(best_i)
}

.addParents <- function(best_hap, p_geno, param_list){
  best_hap <- cbind(best_hap,
                    matrix(rep(p_geno, param_list$n_parents),
                           ncol = param_list$n_parents))
  index <- c(which(param_list$samples_index), param_list$parents_index)
  best_hap[, index] <- best_hap
  return(best_hap)
}

.saveHap <- function(object, best_hap, n_snp, valid_index){
  if(is.array(best_hap)){
    output <- array(0, c(2, n_snp, nscan(object, valid = FALSE)))
    output[, valid_index, getValidScan(object, parents = TRUE)] <- best_hap
    output <- t(matrix(output, n_snp * 2))
    output[is.na(output)] <- 0
    
  } else {
    output <- array(0, c(2, n_snp, nscan(object, valid = FALSE)))
  }
  
  hap_gdsn <- gdsfmt::index.gdsn(object@data@handler, "estimated.haplotype")
  gdsn_dim <- gdsfmt::objdesp.gdsn(hap_gdsn)$dim
  if(gdsn_dim[1] == 0){
    gdsfmt::add.gdsn(node=object@data@handler,
                     name="estimated.haplotype",
                     storage="bit6",
                     compress="LZMA_RA",
                     val = output,
                     replace=TRUE)
  } else {
    gdsfmt::append.gdsn(hap_gdsn, val = output)
  }
}

.parentPat2Geno <- function(best_hap, param_list){
  p_pat <- best_hap$p_geno
  p_geno <- t(param_list$pat$geno_parents[p_pat, ])
  return(p_geno)
}

.hap2geno <- function(hap, p_geno, param_list){
  out_geno <- apply(hap, 3, function(x){
    sapply(1:ncol(x), function(y){
      return(p_geno[x[1, y], y] + p_geno[x[2, y], y])
    })
  })
  sample_order <- c(param_list$parents_index, which(param_list$samples_index))
  out_geno[, sample_order] <- out_geno
  return(out_geno)
}

.pat2hap <- function(best_hap, param_list){
  n_parents_chr <- (param_list$n_parents * param_list$n_ploidy)
  if(param_list$het_parent){
    parent_hap <- rep(1:n_parents_chr,
                      each = param_list$n_snp)
  } else {
    parent_hap <- rep(seq(1, by = 2, length.out = param_list$n_parents),
                      each = param_list$n_ploidy * param_list$n_snp)
  }
  parent_hap <- matrix(parent_hap, nrow = n_parents_chr, byrow = TRUE)
  sample_pat <- best_hap$best_seq
  sample_hap <- apply(sample_pat, 1, function(x){
    return(rbind(param_list$pat$hap_progeny[x, 1],
                 param_list$pat$hap_progeny[x, 2]))
  })
  out_hap <- rbind(parent_hap, sample_hap)
  out_hap <- array(out_hap,
                   c(param_list$n_ploidy,
                     length(param_list$samples_index),
                     param_list$n_snp))
  out_hap <- apply(out_hap, 2, function(x)return(x))
  out_hap <- array(out_hap,
                   c(param_list$n_ploidy,
                     param_list$n_snp,
                     length(param_list$samples_index)))
  return(out_hap)
}

.saveGeno <- function(object, best_geno, n_snp, valid_index){
  if(is.na(best_geno[1])){
    output <- matrix(3, n_snp, nscan(object, valid = FALSE))
    
  } else {
    output <- matrix(3, n_snp, nscan(object, valid = FALSE))
    output[valid_index, getValidScan(object, parents = TRUE)] <- best_geno
    output[is.na(output)] <- 3
  }
  
  out_gdsn <- gdsfmt::index.gdsn(object@data@handler, "corrected.genotype")
  gdsn_dim <- gdsfmt::objdesp.gdsn(out_gdsn)$dim
  if(gdsn_dim[1] == 0){
    gdsfmt::add.gdsn(node=object@data@handler,
                     name="corrected.genotype",
                     storage="bit2",
                     compress="LZMA_RA",
                     val = t(output),
                     replace=TRUE)
  } else {
    gdsfmt::append.gdsn(out_gdsn, val = t(output))
  }
}

.savePGeno <- function(object, p_geno, n_snp, valid_index){
  if(is.na(p_geno[1])){
    n_p <- nrow(getParents(object))
    output <- matrix(3, n_p * 2, n_snp)
  } else {
    output <- matrix(3, nrow(p_geno), n_snp)
    output[, valid_index] <- p_geno
  }
  out_gdsn <- gdsfmt::index.gdsn(object@data@handler, "parents.genotype")
  gdsn_dim <- gdsfmt::objdesp.gdsn(out_gdsn)$dim
  if(gdsn_dim[1] == 0){
    gdsfmt::add.gdsn(node=object@data@handler,
                     name="parents.genotype",
                     storage="bit2",
                     compress="LZMA_RA",
                     val = output,
                     replace=TRUE)
  } else {
    gdsfmt::append.gdsn(out_gdsn, val = output)
  }
}

.solveVit <- function(vit_list, param_list){
  last_path <- matrix(rep(apply(vit_list$vit_s, 2, which.max),
                          each = param_list$pat$n_hap_pat),
                      param_list$pat$n_hap_pat)
  vit_list$vit_path[,,param_list$n_snp] <-last_path
  best_hap <- apply(vit_list$vit_path, 2, .traceback)
  best_hap <- .addParents(best_hap = best_hap,
                          p_geno = vit_list$p_geno,
                          param_list = param_list)
  
  return(best_hap)
}

.getInitProb <- function(prob, n_pat, n_samples){
  ev1 <- eigen(t(prob))$vectors[, 1]
  init <- ev1 / sum(ev1)
  return(init)
}

.getParams <- function(object,
                       index,
                       error_rate,
                       recomb_rate,
                       call_threshold,
                       het_parent){
  
  reads <- .loadReadCounts(object, index)
  
  n_parents <- nrow(reads$p_ref)
  n_samples <- nrow(reads$ref)
  n_alleles <- 2
  n_ploidy <- getPloidy(object, valid = FALSE)[index][1]
  n_origin <- n_parents * 2^het_parent
  pat <- .makePattern(n_parents,
                      n_ploidy,
                      n_alleles,
                      n_samples,
                      het_parent,
                      scheme = object@scheme,
                      n_origin)
  
  pos <- getPosition(object, valid = FALSE)[index]
  n_snp <- length(pos)
  trans_prob <- .transitionProb(pat, pos, recomb_rate,
                                scheme = object@scheme, n_origin, het_parent)
  init_prob <- .getInitProb(trans_prob[,, 1], pat$n_p_pat, n_samples)
  
  return(list(n_parents = n_parents,
              n_samples = n_samples,
              n_alleles = n_alleles,
              n_ploidy = n_ploidy,
              n_snp = n_snp,
              pos = pos,
              samples_index = reads$samples_index,
              parents_index = reads$parents_index,
              het_parent = het_parent,
              error_rate = c(1 - error_rate, error_rate),
              recomb_rate = recomb_rate,
              call_threshold = call_threshold,
              reads = reads[grepl("ref$|alt$", names(reads))],
              pat = pat,
              bias = rep(0.5, n_snp),
              mismap = matrix(0, n_snp, 2),
              trans_prob = log10(trans_prob),
              init_prob = log10(init_prob),
              count = 0,
              p_geno_fix = -1,
              flip = FALSE))
}


.makePattern <- function(n_parents,
                         n_ploidy,
                         n_alleles,
                         n_samples,
                         het_parent,
                         scheme,
                         n_origin){
  alleles <- seq(0, length.out = n_alleles)
  
  geno_pat <- NULL
  for(i in 1:n_ploidy){
    geno_pat <- c(geno_pat, list(alleles))
  }
  geno_pat <- expand.grid(geno_pat)
  geno_pat <- t(apply(geno_pat, 1, sort))
  geno_pat <- rowSums(geno_pat[!duplicated(geno_pat), ])
  
  geno_parents <- NULL
  for(i in 1:(n_parents * 2)){
    geno_parents <- c(geno_parents, list(alleles))
  }
  geno_parents <- expand.grid(geno_parents, KEEP.OUT.ATTRS = FALSE)
  valid_pat <- apply(geno_parents, 1, function(x){length(unique(x)) > 1})
  geno_parents <- as.matrix(geno_parents[valid_pat, ])
  if(!het_parent){
    valid <- apply(geno_parents[, c(T, F)] == geno_parents[, c(F, T)], 1, all)
    geno_parents <- geno_parents[valid, ]
  }
  attributes(geno_parents) <- list(dim = dim(geno_parents))
  
  hap_progeny <- .getValidPat(scheme = scheme, het_parent, n_origin)
  
  hap_vec <- as.vector(t(hap_progeny))
  
  derived_geno <- apply(geno_parents, 1, function(x){
    x <- x[hap_vec]
    x <- x[c(T, F)] + x[c(F, T)]
    return(x)
  })
  derived_geno <- as.vector(derived_geno)
  possiblehap <- rbind(derived_geno, derived_geno, derived_geno) == geno_pat
  possiblehap <- apply(possiblehap, 2, which)
  
  possiblegeno <- apply(geno_parents, 1, function(x){
    x <- x[c(T, F)] + x[c(F, T)]
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

.getBestSeq <- function(param_list, outprob){
  param_list$trans_prob <- matrix(param_list$trans_prob,
                                  nrow = dim(param_list$trans_prob)[1])
  
  
  if(param_list$sep_vit){
    
    out_list <- run_viterbi2(p_ref = param_list$reads$p_ref,
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
                             n_m = param_list$n_snp,
                             possiblehap = param_list$pat$possiblehap - 1,
                             possiblegeno = param_list$pat$possiblegeno - 1,
                             p_geno_fix = param_list$p_geno_fix - 1)
  } else {
    
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
                             n_m = param_list$n_snp,
                             possiblehap = param_list$pat$possiblehap - 1,
                             possiblegeno = param_list$pat$possiblegeno - 1,
                             p_geno_fix = param_list$p_geno_fix - 1)
  }
  
  if(outprob){
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
                   n_m = param_list$n_snp,
                   possiblehap = param_list$pat$possiblehap - 1)
    prob <- array(apply(prob, 3, function(x) return(t(x))),
                  dim = c(3, param_list$n_samples, param_list$n_snp))
    out_list$prob <- prob
  } else {
    out_list$prob <- NULL
  }
  out_list$best_seq <- out_list$best_seq + 1
  out_list$p_geno <- out_list$p_geno + 1
  return(out_list)
}

.flipParam <- function(param_list){
  n_snp <- param_list$n_snp
  param_list$trans_prob <- param_list$trans_prob[,, (n_snp - 1):1]
  param_list$reads$p_ref <- param_list$reads$p_ref[, n_snp:1]
  param_list$reads$p_alt <- param_list$reads$p_alt[, n_snp:1]
  param_list$reads$ref <- param_list$reads$ref[, n_snp:1]
  param_list$reads$alt <- param_list$reads$alt[, n_snp:1]
  param_list$bias <- param_list$bias[n_snp:1]
  param_list$mismap <- param_list$mismap[n_snp:1, ]
  return(param_list)
}

.halfJoint <- function(best_pat_f, best_pat_r, n_snp, n_sample){
  half <- round(n_snp / 2)
  first <- (n_snp:1)[1:half]
  latter <- (half + 1):n_snp
  
  best_seq <- rbind(best_pat_r$best_seq[first, ],
                    best_pat_f$best_seq[latter, ])
  p_geno <- c(best_pat_r$p_geno[first], best_pat_f$p_geno[latter])
  
  prob <- c(best_pat_r$prob[,, first],
            best_pat_f$prob[,, latter])
  prob <- array(prob, dim(best_pat_f$prob))
  prob <- apply(prob, 2, function(x)return(x))
  prob <- array(prob, c(3, n_snp, n_sample))
  out_list <- list(best_seq = best_seq, p_geno = p_geno, prob = prob)
  return(out_list)
}

.getBias <- function(best_seq, type, ref, alt){
  if(type == 1){
    est_het <- best_seq == 1
    ref[!est_het] <- NA
    ref <- rowSums(ref, na.rm = TRUE)
    n_ref <- rowSums(est_het, na.rm = TRUE)
    ref_prop <- ref / n_ref
    alt[!est_het] <- NA
    alt <- rowSums(alt, na.rm = TRUE)
    n_alt <- rowSums(est_het, na.rm = TRUE)
    alt_prop <- alt / n_alt
    bias <- ref_prop / (ref_prop + alt_prop)
    
  } else {
    est_ref <- best_seq == 0
    ref[!est_ref] <- NA
    ref <- rowSums(ref, na.rm = TRUE)
    n_ref <- rowSums(est_ref, na.rm = TRUE)
    ref_prop <- ref / n_ref
    
    est_alt <- best_seq == 2
    alt[!est_alt] <- NA
    alt <- rowSums(alt, na.rm = TRUE)
    n_alt <- rowSums(est_alt, na.rm = TRUE)
    alt_prop <- alt / n_alt
    bias <- ref_prop / (ref_prop + alt_prop)
  }
  return(list(bias = bias, ref = ref, alt = alt, n_ref = n_ref, n_alt = n_alt))
}

.sumUpBias <- function(bias1, bias2){
  ref_prop <- (bias1$ref + bias2$ref) / (bias1$n_ref + bias2$n_ref * 2)
  alt_prop <- (bias1$alt + bias2$alt) / (bias1$n_alt + bias2$n_alt * 2)
  bias <- ref_prop / (ref_prop + alt_prop)
  return(bias)
}

.checkCompleHomo <- function(best_seq){
  apply(best_seq == 1, 1, sum)
}

.calcReadBias <- function(best_seq, param_list){
  indices <- c(param_list$parents_index, which(param_list$samples_index))
  ref <- t(rbind(param_list$reads$p_ref, param_list$reads$ref))
  alt <- t(rbind(param_list$reads$p_alt, param_list$reads$alt))
  ref[, indices] <- ref
  alt[, indices] <- alt
  
  bias1 <- .getBias(best_seq, 1, ref, alt)
  bias2 <- .getBias(best_seq, 2, ref, alt)
  bias_cor <- suppressWarnings(cor(bias1$bias, bias2$bias, "pair"))
  if(!is.na(bias_cor) & bias_cor > 0.7){
    bias <- .sumUpBias(bias1, bias2)
  } else {
    bias <- bias1$bias
  }
  bias[is.na(bias)] <- 0.5
  return(bias)
}

.calcMissmap <- function(best_seq, param_list){
  est <- best_seq[, param_list$samples_index] == 0
  n_ref <- rowSums(est, na.rm = TRUE)
  alt <- t(param_list$reads$alt) > 0
  alt[!est] <- NA
  alt_mis <- rowSums(alt, na.rm = TRUE) / n_ref
  
  est <- best_seq[, param_list$samples_index] == 2
  n_alt <- rowSums(est, na.rm = TRUE)
  ref <- t(param_list$reads$ref) > 0
  ref[!est] <- NA
  ref_mis <- apply(ref, 1, sum, na.rm = TRUE) / n_alt
  return(cbind(alt_mis, ref_mis))
}

.calcErrors <- function(best_seq, param_list){
  bias <- .calcReadBias(best_seq, param_list)
  mismap <- .calcMissmap(best_seq, param_list)
  return(list(bias = bias, mismap = mismap))
}

.bindErrors <- function(error_f, error_r){
  nsnp <- length(error_f$bias)
  bias <- colMeans(rbind(error_f$bias,
                         error_r$bias), na.rm = TRUE)
  bias[is.na(bias)] <- 0.5
  mismap <- cbind(colMeans(rbind(error_f$mismap[, 1],
                                 error_r$mismap[, 1]), na.rm = TRUE),
                  colMeans(rbind(error_f$mismap[, 2],
                                 error_r$mismap[, 2]), na.rm = TRUE))
  mismap[is.na(mismap)] <- 0
  return(list(bias = bias, mismap = mismap))
}

.summarizeEst <- function(best_hap, best_geno, pat_prob, param_list){
  n_snp <- param_list$n_snp
  n_sample <- param_list$n_samples
  sample_geno <- best_geno[, param_list$samples_index] + 1
  sample_geno <- sample_geno + seq(from = 0, by = 3, length.out = n_snp * n_sample)
  geno_prob <- matrix(pat_prob[sample_geno], n_snp, n_sample) < log10(param_list$call_threshold)
  best_geno[, param_list$samples_index][geno_prob] <- NA
  best_geno <- abs(best_geno - 2)
  best_geno[is.na(best_geno)] <- 3
  if(!param_list$het_parent){
    best_hap[best_hap != 1] <- (best_hap[best_hap != 1] + 1) / 2
  }
  best_hap[1, , param_list$samples_index][geno_prob] <- 0
  best_hap[2, , param_list$samples_index][geno_prob] <- 0
  out_list <- list(best_hap = best_hap, best_geno = best_geno)
  return(out_list)
}

.runCycle <- function(param_list, outprob, outgeno){
  param_list$count <- param_list$count + 1
  message("Cycle ", param_list$count, "...")
  cycle <- paste0("Cycle ", param_list$count, ": ")
  message("\r", paste0(cycle, " Forward path..."))
  
  if(param_list$flip){
    best_pat_r <- .getBestSeq(.flipParam(param_list), outprob)
    best_hap_r <- .pat2hap(best_pat_r, param_list)
    p_geno_r <- .parentPat2Geno(best_pat_r, param_list)
    best_geno_r <- .hap2geno(best_hap_r, p_geno_r, param_list)
    param_list$p_geno_fix <- best_pat_r$p_geno[param_list$n_snp]
    
  } else {
    best_pat_f <- .getBestSeq(param_list, outprob)
    best_hap_f <- .pat2hap(best_pat_f, param_list)
    p_geno_f <- .parentPat2Geno(best_pat_f, param_list)
    best_geno_f <- .hap2geno(best_hap_f, p_geno_f, param_list)
    param_list$p_geno_fix <- best_pat_f$p_geno[param_list$n_snp]
  }
  
  message("\r", paste0(cycle, "Backward path..."))
  if(param_list$flip){
    best_pat_f <- .getBestSeq(param_list, outprob)
    best_hap_f <- .pat2hap(best_pat_f, param_list)
    p_geno_f <- .parentPat2Geno(best_pat_f, param_list)
    best_geno_f <- .hap2geno(best_hap_f, p_geno_f, param_list)
    param_list$p_geno_fix <- best_pat_f$p_geno[param_list$n_snp]
    
  } else {
    best_pat_r <- .getBestSeq(.flipParam(param_list), outprob)
    best_hap_r <- .pat2hap(best_pat_r, param_list)
    p_geno_r <- .parentPat2Geno(best_pat_r, param_list)
    best_geno_r <- .hap2geno(best_hap_r, p_geno_r, param_list)
    param_list$p_geno_fix <- best_pat_r$p_geno[param_list$n_snp]
  }
  
  if(!outprob){
    message("\r", paste0(cycle, "Estimating error pattern..."))
    error_f <- .calcErrors(best_geno_f, param_list)
    error_r <- .calcErrors(best_geno_r[param_list$n_snp:1, ], param_list)
    error <- .bindErrors(error_f, error_r)
    param_list$bias <- error$bias
    param_list$mismap <- error$mismap
  }
  
  if(outgeno){
    message("\r", "Summarizing output...")
    best_pat <- .halfJoint(best_pat_f,
                           best_pat_r,
                           param_list$n_snp,
                           param_list$n_samples)
    best_hap <- .pat2hap(best_pat, param_list)
    p_geno <- .parentPat2Geno(best_pat, param_list)
    best_geno <- .hap2geno(best_hap, p_geno, param_list)
    out_list <- .summarizeEst(best_hap, best_geno, best_pat$prob, param_list)
    out_list$p_geno <- p_geno
    
    message("\r", "Done!")
    return(out_list)
  } else {
    return(param_list)
  }
}

.checkPread <- function(param_list){
  param_list$flip = FALSE
  p_read_s <- sum(param_list$reads$p_ref[, 1] == 0 & param_list$reads$p_alt[, 1] == 0)
  p_read_e <- sum(param_list$reads$p_ref[, param_list$n_snp] == 0 & param_list$reads$p_alt[, param_list$n_snp] == 0)
  if(p_read_s < p_read_e){
    param_list$flip = TRUE
  } else if(p_read_s == p_read_e){
    p_read_s <- sum(param_list$reads$p_ref[, 1] + param_list$reads$p_alt[, 1])
    p_read_e <- sum(param_list$reads$p_ref[, param_list$n_snp] + param_list$reads$p_alt[, param_list$n_snp])
    if(p_read_s < p_read_e){
      param_list$flip = TRUE
    }
  }
  return(param_list)
}

.cleanEachChr <- function(object,
                          index,
                          error_rate,
                          recomb_rate,
                          call_threshold,
                          het_parent,
                          optim,
                          iter,
                          sep_vit){
  
  param_list <- .getParams(object,
                           index,
                           error_rate,
                           recomb_rate,
                           call_threshold,
                           het_parent)
  param_list <- .checkPread(param_list)
  param_list$sep_vit <- sep_vit
  
  if(iter == 1){
    optim <- FALSE
  }
  
  if(optim){
    param_list <- .runCycle(param_list, FALSE, FALSE)
    
    for(i in 2:iter){
      
      if(i == iter){
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


# # dir <- "~/02_gbscleanr/ForManuscript/simpop_8way_RIL_noADbias_homoParents"
# # dir <- "~/02_gbscleanr/ForManuscript/simpop_2way_F2_noADbias_homoParents"
# dir <- "~/02_gbscleanr/ForManuscript/simpop_2way_F2_noADbias_hetParents_sibling/"
# 
# ############
# library(GBScleanR)
# files <- list.files(dir, ".vcf")
# files <- grep("_LB.vcf", files, invert = T, value = T)
# file_i<-files[42]
# vcf_fn <- paste0(dir, "/", file_i)
# gds_fn <- sub(".vcf", "_gbsr.gds", vcf_fn)
# # gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
# gds <- loadGDS(gds_fn)
# parents <- grep("Founder", getScanID(gds), value = TRUE)
# gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
# 
# # gds <- initScheme(gds, "pairing", matrix(1:8, 2))
# # gds <- addScheme(gds, "pairing", matrix(9:12, 2))
# # gds <- addScheme(gds, "pairing", matrix(13:14, 2))
# # gds <- addScheme(gds, "self")
# # gds <- addScheme(gds, "self")
# # gds <- addScheme(gds, "self")
# # gds <- addScheme(gds, "self")
# # gds <- addScheme(gds, "self")
# 
# # gds <- initScheme(gds, "pair", matrix(1:2, 2))
# # gds <- addScheme(gds, "self")
# 
# gds <- initScheme(gds, "pairing", matrix(1:2, 2))
# gds <- addScheme(gds, "sibling")
# 
# object = gds
# recomb_rate = 0.04
# call_threshold = 0.9
# error_rate = 0.0025
# iter = 4
# het_parent = TRUE
# # het_parent = FALSE
# optim = TRUE
# n_threads = NULL
# index <- getChromosome(object, valid = FALSE) %in% 1
# valid_index <- index & getValidSnp(object)
# param_list <- .getParams(object,
#                          valid_index,
#                          error_rate,
#                          recomb_rate,
#                          call_threshold,
#                          het_parent)
# 
# param_list <- .checkPread(param_list)
# 
# # param_list <- .flipParam(param_list)
# 
# param_list$trans_prob <- matrix(param_list$trans_prob,
#                                 nrow = dim(param_list$trans_prob)[1])
# 
# Rcpp::sourceCpp("src/GBSR_HMM.cpp")
# out_list <- run_viterbi2(p_ref = param_list$reads$p_ref,
#                         p_alt = param_list$reads$p_alt,
#                         ref = param_list$reads$ref,
#                         alt = param_list$reads$alt,
#                         eseq_in = param_list$error_rate,
#                         bias = param_list$bias,
#                         mismap = param_list$mismap,
#                         trans_prob = param_list$trans_prob,
#                         init_prob = param_list$init_prob,
#                         n_p = param_list$pat$n_p_pat,
#                         n_h = param_list$pat$n_hap_pat,
#                         n_f = param_list$n_parents,
#                         n_o = param_list$n_samples,
#                         n_m = param_list$n_snp,
#                         possiblehap = param_list$pat$possiblehap - 1,
#                         possiblegeno = param_list$pat$possiblegeno - 1,
#                         p_geno_fix = param_list$p_geno_fix - 1)
# out_list$p_geno
# prob <- run_fb(ref = param_list$reads$ref,
#                alt = param_list$reads$alt,
#                eseq_in = param_list$error_rate,
#                bias = param_list$bias,
#                mismap = param_list$mismap,
#                trans_prob = param_list$trans_prob,
#                init_prob = param_list$init_prob,
#                p_geno = out_list$p_geno,
#                n_h = param_list$pat$n_hap_pat,
#                n_o = param_list$n_samples,
#                n_m = param_list$n_snp,
#                possiblehap = param_list$pat$possiblehap - 1)
# out_list$prob <- prob
# 
# 
# out_list <- .runCycle(param_list, TRUE, TRUE)
# geno <- out_list$best_geno[, -(1:2)]
# geno[geno == 3] <- NA
# geno<-abs(geno -2)
# geno <- t(geno)
# df <- compareGeno(geno, true)
# df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
# df_ind_mean
# 
# geno <- best_geno_f[, -(1:2)]
