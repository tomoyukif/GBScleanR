#' @param parents Indices of parental lines.
#' @rdname initScheme
#'
setMethod("initScheme",
          "GbsrScheme",
          function(object, mating, parents) {
              if (!all(parents %in% mating)) {
                  stop('All parents should be listed in the matrix ",
                       "specified as the "mating" argument.',
                       call. = FALSE)
              }
              if (!all(mating %in% parents)) {
                  stop('Invalid member ID(s) specified to mating.',
                       call. = FALSE)
              }
              object@parents <- parents
              object@crosstype <- list(rep("pairing", ncol(mating)))
              object@mating <- list(mating)
              object@progenies <- list(seq(max(parents) + 1, by = 1,
                                           length.out = ncol(mating)))
              object@samples <- numeric()
              return(object)
          })

.validMating <- function(crosstype, mating){
    if(crosstype == "pairing"){
        if(any(mating[1] == mating[2])){
            stop('The same member ID was specified to ',
                 'mating for "crosstype = pairing"',
                 call. = FALSE)
        }
    }
    if(crosstype == "sibling"){
        if(any(mating[1] != mating[2])){
            stop('Different member IDs were specified to ',
                 'mating for "crosstype = sibling"',
                 call. = FALSE)
        }
    }
    if(crosstype == "selfing"){
        if(any(mating[1] != mating[2])){
            stop('Different member IDs were specified to ',
                 'mating for "crosstype = selfing"',
                 call. = FALSE)
        }
    }
}

#'
#' @importFrom utils tail
#' @rdname addScheme
#'
setMethod("addScheme",
          "GbsrScheme",
          function(object, crosstype, mating) {
              crosstype <- sapply(crosstype, match.arg,
                                  choices = c("pairing", "selfing", "sibling"))

              last_gen <- unlist(tail(object@progenies, 1))
              n_crosstype <- length(crosstype)
              if(n_crosstype != 1){
                  stop("length(crosstype) should be 1.", call. = FALSE)
              }
              is_na_mating <- any(is.na(mating))

              if(n_crosstype == 1){
                  if(is_na_mating){
                      if (crosstype == "pairing") {
                          mating <- t(expand.grid(last_gen, last_gen))
                          mating <- unique(t(apply(mating, 2, sort)))
                          mating <- t(mating[mating[, 1] != mating[, 2], ])
                          message("As `mating` was not specified,",
                                  " set the following mating design.")
                          print(mating)
                          crosstype <- rep(crosstype, ncol(mating))

                      } else {
                          mating <- matrix(rep(last_gen, each = 2), 2)
                          message("As `mating` was not specified,",
                                  " set the following mating design.")
                          print(mating)
                          crosstype <- rep(crosstype, ncol(mating))
                      }
                  } else {
                      crosstype <- rep(crosstype, ncol(mating))
                  }
              } else {
                  if(length(crosstype) != ncol(mating)){
                      stop("length(crosstype) and ncol(mating) should be same",
                           " if length(crosstype) > 1.",
                           call. = FALSE)
                  }
              }
              if (!all(mating %in% last_gen)) {
                  stop('Invalid member ID(s) specified to mating.',
                       call. = FALSE)
              }
              sapply(seq_along(crosstype), function(i){
                  .validMating(crosstype[i], mating[, i])
              })
              object@crosstype <- c(object@crosstype, list(crosstype))
              object@mating <- c(object@mating, list(mating))
              object@progenies <- c(object@progenies,
                                    list(seq(max(last_gen) + 1,
                                             by = 1,
                                             length.out = ncol(mating))))
              return(object)
          })


#'
#' @rdname assignScheme
setMethod("assignScheme",
          "GbsrScheme",
          function(object, id) {
              valid_id <- unlist(slot(object, "progenies"))
              if(!all(id %in% valid_id)){
                  stop("Invalid specification of member ID(s) to sample(s)!\n",
                       "You can asign only progenies' member ID(s) to sample(s).\n",
                       "Progenies' member ID(s) are: ", paste(valid_id, collapse = ", "),
                       call. = FALSE)
              }
              slot(object, "samples") <- id
              return(object)
          })

#'
#' @param parents_name A vector of strings to indicate names
#' of parental samples. This argument is used internally by showScheme()
#' for the gbsrGenotypeData object.
#' @param pedigree A integer vector indicating the member
#' ID assignment to samples. This argument is used internally by 
#' showScheme() for the gbsrGenotypeData object.
#' @rdname showScheme
#'
setMethod("showScheme",
          "GbsrScheme",
          function(object, parents_name, pedigree) {
              if(length(object@crosstype) == 0){
                  message("No scheme information.\nRun initScheme().")
              }
              n_gen <- length(object@progenies)
              for (i in seq_len(n_gen)) {
                  if (i == 1) {
                      message('\nGeneration: ', appendLF = FALSE)
                      cat("Parents")
                      parents <- object@parents
                      message('\nThe number of parents: ', appendLF = FALSE)
                      cat(length(unique(parents)))
                      message('\nSample IDs: ', appendLF = FALSE)
                      cat(parents_name)
                      message('\nMember IDs: ', appendLF = FALSE)
                      cat(parents)
                      message('\nCross type: ', appendLF = FALSE)
                      cat(object@crosstype[[i]])
                      message('\nMating: ')
                      mating <- object@mating[[i]]
                      rownames(mating) <-
                          paste("mate", seq_len(2), sep = "")
                      colnames(mating) <-
                          paste("combination", seq_len(ncol(mating)), sep = "")
                      print(mating)
                      message('Member IDs of progenites: ', appendLF = FALSE)
                      cat(object@progenies[[i]])
                  } else {
                      message("")
                      message("-------------------------")
                      message('Generation: ', appendLF = FALSE)
                      cat(i - 1)
                      members <- object@progenies[[i - 1]]
                      message('\nThe number of members: ', appendLF = FALSE)
                      cat(length(members))
                      message('\nMember IDs: ', appendLF = FALSE)
                      cat(members)
                      message('\nCross type: ', appendLF = FALSE)
                      cat(object@crosstype[[i]])
                      message('\nMating: ')
                      mating <- object@mating[[i]]
                      rownames(mating) <-
                          paste("mate", seq_len(2), sep = "")
                      colnames(mating) <-
                          paste("combination", seq_len(ncol(mating)), sep = "")
                      print(mating)
                      message('Member IDs of progenites: ', appendLF = FALSE)
                      cat(object@progenies[[i]])
                  }
                  if(i == n_gen){
                      message('\nAssigned member IDs to samples: ',
                              appendLF = FALSE)
                      if(length(object@samples) == 0){
                          cat("Not assigned.")
                      } else {
                          message('')
                          txt <- apply(pedigree, 1, paste, collapse = ":")
                          n_txt <- length(txt)
                          for(i in seq_len(n_txt)){
                              if(i == n_txt){
                                  cat(txt[i])
                              } else {
                                  if(i %% 6 == 0){
                                      cat(paste0(txt[i], "\n"))
                                  } else {
                                      cat(paste0(txt[i], ", "))
                                  }
                              }
                          }
                      }
                  }
              }
          })
