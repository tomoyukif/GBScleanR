setMethod("initScheme",
          "GbsrScheme",
          function(object, crosstype, mating, parents, ...){
            crosstype <- match.arg(crosstype,
                                   c("pairing", "random"),
                                   FALSE)
            if(length(parents) / 2 != ncol(mating)){
              stop('The number of mating combinations should match with half the number of parents.')
            }
            if(!all(parents %in% mating)){
              stop('All parents should be listed in the matrix specified as the "mating" argument.')
            }
            if(length(parents) != length(mating) & crosstype == "pairing"){
              stop('Invalid mating combinations.')
            }
            object@parents <- parents
            object@crosstype <-  crosstype
            object@pop_size <- length(parents)
            object@mating <- list(mating)
            object@progenies <- list(switch(crosstype,
                                            pairing = seq(max(parents) + 1,
                                                          by = 1,
                                                          length.out = length(parents) / 2),
                                            random = seq(max(parents) + 1,
                                                         by = 1,
                                                         length.out = length(parents))))
            return(object)
          })


setMethod("addScheme",
          "GbsrScheme",
          function(object, crosstype, mating, pop_size, ...){
            crosstype <- match.arg(crosstype,
                                   c("pairing", "selfing", "sibling", "random"),
                                   FALSE)

            last_gen <- unlist(tail(object@progenies, 1))
            last_id <- max(last_gen)
            n_last <- length(last_gen)

            if(crosstype == "random" & is.na(pop_size)){
              stop('"pop_size" is required for crosstype = "random".')
            }
            if(crosstype %in% c("selfing", "sibling")){
              pop_size <- 1
            } else if(crosstype == "pairing"){
              pop_size <- n_last
            }
            if(is.na(mating[1])){
              if(crosstype == "pairing"){
                stop('Need a "mating" matrix for crosstype = "paring".')
              } else {
                mating <- matrix(last_id, 2, 1)
              }
            }

            object@crosstype <- c(object@crosstype, crosstype)
            object@pop_size <- c(object@pop_size, pop_size)
            object@mating <- c(object@mating, list(mating))
            object@progenies <- c(object@progenies,
                                  list(switch(crosstype,
                                              pairing = seq(last_id + 1,
                                                            by = 1,
                                                            length.out = n_last / 2),
                                              random = seq(last_id + 1,
                                                           by = 1,
                                                           length.out = n_last),
                                              selfing = last_id + 1,
                                              sibling = last_id + 1)))
            return(object)
          })


setMethod("showScheme",
          "GbsrScheme",
          function(object, parents_name, ...){
            for(i in 1:length(object@crosstype)){
              if(i == 1){
                message('Generation: Parents')
                parents <- object@parents
                message('The number of parents: ', length(parents))
                message('Scan IDs: ')
                print(parents_name)
                message('Member IDs: ')
                print(parents)
                message('Cross type: ', object@crosstype[i])
                message('Mating: ')
                mating <- object@mating[[i]]
                rownames(mating) <- paste("mate", 1:2, sep = "")
                colnames(mating) <- paste("combination", 1:ncol(mating), sep = "")
                print(mating)
                message('Member IDs of progenites: ')
                print(object@progenies[[i]])
              } else {
                message("")
                message("-------------------------")
                message('Generation: ', i - 1)
                members <- object@progenies[[i - 1]]
                message('The number of members: ', length(members))
                message('Member IDs: ')
                print(members)
                message('Cross type: ', object@crosstype[i])
                message('Mating: ')
                mating <- object@mating[[i]]
                rownames(mating) <- paste("mate", 1:2, sep = "")
                colnames(mating) <- paste("combination", 1:ncol(mating), sep = "")
                print(mating)
                message('Member IDs of progenites: ')
                print(object@progenies[[i]])
              }
            }
          })

