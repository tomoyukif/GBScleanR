#' Calculate smoothed read ratio per marker.
#'
#' Calculate smoothed read ratio per marker.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param window A numeric to indicate the size of window to calculate
#' smoothed read ratio.
#' @param n_threads An integer value to specify the number of
#' threads used for the calculation. The default is 1 and if `n_threads = NULL`,
#' automatically set half the number of available threads on the computer.
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object.
#'
#' @details
#' Calculate smoothed read ratio per marker.
#'
#' @importFrom RcppParallel defaultNumThreads setThreadOptions
#' @importFrom gdsfmt index.gdsn addfolder.gdsn add.gdsn
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' gds <- smoothReadRatio(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("smoothReadRatio", function(object, node = "raw", window = 1e+05,
                                       n_threads = 1, verbose = TRUE, ...)
    standardGeneric("smoothReadRatio"))

## Calculate smoothed read ratio per marker.
#' @rdname smoothReadRatio
setMethod("smoothReadRatio",
          "GbsrGenotypeData",
          function(object, node, window, n_threads, verbose){
              node <- match.arg(node,  c("raw", "filt"))
              if(is.null(slot(object, "sample")[["parents"]])){
                  parents <- FALSE
              } else {
                  parents <- TRUE
              }
              read <- getRead(object, node = node, parents = parents)
              chr <- getChromosome(object)
              pos <- getPosition(object)

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

              ratio <- NULL
              for(i in unique(getChromosome(object))){
                  target_marker <- chr == i
                  ratio <- cbind(ratio,
                                 calc_ratio(read$ref[, target_marker],
                                            read$alt[, target_marker],
                                            pos[target_marker],
                                            window))
              }
              if(verbose){
                  message("Storing data in the GDS file.")
              }

              arr <- addfolder.gdsn(index.gdsn(object, "annotation/format"),
                                    "ARR", replace = TRUE)
              put.attr.gdsn(arr, "Number", "1")
              put.attr.gdsn(arr, "Type", "Float")
              put.attr.gdsn(arr, "Description",
                            "Smoothed allelic read ratio calculated by GBScleanR")

              arr_data <- add.gdsn(arr, "data", ratio, dim(ratio),
                                   compress = "LZMA_RA", storage = "float64",
                                   replace = TRUE)
              readmode.gdsn(arr_data)
              return(object)
          })


#' Draw line plots of proportion of reference allele read counts
#'  per marker per sample.
#'
#' This function calculate a proportion of reference allele read
#' counts per marker per sample and
#' draw line plots of them in facets for each chromosome for each sample.
#'
#' @param x A [GbsrGenotypeData] object.
#' @param coord A vector with two integer specifying the number of rows and
#' columns to draw faceted line plots for chromosomes.
#' @param chr A vector of indexes to specify chromosomes to be drawn.
#' @param ind An index to specify samples to be drawn.
#' @param node Either one of "raw", "filt", and "cor" to output raw
#' genotype data, filtered genotype data, or corrected genotype data,
#' respectively.
#' @param dot_fill A string to indicate the dot color in a plot.
#'
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#'
#' plotReadRatio(gdata, ind = 1)
#'
#' # Close the connection to the GDS file
#' closeGDS(gdata)
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point labs
#' @importFrom ggplot2 ylim xlab ylab facet_wrap theme
plotSmoothedRatio <- function(x,
                              coord = NULL,
                              chr = NULL,
                              ind = 1,
                              dot_fill = "blue",
                              line_color = "magenta") {
    pos <- ad <- NULL

    if (is.null(chr)) {
        chr <- rep(TRUE, nsam(x))
    } else {
        chr <- getChromosome(x) %in% chr
    }

    id <- getSamID(x)[ind]
    ad <- getGenotype(x, node = "ratio")

    df <- data.frame(chr = getChromosome(x)[chr],
                     pos = getPosition(x)[chr],
                     ad = ad[ind, chr])

    p <- ggplot(df, aes(x = pos * 10^-6, y = ad, group = chr, color = chr)) +
        geom_point(size=0.8, alpha=0.2, stroke=0) +
        geom_line(color=line_color) +
        labs(title=paste0("Alternative allele read ratio: ", id)) +
        ylim(0, 1) +
        xlab("Position (Mb)") +
        ylab("Alternative allele read ratio")

    p <- p + facet_wrap(~ chr, coord[1], coord[2], "free_x",
                        dir="v", strip.position="right") +
        theme(legend.position="none")
    return(p)
}



#'
#' Estimate dosage
#'
#'
#' @param object A [GbsrGenotypeData] object.
#' @param n_threads An integer value to specify the number of
#' threads used for the calculation. The default is 1 and if `n_threads = NULL`,
#' automatically set half the number of available threads on the computer.
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object.
#'
#' @details
#' Calculate smoothed read ratio per marker.
#'
#' @importFrom RcppParallel defaultNumThreads setThreadOptions
#' @importFrom gdsfmt index.gdsn addfolder.gdsn add.gdsn exist.gdsn
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' gds <- estDosage(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("estDosage", function(object, recomb_rate = 0.04, error_rate = 0.0025,
                                 n_threads = 1, smooth = TRUE, mindp = 10,
                                 ...)
    standardGeneric("estDosage"))

#' @rdname estDosage
setMethod("estDosage",
          "GbsrGenotypeData",
          function(object,
                   recomb_rate,
                   error_rate,
                   n_threads,
                   smooth,
                   mindp) {

              if(smooth){
                  if(!exist.gdsn(object, "annotation/format/ARR")){
                      stop("No smoothed allelic read ratio data!",
                           "\n Run smoothReadRatio().")
                  }
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
              chr <- unique(getChromosome(object))
              for (chr_i in chr) {

                  n_mar_i <- nmar(object, FALSE, chr_i)
                  message("\nNow cleaning chr ", chr_i, "...")

                  best_dosage <- .cleanEachChr(object, chr_i, error_rate,
                                               recomb_rate, smooth, mindp)
                  valid_mar_i <- validMar(object, chr_i)
                  .saveDosage(object, best_dosage, n_mar_i, valid_mar_i)
              }
              closeGDS(object, verbose = FALSE)
              seqOptimize(object$filename, "by.sample",
                          "EDS", verbose = FALSE)
              object <- reopenGDS(object)
              return(object)
          })

.initGDS <- function(object, het_parent) {
    cgt <- addfolder.gdsn(index.gdsn(object, "annotation/format"), "EDS",
                          replace = TRUE)
    add.gdsn(cgt, "data", storage = "bit6", compress = "", replace = TRUE)
}

.saveDosage <- function(object, best_hap, n_mar, valid_index) {
    output <- array(0, c(nsam(object, FALSE), n_mar))
    i_sample <- c(which(slot(object, "sample")$parents != 0), which(validSam(object)))
    output[i_sample, valid_index] <- best_hap

    hap_gdsn <- index.gdsn(object, "annotation/format/EDS/data")
    gdsn_dim <- objdesp.gdsn(hap_gdsn)$dim
    if (gdsn_dim[1] == 0) {
        add.gdsn(index.gdsn(object, "annotation/format/EDS"), "data", output,
                 "bit6", compress = "", replace = TRUE)
    } else {
        append.gdsn(hap_gdsn, output)
    }
}

.cleanEachChr <- function(object,
                          chr_i,
                          error_rate,
                          recomb_rate,
                          smooth,
                          mindp) {
    param_list <- .getParams(object, chr_i, error_rate, recomb_rate, smooth, mindp)

    param_list$trans_prob <- matrix(param_list$trans_prob,
                                    dim(param_list$trans_prob)[1])
    best_dosage <- dosage_viterbi(ref = param_list$reads$ref,
                                  alt = param_list$reads$alt,
                                  ratio = param_list$reads$ratio,
                                  eseq_in = param_list$error_rate,
                                  trans_prob = param_list$trans_prob,
                                  init_prob = param_list$init_prob,
                                  n_o = param_list$n_samples,
                                  n_m = param_list$n_mar,
                                  ploidy = param_list$n_ploidy,
                                  mindp = param_list$mindp)
    return(t(best_dosage))
}

.getParams <- function(object, chr_i, error_rate, recomb_rate, smooth, mindp) {
    reads <- .loadReadCounts(object, chr_i, smooth)
    n_samples <- nsam(object)
    n_alleles <- 2
    n_ploidy <- attributes(slot(object, "sample"))$ploidy
    pos <- getPosition(object, TRUE, chr_i)
    n_mar <- nmar(object, TRUE, chr_i)
    trans_prob <- .transitionProb(pos, recomb_rate, n_ploidy)
    init_prob <- .getInitProb(trans_prob[,, 1])

    return(list(n_samples = n_samples,
                n_alleles = n_alleles,
                n_ploidy = n_ploidy,
                n_mar = n_mar,
                error_rate = c(1 - error_rate, error_rate),
                recomb_rate = recomb_rate,
                reads = reads,
                trans_prob = log10(trans_prob),
                init_prob = log10(init_prob),
                mindp = mindp))
}

.loadReadCounts <- function(object, chr_i, smooth) {
    if (exist.gdsn(object, "annotation/format/FAD")) {
        ad_node <- "filt"
    } else {
        ad_node <- "raw"
    }
    reads <- getRead(object, ad_node, FALSE, TRUE, chr_i)

    if(smooth){
        ratio <- getGenotype(object, "ratio", FALSE, TRUE, chr_i)
        ratio[is.na(ratio)] <- -1
        ratio[is.nan(ratio)] <- -1
        ratio[is.infinite(ratio)] <- -1
        reads$alt <- reads$ref <- matrix(-1, nrow = nrow(reads$ref),
                                         ncol = ncol(reads$ref))
    } else {
        dp <- reads$ref + reads$alt
        dp[dp == 0] <- NA
        ratio <- reads$ref / dp
        ratio[is.na(ratio)] <- -1
    }
    return(list(ref = reads$ref,
                alt = reads$alt,
                ratio = ratio))
}

.transitionProb <- function(pos, recomb_rate, n_ploidy){

    # Calculate the recombination frequency between each pair of flanking markers
    mar_dist <- abs(diff(pos) * 1e-06)
    rf <- tanh(2 * mar_dist * recomb_rate) / 2

    # Calculate tansition probabilities from each states to each state
    dosage <- seq(0, n_ploidy, 1)
    dosage_cahnge <- sapply(dosage, function(x){
        return(abs(x - dosage))
    })

    prob <- sapply(dosage_cahnge, function(x){
        return(rf^x * (1 - rf)^(n_ploidy - x))
    })

    prob <- array(t(prob), c(n_ploidy + 1, n_ploidy + 1, length(mar_dist)))
    return(prob)
}

.getInitProb <- function(prob) {
    ev1 <- eigen(t(prob))$vectors[, 1]
    init <- ev1 / sum(ev1)
    return(init)
}
