#' GBScleanR: A package to conduct error correction for noisy genotyping by
#'  sequencing (reduced representation sequencing based genotyoing) data.
#'
#' GBScleanR is a package for quality check, filtering, and error correction of
#' genotype data derived from next generation sequcener (NGS) based genotyping
#' platforms. GBScleanR takes Variant Call Format (VCF) file as input. The main
#' function of this package is "clean.geno" which estimates the true genotypes
#' of samples from given read counts for genotype markers using a hidden Markov
#' model with incorporating uneven observation ratio of allelic reads. This
#' implementation gives robust genotype estimation even in noisy genotype data
#' usually observed in Genotyping-By-Sequnencing (GBS) and similar methods, e.g.
#' RADseq. GBScleanR currenly only supports genotype data of
#' biparental populations.
#'
#'
#' @docType package
#' @name GBScleanR
#' @useDynLib GBScleanR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"
