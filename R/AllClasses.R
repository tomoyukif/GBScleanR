#' GbsrScheme class
#'
#' GBScleanR uses breeding scheme information to set the expected
#' number of cross overs in a chromosome which is a required parameter
#' for the genotype error correction with the Hidden Markov model
#' implemented in the "clean" function. This class stores those
#' information including ID of parental samples, type crosses
#' performed at each generation of breeding and population
#' sizes of each generation. This class is not exported.
#'
#'
setClass(
    Class = "GbsrScheme",
    slots = c(
        crosstype = "character",
        pop_size = "numeric",
        mating = "list",
        parents = "numeric",
        progenies = "list"
    )
)

#' GbsrGenotypeData class
#'
#' The `GbsrGenotypeData` class is a container for genotype data from a VCF file
#' together with the metadata associated with the samples and SNPs.
#'
#' @details
#' The `GbsrGenotypeData` class has three slots: data, SNP annotation,
#' and scan annotation. "Samples" or "individuals" subjected to
#' genotyping are called "scan" here with following the way used
#' in the `GWASTools` package. SNP and scan annotation are
#' automatically obtained via the [loadGDS()] function and update
#' when the class methods to calculate SNP wise and scan wise
#'  statistics and set some labels on them.
#'
#' @export
setClass(
    Class = "GbsrGenotypeData",
    contains = "GenotypeData",
    slots = c(scheme = "GbsrScheme")
)
