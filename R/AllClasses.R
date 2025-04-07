#' Class `GbsrScheme`
#'
#' [GBScleanR] uses breeding scheme information to set the expected
#' number of cross overs in a chromosome which is a required parameter
#' for the genotype error correction with the hidden Markov model
#' implemented in the [estGeno()] function. This class stores those
#' information including ID of parental samples, type crosses
#' performed at each generation of breeding and population
#' sizes of each generation. This class is not exported.
#'
#' @slot crosstype A vector of strings indicating the type of
#' crossing done at each generation.
#' @slot mating A list of matrices showing combinations member
#' IDs of samples mated.
#' @slot parents A vector of member IDs of parents.
#' @slot progenies A vector of memeber IDs of progenies produced
#' at each generation.
#' @slot samples A vector of member IDs of samples indicating which samples are
#' derived from which pedigrees.
#'
#' @rdname GbsrScheme-class
#' @aliases GbsrScheme-class GbsrScheme
#' @importFrom methods setClass slot
#' @exportClass GbsrScheme
#' @seealso [GbsrGenotypeData] dnd [loadGDS()].
#'
#' @examples
#' # [loadGDS()] initialize a `GbsrScheme` object internally and
#' # attache it to the shceme slot of a [GbsrGenotypeData] object.
#'
#' # Load data in the GDS file and instantiate
#' # a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Print the information stored in the `GbsrScheme` object.
#' showScheme(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
setClass(Class = "GbsrScheme",
         slots = c(crosstype = "list",
                   mating = "list",
                   parents = "numeric",
                   progenies = "list",
                   samples = "numeric"))


#' Class [GbsrGenotypeData]
#'
#' The [GbsrGenotypeData] class is the main class of [GBScleanR] and
#' user work with this class object.
#'
#' @details
#' The [GbsrGenotypeData] class is an extention of \code{\link[SeqArray]{SeqVarGDSClass-class}} in the
#' [SeqArray] package to store summary data of genotypes and reads and a
#' [GbsrScheme] object that contains mating scheme information of the given
#' population..
#' The slots `marker` and `sample` store a [data.frame] object for variant-wise
#' and sample-wise summary information, respectively. The `scheme` slot holds a
#' [GbsrScheme] object. The function [loadGDS()] initialize the
#' [GbsrGenotypeData] class.
#'
#' @importClassesFrom SeqArray SeqVarGDSClass
#' @aliases  GbsrGenotypeData-class GbsrGenotypeData
#' @slot marker A [data.frame] object.
#' @slot sample A [data.frame] object.
#' @slot scheme A [GbsrScheme] object.
#'
#' @examples
#' # [loadGDS()] initialize the [GbsrGenotypeData] object.
#'
#' # Load a GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Close connection to the GDS file.
#' closeGDS(gds)
#'
#' @exportClass GbsrGenotypeData
#' @importFrom methods setClass slot
#'
setClass(Class = "GbsrGenotypeData",
         contains = "SeqVarGDSClass",
         slots = c(marker = "data.frame",
                   sample = "data.frame",
                   scheme = "GbsrScheme"))
