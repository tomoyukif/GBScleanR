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
#' @slot pop_size A vector  of integers of the population size
#'  of each generation.
#' @slot mating A list of matrices showing combinations member 
#' IDs of samples mated.
#' @slot parents A vector of member IDs of parents.
#' @slot progenies A vector of memeber IDs of progenies produced 
#' at each generation.
#' 
#' @rdname GbsrScheme-class
#' @aliases GbsrScheme-class GbsrScheme
#' @import methods
#' @exportClass GbsrScheme
#' @seealso [GbsrGenotypeData] dnd [loadGDS()].
#' 
#' @examples 
#' # `loadGDS()` initialize a `GbsrScheme` object internally and 
#' # attache it to the shceme slot of a [GbsrGenotypeData] object.
#' 
#' # Load data in the GDS file and instantiate 
#' # a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#' 
#' # Print the information stored in the `GbsrScheme` object.
#' showScheme(gdata)
#' 
#' # Close the connection to the GDS file.
#' closeGDS(gdata)
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


#' Class `GbsrGenotypeData`
#' 
#' The `GbsrGenotypeData` class is the main class of [GBScleanR] and 
#' user work with this class object. 
#'
#' @details
#' The `GbsrGenotypeData` class has slots to store summarized genotype 
#' data, a [GbsrScheme] object, and the connection information
#' to a GDS file which can be generated from a VCF file of your 
#' genotype data via [gbsrVCF2GDS()].
#' The [GbsrGenotypeData] class has four slots: `data`, `snpAnnot`,
#' and `scanAnnot` and `scheme.` "Samples" or "individuals" subjected to
#' genotyping are called "scan" here with following the way used
#' in [GWASTools]. The slots `snpAnnot` and `scanAnnot` store
#' a [SnpAnnotationDataFrame] object and a 
#' [ScanAnnotationDataFrame] object, respectively. 
#' The slot `data` stores a [GdsGenotypeReader] object. 
#' Those three classes are included from [GWASTools]. 
#' The slot `scheme` holds a [GbsrScheme] object.
#' The function [loadGDS()] initialize those class objects internally.
#' 
#' @importClassesFrom GWASTools GenotypeData
#' @aliases  GbsrGenotypeData-class GbsrGenotypeData
#' @slot data A [GdsGenotypeReader] object.
#' @slot snpAnnot A [SnpAnnotationDataFrame] object.
#' @slot scanAnnot A [ScanAnnotationDataFrame] object.
#' @slot scheme A [GbsrScheme] object.
#' 
#' @examples
#' # `loadGDS()` initialize objects in the slots of a `GbsrGenotypeData`
#' # object internally.
#' 
#' # Load data in the GDS file and instantiate 
#' # a `GbsrGenotypeData` object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#' 
#' # Close connection to the GDS file.
#' closeGDS(gdata)
#'
#' @exportClass GbsrGenotypeData
#' 
setClass(
    Class = "GbsrGenotypeData",
    contains = "GenotypeData",
    slots = c(scheme = "GbsrScheme")
)
