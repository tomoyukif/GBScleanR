#' Open the connection to the GDS file.
#'
#'
#' @param object A [GbsrGenotypeData] object.
#'
#' @details The GbsrGenotypeData object stores the file path of the GDS file
#' even after closing the connection the file. This function open again the
#' connection to the GDS file at the file path stored in the GbsrGenotypeData
#' object. If the GbsrGenotypeData object witch has an open connection to
#' the GDS file, this function will reopen the connection. The data stored in
#' the SnpAnnotationDataFrame and ScanAnnotationDataFrame will not be changed.
#' Thus, you can open a connection with the GDS file with keeping information
#' of filtering and summary statistics.
#'
#' @return A [GbsrGenotypeData] object.
#' @param ... Unused.
#'
#' @export
#'
#' @examples
#' # Use a GDS file of example data.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#'
#' # Instantiation of [GbsrGenotypeData]
#' gds <- loadGDS(gds_fn)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' gds <- openGDS(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)

setGeneric("openGDS", function(object, ...) standardGeneric("openGDS"))


#' Check if a GDS file has been opened or not.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param ... Unused.
#'
#' @return `TRUE` if the GDS file linked to the input [GbsrGenotypeData] object
#'  has been opened, while `FALSE` if closed.
#'
#' @examples
#' # Use a GDS file of example data.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#'
#' # Instantiation of [GbsrGenotypeData]
#' gds <- loadGDS(gds_fn)
#'
#' # Check connection to the GDS file
#' isOpenGDS(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#'@export
#'
setGeneric("isOpenGDS", function(object, ...)
    standardGeneric("isOpenGDS"))


#' Close the connection to the GDS file
#'
#' Close the connection to the GDS file linked to the given
#'  [GbsrGenotypeData] object.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param verbose if TRUE, show information.
#' @param ... Unused.
#'
#' @return NULL.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#'@export
#'
setGeneric("closeGDS", function(object, verbose = TRUE, ...)
    standardGeneric("closeGDS"))


#' Write out the information stored in the [SnpAnnotationDataFrame] slot
#'
#' All the data stored in the `SnpAnnotatoinDataFrame` slot of
#' the [GbsrGenotypeData]
#' object can be saved in the GDS file linked to the given
#' [GbsrGenotypeData] object. You can load the saved data
#' using [loadSnpAnnot()].
#'
#' @param object A [GbsrGenotypeData] object
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object.
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds <- loadGDS(gds_fn)
#'
#' # Save data stored in the [SnpAnnotationDataFrame] object
#' # linked in the slot of the `GbsrGenotypeData.`
#' gds <- saveSnpAnnot(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("saveSnpAnnot", function(object, ...)
    standardGeneric("saveSnpAnnot"))


#' Write out the information stored in the [ScanAnnotationDataFrame] slot
#'
#' All the data stored in a [ScanAnnotationDataFrame] slot of
#' the [GbsrGenotypeData] object can be saved in the GDS file linked
#' to the given [GbsrGenotypeData] object. You can load the saved data
#'  using [loadSnpAnnot()].
#'
#' @param object A [GbsrGenotypeData] object
#' @param ... Unused.
#'
#' @return NA [GbsrGenotypeData] object.
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds <- loadGDS(gds_fn)
#'
#' # Save data stored in a [ScanAnnotationDataFrame] object
#' # linked in the slot of the `GbsrGenotypeData.`
#' gds <- saveSnpAnnot(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("saveScanAnnot", function(object, ...)
    standardGeneric("saveScanAnnot"))


#' Load the stored [SnpAnnotationDataFrame] information
#'
#' All the data stored in the `SnpAnnotationnDataFrame` slot
#' of the [GbsrGenotypeData]
#' object can be saved in the GDS file linked to the given [GbsrGenotypeData]
#'  object via [saveSnpAnnot()]. You can load the saved data using this
#' function.
#'
#' @param object A [GbsrGenotypeData] object
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object.
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds <- loadGDS(gds_fn)
#'
#' # Save data stored in a [ScanAnnotationDataFrame] object
#' # linked in the slot of the `GbsrGenotypeData.`
#' gds <- saveSnpAnnot(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' # Open the connection to the GDS file again.
#' gds <- loadGDS(gds)
#'
#' # Load the saved [SnpAnnotationDataFrame] object
#' gds <- loadSnpAnnot(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#'
#' @export
#'
setGeneric("loadSnpAnnot", function(object, ...)
    standardGeneric("loadSnpAnnot"))


#' Load the stored [ScanAnnotationDataFrame] information
#'
#' All the data stored in the [ScanAnnotationDataFrame] slot of
#' the [GbsrGenotypeData] object can be saved in the GDS file
#' linked to the given [GbsrGenotypeData] object via [saveScanAnnot()].
#' You can load the saved data using this function.
#'
#' @param object A [GbsrGenotypeData] object
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object.
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds <- loadGDS(gds_fn)
#'
#' # Save data stored in a [ScanAnnotationDataFrame] object
#' # linked in the slot of the `GbsrGenotypeData.`
#' gds <- saveScanAnnot(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' # Open the connection to the GDS file again.
#' gds <- loadGDS(gds)
#'
#' # Load the saved [ScanAnnotationDataFrame] object
#' gds <- loadScanAnnot(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("loadScanAnnot", function(object, ...)
    standardGeneric("loadScanAnnot"))

#' Write a VCF file based on data in a GDS file
#'
#' Write out a VCF file with raw, filtered, or corrected genotype data
#' stored in a GDS file. The output VCF file contains only the GT filed,
#' while other annotations, AD, DP and other information will be omitted.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param out_fn A string to specify the path to an output VCF file.
#' @param node Either one of "raw", "filt", and "cor" to output raw
#' genotype data, filtered genotype data, or corrected genotype data,
#' respectively.
#' @param incl_parents A logical value to specify whether parental
#' samples should be included in an output VCF file or not.
#' @param ... Unused.
#' @export
#'
#' @return The path to the VCF file.
#'
#' @details Create a VCF file at location specified by out_fn.
#' The connection to the GDS file of the input [GbsrGenotypeData] object will be
#' automatically closed for internal file handling in this function. Please use
#' [openGDS()] to open the connection again. If you use [loadGDS()], summary
#' statistics and filtering information will be discarded.
#'
#' @importFrom SeqArray seqSNP2GDS seqGDS2VCF seqOpen
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Create a VCF file with data from the GDS file
#' #  connected to the [GbsrGenotypeData] oobject.
#' out_fn <- tempfile("sample_out", fileext = ".vcf.gz")
#' gbsrGDS2VCF(gds, out_fn)
#'
setGeneric("gbsrGDS2VCF", function(object,
                                   out_fn,
                                   node = "raw",
                                   incl_parents = TRUE,
                                   ...)
    standardGeneric("gbsrGDS2VCF"))


#' Return the number of SNPs.
#'
#' This function returns the number of SNPs recorded in the GDS file
#' connected to the given [GbsrGenotypeData] object.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param chr A index to spefcify chromosome to get information.
#' @param ... Unused.
#'
#' @details
#' If `valid = TRUE`, the number of SNPs which are labeled `TRUE` in
#' the [SnpAnnotationDataFrame] slot will be returned. You need the number
#' of over all SNPs, set `valid = FALSE`. [getValidSnp()] tells you
#' which markers are valid.
#'
#' @return An integer value to indicate the number of SNP markers.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' nsnp(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @seealso [getValidSnp()]
#'
#' @export
#'
setGeneric("nsnp", function(object, valid = TRUE, chr = NULL, ...)
    standardGeneric("nsnp"))


#' Return the number of scans (samples).
#'
#' This function returns the number of samples recorded in the GDS file
#' connected to the given [GbsrGenotypeData] object.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' If `valid = TRUE`, the number of samples which are labeled `TRUE`
#' in the [ScanAnnotationDataFrame] slot will be returned. You need
#' the number of over all samples, set `valid = FALSE`.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return An integer value to indicate the number of samples.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' nscan(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @seealso [getValidScan()]
#'
#'@export
#'
setGeneric("nscan", function(object, valid = TRUE, ...)
    standardGeneric("nscan"))


#' Return a logical vector indicating which are valid SNP markers.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param chr A index to spefcify chromosome to get information.
#' @param ... Unused.
#'
#' @return A logical vector of the same length with
#' the number of total SNP markers
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getValidSnp(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @seealso [setValidSnp()]
#'
#' @export
#'
setGeneric("getValidSnp", function(object, chr = NULL, ...)
    standardGeneric("getValidSnp"))


#' Return a logical vector indicating which are valid scans (samples).
#'
#' @param object A [GbsrGenotypeData] object.
#' @param parents A logical value to indicate to set FALSE or TRUE
#' to parental samples. If you specify `parents = "only"`, this
#' function returns a logical vector indicating TRUE for only parental samples.
#' @param ... Unused.
#'
#' @return A logical vector of the same length
#' with the number of total samples.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getValidScan(gds)
#'
#' # Close the connection the GDS file.
#' closeGDS(gds)
#'
#' @seealso [setValidScan()]
#'
#' @export
#'
setGeneric("getValidScan", function(object, parents = FALSE, ...)
    standardGeneric("getValidScan"))


#' Manually set valid SNP markers.
#'
#' If you need manually set valid and invalid SNP markers,
#' you can do it via this function,
#' e.g in the case you conducted a filtering on SNP
#' markers manually by your self.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param new A logical vector of the same length with the over
#' all number of the SNP markers.
#' @param update A logical vector of the same length with the
#' currently valid SNP markers.
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object with updated SNP validity information.
#'
#' @details
#' To over write the current validity information,
#' give a logical vector to `new`.
#' On the other hand, a logical vector specified to `update`
#' will be used to update validity information of the currently
#' valid SNP markers. If you gave a vector for both argument,
#' only the vector passed to `new` will be used to
#' over write the validity information.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Randomly remove SNP markers.
#' gds <- setValidSnp(gds,
#'                      update = sample(c(TRUE, FALSE),
#'                                      nsnp(gds),
#'                                      replace = TRUE))
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @seealso [setSnpFilter()] to filter out SNP markers based on
#' some summary statistics.
#'
#' @export
#'
setGeneric("setValidSnp", function(object, new, update, ...)
    standardGeneric("setValidSnp"))


#' Manually set valid scans (samples).
#'
#' If you need manually set valid and invalid samples,
#' you can do it via this function,
#' e.g in the case you conducted a filtering on
#' samples manually by your self.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param new A logical vector of the same length with
#' the over all number of the samples.
#' @param update A logical vector of the same length with
#' the currently valid samples.
#' @param ... Unused.
#'
#' @details
#' To over write the current validity information,
#' give a logical vector to `new`. On the other hand, a logical
#' vector specified to `update` will be used to update validity
#' information of the currently valid samples. If you gave
#' a vector for both argument, only the vector passed to `new`
#' will be used to over write the validity information.
#'
#' @return A [GbsrGenotypeData] object with updated
#' scan(sample) validity information.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Randomly remove samples.
#' gds <- setValidScan(gds,
#'                       update = sample(c(TRUE, FALSE),
#'                                       nscan(gds),
#'                                       replace = TRUE))
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @seealso [setScanFilter()] to filter out samples based on
#' some summary statistics
#'
#' @export
#'
setGeneric("setValidScan", function(object, new, update, ...)
    standardGeneric("setValidScan"))


#' Get a logical vector indicating flipped SNP markers.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value to specify flip alleles
#' only of valid markers. see also[setSnpFilter()].
#' @param chr A index to spefcify chromosome to get information.
#' @param ... Unused.
#'
#' @details
#' Flipped markers are markers where the alleles expected as
#' reference allele are called as alternative allele. If you
#' specify two parents in the `parents` argument of [setParents()]
#' with `flip = TRUE`, `bi = TRUE`, and `homo = TRUE`, the alleles found
#' in the parent specified as the first element to the `parents`
#' argument are supposed as reference alleles of the markers.
#' If the "expected" reference alleles are not actually
#' called as reference alleles but alternative alleles
#' in the given data. [setParents()] will automatically labels
#' those markers "flipped". The SnpAnnotatoinDataFrame slot sores this
#' information and accessible via [getFlipped()] which gives you
#' a logical vector indicating which markers are labeled as flipped
#' `TRUE` or not flipped `FALSE`. [hasFlipped()] just tells you
#' whether the SnpAnnotatoinDataFrame slot has the information of
#' flipped markers or not.
#'
#' @return A logical vector indicating which markers of alleles were flipped.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Find the IDs of parental samples.
#' parents <- grep("Founder", getScanID(gds), value = TRUE)
#'
#' # Set the parents and flip allele information
#' # if the reference sample (Founder1 in our case) has homozygous
#' # alternative genotype at some markers of which alleles will
#' # be swapped to make the reference sample have homozygous
#' # reference genotype.
#' gds <- setParents(gds, parents = parents, flip = TRUE)
#'
#' getFlipped(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @seealso [setParents()] and [hasFlipped()].
#'
#' @export
#'
setGeneric("getFlipped", function(object, valid = TRUE, chr = NULL, ...)
    standardGeneric("getFlipped"))


#' Get a logical value indicating flipped SNP markers
#' whether information exists.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param ... Unused.
#'
#' @details
#' Flipped markers are markers where the alleles expected as reference
#' allele are called as alternative allele. If you specify two parents
#' in the `parents` argument of [setParents()] with `flip = TRUE`,
#' `bi = TRUE`, and `homo = TRUE`, the alleles found in the parent
#' specified as the first element to the `parents` argument are supposed as
#' reference alleles of the markers. If the "expected" reference alleles
#' are not actually called as reference alleles but alternative alleles
#' in the given data. [setParents()] will automatically labels those
#' markers "flipped". The SnpAnnotatoinDataFrame slot sores this information
#' and accessible via [getFlipped()] which gives you a logical vector
#' indicating which markers are labeled as flipped `TRUE` or
#' not flipped `FALSE`. [hasFlipped()] just tells you whether
#' the SnpAnnotatoinDataFrame slot has the information of flipped
#' markers or not.
#'
#' @return A logical value to indicate the GbsrGenotypeData
#' contains markers at which alleles were flipped.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Find the IDs of parental samples.
#' parents <- grep("Founder", getScanID(gds), value = TRUE)
#'
#' # Set the parents and flip allele information
#' # if the reference sample (Founder1 in our case) has homozygous
#' # alternative genotype at some markers of which alleles will
#' # be swapped to make the reference sample have homozygous
#' # reference genotype.
#' gds <- setParents(gds, parents = parents, flip = TRUE)
#'
#' hasFlipped(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds) # Close the connection to the GDS file
#'
#'
#' @seealso [setParents()] and [getFlipped()].
#'
#' @export
#'
setGeneric("hasFlipped", function(object, ...)
    standardGeneric("hasFlipped"))


#' Get genotype call data.
#'
#' Genotype calls are retrieved from the GDS file linked to the given
#' [GbsrGenotypeData] object.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param chr A integer vector of indexes indicating chromosomes
#' to get read count data.
#' @param node Either of "raw", "filt", and "cor. See details.
#' @param parents A logical value or "only" to include data for
#' parents or to get data only for parents.
#' @param ... Unused.
#'
#' @return An integer matirix of genotype data which is represented
#' by the number of reference alleles at each marker of each sample.
#'
#' @details
#' Genotype call data can be obtained from the "genotype" node,
#' the "filt.genotype" node, or the "corrected.genotype" node
#' of the GDS file with `node = "raw"`, `node = "filt"`, or `node = "raw"`,
#' respectively. If `node = "parents`, the data in the "parents.genotype"
#' node will be returned. The "parents.genotype" node stores phased
#' parental genotypes estimated by the [estGeno()] function.
#' The [setCallFilter()] function generate filtered genotype call data in the
#' "filt.genotype" node which can be accessed as mentioned above.
#' On the other hand, the "corrected.genotype" node can be generated
#' via the [estGeno()] function.
#' If `valid = TRUE`, read counts for only valid marker and valid samples will
#' be obtained.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' geno <- getGenotype(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#'
#' @seealso [setCallFilter()] and [estGeno()]
#'
#' @export
#'
setGeneric("getGenotype", function(object,
                                   valid = TRUE,
                                   chr = NULL,
                                   node = "raw",
                                   parents = FALSE,
                                   ...)
    standardGeneric("getGenotype"))


#' Get haplotype call data.
#'
#' Haplotype calls are retrieved from the GDS file linked to the given
#' [GbsrGenotypeData] object.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param chr A integer vector of indexes indicating chromosomes
#' to get read count data.
#' @param parents A logical value or "only" to include data for
#' parents or to get data only for parents.
#' @param ... Unused.
#'
#' @details
#' Haplotype call data can be obtained from the "estimated.haplotype" node of
#' the GDS file which can be generated via the [estGeno()] function.
#' Thus, this function is valid only after having executed [estGeno()].
#' If `valid = TRUE`, read counts for only valid marker and valid samples will
#' be obtained.
#'
#' @return An integer array of haplotype data. The array have 2 x M x N
#' dimensions, where M is the number of markers and N is the number of
#' samples. Each integer values represent the origin of the haplotype.
#' For example, in the population with two inbred founders, values take
#' either 1 or 2 indicating the hapotype descent from founder 1 and 2.
#' If two outbred founders, values take 1, 2, 3, or 4 indicating
#' the first and second haplotype in founder 1 and the first and
#' second haplotype in founder 2.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Find the IDs of parental samples.
#' parents <- grep("Founder", getScanID(gds), value = TRUE)
#'
#' # Set the parents and flip allele information
#' # if the reference sample (Founder1 in our case) has homozygous
#' # alternative genotype at some markers of which alleles will
#' # be swapped to make the reference sample have homozygous
#' # reference genotype.
#' gds <- setParents(gds, parents = parents, flip = TRUE)
#'
#' # Initialize a scheme object stored in the slot of the GbsrGenotypeData.
#' # We chose `crosstype = "pair"` because two inbred founders were mated
#' # in our breeding scheme.
#' # We also need to specify the mating matrix which has two rows and
#' # one column with integers 1 and 2 indicating a sample (founder)
#' # with the memberID 1 and a sample (founder) with the memberID 2
#' # were mated.
#' gds <- initScheme(gds, crosstype = "pair", mating = cbind(c(1:2)))
#'
#' # Add information of the next cross conducted in our scheme.
#' # We chose 'crosstype = "selfing"', which do not require a
#' # mating matrix.
#' gds <- addScheme(gds, crosstype = "selfing")
#'
#' # Execute error correction by estimating genotype and haplotype of
#' # founders and offspring.
#' gds <- estGeno(gds)
#'
#' hap <- getHaplotype(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @seealso [estGeno()]
#'
#' @export
#'
setGeneric("getHaplotype", function(object,
                                    valid = TRUE,
                                    chr = NULL,
                                    parents = FALSE,
                                    ...)
    standardGeneric("getHaplotype"))


#' Get read count data.
#'
#' Read counts for reference allele and alternative allele are retrieved
#' from the GDS file linked to the given [GbsrGenotypeData] object.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param chr A integer vector of indexes indicating
#' chromosomes to get read count data.
#' @param node Either of "raw" and "filt". See details.
#' @param parents A logical value or "only" to include data
#' for parents or to get data only for parents.
#' @param ... Unused.
#'
#' @return A named list with two elements "ref" and "alt" storing
#' a matrix of reference allele read counts and a matrix of
#' alternative read counts for all markers in all samples.
#'
#' @details
#' Read count data can be obtained from the "annotation/format/AD/data" node,
#' the "annotation/format/AD/norm" node or the
#' "annotation/format/AD/filt.data" node of the GDS file with
#' `node = "raw"`,  `node = "norm"`, or `node = "filt"`, respectively.
#' The [calcReadStats()] function generate normalized read count data
#' in the "annotation/format/AD/norm" node, while [setCallFilter()]
#' function generate filtered read count data
#' in the "annotation/format/AD/filt.data" node which can be accessed as
#' mentioned above.
#' If `valid = TRUE`, read counts for only valid marker and valid samples will
#' be obtained.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' read <- getRead(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @seealso [setCallFilter()]
#'
#' @export
#'
setGeneric("getRead", function(object,
                               valid = TRUE,
                               chr = NULL,
                               node = "raw",
                               parents = FALSE,
                               ...)
    standardGeneric("getRead"))


#' Obtain chromosome information of each SNP marker
#'
#' This function returns indexes or names of chromsomes of each SNP or just a
#' set of unique chromosome names.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param levels A logical value. See details.
#' @param name A logical value. See details.
#' @param ... Unused.
#'
#' @return A vector of integers or strings indicating chromosome IDs.
#'
#' @details
#' A GDS file created via [GBScleanR] stores chromosome names as
#' sequential integers from 1 to N, where N is the number of chromosomes.
#' This function returns those indexes as default. If you need actual names
#' of the chromosomes, set `name = TRUE`.
#' `levels = TRUE` gives you only unique chromosome names with length N.
#' If `valid = TRUE`, the chromosome information of markers which
#' are labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getChromosome(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
#'
setGeneric("getChromosome", function(object,
                                     valid = TRUE,
                                     levels = FALSE,
                                     name = FALSE,
                                     ...)
    standardGeneric("getChromosome"))


#' Obtain physical position information of each SNP marker
#'
#' This function returns physical positions of SNP markers.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param chr A index to spefcify chromosome to get information.
#' @param ... Unused.
#'
#' @return A integer vetor indicating the physical positions of markers.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which
#' are labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getPosition(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getPosition", function(object, valid = TRUE, chr = NULL, ...)
    standardGeneric("getPosition"))


#' Obtain reference allele information of each SNP marker
#'
#' This function returns reference alleles, either
#' of A, T, G, and C, of SNP markers.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param chr A index to spefcify chromosome to get information.
#' @param ... Unused.
#' @return A character vector indicating the reference alleles.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers
#' which are labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be
#' returned. [getValidSnp()] tells you which samples are valid.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getAlleleA(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getAlleleA", function(object, valid = TRUE, chr = NULL, ...)
    standardGeneric("getAlleleA"))


#' Obtain alternative allele information of each SNP marker
#'
#' This function returns alternative alleles,
#' either of A, T, G, and C, of SNP markers.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param chr A index to spefcify chromosome to get information.
#' @param ... Unused.
#'
#' @return A character vector indicating the alternative alleles.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getAlleleB(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getAlleleB", function(object, valid = TRUE, chr = NULL, ...)
    standardGeneric("getAlleleB"))


#' Obtain SNP ID
#'
#' This function returns SNP ID of SNP markers.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param chr A index to spefcify chromosome to get information.
#' @param ... Unused.
#'
#' @return A character vector of SNP IDs.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getSnpID(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getSnpID", function(object, valid = TRUE, chr = NULL, ...)
    standardGeneric("getSnpID"))


#' Obtain scan (sample) ID
#'
#' This function returns scan (sample) ID.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @return A character vector of scan(sample) IDs.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getScanID(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getScanID", function(object, valid = TRUE, ...)
    standardGeneric("getScanID"))


#' Obtain ploidy information of each SNP marker
#'
#' This function returns ploidy of each SNP marker. The ploidy of all
#' the markers in a dataset is a same value and the current implementation
#' of [GBScleanR] only works with data having ploidy = 2 for all markers.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param valid A logical value. See details.
#' @param chr A index to spefcify chromosome to get information.
#' @param ... Unused.
#'
#' @return A integer vector indicating the ploidy of each marker.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getPloidy(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getPloidy", function(object, valid = TRUE, chr = NULL, ...)
    standardGeneric("getPloidy"))


#' Obtain information stored in the "annotation/info" node
#'
#' The "annotation/info" node stores annotation infromation of markers obtained
#' via SNP calling tools like bcftools and GATK.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param var A string to indicate which annotation info should be retrieved.
#' @param valid A logical value. See details.
#' @param chr A index to spefcify chromosome to get information.
#' @param ... Unused.
#'
#' @return A numeric vector of data stored in INFO node of the GDS file.
#'
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Get mapping qualities (MQ) of markers.
#' mq <- getInfo(gds, "MQ")
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#' @export
#'
setGeneric("getInfo", function(object, var, valid = TRUE, chr = NULL, ...)
    standardGeneric("getInfo"))


#' Obtain total reference read counts per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param prop A logical value whether to return values as proportions of
#' total reference read counts in total read counts per SNP or not.
#' @param ... Unused.
#'
#' @return A numeric vector of (proportion of) reference read
#' counts per marker.
#'
#' @details
#' You need to execute [countRead()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getCountReadRef(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getCountReadRef", function(object,
                                       target = "snp",
                                       valid = TRUE,
                                       prop = FALSE,
                                       ...)
    standardGeneric("getCountReadRef"))


#' Obtain total alternative read counts per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param prop A logical value whether to return values as proportions
#' of total alternative read counts in total read counts per SNP or not.
#' @param ... Unused.
#'
#' @return A numeric vector of (proportion of) alternative
#' allele read counts per marker.
#'
#' @details
#' You need to execute [countRead()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getCountReadAlt(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getCountReadAlt", function(object,
                                       target = "snp",
                                       valid = TRUE,
                                       prop = FALSE,
                                       ...)
    standardGeneric("getCountReadAlt"))


#' Obtain total read counts per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [countRead()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A integer vector of total read counts
#' (reference allele reads + alternative allele reads) per marker.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' getCountRead(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getCountRead", function(object,
                                    target = "snp",
                                    valid = TRUE,
                                    ...)
    standardGeneric("getCountRead"))


#' Obtain total reference genotype counts per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as
#' proportions of total reference genotype counts to total
#' non missing genotype counts or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of (proportion of) homozygous reference
#' genotype calls per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the genotype count information and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- countGenotype(gds)
#'
#' getCountGenoRef(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getCountGenoRef", function(object,
                                       target = "snp",
                                       valid = TRUE,
                                       prop = FALSE,
                                       ...)
    standardGeneric("getCountGenoRef"))


#' Obtain total heterozygote counts per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of
#' total heterozygote counts to total non missing genotype counts or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of (proportion of) heterozugous
#' genotype calls per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the genotype count information and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- countGenotype(gds)
#'
#' getCountGenoHet(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getCountGenoHet", function(object,
                                       target = "snp",
                                       valid = TRUE,
                                       prop = FALSE,
                                       ...)
    standardGeneric("getCountGenoHet"))


#' Obtain total alternative genotype counts per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as
#' proportions of total alternative genotype counts to total
#' non missing genotype counts or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of (proportion of) homozygous alternative
#' genotype calls per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the genotype count information and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- countGenotype(gds)
#'
#' getCountGenoAlt(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getCountGenoAlt", function(object,
                                       target = "snp",
                                       valid = TRUE,
                                       prop = FALSE,
                                       ...)
    standardGeneric("getCountGenoAlt"))


#' Obtain total missing genotype counts per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions
#' of total missing genotype counts to the total genotype calls or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of (proportion of)
#' missing genotype calls per marker.
#' @examples
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#' gds <- countGenotype(gds)
#' getCountGenoMissing(gds)
#' closeGDS(gds) # Close the connection to the GDS file
#'
#' @export
#'
setGeneric("getCountGenoMissing", function(object,
                                           target = "snp",
                                           valid = TRUE,
                                           prop = FALSE,
                                           ...)
    standardGeneric("getCountGenoMissing"))


#' Obtain total reference allele counts per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions
#' of total reference allele counts to total non missing allele counts or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of (proportion of) reference alleles per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the genotype count information and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- countGenotype(gds)
#'
#' getCountAlleleRef(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getCountAlleleRef", function(object,
                                         target = "snp",
                                         valid = TRUE,
                                         prop = FALSE,
                                         ...)
    standardGeneric("getCountAlleleRef"))


#' Obtain total alternative allele counts per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions
#' of total alternative allele counts to total non missing allele
#' counts or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of (proportion of) alternative alleles per marker.
#' @examples
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#' gds <- countGenotype(gds)
#' getCountAlleleAlt(gds)
#' closeGDS(gds) # Close the connection to the GDS file
#'
#' @export
#'
setGeneric("getCountAlleleAlt", function(object,
                                         target = "snp",
                                         valid = TRUE,
                                         prop = FALSE,
                                         ...)
    standardGeneric("getCountAlleleAlt"))


#' Obtain total missing allele counts per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of
#' total missing allele counts to the total allele number or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of (proportion of) missing alleles per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the genotype count information and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- countGenotype(gds)
#'
#' getCountAlleleMissing(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getCountAlleleMissing", function(object,
                                             target = "snp",
                                             valid = TRUE,
                                             prop = FALSE,
                                             ...)
    standardGeneric("getCountAlleleMissing"))


#' Obtain mean values of total reference read counts per
#' SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of the mean values of
#' reference allele reads per marker.
#' @examples
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#' gds <- calcReadStats(gds)
#' getMeanReadRef(gds)
#' closeGDS(gds) # Close the connection to the GDS file
#'
#' @export
#'
setGeneric("getMeanReadRef", function(object,
                                      target = "snp",
                                      valid = TRUE,
                                      ...)
    standardGeneric("getMeanReadRef"))


#' Obtain mean values of total alternative read counts
#' per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of the mean values of
#' alternative allele reads per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Calculate means, standard deviations, quantiles of read counts
#' # per marker and per sample with or without standardization of
#' # the counts and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- calcReadStats(gds)
#'
#' getMeanReadAlt(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getMeanReadAlt", function(object,
                                      target = "snp",
                                      valid = TRUE,
                                      ...)
    standardGeneric("getMeanReadAlt"))


#' Obtain standard deviations of total reference read counts
#' per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate summary statistics to be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#'
#' @return A numeric vector of the standard deviations of
#' reference allele reads per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Calculate means, standard deviations, quantiles of read counts
#' # per marker and per sample with or without standardization of
#' # the counts and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- calcReadStats(gds)
#'
#' getSDReadRef(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getSDReadRef", function(object,
                                    target = "snp",
                                    valid = TRUE,
                                    ...)
    standardGeneric("getSDReadRef"))


#' Obtain standard deviations of total alternative read counts per
#' SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of the standard deviations of
#' alternative allele reads per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Calculate means, standard deviations, quantiles of read counts
#' # per marker and per sample with or without standardization of
#' # the counts and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- calcReadStats(gds)
#'
#' getSDReadAlt(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getSDReadAlt", function(object,
                                    target = "snp",
                                    valid = TRUE,
                                    ...)
    standardGeneric("getSDReadAlt"))


#' Obtain quantile values of total reference read counts per
#' SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param q A numeric value \[0-1\] to indicate quantile to obtain.
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of the quantile values of
#' alternative allele reads per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Calculate means, standard deviations, quantiles of read counts
#' # per marker and per sample with or without standardization of
#' # the counts and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- calcReadStats(gds)
#'
#' getQtileReadRef(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getQtileReadRef", function(object,
                                       target = "snp",
                                       q = 0.5,
                                       valid = TRUE,
                                       ...)
    standardGeneric("getQtileReadRef"))


#' Obtain quantile values of total alternative read counts per
#' SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param q A numeric value \[0-1\] to indicate quantile to obtain.
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of the quantile values of
#' alternative allele reads per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Calculate means, standard deviations, quantiles of read counts
#' # per marker and per sample with or without standardization of
#' # the counts and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- calcReadStats(gds)
#'
#' getQtileReadAlt(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getQtileReadAlt", function(object,
                                       target = "snp",
                                       q = 0.5,
                                       valid = TRUE,
                                       ...)
    standardGeneric("getQtileReadAlt"))


#' Obtain minor allele frequencies per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of the minor allele frequencies per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the genotype count information and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- countGenotype(gds)
#'
#' getMAF(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getMAF", function(object,
                              target = "snp",
                              valid = TRUE,
                              ...)
    standardGeneric("getMAF"))


#' Obtain minor allele counts per SNP or per scan (sample)
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#'
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are
#' labeled `TRUE` in the [ScanAnnotationDataFrame] slot will be returned.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @return A numeric vector of the minor allele counts per marker.
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the genotype count information and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- countGenotype(gds)
#'
#' getMAC(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("getMAC", function(object,
                              target = "snp",
                              valid = TRUE,
                              ...)
    standardGeneric("getMAC"))


#' Count genotype calls and alleles per sample and per marker.
#'
#' This function calculates several summary statistics of
#' genotype calls and alleles per marker and per sample.
#' Those values will be stored in the SnpAnnotaionDataFrame slot
#' and the [ScanAnnotationDataFrame] slot and obtained via getter
#' functions, e.g.s
#' [getCountGenoRef()], [getCountAlleleRef()], and [getMAF()].
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param node Either of "raw", "filt", and "cor". See details.
#' @param ... Unused.
#'
#' @details
#' #' Genotype call data can be obtained from the "genotype" node,
#' the "filt.genotype" node, or the "corrected.genotype" node of
#' the GDS file with `node = "raw"`, `node = "filt"`, or `node = "raw"`,
#' respectively.
#' The [setCallFilter()] function generate filtered genotype call data in the
#' "filt.genotype" node which can be accessed as mentioned above.
#' On the other hand, the "corrected.genotype" node can be generated
#' via the [estGeno()] function.
#'
#' @return A [GbsrGenotypeData] object with genotype count information.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the genotype count information and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- countGenotype(gds)
#'
#' # Get the proportion of missing genotype per sample.
#' sample_missing_rate <- getCountGenoMissing(gds,
#'                                            target = "scan",
#'                                            prop = TRUE)
#'
#' # Get the minor allele frequency per marker.
#' marker_minor_allele_freq <- getMAF(gds, target = "snp")
#'
#' # Draw histograms of the missing rate per sample and marker.
#' histGBSR(gds, stats = "missing")
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("countGenotype", function(object,
                                     target = "both",
                                     node = "raw",
                                     ...)
    standardGeneric("countGenotype"))


#' Count reads per sample and per marker.
#'
#' This function calculates several summary statistics of read counts
#' per marker and per sample. Those values will be stored
#' in the SnpAnnotaionDataFrame slot and the [ScanAnnotationDataFrame] slot
#' and obtained via getter functions, e.g.
#' [getCountReadRef()] and [getCountReadAlt()].
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param node Either of "raw" and "filt". See details.
#' @param ... Unused.
#'
#' @details
#' Read count data can be obtained from the "annotation/format/AD/data" node
#' or the "annotation/format/AD/filt.data" node of the GDS file
#' with `node = "raw"` or `node = "filt"`, respectively.
#' The [setCallFilter()] function generate filtered read count data
#' in the "annotation/format/AD/filt.data" node which can be accessed as
#' mentioned above.
#'
#' @return A [GbsrGenotypeData] object with read count information.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the read count information and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- countRead(gds)
#'
#' # Get the total read counts per marker
#' read_depth_per_marker <- getCountRead(gds, target = "snp")
#'
#' # Get the proportion of reference allele rads per marker.
#' reference_read_freq <- getCountReadRef(gds, target = "snp", prop = TRUE)
#'
#' # Draw histgrams of reference allele read counts per sample and marker.
#' histGBSR(gds, stats = "ad_ref")
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("countRead", function(object,
                                 target = "both",
                                 node = "raw",
                                 ...)
    standardGeneric("countRead"))


#' Calculate mean, standard deviation, and quantile values of
#' normalized read counts per sample and per marker.
#'
#' This function first calculates normalized allele read counts by dividing
#' allele read counts at each marker in each sample by the total allele read
#' of the sample followed by multiplication by 10^6. In other words, it
#' calculates reads per million (rpm). Then, the function calculates
#' mean, standard deviation, quantile values of rpm per marker and per sample.
#' The results will be stored in the SnpAnnotaionDataFrame slot and the
#' [ScanAnnotationDataFrame] slot and obtained via getter functions, e.g.
#' [getMeanReadRef()] and [getQtileReadAlt()].
#'
#' @param object A [GbsrGenotypeData] object.
#' @param target Either of "snp" and "scan".
#' @param q A numeric value \[0-1\] to indicate quantile to obtain.
#' @param ... Unused.
#'
#' @details
#' Read count data can be obtained from the "annotation/format/AD/data" node
#' or the "annotation/format/AD/filt.data" node of the GDS file
#' with `node = "raw"` or `node = "filt"`, respectively.
#' The [setCallFilter()] function generate filtered read count data
#' in the "annotation/format/AD/filt.data" node which can be accessed as
#' mentioned above. The calculation of mean, standard deviation, and quantile
#' values omits 0 in the read count data.
#'
#' @importFrom stats sd quantile
#'
#' @return A [GbsrGenotypeData] object with read statistics information.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Calculate means, standard deviations, quantiles of read counts
#' # per marker and per sample with or without standardization of
#' # the counts and store them in the
#' # [SnpAnnotationDataFrame] and [ScanAnnotationDataFrame] objects
#' # linked at the slots of the [GbsrGenotypeData] object.
#' gds <- calcReadStats(gds, q = 0.5)
#'
#' # Get the means of reference allele read counts (rpm) per marker.
#' mean_reference_read_depth <- getMeanReadRef(gds, target = "snp")
#'
#' # Get the 0.5 percentiles (medians) of the alternative allele
#' # read counts (rpm) per marker.
#' median_reference_read_depth <- getQtileReadAlt(gds,
#'                                                target = "snp",
#'                                                q = 0.5)
#' # Draw histograms of the means of reference allele read counts (rpm)
#' # per sample and marker.
#' histGBSR(gds, stats = "mean_ref")
#'
#' # Draw histograms of the 0.5 percentiles (medians) of
#' # alternative allele read counts (rpm) per sample and marker.
#' histGBSR(gds, stats = "qtile_alt", q = 0.5)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("calcReadStats", function(object,
                                     target = "both",
                                     q = NULL,
                                     ...)
    standardGeneric("calcReadStats"))


#' Set labels to samples which should be recognized as
#' parents of the population to be subjected to error correction.
#'
#' Specify two or more samples in the dataset as parents
#' of the population. Markers will be filtered out up on your specification.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param parents A vector of strings with at least length two.
#' The specified strings should match with the samples
#' ID available via [getScanID()].
#' @param flip A logical value to indicate whether markers should
#' be checked for "flip". See details.
#' @param mono A logical value whether to filter out markers which
#' are not monomorphic in parents.
#' @param bi A logical value whether to filter out marekrs which
#' are not biallelic between parents.
#' @param ... Unused.
#'
#' @details
#' The `clean` function of [GBScleanR] uses read count information of
#' samples and their parents separately to estimate most probable
#' genotype calls of them. Therefore, you must specify proper samples
#' as parents via this function. If you would like to remove SNP markers
#' which are not biallelic and/or not monomorphic in each parent,
#' set `mono = TRUE` and `bi = TRUE`. `flip = TRUE` flips alleles of
#' markers where the alleles expected as reference allele are called
#' as alternative allele. The alleles found in the parent specified as
#' the first element to the `parents` argument are supposed as reference
#' alleles of the markers. If the "expected" reference alleles are not
#' actually called as reference alleles but alternative alleles
#' in the given data. setParents() will automatically labels
#' those markers "flipped".
#' The SnpAnnotatoinDataFrame slot sores this information and accessible
#' via [getFlipped()] which gives you a logical vector
#' indicating which markers are labeled as flipped `TRUE` or
#' not flipped `FALSE`. [hasFlipped()] just tells you whether
#' the SnpAnnotatoinDataFrame slot has the information of flipped
#' markers or not.
#'
#' @return A [GbsrGenotypeData] object with parents information.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Find the IDs of parental samples.
#' parents <- grep("Founder", getScanID(gds), value = TRUE)
#'
#' # Set the parents and flip allele information
#' # if the reference sample (Founder1 in our case) has homozygous
#' # alternative genotype at some markers of which alleles will
#' # be swapped to make the reference sample have homozygous
#' # reference genotype.
#' gds <- setParents(gds, parents = parents, flip = TRUE)
#'
#' # Initialize a scheme object stored in the slot of the GbsrGenotypeData.
#' # We chose `crosstype = "pair"` because two inbred founders were mated
#' # in our breeding scheme.
#' # We also need to specify the mating matrix which has two rows and
#' # one column with integers 1 and 2 indicating a sample (founder)
#' # with the memberID 1 and a sample (founder) with the memberID 2
#' # were mated.
#' gds <- initScheme(gds, crosstype = "pair", mating = cbind(c(1:2)))
#'
#' # Add information of the next cross conducted in our scheme.
#' # We chose 'crosstype = "selfing"', which do not require a
#' # mating matrix.
#' gds <- addScheme(gds, crosstype = "selfing")
#'
#' # Execute error correction by estimating genotype and haplotype of
#' # founders and offspring.
#' gds <- estGeno(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("setParents", function(object,
                                  parents,
                                  flip = FALSE,
                                  mono = FALSE,
                                  bi = FALSE,
                                  ...)
    standardGeneric("setParents"))


#' Get parental sample information
#'
#' This function returns scan IDs, member IDs and indexes of parental samples
#' set via [setParents()]. Scan IDs are IDs given by user or obtained from the
#' original VCF file. Member IDs are serial numbers assigned by [setParents()].
#'
#' @param object A [GbsrGenotypeData] object.
#' @param bool If TRUE, the function returns a logical vector indicating
#' which scans (samples) have been set as parents.
#' @param ... Unused.
#'
#' @export
#'
#' @return A data frame of parents information indicating scanIDs, memberIDs
#' and indexes of parental lines assigned via [setParents()].
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Find the IDs of parental samples.
#' parents <- grep("Founder", getScanID(gds), value = TRUE)
#'
#' # Set the parents.
#' gds <- setParents(gds, parents = parents, flip = TRUE)
#'
#' # Get the information of parents.
#' getParents(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#'
setGeneric("getParents", function(object,
                                  bool=FALSE,
                                  ...)
    standardGeneric("getParents"))


#' Filter out each genotype call meeting criteria
#'
#' Perform filtering of each genotype call, neither markers
#' nor samples. Each genotype call is supported by its read counts
#' for the reference allele and the alternative allele of
#' a marker of a sample. `setCallFilter()` set missing to
#' the genotype calls which are not reliable enough and set zero
#' to reference and alternative read counts of the genotype calls.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param dp_count A numeric vector with length two specifying lower and
#' upper limit of total read counts (reference reads + alternative reads).
#' @param ref_count A numeric vector with length two specifying lower and
#' upper limit of reference read counts.
#' @param alt_count A numeric vector with length two specifying lower and
#' upper limit of alternative read counts.
#' @param norm_dp_count A numeric vector with length two specifying lower
#' and upper limit of normalized total read counts
#' (normalized reference reads + normalized alternaitve reads).
#' @param norm_ref_count A numeric vector with length two specifying lower
#' and upper limit of normalized reference read counts.
#' @param norm_alt_count A numeric vector with length two specifying lower
#' and upper limit of normalized alternative read counts
#' @param scan_dp_qtile A numeric vector with length two specifying lower
#' and upper limit of quantile of total read counts in each scan (sample).
#'@param scan_ref_qtile A numeric vector with length two specifying lower
#'and upper limit of quantile of reference read counts in each scan (sample).
#' @param scan_alt_qtile A numeric vector with length two specifying lower
#' and upper limit of quantile of alternative read counts in each
#' scan (sample).
#' @param snp_dp_qtile A numeric vector with length two specifying lower
#' and upper limit of quantile of total read counts in each SNP marker
#'@param snp_ref_qtile A numeric vector with length two specifying lower
#'and upper limit of quantile of reference read counts in each SNP marker.
#' @param snp_alt_qtile A numeric vector with length two specifying lower
#' and upper limit of quantile of alternative read counts in each SNP marker.
#' @param ... Unused.
#'
#' @details
#' `norm_dp_count`, `norm_ref_count`, and `norm_alt_count` use
#' normalized read counts which are obtained by dividing each read count
#' by the total read count of each sample. `scan_dp_qtile`, `scan_ref_qtile`,
#' and `scan_alt_qtile` work similarly but use quantile values of read counts
#' of each sample to decide the lower and upper limit of read counts.
#' This function generate two new nodes in the GDS file linked with
#' the given [GbsrGenotypeData] object. The new nodes "filt.data"
#' in the AD node
#' and "filt.genotype" contains read count data and genotype data
#' after filtering, respectively.
#' To reset the filter applied by setCallFilter(), run [setRawGenotype()].
#'
#' @return A [GbsrGenotypeData] object with filters on genotype calls.
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds <- loadGDS(gds_fn)
#'
#' # Filter out genotype calls supported by less than 5 reads.
#' gds <- setCallFilter(gds, dp_count = c(5, Inf))
#'
#' # Filter out genotype calls supported by reads less than
#' # the 20 percentile of read counts per marker in each sample.
#' gds <- setCallFilter(gds, scan_dp_qtile = c(0.2, 1))
#'
#' # Reset the filter
#' gds <- setRawGenotype(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("setCallFilter", function(object,
                                     dp_count = c(0, Inf),
                                     ref_count = c(0, Inf),
                                     alt_count = c(0, Inf),
                                     norm_dp_count = c(0, Inf),
                                     norm_ref_count = c(0, Inf),
                                     norm_alt_count = c(0, Inf),
                                     scan_dp_qtile = c(0, 1),
                                     scan_ref_qtile = c(0, 1),
                                     scan_alt_qtile = c(0, 1),
                                     snp_dp_qtile = c(0, 1),
                                     snp_ref_qtile = c(0, 1),
                                     snp_alt_qtile = c(0, 1),
                                     ...)
           standardGeneric("setCallFilter"))


#' Filter out scans (samples)
#'
#' Search samples which do not meet the criteria and label them as "invalid".
#'
#' @param object A [GbsrGenotypeData] object.
#' @param id A vector of strings matching with scan ID which can
#' be retrieve by `getScanID()`. The samples with the specified IDs
#' will be filtered out.
#' @param missing A numeric value \[0-1\] to specify the maximum
#' missing genotype call rate per sample.
#' @param het A vector of two numeric values \[0-1\] to specify
#' the minimum and maximum heterozygous genotype call rate per sample.
#' @param mac A integer value to specify the minimum minor
#' allele count per sample.
#' @param maf A numeric value to specify the minimum minor
#' allele frequency per sample.
#' @param ad_ref A numeric vector with length two specifying lower
#' and upper limit of reference read counts per sample.
#' @param ad_alt A numeric vector with length two specifying lower
#' and upper limit of alternative read counts per sample.
#' @param dp A numeric vector with length two specifying lower
#' and upper limit of total read counts per sample.
#' @param mean_ref A numeric vector with length two specifying lower
#' and upper limit of mean of reference read counts per sample.
#' @param mean_alt A numeric vector with length two specifying lower
#' and upper limit of mean of alternative read counts per sample.
#' @param sd_ref A numeric value specifying the upper limit of
#' standard deviation of reference read counts per sample.
#' @param sd_alt A numeric value specifying the upper limit of
#' standard deviation of alternative read counts per sample.
#' @param ... Unused.
#'
#' @details
#' For `mean_ref`, `mean_alt`, `sd_ref`, and `sd_alt`,
#' this function calculate mean and standard deviation of reads
#' obtained at SNP markers of each sample. If a mean read counts
#' of a sample was smaller than the specified lower limit or larger
#' than the upper limit, this function labels the sample as "invalid".
#' In the case of `sd_ref` and `sd_alt`, standard deviations of read counts
#' of each sample are checked and the samples having
#' a larger standard deviation will be labeled as "invalid".
#' To check valid and invalid samples, run [getValidScan()].
#'
#' @return A [GbsrGenotypeData] object with filters on scan(samples).
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the information needed for filtering.
#' gds <- countGenotype(gds)
#' gds <- countRead(gds)
#'
#' gds <- setScanFilter(gds,
#'                        id = getScanID(gds)[1:10],
#'                        missing = 0.2,
#'                        dp = c(5, Inf))
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("setScanFilter", function(object,
                                     id = "",
                                     missing = 1,
                                     het = c(0, 1),
                                     mac = 0,
                                     maf = 0,
                                     ad_ref = c(0, Inf),
                                     ad_alt = c(0, Inf),
                                     dp = c(0, Inf),
                                     mean_ref = c(0, Inf),
                                     mean_alt = c(0, Inf),
                                     sd_ref = Inf,
                                     sd_alt = Inf,
                                     ...)
           standardGeneric("setScanFilter"))


#' Filter out markers
#'
#' Search markers which do not meet the criteria and label them as "invalid".
#'
#' @param object A [GbsrGenotypeData] object.
#' @param id A vector of integers matching with snp ID which can
#' be retrieve by `getSnpID()`. The markers with the specified IDs
#' will be filtered out.
#' @param missing A numeric value \[0-1\] to specify
#' the maximum missing genotype call rate per marker
#' @param het A numeric vector with length two \[0-1\] to specify
#' the minimum and maximum heterozygous genotype call rate per marker
#' @param mac A integer value to specify the minimum minor allele
#' count per marker
#' @param maf A numeric value to specify the minimum minor allele
#'  frequency per marker.
#' @param ad_ref A numeric vector with length two specifying lower
#' and upper limit of reference read counts per marker.
#' @param ad_alt A numeric vector with length two specifying lower
#' and upper limit of alternative read counts per marker.
#' @param dp A numeric vector with length two specifying lower
#' and upper limit of total read counts per marker.
#' @param mean_ref A numeric vector with length two specifying lower
#' and upper limit of mean of reference read counts per marker.
#' @param mean_alt A numeric vector with length two specifying lower
#' and upper limit of mean of alternative read counts per marker.
#' @param sd_ref A numeric value specifying the upper limit of
#' standard deviation of reference read counts per marker.
#' @param sd_alt A numeric value specifying the upper limit of
#' standard deviation of alternative read counts per marker.
#' @param ... Unused.
#'
#' @details
#' For `mean_ref`, `mean_alt`, `sd_ref`, and `sd_alt`, this function
#' calculate mean and standard deviation of reads obtained for samples
#' at each SNP marker. If a mean read counts of a marker was smaller
#' than the specified lower limit or larger than the upper limit,
#' this function labels the marker as "invalid". In the case of `sd_ref`
#' and `sd_alt`, standard deviations of read counts of each marker are
#' checked and the markers having a larger standard deviation will be
#' labeled as "invalid". To check valid and invalid
#' markers, run [getValidSnp()].
#'
#' @return A [GbsrGenotypeData] object with filters on markers.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the information needed for filtering.
#' gds <- countGenotype(gds)
#' gds <- countRead(gds)
#'
#' gds <- setSnpFilter(gds,
#'                       id = getSnpID(gds)[1:100],
#'                       missing = 0.2,
#'                       dp = c(5, Inf))
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("setSnpFilter", function(object,
                                    id = NA_integer_,
                                    missing = 1,
                                    het = c(0, 1),
                                    mac = 0,
                                    maf = 0,
                                    ad_ref = c(0, Inf),
                                    ad_alt = c(0, Inf),
                                    dp = c(0, Inf),
                                    mean_ref = c(0, Inf),
                                    mean_alt = c(0, Inf),
                                    sd_ref = Inf,
                                    sd_alt = Inf,
                                    ...)
           standardGeneric("setSnpFilter"))


#' Filter out markers based on marker quality metrics
#'
#' A VCF file usually has marker quality metrics in
#' the INFO filed and those are stored in
#' a GDS file created via [GBScleanR]. This function filter
#' out markers based on those marker
#' quality metrics.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param mq A numeric value to specify minimum mapping
#' quality (shown as MQ in the VCF format).
#' @param fs A numeric value to specify maximum Phred-scaled
#' p-value (strand bias) (shown as FS in the VCF format).
#' @param qd A numeric value to specify minimum Variant
#' Quality by Depth (shown as QD in the VCF format).
#' @param sor A numeric value to specify maximum
#' Symmetric Odds Ratio (strand bias) (shown as SOR in the VCF format).
#' @param mqranksum A numeric values to specify the lower
#' and upper limit of Alt vs. Ref read mapping
#' qualities (shown as MQRankSum in the VCF format).
#' @param readposranksum A numeric values to specify the lower
#' and upper limit of Alt vs. Ref read position
#' bias (shown as ReadPosRankSum in the VCF format).
#' @param baseqranksum A numeric values to specify the lower
#' and upper limit of Alt Vs. Ref base
#' qualities (shown as BaseQRankSum in the VCF format).
#' @param ... Unused.
#'
#' @details
#' Detailed explanation of each metric can be found
#' in [GATK's web site](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @return A [GbsrGenotypeData] object with filters on markers.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' gds <- setInfoFilter(gds, mq = 40, qd = 20)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#' @export
#'
setGeneric("setInfoFilter", function(object,
                                     mq = 0 ,
                                     fs = Inf,
                                     qd = 0,
                                     sor = Inf,
                                     mqranksum = c(-Inf, Inf),
                                     readposranksum = c(-Inf, Inf),
                                     baseqranksum = c(-Inf, Inf),
                                     ...)
    standardGeneric("setInfoFilter"))


#' Reset the filter made by [setScanFilter()]
#'
#' Remove "invalid" labels put on samples and make all samples valid.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object after removing
#' all filters on scan(samples).
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize the information needed for filtering.
#' gds <- countGenotype(gds)
#' gds <- countRead(gds)
#'
#' gds <- setScanFilter(gds,
#'                        id = getScanID(gds)[1:10],
#'                        missing = 0.2,
#'                        dp = c(5, Inf))
#'
#' # Reset all filters applied above.
#' gds <- resetScanFilters(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("resetScanFilters", function(object, ...)
    standardGeneric("resetScanFilters"))


#' Reset the filter made by [setSnpFilter()]
#'
#' Remove "invalid" labels put on markers and make all markers valid.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object after removing all filters on markers.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Check the number of markers.
#' nsnp(gds)
#'
#' # Summarize the information needed for filtering.
#' gds <- countGenotype(gds)
#' gds <- countRead(gds)
#'
#' # filter out some markers meeting the criteria.
#' gds <- setSnpFilter(gds,
#'                       id = getSnpID(gds)[1:100],
#'                       missing = 0.2,
#'                       dp = c(5, Inf))
#'
#' # Check the number of the retained markers.
#' nsnp(gds)
#'
#' # Reset all filters applied above.
#' gds <- resetSnpFilters(gds)
#'
#' # Check the number of the markers again.
#' nsnp(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#' @export
#'
setGeneric("resetSnpFilters", function(object, ...)
    standardGeneric("resetSnpFilters"))


#' Reset all filters made by [setScanFilter()], [setSnpFilter()],
#' and [setCallFilter()].
#'
#' Return all data intact.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object after removing all filters.
#'
#'
#' @return A [GbsrGenotypeData] object after removing all filters on markers.
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds <- loadGDS(gds_fn)
#'
#' # `setCallFilter()` do not require summarized information of
#' # genotype counts and read counts.
#' gds <- setCallFilter(gds, dp_count = c(5, Inf))
#'
#' # `setScanFilter()` and `setSnpFilter()` needs information of
#' # the genotype count summary and the read count summary.
#' gds <- countGenotype(gds)
#' gds <- countRead(gds)
#'
#' gds <- setScanFilter(gds,
#'                        id = getScanID(gds)[1:10],
#'                        missing = 0.2,
#'                        dp = c(5, Inf))
#'
#' gds <- setSnpFilter(gds,
#'                       id = getSnpID(gds)[1:100],
#'                       missing = 0.2,
#'                       dp = c(5, Inf))
#'
#' gds <- setInfoFilter(gds, mq = 40, qd = 20)
#'
#' # Reset all filters applied above.
#' gds <- resetFilters(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#' @export
#'
setGeneric("resetFilters", function(object, ...)
    standardGeneric("resetFilters"))


#' Remove markers potentially having redundant information.
#'
#' Markers within the length of the sequenced reads
#' (usually ~ 150 bp, up to your sequencer)
#' potentially have redundant information and
#' those will cause unexpected errors
#' in error correction which assumes
#' independency of markers each other.
#' This function only retains the first marker or
#' the least missing rate marker
#' from the markers locating within the specified stretch.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param range A integer value to indicate the stretch to search markers.
#' @param ... Unused.
#'
#' @details
#' This function search valid markers from the first marker
#' of each chromosome and
#' compare its physical position with a neighbor marker.
#' If the distance between those
#' markers are equal or less then `range`, one of them
#' which has a larger missing rate
#' will be removed (labeled as invalid marker).
#' When the first marker was retained and
#' the second marker was removed as invalid marker,
#' next the distance between the first marker
#' and the third marker will be checked and
#' this cycle is repeated until reaching the
#' end of each chromosome. Run [getValidSnp()]
#' to check the valid SNP markers.
#'
#' @return A [GbsrGenotypeData] object with filters on markers.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize genotype count information to be used in thinMarker().
#' gds <- countGenotype(gds)
#' gds <- thinMarker(gds, range = 150)
#'
#' closeGDS(gds) # Close the connection to the GDS file
#' @export
#'
setGeneric("thinMarker", function(object, range = 150, ...)
    standardGeneric("thinMarker"))


#' Create a GDS file with subset data of the current GDS file
#'
#' Create a new GDS file storing the subset
#' data from the current GDS file linked to
#' the given [GbsrGenotypeData] object with
#' keeping (or removing) information based on
#' valid markers and samples information.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param out_fn A string to specify the path to an output GDS file.
#' @param snp_incl A logical vector having the same length
#' with the total number of markers.
#' The values obtained via [getValidSnp()] are used.
#' @param scan_incl A logical vector having the same length
#' with the total number of scans (samples).
#' The values obtained via [getValidScan()] are used.
#' @param incl_parents A logical value to specify whether parental
#' samples should be included in a subset data or not.
#' @param verbose if TRUE, show information.
#' @param ... Unused.
#'
#' @details
#' A resultant subset data in a new GDS file includes subsets
#' of each category of data, e.g.
#' genotype, SNP ID, scan ID, read counts,
#' and quality metrics of SNP markers.
#' The connection to the GDS file of the input [GbsrGenotypeData] object will
#' be automatically closed for internal file handling in this function. Please
#' use [openGDS()] to open the connection again. If you use [loadGDS()],
#' summary statistics and filtering information will be discarded.
#'
#' @return A [GbsrGenotypeData] object linking to
#' the new GDS file storing subset data.
#'
#' @examples
#' #' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize genotype count information to be used in `setSnpFilter()`
#' gds <- countGenotype(gds)
#'
#' # Filter out markers meeting the criteia.
#' gds <- setSnpFilter(gds, missing = 0.2, het = c(0.1, 0.9), maf = 0.05)
#'
#' # Create a new GDS file with the subset data obtained by applying
#' # the filter maed via `setSnpFilter()`.
#' subsetgds_fn <- tempfile("sample_subset", fileext = ".gds")
#' subset_gds <- subsetGDS(gds, out_fn = subsetgds_fn)
#'
#' # Close the connection to the GDS files.
#' closeGDS(subset_gds)
#'
#' @export
#'
setGeneric("subsetGDS", function(object,
                                 out_fn = "./susbet.gds",
                                 snp_incl,
                                 scan_incl,
                                 incl_parents = TRUE,
                                 verbose = TRUE,
                                 ...)
    standardGeneric("subsetGDS"))


#' Set the filtered data to be used in GBScleanR's functions
#'
#' Set the "filt.genotype" node and the "filt.data" node as
#' primary nodes for genotype
#' data and read count data. The data stored
#' in the primary nodes are used in the
#' functions of GBScleanR.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param ... Unused.
#'
#' @details
#' A [GbsrGenotypeData] object storing information of
#' the primary node of genotype data and
#' read count data. All of the functions implemented
#' in [GBScleanR] check the primary nodes
#' and use data stored in those nodes. [setCallFilter()]
#' create new nodes storing
#' filtered genotype calls and read counts in a GDS file
#' and change the primary nodes to
#' "filt.genotype" and "filt.data" for genotype and
#' read count data, respectively.
#' [setRawGenotype()] set back the nodes to the original,
#' those are "genotype" and "data" for
#' genotype and read count data, respectively.
#' You can set the filtered data again by [setFiltGenotype()].
#'
#' @return A [GbsrGenotypeData] object.
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds <- loadGDS(gds_fn)
#'
#' # Filter out set zero to read counts and
#' # missing to genotype calls of which meet the criteria.
#' gds <- setCallFilter(gds, dp_count = c(5, Inf))
#'
#' # Now any functions of [GBScleanR] reference the genotype data
#' # stored in the "filt.genotype" node of the GDS file.
#'
#' # Close and reopen the connection to the GDS file.
#' closeGDS(gds)
#' gds <- loadGDS(gds)
#'
#' # If the connection was closed once and re-loaded,
#' # raw genotype data is referenced by the functions.
#'
#' # To set the "filt.genotype" node as genotype data again,
#' # run the setFiltGenotype().
#' gds <- setFiltGenotype(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("setFiltGenotype", function(object, ...)
    standardGeneric("setFiltGenotype"))


#' Set the origina; data to be used in GBScleanR's functions
#'
#' Set the "genotype" node and the "data" node
#' as primary nodes for genotype
#' data and read count data. The data stored
#' in the primary nodes are used in the
#' functions of GBScleanR.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param ... Unused.
#'
#' @details
#' A [GbsrGenotypeData] object storing information of
#' the primary node of genotype data and
#' read count data. All of the functions implemented
#' in [GBScleanR] check the primary nodes
#' and use data stored in those nodes.
#' [setCallFilter()] create new nodes storing
#' filtered genotype calls and read counts in
#' a GDS file and change the primary nodes to
#' "filt.genotype" and "filt.data" for genotype and
#' read count data, respectively.
#' [setRawGenotype()] set back the nodes to
#' the original, those are "genotype" and "data" for
#' genotype and read count data, respectively.
#' You can set the filtered data again by [setFiltGenotype()].
#'
#' @return A [GbsrGenotypeData] object.
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds <- loadGDS(gds_fn)
#'
#' # Filter out set zero to read counts and
#' # missing to genotype calls of which meet the criteria.
#' gds <- setCallFilter(gds, dp_count = c(5, Inf))
#'
#' # Now any functions of [GBScleanR] reference the genotype data
#' # stored in the "filt.genotype" node of the GDS file.
#'
#' # If you need to set the "genotype" node, where store the raw genotype data
#' # as genotype to be referenced by the functions of GBScleanR,
#' # run the following.
#' gds <- setRawGenotype(gds)
#'
#' # Reopening the connection to the GDS file also set the raw genotype again.
#' gds <- loadGDS(gds)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#' @export
#'
setGeneric("setRawGenotype", function(object, ...)
    standardGeneric("setRawGenotype"))


#' Genotype estimation using a hiden Morkov model
#'
#' Clean up genotype data by error correction based on
#' genotype estimation using a hidden Markov model.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param chr An integer vector of chromosome indices to
#' be analyzed. All chromosomes will be analyzed if you left it default.
#' @param recomb_rate A numeric value to indicate the expected
#' recombination frequency per chromosome per megabase pairs.
#' @param error_rate A numeric value of the expected sequence error rate.
#' @param call_threshold A numeric value of the probability threshold
#' to accept estimated genotype calls.
#' @param het_parent A logical value to indicate whether parental
#' samples are outbred or inbred. If FALSE, this function assume all
#' true genotype of markers in parents are homozygotes.
#' @param optim A logical value to specify whether to
#' conduct parameter optimization for
#' error correction.
#' @param iter An integer value to specify the number of
#' iterative parameter updates.
#' @param n_threads An integer value to specify the number of
#' threads used for the calculation. The default is `n_threads = NULL`
#' and automatically set half the number of available threads on the computer.
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object in which the "estimated.haplotype",
#' "corrected.genotype" and
#' "parents.genotype" nodes were added.
#'
#' @importFrom RcppParallel setThreadOptions defaultNumThreads
#' @importFrom expm expm.Higham08
#' @importFrom stats cor
#'
#' @examples

#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Find the IDs of parental samples.
#' parents <- grep("Founder", getScanID(gds), value = TRUE)
#'
#' # Set the parents and flip allele information
#' # if the reference sample (Founder1 in our case) has homozygous
#' # alternative genotype at some markers of which alleles will
#' # be swapped to make the reference sample have homozygous
#' # reference genotype.
#' gds <- setParents(gds, parents = parents, flip = TRUE)
#'
#' # Initialize a scheme object stored in the slot of the GbsrGenotypeData.
#' # We chose `crosstype = "pair"` because two inbred founders were mated
#' # in our breeding scheme.
#' # We also need to specify the mating matrix which has two rows and
#' # one column with integers 1 and 2 indicating a sample (founder)
#' # with the memberID 1 and a sample (founder) with the memberID 2
#' # were mated.
#' gds <- initScheme(gds, crosstype = "pair", mating = cbind(c(1:2)))
#'
#' # Add information of the next cross conducted in our scheme.
#' # We chose 'crosstype = "selfing"', which do not require a
#' # mating matrix.
#' gds <- addScheme(gds, crosstype = "selfing")
#'
#' # Execute error correction by estimating genotype and haplotype of
#' # founders and offspring.
#' gds <- estGeno(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
setGeneric("estGeno", function(object,
                               chr,
                               recomb_rate = 0.04,
                               error_rate = 0.0025,
                               call_threshold = 0.9,
                               het_parent = FALSE,
                               optim = TRUE,
                               iter = 2,
                               n_threads = NULL,
                               ...)
    standardGeneric("estGeno"))

#' Add genotype data into the GDS file
#'
#' @param object A [GbsrGenotypeData] object.
#' @param id A character vector.
#' @param genotype A numeric vector or matrix indicating genotypes at markers
#' of given samples. The length or the number of columns should match with the
#' number of markers recorded in the GDS file, which can be obtained via
#' nsnp() with `valid = FALSE`.
#' @param reads A numeric vector of matrix indicating read counts at markers of
#' given samples. The length or the number of columns should match with twice
#' the number of markers recorded in the GDS file, which can be obtained via
#' nsnp() with `valid = FALSE`.
#' @param ... Unused.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' id <- "Dummy1"
#' genotype <- sample(c(0:3), nsnp(gds), replace = TRUE)
#' reads <- round(rexp(nsnp(gds) * 2, 1/5))
#' missing_pos <- which(genotype == 3)
#' reads[c(missing_pos, missing_pos + 1)] <- 0
#' ref_pos <- which(genotype == 2)
#' reads[ref_pos + 1] <- 0
#' alt_pos <- which(genotype == 0)
#' reads[alt_pos] <- 0
#'
#' addScan(gds, id, genotype, reads)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' @export
setGeneric("addScan", function(object, id, genotype, reads, ...)
    standardGeneric("addScan"))


#' Build a [GbsrScheme] object
#'
#' [GBScleanR] uses breeding scheme information to set the expected
#' number of cross overs in a chromosome which is a required parameter
#' for the genotype error correction with the hidden Markov model
#' implemented in the [estGeno()] function.
#' This function build the object storing
#' type crosses performed at each generation of breeding and population sizes.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param crosstype A string to indicate the type of
#' cross conducted with a given generation.
#' @param mating An integer matrix to indicate mating combinations.
#' The each element should match with IDs of
#' parental samples which are 1 to N. see Details.
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object storing
#' a [GbsrScheme] object in the "scheme" slot.
#'
#' @details
#' A [GbsrScheme] object stores information of a population size,
#' mating combinations and
#' a type of cross applied to each generation of the breeding process
#' to generate the population which you are going to subject
#' to the [estGeno()] function.
#' The first generation should be parents of the population.
#' It is supposed that
#' [setParents()] has been already executed and parents
#' are labeled in the
#' [GbsrGenotypeData] object. The number of parents
#' are automatically recognized.
#' The "crosstype" of the first generation can be
#' "pairing" or "random" with
#' `pop_size = N`, where N is the number of parents.
#' You need to specify a matrix indicating combinations
#' of `mating`, in which each column shows
#' a pair of parental samples. For example, if you have
#' only two parents, the `mating` matrix
#' is `mating = cbind(c(1:2))`. The indices used in the matrix
#' should match with the IDs labeled to parental samples by [setParents()].
#' The created [GbsrScheme] object is set
#' in the `scheme` slot of the [GbsrGenotypeData] object.
#'
#' @export
#'
#' @seealso [addScheme()] and [showScheme()]
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Biparental F2 population.
#' gds <- setParents(gds, parents = c("Founder1", "Founder2"))
#'
#' # setParents gave member ID 1 and 2 to Founder1 and Founder2, respectively.
#' gds <- initScheme(gds, crosstype = "pair", mating = cbind(c(1:2)))
#'
#' # Now the progenies of the cross above have member ID 3.
#' # If `crosstype = "selfing"` or `"sibling"`, you can omit a `mating` matrix.
#' gds <- addScheme(gds, crosstype = "self")
#'
#' # Now you can execute `estGeno()` which requires a [GbsrScheme] object.
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
setGeneric("initScheme", function(object, crosstype, mating, ...)
    standardGeneric("initScheme"))


#' #' Build a [GbsrScheme] object
#'
#' [GBScleanR] uses breeding scheme information to set the expected
#' number of cross overs in a chromosome which is a required parameter
#' for the genotype error correction with the Hidden Markov model
#' implemented in the `estGeno()` function.
#' This function build the object storing
#' type crosses performed at each generation of breeding and population sizes.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param crosstype A string to indicate the type of
#' cross conducted with a given generation.
#' @param mating An integer matrix to indicate mating combinations.
#' The each element should match with member IDs of the last generation.
#' @param pop_size An integer of the number of
#' individuals in a given generation.
#' @param ... Unused.
#'
#' @return A [GbsrGenotypeData] object storing
#' a [GbsrScheme] object in the "scheme" slot.
#'
#' @details
#' A scheme object is just a data.frame indicating a population size and
#' a type of cross applied to each generation of the breeding process
#' to generate the population which you are going to subject
#' to the [estGeno()] function.
#' The `crosstype` can take either of "selfing", "sibling",
#' "pairing", and "random".
#' When you set `crosstype = "random"`, you need to
#' specify `pop_size` to indicate how many
#' individuals were crossed in the random mating.
#' You also need to specify a matrix indicating
#' combinations of `mating`, in which
#' each column shows a pair of member IDs indicating
#' parental samples of the cross.
#' Member IDs are serial numbers starts from 1 and
#' automatically assigned by
#' [initScheme()] and [addScheme()]. To check the member IDs,
#' run [showScheme()].
#' Please see the examples section for more details of
#' specifying a `mating` matrix.
#' The created [GbsrScheme] object is set in the `scheme`
#' slot of the [GbsrGenotypeData] object.
#'
#' @export
#'
#' @seealso [addScheme()] and [showScheme()]
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Biparental F2 population.
#' gds <- setParents(gds, parents = c("Founder1", "Founder2"))
#'
#' # setParents gave member ID 1 and 2 to Founder1 and Founder2, respectively.
#' gds <- initScheme(gds, crosstype = "pair", mating = cbind(c(1:2)))
#'
#' # Now the progenies of the cross above have member ID 3.
#' # If `crosstype = "selfing"` or `"sibling"`, you can omit a `mating` matrix.
#' gds <- addScheme(gds, crosstype = "self")
#'
#' ############################################################################
#' # Now you can execute `estGeno()` which requires a [GbsrScheme] object.
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
setGeneric("addScheme", function(object, crosstype, mating, pop_size, ...)
    standardGeneric("addScheme"))


#' Show the information stored in a [GbsrScheme] object
#'
#' Print the information of each generation in
#' a [GbsrScheme] object in the scheme
#' slot of a [GbsrGenotypeData] object.
#' A [GbsrScheme] object stores information of a population size,
#' mating combinations and
#' a type of cross applied to each generation of the breeding process
#' to generate the population which you are going to
#' subject to the `estGeno()` function.
#'
#' @param object A [GbsrGenotypeData] object.
#' @param ... Unused.
#'
#' @return NULL. Print the scheme information on the R console.
#'
#' @export
#'
#' @seealso [initScheme()] and [addScheme()]
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Biparental F2 population.
#' gds <- setParents(gds, parents = c("Founder1", "Founder2"))
#'
#' # setParents gave member ID 1 and 2 to Founder1 and Founder2, respectively.
#' gds <- initScheme(gds, crosstype = "pair", mating = cbind(c(1:2)))
#'
#' # Now the progenies of the cross above have member ID 3.
#' # If `crosstype = "selfing"` or `"sibling"`, you can omit a `mating` matrix.
#' gds <- addScheme(gds, crosstype = "self")
#'
#' # Now you can execute `estGeno()` which requires a [GbsrScheme] object.
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
setGeneric("showScheme", function(object, ...)
    standardGeneric("showScheme"))
