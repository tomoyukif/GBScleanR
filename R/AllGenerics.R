setGeneric("openGDS", function(object, ...)
  standardGeneric("openGDS"))


#' Check if a GDS file has been opened or not.
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#'
#' @return
#' `TRUE` if the GDS file linked to the input GbsrGenotypeData object has been opened, while `FALSE` if closed.
#'
#'@export
#'
setGeneric("isOpenGDS", function(object, ...)
  standardGeneric("isOpenGDS"))


#' Close the connection to the GDS file
#'
#' Close the connection to the GDS file linked to the given GbsrGenotypeData object.
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#'@export
#'
setGeneric("closeGDS", function(object, ...)
  standardGeneric("closeGDS"))


#' Write out the information stored in the SnpAnnotationDataSet slot
#'
#' All the data stored in the SnpAnnotatoinDataSet slot of the GbsrGenotypeData
#' object can be saved in the GDS file linked to the given GbsrGenotypeData object.
#' You can load the saved data using [loadSnpAnnot()].
#'
#' @param object A GbsrGenotypeData object
#' @param ... Unused.
#'
#' @export
#'
setGeneric("saveSnpAnnot", function(object, ...)
  standardGeneric("saveSnpAnnot"))


#' Write out the information stored in the ScanAnnotationDataSet slot
#'
#' All the data stored in the ScanAnnotationDataSet slot of the GbsrGenotypeData
#' object can be saved in the GDS file linked to the given GbsrGenotypeData object.
#' You can load the saved data using [loadSnpAnnot()].
#'
#' @param object A GbsrGenotypeData object
#' @param ... Unused.
#' 
#' @export
#'
setGeneric("saveScanAnnot", function(object, ...)
  standardGeneric("saveScanAnnot"))


#' Load the stored SnpAnnotationDataSet information
#'
#' All the data stored in the SnpAnnotatoinDataSet slot of the GbsrGenotypeData
#' object can be saved in the GDS file linked to the given GbsrGenotypeData object via [saveSnpAnnot()].
#' You can load the saved data using this function.
#'
#' @param object A GbsrGenotypeData object
#' @param ... Unused.
#' 
#' @export
#'
setGeneric("loadSnpAnnot", function(object, ...)
  standardGeneric("loadSnpAnnot"))


#' Load the stored ScanAnnotationDataSet information
#'
#' All the data stored in the ScanAnnotationDataSet slot of the GbsrGenotypeData
#' object can be saved in the GDS file linked to the given GbsrGenotypeData object via [saveScanAnnot()].
#' You can load the saved data using this function.
#'
#' @param object A GbsrGenotypeData object
#' @param ... Unused.
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
#' @param object A GbsrGenotypeData object.
#' @param out_fn A string to specify the path to an output VCF file.
#' @param node Either one of "raw", "filt", and "cor" to output raw genotype data, filtered genotype data, or corrected genotype data, respectively.
#' @param valid A logical value to specify whether to output valid markers and samples only or all.
#' @param out_fmt A character vector to specify which variables in the annotation/format node should be output.
#' @param out_info A character vector to specify which variables in the annotation/info node should be output.
#' @param ... Unused.
#' @export
#'
#' @importFrom SeqArray seqSNP2GDS seqGDS2VCF
#'
#' @examples
#' \dontrun{
#' gdata <- loadGDS("/path/to/GDS.gds")
#' gdata <- clean(gdata)
#' gbsrGDS2VCF(gdata, "/path/to/output.vcf", node = "cor")
#' }
#'
setGeneric("gbsrGDS2VCF", function(object, out_fn, node = "raw", valid = TRUE,
                                   out_fmt = NULL, out_info = NULL, ...)
  standardGeneric("gbsrGDS2VCF"))


#' Return the number of SNPs.
#'
#' This function returns the number of SNPs recorded in the GDS file
#' connected to the given GbsrGenotypeData object.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' If `valid = TRUE`, the number of SNPs which are labeled `TRUE` in
#' the SnpAnnotationDataSet slot will be returned. You need the number
#' of over all SNPs, set `valid = FALSE`. [getValidSnp()] tells you
#' which markers are valid.
#'
#' @seealso [getValidSnp()]
#'
#' @export
#'
setGeneric("nsnp", function(object, valid = TRUE, ...)
  standardGeneric("nsnp"))


#' Return the number of scans (samples).
#'
#' This function returns the number of samples recorded in the GDS file
#' connected to the given GbsrGenotypeData object.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' If `valid = TRUE`, the number of samples which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. You need
#' the number of over all samples, set `valid = FALSE`.
#' [getValidSnp()] tells you which samples are valid.
#'
#' @seealso [getValidSnp()]
#'
#'@export
#'
setGeneric("nscan", function(object, valid = TRUE, ...)
  standardGeneric("nscan"))


#' Return a logical vector indicating which are valid SNP markers.
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#' @seealso [setValidSnp()]
#'
#' @export
#'
setGeneric("getValidSnp", function(object, ...)
  standardGeneric("getValidSnp"))


#' Return a logical vector indicating which are valid scans (samples).
#'
#' @param object A GbsrGenotypeData object.
#' @param parents A logical value to indicate to set FALSE or TRUE to parental samples. If you specify `parents = "only"`, this function returns a logical vector indicating TRUE for only parental samples.
#' @param ... Unused.
#' 
#' @seealso [setValidScan()]
#'
#' @export
#'
setGeneric("getValidScan", function(object, parents = FALSE, ...)
  standardGeneric("getValidScan"))


#' Manually set valid SNP markers.
#'
#' If you need manually set valid and invalid SNP markers, you can do it via this function,
#' e.g in the case you conducted a filtering on SNP markers manually by your self.
#'
#' @param object A GbsrGenotypeData object.
#' @param new A logical vector of the same length with the over all number of the SNP markers.
#' @param update A logical vector of the same length with the currently valid SNP markers.
#' @param ... Unused.
#' 
#' @details
#' To over write the current validity information, give a logical vector to `new`.
#' On the other hand, a logical vector specified to `update` will be used to
#' update validity information of the currently valid SNP markers. If you gave
#' a vector for both argument, only the vector passed to `new` will be used to
#' over write the validity information.
#'
#' @seealso [setSnpFilter()] to filter out SNP markers based on some summary statistics.
#'
#' @export
#'
setGeneric("setValidSnp", function(object, new, update, ...)
  standardGeneric("setValidSnp"))


#' Manually set valid scans (samples).
#'
#' If you need manually set valid and invalid samples, you can do it via this function,
#' e.g in the case you conducted a filtering on samples manually by your self.
#'
#' @param object A GbsrGenotypeData object.
#' @param new A logical vector of the same length with the over all number of the samples.
#' @param update A logical vector of the same length with the currently valid samples.
#' @param ... Unused.
#' 
#' @details
#' To over write the current validity information, give a logical vector to `new`.
#' On the other hand, a logical vector specified to `update` will be used to
#' update validity information of the currently valid samples. If you gave
#' a vector for both argument, only the vector passed to `new` will be used to
#' over write the validity information.
#'
#' @seealso [setScanFilter()] to filter out samples based on some summary statistics
#'
#' @export
#'
setGeneric("setValidScan", function(object, new, update, ...)
  standardGeneric("setValidScan"))


#' Get a logical vector indicating flipped SNP markers.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value to specify flip alleles only of valid markers. see also[setSnpFilter()].
#' @param ... Unused.
#' 
#' @details
#' Flipped markers are markers where the alleles expected as reference allele are called as
#' alternative allele. If you specify two parents in the `parents` argument of
#' [setParents()] with `flip = TRUE`, `bi = TRUE`, and `homo = TRUE`, the alleles found
#' in the parent specified as the first element to the `parents` argument are supposed as
#' reference alleles of the markers. If the "expected" reference alleles are not actually
#' called as reference alleles but alternative alleles in the given data. [setParents()] will
#' automatically labels those markers "flipped". The SnpAnnotatoinDataSet slot sores this
#' information and accessible via [getFlipped()] which gives you a logical vector
#' indicating which markers are labeled as flipped `TRUE` or not flipped `FALSE`.
#' [haveFlipped()] just tells you whether the SnpAnnotatoinDataSet slot has
#' the information of flipped markers or not.
#'
#' @seealso [setParents()] and [haveFlipped()].
#'
#' @export
#'
setGeneric("getFlipped", function(object, valid = TRUE, ...)
  standardGeneric("getFlipped"))


#' Get a logical value indicating flipped SNP markers whether information exists.
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#' @details
#' Flipped markers are markers where the alleles expected as reference allele are called as
#' alternative allele. If you specify two parents in the `parents` argument of
#' [setParents()] with `flip = TRUE`, `bi = TRUE`, and `homo = TRUE`, the alleles found
#' in the parent specified as the first element to the `parents` argument are supposed as
#' reference alleles of the markers. If the "expected" reference alleles are not actually
#' called as reference alleles but alternative alleles in the given data. [setParents()] will
#' automatically labels those markers "flipped". The SnpAnnotatoinDataSet slot sores this
#' information and accessible via [getFlipped()] which gives you a logical vector
#' indicating which markers are labeled as flipped `TRUE` or not flipped `FALSE`.
#' [haveFlipped()] just tells you whether the SnpAnnotatoinDataSet slot has
#' the information of flipped markers or not.
#'
#' @seealso [setParents()] and [getFlipped()].
#'
#' @export
#'
setGeneric("haveFlipped", function(object, ...)
  standardGeneric("haveFlipped"))


#' Flip genotype, allele information, and allele read counts of the flipped SNP markers.
#'
#' Genotype data, allele information, and allele read counts to be the expected reference
#' allele as actually reference allele.
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#' @details
#' Flipped markers are markers where the alleles expected as reference allele are called as
#' alternative allele. If you specify two parents in the `parents` argument of
#' [setParents()] with `flip = TRUE`, `bi = TRUE`, and `homo = TRUE`, the alleles found
#' in the parent specified as the first element to the `parents` argument are supposed as
#' reference alleles of the markers. If the "expected" reference alleles are not actually
#' called as reference alleles but alternative alleles in the given data. [setParents()] will
#' automatically labels those markers "flipped". The SnpAnnotatoinDataSet slot sores this
#' information and accessible via [getFlipped()] which gives you a logical vector
#' indicating which markers are labeled as flipped `TRUE` or not flipped `FALSE`.
#' [haveFlipped()] just tells you whether the SnpAnnotatoinDataSet slot has
#' the information of flipped markers or not.
#'
#' @seealso [setParents()], [getFlipped()], and [haveFlipped()].
#'
#' @export
#'
setGeneric("flipData", function(object, ...)
  standardGeneric("flipData"))


#' Get genotype call data.
#'
#' Genotype calls are retrieved from the GDS file linked to the given
#' GbsrGenotypeData object.
#'
#' @param object A gbsrGenotypeData object.
#' @param chr A integer vector of indexes indicating chromosomes to get read count data.
#' @param node Either of "raw", "filt", and "cor. See details.
#' @param parents A logical value or "only" to include data for parents or to get data only for parents.
#' @param ... Unused.
#' 
#' @details
#' Genotype call data can be obtained from the "genotype" node, the "filt.genotype"
#' node, or the "corrected.genotype" node of the GDS file with `node = "raw"`,
#' `node = "filt"`, or `node = "raw"`, respectively. If `node = "parents`, the data in the "parents.genotype" node will be returned. The "parents.genotype" node stores phased parental genotypes estimated by the [clean()] function.
#' The [setCallFilter()] function generate filtered genotype call data in the
#' "filt.genotype" node which can be accessed as mentioned above. On the other hand, the
#' "corrected.genotype" node can be generated via the [clean()] function.
#'
#' @seealso [setCallFilter()] and [clean()]
#'
#' @export
#'
setGeneric("getGenotype", function(object,
                                   chr = NULL,
                                   node = "raw",
                                   parents = FALSE,
                                   ...)
  standardGeneric("getGenotype"))


#' Get haplotype call data.
#'
#' Haplotype calls are retrieved from the GDS file linked to the given
#' GbsrGenotypeData object.
#'
#' @param object A gbsrGenotypeData object.
#' @param chr A integer vector of indexes indicating chromosomes to get read count data.
#' @param parents A logical value or "only" to include data for parents or to get data only for parents.
#' @param ... Unused.
#' 
#' @details
#' Haplotype call data can be obtained from the "estimated.haplotype" node of
#' the GDS file which can be generated via the [clean()] function. Thus, this function
#' is valid only after having executed [clean()].
#'
#' @seealso [clean()]
#'
#' @export
#'
setGeneric("getHaplotype", function(object,
                                    chr = NULL,
                                    parents = FALSE,
                                    ...)
  standardGeneric("getHaplotype"))


#' Get read count data.
#'
#' Read counts for reference allele and alternative allele are retrieved
#' from the GDS file linked to the given GbsrGenotypeData object.
#'
#' @param object A gbsrGenotypeData object.
#' @param chr A integer vector of indexes indicating chromosomes to get read count data.
#' @param node Either of "raw" and "filt". See details.
#' @param parents A logical value or "only" to include data for parents or to get data only for parents.
#' @param ... Unused.
#' 
#' @details
#' Read count data can be obtained from the "annotation/format/AD/data" node or the
#' "annotation/format/AD/filt.data" node of the GDS file with `node = "raw"` or
#' `node = "filt"`, respectively. The [setCallFilter()] function generate filtered
#' read count data in the "annotation/format/AD/filt.data" node which can be accessed as
#' mentioned above.
#'
#' @seealso [setCallFilter()]
#'
#' @export
#'
setGeneric("getRead", function(object,
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
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#' @param levels A logical value. See details.
#' @param name A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' A GDS file created via GBScleanR stores chromosome names as sequential integers
#' from 1 to N, where N is the number of chromosomes. This function returns those
#' indexes as default. If you need actual names of the chromosomes, set `name = TRUE`.
#' `levels = TRUE` gives you only unique chromosome names with length N.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getPosition", function(object, valid = TRUE, ...)
  standardGeneric("getPosition"))


#' Obtain reference allele information of each SNP marker
#'
#' This function returns reference alleles, either of A, T, G, and C, of SNP markers.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getAlleleA", function(object, valid = TRUE, ...)
  standardGeneric("getAlleleA"))


#' Obtain alternative allele information of each SNP marker
#'
#' This function returns alternative alleles, either of A, T, G, and C, of SNP markers.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getAlleleB", function(object, valid = TRUE, ...)
  standardGeneric("getAlleleB"))


#' Obtain SNP ID
#'
#' This function returns SNP ID of SNP markers.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getSnpID", function(object, valid = TRUE, ...)
  standardGeneric("getSnpID"))


#' Obtain scan (sample) ID
#'
#' This function returns scan (sample) ID.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getScanID", function(object, valid = TRUE, ...)
  standardGeneric("getScanID"))


#' Obtain ploidy information of each SNP marker
#'
#' This function returns ploidy of each SNP marker. The ploidy of all the markers in
#' a dataset is a same value and the current implementation of GBScleanR only works
#' with data having ploidy = 2 for all markers.
#'
#' @param object A GbsrGenotypeData object.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getPloidy", function(object, valid = TRUE, ...)
  standardGeneric("getPloidy"))


#' Obtain information stored in the "annotation/info" node
#'
#' The "annotation/info" node stores annotation infromation of markers obtained
#' via SNP calling tools like bcftools and GATK.
#'
#' @param object A GbsrGenotypeData object.
#' @param var A string to indicate which annotation info should be retrieved.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getInfo", function(object, var, valid = TRUE, ...)
  standardGeneric("getInfo"))


#' Obtain total reference read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param prop A logical value whether to return values as proportions of total reference read counts in total read counts per SNP or not. 
#' @param ... Unused.
#' 
#'
#' @details
#' You need to execute [countRead()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param prop A logical value whether to return values as proportions of total alternative read counts in total read counts per SNP or not.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [countRead()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [countRead()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total reference genotype counts to total non missing genotype counts or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total heterozygote counts to total non missing genotype counts or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total alternative genotype counts to total non missing genotype counts or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total missing genotype counts to the total genotype calls or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total reference allele counts to total non missing allele counts or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total alternative allele counts to total non missing allele counts or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param prop A logical value whether to return values as proportions of total missing allele counts to the total allele number or not.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getCountAlleleMissing", function(object,
                                             target = "snp",
                                             valid = TRUE,
                                             prop = FALSE,
                                             ...)
  standardGeneric("getCountAlleleMissing"))


#' Obtain mean values of total reference read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getMeanReadRef", function(object,
                                      target = "snp",
                                      valid = TRUE,
                                      ...)
  standardGeneric("getMeanReadRef"))


#' Obtain mean values of total alternative read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getMeanReadAlt", function(object,
                                      target = "snp",
                                      valid = TRUE,
                                      ...)
  standardGeneric("getMeanReadAlt"))


#' Obtain standard deviations of total reference read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getSDReadRef", function(object,
                                    target = "snp",
                                    valid = TRUE,
                                    ...)
  standardGeneric("getSDReadRef"))


#' Obtain standard deviations of total alternative read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getSDReadAlt", function(object,
                                    target = "snp",
                                    valid = TRUE,
                                    ...)
  standardGeneric("getSDReadAlt"))


#' Obtain quantile values of total reference read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param q A numeric value [0-1] to indicate quantile to obtain.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
#'
#' @export
#'
setGeneric("getQtileReadRef", function(object,
                                       target = "snp",
                                       q = 0.5,
                                       valid = TRUE,
                                       ...)
  standardGeneric("getQtileReadRef"))


#' Obtain quantile values of total alternative read counts per SNP or per scan (sample)
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param q A numeric value [0-1] to indicate quantile to obtain.
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [calcReadStats()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param valid A logical value. See details.
#' @param ... Unused.
#' 
#' @details
#' You need to execute [countGenotype()] to calculate sumaary statisticsto be
#' obtained via this function.
#' If `valid = TRUE`, the chromosome information of markers which are labeled `TRUE`
#' in the ScanAnnotationDataSet slot will be returned. [getValidSnp()] tells you
#' which samples are valid.
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
#' This function calculates several summary statistics of genotype calls and alleles
#' per marker and per sample. Those values will be stored in the SnpAnnotaionDataSet slot
#' and the ScanAnnotationDataSet slot and obtained via getter functions, e.g.
#' [getCountGenoRef()], [getCountAlleleRef()], and [getMAF()].
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param node Either of "raw", "filt", and "cor". See details.
#' @param ... Unused.
#' 
#' @details
#' #' Genotype call data can be obtained from the "genotype" node, the "filt.genotype"
#' node, or the "corrected.genotype" node of the GDS file with `node = "raw"`,
#' `node = "filt"`, or `node = "raw"`, respectively.
#' The [setCallFilter()] function generate filtered genotype call data in the
#' "filt.genotype" node which can be accessed as mentioned above. On the other hand, the
#' "corrected.genotype" node can be generated via the [clean()] function.
#'
#' @examples
#' \dontrun{
#' gdata <- loadGDS("/path/to/GDS.gds")
#' gdata <- countGenotype(gdata)
#' sample_missing_rate <- getCountGenoMissing(gdata, target = "scan", prop = TRUE)
#' marker_minor_allele_freq <- getMAF(gdata, target = "snp")
#' hist(gdata, stats = "missing")
#' }
#'
#' @export
#'
setGeneric("countGenotype", function(object,
                                     target = "both",
                                     node = "",
                                     ...)
  standardGeneric("countGenotype"))


#' Count reads per sample and per marker.
#'
#' This function calculates several summary statistics of read counts
#' per marker and per sample. Those values will be stored in the SnpAnnotaionDataSet slot
#' and the ScanAnnotationDataSet slot and obtained via getter functions, e.g.
#' [getCountReadRef()] and [getCountReadAlt()].
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param node Either of "raw" and "filt". See details.
#' @param ... Unused.
#' 
#' @details
#' Read count data can be obtained from the "annotation/format/AD/data" node or the
#' "annotation/format/AD/filt.data" node of the GDS file with `node = "raw"` or
#' `node = "filt"`, respectively. The [setCallFilter()] function generate filtered
#' read count data in the "annotation/format/AD/filt.data" node which can be accessed as
#' mentioned above.
#'
#' @examples
#' \dontrun{
#' gdata <- loadGDS("/path/to/GDS.gds")
#' gdata <- countRead(gdata)
#' read_depth_per_marker <- getCountRead(gdata, target = "snp")
#' reference_read_freq <- getCountReadRef(gdata, target = "snp", prop = TRUE)
#' hist(gdata, stats = "ad_ref")
#' }
#'
#' @export
#'
setGeneric("countRead", function(object,
                                 target = "both",
                                 node = "",
                                 ...)
  standardGeneric("countRead"))


#' Calculate mean, standard deviation, and quantile values of reads per sample and per marker.
#'
#' This function calculates several summary statistics of read counts
#' per marker and per sample. Those values will be stored in the SnpAnnotaionDataSet slot
#' and the ScanAnnotationDataSet slot and obtained via getter functions, e.g.
#' [getMeanReadRef()] and [getQtileReadAlt()].
#'
#' @param object A GbsrGenotypeData object.
#' @param target Either of "snp" and "scan".
#' @param q A numeric value [0-1] to indicate quantile to obtain.
#' @param ... Unused.
#' 
#' @details
#' Read count data can be obtained from the "annotation/format/AD/data" node or the
#' "annotation/format/AD/filt.data" node of the GDS file with `node = "raw"` or
#' `node = "filt"`, respectively. The [setCallFilter()] function generate filtered
#' read count data in the "annotation/format/AD/filt.data" node which can be accessed as
#' mentioned above.
#'
#' @importFrom stats sd quantile
#' 
#' @examples
#' \dontrun{
#' gdata <- loadGDS("/path/to/GDS.gds")
#' gdata <- calcReadStats(gdata, q = 0.5)
#' mean_reference_read_depth <- getMeanReadRef(gdata, target = "snp")
#' median_reference_read_depth <- getQtileReadAlt(gdata, target = "snp", q = 0.5)
#' hist(gdata, stats = "mean_ref")
#' hist(gdata, stats = "qtile_alt", q = 0.5)
#' }
#' 
#' @export
#'
setGeneric("calcReadStats", function(object,
                                     target = "both",
                                     q = NULL,
                                     ...)
  standardGeneric("calcReadStats"))


#' Set labels to samples which should be recognized as parents of the population to be subjected to error correction.
#'
#' Specify two or more samples in the dataset as parents of the population. Markers will be filtered out up on your specification.
#'
#' @param object A GbsrGenotypeData object.
#' @param parents A vector of strings with at least length two. The specified strings should match with the samples ID available via [getScanID()].
#' @param flip A logical value to indicate whether markers should be checked for "flip". See details.
#' @param mono A logical value whether to filter out markers which are not monomorphic in parents.
#' @param bi A logical value whether to filter out marekrs which are not biallelic between parents.
#' @param ... Unused.
#' 
#' @details
#' The `clean` function of `GBScleanR` uses read count information of samples and
#' their parents separately to estimate most probable genotype calls of them.
#' Therefore, you must specify proper samples as parents via this function.
#' If you would like to remove SNP markers which are not biallelic and/or
#' not monomorphic in each parent, set `mono = TRUE` and `bi = TRUE`.
#' `flip = TRUE` flips alleles of markers where the alleles expected as reference
#' allele are called as alternative allele. The alleles found in the parent specified as
#' the first element to the `parents` argument are supposed as reference alleles
#' of the markers. If the "expected" reference alleles are not actually called
#' as reference alleles but alternative alleles in the given data. setParents()
#' will automatically labels those markers "flipped".
#' The SnpAnnotatoinDataSet slot sores this information and accessible
#' via [getFlipped()] which gives you a logical vector
#' indicating which markers are labeled as flipped `TRUE` or not flipped `FALSE`.
#' [haveFlipped()] just tells you whether the SnpAnnotatoinDataSet slot has
#' the information of flipped markers or not.
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' \dontrun{
#' gds <- loadGDS("/path/to/GDS.gds")
#' parents <- grep("parent", getScanID(gds), value = TRUE)
#' gds <- setParents(gds, parents = parents, mono = TRUE, bi = TRUE, flip = TRUE)
#' gds <- clean(gds)
#' }
#'
#' @export
#'
setGeneric("setParents", function(object,
                                  parents,
                                  flip = FALSE,
                                  mono = FALSE,
                                  bi = TRUE,
                                  ...)
  standardGeneric("setParents"))


#' Get parental sample information
#'
#' This function returns scan IDs, member IDs and indexes of parental samples
#' set via [setParents()]. Scan IDs are IDs given by user or obtained from the
#' original VCF file. Member IDs are serial numbers assigned by [setParents()].
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = c("parent1", "parent2"))
#' getParents(gds)
#' }
#' 
#'
setGeneric("getParents", function(object,
                                  ...)
  standardGeneric("getParents"))


#' Swap the alleles recorded in a GDS file linked to the given GbsrGenotypeData object.
#'
#' The alleles of each marker are automatically obtained to match with those
#' recorded in an input VCF file when it was converted to a GDS file. This function swap
#' those alleles.
#'
#' @param object A GbsrGenotypeData object.
#' @param allele A vector or matrix of characters each of which is either of "A", "T", "G", and "C". The length and the number of rows should be same with the number of "valid" markers. See details.
#' @param ... Unused.
#' 
#' @details
#' The `allele` argument can take a vector or a two-column matrix of characters
#' indicating reference alleles or both alleles. If A vector was given, this
#' function check the current reference and alternative allele of each marker is
#' same with the specified allele for each marker in the vector.
#' If a marker showed that the current alternative allele matched with the allele
#' in the vector, the reference allele and the alternative allele will be swapped
#' each other. In the case of a matrix, the alleles specified in the first column are
#' supposed to be reference alleles while the second column is for alternative alleles.
#' This function compares both alleles between the current record and the specified in
#' the matrix. If a marker showed one of the alleles specified in the allele matrix
#' do not exist in the current allele of the marker, this marker will be labeled as "invalid"
#' marker. In the case of that both of the specified alleles exist but swapped
#' in the current record, the current alleles will be swapped to match with those specified in
#' the allele matrix.
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' \dontrun{
#' # In the case of that you have a reference genome data but it was not used 
#' # for the SNP call, or reference alleles in the genotype data do not match 
#' # with the alleles in the reference genome, e.g. TASSEL-GBS do.
#' gds <- loadGDS("/path/to/GDS.gds")
#' ref_genome <- Biostrings::readDNAStringSet("/path/to/genome.fasta")
#' chr_names <- getChromosome(object, name = TRUE)
#' snp_pos <- GenomicRanges::GRanges(seqnames = chr_names,
#'                                   ranges = IRanges::IRanges(start = getPosition(object),
#'                                                             width = 1))
#' ref_allele <- as.character(genome[snp_pos])
#' gds <- swapAlleles(gds, allele = ref_allele)
#' }
#'
#' @export
#'
setGeneric("swapAlleles", function(object, allele, ...)
  standardGeneric("swapAlleles"))


#' Filter out each genotype call meeting criteria
#'
#' Perform filtering of each genotype call, neither markers nor samples. Each genotype call
#' is supported by its read counts for the reference allele and the alternative allele of
#' a marker of a sample. `setCallFilter()` set missing to the genotype calls which are
#' not reliable enough and set zero to reference and alternative read counts of
#' the genotype calls.
#'
#' @param object A GbsrGenotypeData object.
#' @param dp_count A numeric vector with length two specifying lower and upper limit of total read counts (reference reads + alternative reads).
#' @param ref_count A numeric vector with length two specifying lower and upper limit of reference read counts.
#' @param alt_count A numeric vector with length two specifying lower and upper limit of alternative read counts.
#' @param norm_dp_count A numeric vector with length two specifying lower and upper limit of normalized total read counts (normalized reference reads + normalized alternaitve reads).
#' @param norm_ref_count A numeric vector with length two specifying lower and upper limit of normalized reference read counts.
#' @param norm_alt_count A numeric vector with length two specifying lower and upper limit of normalized alternative read counts
#' @param scan_dp_qtile A numeric vector with length two specifying lower and upper limit of quantile of total read counts in each scan (sample).
#'@param scan_ref_qtile A numeric vector with length two specifying lower and upper limit of quantile of reference read counts in each scan (sample).
#' @param scan_alt_qtile A numeric vector with length two specifying lower and upper limit of quantile of alternative read counts in each scan (sample).
#' @param snp_dp_qtile A numeric vector with length two specifying lower and upper limit of quantile of total read counts in each SNP marker
#'@param snp_ref_qtile A numeric vector with length two specifying lower and upper limit of quantile of reference read counts in each SNP marker.
#' @param snp_alt_qtile A numeric vector with length two specifying lower and upper limit of quantile of alternative read counts in each SNP marker.
#' @param omit_geno A vector of string with combinations of "ref", "het", and "alt" to remove specified genotype calls.
#' @param ... Unused.
#' 
#' @details
#' `norm_dp_count`, `norm_ref_count`, and `norm_alt_count` use normalized read counts which
#' are obtained by dividing each read count by the total read count of each sample.
#' `scan_dp_qtile`, `scan_ref_qtile`, and `scan_alt_qtile` work similarly but use quantile
#' values of read counts of each sample to decide the lower and upper limit of read counts.
#' This function generate two new nodes in the GDS file linked with the given GbsrGenotypeData
#' object. The new nodes "filt.data" in the AD node and "filt.genotype" contains read count
#' data and genotype data after filtering, respectively.
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' 
#' \dontrun{
#' # Filter out genotype calls supported by less than 5 reads.
#' gds <- setCallFilter(gds, dp_count = c(5, Inf))
#'
#' # Filter out genotype calls supported by reads less than 
#' # the 20 percentile of read counts per marker in each sample.
#' gds <- setCallFilter(gds, scan_dp_qtile = c(0.2, 1))
#'
#' # Filter out all reference homozygote genotype calls.
#' gds <- setCallFilter(gds, omit_geno = "ref")
#' }
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
                                     omit_geno = NULL,
                                     ...)
           standardGeneric("setCallFilter"))


#' Filter out scans (samples)
#'
#' Search samples which do not meet the criteria and label them as "invalid".
#'
#' @param object A GbsrGenotypeData object.
#' @param id A vector of strings match with scan ID which can be retrieve by `getScanID()`.
#' @param missing A numeric value [0-1] to specify the maximum missing genotype call rate per sample.
#' @param het A numeric value [0-1] to specify the maximum heterozygous genotype call rate per sample.
#' @param mac A integer value to specify the minimum minor allele count per sample.
#' @param maf A numeric value to specify the minimum minor allele frequency per sample.
#' @param ad_ref A numeric vector with length two specifying lower and upper limit of reference read counts per sample.
#' @param ad_alt A numeric vector with length two specifying lower and upper limit of alternative read counts per sample.
#' @param dp A numeric vector with length two specifying lower and upper limit of total read counts per sample.
#' @param mean_ref A numeric vector with length two specifying lower and upper limit of mean of reference read counts per sample.
#' @param mean_alt A numeric vector with length two specifying lower and upper limit of mean of alternative read counts per sample.
#' @param sd_ref A numeric value specifying the upper limit of standard deviation of reference read counts per sample.
#' @param sd_alt A numeric value specifying the upper limit of standard deviation of alternative read counts per sample.
#' @param ... Unused.
#' 
#' @details
#' For `mean_ref`, `mean_alt`, `sd_ref`, and `sd_alt`, this function calculate mean and
#' standard deviation of reads obtained at SNP markers of each sample. If a mean read counts
#' of a sample was smaller than the specified lower limit or larger than the upper limit,
#' this function labels the sample as "invalid". In the case of `sd_ref` and `sd_alt`,
#' standard deviations of read counts of each sample are checked and the samples having a
#' larger standard deviation will be labeled as "invalid". To check valid and invalid
#' samples, run [getValidScan()].
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' \dontrun{
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setScanFilter(gds, id = getScanID(gds)[1:10], missing = 0.2, dp = c(5, Inf))
#' }
#' 
#' @export
#'
setGeneric("setScanFilter", function(object,
                                     id,
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
#' @param object A GbsrGenotypeData object.
#' @param id A vector of strings match with scan ID which can be retrieve by `getScanID()`.
#' @param missing A numeric value [0-1] to specify the maximum missing genotype call rate per marker
#' @param het A numeric value [0-1] to specify the maximum heterozygous genotype call rate per marker
#' @param mac A integer value to specify the minimum minor allele count per marker
#' @param maf A numeric value to specify the minimum minor allele frequency per marker.
#' @param ad_ref A numeric vector with length two specifying lower and upper limit of reference read counts per marker.
#' @param ad_alt A numeric vector with length two specifying lower and upper limit of alternative read counts per marker.
#' @param dp A numeric vector with length two specifying lower and upper limit of total read counts per marker.
#' @param mean_ref A numeric vector with length two specifying lower and upper limit of mean of reference read counts per marker.
#' @param mean_alt A numeric vector with length two specifying lower and upper limit of mean of alternative read counts per marker.
#' @param sd_ref A numeric value specifying the upper limit of standard deviation of reference read counts per marker.
#' @param sd_alt A numeric value specifying the upper limit of standard deviation of alternative read counts per marker.
#' @param ... Unused.
#' 
#' @details
#' For `mean_ref`, `mean_alt`, `sd_ref`, and `sd_alt`, this function calculate mean and
#' standard deviation of reads obtained for samples at each SNP marker. If a mean read counts
#' of a marker was smaller than the specified lower limit or larger than the upper limit,
#' this function labels the marker as "invalid". In the case of `sd_ref` and `sd_alt`,
#' standard deviations of read counts of each marker are checked and the markers having a
#' larger standard deviation will be labeled as "invalid". To check valid and invalid
#' markers, run [getValidSnp()].
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' \dontrun{
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setSnpFilter(gds, id = getSnpID(gds)[1:1000], missing = 0.2, dp = c(5, Inf))
#' }
#'
#' @export
#'
setGeneric("setSnpFilter", function(object,
                                    id,
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
#' A VCF file usually has marker quality metrics in the INFO filed and those are stored in
#' a GDS file created via `GBScleanR`. This function filter out markers based on those marker
#' quality metrics.
#'
#' @param object A GbsrGenotypeData object.
#' @param mq A numeric value to specify minimum mapping quality (shown as MQ in the VCF format).
#' @param fs A numeric value to specify maximum Phred-scaled p-value (strand bias) (shown as FS in the VCF format).
#' @param qd A numeric value to specify minimum Variant Quality by Depth (shown as QD in the VCF format).
#' @param sor A numeric value to specify maximum Symmetric Odds Ratio (strand bias) (shown as SOR in the VCF format).
#' @param mqranksum A numeric values to specify the lower and upper limit of Alt vs. Ref read mapping qualities (shown as MQRankSum in the VCF format).
#' @param readposranksum A numeric values to specify the lower and upper limit of Alt vs. Ref read position bias (shown as ReadPosRankSum in the VCF format).
#' @param baseqranksum A numeric values to specify the lower and upper limit of Alt Vs. Ref base qualities (shown as BaseQRankSum in the VCF format).
#' @param ... Unused.
#' 
#' @details
#' Detailed explanation of each metric can be found in [GATK's web site](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' \dontrun{
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setInfoFilter(gds, mq = 40, qd = 20)
#' }
#' 
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


#' Reset the filter made by [setCallFiler()]
#'
#' Return genotype calls and read count data to the original data which are same with
#' those data before running [setCallFilter()].
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setGeneric("resetCallFilters", function(object, ...)
  standardGeneric("resetCallFilters"))


#' Reset the filter made by [setScanFiler()]
#'
#' Remove "invalid" labels put on samples and make all samples valid.
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setGeneric("resetScanFilters", function(object, ...)
  standardGeneric("resetScanFilters"))


#' Reset the filter made by [setSnpFiler()]
#'
#' Remove "invalid" labels put on markers and make all markers valid.
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setGeneric("resetSnpFilters", function(object, ...)
  standardGeneric("resetSnpFilters"))


#' Reset all filters made by [setScanFiler()], [setSnpFiler()], and [setCallFiler()].
#'
#' Return all data intact.
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setGeneric("resetFilters", function(object, ...)
  standardGeneric("resetFilters"))


#' Remove markers potentially having redundant information.
#'
#' Markers within the length of the sequenced reads (usually ~ 150 bp, up to your sequencer)
#' potentially have redundant information and those will cause unexpected errors
#' in error correction which assumes independency of markers each other.
#' This function only retains the first marker or the least missing rate marker
#' from the markers locating within the specified stretch.
#'
#' @param object A GbsrGenotypeData object.
#' @param range A integer value to indicate the stretch to search markers.
#' @param ... Unused.
#' 
#' @details
#' This function search valid markers from the first marker of each chromosome and
#' compare its physical position with a neighbor marker. If the distance between those
#' markers are equal or less then `range`, one of them which has a larger missing rate
#' will be removed (labeled as invalid marker). When the first marker was retained and
#' the second marker was removed as invalid marker, next the distance between the first marker
#' and the third marker will be checked and this cycle is repeated until reaching the
#' end of each chromosome. Run [getValidSnp()] to check the valid SNP markers.
#'
#' @return A GbsrGenotypeData object.
#'
#' @examples
#' \dontrun{
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- thinMarker(gds, range = 150)
#' }
#' 
#' @export
#'
setGeneric("thinMarker", function(object, range = 150, ...)
  standardGeneric("thinMarker"))


#' Create a GDS file with subset data of the current GDS file
#'
#' Create a new GDS file storing the subset data from the current GDS file linked to
#' the given GbsrGenotypeData object with keeping (or removing) information based on
#' valid markers and samples information.
#'
#' @param object A GbsrGenotypeData object.
#' @param out_fn A string to specify the path to an output GDS file.
#' @param snp_incl A logical vector having the same length with the total number of markers. The values obtained via [getValidSnp()] are used.
#' @param scan_incl A logical vector having the same length with the total number of scans (samples). The values obtained via [getValidScan()] are used.
#' @param incl_parents A logical value to specify whether parental samples should be included in a subset data or not.
#' @param ... Unused.
#' 
#' @details
#' A resultant subset data in a new GDS file includes subsets of each category of data, e.g.
#' genotype, SNP ID, scan ID, read counts, and quality metrics of SNP markers.
#'
#' @return A GbsrGenotypeData object linking to the new GDS file storing subset data.
#'
#' @examples
#' \dontrun{
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setSnpFilter(gds, missing = 0.2, het = c(0.1, 0.9), maf = 0.05)
#' new_gds <- subsetGDS(gds, out_fn = "/path/to/newGDS.gds")
#' }
#'
#' @export
#'
setGeneric("subsetGDS", function(object,
                                 out_fn,
                                 snp_incl,
                                 scan_incl,
                                 incl_parents = TRUE,
                                 ...)
  standardGeneric("subsetGDS"))


#' Set the filtered data to be used in GBScleanR's functions
#'
#' Set the "filt.genotype" node and the "filt.data" node as primary nodes for genotype
#' data and read count data. The data stored in the primary nodes are used in the
#' functions of GBScleanR.
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#' @details
#' A GbsrGenotypeData object storing information of the primary node of genotype data and
#' read count data. All of the functions implemented in `GBScleanR` check the primary nodes
#' and use data stored in those nodes. [setCallFilter()] create new nodes storing
#' filtered genotype calls and read counts in a GDS file and change the primary nodes to
#' "filt.genotype" and "filt.data" for genotype and read count data, respectively.
#' [SetRawGenotype()] set back the nodes to the original, those are "genotype" and "data" for
#' genotype and read count data, respectively. You can set the filtered data again by
#' `SetFiltGenotype()`.
#'
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setGeneric("setFiltGenotype", function(object, ...)
  standardGeneric("setFiltGenotype"))


#' Set the origina; data to be used in GBScleanR's functions
#'
#' Set the "genotype" node and the "data" node as primary nodes for genotype
#' data and read count data. The data stored in the primary nodes are used in the
#' functions of GBScleanR.
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#' @details
#' A GbsrGenotypeData object storing information of the primary node of genotype data and
#' read count data. All of the functions implemented in `GBScleanR` check the primary nodes
#' and use data stored in those nodes. [setCallFilter()] create new nodes storing
#' filtered genotype calls and read counts in a GDS file and change the primary nodes to
#' "filt.genotype" and "filt.data" for genotype and read count data, respectively.
#' `SetRawGenotype()` set back the nodes to the original, those are "genotype" and "data" for
#' genotype and read count data, respectively. You can set the filtered data again by
#' [SetFiltGenotype()].
#'
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setGeneric("setRawGenotype", function(object, ...)
  standardGeneric("setRawGenotype"))


#' Replace data stored in the GDS file
#'
#' Replace the original data with the filtered data or replace sample IDs in the GDS file linked
#' to an input GbsrGenotypeData object.
#'
#' @param object A GbsrGenotypeData object.
#' @param target A vector of combinations of "sample.id", "genotype", and "ad".
#' @param node Either one of "filt" and "cor" to replace raw genotype data with filtered genotype data or corrected genotype data, respectively.
#' @param ... Unused.
#' 
#' @details
#' If `target = "genotype"`, replace the data stored in the "genotype" node with the
#' data in the "filt.genotype" node for genotype call data. If `target = "ad`,
#' The data stored in the "data" node is also replaced with the data in the
#' "filt.data" node for read count data. `target = "sample.id"` makes this function to
#' replace data in the "sample.id" node with the sample IDs stored in the
#' ScanAnnotationDataSet slot of the input GbsrGenotypeData object.
#'
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
setGeneric("replaceGDSdata", function(object,
                                      target,
                                      node = "filt",
                                      ...)
  standardGeneric("replaceGDSdata"))


#' Genotype estimation using a hiden Morkov model
#'
#' Clean up genotype data by error correction based on genotype estimation using a hidden Markov model.
#'
#' @param object A GbsrGenotypeData object.
#' @param chr An integer vector of chromosome indices to be analyzed. All chromosomes will be analyzed if you left it default.
#' @param recomb_rate A numeric value to indicate the expected recombination frequency per chromosome per megabase pairs.
#' @param error_rate A numeric value of the expected sequence error rate.
#' @param call_threshold A numeric value of the probability threhosld to accept estimated genotype calls.
#' @param het_parent A logical value to indicate whether parental samples are outbred or inbred. If FALSE, this function assume all true genotype of markers in parents are homozygotes.
#' @param optim A logical value to specify whether to conduct parameter optimization for
#' error correction.
#' @param iter An integer value to specify the number of iterative parameter updates.
#' @param n_threads An integer value to specify the number of threads used for the calculation. The default is `n_threads = NULL` and automatically set half the number of available threads on the computer.
#' @param ... Unused.
#' 
#' @return A GbsrGenotypeData object.
#' 
#' @importFrom RcppParallel setThreadOptions defaultNumThreads
#' @importFrom expm expm.Higham08
#' @importFrom stats cor
#'
#' @examples
#' \dontrun{
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = getScanID(gds)[1:4])
#' gds <- initScheme(gds, crosstype = "pairing", mating = matrix(1:4, 2))
#' gds <- addScheme(gds, "pairing", mating = matrix(5:6, 2))
#' gds <- addScheme(gds, "selfing")
#' gds <- estGeno(gds, recomb_rate = 0.03, het_parent = FALSE, iter = 3)
#' }
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


#' Build a GbsrScheme object
#'
#' GBScleanR uses breeding scheme information to set the expected
#' number of cross overs in a chromosome which is a required parameter
#' for the genotype error correction with the Hidden Markov model
#' implemented in the "clean" function. This function build the object storing
#' type crosses performed at each generation of breeding and population sizes.
#'
#' @param object A GbsrGenotypeData object.
#' @param crosstype A string to indicate the type of cross conducted with a given generation.
#' @param mating An integer matrix to indicate mating combinations. The each element should match with IDs of parental samples which are 1 to N. see Details.
#' @param ... Unused.
#' 
#' @return A GbsrScheme object.
#'
#' @details
#' A GbsrScheme object stores information of a population size, mating combinations and
#' a type of cross applied to each generation of the breeding process
#' to generate the population which you are going to subject to the "clean" function.
#' The first generation should be parents of the population. It is supposed that
#' [setParents()] has been already executed and parents are labeled in the
#' GbsrGenotypeData object. The number of parents are automatically recognized.
#' The "crosstype" of the first generation can be "pairing" or "random" with
#' `pop_size = N`, where N is the number of parents.
#' You need to specify a matrix indicating combinations of `mating`, in which each column shows
#' a pair of parental samples. For example, if you have only two parents, the `mating` matirix
#' should be `mating = matrix(1:2, nrow = 1, ncol = 2)`. The indices used in the matrix
#' should match with the IDs labeled to parental samples by [setParents()].
#' The created GbsrScheme object is set in the `scheme` slot of the GbsrGenotypeData object.
#'
#' @export
#'
#' @seealso [addScheme()] and [showScheme()]
#'
#' @examples
#' \dontrun{
#' # Biparental F2 population.
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = c("parent1", "parent2"))
#' # setParents gave member ID 1 and 2 to parent1 and parent2, respectively.
#' gds <- initScheme(gds, crosstype = "pair", mating = matrix(1:2, nrow = 1, ncol = 2))
#' # Now the progenies of the cross above have member ID 3.
#' # If `crosstype = "selfing"` or `"sibling"`, you can omit a `mating` matrix.
#' gds <- addScheme(gds, crosstype = "self")
#'
#' # 8-way RILs with sibling mating.
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = paste("parent", 1:8, sep = ""))
#' # setParents set member ID 1 to 8 to parent1 to parent8, respectively.
#'
#' # If you made crosses of parent1 x parent2, parent3 x parent4,
#' # parent5 x parent6, and parent7 x parent8, run the following.
#' gds <- initScheme(gds, crosstype = "pair", mating = matrix(1:8, nrow = 4, ncol = 2))
#'
#' # Now the progenies of the crosses above have member ID 9, 10, 11, 
#' # and 12 for each combination of mating.You can check IDs with showScheme().
#' showScheme(gds)
#'
#' # Then, produce 4-way crosses.
#' gds <- addScheme(gds, crosstype = "pair", mating = matrix(9:12, nrow = 2, ncol = 2))
#' # 8-way crosses.
#' gds <- addScheme(gds, crosstype = "pair", mating = matrix(13:14, nrow = 1, ncol = 2))
#' # Inbreeding by 4 times selfing.
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#'
#' # Execute error correction
#' gds <- clean(gds)
#' }
#'
setGeneric("initScheme", function(object, crosstype, mating, ...)
  standardGeneric("initScheme"))


#' #' Build a GbsrScheme object
#'
#' GBScleanR uses breeding scheme information to set the expected
#' number of cross overs in a chromosome which is a required parameter
#' for the genotype error correction with the Hidden Markov model
#' implemented in the "clean" function. This function build the object storing
#' type crosses performed at each generation of breeding and population sizes.
#'
#' @param object A GbsrGenotypeData object.
#' @param crosstype A string to indicate the type of cross conducted with a given generation.
#' @param mating An integer matrix to indicate mating combinations. The each element should match with member IDs of the last generation.
#' @param pop_size An integer of the number of individuals in a given generation.
#' @param ... Unused.
#' 
#' @return A GbsrScheme object.
#'
#' @details
#' A scheme object is just a data.frame indicating a population size and
#' a type of cross applied to each generation of the breeding process
#' to generate the population which you are going to subject to the "clean" function.
#' The `crosstype` can take either of "selfing", "sibling", "pairing", and "random".
#' When you set `crosstype = "random"`, you need to specify `pop_size` to indicate how many
#' individuals were crossed in the random mating.
#' You also need to specify a matrix indicating combinations of `mating`, in which
#' each column shows a pair of member IDs indicating parental samples of the cross.
#' Member IDs are serial numbers starts from 1 and automatically assigned by
#' [initScheme()] and [addScheme()]. To check the member IDs, run [showScheme()].
#' Please see the examples section for more details of specifying a `mating` matrix.
#' The created GbsrScheme object is set in the `scheme` slot of the GbsrGenotypeData object.
#'
#' @export
#'
#' @seealso [addScheme()] and [showScheme()]
#'
#' @examples
#' \dontrun{
#' # Biparental F2 population.
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = c("parent1", "parent2"))
#' # setParents gave member ID 1 and 2 to parent1 and parent2, respectively.
#' gds <- initScheme(gds, crosstype = "pair", mating = matrix(1:2, nrow = 1, ncol = 2))
#' # Now the progenies of the cross above have member ID 3.
#' # If `crosstype = "selfing"` or `"sibling"`, you can omit a `mating` matrix.
#' gds <- addScheme(gds, crosstype = "self")
#'
#' # 8-way RILs with sibling mating.
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = paste("parent", 1:8, sep = ""))
#' # setParents set member ID 1 to 8 to parent1 to parent8, respectively.
#'
#' # If you made crosses of parent1 x parent2, parent3 x parent4, 
#' # parent5 x parent6, and parent7 x parent8, run the following.
#' gds <- initScheme(gds, crosstype = "pair", mating = matrix(1:8, nrow = 4, ncol = 2))
#'
#' # Now the progenies of the crosses above have member ID 9, 10, 11, 
#' # and 12 for each combination of mating.You can check IDs with showScheme().
#' showScheme(gds)
#'
#' # Then, produce 4-way crosses.
#' gds <- addScheme(gds, crosstype = "pair", mating = matrix(9:12, nrow = 2, ncol = 2))
#' # 8-way crosses.
#' gds <- addScheme(gds, crosstype = "pair", mating = matrix(13:14, nrow = 1, ncol = 2))
#' # Inbreeding by 4 times selfing.
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#' gds <- addScheme(gds, crosstype = "self")
#'
#' # Execute error correction
#' gds <- estGeno(gds)
#'
#'}
#'
setGeneric("addScheme", function(object, crosstype, mating, pop_size, ...)
  standardGeneric("addScheme"))


#' Show the information stored in a GbsrScheme object
#'
#' Print the information of each generation in a GbsrScheme object in the scheme
#' slot of a GbsrGenotypeData object.
#' A GbsrScheme object stores information of a population size, mating combinations and
#' a type of cross applied to each generation of the breeding process
#' to generate the population which you are going to subject to the "clean" function.
#'
#' @param object A GbsrGenotypeData object.
#' @param ... Unused.
#' 
#' @export
#'
#' @seealso [initScheme()] and [addScheme()]
#'
#' @examples
#' \dontrun{
#' gds <- loadGDS("/path/to/GDS.gds")
#' gds <- setParents(gds, parents = paste("parent", 1:8, sep = ""))
#' # setParents set member ID 1 to 8 to parent1 to parent8, respectively.
#' # If you made crosses of parent1 x parent2, parent3 x parent4, 
#' # parent5 x parent6, and parent7 x parent8, run the following.
#' gds <- initScheme(gds, crosstype = "pair", mating = matrix(1:8, nrow = 4, ncol = 2))
#'
#' # Now the progenies of the crosses above have member ID 9, 10, 11,
#' # and 12 for each combination of mating. You can check IDs with showScheme().
#' showScheme(gds)
#' }
#'
setGeneric("showScheme", function(object, ...)
  standardGeneric("showScheme"))
