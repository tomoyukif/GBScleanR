print.GbsrGenotypeData <- function(x) {
  show(x)
}

#' Convert a VCF file to a GDS file
#'
#' This function converts a variant call data in the VCF format.
#' The current implementation only accepts biallelic
#' single nucleotide polymorphisms.
#' Please filter out variants which are insertions and
#' deletions or multiallelic.
#' You may use "bcftools" or "vcftools" for filtering.
#'
#' @param vcf_fn A string to indicate path to an input VCF file.
#' @param out_fn A string to indicate path to an output GDS file.
#' @param force A logical value to overwrite a GDS file even if
#' the file specified in "out_fn" exists.
#' @param verbose if TRUE, show information.
#'
#' @details
#' gbsrVCF2GDS converts a VCF file to a GDS file.
#' The data structure of the GDS file created via this functions is same with
#' those created by `snpgdsVCF2GDS` of `SNPRelate`.
#' Unlike the GDS files created by `SNPRelate` `GBScleanR`'s GDS files: \cr
#' \itemize{
#' \item{"Include annotation information of variants
#' recorded in the `INFO` filed of the input VCF."}
#' \item{"Include allelic read count data (indicated as AD)
#' recoreded in the `FORMAT` filed of the input VCF."}
#' }
#' @return The output GDS file path.
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#' 
#' # Load data in the GDS file and instantiate a `GbsrGenotypeData` object.
#' gdata <- loadGDS(gds_fn)
#' 
#' # Close the connection to the GDS file.
#' closeGDS(gdata)
#'
#' @export
#'
#' @importFrom SeqArray seqVCF2GDS seqGDS2SNP seqStorageOption
#' @import gdsfmt
#'
gbsrVCF2GDS <- function(vcf_fn,
                        out_fn,
                        force = FALSE,
                        verbose = TRUE) {
  if (file.exists(out_fn)) {
    msg <- paste0(out_fn, ' exists!')
    message(msg)
    while (TRUE) {
      if (force) {
        file.remove(out_fn)
        break()
      }
      user <- readline("Are you sure to overwirte?(y/n)")
      if (user == "y") {
        file.remove(out_fn)
        break()
      } else if (user == "n") {
        message('Stopped by user.')
        return(out_fn)
      }
    }
  }
  
  # Create GDS file formatted in the SeqArray style.
  str_opt <- seqStorageOption(mode=c('annotation/format/AD'="uint32"))
  tmp_fn <- tempfile(fileext="gds")
  tmp_fn <- seqVCF2GDS(vcf_fn, tmp_fn,
                       fmt.import="AD", storage.option=str_opt, verbose=verbose)
  on.exit({unlink(tmp_fn)})
  .convertGDS(out_fn, tmp_fn, verbose)
  return(out_fn)
}

# This is the function internally used to convert the format of
# a GDS file from that in SeqArray's GDS to SNPRelate's
# GDS to make the GDS valid for "GWASTools".
.convertGDS <- function(out_fn, tmp_fn, verbose) {
  if(verbose){
    message('Reformatting genotype data.')
  }
  
  out_fn <- seqGDS2SNP(tmp_fn, out_fn, verbose=verbose)
  tmp_gds <- openfn.gds(tmp_fn, FALSE)
  on.exit({closefn.gds(tmp_gds)})
  out_gds <- openfn.gds(out_fn, FALSE)
  on.exit({closefn.gds(out_gds)}, TRUE)
  .gds_decomp(tmp_gds)
  .gds_decomp(out_gds)
  
  snp_ano <- addfolder.gdsn(out_gds, "annotation")
  addfolder.gdsn(snp_ano, "info")
  addfolder.gdsn(snp_ano, "format")
  
  .insertAnnot(tmp_gds, out_gds)
  
  chr_node <-.getNodeIndex(out_gds, "snp.chromosome")
  rename.gdsn(chr_node, "snp.chromosome.name")
  
  chr <- read.gdsn(chr_node)
  chr <- factor(chr, levels = unique(chr))
  add.gdsn(out_gds, "snp.chromosome", as.integer(chr), "int8", replace=TRUE)
}

## Get the list of the paths including `path` in its path.
.getNodeList <- function(gds_class, path){
  gdsn_i <- ls.gdsn(gds_class, recursive=TRUE, include.dirs=FALSE)
  out <- gdsn_i[grep(path, gdsn_i)]
  names(out) <- sub(path, "", out)
  return(out)
}

.insertAnnot <- function(tmp_gds, out_gds){
  info_node_list <- .getNodeList(tmp_gds, "annotation/info")
  
  # Copy INFO data to new GDS file.
  out_info <- index.gdsn(out_gds, "annotation/info")
  for (gdsn_i in info_node_list){
    i_node <- index.gdsn(tmp_gds, gdsn_i)
    copyto.gdsn(out_info, i_node)
  }
  
  # Copy FORMAT data to new GDS file.
  format_node_list <- .getNodeList(tmp_gds, "annotation/format")
  format_node_list <- grep("/data$", format_node_list, value=TRUE)
  out_format <- index.gdsn(out_gds, "annotation/format")
  for (gdsn_i in format_node_list){
    fld_name <- sub("/data$", "", gdsn_i)
    fld <- addfolder.gdsn(out_format, sub(".*/", "", fld_name))
    node_data <- index.gdsn(tmp_gds, gdsn_i)
    copyto.gdsn(fld, node_data)
    att <- get.attr.gdsn(index.gdsn(tmp_gds, fld_name))
    for(i in seq_along(att)){
      put.attr.gdsn(fld, names(att)[i], att[[i]])
    }
  }
}

# Internal function to build data for slots of
# the GenotypData class in GWASTools.
#' @importFrom GWASTools getVariable
.buildSnpAnnot <- function(gds) {
  snpID <- getVariable(gds, "snp.id")
  chromosome.name <- getVariable(gds, "snp.chromosome.name")
  chromosome <- getVariable(gds, "snp.chromosome")
  position <- getVariable(gds, "snp.position")
  allele <- getVariable(gds, "snp.allele")
  alleleA <- sub("/.*", "", allele)
  alleleB <- sub(".*/", "", allele)
  validMarker <- rep(TRUE, length(snpID))
  data <- data.frame(snpID, chromosome, chromosome.name, position, 
                     alleleA, alleleB, validMarker,
                     ploidy=2L, stringsAsFactors=FALSE)
  
  return(SnpAnnotationDataFrame(data))
}

# Internal function to build data for slots of
# the GenotypData class in GWASTools.
#' @importFrom GWASTools getVariable
.buildScanAnnot <- function(gds) {
  scanID <- getVariable(gds, "sample.id")
  validScan <- rep(TRUE, length(scanID))
  data <- data.frame(scanID, validScan, stringsAsFactors=FALSE)
  return(ScanAnnotationDataFrame(data))
}

#' Load a GDS file and construct a `GbsrGenotypeData` object.
#'
#' Load data stored in an input GDS file to R environment and
#' create a `GbsrGenotypeData` instance.
#' GBScleanR handles only one class `GbsrGenotypeData` and
#' conducts all data manipulation via class methods for it.
#'
#' @param x A string of the path to an input GDS file or 
#' a `GbsrGenotypeData` object to reload.
#' @param verbose if TRUE, show information.
#'
#' @return A `GbsrGenotypeData` object.
#'
#' @export
#'
#' @importClassesFrom GWASTools GenotypeData ScanAnnotationDataFrame
#' @importClassesFrom GWASTools SnpAnnotationDataFrame GdsGenotypeReader
#' @importFrom GWASTools GdsGenotypeReader GenotypeData
#' @importFrom methods new
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#' 
#' # Load data in the GDS file and instantiate a `GbsrGenotypeData` object.
#' gdata <- loadGDS(gds_fn)
#'
#' # Reload data from the GDS file.
#' gdata <- loadGDS(gdata)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gdata)
#' 
loadGDS <- function(x, 
                    verbose = TRUE) {
  if(verbose){ message('Loading GDS file.') }
  if(inherits(x, "GbsrGenotypeData")){
    if(isOpenGDS(x)){
        closeGDS(x)
    }
    x <- openfn.gds(.getGDSFileName(x), FALSE)
    
  } else if(is.character(x)){
    x <- openfn.gds(x, FALSE)
    
  } else {
    stop("x must be a file path or a GbsrGenotypeData object",
         call. = FALSE)
  }
  
  gds <- GdsGenotypeReader(x, "scan,snp")
  snp_ann <- .buildSnpAnnot(gds)
  scan_ann <- .buildScanAnnot(gds)
  gds <- GenotypeData(gds, snp_ann, scan_ann)
  gds <- new("GbsrGenotypeData", gds)
  .gds_comp(gds)
  return(gds)
}