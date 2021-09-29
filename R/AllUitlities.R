print.GbsrGenotypeData <- function(x){
  message('Data in GDS file...')
  message('GDS file name')
  message(x@data@handler$filename)
  if(isOpenGDS(x)){
    print(x@data)
  } else {
    message('Connection to GDS file has been closed.')
  }
  message('-----------------------------------------------------')
  message('SnpAnnotationDataSet')
  print(x@snpAnnot)
  message('-----------------------------------------------------')
  message('ScanAnnotationDataSet')
  print(x@scanAnnot)
}

#' Convert a VCF file to a GDS file
#'
#' This function converts a variant call data in the VCF format.
#' The current implementation only accepts biallelic single nucleotide polymorphisms.
#' Please filter out variants which are insertions and deletions or multiallelic.
#' You may use "bcftools" or "vcftools" for filtering.
#'
#' @param vcf_fn A string to indicate path to an input VCF file.
#' @param out_fn A string to indicate path to an output GDS file.
#' @param force A logical value to overwrite a GDS file even if the file specified in "out_fn" exists.
#'
#' @details
#' gbsrVCF2GDS converts a VCF file to a GDS file.
#' The data structure of the GDS file created via this functions is same with
#' those created by snpgdsVCF2GDS of
#' [SNPRelate](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html).
#'  Unlike the GDS files created by SNPRelate GBScleanR's GDS files: \cr
#' \itemize{
#' \item{"Include annotation information of variants recorded in the INFO filed of the input VCF."}
#' \item{"Include allelic read count data (indicated as AD) recoreded in the FORMAT filed of the input VCF."}
#' }
#' @return The output GDS file path.
#'
#' @export
#'
#' @importFrom SeqArray seqVCF2GDS seqGDS2SNP
#' @import gdsfmt
#'
gbsrVCF2GDS <- function(vcf_fn,
                       out_fn,
                       force=FALSE){

  if(file.exists(out_fn)){
    message(paste0(out_fn, ' exists!'))
    while(TRUE){
      if(force){
        file.remove(out_fn)
        break()
      }
      user <- readline(prompt = "Are you sure to overwirte?(y/n)")
      if(user == "y"){
        file.remove(out_fn)
        break()
      } else if(user == "n"){
        message('Stopped by user.')
        return(out_fn)
      }
    }
  }

  # Create GDS file formatted in the SeqArray style.
  out <- SeqArray::seqVCF2GDS(vcf.fn = vcf_fn,
                              out.fn = sub("\\.gds", ".tmp.gds", out_fn))
  out <- .convertGDS(gds_fn = out)
  return(out)
}

# This is the function internally used to convert the format of a GDS file from that in SeqArray's GDS to SNPRelate's GDS to make the GDS valid for "GWASTools".

.convertGDS <- function(gds_fn){
  message('Reformatting genotype data.')

  new_gds_fn <- sub("\\.tmp\\.gds", ".gds", gds_fn)
  SeqArray::seqGDS2SNP(gdsfile = gds_fn, out.gdsfn = new_gds_fn)
  old_gds <- gdsfmt::openfn.gds(gds_fn, readonly = FALSE)
  new_gds <- gdsfmt::openfn.gds(new_gds_fn, readonly = FALSE)
  on.exit(expr = {
    gdsfmt::closefn.gds(old_gds)
    gdsfmt::closefn.gds(new_gds)
    file.remove(gds_fn)
    gdsfmt::cleanup.gds(filename=new_gds_fn)
    })
  new_ano <- gdsfmt::addfolder.gdsn(node=new_gds, name="annotation")
  gdsfmt::copyto.gdsn(node=new_ano,
                      source=gdsfmt::index.gdsn(node=old_gds, path="annotation/info"))
  ls_format <- try(gdsfmt::ls.gdsn(node=old_gds,
                               recursive=TRUE,
                               include.dirs=TRUE), silent=TRUE)

  if(!class(ls_format) %in% "try-error"){
    ls_format <- grep("format/.+/data", ls_format, value=TRUE)
    if(length(ls_format) != 0){
      new_format <- gdsfmt::addfolder.gdsn(node=new_ano, name="format", replace=TRUE)
      for(i in ls_format){
        label <- gsub(".+format/|/data", "", i)
        srcnode <- gdsfmt::index.gdsn(node=old_gds, path=i)
        desp <- gdsfmt::objdesp.gdsn(node=srcnode)

        new_sub <- gdsfmt::addfolder.gdsn(node=new_format, name=label, replace=TRUE)
        newdata <- gdsfmt::add.gdsn(node=new_sub,
                                    name="data",
                                    replace=TRUE,
                                    storage=desp$storage,
                                    compress=desp$compress)
        gdsfmt::assign.gdsn(node=newdata, src.node=srcnode)


        old_attr <- gdsfmt::get.attr.gdsn(gdsfmt::index.gdsn(old_gds,
                                                             gsub("/data", "", i)))
        if(!is.null(old_attr)){
          for(i in 1:length(old_attr))
            gdsfmt::put.attr.gdsn(new_sub, name = names(old_attr)[i], val = old_attr[[i]])
        }
        old_attr <- gdsfmt::get.attr.gdsn(srcnode)
        if(!is.null(old_attr)){
          for(i in 1:length(old_attr))
            gdsfmt::put.attr.gdsn(newdata, name = names(old_attr)[i], val = old_attr[[i]])
        }
      }
    }
  }
  chr_node <- gdsfmt::index.gdsn(node=new_gds, path="snp.chromosome")
  new_chr_node <- gdsfmt::add.gdsn(node=new_gds,
                                   name="snp.chromosome.name",
                                   storage = "string",
                                   compress = "LZMA_ra",
                                   replace=TRUE)
  gdsfmt::moveto.gdsn(node=chr_node,
                      loc.node=new_chr_node,
                      relpos="replace+rename")
  chr <- gdsfmt::read.gdsn(node=new_chr_node)
  chr <- factor(chr, levels = unique(chr))
  gdsfmt::add.gdsn(node=new_gds,
                   name="snp.chromosome",
                   val=as.integer(chr),
                   storage="int8",
                   compress="LZMA_ra",
                   replace=TRUE)
  gdsfmt::readmode.gdsn(gdsfmt::index.gdsn(new_gds, "snp.chromosome"))
  return(new_gds_fn)
}

# Internal function to build data for slots of the GenotypData class in GWASTools.

.buildSnpAnnot <- function(gds){
  snpID <- GWASTools::getSnpID(gds)
  chromosome.name <- factor(GWASTools::getChromosome(gds))
  chromosome <- as.integer(factor(GWASTools::getChromosome(gds)))
  position <- GWASTools::getPosition(gds)
  alleleA <- GWASTools::getAlleleA(gds)
  alleleB <- GWASTools::getAlleleB(gds)
  validMarker <- rep(TRUE, length(snpID))
  pdata <- data.frame(snpID, chromosome, chromosome.name, position,
             alleleA, alleleB, validMarker,
             ploidy=2L,
             stringsAsFactors=FALSE)

  ls_gdsn <- gdsfmt::ls.gdsn(node=gds@handler, recursive=TRUE, include.dirs=TRUE)
  if("info/QD" %in% ls_gdsn){
    pdata <- cbind(pdata,
                   QD = gdsfmt::read.gdsn(node=gdsfmt::index.gdsn(node=gds@handler,
                                                                  path="info/QD")))
  }
  if("info/FS" %in% ls_gdsn){
    pdata <- cbind(pdata,
                   FS = gdsfmt::read.gdsn(node=gdsfmt::index.gdsn(node=gds@handler,
                                                                  path="info/FS")))
  }
  if("info/SOR" %in% ls_gdsn){
    pdata <- cbind(pdata,
                   SOR = gdsfmt::read.gdsn(node=gdsfmt::index.gdsn(node=gds@handler,
                                                                   path="info/SOR")))
  }
  if("info/MQ" %in% ls_gdsn){
    pdata <- cbind(pdata,
                   MQ = gdsfmt::read.gdsn(node=gdsfmt::index.gdsn(node=gds@handler,
                                                                  path="info/MQ")))
  }
  if("info/MQRankSum" %in% ls_gdsn){
    pdata <- cbind(pdata,
                   MQRankSum = gdsfmt::read.gdsn(node=gdsfmt::index.gdsn(node=gds@handler,
                                                                         path="info/MQRankSum")))
  }
  if("info/ReadPosRankSum" %in% ls_gdsn){
    pdata <- cbind(pdata,
                   ReadPosRankSum = gdsfmt::read.gdsn(node=gdsfmt::index.gdsn(node=gds@handler,
                                                                              path="info/ReadPosRankSum")))
  }

  return(
    GWASTools::SnpAnnotationDataFrame(pdata)
  )
}

# Internal function to build data for slots of the GenotypData class in GWASTools.

.buildScanAnnot <- function(gds){
  scanID <- GWASTools::getScanID(gds)
  validScan <- rep(TRUE, length(scanID))

  return(
    GWASTools::ScanAnnotationDataFrame(data.frame(scanID, validScan, stringsAsFactors=FALSE))
  )
}

#' Load a GDS file and construct a GbsrGenotypeData object.
#'
#' Load data stored in an input GDS file to R environment and
#' create a GbsrGenotypeData instance.
#' GBScleanR handles only one class GbsrGenotypeData and
#' conducts all data manipulation via class methods for it.
#'
#' @param gds_fn A string of the path to an input GDS file.
#' @param non_autosomes A list with elements X, Y, XY, and/or M. See details.
#' @param genotypeData A GbsrGenotypeData object to reload.
#'
#' @return A GbsrGenotypeData object.
#'
#' @export
#'
#' @importClassesFrom GWASTools GenotypeData ScanAnnotationDataFrame SnpAnnotationDataFrame GdsGenotypeReader 
#' @importMethodsFrom GWASTools getSnpVariable getSnpAnnotation getScanAnnotation
#' @importFrom methods new
#'
#' @examples
#' \dontrun{
#' gds_fn <- gbsrVCF2GDS(vcf_fn = "/path/to/vcf.vcf", out_fn = "/path/to/gds.gds/", force = TRUE)
#' gdata <- loadGDS(gds_fn = gds_fn, non_autosomes = NULL)
#'
#' # If you would like to reload a GDS file.
#' gdata <- loadGDS(genotypeData = gdata)
#' }
#'
loadGDS <- function(gds_fn, non_autosomes=NULL, genotypeData){
  if(missing(genotypeData)){
    message('Loading GDS file.')

    gds <- GWASTools::GdsGenotypeReader(filename=gds_fn, genotypeDim="scan,snp", allow.fork = TRUE)

    message('Building SnpAnnotationDataFrame and ScanAnnotationDataFrame.')

    out <- new("GbsrGenotypeData", data=gds,
               snpAnnot=.buildSnpAnnot(gds),
               scanAnnot=.buildScanAnnot(gds))

    if(is.null(non_autosomes)){
      message('Set no non-autosomes.')
      n_chr <- length(unique(GWASTools::getSnpVariable(object=out, varname="chromosome")))
      out@data@XchromCode <- n_chr + 1L
      out@data@XYchromCode <- n_chr + 2L
      out@data@YchromCode <- n_chr + 3L
      out@data@MchromCode <- n_chr + 4L
    } else {
      message('Set the indices of non-autosomes.')
      out@XchromCode <- non_autosomes$X
      out@XYchromCode <-non_autosomes$XY
      out@YchromCode <- non_autosomes$Y
      out@MchromCode <- non_autosomes$M
    }

  } else {
    message('Reload GDS file.')
    gds_fn <- genotypeData@data@filename
    if(isOpenGDS(genotypeData)){
      closeGDS(genotypeData)
    }
    gds <- GWASTools::GdsGenotypeReader(filename=gds_fn, genotypeDim="scan,snp")

    out <- new("GbsrGenotypeData", data=gds,
               snpAnnot=getSnpAnnotation(genotypeData),
               scanAnnot=getScanAnnotation(genotypeData))

    n_chr <- length(unique(GWASTools::getSnpVariable(object=out, varname="chromosome")))
    out@data@XchromCode <- genotypeData@data@XchromCode
    out@data@XYchromCode <- genotypeData@data@XYchromCode
    out@data@YchromCode <- genotypeData@data@YchromCode
    out@data@MchromCode <- genotypeData@data@MchromCode
  }

  gdsfmt::closefn.gds(out@data@handler)
  out@data@handler <- gdsfmt::openfn.gds(filename=out@data@filename,
                                         readonly=FALSE,
                                         allow.duplicate=TRUE)
  return(out)
}
