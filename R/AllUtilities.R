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
#' @param gt the ID for genotypic data in the FORMAT column; "GT" by default,
#'  VCFv4.0.
#' @param info.import characters, the variable name(s) in the INFO field for
#' import; or NULL for all variables. If you specify `character(0)`, nothing
#' will be retrieved from the INFO filed.
#' @param fmt.import characters, the variable name(s) in the FORMAT field for
#' import; or NULL for all variables. If you specify `character(0)`, nothing
#' will be retrieved from the FORMAT filed, except for GT.
#' @param force A logical value to overwrite a GDS file even if
#' the file specified in "out_fn" exists.
#' @param verbose if TRUE, show information.
#'
#' @details
#' gbsrVCF2GDS converts a VCF file to a GDS file.
#' The data structure of the GDS file created via this functions is same with
#' those created by `seqVCF2GDS()` of `SeqArray`.
#'
#' @return The output GDS file path.
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a `GbsrGenotypeData` object.
#' gds <- loadGDS(gds_fn)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
#' @export
#'
#' @importFrom SeqArray seqVCF2GDS seqOptimize seqStorageOption seqExport seqNumAllele
#'
gbsrVCF2GDS <- function(vcf_fn,
                        out_fn,
                        gt = "GT",
                        info.import = NULL,
                        fmt.import = NULL,
                        force = FALSE,
                        verbose = TRUE) {
    if (file.exists(out_fn)) {
        message(out_fn, ' exists!')
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
    stropt <- seqStorageOption(mode = c ('annotation/format/AD' = "int64"))
    out_fn <- seqVCF2GDS(vcf_fn, out_fn, ignore.chr.prefix = "",
                         storage.option = stropt,
                         info.import = info.import,
                         fmt.import = fmt.import,
                         genotype.var.name = gt,
                         verbose = verbose)
    gds <- seqOpen(out_fn)
    n_allele <- seqNumAllele(gds)
    if(any(n_allele != 2)){
        message("Some markers are not bi-allelic.",
                "\nGBScleanR filters out those markers.")
        message("Filtering non-biallelic markers.")
        temp_fn <- tempfile("biallelic", tempdir(), ".gds")
        seqSetFilter(gds, variant.sel = n_allele == 2)
        seqExport(gds, temp_fn)
        seqClose(gds)
        file.copy(temp_fn, out_fn, overwrite = TRUE)
    } else {
        seqClose(gds)
    }
    return(out_fn)
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
#' @param load_filter A logical whether to load the filtering information made
#'  via [setSamFilter()] and [setMarFilter()] and saved in the GDS file via
#'  [closeGDS()] with `save_filter = TRUE`.
#' @param ploidy Set the ploidy of the population.
#' @param verbose if TRUE, show information.
#'
#' @return A `GbsrGenotypeData` object.
#'
#' @export
#'
#' @importFrom methods new
#' @importFrom SeqArray seqOpen seqSummary seqOptimize
#' @importFrom gdsfmt exist.gdsn delete.gdsn objdesp.gdsn setdim.gdsn add.gdsn
#'
#' @details
#' The first time to load a newly produced GDS file will take long time due to
#' data reformatting for quick access. The `GbsrGenotypeData` object returned
#' from `loadGDS()` can be also handled as `SeqVarGDSClass` of the [SeqArray]
#' package.
#'
#' @examples
#' # Create a GDS file from a sample VCF file.
#' vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
#' gds_fn <- tempfile("sample", fileext = ".gds")
#' gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)
#'
#' # Load data in the GDS file and instantiate a `GbsrGenotypeData` object.
#' gds <- loadGDS(gds_fn)
#'
#' # Reload data from the GDS file.
#' gds <- loadGDS(gds)
#'
#' # Close the connection to the GDS file.
#' closeGDS(gds)
#'
loadGDS <- function(x, load_filter = FALSE, ploidy = 2, verbose = TRUE) {
    if(verbose){ message('Loading GDS file.') }
    if(inherits(x = x, what = "GbsrGenotypeData")){
        if(isOpenGDS(object = x)){
            closeGDS(object = x, save_filter = FALSE)
        }
        # Leave the following connection open to build the GdsGenotypeReader object
        # with `readonly=FALSE` mode.
        gds <- seqOpen(gds.fn = x$filename, readonly = FALSE)

    } else if(is.character(x)){
        # Leave the following connection open to build the GdsGenotypeReader object
        # with `readonly=FALSE` mode.
        gds <- seqOpen(gds.fn = x, readonly = FALSE)

    } else {
        stop("x must be a file path or a GbsrGenotypeData object",
             call. = FALSE)
    }

    for(node in c("FGT", "CGT", "HAP")){
        gdsn <- paste0("annotation/format/", node)
        if(exist.gdsn(gds, gdsn)){
            message("Reformatting ", node)
            i_gdsn <- index.gdsn(gds, paste0(gdsn, "/data"))
            gdsn_dim <- objdesp.gdsn(i_gdsn)
            if(length(gdsn_dim$dim) != 3){
                if(node == "FGT"){storage <- "bit2"; na_val <- 3}
                if(node == "CGT"){storage <- "bit2"; na_val <- 3}
                if(node == "HAP"){storage <- "bit6"; na_val <- 0}

                tmp_gdsn <- add.gdsn(index.gdsn(gds, gdsn),
                                     "tmp", NULL, storage, NULL,
                                     replace = TRUE)
                apply.gdsn(i_gdsn, margin = 1,
                           function(x){
                               x <- unlist(strsplit(x, "\\||/"))
                               x[x == "."] <- na_val
                               x <- as.integer(x)
                               return(x)
                           },
                           as.is = "gdsnode",
                           target.node = tmp_gdsn)
                setdim.gdsn(tmp_gdsn, c(2, gdsn_dim$dim[2], gdsn_dim$dim[1]))
                i_gdsn <- add.gdsn(index.gdsn(gds, gdsn),
                                   "data", NULL, storage, NULL,
                                   replace = TRUE)
                apply.gdsn(tmp_gdsn, margin = 2, c, as.is = "gdsnode",
                           target.node = i_gdsn)
                setdim.gdsn(i_gdsn, c(2, gdsn_dim$dim[1], gdsn_dim$dim[2]))
                delete.gdsn(tmp_gdsn)
            }
        }
    }

    opt <- exist.gdsn(node = gds, path = "genotype/~data") &
        exist.gdsn(node = gds, path = "annotation/format/AD/~data")
    if(!opt){
        seqClose(object = gds)
        seqOptimize(gdsfn = gds$filename, target = "by.sample",
                    format.var = c("AD", "FAD", "HAP", "CGT", "FGT"),
                    verbose = verbose)
        gds <- seqOpen(gds.fn = gds$filename, readonly = FALSE)
    }

    d <- seqSummary(gdsfile = gds, varname = "genotype", verbose = FALSE)
    marker <- data.frame(valid = rep(TRUE, d$dim[3]))
    sample <- data.frame(valid = rep(TRUE, d$dim[2]))
    attributes(sample) <- c(attributes(sample), list(ploidy = ploidy))
    gds <- new(Class = "GbsrGenotypeData", gds, sample = sample, marker = marker)

    if(load_filter){
        if(exist.gdsn(gds, "sample.annotation/GFL")){
            validSam(gds) <- read.gdsn(index.gdsn(gds, "sample.annotation/GFL"))
            validMar(gds) <- read.gdsn(index.gdsn(gds, "annotation/info/GFL"))

        } else {
            warnings("No filtering information stored in the GDS.")
        }
    }
    return(gds)
}

