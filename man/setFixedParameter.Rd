% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/Methods-GbsrGenotypeData.R
\name{setFixedParameter}
\alias{setFixedParameter}
\alias{setFixedParameter,GbsrGenotypeData-method}
\title{Set fixed allele read biases and mismapping rate}
\usage{
setFixedParameter(object, bias = NULL, mismap = NULL, parent_geno = FALSE, ...)

\S4method{setFixedParameter}{GbsrGenotypeData}(object, bias, mismap, parent_geno)
}
\arguments{
\item{object}{A \link{GbsrGenotypeData} object.}

\item{bias}{A numeric vector of fixed allele read biases to be assigned to
valid markers. The length of \code{bias} vector should match the number of valid
markers.}

\item{mismap}{A numeric matrix of fixed reference and alternative read
mismapping rates to be assigned to valid markers. The number of rows of the
given matrix should match the number of valid markers and should have two
columns that are for reference and alternative read mismapping rates,
respectively.}

\item{parent_geno}{A logical value indicating whether to use fixed parental
genotypes in the genoype estimation by \code{\link[=estGeno]{estGeno()}}. This mode requires the
estimated genotypes for parental samples that were estimated by \code{\link[=estGeno]{estGeno()}}
and stored in the GDS file linked to the input \link{GbsrGenotypeData} object.}

\item{...}{Unused.}
}
\value{
A \link{GbsrGenotypeData} object after adding dominant marker information
}
\description{
Set fixed allele read biases and mismapping rates of markers
}
\details{
If you have already executed genotype estimation and want to reuse the
marker-wise allele read biases and mismapping rates estimated in the
completed run of estGeno(), you can use them in the next genotype estimation
run.
For example, if you want to estimate genotypes with different argument
settings of setGeno(), it is worth to set fixed parameters and run
estGeno() with setting \code{optim = FALSE} to skip time-consuming iterative
parameter optimization steps but use the estimated parameters from the first
run to incorporate the marker-wise error parameters.
Since the bias set by \code{\link[=setFixedParameter]{setFixedParameter()}} function is the reference allele read
bias, the values 0 and 1 mean that the marker only gives alternative
and reference allele reads, respectively.
The values in the \code{bias} vector are assigned to the valid markers.
Similarly, the values in the \code{mismap} matrix are
assigned to the valid markers in the order they appear in the rows.
}
\examples{
# Create a GDS file from a sample VCF file.
vcf_fn <- system.file("extdata", "sample.vcf", package = "GBScleanR")
gds_fn <- tempfile("sample", fileext = ".gds")
gbsrVCF2GDS(vcf_fn = vcf_fn, out_fn = gds_fn, force = TRUE)

# Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
gds <- loadGDS(gds_fn)

# Not run.
# Run estGeno() and reuse the estimated parameters in the second run.
# gds <- makeScheme(gds, generation = 2, crosstype = "self")
# gds <- estGeno(gds)
# fixed_param <- getFixedParameter(gds)
# gds <- setFixedParameter(gds,
#                          bias = fixed_param$bias,
#                          mismap = fixed_param$mismap)
# gds <- estGeno(gds, optim = FALSE, call_threshold = 0.5)

# You can also set arbitrary values.
bias <- sample(seq(0, 1, 0.01), nmar(gds), replace = TRUE)
mismap <- cbind(sample(seq(0, 0.2, 0.01), nmar(gds), replace = TRUE),
                sample(seq(0, 0.2, 0.01), nmar(gds), replace = TRUE))
gds <- setFixedParameter(gds, bias = bias, mismap = mismap)

# Close the connection to the GDS file
closeGDS(gds)


}
