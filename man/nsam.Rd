% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/Methods-GbsrGenotypeData.R
\name{nsam}
\alias{nsam}
\alias{nsam,GbsrGenotypeData-method}
\title{Return the number of samples.}
\usage{
nsam(object, valid = TRUE, parents = FALSE, ...)

\S4method{nsam}{GbsrGenotypeData}(object, valid, parents)
}
\arguments{
\item{object}{A \link{GbsrGenotypeData} object.}

\item{valid}{A logical value. See details.}

\item{parents}{A logical value whether to include to parental samples or not.}

\item{...}{Unused.}
}
\value{
An integer value to indicate the number of samples.
}
\description{
This function returns the number of samples recorded in the GDS file
connected to the given \link{GbsrGenotypeData} object.
}
\details{
If \code{valid = TRUE}, the number of the samples which are labeled \code{TRUE}
in the "valid" column of the "sample" slot will be returned. If you need
the number of over all samples, set \code{valid = FALSE}.
\code{\link[=validSam]{validSam()}} tells you which samples are valid.
}
\examples{
# Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
gds <- loadGDS(gds_fn)

nsam(gds)

# Close the connection to the GDS file.
closeGDS(gds)

}
\seealso{
\code{\link[=validSam]{validSam()}}
}
