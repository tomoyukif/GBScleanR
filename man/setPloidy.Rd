% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/Methods-GbsrGenotypeData.R
\name{setPloidy}
\alias{setPloidy}
\alias{setPloidy,GbsrGenotypeData-method}
\title{Set ploidy}
\usage{
setPloidy(object, ploidy = 2, ...)

\S4method{setPloidy}{GbsrGenotypeData}(object, ploidy = 2)
}
\arguments{
\item{object}{A \link{GbsrGenotypeData} object.}

\item{ploidy}{A integer value to specify the ploidy of the given population.}

\item{...}{Unused.}
}
\value{
A \link{GbsrGenotypeData} object with filters on markers.
}
\description{
Set the ploidy of a given population in the input \link{GbsrGenotypeData} object.
}
\details{
The genotype estimation by \code{\link[=estGeno]{estGeno()}} would be performed with the assumption
of the ploidy specified through this function or the \code{ploidy} argument of
\code{\link[=loadGDS]{loadGDS()}}.
When an odd number was specified as \code{ploidy}, the ploidy of intermediate
generations would be treated as \code{ploidy} + 1 to properly list up possible
descendent haplotype patterns in the process by \code{\link[=estGeno]{estGeno()}}.
}
\examples{
# Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
gds <- loadGDS(gds_fn)

gds <- setPloidy(gds, ploidy = 4)

# Close the connection to the GDS file.
closeGDS(gds)

}
\seealso{
\code{\link[=getPloidy]{getPloidy()}}
}
