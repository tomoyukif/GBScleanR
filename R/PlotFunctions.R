#' Draw histograms of specified statistics
#'
#' @param x A [GbsrGenotypeData] object.
#' @param stats A string to specify statistics to be drawn.
#' @param target Either or both of "marker" and "sample", e.g. `target = "marker"`
#' to draw a histogram only for SNPs.
#' @param binwidth An integer to specify bin width of the histogram.
#' This value is passed to the ggplot function.
#' @param color A named vector "Marker" and "Sample" to specify border
#' color of bins in the histograms.
#' @param fill A named vector "Marker" and "Sample" to specify fill
#' color of bins in the histograms.
#'
#' @details
#' You can draw histograms of several summary statistics of genotype counts
#' and read counts per sample and per marker. The "stats" argument can take
#' the following values:
#' \describe{
#'   \item{missing}{Proportion of missing genotype calls.}
#'   \item{het}{Proportion of heterozygote calls.}
#'   \item{raf}{Reference allele frequency.}
#'   \item{dp}{Total read counts.}
#'   \item{ad_ref}{Reference allele read counts.}
#'   \item{ad_alt}{Alternative allele read counts.}
#'   \item{rrf}{Reference allele read frequency.}
#'   \item{mean_ref}{Mean of reference allele read counts.}
#'   \item{sd_ref}{Standard deviation of reference allele read counts.}
#'   \item{median_ref}{Quantile of reference allele read counts.}
#'   \item{mean_alt}{Mean of alternative allele read counts.}
#'   \item{sd_alt}{Standard deviation of alternative allele read counts.}
#'   \item{median_alt}{Quantile of alternative allele read counts.}
#'   \item{mq}{Mapping quality.}
#'   \item{fs}{Phred-scaled p-value (strand bias)}
#'   \item{qd}{Variant Quality by Depth}
#'   \item{sor}{Symmetric Odds Ratio (strand bias)}
#'   \item{mqranksum}{Alt vs. Ref read mapping qualities}
#'   \item{readposranksum}{Alt vs. Ref read position bias}
#'   \item{baseqranksum}{Alt Vs. Ref base qualities}
#' }
#'
#' To draw histograms for "missing", "het", "raf", you need to run
#' [countGenotype()] first to obtain statistics. Similary, "dp",
#'  "ad_ref", "ad_alt", "rrf" requires values obtained via
#'  [countRead()]. "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' and "baseqranksum" only work with `target = "marker"`, if your data
#'  contains those values supplied via SNP calling tools like
#' [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @export
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes xlab ylab xlim scale_fill_manual
#' @importFrom ggplot2 scale_color_manual theme element_text rel
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize genotype count information to be used in `histGBSR()`
#' gds <- countGenotype(gds)
#'
#' # Draw histograms of missing rate, heterozygosity, and reference
#' # allele frequency per SNP and per sample.
#' histGBSR(gds, stats = "missing")
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
histGBSR  <- function(x,
                      stats = c("dp", "missing", "het"),
                      target = c("marker", "sample"),
                      binwidth = NULL,
                      color = c(Marker = "darkblue", Sample = "darkblue"),
                      fill = c(Marker = "skyblue", Sample = "skyblue")) {
    stats_list <- c("dp",
                    "ad_ref",
                    "ad_alt",
                    "missing",
                    "het",
                    "raf",
                    "rrf",
                    "mean_ref",
                    "sd_ref",
                    "median_ref",
                    "mean_alt",
                    "sd_alt",
                    "median_alt",
                    "mq",
                    "fs",
                    "qd",
                    "sor",
                    "mqranksum",
                    "readposranksum",
                    "baseqranksum")

    if (missing(x)) {
        return(data.frame(
            stats = stats_list,
            description = c(
                "Total read depth",
                "Reference allele read depth",
                "Alternative allele read depth",
                "Missing rate",
                "Heterozygosity",
                "Referernce allele frequency",
                "Reference allele read frequency",
                "Mean of reference read depth",
                "SD of reference read depth",
                "Median of reference read depth",
                "Mean of alternative read depth",
                "SD of alternative read depth",
                "Median of alternative read depth",
                "Phred-scaled p-value (strand bias)",
                "Variant Quality by Depth",
                "Symmetric Odds Ratio (strand bias)",
                "Alt vs. Ref read mapping qualities",
                "Alt vs. Ref read position bias",
                "Alt Vs. Ref base qualities"
            )
        ))
    }

    stats <- match.arg(stats, stats_list)

    if (is.null(binwidth)) {
        if (stats %in% c("missing", "het", "raf", "rrf")) {
            binwidth <- 0.025
        }
    }

    df <- .df.maker(x, stats, target)
    p <- ggplot(df)
    p <- .hist.maker(p, stats, binwidth)
    p <- p + xlab(.lab.maker(stats)) +
        ylab("Count") +
        xlim(.limit.maker(stats)) +
        scale_fill_manual(values=fill, breaks=names(fill)) +
        scale_color_manual(values=color, breaks=names(color)) +
        theme(axis.title=element_text(face="bold", size=rel(1)),
              axis.text=element_text(size=rel(0.5)),
              strip.text=element_text(size=rel(0.8), face="bold"),
              legend.position="none")
    return(p)
}

#' Draw boxplots of specified statistics
#'
#' @param x A [GbsrGenotypeData] object.
#' @param stats A string to specify statistics to be drawn.
#' @param target Either or both of "marker" and "sample", e.g. `target = "marker"`
#' to draw a histogram only for SNPs.
#' @param color A named vector "Marker" and "Sample" to specify
#' border color of bins in the histograms.
#' @param fill A named vector "Marker" and "Sample" to specify
#' fill color of bins in the histograms.
#'
#' @details
#' You can draw boxplots of several summary statistics of genotype counts
#' and read counts per sample and per marker. The "stats" argument can take
#' the following values:
#'
#' \describe{
#'   \item{missing}{Proportion of missing genotype calls.}
#'   \item{het}{Proportion of heterozygote calls.}
#'   \item{raf}{Reference allele frequency.}
#'   \item{dp}{Total read counts.}
#'   \item{ad_ref}{Reference allele read counts.}
#'   \item{ad_alt}{Alternative allele read counts.}
#'   \item{rrf}{Reference allele read frequency.}
#'   \item{mean_ref}{Mean of reference allele read counts.}
#'   \item{sd_ref}{Standard deviation of reference allele read counts.}
#'   \item{median_ref}{Quantile of reference allele read counts.}
#'   \item{mean_alt}{Mean of alternative allele read counts.}
#'   \item{sd_alt}{Standard deviation of alternative allele read counts.}
#'   \item{median_alt}{Quantile of alternative allele read counts.}
#'   \item{mq}{Mapping quality.}
#'   \item{fs}{Phred-scaled p-value (strand bias)}
#'   \item{qd}{Variant Quality by Depth}
#'   \item{sor}{Symmetric Odds Ratio (strand bias)}
#'   \item{mqranksum}{Alt vs. Ref read mapping qualities}
#'   \item{readposranksum}{Alt vs. Ref read position bias}
#'   \item{baseqranksum}{Alt Vs. Ref base qualities}
#' }
#'
#' To draw boxplots for "missing", "het", "raf", you need to run
#' [countGenotype()] first to obtain statistics. Similary, "dp",
#' "ad_ref", "ad_alt", "rrf" requires values obtained via [countRead()].
#' "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' and "baseqranksum" only work with `target = "marker"`, if your data
#' contains those values supplied via SNP calling tools like
#' [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @export
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes xlab ylab xlim scale_fill_manual
#' @importFrom ggplot2 scale_color_manual theme element_text rel
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize genotype count information to be used in `boxplotGBSR()`
#' gds <- countGenotype(gds)
#'
#' boxplotGBSR(gds, stats = "missing")
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
boxplotGBSR <- function(x,
                        stats = "missing",
                        target = c("marker", "sample"),
                        color = c(Marker = "darkblue", Sample = "darkblue"),
                        fill = c(Marker = "skyblue", Sample = "skyblue")) {
    stats_list <- c("dp",
                    "ad_ref",
                    "ad_alt",
                    "missing",
                    "het",
                    "raf",
                    "rrf",
                    "mean_ref",
                    "sd_ref",
                    "median_ref",
                    "mean_alt",
                    "sd_alt",
                    "median_alt",
                    "mq",
                    "fs",
                    "qd",
                    "sor",
                    "mqranksum",
                    "readposranksum",
                    "baseqranksum")

    if (missing(x)) {
        return(data.frame(
            stats = stats_list,
            description = c(
                "Total read depth",
                "Reference allele read depth",
                "Alternative allele read depth",
                "Missing rate",
                "Heterozygosity",
                "Referernce allele frequency",
                "Reference allele read frequency",
                "Mean of reference read depth",
                "SD of reference read depth",
                "Median of reference read depth",
                "Mean of alternative read depth",
                "SD of alternative read depth",
                "Median of alternative read depth",
                "Phred-scaled p-value (strand bias)",
                "Variant Quality by Depth",
                "Symmetric Odds Ratio (strand bias)",
                "Alt vs. Ref read mapping qualities",
                "Alt vs. Ref read position bias",
                "Alt Vs. Ref base qualities"
            )
        ))
    }

    stats <- match.arg(stats, stats_list)

    df <- .df.maker(x, stats, target)
    p <- ggplot(df)
    p <- .boxplot.maker(p, stats)
    p <- p + ylab("") +
        xlab(.lab.maker(stats)) +
        xlim(.limit.maker(stats)) +
        scale_fill_manual(values=fill, breaks=names(fill)) +
        scale_color_manual(values=color, breaks=names(color)) +
        theme(axis.title=element_text(face="bold", size=rel(1)),
              axis.text=element_text(size=rel(0.5)),
              strip.text=element_text(size=rel(0.8), face="bold"),
              legend.position="none")
    return(p)
}


#' Draw line plots of specified statistics
#'
#' @param x A [GbsrGenotypeData] object.
#' @param stats A string to specify statistics to be drawn.
#' @param coord A vector with two integer specifying the number of rows
#' and columns to draw faceted line plots for chromosomes.
#' @param lwd A numeric value to specify the line width in plots.
#' @param binwidth An integer to specify bin width of the histogram.
#' This argument only work with `stats = "marker"` and is passed to the ggplot
#' function.
#' @param color A strings vector named "Marker", "Ref", "Het", "Alt"
#' to specify line colors. `stats = "geno` only requires "Ref", "Het" and "Alt",
#' while others uses the value named "Marker".
#'
#' @details
#' You can draw line plots of several summary statistics of genotype counts
#' and read counts per sample and per marker. The "stats" argument can take
#' the following values:
#' \describe{
#'   \item{marker}{Marker density.}
#'   \item{geno}{Proportion of missing genotype calls.}
#'   \item{missing}{Proportion of missing genotype calls.}
#'   \item{het}{Proportion of heterozygote calls.}
#'   \item{raf}{Reference allele frequency.}
#'   \item{dp}{Total read counts.}
#'   \item{ad_ref}{Reference allele read counts.}
#'   \item{ad_alt}{Alternative allele read counts.}
#'   \item{rrf}{Reference allele read frequency.}
#'   \item{mean_ref}{Mean of reference allele read counts.}
#'   \item{sd_ref}{Standard deviation of reference allele read counts.}
#'   \item{median_ref}{Quantile of reference allele read counts.}
#'   \item{mean_alt}{Mean of alternative allele read counts.}
#'   \item{sd_alt}{Standard deviation of alternative allele read counts.}
#'   \item{median_alt}{Quantile of alternative allele read counts.}
#'   \item{mq}{Mapping quality.}
#'   \item{fs}{Phred-scaled p-value (strand bias)}
#'   \item{qd}{Variant Quality by Depth}
#'   \item{sor}{Symmetric Odds Ratio (strand bias)}
#'   \item{mqranksum}{Alt vs. Ref read mapping qualities}
#'   \item{readposranksum}{Alt vs. Ref read position bias}
#'   \item{baseqranksum}{Alt Vs. Ref base qualities}
#' }
#'
#' To draw line plots for "missing", "het", "raf", you need to run
#' [countGenotype()] first to obtain statistics. Similary, "dp",
#' "ad_ref", "ad_alt", "rrf" requires values obtained via [countRead()].
#' "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' #' and "baseqranksum" only work with `target = "marker"`, if your data
#' contains those values supplied via SNP calling tools like
#' [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @export
#' @importFrom ggplot2 ggplot aes xlab ylab labs ylim scale_x_continuous
#' @importFrom ggplot2 scale_color_manual theme element_text rel expansion
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize genotype count information to be used in `plotGBSR()`
#' gds <- countGenotype(gds)
#'
#' # Draw line plots of missing rate, heterozygosity, proportion of genotype
#' # calls per SNP.
#' plotGBSR(gds, stats = "missing")
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
plotGBSR  <- function(x,
                      stats = c("dp", "missing", "het"),
                      coord = NULL,
                      lwd = 0.5,
                      binwidth = NULL,
                      color = c(
                          Marker = "darkblue",
                          Ref = "darkgreen",
                          Het = "magenta",
                          Alt = "blue")) {
    stats_list <- c("marker",
                    "geno",
                    "dp",
                    "ad_ref",
                    "ad_alt",
                    "missing",
                    "het",
                    "raf",
                    "rrf",
                    "mean_ref",
                    "sd_ref",
                    "median_ref",
                    "mean_alt",
                    "sd_alt",
                    "median_alt",
                    "mq",
                    "fs",
                    "qd",
                    "sor",
                    "mqranksum",
                    "readposranksum",
                    "baseqranksum"
    )
    if (missing(x)) {
        return(data.frame(
            stats = stats_list,
            description = c(
                "Marker density",
                "genotype ratio",
                "Total read depth",
                "Reference allele read depth",
                "Alternative allele read depth",
                "Missing rate",
                "Heterozygosity",
                "Referernce allele frequency",
                "Reference allele read frequency",
                "Mean of reference read depth",
                "SD of reference read depth",
                "Median of reference read depth",
                "Mean of alternative read depth",
                "SD of alternative read depth",
                "Median of alternative read depth",
                "Phred-scaled p-value (strand bias)",
                "Variant Quality by Depth",
                "Symmetric Odds Ratio (strand bias)",
                "Alt vs. Ref read mapping qualities",
                "Alt vs. Ref read position bias",
                "Alt Vs. Ref base qualities"
            )
        ))
    }
    stats <- match.arg(stats, stats_list)

    if(stats == "geno"){
        color <- color[2:4]
    }
    if(stats == "marker"){
        color <- color[1]
    }

    df <- .df.maker(x, stats, "marker", TRUE)
    p <- ggplot(df)
    p <- .plot.maker(p, stats, binwidth, lwd, coord)
    p <- p + xlab("Physical position (Mb)") +
        ylab("") +
        labs(title = .lab.maker(stats)) +
        ylim(.limit.maker(stats)) +
        scale_x_continuous(expand=expansion(0, 0.5), limits=c(0, NA)) +
        scale_color_manual(values=color, breaks=names(color)) +
        theme(
            plot.title=element_text(face="bold", size=rel(1)),
            axis.title=element_text(face="bold", size=rel(0.7)),
            axis.text=element_text(size=rel(0.5)),
            strip.text.y=element_text(size=rel(0.8), face="bold", angle=0)
        )
    if (stats != "geno") {
        p <- p + theme(legend.position="none")
    }
    return(p)
}



#' Draw a scatter plot of a pair of specified statistics
#'
#' @param x A [GbsrGenotypeData] object.
#' @param stats1 A string to specify statistics to be drawn.
#' @param stats2 A string to specify statistics to be drawn.
#' @param target Either or both of "marker" and "sample", e.g.
#' `target = "marker"` to draw a histogram only for SNPs.
#' @param size A numeric value to specify the dot size of a scatter plot.
#' @param alpha A numeric value \[0-1\] to specify the transparency of
#' dots in a scatter plot.
#' @param color A named vector "Marker" and "Sample" to specify border
#' color of bins in the histograms.
#' @param fill A named vector "Marker" and "Sample" to specify fill color
#' of bins in the histograms.`stats = "geno` only requires "Ref", "Het"
#' and "Alt", while others uses the value named "Marker".
#' @param smooth A logical value to indicate whether draw a smooth line for
#' data points. See also [ggplot2::stat_smooth()].
#'
#' @details
#' You can draw a scatter plot of per-marker and/or per-sample summary
#' statistics specified at `stats1` and `stats2`. The "stats1" and "stats2"
#' arguments can take the following values:
#' \describe{
#'   \item{missing}{Proportion of missing genotype calls.}
#'   \item{het}{Proportion of heterozygote calls.}
#'   \item{raf}{Reference allele frequency.}
#'   \item{dp}{Total read counts.}
#'   \item{ad_ref}{Reference allele read counts.}
#'   \item{ad_alt}{Alternative allele read counts.}
#'   \item{rrf}{Reference allele read frequency.}
#'   \item{mean_ref}{Mean of reference allele read counts.}
#'   \item{sd_ref}{Standard deviation of reference allele read counts.}
#'   \item{median_ref}{Quantile of reference allele read counts.}
#'   \item{mean_alt}{Mean of alternative allele read counts.}
#'   \item{sd_alt}{Standard deviation of alternative allele read counts.}
#'   \item{median_alt}{Quantile of alternative allele read counts.}
#'   \item{mq}{Mapping quality.}
#'   \item{fs}{Phred-scaled p-value (strand bias)}
#'   \item{qd}{Variant Quality by Depth}
#'   \item{sor}{Symmetric Odds Ratio (strand bias)}
#'   \item{mqranksum}{Alt vs. Ref read mapping qualities}
#'   \item{readposranksum}{Alt vs. Ref read position bias}
#'   \item{baseqranksum}{Alt Vs. Ref base qualities}
#' }
#'
#' To draw scatter plots for "missing", "het", "raf", you need to run
#' [countGenotype()] first to obtain statistics. Similary, "dp",
#' "ad_ref", "ad_alt", "rrf" requires values obtained via [countRead()].
#' "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' and "baseqranksum" only work with `target = "marker"`, if your data
#' contains those values supplied via SNP calling tools like
#' [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @return A ggplot object.
#'
#' @export
#' @importFrom ggplot2 ggplot aes xlab ylab xlim ylim scale_fill_manual
#' @importFrom ggplot2 scale_color_manual theme stat_smooth element_text rel
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' # Summarize genotype count information to be used in `pairsGBSR()`
#' gds <- countGenotype(gds)
#'
#' # Draw scatter plots of missing rate vs heterozygosity.
#' pairsGBSR(gds, stats1 = "missing", stats2 = "het")
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#'
pairsGBSR  <- function(x,
                       stats1 = "dp",
                       stats2 = "missing",
                       target = "marker",
                       size = 0.5,
                       alpha = 0.8,
                       color = c(Marker = "darkblue", Sample = "darkblue"),
                       fill = c(Marker = "skyblue", Sample = "skyblue"),
                       smooth = FALSE) {
    stats_list <- c("dp",
                    "ad_ref",
                    "ad_alt",
                    "missing",
                    "het",
                    "raf",
                    "rrf",
                    "mean_ref",
                    "sd_ref",
                    "median_ref",
                    "mean_alt",
                    "sd_alt",
                    "median_alt",
                    "mq",
                    "fs",
                    "qd",
                    "sor",
                    "mqranksum",
                    "readposranksum",
                    "baseqranksum")

    if (missing(x)) {
        return(data.frame(
            stats = stats_list,
            description = c(
                "Total read depth",
                "Reference allele read depth",
                "Alternative allele read depth",
                "Missing rate",
                "Heterozygosity",
                "Referernce allele frequency",
                "Reference allele read frequency",
                "Mean of reference read depth",
                "SD of reference read depth",
                "Median of reference read depth",
                "Mean of alternative read depth",
                "SD of alternative read depth",
                "Median of alternative read depth",
                "Mapping quality",
                "Phred-scaled p-value (strand bias)",
                "Variant Quality by Depth",
                "Symmetric Odds Ratio (strand bias)",
                "Alt vs. Ref read mapping qualities",
                "Alt vs. Ref read position bias",
                "Alt Vs. Ref base qualities"
            )
        ))
    }

    stats1 <-
        match.arg(stats1, stats_list)
    stats2 <-
        match.arg(stats2, stats_list)

    df1 <- .df.maker(x, stats1, target)
    df2 <- .df.maker(x, stats2, target)
    df <-
        data.frame(val1 = df1$val,
                   val2 = df2$val,
                   target = df1$target)

    val1 <- val2 <- NULL
    p <- ggplot(df)
    p <- .pairs.maker(p, size, alpha)
    p <- p + xlab(.lab.maker(stats1)) +
        ylab(.lab.maker(stats2)) +
        xlim(.limit.maker(stats1)) +
        ylim(.limit.maker(stats2)) +
        scale_fill_manual(values=fill, breaks=names(fill)) +
        scale_color_manual(values=color, breaks=names(color)) +
        theme(
            plot.title=element_text(face="bold", size=rel(1)),
            axis.title=element_text(face="bold", size=rel(0.7)),
            axis.text=element_text(size=rel(0.5)),
            strip.text.y=element_text(size=rel(0.8), face="bold", angle=0),
            legend.position="none"
        )
    if (smooth) {
        p <- p + stat_smooth(aes(x=val1, y=val2),)
    }
    return(p)
}


# Internal function to build a data.frame passed to ggplot().
#' @importFrom tidyr pivot_longer
.df.maker <- function(x, stats, target, pos = FALSE) {
    Ref <- NULL
    Alt <- NULL
    if (stats == "geno") {
        snp <- data.frame(Ref = getCountGenoRef(x, "marker", TRUE, TRUE),
                          Het = getCountGenoHet(x, "marker", TRUE, TRUE),
                          Alt = getCountGenoAlt(x, "marker", TRUE, TRUE),
                          chr = getChromosome(x),
                          pos = getPosition(x) * 10^-6,
                          stringsAsFactors=FALSE)
        if (is.null(snp)) {
            msg <- paste0('No data for the statistic: ', stats)
            stop(msg)
        }
        snp <- pivot_longer(snp, Ref:Alt,
                            names_to="genotype", values_to="val")
        scan <- NULL

    } else if (stats == "marker") {
        snp <- data.frame(chr = getChromosome(x),
                          pos = getPosition(x) * 10^-6,
                          target = "Marker",
                          stringsAsFactors=FALSE)
        scan <- NULL

    } else {
        if ("marker" %in% target) {
            snp <- switch(
                stats,
                "dp" = getCountRead(x, "marker"),
                "missing" = getCountGenoMissing(x, "marker", TRUE, TRUE),
                "het" = getCountGenoHet(x, "marker", TRUE, TRUE),
                "raf" = getCountAlleleRef(x, "marker", TRUE, TRUE),
                "rrf" = getCountReadRef(x, "marker", TRUE, TRUE),
                "ad_ref" = getCountReadRef(x, "marker"),
                "mean_ref" = getMeanReadRef(x, "marker"),
                "sd_ref" = getSDReadRef(x, "marker"),
                "median_ref" = getMedianReadRef(x, "marker"),
                "ad_alt" = getCountReadAlt(x, "marker"),
                "mean_alt" = getMeanReadAlt(x, "marker"),
                "sd_alt" = getSDReadAlt(x, "marker"),
                "median_alt" = getMedianReadAlt(x, "marker"),
                "mq" = getInfo(x, "MQ"),
                "fs" = getInfo(x, "FS"),
                "qd" = getInfo(x, "QD"),
                "sor" = getInfo(x, "SOR"),
                "mqranksum" = getInfo(x, "MQRankSum"),
                "readposranksum" = getInfo(x, "ReadPosRankSum"),
                "baseqranksum" = getInfo(x, "BaseQRankSum")
            )
            if (is.null(snp)) {
                msg <- paste0('No SNP summary for the statistic: ', stats)
                warning(msg)
            } else {
                if (pos) {
                    snp <- data.frame(val = snp,
                                      target = "Marker",
                                      chr = getChromosome(x),
                                      pos = getPosition(x) * 10^-6,
                                      stringsAsFactors=FALSE)
                } else {
                    snp <- data.frame(val = snp, target = "Marker")
                }
            }
        } else {
            snp <- NULL
        }
        if ("sample" %in% target) {
            scan <- switch(
                stats,
                "dp" = getCountRead(x, "sample"),
                "missing" = getCountGenoMissing(x, "sample", TRUE, TRUE),
                "het" = getCountGenoHet(x, "sample", TRUE, TRUE),
                "raf" = getCountAlleleRef(x, "sample", TRUE, TRUE),
                "rrf" = getCountReadRef(x, "sample", TRUE, TRUE),
                "ad_ref" = getCountReadRef(x, "sample"),
                "mean_ref" = getMeanReadRef(x, "sample"),
                "sd_ref" = getSDReadRef(x, "sample"),
                "median_ref" = getMedianReadRef(x, "sample"),
                "ad_alt" = getCountReadAlt(x, "sample"),
                "mean_alt" = getMeanReadAlt(x, "sample"),
                "sd_alt" = getSDReadAlt(x, "sample"),
                "median_alt" = getMedianReadAlt(x, "sample")
            )
            if (is.null(scan)) {
                msg <- paste0('No sample summary for the statistic: ', stats)
                warning(msg)
            } else {
                scan <- data.frame(val = scan, target = "Sample")
            }
        } else {
            scan <- NULL
        }
    }
    df <- rbind(snp, scan)
    return(df)
}

# Internal function to build label information passed to ggplot().
.lab.maker <- function(stats) {
    lab <- switch(
        stats,
        "marker" = "Marker density",
        "geno" = "Genotype ratio",
        "dp" = "Total read depth",
        "ad_ref" = "Reference allelic read depth",
        "ad_alt" = "Alternative allelic read depth",
        "missing" = "Missing rate",
        "het" = "Heterozygosity",
        "raf" = "Referernce allele frequency",
        "rrf" = "Reference allele read frequency",
        "mean_ref" = "Mean of reference read depth",
        "sd_ref" = "SD of reference read depth",
        "median_ref" = "Median of reference read depth",
        "mean_alt" = "Mean of alternative read depth",
        "sd_alt" = "SD of alternative read depth",
        "median_alt" = "Median of alternative read depth",
        "mq" = "Mapping quality (MQ)",
        "fs" = "Phred-scaled p-value (strand bias) (FS)",
        "qd" = "Variant Quality by Depth (QD)",
        "sor" = "Symmetric Odds Ratio (strand bias) (SOR)",
        "mqranksum" = "Alt vs. Ref read mapping qualities (MQRankSum)",
        "readposranksum" = "Alt vs. Ref read position bias (ReadPosRankSum)",
        "baseqranksum" = "Alt Vs. Ref base qualities (BaseQRankSum)"
    )
    return(lab)
}

# Internal function to build plot limit information passed to ggplot().
.limit.maker <- function(stats) {
    lim <- switch(
        stats,
        "marker" = c(0, NA),
        "geno" = c(0, 1),
        "dp" = c(0, NA),
        "ad_ref" = c(0, NA),
        "ad_alt" = c(0, NA),
        "missing" = c(0, 1),
        "het" = c(0, 1),
        "raf" = c(0, 1),
        "rrf" = c(0, 1),
        "mean_ref" = c(0, NA),
        "sd_ref" = c(0, NA),
        "median_ref" = c(0, NA),
        "mean_alt" = c(0, NA),
        "sd_alt" = c(0, NA),
        "median_alt" = c(0, NA),
        "mq" = c(0, NA),
        "fs" = c(0, NA),
        "qd" = c(0, NA),
        "sor" = c(0, NA),
        "mqranksum" = as.numeric(c(NA, NA)),
        "readposranksum" = as.numeric(c(NA, NA)),
        "baseqranksum" = as.numeric(c(NA, NA))
    )
    return(lim)
}

# Internal function to draw a histogram.
#' @importFrom ggplot2 geom_histogram facet_wrap
.hist.maker <- function(p, stats, binwidth, color, fill) {
    val <- target <- NULL
    p <- p +
        geom_histogram(aes(x=val, color=target, fill=target),
                       binwidth=binwidth,
                       boundary=0) +
        facet_wrap( ~ target, scales="free")
    return(p)
}

# Internal function to draw a boxplot.
#' @importFrom ggplot2 geom_histogram facet_wrap
.boxplot.maker <- function(p, stats) {
    val <- NULL
    target <- NULL
    p <- p + geom_histogram(aes(x=val, color=target)) +
        facet_wrap( ~ target, scales="free")
    return(p)
}

# Internal function to draw a line plot.
#' @importFrom ggplot2 geom_line geom_histogram facet_wrap
.plot.maker <- function(p, stats, binwidth, lwd, coord) {
    val <- genotype <- target <- pos <- NULL
    if ("geno" %in% stats) {
        p <-p +
            geom_line(aes(x=pos, y=val, color=genotype),
                      size=lwd) +
            facet_wrap(~ chr, coord[1], coord[2], "free_x",
                       dir="v", strip.position="right")
    } else if ("marker" %in% stats) {
        p <- p +
            geom_histogram(aes(x=pos, fill=target, color=target),
                           binwidth=binwidth, boundary=0) +
            facet_wrap(~ chr, coord[1], coord[2], "free_x",
                       dir="v", strip.position="right")
    } else {
        p <- p +
            geom_line(aes(x=pos, y=val, color=target), size=lwd) +
            facet_wrap(~ chr, coord[1], coord[2], "free_x",
                       dir="v", strip.position="right")
    }
    return(p)
}

# Internal function to draw a scatter plot.
#' @importFrom ggplot2 geom_point facet_wrap
.pairs.maker <- function(p, size, alpha) {
    val1 <- val2 <- target <- NULL
    p <- p +
        geom_point(aes(x=val1, y=val2, color=target, fill=target),
                   size=size, alpha=alpha) +
        facet_wrap( ~ target, scales="free")
    return(p)
}

#' Draw line plots of allele dosage per marker per sample.
#'
#' This function counts a reference allele dosage per marker per sample and
#' draw line plots of them in facets for each chromosome for each sample.
#'
#' @param x A [GbsrGenotypeData] object.
#' @param coord A vector with two integer specifying the number of rows
#' and columns to draw faceted line plots for chromosomes.
#' @param chr A vector of indexes to specify chromosomes to be drawn.
#' @param ind An index to specify samples to be drawn.
#' @param node Either one of "raw" or "filt" to output raw read data, or filtered read data,
#' respectively.
#' @param showratio If `TRUE`, draw dots indicating read ratio.
#' @param line_color A string to indicate the line color in the plot.
#' @param dot_fill A vector of two strings to indicate the dot colors in the plot. The first and second elements of the vector are set as the colors for the lowest and highest values in the gradient coloring of the dots indicating total read counts par marker.
##' @param size A positive number to indicate the dot size in a plot.
#' @param alpha A positive number in 0-1 to indicate the dot opacity in a plot.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' plotDosage(gds, ind = 1)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs
#' @importFrom ggplot2 ylim xlab ylab facet_wrap theme
#' @importFrom ggplot2 scale_colour_gradient scale_size_manual
#' @importFrom ggplot2 scale_linewidth_manual
#'
plotDosage <- function(x,
                       coord = NULL,
                       chr = NULL,
                       ind = 1,
                       node = "raw",
                       showratio = TRUE,
                       dot_fill = c("green", "darkblue"),
                       size = 0.8,
                       alpha = 0.8,
                       line_color = "magenta") {
    pos <- na <- NULL

    if (is.null(chr)) {
        chr <- rep(TRUE, nmar(x))
    } else {
        chr <- getChromosome(x) %in% chr
    }

    if(is.character(ind)){
        parents <- TRUE
        ind <- which(getSamID(x, valid = FALSE) %in% ind)
        id <- getSamID(x, valid = FALSE)[ind]

    } else {
        if(is.null(slot(x, "sample")[["parents"]])){
            parents <- FALSE
        } else {
            parents <- TRUE
        }
        id <- getSamID(x, valid = FALSE)[validSam(x, parents = parents)][ind]
    }

    if(showratio){
        ploidy <- attributes(slot(x, "sample"))[["ploidy"]]
        read <- getRead(x, node = node, valid = TRUE, parents = parents)
        ref <- read$ref[ind, chr]
        alt <- read$alt[ind, chr]
        dp <- ref + alt
        ad <- alt / dp * ploidy
    }

    geno <- getGenotype(x, node = "dosage", valid = TRUE, parents = parents)[ind, chr]
    if(all(is.na(geno))){
        warning("Missing values at all markers./n",
                "This sample might have no read at all markers.")
    }

    pos <- getPosition(x)[chr]
    chr <- getChromosome(x)[chr]
    max_pos <- tapply(pos, chr, max)
    min_pos <- tapply(pos, chr, min)
    uniq_chr <- sort(unique(chr))
    dummy <- rbind(data.frame(chr = uniq_chr, pos = max_pos, ad = 0, dp = 0),
                   data.frame(chr = uniq_chr, pos = min_pos, ad = 0, dp = 0))

    df <- data.frame(chr = chr, pos = pos, geno = geno, ad = ad, dp = dp)
    df <- df[!is.na(df$geno), ]
    ploidy <- attributes(slot(x, "sample"))[["ploidy"]]
    p <- ggplot(df) +
        geom_point(data = dummy, mapping = aes(x = pos * 10^-6, y = ad),
                   size = 0)

    if(showratio){
        p <- p + geom_point(mapping = aes(x = pos * 10^-6, y = ad, colour = dp),
                            size = size, alpha = alpha, stroke = 0) +
            scale_colour_gradient(low = dot_fill[1], high = dot_fill[2])
    }
    p <- p + geom_line(mapping = aes(x = pos * 10^-6, y = geno, group = chr),
                       color = line_color) +
        labs(title=paste0("Alternative allele dosage: ", id)) +
        ylim(0, ploidy) +
        xlab("Physical position (Mb)") +
        ylab("Alternative allele dosage")

    p <- p + facet_wrap(~ chr, coord[1], coord[2], "free_x",
                        dir="v", strip.position="right") +
        theme(legend.position="none")
    return(p)
}

#' Draw line plots of proportion of reference allele read counts
#'  per marker per sample.
#'
#' This function calculate a proportion of reference allele read
#' counts per marker per sample and
#' draw line plots of them in facets for each chromosome for each sample.
#'
#' @param x A [GbsrGenotypeData] object.
#' @param coord A vector with two integer specifying the number of rows and
#' columns to draw faceted line plots for chromosomes.
#' @param chr A vector of indexes to specify chromosomes to be drawn.
#' @param ind A string of sample id or an index to specify the sample to be drawn.
#' @param node Either one of "raw" or "filt" to output raw read data, or filtered read data,
#' respectively.
#' @param dot_fill A vector of two strings to indicate the dot colors in the plot. The first and second elements of the vector are set as the colors for the lowest and highest values in the gradient coloring of the dots indicating total read counts par marker.
#' @param size A positive number to indicate the dot size in the plot.
#' @param alpha A positive number in 0-1 to indicate the dot opacity in the plot.
#'
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gds <- loadGDS(gds_fn)
#'
#' plotReadRatio(gds, ind = 1)
#'
#' # Close the connection to the GDS file
#' closeGDS(gds)
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point labs
#' @importFrom ggplot2 ylim xlab ylab facet_wrap theme
#' @importFrom ggplot2 scale_colour_gradient scale_size_manual
#'
plotReadRatio <- function(x,
                          coord = NULL,
                          chr = NULL,
                          ind = 1,
                          node = "raw",
                          dot_fill = c("green", "darkblue"),
                          size = 0.8,
                          alpha = 0.8) {
    pos <- ad <- NULL

    if (is.null(chr)) {
        chr <- rep(TRUE, nmar(x))
    } else {
        chr <- getChromosome(x) %in% chr
    }


    if(is.character(ind)){
        parents <- TRUE
        ind <- which(getSamID(x, valid = FALSE) %in% ind)
        id <- getSamID(x, valid = FALSE)[ind]

    } else {
        if(is.null(slot(x, "sample")[["parents"]])){
            parents <- FALSE
        } else {
            parents <- TRUE
        }
        id <- getSamID(x)[ind]
    }

    read <- getRead(x, node = node, valid = TRUE, parents = parents)
    ref <- read$ref[ind, chr]
    alt <- read$alt[ind, chr]
    dp <- ref + alt
    ad <- alt / dp

    if(all(is.na(ad))){
        warning("Missing values at all markers./n",
                "This sample might have no read at all markers.")
    }

    pos <- getPosition(x)[chr]
    chr <- getChromosome(x)[chr]
    max_pos <- tapply(pos, chr, max)
    min_pos <- tapply(pos, chr, min)
    uniq_chr <- sort(unique(chr))
    dummy <- rbind(data.frame(chr = uniq_chr, pos = max_pos, ad = 0, dp = 0),
                   data.frame(chr = uniq_chr, pos = min_pos, ad = 0, dp = 0))

    df <- data.frame(chr = chr, pos = pos, ad = ad, dp = dp)

    p <- ggplot(df, aes(x = pos * 10^-6, y = ad, colour = dp)) +
        geom_point(data = dummy, mapping = aes(x = pos * 10^-6, y = ad),
                   size = 0) +
        geom_point(size = size, alpha = alpha, stroke = 0) +
        scale_colour_gradient(low = dot_fill[1], high = dot_fill[2]) +
        labs(title=paste0("Alternative allele read ratio: ", id)) +
        ylim(0, 1) +
        xlab("Position (Mb)") +
        ylab("Alternative allele read ratio")

    p <- p + facet_wrap(~ chr, coord[1], coord[2], "free_x",
                        dir="v", strip.position="right") +
        theme(legend.position="none")
    return(p)
}
