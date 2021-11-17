#' Draw histograms of specified statistics
#'
#' @param x A [GbsrGenotypeData] object.
#' @param stats A string to specify statistics to be drawn.
#' @param target Either or both of "snp" and "scan", e.g. `target = "snp"`
#' to draw a histogram only for SNPs.
#' @param q An integer to specify a quantile calculated via [calcReadStats()].
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
#' \itemize{
#' \item{"missing"}{"Proportion of missing genotype calls."},
#' \item{"het"}{"Proportion of heterozygote calls."},
#' \item{"raf"}{"Reference allele frequency."},
#' \item{"dp"}{"Total read counts."},
#' \item{"ad_ref"}{"Reference allele read counts."},
#' \item{"ad_alt"}{"Alternative allele read counts."},
#' \item{"rrf"}{"Reference allele read frequency."},
#' \item{"mean_ref"}{"Mean of reference allele read counts."},
#' \item{"sd_ref"}{"Standard deviation of reference allele read counts."},
#' \item{"qtile_ref"}{"Quantile of reference allele read counts."},
#' \item{"mean_alt"}{"Mean of alternative allele read counts."},
#' \item{"sd_alt"}{"Standard deviation of alternative allele read counts."},
#' \item{"qtile_alt"}{"Quantile of alternative allele read counts."},
#' \item{"mq"}{"Mapping quality."},
#' \item{"fs"}{"Phred-scaled p-value (strand bias)"},
#' \item{"qd"}{"Variant Quality by Depth"},
#' \item{"sor"}{"Symmetric Odds Ratio (strand bias)"},
#' \item{"mqranksum"}{"Alt vs. Ref read mapping qualities"},
#' \item{"readposranksum"}{"Alt vs. Ref read position bias"},
#' \item{"baseqranksum"}{"Alt Vs. Ref base qualities"},
#' }
#'
#' To draw histograms for "missing", "het", "raf", you need to run
#' [countGenotype()] first to obtain statistics. Similary, "dp",
#'  "ad_ref", "ad_alt", "rrf" requires values obtained via
#'  [countRead()]. [calcReadStats()] should be executed before drawing
#'  histograms of "mean_ref", "sd_ref", "qtile_ref", "mean_alt", "sd_alt",
#' and "qtile_alt". "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' and "baseqranksum" only work with `target = "snp"`, if your data
#'  contains those values supplied via SNP calling tools like
#' [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @export
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @import graphics
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#'
#' # Summarize genotype count information to be used in `histGBSR()`
#' gdata <- countGenotype(gdata)
#'
#' # Draw histograms of missing rate, heterozygosity, and reference
#' # allele frequency per SNP and per sample.
#' histGBSR(gdata, stats = "missing")
#'
#'
#' # Calculate means, standard deviations, quantile values of read counts
#' # to be used in `histGBSR()`
#' gdata <- calcReadStats(gdata, q = 0.9)
#'
#' # Draw histograms of 90 percentile values of reference read counts
#' # and alternative read counts per SNP and per sample.
#' histGBSR(gdata, stats = "qtile_ref", q = 0.9)
#'
#' # Close the connection to the GDS file
#' closeGDS(gdata)
histGBSR  <- function(x,
                      stats = c("dp", "missing", "het"),
                      target = c("snp", "scan"),
                      q = 0.5,
                      binwidth = NULL,
                      color = c(Marker = "darkblue", Sample = "darkblue"),
                      fill = c(Marker = "skyblue", Sample = "skyblue")) {
    stats_list = c("dp",
                   "ad_ref",
                   "ad_alt",
                   "missing",
                   "het",
                   "raf",
                   "rrf",
                   "mean_ref",
                   "sd_ref",
                   "qtile_ref",
                   "mean_alt",
                   "sd_alt",
                   "qtile_alt",
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
                paste0("Quantile of reference read depth (q=", q, ")"),
                "Mean of alternative read depth",
                "SD of alternative read depth",
                paste0("Quantile of alternative read depth (q=", q, ")"),
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

    df <- .df.maker(x, stats, q, target)
    p <- ggplot(df)
    p <- .hist.maker(p, stats, binwidth)
    p <- p + xlab(.lab.maker(stats, q)) +
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
#' @param target Either or both of "snp" and "scan", e.g. `target = "snp"`
#' to draw a histogram only for SNPs.
#' @param q An integer to specify a quantile calculated via
#' [calcReadStats()].
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
#' \itemize{
#' \item{"missing"}{"Proportion of missing genotype calls."},
#' \item{"het"}{"Proportion of heterozygote calls."},
#' \item{"raf"}{"Reference allele frequency."},
#' \item{"dp"}{"Total read counts."},
#' \item{"ad_ref"}{"Reference allele read counts."},
#' \item{"ad_alt"}{"Alternative allele read counts."},
#' \item{"rrf"}{"Reference allele read frequency."},
#' \item{"mean_ref"}{"Mean of reference allele read counts."},
#' \item{"sd_ref"}{"Standard deviation of reference allele read counts."},
#' \item{"qtile_ref"}{"Quantile of reference allele read counts."},
#' \item{"mean_alt"}{"Mean of alternative allele read counts."},
#' \item{"sd_alt"}{"Standard deviation of alternative allele read counts."},
#' \item{"qtile_alt"}{"Quantile of alternative allele read counts."},
#' \item{"mq"}{"Mapping quality."},
#' \item{"fs"}{"Phred-scaled p-value (strand bias)"},
#' \item{"qd"}{"Variant Quality by Depth"},
#' \item{"sor"}{"Symmetric Odds Ratio (strand bias)"},
#' \item{"mqranksum"}{"Alt vs. Ref read mapping qualities"},
#' \item{"readposranksum"}{"Alt vs. Ref read position bias"},
#' \item{"baseqranksum"}{"Alt Vs. Ref base qualities"},
#' }
#'
#' To draw boxplots for "missing", "het", "raf", you need to run
#' [countGenotype()] first to obtain statistics. Similary, "dp",
#' "ad_ref", "ad_alt", "rrf" requires values obtained via [countRead()].
#' [calcReadStats()] should be executed before drawing boxplots of
#' "mean_ref", "sd_ref", "qtile_ref", "mean_alt", "sd_alt", and
#' "qtile_alt". "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' and "baseqranksum" only work with `target = "snp"`, if your data
#' contains those values supplied via SNP calling tools like
#' [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @export
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @import graphics
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#'
#' # Summarize genotype count information to be used in `boxplotGBSR()`
#' gdata <- countGenotype(gdata)
#'
#' boxplotGBSR(gdata, stats = "missing")
#'
#' # Calculate means, standard deviations, quantile values of read counts
#' # to be used in `boxplotGBSR()`
#' gdata <- calcReadStats(gdata, q = 0.9)
#'
#' # Draw boxplots of 90 percentile values of reference read counts and
#' # alternative read counts per SNP and per sample.
#' boxplotGBSR(gdata, stats = "qtile_ref", q = 0.9)
#'
#' # Close the connection to the GDS file
#' closeGDS(gdata)
#'
boxplotGBSR <- function(x,
                        stats = "missing",
                        target = c("snp", "scan"),
                        q = 0.5,
                        color = c(Marker = "darkblue", Sample = "darkblue"),
                        fill = c(Marker = "skyblue", Sample = "skyblue")) {
    stats_list = c("dp",
                   "ad_ref",
                   "ad_alt",
                   "missing",
                   "het",
                   "raf",
                   "rrf",
                   "mean_ref",
                   "sd_ref",
                   "qtile_ref",
                   "mean_alt",
                   "sd_alt",
                   "qtile_alt",
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
                paste0("Quantile of reference read depth (q=", q, ")"),
                "Mean of alternative read depth",
                "SD of alternative read depth",
                paste0("Quantile of alternative read depth (q=", q, ")"),
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

    df <- .df.maker(x, stats, q, target)
    p <- ggplot(df)
    p <- .boxplot.maker(p, stats)
    p <- p + ylab("") +
        xlab(.lab.maker(stats, q)) +
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
#' @param q An integer to specify a quantile calculated via
#' [calcReadStats()].
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
#' \itemize{
#' \item{"marker"}{"Marker density."},
#' \item{"geno"}{"Proportion of missing genotype calls."},
#' \item{"missing"}{"Proportion of missing genotype calls."},
#' \item{"het"}{"Proportion of heterozygote calls."},
#' \item{"raf"}{"Reference allele frequency."},
#' \item{"dp"}{"Total read counts."},
#' \item{"ad_ref"}{"Reference allele read counts."},
#' \item{"ad_alt"}{"Alternative allele read counts."},
#' \item{"rrf"}{"Reference allele read frequency."},
#' \item{"mean_ref"}{"Mean of reference allele read counts."},
#' \item{"sd_ref"}{"Standard deviation of reference allele read counts."},
#' \item{"qtile_ref"}{"Quantile of reference allele read counts."},
#' \item{"mean_alt"}{"Mean of alternative allele read counts."},
#' \item{"sd_alt"}{"Standard deviation of alternative allele read counts."},
#' \item{"qtile_alt"}{"Quantile of alternative allele read counts."},
#' \item{"mq"}{"Mapping quality."},
#' \item{"fs"}{"Phred-scaled p-value (strand bias)"},
#' \item{"qd"}{"Variant Quality by Depth"},
#' \item{"sor"}{"Symmetric Odds Ratio (strand bias)"},
#' \item{"mqranksum"}{"Alt vs. Ref read mapping qualities"},
#' \item{"readposranksum"}{"Alt vs. Ref read position bias"},
#' \item{"baseqranksum"}{"Alt Vs. Ref base qualities"},
#' }
#'
#' To draw line plots for "missing", "het", "raf", you need to run
#' [countGenotype()] first to obtain statistics. Similary, "dp",
#' "ad_ref", "ad_alt", "rrf" requires values obtained via [countRead()].
#' [calcReadStats()] should be executed before drawing line plots of
#' "mean_ref", "sd_ref", "qtile_ref", "mean_alt", "sd_alt", and
#' "qtile_alt". "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' #' and "baseqranksum" only work with `target = "snp"`, if your data
#' contains those values supplied via SNP calling tools like
#' [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @export
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @import graphics
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#'
#' # Summarize genotype count information to be used in `plotGBSR()`
#' gdata <- countGenotype(gdata)
#'
#' # Draw line plots of missing rate, heterozygosity, proportion of genotype
#' # calls per SNP.
#' plotGBSR(gdata, stats = "missing")
#'
#' # Calculate means, standard deviations, quantile values of read counts
#' # to be used in `plotGBSR()`
#' gdata <- calcReadStats(gdata, q = 0.9)
#'
#' # Draw line plots of 90 percentile values of reference read counts and
#' # alternative read counts per SNP and per sample.
#' plotGBSR(gdata, stats = "qtile_ref", q = 0.9)
#'
#' # Close the connection to the GDS file
#' closeGDS(gdata)
#'
plotGBSR  <- function(x,
                      stats = c("dp", "missing", "het"),
                      coord = NULL,
                      q = 0.5,
                      lwd = 0.5,
                      binwidth = NULL,
                      color = c(
                          Marker = "darkblue",
                          Ref = "darkgreen",
                          Het = "magenta",
                          Alt = "blue")) {
    stats_list = c("marker",
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
                   "qtile_ref",
                   "mean_alt",
                   "sd_alt",
                   "qtile_alt",
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
                paste0("Quantile of reference read depth (q=", q, ")"),
                "Mean of alternative read depth",
                "SD of alternative read depth",
                paste0("Quantile of alternative read depth (q=", q, ")"),
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

    df <- .df.maker(x, stats, q, "snp", TRUE)
    p <- ggplot(df)
    p <- .plot.maker(p, stats, binwidth, lwd, coord)
    p <- p + xlab("Physical position (Mb)") +
        ylab("") +
        labs(title = .lab.maker(stats, q)) +
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
#' @param target Either or both of "snp" and "scan", e.g.
#' `target = "snp"` to draw a histogram only for SNPs.
#' @param q An integer to specify a quantile calculated via
#' [calcReadStats()].
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
#' \itemize{
#' \item{"missing"}{"Proportion of missing genotype calls."},
#' \item{"het"}{"Proportion of heterozygote calls."},
#' \item{"raf"}{"Reference allele frequency."},
#' \item{"dp"}{"Total read counts."},
#' \item{"ad_ref"}{"Reference allele read counts."},
#' \item{"ad_alt"}{"Alternative allele read counts."},
#' \item{"rrf"}{"Reference allele read frequency."},
#' \item{"mean_ref"}{"Mean of reference allele read counts."},
#' \item{"sd_ref"}{"Standard deviation of reference allele read counts."},
#' \item{"qtile_ref"}{"Quantile of reference allele read counts."},
#' \item{"mean_alt"}{"Mean of alternative allele read counts."},
#' \item{"sd_alt"}{"Standard deviation of alternative allele read counts."},
#' \item{"qtile_alt"}{"Quantile of alternative allele read counts."},
#' \item{"mq"}{"Mapping quality."},
#' \item{"fs"}{"Phred-scaled p-value (strand bias)"},
#' \item{"qd"}{"Variant Quality by Depth"},
#' \item{"sor"}{"Symmetric Odds Ratio (strand bias)"},
#' \item{"mqranksum"}{"Alt vs. Ref read mapping qualities"},
#' \item{"readposranksum"}{"Alt vs. Ref read position bias"},
#' \item{"baseqranksum"}{"Alt Vs. Ref base qualities"},
#' }
#'
#' To draw scatter plots for "missing", "het", "raf", you need to run
#' [countGenotype()] first to obtain statistics. Similary, "dp",
#' "ad_ref", "ad_alt", "rrf" requires values obtained via [countRead()].
#' [calcReadStats()] should be executed before drawing line plots of
#' "mean_ref", "sd_ref", "qtile_ref", "mean_alt", "sd_alt", and
#' "qtile_alt". "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' and "baseqranksum" only work with `target = "snp"`, if your data
#' contains those values supplied via SNP calling tools like
#' [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#'
#' @return NULL.
#'
#' @export
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @import graphics
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#'
#' # Summarize genotype count information to be used in `pairsGBSR()`
#' gdata <- countGenotype(gdata)
#'
#' # Draw scatter plots of missing rate vs heterozygosity.
#' pairsGBSR(gdata, stats1 = "missing", stats2 = "het")
#'
#' # Close the connection to the GDS file
#' closeGDS(gdata)
#'
#'
pairsGBSR  <- function(x,
                       stats1 = "dp",
                       stats2 = "missing",
                       target = "snp",
                       q = 0.5,
                       size = 0.5,
                       alpha = 0.8,
                       color = c(Marker = "darkblue", Sample = "darkblue"),
                       fill = c(Marker = "skyblue", Sample = "skyblue"),
                       smooth = FALSE) {
    stats_list = c("dp",
                   "ad_ref",
                   "ad_alt",
                   "missing",
                   "het",
                   "raf",
                   "rrf",
                   "mean_ref",
                   "sd_ref",
                   "qtile_ref",
                   "mean_alt",
                   "sd_alt",
                   "qtile_alt",
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
                paste0("Quantile of reference read depth (q=", q, ")"),
                "Mean of alternative read depth",
                "SD of alternative read depth",
                paste0("Quantile of alternative read depth (q=", q, ")"),
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

    df1 <- .df.maker(x, stats1, q, target)
    df2 <- .df.maker(x, stats2, q, target)
    df <-
        data.frame(val1 = df1$val,
                   val2 = df2$val,
                   target = df1$target)

    val1 <- val2 <- NULL
    p <- ggplot(df)
    p <- .pairs.maker(p, size, alpha)
    p <- p + xlab(.lab.maker(stats1, q)) +
        ylab(.lab.maker(stats2, q)) +
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
.df.maker <- function(x, stats, q, target, pos = FALSE) {
    Ref <- NULL
    Missing <- NULL
    if (stats == "geno") {
        snp <-data.frame(Ref = getCountGenoRef(x, "snp", TRUE, TRUE),
                         Het = getCountGenoHet(x, "snp", TRUE, TRUE),
                         Alt = getCountGenoAlt(x, "snp", TRUE, TRUE),
                         Missing = getCountGenoMissing(x, "snp", TRUE, TRUE),
                         chr = getChromosome(x),
                         pos = getPosition(x) * 10^-6,
                         stringsAsFactors=FALSE)
        if (is.null(snp)) {
            msg <- paste0('No data for the statistic: ', stats)
            stop(msg)
        }
        snp <- pivot_longer(snp, Ref:Missing,
                            names_to="genotype", values_to="val")
        scan <- NULL

    } else if (stats == "marker") {
        snp <- data.frame(chr = getChromosome(x),
                          pos = getPosition(x) * 10^-6,
                          target = "Marker",
                          stringsAsFactors=FALSE)
        scan <- NULL

    } else {
        if ("snp" %in% target) {
            snp <- switch(
                stats,
                "dp" = getCountRead(x, "snp"),
                "missing" = getCountGenoMissing(x, "snp", TRUE, TRUE),
                "het" = getCountGenoHet(x, "snp", TRUE, TRUE),
                "raf" = getCountAlleleRef(x, "snp", TRUE, TRUE),
                "rrf" = getCountReadRef(x, "snp", TRUE, TRUE),
                "ad_ref" = getCountReadRef(x, "snp"),
                "mean_ref" = getMeanReadRef(x, "snp"),
                "sd_ref" = getSDReadRef(x, "snp"),
                "qtile_ref" = getQtileReadRef(x, "snp", q),
                "ad_alt" = getCountReadAlt(x, "snp"),
                "mean_alt" = getMeanReadAlt(x, "snp"),
                "sd_alt" = getSDReadAlt(x, "snp"),
                "qtile_alt" = getQtileReadAlt(x, "snp", q),
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
        if ("scan" %in% target) {
            scan <- switch(
                stats,
                "dp" = getCountRead(x, "scan"),
                "missing" = getCountGenoMissing(x, "scan", TRUE, TRUE),
                "het" = getCountGenoHet(x, "scan", TRUE, TRUE),
                "raf" = getCountAlleleRef(x, "scan", TRUE, TRUE),
                "rrf" = getCountReadRef(x, "scan", TRUE, TRUE),
                "ad_ref" = getCountReadRef(x, "scan"),
                "mean_ref" = getMeanReadRef(x, "scan"),
                "sd_ref" = getSDReadRef(x, "scan"),
                "qtile_ref" = getQtileReadRef(x, "scan", q),
                "ad_alt" = getCountReadAlt(x, "scan"),
                "mean_alt" = getMeanReadAlt(x, "scan"),
                "sd_alt" = getSDReadAlt(x, "scan"),
                "qtile_alt" = getQtileReadAlt(x, "scan", q)
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
.lab.maker <- function(stats, q) {
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
        "qtile_ref" = paste0("Quantile of reference read depth (q=", q, ")"),
        "mean_alt" = "Mean of alternative read depth",
        "sd_alt" = "SD of alternative read depth",
        "qtile_alt" = paste0("Quantile of alternative read depth (q=", q, ")"),
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
        "qtile_ref" = c(0, NA),
        "mean_alt" = c(0, NA),
        "sd_alt" = c(0, NA),
        "qtile_alt" = c(0, NA),
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
.boxplot.maker <- function(p, stats) {
    val <- NULL
    target <- NULL
    p <- p + geom_histogram(aes(x=val, color=target)) +
        facet_wrap( ~ target, scales="free")
    return(p)
}

# Internal function to draw a line plot.
.plot.maker <- function(p, stats, binwidth, lwd, coord) {
    val <- genotype <- target <- pos <- NULL
    if ("geno" %in% stats) {
        p <-p +
            geom_line(aes(x=pos, y=val, color=genotype, group=genotype),
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
#' @param node Either one of "raw", "filt", and "cor" to output raw
#' genotype data, filtered genotype data, or corrected genotype data,
#' respectively.
#' @param dot_fill A string to indicate the dot color in the plot.
#' @param line_color A string to indicate the line color in the plot.
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#'
#' plotDosage(gdata, ind = 1)
#'
#' # Close the connection to the GDS file
#' closeGDS(gdata)
#'
#' @return Invisibly returns a ggplot object.
#' @export
#' @import ggplot2
plotDosage <- function(x,
                       coord = NULL,
                       chr = NULL,
                       ind = 1,
                       node = "raw",
                       dot_fill = "green",
                       line_color = "magenta") {
    pos <- na <- NULL

    if (is.null(chr)) {
        chr <- rep(TRUE, nscan(x))
    } else {
        chr <- getChromosome(x) %in% chr
    }

    geno <- getGenotype(x, node=node, parents=TRUE)[ind, chr]
    id <- getScanID(x, valid = FALSE)[ind]
    df <- data.frame(chr = getChromosome(x)[chr],
                     pos = getPosition(x)[chr],
                     geno = geno)
    df <- df[!is.na(df$geno), ]

    p <- ggplot(df, aes( x = pos * 10^-6, y=geno, group=chr)) +
        geom_point(size=0.8, stroke=0, color=dot_fill) +
        geom_line(color=line_color) +
        labs(title=paste0("Reference allele dosage: ", id)) +
        ylim(0, max(getPloidy(x))) +
        xlab("Physical position (Mb)") +
        ylab("Reference allele dosage")

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
#' @param ind An index to specify samples to be drawn.
#' @param node Either one of "raw", "filt", and "cor" to output raw
#' genotype data, filtered genotype data, or corrected genotype data,
#' respectively.
#' @param dot_fill A string to indicate the dot color in a plot.
#'
#'
#' @examples
#' # Load data in the GDS file and instantiate a [GbsrGenotypeData] object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#'
#' plotReadRatio(gdata, ind = 1)
#'
#' # Close the connection to the GDS file
#' closeGDS(gdata)
#'
#' @return Invisibly returns a ggplot object.
#' @export
#' @import ggplot2
plotReadRatio <- function(x,
                          coord = NULL,
                          chr = NULL,
                          ind = 1,
                          node = "raw",
                          dot_fill = "blue") {
    pos <- ad <- NULL
    validmarker <- getValidSnp(x)

    if (is.null(chr)) {
        chr <- rep(TRUE, nscan(x))
    } else {
        chr <- getChromosome(x) %in% chr
    }

    id <- getScanID(x, FALSE)[ind]
    read <- getRead(x, node=node, parents=TRUE)
    ref <- read$ref[ind, chr]
    alt <- read$alt[ind, chr]

    df <- data.frame(chr = getChromosome(x)[validmarker],
                     pos = getPosition(x)[validmarker],
                     ad = ref / (ref + alt))

    p <- ggplot(df, aes(x = pos * 10^-6, y = ad, group = chr, color = chr)) +
        geom_point(size=0.8, alpha=0.2, stroke=0) +
        labs(title=paste0("Reference allele read ratio: ", id)) +
        ylim(0, 1) +
        xlab("Position (Mb)") +
        ylab("Reference allele read ratio")

    p <- p + facet_wrap(~ chr, coord[1], coord[2], "free_x",
                        dir="v", strip.position="right") +
        theme(legend.position="none")
    return(p)
}
