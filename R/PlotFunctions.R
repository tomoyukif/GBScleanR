#' Draw histograms of specified statistics
#'
#' @param x A GbsrGenotypeData object.
#' @param stats A vector of strings to specify statistics to be drawn.
#' @param target Either or both of "snp" and "scan", e.g. `target = "snp"` to draw a histogram only for SNPs.
#' @param q An integer to specify a quantile calculated via [calcReadStats()].
#' @param binwidth An integer to specify bin width of the histogram. This value is passed to the ggplot function.
#' @param color A named vector "Marker" and "Sample" to specify border color of bins in the histograms.
#' @param fill A named vector "Marker" and "Sample" to specify fill color of bins in the histograms.
#' @param ggargs An expression to be internally passed to a ggplot object to be drawn.
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
#' To draw histograms for "missing", "het", "raf", you need to run [countGenotype()]
#' first to obtain statistics. Similary, "dp", "ad_ref", "ad_alt", "rrf" requires
#' values obtained via [countRead()]. [calcReadStats()] should be executed before
#' drawing histograms of "mean_ref", "sd_ref", "qtile_ref", "mean_alt", "sd_alt",
#' and "qtile_alt". "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' and "baseqranksum" only work with `target = "snp"`, if your data contains those
#' values supplied via SNP calling tools like [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @export
#'
#' @return NULL.
#'
#' @import ggplot2
#' @import graphics
#'
#' @examples
#' # Draw histograms of missing rate, heterozygosity, and reference
#' # allele frequency per SNP and per sample.
#' gds_fn <- system.file("extdata", "simpop.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#' gdata <- countGenotype(gdata)
#' histGBSR(gdata, stats = c("missing", "het", "raf"))
#'
#' # Draw histograms of 90 percentile values of reference read counts
#' # and alternative read counts per SNP and per sample.
#' gdata <- calcReadStats(gdata, q = 0.9)
#' histGBSR(gdata, stats = c("qtile_ref", "qtile_alt"), q = 0.9)
#'
histGBSR  <- function(x,
                      stats = c("dp", "missing", "het"),
                      target = c("snp", "scan"),
                      q = 0.5,
                      binwidth = NULL,
                      color = c(Marker = "darkblue", Sample = "darkblue"),
                      fill = c(Marker = "skyblue", Sample = "skyblue"),
                      ggargs = NULL) {
  stats_list = c(
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
  if ("all" %in% stats) {
    stats <- stats_list
  } else {
    stats <- match.arg(arg = stats,
                       choices = stats_list,
                       several.ok = TRUE)
  }
  
  for (i in stats) {
    if (is.null(binwidth)) {
      if (i %in% c("missing", "het", "raf", "rrf")) {
        binwidth <- 0.025
      }
    }
    
    df <- .df.maker(x, i, q, target)
    p <- ggplot(df)
    p <- .hist.maker(p, i, binwidth)
    p <- p + xlab(.lab.maker(i, q)) +
      ylab("Count") +
      xlim(.limit.maker(i)) +
      scale_fill_manual(values = fill, breaks = names(fill)) +
      scale_color_manual(values = color, breaks = names(color)) +
      theme(
        axis.title = element_text(face = "bold", size = rel(1.3)),
        axis.text = element_text(size = rel(1.2)),
        strip.text = element_text(size = rel(1.3), face = "bold"),
        legend.position = "none"
      )
    if (!is.null(ggargs)) {
      p <- eval(parse(text = paste0("p + ", ggargs)))
    }
    print(p)
  }
}

#' Draw boxplots of specified statistics
#'
#' @param x A GbsrGenotypeData object.
#' @param stats A vector of strings to specify statistics to be drawn.
#' @param target Either or both of "snp" and "scan", e.g. `target = "snp"` to draw a histogram only for SNPs.
#' @param q An integer to specify a quantile calculated via [calcReadStats()].
#' @param color A named vector "Marker" and "Sample" to specify border color of bins in the histograms.
#' @param fill A named vector "Marker" and "Sample" to specify fill color of bins in the histograms.
#' @param ggargs An expression to be internally passed to a ggplot object to be drawn.
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
#' To draw boxplots for "missing", "het", "raf", you need to run [countGenotype()]
#' first to obtain statistics. Similary, "dp", "ad_ref", "ad_alt", "rrf" requires
#' values obtained via [countRead()]. [calcReadStats()] should be executed before
#' drawing boxplots of "mean_ref", "sd_ref", "qtile_ref", "mean_alt", "sd_alt",
#' and "qtile_alt". "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' and "baseqranksum" only work with `target = "snp"`, if your data contains those
#' values supplied via SNP calling tools like [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @export
#'
#' @return NULL.
#'
#' @import ggplot2
#' @import graphics
#'
#' @examples
#' # Draw boxplots of missing rate, heterozygosity, and reference
#' # allele frequency per SNP and per sample.
#' gds_fn <- system.file("extdata", "simpop.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#' gdata <- countGenotype(gdata)
#' boxplotGBSR(gdata, stats = c("missing", "het", "raf"))
#'
#' # Draw boxplots of 90 percentile values of reference read counts and
#' # alternative read counts per SNP and per sample.
#' gdata <- calcReadStats(gdata, q = 0.9)
#' boxplotGBSR(gdata, stats = c("qtile_ref", "qtile_alt"), q = 0.9)
#'
boxplotGBSR <- function(x,
                        stats = "missing",
                        target = c("snp", "scan"),
                        q = 0.5,
                        color = c(Marker = "darkblue", Sample = "darkblue"),
                        fill = c(Marker = "skyblue", Sample = "skyblue"),
                        ggargs = NULL) {
  stats_list = c(
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
  if ("all" %in% stats) {
    stats <- stats_list
  } else {
    stats <- match.arg(arg = stats,
                       choices = stats_list,
                       several.ok = TRUE)
  }
  
  for (i in stats) {
    df <- .df.maker(x, i, q, target)
    p <- ggplot(df)
    p <- .boxplot.maker(p, i)
    p <- p + ylab("") +
      xlab(.lab.maker(i, q)) +
      xlim(.limit.maker(i)) +
      scale_fill_manual(values = fill, breaks = names(fill)) +
      scale_color_manual(values = color, breaks = names(color)) +
      theme(
        axis.title = element_text(face = "bold", size = rel(1.3)),
        axis.text = element_text(size = rel(1.2)),
        strip.text = element_text(size = rel(1.3), face = "bold"),
        legend.position = "none"
      )
    if (!is.null(ggargs)) {
      p <- eval(parse(text = paste0("p + ", ggargs)))
    }
    print(p)
  }
}


#' Draw line plots of specified statistics
#'
#' @param x A GbsrGenotypeData object.
#' @param stats A vector of strings to specify statistics to be drawn.
#' @param coord A vector with two integer specifying the number of rows and columns to draw faceted line plots for chromosomes.
#' @param q An integer to specify a quantile calculated via [calcReadStats()].
#' @param lwd A numeric value to specify the line width in plots.
#' @param binwidth An integer to specify bin width of the histogram. This argument only work with `stats = "marker"` and is passed to the ggplot function.
#' @param color A strings vector named "Marker", "Ref", "Het", "Alt" to specify line colors. `stats = "geno` only requires "Ref", "Het" and "Alt", while others uses the value named "Marker".
#' @param ggargs An expression to be internally passed to a ggplot object to be drawn.
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
#' To draw line plots for "missing", "het", "raf", you need to run [countGenotype()]
#' first to obtain statistics. Similary, "dp", "ad_ref", "ad_alt", "rrf" requires
#' values obtained via [countRead()]. [calcReadStats()] should be executed before
#' drawing line plots of "mean_ref", "sd_ref", "qtile_ref", "mean_alt", "sd_alt",
#' and "qtile_alt". "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' and "baseqranksum" only work with `target = "snp"`, if your data contains those
#' values supplied via SNP calling tools like [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#' @export
#' @import ggplot2
#' @import graphics
#'
#' @return NULL.
#'
#' @examples
#' # Draw line plots of missing rate, heterozygosity, proportion of genotype calls per SNP.
#' gds_fn <- system.file("extdata", "simpop.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#' gdata <- countGenotype(gdata)
#' plotGBSR(gdata, stats = c("missing", "het", "raf", "geno"))
#'
#' # Draw line plots of 90 percentile values of reference read counts and
#' # alternative read counts per SNP and per sample.
#' gdata <- calcReadStats(gdata, q = 0.9)
#' plotGBSR(gdata, stats = c("qtile_ref", "qtile_alt"), q = 0.9)
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
                        Alt = "blue"
                      ),
                      ggargs = NULL) {
  stats_list = c(
    "marker",
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
  if ("all" %in% stats) {
    stats <- stats_list
  } else {
    stats <- match.arg(arg = stats,
                       choices = stats_list,
                       several.ok = TRUE)
  }
  
  for (i in stats) {
    df <- .df.maker(x, i, q, target = "snp", TRUE)
    p <- ggplot(df)
    p <- .plot.maker(p, i, binwidth, lwd, coord)
    p <- p + xlab("Physical position (Mb)") +
      ylab("") +
      labs(title = .lab.maker(i, q)) +
      ylim(.limit.maker(i)) +
      scale_x_continuous(expand = expansion(0, 0.5), limits = c(0, NA)) +
      scale_color_manual(values = color, breaks = names(color)) +
      theme(
        plot.title = element_text(face = "bold", size = rel(1.5)),
        axis.title = element_text(face = "bold", size = rel(1.3)),
        axis.text = element_text(size = rel(1)),
        strip.text.y = element_text(
          size = rel(1.2),
          face = "bold",
          angle = 0
        )
      )
    if (i != "geno") {
      p <- p + theme(legend.position = "none")
    }
    if (!is.null(ggargs)) {
      p <- eval(parse(text = paste0("p + ", ggargs)))
    }
    print(p)
  }
}



#' Draw a scatter plot of a pair of specified statistics
#'
#' @param x A GbsrGenotypeData object.
#' @param stats1 A string to specify statistics to be drawn.
#' @param stats2 A string to specify statistics to be drawn.
#' @param target Either or both of "snp" and "scan", e.g. `target = "snp"` to draw a histogram only for SNPs.
#' @param q An integer to specify a quantile calculated via [calcReadStats()].
#' @param size A numeric value to specify the dot size of a scatter plot.
#' @param alpha A numeric value \[0-1\] to specify the transparency of dots in a scatter plot.
#' @param color A named vector "Marker" and "Sample" to specify border color of bins in the histograms.
#' @param fill A named vector "Marker" and "Sample" to specify fill color of bins in the histograms.`stats = "geno` only requires "Ref", "Het" and "Alt", while others uses the value named "Marker".
#' @param smooth A logical value to indicate whether draw a smooth line for data points. See also [ggplot2::stat_smooth()].
#' @param ggargs An expression to be internally passed to a ggplot object to be drawn.
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
#' To draw scatter plots for "missing", "het", "raf", you need to run [countGenotype()]
#' first to obtain statistics. Similary, "dp", "ad_ref", "ad_alt", "rrf" requires
#' values obtained via [countRead()]. [calcReadStats()] should be executed before
#' drawing line plots of "mean_ref", "sd_ref", "qtile_ref", "mean_alt", "sd_alt",
#' and "qtile_alt". "mq", "fs", "qd", "sor", "mqranksum", "readposranksum",
#' and "baseqranksum" only work with `target = "snp"`, if your data contains those
#' values supplied via SNP calling tools like [GATK](https://gatk.broadinstitute.org/hc/en-us).
#'
#'
#' @return NULL.
#'
#' @export
#' @import ggplot2
#' @import graphics
#'
#' @examples
#' # Draw scatter plots of missing rate vs heterozygosity.
#' gds_fn <- system.file("extdata", "simpop.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#' gdata <- countGenotype(gdata)
#' plotGBSR(gdata, stats1 = "missing", stats2 = "het")
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
                       smooth = FALSE,
                       ggargs = NULL) {
  stats_list = c(
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
    match.arg(arg = stats1,
              choices = stats_list,
              several.ok = TRUE)
  stats2 <-
    match.arg(arg = stats2,
              choices = stats_list,
              several.ok = TRUE)
  
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
    scale_fill_manual(values = fill, breaks = names(fill)) +
    scale_color_manual(values = color, breaks = names(color)) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.5)),
      axis.title = element_text(face = "bold", size = rel(1.3)),
      axis.text = element_text(size = rel(1)),
      strip.text.y = element_text(
        size = rel(1.2),
        face = "bold",
        angle = 0
      ),
      legend.position = "none"
    )
  if (smooth) {
    p <- p + stat_smooth(aes(x = val1, y = val2),)
  }
  if (!is.null(ggargs)) {
    p <- eval(parse(text = paste0("p + ", ggargs)))
  }
  print(p)
}


# Internal function to build a data.frame passed to ggplot().
.df.maker <- function(x, stats, q, target, pos = FALSE) {
  Ref <- NULL
  Missing <- NULL
  if (stats == "geno") {
    snp <-
      data.frame(
        Ref = getCountGenoRef(x, target = "snp", prop = TRUE),
        Het = getCountGenoHet(x, target = "snp", prop = TRUE),
        Alt = getCountGenoAlt(x, target = "snp", prop = TRUE),
        Missing = getCountGenoMissing(x, target = "snp", prop =
                                        TRUE),
        chr = getChromosome(x),
        pos = getPosition(x) * 10 ^ -6,
        stringsAsFactors = FALSE
      )
    if (is.null(snp)) {
      msg <- paste0('No data for the statistic: ', stats)
      stop(msg)
    }
    snp <-
      tidyr::pivot_longer(
        snp,
        cols = Ref:Missing,
        names_to = "genotype",
        values_to = "val"
      )
    scan <- NULL
    
  } else if (stats == "marker") {
    snp <- data.frame(
      chr = getChromosome(x),
      pos = getPosition(x) * 10 ^ -6,
      target = "Marker",
      stringsAsFactors = FALSE
    )
    scan <- NULL
    
  } else {
    if ("snp" %in% target) {
      snp <- switch(
        stats,
        "dp" = getCountRead(x, target = "snp"),
        "missing" = getCountGenoMissing(x, target = "snp", prop = TRUE),
        "het" = getCountGenoHet(x, target = "snp", prop = TRUE),
        "raf" = getCountAlleleRef(x, target = "snp", prop = TRUE),
        "rrf" = getCountReadRef(x, target = "snp", prop = TRUE),
        "ad_ref" = getCountReadRef(x, target = "snp"),
        "mean_ref" = getMeanReadRef(x, target = "snp"),
        "sd_ref" = getSDReadRef(x, target = "snp"),
        "qtile_ref" = getQtileReadRef(x, target = "snp", q = q),
        "ad_alt" = getCountReadAlt(x, target = "snp"),
        "mean_alt" = getMeanReadAlt(x, target = "snp"),
        "sd_alt" = getSDReadAlt(x, target = "snp"),
        "qtile_alt" = getQtileReadAlt(x, target = "snp", q = q),
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
          snp <- data.frame(
            val = snp,
            target = "Marker",
            chr = getChromosome(x),
            pos = getPosition(x) * 10 ^ -6,
            stringsAsFactors = FALSE
          )
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
        "dp" = getCountRead(x, target = "scan"),
        "missing" = getCountGenoMissing(x, target = "scan", prop = TRUE),
        "het" = getCountGenoHet(x, target = "scan", prop = TRUE),
        "raf" = getCountAlleleRef(x, target = "scan", prop = TRUE),
        "rrf" = getCountReadRef(x, target = "scan", prop = TRUE),
        "ad_ref" = getCountReadRef(x, target = "scan"),
        "mean_ref" = getMeanReadRef(x, target = "scan"),
        "sd_ref" = getSDReadRef(x, target = "scan"),
        "qtile_ref" = getQtileReadRef(x, target = "scan", q = q),
        "ad_alt" = getCountReadAlt(x, target = "scan"),
        "mean_alt" = getMeanReadAlt(x, target = "scan"),
        "sd_alt" = getSDReadAlt(x, target = "scan"),
        "qtile_alt" = getQtileReadAlt(x, target = "scan", q = q)
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
    "mqranksum" = c(NA, NA),
    "readposranksum" = c(NA, NA),
    "baseqranksum" = c(NA, NA)
  )
  return(lim)
}

# Internal function to draw a histogram.
.hist.maker <- function(p, stats, binwidth, color, fill) {
  val <- target <- NULL
  p <-
    p + geom_histogram(
      mapping = aes(x = val, color = target, fill = target),
      binwidth = binwidth,
      boundary = 0
    ) +
    facet_wrap( ~ target, scales = "free")
  return(p)
}

# Internal function to draw a boxplot.
.boxplot.maker <- function(p, stats) {
  val <- NULL
  target <- NULL
  p <- p + geom_histogram(mapping = aes(x = val, color = target)) +
    facet_wrap( ~ target, scales = "free")
  return(p)
}

# Internal function to draw a line plot.
.plot.maker <- function(p, stats, binwidth, lwd, coord) {
  val <- genotype <- target <- pos <- NULL
  if ("geno" %in% stats) {
    p <-
      p + geom_line(
        mapping = aes(
          x = pos,
          y = val,
          color = genotype,
          group = genotype
        ),
        size = lwd
      ) +
      facet_wrap(
        facets = ~ chr,
        nrow = coord[1],
        ncol = coord[2],
        scales = "free_x",
        dir = "v",
        strip.position = "right"
      )
  } else if ("marker" %in% stats) {
    p <-
      p + geom_histogram(
        mapping = aes(
          x = pos,
          fill = target,
          color = target
        ),
        binwidth = binwidth,
        boundary = 0
      ) +
      facet_wrap(
        facets = ~ chr,
        nrow = coord[1],
        ncol = coord[2],
        scales = "free_x",
        dir = "v",
        strip.position = "right"
      )
  } else {
    p <-
      p + geom_line(mapping = aes(x = pos, y = val, color = target),
                    size = lwd) +
      facet_wrap(
        facets = ~ chr,
        nrow = coord[1],
        ncol = coord[2],
        scales = "free_x",
        dir = "v",
        strip.position = "right"
      )
  }
  return(p)
}

# Internal function to draw a scatter plot.
.pairs.maker <- function(p, size, alpha) {
  val1 <- val2 <- target <- NULL
  p <-
    p + geom_point(
      mapping = aes(
        x = val1,
        y = val2,
        color = target,
        fill = target
      ),
      size = size,
      alpha = alpha
    ) +
    facet_wrap( ~ target, scales = "free")
  return(p)
}

#' Draw line plots of allele dosage per marker per sample.
#'
#' This function counts a reference allele dosage per marker per sample and
#' draw line plots of them in facets for each chromosome for each sample.
#'
#' @param x A GbsrGenotypeData object.
#' @param coord A vector with two integer specifying the number of rows and columns to draw faceted line plots for chromosomes.
#' @param chr A vector of indexes to specify chromosomes to be drawn.
#' @param ind A vector of indexes to specify samples to be drawn.
#' @param valid_only A logical value whether to draw a plot only for valid markers and samples. You can get validity with `getValidSnp()` and `getValidScan`.
#' @param dot_fill A string to indicate the dot color in a plot.
#'
#' @examples
#' gds_fn <- system.file("extdata", "simpop.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#' gdata <- countGenotype(gdata)
#' plotDosage(gdata, ind = 1)
#'
#' @return NULL.
#' @export
#' @import ggplot2
plotDosage <- function(x,
                       coord = NULL,
                       chr = NULL,
                       ind = NULL,
                       valid_only = TRUE,
                       dot_fill = "blue") {
  pos <- na <- NULL
  validscan <- getValidScan(x)
  validmarker <- getValidSnp(x)
  if (is.null(ind)) {
    ind <- seq_len(validscan)
  }
  
  if (is.null(chr)) {
    chr <- rep(TRUE, sum(validmarker))
  } else {
    chr <- getChromosome(x) %in% chr
  }
  sel_marker <- validmarker & chr
  
  parents <- x@scanAnnot$parents[x@scanAnnot$parents != 0]
  if (!is.null(parents)) {
    ind <- c(ind[ind %in% parents], ind[!ind %in% parents])
  }
  
  genotype_node <-
    gdsfmt::index.gdsn(node = x@data@handler,
                       path = x@data@genotypeVar)
  
  for (i in ind) {
    if (valid_only & !validscan[i]) {
      next()
    }
    id <- getScanID(x, FALSE)[i]
    sel <- list(seq_len(nscan(x)) %in% i, sel_marker)
    geno <- gdsfmt::readex.gdsn(genotype_node, sel = sel)
    
    df <- data.frame(chr = getChromosome(x)[sel_marker],
                     pos = getPosition(x)[sel_marker],
                     geno = geno)
    df$na <- is.na(geno)
    
    p <- ggplot(
      data = df,
      mapping = aes(
        x = pos * 10 ^ -6,
        y = geno,
        group = chr,
        color = chr,
        fill = na,
        alpha = na
      )
    ) +
      geom_point() +
      geom_line() +
      labs(title = paste0("Reference allele dosage: ", id)) +
      xlab("Physical position (Mb)") +
      ylab("Reference allele dosage") +
      scale_alpha_manual(values = c(0.3, 0.5),
                         breaks = c(TRUE, FALSE)) +
      scale_fill_manual(values = c("darkgray", dot_fill),
                        breaks = c(TRUE, FALSE))
    
    p <- p +
      facet_wrap(
        facets = ~ chr,
        nrow = coord[1],
        ncol = coord[2],
        scales = "free_x",
        dir = "v",
        strip.position = "right"
      )
    print(p)
  }
}

#' Draw line plots of proportion of reference allele read counts per marker per sample.
#'
#' This function calculate a proportion of reference allele read counts per marker per sample and
#' draw line plots of them in facets for each chromosome for each sample.
#'
#' @param x A GbsrGenotypeData object.
#' @param coord A vector with two integer specifying the number of rows and columns to draw faceted line plots for chromosomes.
#' @param chr A vector of indexes to specify chromosomes to be drawn.
#' @param ind A vector of indexes to specify samples to be drawn.
#' @param valid_only A logical value whether to draw a plot only for valid markers and samples. You can get validity with `getValidSnp()` and `getValidScan`.
#' @param dot_fill A string to indicate the dot color in a plot.
#'
#'
#' @examples
#' gds_fn <- system.file("extdata", "simpop.gds", package = "GBScleanR")
#' gdata <- loadGDS(gds_fn)
#' gdata <- countGenotype(gdata)
#' plotReadRatio(gdata, ind = 1)
#'
#' @return NULL.
#' @export
#' @import ggplot2
plotReadRatio <- function(x,
                          coord = NULL,
                          chr = NULL,
                          ind = NULL,
                          valid_only = TRUE,
                          dot_fill = "blue") {
  pos <- NULL
  validscan <- getValidScan(x)
  validmarker <- getValidSnp(x)
  if (is.null(ind)) {
    ind <- seq_len(validscan)
  }
  
  if (is.null(chr)) {
    chr <- rep(TRUE, sum(validmarker))
  } else {
    chr <- getChromosome(x) %in% chr
  }
  validmarker <- validmarker & chr
  sel_marker <- rep(validmarker, each = 2)
  
  parents <- x@scanAnnot$parents[x@scanAnnot$parents != 0]
  if (!is.null(parents)) {
    ind <- c(ind[ind %in% parents], ind[!ind %in% parents])
  }
  
  if (x@data@genotypeVar == "genotype.filt") {
    ad_data_node <-
      gdsfmt::index.gdsn(node = x@data@handler, path = "/annotation/format/AD/filt.data")
  } else {
    ad_data_node <-
      gdsfmt::index.gdsn(node = x@data@handler, path = "/annotation/format/AD/data")
  }
  
  for (i in ind) {
    if (valid_only & !validscan[i]) {
      next()
    }
    id <- getScanID(x, FALSE)[i]
    sel <- list(seq_len(nscan(x)) %in% i, sel_marker)
    ad <- gdsfmt::readex.gdsn(ad_data_node, sel = sel)
    ref <- ad[c(TRUE, FALSE)]
    alt <- ad[c(FALSE, TRUE)]
    
    if (!is.null(x@snpAnnot$flipped)) {
      flipped <- ref[x@snpAnnot$flipped]
      ref[x@snpAnnot$flipped] <- alt[x@snpAnnot$flipped]
      alt[x@snpAnnot$flipped] <- flipped
    }
    
    df <- data.frame(chr = getChromosome(x)[validmarker],
                     pos = getPosition(x)[validmarker],
                     ad = ref / (ref + alt))
    
    p <- ggplot(data = df,
                mapping = aes(
                  x = pos * 10 ^ -6,
                  y = ad,
                  group = chr,
                  color = chr
                )) +
      geom_point(alpha = 0.2) +
      # geom_line() +
      labs(title = paste0("Reference allele read ratio: ", id)) +
      ylim(0, 1) +
      xlab("Position (Mb)") +
      ylab("Reference allele read ratio")
    
    p <- p +
      facet_wrap(
        facets = ~ chr,
        nrow = coord[1],
        ncol = coord[2],
        scales = "free_x",
        dir = "v",
        strip.position = "right"
      )
    print(p)
  }
}
