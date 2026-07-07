object <- loadGDS("~/hdd5/ngsData/gbs/20230309_yamamotoTP/mcptaggr/genotype.gds", ploidy = 4, load_filter = TRUE)

recomb_rate <- 0.04
error_rate <- 0.0025
call_threshold <- 0.9
het_parent <- FALSE
optim <- TRUE
iter <- 4
n_threads <- 5
dummy_reads <- 10
fix_bias <- NULL
fix_mismap <- NULL
chr_i <- "chr01"

object <- initScheme(object = object)
object <- addScheme(object = object, crosstype = "self")
object <- addScheme(object = object, crosstype = "self")
object <- addScheme(object = object, crosstype = "self")
object <- addScheme(object = object, crosstype = "self")
object <- addScheme(object = object, crosstype = "self")
object <- addScheme(object = object, crosstype = "self")
object <- addScheme(object = object, crosstype = "self")
object <- addScheme(object = object, crosstype = "self")

closeGDS(object = object)
w1 <- 0.5
ploidy <- 4
ref = 10
alt = 1
for(g in 0:4){
    denomi = w1 * (ploidy - g) + (1 - w1) * g
    prob <- w1 * (ploidy - g) / denomi
    # print(prob^ref * (1 - prob)^alt)
    print(prob)
}

calcGenoprobR(ref = 10, alt = 5, eseq0 = 0.9975, eseq1 = 0.0025, w1 = 0.5, het = TRUE, ploidy = 4)

c(4.9414e-69,
  3.2589e-10,
  0.00531277,
  0.994687,
  5.16783e-17) + 0.0005/
sum(c(4.9414e-69,
      3.2589e-10,
      0.00531277,
      0.994687,
      5.16783e-17)) + 0.0005 * 5
