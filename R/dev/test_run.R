source("~/11_wd2/gbscleanr/v1/ForManuscript/script/compareGeno.R")
devtools::load_all(path = "~/01_wd/softDevel/SimPop/")
devtools::load_all(path = "~/01_wd/softDevel/GBScleanR/")

dirs <- grep("simpop", list.dirs("~/11_wd2/gbscleanr/v1/ForManuscript", recursive = F), value = T)

# homoP2 F2
rec1 <- NULL
vcf_files <- list.files(dirs[3], ".vcf", full.names = TRUE)
vcf_files <- grep("_LB", vcf_files, invert = TRUE, value = TRUE)
gds_files <- list.files(dirs[3], "_gbsr.gds", full.names = TRUE)
gbsrVCF2GDS(vcf_fn = vcf_files[51], out_fn = gds_files[51], force = TRUE)
gds <- loadGDS(gds_files[51])
parents <- grep("Founder", getSamID(gds), value = TRUE)
gds <- setParents(gds, parents)
gds <- initScheme(gds, matrix(1:2, 2))
gds <- addScheme(gds, "self")
runtime1 <- system.time({gds <- estGeno(object = gds,
                                        iter = 4,
                                        het_parent = FALSE,
                                        n_threads = 50)
})
parents <- grep("Founder", getSamID(gds), value = TRUE)
geno <- getGenotype(gds, node = "cor")
pop_file <- sub("_gbsr.gds", "_SimPop.Rdata", gds_files[51])
load(pop_file)
true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
true <- t(sapply(true, colSums))
df <- compareGeno(geno, true)
apply(df$ind, 2, mean, na.rm = TRUE)
# correct     miscall     missing
# 0.983185484 0.003167742 0.013646774
closeGDS(gds)
print(runtime1)
# user  system elapsed
# 220.487  25.499  13.710

# hetP2 F1
rec2 <- NULL
vcf_files <- list.files(dirs[1], ".vcf", full.names = TRUE)
gds_files <- list.files(dirs[1], "_gbsr.gds", full.names = TRUE)
gbsrVCF2GDS(vcf_fn = vcf_files[51], out_fn = gds_files[51], force = TRUE)
gds <- loadGDS(gds_files[51])
parents <- grep("Founder", getSamID(gds), value = TRUE)
gds <- setParents(gds, parents[1])
gds <- initScheme(gds, matrix(1:2, 2))
runtime2 <- system.time({gds <- estGeno(object = gds,
                                        iter = 4,
                                        het_parent = TRUE,
                                        n_threads = 50)
})
parents <- grep("Founder", getSamID(gds), value = TRUE)
geno <- getGenotype(gds, node = "cor")
pop_file <- sub("_gbsr.gds", "_SimPop.Rdata", gds_files[51])
load(pop_file)
true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
true <- t(sapply(true, colSums))
df <- compareGeno(geno, true)
apply(df$ind, 2, mean, na.rm = TRUE)
# correct    miscall    missing
# 0.92197419 0.06662097 0.01140484
closeGDS(gds)
print(runtime2)
# user   system  elapsed
# 1688.315   18.408   48.784

# hetP8 RIL
vcf_files <- list.files(dirs[5], ".vcf", full.names = TRUE)
gds_files <- list.files(dirs[5], "_gbsr.gds", full.names = TRUE)
gbsrVCF2GDS(vcf_fn = vcf_files[51], out_fn = gds_files[51], force = TRUE)
gds <- loadGDS(gds_files[51])
parents <- grep("Founder", getSamID(gds), value = TRUE)
gds <- setParents(gds, parents)
gds <- initScheme(gds, matrix(1:8, 2))
gds <- addScheme(gds, "pairing", matrix(9:12, 2))
gds <- addScheme(gds, "pairing", matrix(13:14, 2))
gds <- addScheme(gds, "self")
gds <- addScheme(gds, "self")
gds <- addScheme(gds, "self")
gds <- addScheme(gds, "self")
gds <- addScheme(gds, "self")
runtime3 <- system.time({gds <- estGeno(object = gds,
                                        iter = 1,
                                        het_parent = FALSE,
                                        n_threads = 50)
})
parents <- grep("Founder", getSamID(gds), value = TRUE)
geno <- getGenotype(gds, node = "cor")
pop_file <- sub("_gbsr.gds", "_SimPop.Rdata", gds_files[51])
load(pop_file)
true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
true <- t(sapply(true, colSums))
df <- compareGeno(geno, true)
apply(df$ind, 2, mean, na.rm = TRUE)
# correct    miscall    missing
# 0.97278548 0.01157903 0.01563548
closeGDS(gds)
print(runtime3)
# user  system elapsed
# 382.115  35.309  15.115
