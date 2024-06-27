setwd("/home/wolffha/CC-GWAS/")

library(ggplot2, lib = "/home/math/wolffha/R/myLibs/")
library(tidyr, lib = "/home/math/wolffha/R/myLibs/")
library(readr)
library(tidyverse, lib = "/home/math/wolffha/R/myLibs/")
library(PhenotypeSimulator, lib = "/home/math/wolffha/R/myLibs/")
library(CaseControlAF, lib = "/home/math/wolffha/R/myLibs/")
library(DescTools, lib = "/home/math/wolffha/R/myLibs/")

source("simulation_source.R")
source("AFAnalysis_server.R")

# Want to do the following:
# Sample size:  1000, 10000, 100000
# Number Covariates: 0, 1, 3

# first all three sample sizes with 3 covariates
print("Doing 100,000 Inds, 3 confounders")
test <- createSimulation_noPopStruc(seed = 2020, NInds = 100000, prop_case = 1, NSNPs = 10000,
                         NPhenos = 1, NConfounders = 3, confounder_dist = c("bin", "cat_norm", "norm"),
                         cat_confounders = c(5), prop_independent = c(0, 1, 0), prob_binary = c(0.5),
                         var_geno = 0.6, var_noise = 0.1, var_confounders = 0.9, totalSNPeffect = 0.01,
                         prop_shared = 0.6, mConfounders = c(30), sdConfounders = c(20), N_causalSNPs = 100)
p <- runAnalysis_AFComps(test, title = "100,000 individuals | 3 covariates | 10,000 SNPs (100 causal)")
ggsave(plot = p, filename = "100000obs_3covariate.png", device = "png",
       height = 9, width = 12, units = "in")
saveRDS(test, "simDat_100000obs_3covariate.rds")

print("Doing 10,000 Inds, 3 confounders")
test <- createSimulation_noPopStruc(seed = 2021, NInds = 10000, prop_case = 1, NSNPs = 10000,
                                    NPhenos = 1, NConfounders = 3, 
                                    confounder_dist = c("bin", "cat_norm", "norm"),
                                    cat_confounders = c(5), prop_independent = c(0, 1, 0), prob_binary = c(0.5),
                                    var_geno = 0.6, var_noise = 0.1, var_confounders = 0.9, 
                                    totalSNPeffect = 0.01, prop_shared = 0.6, mConfounders = c(30),
                                    sdConfounders = c(20), N_causalSNPs = 100)
p <- runAnalysis_AFComps(test, title = "10,000 individuals | 3 covariates | 10,000 SNPs (100 causal)")
ggsave(plot = p, filename = "10000obs_3covariate.png", device = "png",
       height = 9, width = 12, units = "in")
saveRDS(test, "simDat_10000obs_3covariate.rds")

print("Doing 1,000 Inds, 3 confounders")
test <- createSimulation_noPopStruc(seed = 2022, NInds = 1000, prop_case = 1, NSNPs = 10000,
                                    NPhenos = 1, NConfounders = 3, 
                                    confounder_dist = c("bin", "cat_norm", "norm"),
                                    cat_confounders = c(5), prop_independent = c(0, 1, 0), prob_binary = c(0.5),
                                    var_geno = 0.6, var_noise = 0.1, var_confounders = 0.9, 
                                    totalSNPeffect = 0.01, prop_shared = 0.6, mConfounders = c(30),
                                    sdConfounders = c(20), N_causalSNPs = 100)
p <- runAnalysis_AFComps(test, title = "1,000 individuals | 3 covariates | 10,000 SNPs (100 causal)")
ggsave(plot = p, filename = "1000obs_3covariate.png", device = "png",
       height = 9, width = 12, units = "in")
saveRDS(test, "simDat_1000obs_3covariate.rds")

# now do 1 covariate (will do binary for sex)
print("Doing 100,000 Inds, 1 confounder")
test <- createSimulation_noPopStruc(seed = 2023, NInds = 100000, prop_case = 1, NSNPs = 10000,
                                    NPhenos = 1, NConfounders = 1, 
                                    confounder_dist = c("bin"),
                                    prop_independent = c(.5), prob_binary = c(0.5),
                                    var_geno = 0.6, var_noise = 0.1, var_confounders = 0.9, 
                                    totalSNPeffect = 0.01, prop_shared = 0.6, N_causalSNPs = 100)
p <- runAnalysis_AFComps(test, title = "100,000 individuals | 1 covariate | 10,000 SNPs (100 causal)")
ggsave(plot = p, filename = "100000obs_1covariate.png", device = "png",
       height = 9, width = 12, units = "in")
saveRDS(test, "simDat_100000obs_1covariate.rds")

print("Doing 10,000 Inds, 1 confounder")
test <- createSimulation_noPopStruc(seed = 2024, NInds = 10000, prop_case = 1, NSNPs = 10000,
                                    NPhenos = 1, NConfounders = 1, 
                                    confounder_dist = c("bin"),
                                    prop_independent = c(1), prob_binary = c(0.5),
                                    var_geno = 0.6, var_noise = 0.1, var_confounders = 0.9, 
                                    totalSNPeffect = 0.01, prop_shared = 0.6, N_causalSNPs = 100)
p <- runAnalysis_AFComps(test, title = "10,000 individuals | 1 covariate | 10,000 SNPs (100 causal)")
ggsave(plot = p, filename = "10000obs_1covariate.png", device = "png",
       height = 9, width = 12, units = "in")
saveRDS(test, "simDat_10000obs_1covariate.rds")

print("Doing 1,000 Inds, 1 confounder")
test <- createSimulation_noPopStruc(seed = 2025, NInds = 1000, prop_case = 1, NSNPs = 10000,
                                    NPhenos = 1, NConfounders = 1, 
                                    confounder_dist = c("bin"),
                                    prop_independent = c(0), prob_binary = c(0.5),
                                    var_geno = 0.6, var_noise = 0.1, var_confounders = 0.9, 
                                    totalSNPeffect = 0.01, prop_shared = 0.6, N_causalSNPs = 100)
p <- runAnalysis_AFComps(test, title = "1,000 individuals | 1 covariate | 10,000 SNPs (100 causal)")
ggsave(plot = p, filename = "1000obs_1covariate.png", device = "png",
       height = 9, width = 12, units = "in")
saveRDS(test, "simDat_1000obs_1covariate.rds")

# now with no covariates
print("Doing 100,000 Inds, 0 confounder")
test <- createSimulation_noPopStruc(seed = 2026, NInds = 100000, prop_case = 1, NSNPs = 10000,
                                    NPhenos = 1, NConfounders = 1, 
                                    confounder_dist = c("bin"),
                                    prop_independent = c(0), prob_binary = c(0.5),
                                    var_geno = 0.6, var_noise = 0.99999, var_confounders = 0.00001, 
                                    totalSNPeffect = 0.01, prop_shared = 0.6, N_causalSNPs = 100)
p <- runAnalysis_AFComps(test, title = "100,000 individuals | 0 covariate | 10,000 SNPs (100 causal)")
ggsave(plot = p, filename = "100000obs_0covariate.png", device = "png",
       height = 9, width = 12, units = "in")
saveRDS(test, "simDat_100000obs_0covariate.rds")

print("Doing 10,000 Inds, 0 confounder")
test <- createSimulation_noPopStruc(seed = 2027, NInds = 10000, prop_case = 1, NSNPs = 10000,
                                    NPhenos = 1, NConfounders = 1, 
                                    confounder_dist = c("bin"),
                                    prop_independent = c(0), prob_binary = c(0.5),
                                    var_geno = 0.6, var_noise = 0.99999, var_confounders = 0.00001, 
                                    totalSNPeffect = 0.01, prop_shared = 0.6, N_causalSNPs = 100)
p <- runAnalysis_AFComps(test, title = "10,000 individuals | 0 covariate | 10,000 SNPs (100 causal)")
ggsave(plot = p, filename = "10000obs_0covariate.png", device = "png",
       height = 9, width = 12, units = "in")
saveRDS(test, "simDat_10000obs_0covariate.rds")

print("Doing 1,000 Inds, 0 confounder")
test <- createSimulation_noPopStruc(seed = 2028, NInds = 1000, prop_case = 1, NSNPs = 10000,
                                    NPhenos = 1, NConfounders = 1, 
                                    confounder_dist = c("bin"),
                                    prop_independent = c(0), prob_binary = c(0.5),
                                    var_geno = 0.6, var_noise = 0.99999, var_confounders = 0.00001, 
                                    totalSNPeffect = 0.01, prop_shared = 0.6, N_causalSNPs = 100)
p <- runAnalysis_AFComps(test, title = "1,000 individuals | 0 covariate | 10,000 SNPs (100 causal)")
ggsave(plot = p, filename = "1000obs_0covariate.png", device = "png",
       height = 9, width = 12, units = "in")
saveRDS(test, "simDat_10000obs_0covariate.rds")

