createSimulation <- function(seed = 2024, NInds, prop_case = 0.5, NSNPs, NPhenos = 1,
                             NConfounders, confounder_dist, cat_confounders = NULL,
                             prob_binary = NULL, prop_independent = NULL,
                             var_geno, var_noise, var_confounders, totalSNPeffect = 0.01,
                             prop_shared = 0.6, mConfounders = 0, sdConfounders = 1) {
  
  #require(PhenotypeSimulator)
  #require(tidyverse)
  
  set.seed(seed)
  # first simulate and standardize genotypes
  genos <- simulateGenotypes(N = NInds, # number of individuals
                             NrSNP = NSNPs, # number of SNPs
                             frequencies = seq(0.01, 1, by = 0.01), # allele frequencies to sample from
                             verbose = F)
  genos_std <- standardiseGenotypes(genos$genotypes)
  
  # use the estimated kindship to get genetic effects
  kinship <- getKinship(N = NInds, X=genos_std)
  genBg <- geneticBgEffects(N=NInds,
                            kinship = kinship,
                            P = NPhenos)
  
  # simulate confounders
  noiseFixed <- noiseFixedEffects(N = NInds, # number of individuals
                                  P = NPhenos, # number of phenotypes
                                  NrFixedEffects = NConfounders, # number of confounders
                                  NrConfounders = rep(1, NConfounders), # # confounders to simulate for each dist
                                  distConfounders = confounder_dist, # distributions 
                                  catConfounders = cat_confounders, # number categories for each "cat_norm or cat_unif"
                                  probConfounders = prob_binary,# probability for each "bin"
                                  pIndependentConfounders = prop_independent,
                                  pTraitIndependentConfounders = .5,
                                  mConfounders = mConfounders,
                                  sdConfounders = sdConfounders) 
  
  # simulate random noise to add
  noiseBg <- noiseBgEffects(N = NInds, P = NPhenos)
  
  # Rescale everything
  # set parameters
  genVar <- var_geno # variance explained by genotypes
  noiseVar <- 1 - genVar # variance explained bynoise
  
  phi <- var_noise # proportion of variance for random noise
  delta <- var_confounders # proportion of variance for covariates
  shared <- prop_shared
  independent <- 1 - shared
  
  noiseBg_shared_scaled <- rescaleVariance(component = noiseBg$shared,
                                           propvar = shared*phi*noiseVar)
  noiseBg_independent_scaled <- rescaleVariance(component = noiseBg$independent,
                                                propvar = independent*phi*noiseVar)
  
  noiseFixed_shared_scaled <- rescaleVariance(component = noiseFixed$shared,
                                              propvar = shared*delta*noiseVar)
  noiseFixed_independent_scaled <- rescaleVariance(component = noiseFixed$independent,
                                                   propvar = independent*delta*noiseVar)
  
  genBg_shared_scaled <- rescaleVariance(component = genBg$shared,
                                         propvar = shared*1*genVar)
  genBg_independent_scaled <- rescaleVariance(component = genBg$independent,
                                              propvar = shared*1*genVar)
  
  # ensure total adds up to 1
  total <- shared*phi*noiseVar + independent*phi*noiseVar +
    shared*delta*noiseVar + independent*delta*noiseVar +
    shared*1*genVar + independent*1*genVar
  if(total != 1) {
    return("ERROR: Total variance != 1")
  }
  
  pheno <- scale(genBg_shared_scaled$component + genBg_independent_scaled$component +
                   noiseBg_shared_scaled$component + noiseBg_independent_scaled$component +
                   noiseFixed_shared_scaled$component + noiseFixed_independent_scaled$component)
  
  # calculate the number of cases and controls
  nCase <- length(pheno[pheno > 0])
  nControl <- length(pheno[pheno <= 0])
  
  print(paste0("Original nCase: ", nCase, " nControl: ", nControl))
  
  # keep track of which indices were kept
  original_case <- which(pheno > 0)
  original_control <- which(pheno <= 0)
  
  if(prop_case < 1) {
    desired_nCase = round(prop_case*nControl)
    desired_nControl = nControl
  } else if(prop_case > 1) {
    desired_nControl = round(nControl/prop_case)
    desired_nCase = nCase
  } else {
    desired_nCase = nCase
    desired_nControl = nControl
  }
  
  kept_case <- original_case[1:desired_nCase]
  kept_control <- original_control[1:desired_nControl]
  indices_keep <- c(kept_case, kept_control)
  indices_keep <- sort(indices_keep)
  
  print(paste0("Final nCase: ", desired_nCase, " nControl: ", desired_nControl))
  
  sampleSizes <- data.frame(N_case = desired_nCase,
                            N_control = desired_nCase)
  
  # get the case and control AFs
  # first subset the genotype matrix to just cases
  geno_sub <- genos$genotypes[which(pheno > 0),][1:desired_nCase,]
  
  case_afs <- apply(geno_sub, 2, function(x) getAlleleFrequencies(x))
  case_afs <- case_afs[1,]
  
  # just controls
  geno_sub <- genos$genotypes[which(pheno <= 0),][1:desired_nControl,]
  control_afs <- apply(geno_sub, 2, function(x) getAlleleFrequencies(x))
  control_afs <- control_afs[1,]
  
  # function to get the OR and SE using regression
  getORandSE <- function(genotype_matrix, phenotype_vector, covariate) {
    # Assuming your genotype matrix is named "genotype_matrix" where rows are individuals and columns are SNPs
    # Assuming your phenotype vector is named "phenotype_vector"
    
    # Fit logistic regression models for each SNP
    results <- lapply(1:ncol(genotype_matrix), function(i) {
      # Extract SNP genotype data
      snp_genotype <- genotype_matrix[, i]
      simdat <- data.frame(phenotype = phenotype_vector, 
                           genotype = snp_genotype, 
                           covariate$shared[indices_keep],
                           covariate$independent[indices_keep])
      #print(head(simdat))
      
      # Fit logistic regression model
      model <- glm(phenotype ~., data = simdat, family = binomial())
      
      # Extract coefficients and standard errors
      odds_ratio <- exp(coef(model)[2])
      se <- summary(model)$coefficients[2,2] # gets the std error for the beta estimate for genotype
      
      # Return results
      return(list(se = se, odds_ratio = odds_ratio))
    })
    
    # Extract coefficients, standard errors, and odds ratios
    #coefficients <- sapply(results, function(x) x$coef[2])
    standard_errors <- sapply(results, function(x) x$se)
    odds_ratios <- sapply(results, function(x) x$odds_ratio)
    
    # Create a data frame to store results
    result_df <- data.frame(SE = standard_errors,
                            OR = odds_ratios)
    return(result_df)
    
  }
  
  geno_matrix <- genos$genotypes[indices_keep, ]
  pheno_vector <- ifelse(pheno > 0, 1, 0)[indices_keep]
  sumStats <- getORandSE(genotype_matrix = geno_matrix, phenotype_vector = pheno_vector,
                         covariate = noiseFixed)
  
  # now calculate OR and SE(log(OR))
  # first calculate the 2x2 tables of allele counts
  # a = case_afs * 2 * nCase
  # b = (1-case_afs) * 2 * nCase
  # c = control_afs * 2 * nControl
  # d = (1-control_afs) * 2 * nControl
  # 
  # OR <- (a*d)/(b*c)
  # SE <- sqrt(1/a + 1/b + 1/c + 1/d)
  # SE <- round(SE, 5)
  # 
  # sumStats <- data.frame(OR = OR, SE = SE)
  
  sumStats$AF_pop <- ((case_afs*nCase) + (control_afs*nControl))/(nCase + nControl)
  sumStats$AF_case <- case_afs
  sumStats$AF_control <- control_afs
  
  sumStats <- sumStats %>% rowwise() %>%
    mutate(MAF_case = ifelse(AF_case > 0.5, 1-AF_case, AF_case),
           MAF_control = ifelse(AF_control > 0.5, 1-AF_control, AF_control),
           MAF_pop = ifelse(AF_pop > 0.5, 1-AF_pop, AF_pop))
  
  return(list(SumStats = sumStats,
              SampleSize = sampleSizes))
}

createSimulation_noPopStruc <- function(seed = 2024, NInds, prop_case = 0.5, NSNPs, NPhenos = 1,
                             NConfounders, confounder_dist, cat_confounders = NULL,
                             prob_binary = NULL, prop_independent = NULL,
                             var_geno, var_noise, var_confounders, totalSNPeffect = 0.01,
                             prop_shared = 0.6, mConfounders = 0, sdConfounders = 1,
                             N_causalSNPs = 0) {

  require(PhenotypeSimulator)
  require(tidyverse)
  
  set.seed(seed)
  # first simulate and standardize genotypes
  genos <- simulateGenotypes(N = NInds, # number of individuals
                           NrSNP = NSNPs, # number of SNPs
                           frequencies = seq(0.01, 1, by = 0.01), # allele frequencies to sample from
                           verbose = F)
  genos_std <- standardiseGenotypes(genos$genotypes)

  # this section is for infinitesimal genetic effects (aka population structure)
  # I think this section is very slow so I will try without it
  # use the estimated kindship to get genetic effects
  #kinship <- getKinship(N = NInds, X=genos_std)
  # genBg <- geneticBgEffects(N=NInds,
  #                         kinship = kinship,
  #                         P = NPhenos)
  
  # get genetic variant effects for causal SNPs
  if(N_causalSNPs > 0) {
    causalSNPs <- getCausalSNPs(N = NInds,
                                genotypes = genos$genotypes,
                                NrCausalSNPs = N_causalSNPs,
                                verbose = F)
    genFixed = geneticFixedEffects(N = NInds,
                                   P = NPhenos,
                                   X_causal = causalSNPs)
  }
  
  # simulate confounders
  noiseFixed <- noiseFixedEffects(N = NInds, # number of individuals
                                P = NPhenos, # number of phenotypes
                                NrFixedEffects = NConfounders, # number of confounders
                                NrConfounders = rep(1, NConfounders), # # confounders to simulate for each dist
                                distConfounders = confounder_dist, # distributions 
                                catConfounders = cat_confounders, # number categories for each "cat_norm or cat_unif"
                                probConfounders = prob_binary,# probability for each "bin"
                                pIndependentConfounders = prop_independent,
                                pTraitIndependentConfounders = .5,
                                mConfounders = mConfounders,
                                sdConfounders = sdConfounders) 
  
  # simulate random noise to add
  noiseBg <- noiseBgEffects(N = NInds, P = NPhenos)
  
  # Rescale everything
  # set parameters
  genVar <- var_geno # variance explained by genotypes
  noiseVar <- 1 - genVar # variance explained bynoise

  phi <- var_noise # proportion of variance for random noise
  delta <- var_confounders # proportion of variance for covariates
  shared <- prop_shared
  independent <- 1 - shared
  
  noiseBg_shared_scaled <- rescaleVariance(component = noiseBg$shared,
                                         propvar = shared*phi*noiseVar)
  noiseBg_independent_scaled <- rescaleVariance(component = noiseBg$independent,
                                                propvar = independent*phi*noiseVar)
  
  noiseFixed_shared_scaled <- rescaleVariance(component = noiseFixed$shared,
                                              propvar = shared*delta*noiseVar)
  noiseFixed_independent_scaled <- rescaleVariance(component = noiseFixed$independent,
                                                   propvar = independent*delta*noiseVar)
  
  # genBg_shared_scaled <- rescaleVariance(component = genBg$shared,
  #                                        propvar = shared*1*genVar)
  # genBg_independent_scaled <- rescaleVariance(component = genBg$independent,
  #                                             propvar = shared*1*genVar)
  
  if(N_causalSNPs > 0) {
    genFixed_shared_scaled <- rescaleVariance(component = genFixed$shared,
                                           propvar = shared*1*genVar)
    genFixed_independent_scaled <- rescaleVariance(component = genFixed$independent,
                                                propvar = shared*1*genVar)
  }
  
  # ensure total adds up to 1
  total <- shared*phi*noiseVar + independent*phi*noiseVar +
    shared*delta*noiseVar + independent*delta*noiseVar +
    shared*1*genVar + independent*1*genVar
  if(total != 1) {
    return("ERROR: Total variance != 1")
  }
  
  if(N_causalSNPs > 0) {
    x <- list(genFixed_shared_scaled$component,
              genFixed_independent_scaled$component,
              noiseBg_shared_scaled$component,
              noiseBg_independent_scaled$component,
              noiseFixed_shared_scaled$component,
              noiseFixed_independent_scaled$component)
    pheno = -100
    for(i in 1:length(x)) {
      if(!is.null(x[i][[1]])[1] & pheno[1] == -100) {
        pheno = x[i][[1]]
      } else if(!is.null(x[i][[1]])[1] & pheno[1] != -100) {
        pheno = pheno + x[i][[1]]
      }
    }
    pheno = scale(pheno)
    # pheno <- scale(genFixed_shared_scaled$component + genFixed_independent_scaled$component +
    #                noiseBg_shared_scaled$component + noiseBg_independent_scaled$component +
    #                noiseFixed_shared_scaled$component + noiseFixed_independent_scaled$component)
  } else {
    pheno <- scale(noiseBg_shared_scaled$component + noiseBg_independent_scaled$component +
                   noiseFixed_shared_scaled$component + noiseFixed_independent_scaled$component)
  }
  
  # calculate the number of cases and controls
  nCase <- length(pheno[pheno > 0])
  nControl <- length(pheno[pheno <= 0])
  
  print(paste0("Original nCase: ", nCase, " nControl: ", nControl))
  
  # keep track of which indices were kept
  original_case <- which(pheno > 0)
  original_control <- which(pheno <= 0)
  
  if(prop_case < 1) {
    desired_nCase = round(prop_case*nControl)
    desired_nControl = nControl
  } else if(prop_case > 1) {
    desired_nControl = round(nControl/prop_case)
    desired_nCase = nCase
  } else {
    desired_nCase = nCase
    desired_nControl = nControl
  }
  
  kept_case <- original_case[1:desired_nCase]
  kept_control <- original_control[1:desired_nControl]
  indices_keep <- c(kept_case, kept_control)
  indices_keep <- sort(indices_keep)
  
  print(paste0("Final nCase: ", desired_nCase, " nControl: ", desired_nControl))
  
  sampleSizes <- data.frame(N_case = desired_nCase,
                            N_control = desired_nControl)
  
  # get the case and control AFs
  # first subset the genotype matrix to just cases
  geno_sub <- genos$genotypes[which(pheno > 0),][1:desired_nCase,]

  case_afs <- apply(geno_sub, 2, function(x) getAlleleFrequencies(x))
  case_afs <- case_afs[1,]
  
  # just controls
  geno_sub <- genos$genotypes[which(pheno <= 0),][1:desired_nControl,]
  control_afs <- apply(geno_sub, 2, function(x) getAlleleFrequencies(x))
  control_afs <- control_afs[1,]
  
  # function to get the OR and SE using regression
  getORandSE <- function(genotype_matrix, phenotype_vector, covariate) {
    # Assuming your genotype matrix is named "genotype_matrix" where rows are individuals and columns are SNPs
    # Assuming your phenotype vector is named "phenotype_vector"

    # Fit logistic regression models for each SNP
    results <- lapply(1:ncol(genotype_matrix), function(i) {
      # Extract SNP genotype data
      snp_genotype <- genotype_matrix[, i]
      if(is.null(covariate$shared) & is.null(covariate$independent)) {
        simdat <- data.frame(phenotype = phenotype_vector, 
                             genotype = snp_genotype)
      } else if(is.null(covariate$shared)) {
        simdat <- data.frame(phenotype = phenotype_vector, 
                             genotype = snp_genotype,
                             covariate$independent[indices_keep])
      } else if(is.null(covariate$independent)) {
        simdat <- data.frame(phenotype = phenotype_vector, 
                             genotype = snp_genotype, 
                             covariate$shared[indices_keep])
      } else {
        simdat <- data.frame(phenotype = phenotype_vector, 
                             genotype = snp_genotype, 
                             covariate$shared[indices_keep],
                             covariate$independent[indices_keep])
      }
      
      #print(head(simdat))

      # Fit logistic regression model
      model <- glm(phenotype ~., data = simdat, family = binomial())

      # Extract coefficients and standard errors
      odds_ratio <- exp(coef(model)[2])
      se <- summary(model)$coefficients[2,2] # gets the std error for the beta estimate for genotype

      # Return results
      return(list(se = se, odds_ratio = odds_ratio))
    })

    # Extract coefficients, standard errors, and odds ratios
    #coefficients <- sapply(results, function(x) x$coef[2])
    standard_errors <- sapply(results, function(x) x$se)
    odds_ratios <- sapply(results, function(x) x$odds_ratio)

    # Create a data frame to store results
    result_df <- data.frame(SE = standard_errors,
                            OR = odds_ratios)
    return(result_df)

  }
  
  geno_matrix <- genos$genotypes[indices_keep, ]
  pheno_vector <- ifelse(pheno > 0, 1, 0)[indices_keep]
  sumStats <- getORandSE(genotype_matrix = geno_matrix, phenotype_vector = pheno_vector,
                         covariate = noiseFixed)

  # now calculate OR and SE(log(OR))
  # first calculate the 2x2 tables of allele counts
  # a = case_afs * 2 * nCase
  # b = (1-case_afs) * 2 * nCase
  # c = control_afs * 2 * nControl
  # d = (1-control_afs) * 2 * nControl
  # 
  # OR <- (a*d)/(b*c)
  # SE <- sqrt(1/a + 1/b + 1/c + 1/d)
  # SE <- round(SE, 5)
  # 
  # sumStats <- data.frame(OR = OR, SE = SE)
  
  sumStats$AF_pop <- ((case_afs*nCase) + (control_afs*nControl))/(nCase + nControl)
  sumStats$AF_case <- case_afs
  sumStats$AF_control <- control_afs
  
  sumStats <- sumStats %>% rowwise() %>%
    mutate(MAF_case = ifelse(AF_case > 0.5, 1-AF_case, AF_case),
           MAF_control = ifelse(AF_control > 0.5, 1-AF_control, AF_control),
           MAF_pop = ifelse(AF_pop > 0.5, 1-AF_pop, AF_pop))

  return(list(SumStats = sumStats,
              SampleSize = sampleSizes))
}
