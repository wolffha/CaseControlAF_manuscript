# need to create simulated AFs from autosomes and sex chromosomes

# simulate 100 variants from each chromosome

# will simulate from a sample of 5000 XX and 5000 XY case individuals
# then add 5000 of each in control

set.seed(20200)

N_XX <- 5000
N_XY <- 5000
N <- N_XX + N_XY

for(chr in 1:24) {
  refs <- runif(n = 100, min = 0, max = 1)
  if(chr < 23) {
    genos <- t(sapply(refs, function(x) {
      rmultinom(1, N, c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    afs_toadd <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*N))
  } else {
    if(chr == 23) { # simulate X chromosome variants
      XX_alleles <- (sapply(refs, function(x) {
        sum(rbinom(N_XX*2, 1, x))
      }))
      XY_alleles <- (sapply(refs, function(x) {
        sum(rbinom(N_XY, 1, x))
      }))
      
      total_alleles <- XX_alleles + XY_alleles
      afs_toadd <- total_alleles/(2*N_XX + N_XY)
    } else { # simulate Y chromosome variants
      XY_alleles <- (sapply(refs, function(x) {
        sum(rbinom(N_XY, 1, x))
      }))
      
      afs_toadd <- XY_alleles/(N_XY)
    }
  }
  if(chr == 1) {
    simDat <- data.frame(chr = chr, pos = c(1:100), case_sim_af = afs_toadd)
  } else {
    simDat <- rbind(simDat, data.frame(chr = chr, pos = c(1:100), case_sim_af = afs_toadd))
  }
}

simDat$control_sim_af <- 0

# simulate control AFs

for(chr in 1:24) {
  refs <- simDat[simDat$chr == chr,]$case_sim_af
  refs <- refs - runif(n = length(refs), min = 0, max = .025)
  refs[refs < 0] <- 0
  refs[refs > 1] <- 1
  
  if(chr < 23) {
    genos <- t(sapply(refs, function(x) {
      rmultinom(1, N, c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    afs_toadd <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*N))
  } else {
    if(chr == 23) { # simulate X chromosome variants
      XX_alleles <- (sapply(refs, function(x) {
        sum(rbinom(N_XX*2, 1, x))
      }))
      XY_alleles <- (sapply(refs, function(x) {
        sum(rbinom(N_XY, 1, x))
      }))
      
      total_alleles <- XX_alleles + XY_alleles
      afs_toadd <- total_alleles/(2*N_XX + N_XY)
    } else { # simulate Y chromosome variants
      XY_alleles <- (sapply(refs, function(x) {
        sum(rbinom(N_XY, 1, x))
      }))
      
      afs_toadd <- XY_alleles/(N_XY)
    }
  }
  simDat[simDat$chr == chr, ]$control_sim_af <- afs_toadd
}

simDat$total_sim_af <- (simDat$case_sim_af*10000 + simDat$control_sim_af*10000)/(10000 + 10000)

# need to add OR and SE
OR <- rep(0, nrow(simDat))
SE <- rep(0, nrow(simDat))

for(i in 1:nrow(simDat)) {
  AF1 = simDat[i, ]$case_sim_af
  AF2 = simDat[i,]$control_sim_af
  if(simDat[i,]$chr < 23) {
    # first calculate the 2x2 tables of allele counts
    a = AF1 * 2 * 10000 
    b = (1-AF1) * 2 * 10000
    c = AF2 * 2 * 10000
    d = (1-AF2) * 2 * 10000
  } else {
    if(simDat[i,]$chr == 23) {
      a = AF1 * (2*5000 + 5000)
      b = (1-AF1) * (2*5000 + 5000)
      c = AF2 * (2*5000 + 5000)
      d = (1-AF2) * (2*5000 + 5000)
    } else {
      a = AF1 * 5000
      b = (1-AF1) * 5000
      c = AF2 * 5000
      d = (1-AF2) * 5000
    }
  }
  
  OR[i] <- (a*d)/(b*c)
  SE[i] <- sqrt(1/a + 1/b + 1/c + 1/d)
}

simDat$OR <- OR
simDat$SE <- SE

simDat[simDat$chr == 23, ]$chr <- "X"
simDat[simDat$chr == 24, ]$chr <- "Y"

head(simDat)
summary(simDat)

ggplot(simDat, aes(x = chr, y = SE)) + geom_boxplot()
ggplot(simDat, aes(x = chr, y = OR)) + geom_boxplot() + coord_cartesian(ylim = c(0,3))
ggplot(simDat, aes(x = chr, y = case_sim_af)) + geom_boxplot()
ggplot(simDat, aes(x = chr, y = control_sim_af)) + geom_boxplot()

saveRDS(simDat, "C:/Users/HayBa/Documents/Graduate School/CaseControlAF/WholeGenomeSimDat.rds")

