#require(tidyverse)
#require(DescTools)
#require(CaseControlAF)

# function to correct the SE AFs using formulas above
correct_bias_SE <- function(MAF_df, SE, N_case, N_control) {
  bias <- 1/(min(N_case, N_control)*(SE-(0.0001/(min(N_case,N_control)/max(N_case,N_control)))))
  corrected <- data.frame(MAF_case_corr = MAF_df$MAF_case + bias,
                          MAF_control_corr = MAF_df$MAF_control + bias,
                          MAF_pop_corr = MAF_df$MAF_pop + bias)
  return(corrected)
}

# function to get Case/control AFs using the different methods
getAFs <- function(data) {
  data$SumStats <- data$SumStats %>% filter(is.finite(OR))
  data$SumStats <- data$SumStats %>% filter(is.finite(SE))
  res_AF <- CaseControl_AF(N_case = data$SampleSize$N_case,
                           N_control = data$SampleSize$N_control,
                           AF_population = data$SumStats$AF_pop,
                           OR = data$SumStats$OR)
  res_SE <- CaseControl_SE(N_case = data$SampleSize$N_case,
                           N_control = data$SampleSize$N_control,
                           SE = data$SumStats$SE,
                           OR = data$SumStats$OR)
  res_SE_corr <- correct_bias_SE(MAF_df = res_SE,
                                 SE = data$SumStats$SE,
                                 N_case = data$SampleSize$N_case,
                                 N_control = data$SampleSize$N_control)
  allResults <- cbind(res_AF, res_SE, res_SE_corr)
  allResults$true_AF_case <- data$SumStats$AF_case
  allResults$true_AF_control <- data$SumStats$AF_control
  allResults$true_AF_pop <- data$SumStats$AF_pop
  allResults$SE <- data$SumStats$SE
  allResults$OR <- data$SumStats$OR
  allResults <- allResults %>% rowwise() %>%
    mutate(true_MAF_case = ifelse(true_AF_case > 0.5, 1-true_AF_case, true_AF_case),
           true_MAF_control = ifelse(true_AF_control > 0.5, 1-true_AF_control, true_AF_control),
           true_MAF_pop = ifelse(true_AF_pop > 0.5, 1-true_AF_pop, true_AF_pop),
           MAF_case_AF = ifelse(AF_case > 0.5, 1-AF_case, AF_case),
           MAF_control_AF = ifelse(AF_control > 0.5, 1-AF_control, AF_control))
  allResults$MAF_pop_AF <- (allResults$MAF_case_AF*data$SampleSize$N_case + 
                              allResults$MAF_control_AF * data$SampleSize$N_control) / 
    (data$SampleSize$N_case + data$SampleSize$N_control)
  return(allResults)
}

# function to convert to long format data
getLongData <- function(data) {
  long <- data.frame(method = c(rep("CaseControl_AF", nrow(data)), 
                                rep("CaseControl_SE", nrow(data)), 
                                rep("SE+corr", nrow(data))),
                     status = rep("Case", 3*nrow(data)),
                     estimated = c(data$MAF_case_AF, data$MAF_case, data$MAF_case_corr),
                     true = c(data$true_MAF_case, data$true_MAF_case, data$true_MAF_case))
  long <- rbind(long, 
                data.frame(method = c(rep("CaseControl_AF", nrow(data)), 
                                      rep("CaseControl_SE", nrow(data)), 
                                      rep("SE+corr", nrow(data))),
                           status = rep("Control", 3*nrow(data)),
                           estimated = c(data$MAF_control_AF, data$MAF_control, data$MAF_control_corr),
                           true = c(data$true_MAF_control, data$true_MAF_control, data$true_MAF_control)))
  return(long)
}

getCCC <- function(data) {
  ccc <- data.frame(method = c("CaseControl_AF", "CaseControl_SE", "SE+corr"),
                    status = rep("Case", 3),
                    CCC = c(CCC(data$true_MAF_case, data$MAF_case_AF, na.rm = T)$rho.c$est,
                            CCC(data$true_MAF_case, data$MAF_case, na.rm = T)$rho.c$est,
                            CCC(data$true_MAF_case, data$MAF_case_corr, na.rm = T)$rho.c$est))
  ccc <- rbind(ccc,
               data.frame(method = c("CaseControl_AF", "CaseControl_SE", "SE+corr"),
                          status = rep("Control", 3),
                          CCC = c(CCC(data$true_MAF_control, data$MAF_control_AF, na.rm = T)$rho.c$est,
                                  CCC(data$true_MAF_control, data$MAF_control, na.rm = T)$rho.c$est,
                                  CCC(data$true_MAF_control, data$MAF_control_corr, na.rm = T)$rho.c$est)))
  ccc$label <- paste0("CCC=", round(ccc$CCC, 4))
  return(ccc)
}

# final function to run whole analysis on simulated data output from simulation function
runAnalysis_AFComps <- function(dat, title = "AF Derivation Comparison") {
  res <- getAFs(dat)
  res_long <- getLongData(res)
  ccc <- getCCC(res)
  
  p <- ggplot(res_long %>% filter(method != "SE+corr"), aes(x = true, y = estimated)) +
    geom_point() +
    #scale_color_manual(values = c("grey", "black")) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
    geom_text(mapping = aes(x = .15, y = .47, label = label), 
              data = ccc %>% filter(method != "SE+corr"), size = 5, color = "black") +
    xlab("True MAF") + ylab("Estimated MAF") +
    theme_bw() +
    facet_wrap(~status + method) +
    coord_cartesian(ylim = c(0,.5)) +
    ggtitle(title) +
    theme(strip.text = element_text(size = 15),
          plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 14))
  return(p)
}
