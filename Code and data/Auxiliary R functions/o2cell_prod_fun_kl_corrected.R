# The function o2_cell_fun_Kl_corrected is defined to compute algae productivity
# and respiration rate as measured with the O2 cell set-up. This calculation takes
# into account the reaeration rate measured in the O2 cell.

# INPUTS:
# data_set: direct reading of the csv file obtain when exporting the data from the O2 cell sofware (in csv quite oviously)
# str_results_name: a string used to identify the lines in the data_frame obtained as outputs
# n_discard_start: number of data points discarded from the start of the vectors created with each replicate
# n_discard_end: number of data points discarded from the end of the vectors created with each replicate
# T_test: temperature at which the test was performed (important for Kl calculation, °C)
# O2_sat: O2 saturation concentration as input of the calibration (mg/L)

# OUPUTS:
  # light_intensity
  # productivity
  # productivity_standardError,
  # N_light
  # R2_light
  # pValue_light,
  # respiration (mg O2/L/min)
  # respiration_standardError
  # N_dark
  # R2_dark
  # pValue_dark 

o2cell_prod_fun_kl_corrected <- function(data_set, str_results_name, n_discard_start, n_discard_end,T_test,O2_sat){
  
  
  time <- unlist(data.frame(data_set[,1]))
  o2 <- unlist(data.frame(data_set[,2]))
  light <- unlist(data.frame(data_set[,3]))
  temperature <- unlist(data.frame(data_set[,5]))
  
  n <- dim(data_set)[1]
  
  # Correction of O2 for reaeration rate
  # kl_test = 1.06*10^(-3)*exp(5.24*10^(-2)*T_test)
  # kl_test = 4.86*10^(-3)*1.08^(T_test - 25)
  kl_test = 2.823*10^(-3)
  aeration_rate <- rep(NA,times = length(o2))
  for (ikl in 1:length(o2)){
    aeration_rate[ikl] <- kl_test*(O2_sat - o2[ikl]) #mg O2/L/min
  }
  
  # DATA FORMATTING
  
  index_light = which(light > 0)
  index_dark = which(light < 0)
  index_light <- c(split(index_light, cumsum(c(1, diff(index_light) != 1)))) # separating the different light periods
  index_dark <- c(split(index_dark, cumsum(c(1, diff(index_dark) != 1)))) # separating the different dark periods
  index_light <- Filter(function(x) length(x) > n_discard_start + n_discard_end + 1 , index_light)
  index_dark <- Filter(function(x) length(x) > n_discard_start + n_discard_end + 1 , index_dark)
  
  N_exp_light <- length(index_light)
  light_exp <- vector("list",N_exp_light)
  n_exp_light = 1 
  for (i in 1:N_exp_light){
    if (length(index_light[[i]]) >= n_discard_start + n_discard_end + 1){
      light_exp[[n_exp_light]] <- list(time_exp=time[index_light[[i]][(n_discard_start + 1):(length(index_light[[i]])-n_discard_end)]],
                                       o2_exp=o2[index_light[[i]][(n_discard_start + 1):(length(index_light[[i]])-n_discard_end)]],
                                       aeration_rate_exp=aeration_rate[index_light[[i]][(n_discard_start + 1):(length(index_light[[i]])-n_discard_end)]])
      integral_aeration_rate <- rep(NA,times = length(light_exp[[n_exp_light]]$time_exp))
      for (i_iar in 1:length(integral_aeration_rate)){
        integral_aeration_rate[i_iar] <- light_exp[[n_exp_light]]$o2_exp[i_iar] - light_exp[[n_exp_light]]$o2_exp[1] -
                                         sum(light_exp[[n_exp_light]]$aeration_rate_exp[1:i_iar]*c(0,diff(light_exp[[n_exp_light]]$time_exp[1:i_iar])))
      }
      light_exp[[n_exp_light]]$integral_aeration_rate_exp <- integral_aeration_rate
      n_exp_light = n_exp_light + 1
    }
    else {
    }
  }
  
  n_exp_light <- n_exp_light - 1

  
  N_exp_dark <- length(index_dark)
  dark_exp <- list()
  n_exp_dark = 1 
  for (i in 1:N_exp_dark){
    if (length(index_dark[[i]]) >= n_discard_start + n_discard_end + 1){
      dark_exp[[n_exp_dark]] <- list(time_exp=time[index_dark[[i]][(n_discard_start + 1):(length(index_dark[[i]])-n_discard_end)]],
                                     o2_exp=o2[index_dark[[i]][(n_discard_start + 1):(length(index_dark[[i]])-n_discard_end)]],
                                     aeration_rate_exp=aeration_rate[index_dark[[i]][(n_discard_start + 1):(length(index_dark[[i]])-n_discard_end)]])
      integral_aeration_rate <- rep(NA,times = length(dark_exp[[n_exp_dark]]$time_exp))
      for (i_iar in 1:length(integral_aeration_rate)){
        integral_aeration_rate[i_iar] <- dark_exp[[n_exp_dark]]$o2_exp[i_iar] - dark_exp[[n_exp_dark]]$o2_exp[1] -
                                         sum(dark_exp[[n_exp_dark]]$aeration_rate_exp[1:i_iar]*c(0,diff(dark_exp[[n_exp_dark]]$time_exp[1:i_iar])))
      }
      dark_exp[[n_exp_dark]]$integral_aeration_rate_exp <- integral_aeration_rate
      n_exp_dark = n_exp_dark + 1
    }
    else {
    }
  }
  n_exp_dark <- n_exp_dark - 1
  
  
  # Linear regression
  
  lm_light <- vector("list",max(n_exp_light,n_exp_dark))
  lm_dark <- vector("list",max(n_exp_light,n_exp_dark))
  
  
  for (i in 1:n_exp_light){
    
    lm_light[[i]] <- lm( integral_aeration_rate_exp  ~ unlist(time_exp), light_exp[[i]])
    plot(unlist(light_exp[[i]]$time_exp),unlist(light_exp[[i]]$o2_exp))
    plot(unlist(light_exp[[i]]$time_exp),unlist(light_exp[[i]]$integral_aeration_rate_exp))
    abline(coef = lm_light[[i]]$coefficients)
  }
  
  for (i in 1:n_exp_dark){
    lm_dark[[i]] <- lm(integral_aeration_rate_exp  ~ unlist(time_exp), dark_exp[[i]])
    plot(unlist(dark_exp[[i]]$time_exp),unlist(dark_exp[[i]]$o2_exp))
    plot(unlist(dark_exp[[i]]$time_exp),unlist(dark_exp[[i]]$integral_aeration_rate_exp))
    abline(coef = lm_dark[[i]]$coefficients)
    
  }
  
  lm_light_summary <- vector("list",max(n_exp_light,n_exp_dark));
  lm_dark_summary <- vector("list",max(n_exp_light,n_exp_dark));
  
  for (i in 1:n_exp_light){
    lm_light_summary[[i]] <- cbind("summary" = summary(lm_light[[i]])$coefficients,
                                   "Rsquare" = c(summary(lm_light[[i]])$adj.r.squared,NA),
                                   "N_data" = c(length(summary(lm_light[[i]])$residuals),NA))
  }
  
  for (i in 1:n_exp_dark){
    lm_dark_summary[[i]] <- cbind("summary" = summary(lm_dark[[i]])$coefficients,
                                  "Rsquare" = c(summary(lm_dark[[i]])$adj.r.squared,NA),
                                  "N_data" = c(length(summary(lm_dark[[i]])$residuals),NA))
  }
  
  
  # There is an assumption to make here that there are as many light than dark experiments....  
  
  light_intensity = vector(mode = "double", length = max(n_exp_light,n_exp_dark))
  productivity = vector(mode = "double", length = max(n_exp_light,n_exp_dark))
  productivity_standardError =  vector(mode = "double", length = max(n_exp_light,n_exp_dark))
  N_light = vector(mode = "double", length = max(n_exp_light,n_exp_dark))
  R2_light =  vector(mode = "double", length = max(n_exp_light,n_exp_dark))
  pValue_light = vector(mode = "double", length = max(n_exp_light,n_exp_dark))
  
  respiration = vector(mode = "double", length = max(n_exp_light,n_exp_dark))
  respiration_standardError = vector(mode = "double", length = max(n_exp_light,n_exp_dark))
  N_dark = vector(mode = "double", length = max(n_exp_light,n_exp_dark))
  R2_dark = vector(mode = "double", length = max(n_exp_light,n_exp_dark))
  pValue_dark = vector(mode = "double", length = max(n_exp_light,n_exp_dark)) 
  
  
  
  
  
  for (i in 1:n_exp_light){
    light_intensity[i] = mean(unlist(light[index_light[[i]][(n_discard_start + 1):(length(index_light[[i]])-n_discard_end)]]))
    
    productivity[i] = lm_light_summary[[i]][[2,1]]
    productivity_standardError[i] = lm_light_summary[[i]][[2,2]]
    N_light[i] = lm_light_summary[[i]][[1,6]]
    R2_light[i] = lm_light_summary[[i]][[1,5]]
    pValue_light[i] = lm_light_summary[[i]][[2,4]]
  }
  
  for (i in 1:n_exp_dark){      
    respiration[i] = lm_dark_summary[[i]][[2,1]]
    respiration_standardError[i] = lm_dark_summary[[i]][[2,2]]
    N_dark[i] = lm_dark_summary[[i]][[1,6]]
    R2_dark[i] = lm_dark_summary[[i]][[1,5]]
    pValue_dark[i] = lm_dark_summary[[i]][[2,4]]
  } 
  
  results <- data.frame(light_intensity,productivity,productivity_standardError,N_light,R2_light,pValue_light,
                    respiration,respiration_standardError,N_dark,R2_dark,pValue_dark)
  names = paste(str_results_name,"exp",1:10)
  
  rownames(results) <- names[1:dim(results)[1]]
  print(results)

}
  