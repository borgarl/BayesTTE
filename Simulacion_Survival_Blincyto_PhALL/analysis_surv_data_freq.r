
analysis_surv_data_freq <- function(simul_data, alpha, 
                                    IA, method_IA, 
                                    test.type, 
                                    sample_size, 
                                    n_exp_events,
                                    parameters,
                                    shape_parameter){
  
  # En el IA tendria que dividir el alfa entre el nº de comparaciones en caso de que haya
  # + de 2 brazos
  
  if(is.null(IA)){ 
    
    boundaries <- 1 
    
  }else{
    
    boundaries <- get_boundaries_IA(IA, alpha, method_IA = method_IA, 
                                    test.type = test.type,
                                    sample_size, n_exp_events)
  }
  
  ##################################################
  
  n_analyses <- ifelse(length(boundaries) > 1, length(boundaries), 1)
  
  # En caso de que haya 3 o más brazos se considerara el IA como la suma de eventos
  # de los tres brazos. Si hace falta se cambiara
  
  sum_events <- sum(simul_data$status)
  
  n_arms <- length(unique(simul_data$arm))
  treatment_arm <- n_arms
  
  res_df <- as.data.frame(matrix(c(0), nrow = n_arms - 1, ncol = 7))
  #res_df <- as.data.frame(matrix(c(0), nrow = n_arms - 1, ncol = 8))
  
  True_HR <- c()
  
  # Este bucle hay que optimizarlo (es lo que mas tiempo lleva de la simulacion)
  
  for(i in 1:(n_arms - 1)){
    
    control_arm <- i
    control_vs_treatment <- subset(simul_data, arm %in% c(control_arm, treatment_arm))
    
    reject <- "Not significant"
    j <- 1
    
    if(n_analyses == 1){
      
      fit <- coxph(Surv(time, status) ~ arm, data = control_vs_treatment)
      
      # Este codigo es para calcular el p-valor de la weibull solo para comparar con los resultados
      # del bayesiano. Cuando se deje de usar se comentara para quitarlo.
      
      # fit_weibull <- survreg(Surv(time, status) ~ arm, data = control_vs_treatment, dist="weibull")
      # coef_arm <- fit_weibull$coefficients["arm"]
      # std_error_arm <- sqrt(diag(vcov(fit_weibull)))[["arm"]]
      # z_value_arm <- coef_arm / std_error_arm
      # p_value_weibull <- as.numeric(2 * (1 - pnorm(abs(z_value_arm))))
      
      z <- summary(fit)$coefficients[4]
      p_value <- round(test.type * (1 - pnorm(abs(z))),5)
      reject <- ifelse((p_value < alpha) & (z < 0), 1, 0)
      Estimate <- c(summary(fit)$coefficients[2], exp(confint(fit))[1], exp(confint(fit))[2])
      res <- c(Estimate, reject, p_value)
      #res <- c(Estimate, reject, p_value,p_value_weibull)
      
      if(p_value < alpha){
        
        Test <- as.character(paste("FA"))
        res <- c(res, Test)
        
      }else{
        
        Test <- as.character(paste("Not significant"))
        res <- c(Estimate, reject, p_value, Test)
        #res <- c(Estimate, reject, p_value, p_value_weibull,Test)
        
      }
      
    }else if (n_analyses > 1){
      
      for(j in 1:n_analyses){
        
        if(j == length(boundaries)){
          
          fit <- coxph(Surv(time, status) ~ arm, data = control_vs_treatment)
          z <- summary(fit)$coefficients[4]
          p_value <- round(test.type * (1 - pnorm(abs(z))),5)
          Estimate <- c(summary(fit)$coefficients[2], exp(confint(fit))[1], exp(confint(fit))[2])
          
          if(p_value < boundaries[j] && (z < 0)){
            
            reject <- ifelse((p_value < boundaries[j]), 1, 0)
            Estimate <- c(summary(fit)$coefficients[2], exp(confint(fit))[1], exp(confint(fit))[2])
            Test <- as.character(paste("FA"))
            res <- c(Estimate, reject, p_value, Test)
            break
            
          }else{
            reject <- 0
          }
          
          Test <- as.character(paste("Not significant"))
          res <- c(Estimate, reject, p_value, Test)
          break
        } 
        
        events_IA <- round(IA[j] * sum_events)
        IA_data <- control_vs_treatment[cumsum(control_vs_treatment$status) <= events_IA,]
        fit <- coxph(Surv(time, status) ~ arm, data = IA_data)
        z <- summary(fit)$coefficients[4]
        p_value <- round(test.type * (1 - pnorm(abs(z))),5)
        
        if(p_value < boundaries[j] && (z < 0)){
          
          reject <- 1
          Estimate <- c(summary(fit)$coefficients[2], exp(confint(fit))[1], exp(confint(fit))[2])
          Test <- as.character(paste("IA", j))
          res <- c(Estimate, reject, p_value, Test)
          break
          
        }else{
          
          reject <- 0
          
        }
      }
    }
    
    Comparison <- as.character(paste0("Control",i, " vs Treatment"))
    res <- c(res, Comparison)
    res_df[i,] = res
    True_HR[i] <- (parameters[i]/parameters[n_arms])^shape_parameter 
    res = c()
  }
  
  res_df[, 1:5] <- lapply(res_df[, 1:5], as.numeric)
  #res_df[, 1:6] <- lapply(res_df[, 1:6], as.numeric)

  # Calculamos diferencias para medidas de rendimiento
  # 1) MSE
  # 2) Coverage probability
  
  MSE <- (res_df[,1]-True_HR)^2
  
  # 1: Incluye el True_HR dentro del intervalo y 0 no lo incluye
  
  Coverage_Probability <- ifelse((res_df[,2] < True_HR) 
                                 & (True_HR < res_df[,3]), 1, 0)
  
  res_df$MSE <- MSE
  res_df$Coverage_Probability <- Coverage_Probability

  colnames(res_df) <- c("Estimate", "Lower_IC", "Upper_IC", "Reject", "p_value", "Test", "Comparison", "MSE", "Coverage_Probability")
  #colnames(res_df) <- c("Estimate", "Lower_IC", "Upper_IC", "Reject", "p_value", "p_value_weibull", "Test", "Comparison", "MSE", "Coverage_Probability")
  
  return(res_df)
  
}
