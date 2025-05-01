
# Simulamos datos de tiempo hasta con distribucion Weibull. Aqui es donde se van a hacer las simulaciones


sim_bayes <- function(n_sim = 1, 
                      sample_size,
                      ratio,
                      rand_type,
                      block_size = NULL,
                      Tmax, 
                      scenarios_eff,
                      shape_parameter = 1,
                      censor,
                      alpha = 0.05,
                      test.type = 2,
                      IA = NULL,
                      method_IA = NULL, 
                      n_exp_events = NULL,
                      HR_1 = FALSE,
                      Plot_Power = FALSE,
                      prior = prior,
                      modelo = modelo,
                      modelo_bayes_test = modelo_bayes_test,
                      prior_type = prior_type,
                      P_HR_data_Boundary = P_HR_data_Boundary,
                      prior_gamma,
                      Plot_Power_scenarios = FALSE,
                      plot_pvalues = FALSE,
                      Plot_Control_Scenarios = FALSE,
                      desired_HRs = NULL,
                      w = NULL,
                      external_data = NULL, 
                      a0 = NULL, 
                      Dynamic_Borrowing_PP = FALSE){
  
  # "Bucle" para evaluar cada una de las combinaciones para los valores dados de
  # estimacion de efecto y diferentes tamaños muestrales.

  if(HR_1 == TRUE){
    
    hr_1 <- replicate(n = ncol(scenarios_eff), scenarios_eff[, 1])
    scenarios_eff <- rbind(scenarios_eff, hr_1)
  }
  
  scenarios <- 1:nrow(scenarios_eff)
  
  
  # Crear plantilla de todas las combinaciones de tamaño muestral, escenarios y nº de simulaciones
  
  grid <- expand.grid(scenario = scenarios,sample_size = sample_size, n_sim = 1:n_sim)
  
  # Se aplica la funcion generadora de datos para cada una de las filas de la plantilla
  # siendo el resultado una lista.
  
  start_time <- Sys.time()
  iter <<- 1
  total <- length(sample_size)*dim(scenarios_eff)[1]*n_sim
  
  result_list <- apply(grid, 1, function(x) {
    
    scale_params <- scenarios_eff[x["scenario"], ]
    scenario <- as.integer(x["scenario"])
    sample_size <- as.integer(x["sample_size"])
    n_sim <- as.integer(x["n_sim"])
    medians <- scenarios_eff[scenario, ]
    
    # Calculamos los parametros necesarios para la Weibull
    
    parameters <- medians / (log(2)^(1/shape_parameter))

    df <-  surv_analysis_bayes(sample_size = sample_size,
                              parameters = parameters,
                              ratio = ratio,
                              rand_type = rand_type,
                              Tmax = Tmax,
                              shape_parameter = shape_parameter,
                              censor = censor,
                              alpha = alpha,
                              test.type = test.type,
                              IA = IA,
                              method_IA = method_IA,
                              n_exp_events = n_exp_events,
                              prior = prior,
                              modelo = modelo,
                              modelo_bayes_test = modelo_bayes_test,
                              prior_type = prior_type,
                              P_HR_data_Boundary = P_HR_data_Boundary,
                              prior_gamma = prior_gamma,
                              w = w,
                              external_data = external_data,
                              a0 = a0,
                              Dynamic_Borrowing_PP = Dynamic_Borrowing_PP) 
    
    # colnames(df) <- c("HR", "Upper_IC", "Lower_IC", "Reject", "reject_Cred_Int", "Time_significant", "Comparison",
    #                   "MSE", "Coverage_Probability")
    
    # df[3:8] <- lapply(df[3:8], as.numeric)
    # df[10:14] <- lapply(df[10:14], as.numeric)
    df$sample_size <- x["sample_size"]
    df$n_sim <- x["n_sim"]
    df$scenario <- x["scenario"]
    
    cat("Simulacion", iter, "de", total, "\n")
    cat("Mean_Posterior:", round(df$Mean_Posterior[1], 2), 
        "Lower_Cred_Int:", round(df$Lower_Cred_Int[1], 2), 
        "Upper_Cred_Int:", round(df$Upper_Cred_Int[1], 2), "\n")
    
    # Damos forma al resultado para cada una de las filas en formato dataframe
    
    df$scenario <- scenario
    df$sample_size <- sample_size
    df$n_sim <- n_sim
    
    control_count <- length(medians) - 1
    for (i in 1:control_count){
      
      colname <- paste0("Median control", i)
      df[[colname]] <- medians[i]
      
    }
    
    colname <- paste0("Median Treatment")
    df[[colname]] <- medians[length(medians)]
    
     # Ordenamos las columnas del dataframe
    
     cols_order <- c("scenario",
                     paste0("Median control", 1:control_count),
                     "Median Treatment",
                     "sample_size",
                     # "Comparison",
                     "n_sim")

     remaining_cols <- setdiff(colnames(df), cols_order)
     df <- df[,c(cols_order, remaining_cols)]
     
     end <- Sys.time()
     time <- end - start_time
     seconds <- as.numeric(time, units = "secs")
     time <- sprintf("%02d:%02d:%02d", as.integer(seconds %/% 3600),
                     as.integer((seconds %% 3600) %/% 60),
                     as.integer(seconds %% 60))
     cat("Tiempo acumulado de simulacion:", time, "\n\n")
     iter <<- iter+1

    return(df)
     
  })
  
  # Combinar las listas de dataframe en un df solo
  analysis <- do.call(rbind, result_list)
  
  rm(grid)
  rm(result_list)
  
  # Ordenar el df en funcion del escenario de effecto, tamaño muestral y simulacion
  
  analysis <- analysis[order(analysis$scenario, analysis$Prior,
                             analysis$sample_size, analysis$n_sim),]
  
 # Crear un resumen de resultados

  summary_data <- analysis %>%
    group_by(scenario, sample_size, Prior) %>%
    summarise(
      mean_HR = mean(Mean_Posterior, na.rm = TRUE),
      mean_Upper_CrI = mean(Upper_Cred_Int, na.rm = TRUE),
      mean_Lower_CrI = mean(Lower_Cred_Int, na.rm = TRUE),
      P_HR_Data = mean(P_HR_Data, na.rm = TRUE),
      count_significant = sum(Reject),
      prop_significant = count_significant / n(),
      count_not_significant = n() - count_significant,
      prop_not_significant = count_not_significant / n(),
      count_significant_Cred_Int = sum(reject_Cred_Int),
      prop_significant_Cred_Int = count_significant_Cred_Int / n(),
      count_not_significant_Cred_Int = n() - count_significant_Cred_Int,
      prop_not_significant_Cred_Int = count_not_significant_Cred_Int / n(),
      coverage_prob = mean(Coverage_Probability, na.rm = TRUE),
      MSE = mean(MSE, na.rm = TRUE),
      WAIC_Est = mean(WAIC_Est, na.rm = TRUE),
      WAIC_SE = mean(WAIC_SE, na.rm = TRUE),
      mean_of_the_medians_HR = mean(Median_Posterior, na.rm = TRUE),
      median_HR = median(Mean_Posterior, na.rm = TRUE),
      median_Lower_CrI = median(Lower_Cred_Int, na.rm = TRUE),
      median_Upper_CrI = median(Upper_Cred_Int, na.rm = TRUE),
      .groups = "drop") %>%
    rename(!!paste("P(HR <", sprintf("%.2f", P_HR_data_Boundary), ")| data") := P_HR_Data)
  
  
  # En caso de que haya un IA, incorporar esta info al resumen
  
  if (!is.null(IA)){
    
  IA_counts <- analysis %>%
    group_by(scenario, sample_size, Prior, Test) %>%
    summarise(
      count = n(),
      prop = count / sum(count),
      .groups = "drop"
    ) %>%
    mutate(Test = paste0("IA_", Test)) %>% 
    pivot_wider(names_from = Test, values_from = count, values_fill = 0) %>% 
    rename_with(~ gsub("IA_", "", .x), starts_with("IA_")) 
  
  IA_cols <- grep("IA", colnames(IA_counts), value = TRUE)
  FA_col <- grep("FA", colnames(IA_counts), value = TRUE)
  IA_cols_sorted <- IA_cols[order(as.numeric(gsub("IA ", "", IA_cols)))]
  IA_counts <- IA_counts %>%
    select(scenario, sample_size, Prior, all_of(IA_cols_sorted), all_of(FA_col))
  
  
  # Unir la info general con la de los IAs
  
  summary_data <- left_join(summary_data, IA_counts, by = c("scenario","sample_size","Prior"))
  
  }
  
  
  # Para hacer el plot relacionando el p-valor de cada una de las simulaciones con la probabilidad de la dist.
  # a posteriori del bayesiano con el mismo dataset (se van a usar las mismas semillas), se va elaborar otra tabla
  # al margen poniendo en una columna el nº de simulacion y en otra el p-valor.
  
  if(plot_pvalues == TRUE){
    
    pvalues <- analysis %>%
      select(scenario, sample_size, Prior, n_sim, Mean_Posterior, Upper_Cred_Int, Lower_Cred_Int, Coverage_Probability, `1-post_prob`)
    
  }
  
 #Aunque quizas se pueda dejar como opcion de criterio de parada como pensaba
  
  # En caso de muchos tamaños muestrales, un solo escenario y que no haya IAs
  # Dibujar los graficos para el poder y poder&T1E
  
  if(Plot_Power == TRUE & is.null(IA) & dim(scenarios_eff)[1] == 2){
    
    data <- subset(summary_data, scenario == 1)
    
    desired_powers <- c(0.8, 0.9)
    
    # sample_sizes_for_desired_powers <- sapply(desired_powers, function(x) 
    #   min(data$sample_size[data$prop_significant_Cred_Int >= x]))
    
    sample_sizes_for_desired_powers <- sapply(desired_powers, function(x)
      min(data$sample_size[data$prop_significant >= x]))
    
    # Plot poder con valores de poder interesantes (0.8 y 0.9)
    
      p1 <- ggplot(data, aes(x = sample_size, y = prop_significant)) +
      #p1 <- ggplot(data, aes(x = sample_size, y = prop_significant_Cred_Int)) +
      geom_point() +
      geom_line() +
      geom_hline(yintercept = desired_powers, linetype = "dashed", color = "red") +
      scale_x_continuous(breaks = seq(min(data$sample_size)-5, 
                                      max(data$sample_size)+5, 25)) +
      labs(title = "Power vs. Sample Size",
           x = "Sample Size",
           y = "Power")
    
    # Para marcar con etiquetas el poder deseado (0.8 y 0.9) en el tamaño muestral correspondiente
    for (i in 1:length(desired_powers)) {
      p1 <- p1 + annotate("text", x = sample_sizes_for_desired_powers[i], y = desired_powers[i], 
                        label = sample_sizes_for_desired_powers[i], vjust = -0.5, hjust = 1, color = "blue", size = 3.5)
    }
    
    # Plot2: Poder y error de tipo I
 
    df <- data.frame(
      sample_size = summary_data$sample_size,
      HR = summary_data$scenario,
      #power_and_T1E = summary_data$prop_significant_Cred_Int)
      power_and_T1E = summary_data$prop_significant)
    # Nueva columna para indicar a partir de cuándo se hace el zoom para el ET1
    
    # df$zoom <- ifelse(df$power_and_T1E <= 0.1, "Zoomed", "Regular")
    
    x = scenarios_eff / (log(2)^(1/shape_parameter))
    
    # suppressWarnings(
      p2 <- ggplot(data = df) +
      geom_line(aes(x = sample_size, y = power_and_T1E, color = as.factor(HR))) +
      labs(title = "Power and Type I Error by Sample Size",
           x = "Sample Size",
           y = "%",
           color = "scenarios") +
      scale_color_discrete(labels = c(paste("HR = ", round((x[1,1]/x[1,2])^shape_parameter, 2)), "HR = 1")) +
      scale_linetype_discrete(labels = c("Power", "Type I Error")) +
      theme_minimal() 

      
      p4 <- ggplot(data, aes(x = sample_size, y = MSE)) +
        geom_point() +
        geom_line() +  labs(title = "MSE vs. Sample Size",
                            x = "Sample Size",
                            y = "MSE")
      
    output <- list(summary_data,p1,p2,p4)
    
    if(plot_pvalues == TRUE){
      
      output <- list(summary_data,pvalues,p1,p2,p4)
      return(output)
      
    }
    
    return(output)
    
  }
  
  if (Plot_Power_scenarios == TRUE & is.null(IA)) {
    
    data <- summary_data %>%
      group_by(scenario) %>%
      summarize(prop_significant = sum(count_significant) / (sum(count_significant) + sum(count_not_significant)),
                mean_HR = mean(mean_HR), mean_MSE = mean(MSE))
    
    data$desired_HRs <- desired_HRs[1:nrow(data)]
    
    p3 <- ggplot(data, aes(x = desired_HRs, y = prop_significant)) +
      geom_point(size = 3, color = "steelblue") +
      geom_line(linewidth = 1, color = "steelblue") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
      geom_hline(yintercept = c(0.8, 0.9), linetype = "dotted", color = c("purple", "green"), size = 1) +
      scale_x_continuous(breaks = seq(min(desired_HRs), max(desired_HRs), 0.1), labels = sprintf("%.2f", seq(min(desired_HRs), max(desired_HRs), 0.1))) +
      scale_y_continuous(labels = scales::percent) +
      labs(title = "Power vs Hazard Ratio",
           x = "Theoretical Hazard Ratio",
           y = "Power") +
      theme_minimal() +
      theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12))
    
    p3 <- p3 +
      annotate("text", x = max(desired_HRs), y = 0.8, label = "80% power", hjust = 1, vjust = -0.5, size = 3.5, color = "purple") +
      annotate("text", x = max(desired_HRs), y = 0.9, label = "90% power", hjust = 1, vjust = -0.5, size = 3.5, color = "green")
    
    for (i in 1:nrow(data)) {
      p3 <- p3 + annotate("text", x = data$desired_HRs[i], y = data$prop_significant[i],
                          label = sprintf("%.3f", data$mean_HR[i]), hjust = 1, vjust = 1.5, size = 3, color = "black")
    }
    
    summary_data <- summary_data %>%
      arrange(scenario)
    
    # Dibujar el MSE para cada uno de los escenarios
    
    p5 <- ggplot(data, aes(x = desired_HRs, y = mean_MSE)) +
      geom_point() +
      geom_line() +  labs(title = "MSE vs. Scenarios",
                          x = "Scenarios",
                          y = "MSE")
    
    output <- list(summary_data,p3,p5)
    
    if(plot_pvalues == TRUE){
      
      output <- list(summary_data,pvalues,p3,p5)
      return(output)
      
    }
    
    return(output)
  }
  
  if (Plot_Control_Scenarios == TRUE & is.null(IA)) {
    
    scenarios_eff_df <- as.data.frame(scenarios_eff)
    colnames(scenarios_eff_df) <- c("median_control", "median_treatment")
    
    scenarios_eff_df$scenario <- seq_len(nrow(scenarios_eff_df))
    
    data <- summary_data %>%
      group_by(scenario) %>%
      summarize(prop_significant = sum(count_significant) / (sum(count_significant) + sum(count_not_significant)),
                mean_HR = mean(mean_HR), mean_MSE = mean(MSE)) %>%
      left_join(scenarios_eff_df, by = "scenario")
    
    # Dibujar los diferentes efectos del brazo control en terminos de medianas con los diferentes valores de MSE.
    
    data_filtered_diff_not_zero <- data %>%
      filter(median_control != median_treatment)
    
    data_filtered_diff_zero <- data %>%
      filter(median_control == median_treatment)
    
    median_control_center <- median(data_filtered_diff_not_zero$median_control)
    median_control_center_zero_diff <- median(data_filtered_diff_zero$median_control)
    
    p6 <-  ggplot(data_filtered_diff_not_zero, aes(x = median_control, y = mean_MSE, color = mean_MSE)) +
      geom_line() +
      geom_vline(xintercept = median_control_center, linetype = "dashed") +
      labs(x = "Median Control", y = "Mean MSE") +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red") +
      scale_x_continuous(breaks = seq(min(data_filtered_diff_not_zero$median_control), max(data_filtered_diff_not_zero$median_control), by = 0.5)) +
      scale_y_continuous(limits = c(0, 0.075))
    
    # Dibujar el poder en funcion de diferentes efectos en mediana.
    
    p7 <- ggplot(data_filtered_diff_not_zero, aes(x = median_control, y = prop_significant, color = prop_significant)) +
      geom_line() +
      geom_vline(xintercept = median_control_center, linetype = "dashed") +
      labs(x = "Median Control", y = "Prop. Significant") +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red") +
      scale_x_continuous(breaks = seq(min(data_filtered_diff_not_zero$median_control), max(data_filtered_diff_not_zero$median_control), by = 0.5)) +
      scale_y_continuous(limits = c(0, 1))
    
    # Dibujar los diferentes efectos del brazo control en terminos de medianas con el ET1.
    
    p8 <-  ggplot(data_filtered_diff_zero, aes(x = median_control, y = prop_significant, color = prop_significant)) +
      geom_line() +
      geom_vline(xintercept = median_control_center_zero_diff, linetype = "dashed") +
      geom_hline(yintercept = 0.025, linetype = "dashed", color = "green") +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "purple") +
      labs(x = "Median Control", y = "Prop. Significant") +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red") +
      scale_x_continuous(breaks = seq(min(data_filtered_diff_zero$median_control), max(data_filtered_diff_zero$median_control), by = 0.5)) +
      scale_y_continuous(limits = c(0, 0.25))
    
    
    
    output <- list(summary_data,p6,p7,p8)
    
    if(plot_pvalues == TRUE){
      
      output <- list(summary_data,pvalues,p6,p7,p8)
      return(output)
      
    }else{
      output <- list(summary_data,p6,p7,p8)
      return(output)}
    
    
  }
  
  
  
  if(plot_pvalues == TRUE){
    
    output <- list(summary_data, pvalues)
    return(output)
    
  }
  
  if(plot_pvalues == FALSE & Plot_Power == FALSE & Plot_Power_scenarios == FALSE & Plot_Control_Scenarios == FALSE){
    
    output <- summary_data
    return(output)
    
  }
  
}