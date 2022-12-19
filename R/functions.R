
cortest_cols <- function(data, dv, cols){
  
  data_cov <- data %>% dplyr::select({{cols}})
  data_dv <- data %>% dplyr::select({{dv}}) %>% pull()
  
  final_data <- data.frame()
  
  for(i in 1:length(data_cov)){
    
    y <- data_cov[i] %>% pull()
    corT <- broom::tidy(cor.test(as.numeric(data_dv), as.numeric(y)))
    corT$estimate <- corT$estimate
    corT$pvalue <- corT$p.value
    corT <- corT %>% mutate(var = names(data_cov[i]))
    final_data <- plyr::rbind.fill(final_data, corT)
    
  }
  
  final_data
}

collinearity <- function(data, cols){
  
  collinearity <- data.frame()
  for(i in cols){
    cor <- cortest_cols(data, i, cols )
    cor <- cor %>% mutate(voi = paste(i))
    
    collinearity <- bind_rows(collinearity, cor)
    
  }
  
  collinearity
  
}

filter_collinearity <- function(data) {
  data %>% dplyr::filter(pvalue != 0) %>%
    filter(estimate >= abs(0.7))
}


multi_collin <- function(data, type){
  lm_results <- switch(type,
                       'bf' = {lm(log(mean_bf)~log(fac_taudem_17all_int)+
                                    log(us_precip_1981_2010_cpg_all)+
                                    log(us_tmax_1981_2010_int_cpg_all) +
                                    log(swe_avg_0301_cpg_all) +
                                    log(clay_100_cpg_all)+
                                    log(slope_percent_int_cpg_all)+
                                    MAP_UNIT_N +
                                    mgmt, data = data)},
                       'ba' = {lm(log(mean_bank_an)~log(fac_taudem_17all_int)+
                                    log(ave_basin_elev)+
                                    log(streamslope_percent_int_cpg_all) +
                                    log(silt_100_cpg_all) +
                                    log(twi_100_int_cpg_all) +
                                    MAP_UNIT_N +
                                    mgmt, data = data)},
                       'wd' = {lm(log(mean_wd_2009)~log(fac_taudem_17all_int)+
                                    log(us_precip_1981_2010_cpg_all)+
                                    log(ave_ann_cwd) +
                                    log(clay_100_cpg_all) +
                                    log(swe_avg_0301_cpg_all) +
                                    log(twi_100_int_cpg_all) +
                                    log(slope_percent_int_cpg_all) +
                                    MAP_UNIT_N +
                                    mgmt, data = data)}
  )
  
  car::vif(lm_results)
  
}

linear_plot <- function(data, x, y, iv = NULL, log_y = NULL, log_x = NULL){
  
  if(is.null(iv)){
    mapping <- ggplot2::aes(.data[[x]], .data[[y]])
    r2 <- list(
      ggpmisc::stat_poly_eq(aes(label =  paste(stat(eq.label),
                                               stat(rr.label),
                                               sep = "~~~~")),
                            formula = y~x,
                            rr.digits = 2,
                            # coef.digits = 2,
                            parse = TRUE),
      ggplot2::geom_smooth(formula = y~x,
                           method = 'lm',
                           linetype = "dashed",
                           se = F)
    )
    
    
  } else {
    mapping <- ggplot2::aes(.data[[x]], .data[[y]], color = .data[[iv]])
    r2 <- list(
      ggpmisc::stat_poly_eq(aes(label =  paste(stat(eq.label),
                                               stat(rr.label),
                                               sep = "~~~~")),
                            formula = y~x,
                            rr.digits = 2,
                            # coef.digits = 2,
                            parse = TRUE),
      ggplot2::geom_smooth(formula = y~x,
                           method = 'lm',
                           linetype = "dashed",
                           se = F)
    )
    
      mapping <- ggplot2::aes(.data[[x]], .data[[y]], color = .data[[iv]])
    r2 <- Add_R2()
    
  } 
  
  if(!is.null(log_y)){
    logy <- ggplot2::scale_y_log10()
  } else {
    logy <- NULL
  }
  
  if (!is.null(log_x)){
    logx <- ggplot2::scale_x_log10()
  } else {
    logx <- NULL
  }
  
  
  data %>%
    ggplot(mapping) +
    geom_point() +
    r2 +
    custom_theme() +
    logy +
    logx
}


Add_R2 <- function(){
  
  list(
    ggpmisc::stat_poly_eq(aes(label =  paste(stat(eq.label),
                                             stat(rr.label),
                                             sep = "~~~~")),
                          formula = y~x,
                          rr.digits = 2,
                          # coef.digits = 2,
                          parse = TRUE),
    ggplot2::geom_smooth(formula = y~x,
                         method = 'lm',
                         linetype = 'dashed',
                         se = F)
  )
}

feature_selection <- function(data, type){
  
  data <- switch(type,
                 'bf' = data %>% select(mean_bf, fac_taudem_17all_int,
                                        us_precip_1981_2010_cpg_all, us_tmax_1981_2010_int_cpg_all,
                                        swe_avg_0301_cpg_all, clay_100_cpg_all, slope_percent_int_cpg_all,
                                        MAP_UNIT_N) %>%
                   mutate(across(1:7, ~log(.))),
                 'ba' = data %>% select(mean_bank_an, fac_taudem_17all_int,
                                        ave_basin_elev, streamslope_percent_int_cpg_all,
                                        silt_100_cpg_all, twi_100_int_cpg_all,
                                        MAP_UNIT_N) %>%
                   mutate(across(1:6, ~log(.))),
                 'wd' = data %>% select(mean_wd_2009,'fac_taudem_17all_int','us_precip_1981_2010_cpg_all', 'ave_ann_cwd', 'swe_avg_0301_cpg_all', 'clay_100_cpg_all', 'slope_percent_int_cpg_all',
                                        MAP_UNIT_N) %>%
                   mutate(across(1:7, ~log(.))) %>% na.omit()
  )
  
  
  set.seed(1234)
  indicesBlock <- CreateSpacetimeFolds(data, spacevar = "MAP_UNIT_N", k = 6)
  
  
  # Set up repeated k-fold cross-validation
  train.control <- caret::trainControl(method = "cv",
                                       index = indicesBlock$index,
                                       returnResamp = "all",
                                       savePredictions = 'all')
  
  # Train the model
  set.seed(136564)
  
  lasso_model <- switch(type,
                        'bf' = train(data[,2:7],
                                     data$mean_bf,
                                     method = "glmnet",metric = "RMSE",
                                     tuneGrid = expand.grid(alpha = seq(0.1,.9,by = 0.1),
                                                            lambda = seq(0.001,0.1,by = 0.001)),
                                     trControl = train.control),
                        'ba' = train(data[,2:6],
                                     data$mean_bank_an,
                                     method = "glmnet",metric = "RMSE",
                                     tuneGrid = expand.grid(alpha = seq(0.1,.9,by = 0.1),
                                                            lambda = seq(0.001,0.1,by = 0.001)),
                                     trControl = train.control
                        ),
                        'wd' = train(data[,2:7],
                                     data$mean_wd_2009,
                                     method = "glmnet",metric = "RMSE",
                                     tuneGrid = expand.grid(alpha = seq(0.1,.9,by = 0.1),
                                                            lambda = seq(0.001,0.1,by = 0.001)),
                                     trControl = train.control
                        )
  )
  lasso_model
}


## getting varImportance via stackoverlow answer https://stackoverflow.com/questions/63989057/discripencies-in-variable-importance-calculation-for-glmnet-model-in-r
varImp <- function(object, lambda = NULL, ...) {
  beta <- predict(object, s = lambda, type = "coef")
  if(is.list(beta)) {
    out <- do.call("cbind", lapply(beta, function(x) x[,1]))
    out <- as.data.frame(out)
  } else out <- data.frame(Overall = beta[,1])
  out <- abs(out[rownames(out) != "(Intercept)",,drop = FALSE])
  out <- out/max(out)
  out[order(out$Overall, decreasing = TRUE),,drop=FALSE]
}

check_if_normal_hetero <- function(data, type, var_imp){
  
  var_names <- rownames(var_imp)
  
  vars_bf <- paste0('log(mean_bf)~',paste0('log(', var_names, ')', collapse = ' + '), ' + MAP_UNIT_N + mgmt')
  vars_ba <- paste0('log(mean_bank_an)~',paste0('log(', var_names, ')', collapse = ' + '), ' + MAP_UNIT_N + mgmt')
  vars_wd <- paste0('log(mean_wd_2009)~',paste0('log(', var_names, ')', collapse = ' + '), ' + MAP_UNIT_N + mgmt')
  
  lm_results <- switch(type,
                       'bf' = {lm(formula(noquote(vars_bf)), data = data)},
                       'ba' = {lm(formula(noquote(vars_ba)), data = data)},
                       'wd' = {lm(formula(noquote(vars_wd)), data = data)}
  )
  
  list(norm = check_normality(lm_results), hetero = check_heteroscedasticity(lm_results), model = lm_results)
  
}















