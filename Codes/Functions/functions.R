eval_with_timeout <- function(expr, envir = parent.frame(), timeout, on_timeout = c("error", "warning", "silent")) {
  # substitute expression so it is not executed as soon it is used
  expr <- substitute(expr)
  
  # match on_timeout
  on_timeout <- match.arg(on_timeout)
  
  # execute expr in separate fork
  myfork <- parallel::mcparallel({
    eval(expr, envir = envir)
  }, silent = FALSE)
  
  # wait max n seconds for a result.
  myresult <- parallel::mccollect(myfork, wait = FALSE, timeout = timeout)
  # kill fork after collect has returned
  tools::pskill(myfork$pid, tools::SIGKILL)
  tools::pskill(-1 * myfork$pid, tools::SIGKILL)
  
  # clean up:
  parallel::mccollect(myfork, wait = FALSE)
  
  # timeout?
  if (is.null(myresult)) {
    if (on_timeout == "error") {
      stop("reached elapsed time limit")
    } else if (on_timeout == "warning") {
      warning("reached elapsed time limit")
    } else if (on_timeout == "silent") {
      myresult <- NA
    }
  }
  
  # move this to distinguish between timeout and NULL returns
  myresult <- myresult[[1]]
  
  if ("try-error" %in% class(myresult)) {
    stop(attr(myresult, "condition"))
  }
  
  # send the buffered response
  return(myresult)
}

FitFlextableToPage <- function(ft, pgwidth = 6.5){
  
  ft_out <- ft 
  
  ft_out <- width(ft_out, width = dim(ft_out)$widths*pgwidth /(flextable_dim(ft_out)$widths))
  return(ft_out)
}

run_se <- function(data){
  
  #--------------------------
  # Spatial Error Model estimation
  #--------------------------
  cent <- getSpPPolygonsLabptSlots(as(data,"Spatial"))
  dnb <- dnearneigh(cent,0,15)
  dsts <- nbdists(dnb,cent)
  idw <- lapply(dsts, function(x) 1/x)
  wq <- nb2listw(dnb,glist=idw,style = "W",zero.policy = TRUE)
  se_res<-errorsarlm(yield~n+I(n*n), 
                     data = data, listw=wq, zero.policy=TRUE, na.action=na.omit)
  
  return(se_res)
}

run_qp <- function(data){
  
  #--------------------------
  # Quadratic Plateau Model estimation
  #--------------------------
  data_temp <- as.data.frame(data)[,c("n", "yield")]
  qp_res<-nlsfit(data_temp, model=4)
  return(qp_res$Parameters)
  
}

read_rmd <- function(file_name) {
  file_name <-
    here(file_name)
  
  rmd_file <- readLines(file_name)
  
  return(rmd_file)
}

data_summary <- function(ffy){

    data_temp <- here("Data", ffy, "analysis_data.shp") %>% 
      read_sf()%>% 
      filter(!is.na(yild_vl))%>%
      data.table()%>%
      dplyr::mutate(farm = strsplit(ffy, "_")[[1]][1], 
                    field = strsplit(ffy, "_")[[1]][2], 
                    year = strsplit(ffy, "_")[[1]][3])%>%
      dplyr::select('farm', 'field', 'year', 'yild_vl', 's_rate', 'n_rate')%>%
      mutate_if(is.numeric, round, digits = 1)
    
  return(data_temp)
}

analysis_make_report <- function(ffy, rerun = TRUE){
  
  analysis_temp_rmd <- read_rmd(
    "Code/analysis_template.Rmd"
  )
  
  analysis <- read_rmd(
    "Code/1_analysis.Rmd"
  )
  
  weather <- read_rmd(here("Code/0_1_weather.Rmd"))
  
  analysis_rmd_y <- c(analysis_temp_rmd, weather, analysis)%>%
    gsub("field-year-here", ffy, .) %>% 
    gsub("title-here", "Analysis Report", .)
  
  
  #/*=================================================*/
  #' # Write out the rmd and render
  #/*=================================================*/
  analysis_report_rmd_file_name <- here(
    "Data", 
    ffy, 
    "analysis_report_exp.Rmd"
  )
  
  analysis_report_r_file_name <- here(
    "Data", 
    ffy, 
    "for_analysis_debug.R"
  )
  
  writeLines(analysis_rmd_y, con = analysis_report_rmd_file_name)
  
  purl(analysis_report_rmd_file_name, output = analysis_report_r_file_name)
  
}

make_var_name_consistent <- function(data, dictionary) {
  col_list <- dictionary[, column]
  
  for (col in col_list) {
    temp_names_ls <- dictionary[column == col, names][[1]]
    
    matches <- temp_names_ls %in% names(data)
    
    if (any(matches)) { # if there is a match
      
      data <- setnames(data, temp_names_ls[matches][1], col)
    } else {
      data <- mutate(data, !!col := NA)
    }
  }
  
  return(data)
}

read_data <- function(ffy, var_ls){
  
  data_temp <- here("Data", ffy, "analysis_data.shp") %>% 
                read_sf() 
  select_set <- as.vector(intersect(names(data_temp),var_ls))
  data_temp <- data_temp %>%         
                dplyr::select(select_set) %>%
                filter(!is.na(yild_vl)) %>%
                cbind(., st_coordinates(st_centroid(.))) %>%
                mutate(polygon_area := st_area(.)) %>%
                make_var_name_consistent(., dictionary[type == "final", ])
  return(data_temp)
  
}

find_opt_u <- function(data, res, crop_price, n_price) {
  
  data_dt <- data.table(data)
  
  n_ls <- seq(
    quantile(data_dt$n, prob = 0.025), 
    quantile(data_dt$n, prob = 0.975), 
    length = 100
  )
  
  opt_input_u <- data_dt %>% 
    # this is okay because each zone has the same
    # number of observations
    # unique(by = "zone_txt") %>% 
    .[rep(1:nrow(.), length(n_ls)), ] %>% 
    .[, n := rep(n_ls, each = nrow(.)/length(n_ls))] %>% 
    .[, yield_hat := predict(res, newdata = .)] %>% 
    .[, profit_hat := crop_price * yield_hat - n_price * n] %>% 
    .[, .(profit_hat = mean(profit_hat)), by = n] %>% 
    .[order(profit_hat), ] %>% 
    .[.N, n]
  
  return(opt_input_u)
  
}

run_scam_gam <- function(data, field_vars){
  
  results <- NULL
  
  results <- tryCatch(
    {
      eval_with_timeout(
        {
          formula <- paste0(
            "yield ~ s(n, bs = \"micv\") + s(X, k = 5) + s(Y, k = 5) + te(X, Y, k = c(5, 5))",
            ifelse(
              length(field_vars) > 0,
              paste0(" + ", paste0(field_vars, collapse = " + ")),
              ""
            )
          ) %>% as.formula()
          
          scam(formula, data = data)
          
        },
        timeout = 60, # 60 seconds,
        on_timeout = "silent"
      )
    },
    error = function(cond){
      return(NULL)
    }
  )
  
  
  if (is.null(results)) {
    
    formula <- paste0(
      "yield ~ s(n, k = 3) + s(X, k = 4) + s(Y, k = 4) + te(X, Y, k = c(5, 5))",
      ifelse(
        length(field_vars) > 0,
        paste0(" + ", paste0(field_vars, collapse = " + ")),
        ""
      )
    ) %>% formula()
    
    results <- gam(formula, data = data)
    
  }
  
  return(results) 
  
}

run_gam <- function(data, field_vars){
  
  results <- NULL
  
  formula <- paste0(
    "yield ~ s(n, k = 3) + s(X, k = 4) + s(Y, k = 4) + te(X, Y, k = c(5, 5))",
    ifelse(
      length(field_vars) > 0,
      paste0(" + ", paste0(field_vars, collapse = " + ")),
      ""
    )
  ) %>% formula()
  
  results <- gam(formula, data = data)
  
  return(results) 
  
}

assign_gc_rate <- function(data, input_type, gc_type, gc_rate_n, mrtn, mrtn_min, mrtn_max) {
  
  data_temp <- tryCatch(
    {
      if (gc_type == "Rx") {
        #--------------------------
        # Read Rx data
        #--------------------------
        Rx <- st_read(gc_rate_n) %>% 
          st_set_crs(4326) %>% 
          st_transform(st_crs(data)) %>%
          # st_make_valid() %>%
          setnames(names(.), tolower(names(.)))
        
        dict_input <- dictionary[type == paste0("Rx-", tolower(input_type)), ]
        col_list <- dict_input[, column]
        
        Rx <- make_var_name_consistent(
          Rx, 
          dict_input 
        )
        
        #/*----------------------------------*/
        #' ## Unit conversion
        #/*----------------------------------*/
        if (input_type == "N") {
          Rx <- mutate(Rx, 
                       tgti = convert_N_unit(
                         input_data_n$form, 
                         input_data_n$unit, 
                         tgti, 
                         field_data$reporting_unit
                       ) 
                       # + n_base_rate # add base N rate
          )
        } else if (input_type == "S") {
          #--- seed rate conversion ---#
          if (any(Rx$tgti > 10000)){
            #--- convert to K ---#
            Rx <- mutate(Rx, tgti = tgti / 1000)
          }
        }
        
        #--------------------------
        # Identify grower-chosen rate by observation
        #--------------------------
        obs_tgti <- st_intersection(data, Rx) %>% 
          mutate(area = as.numeric(st_area(.))) %>% 
          data.table() %>% 
          .[, .SD[area == max(area)], by = obs_id] %>% 
          .[, num_obs_per_zone := .N, tgti] %>% 
          .[, analyze := FALSE] %>% 
          .[num_obs_per_zone >= 200, analyze := TRUE] %>% 
          .[, .(obs_id, tgti, analyze)] 
        
        data <- left_join(data, obs_tgti, by = "obs_id") %>% 
          rename(gc_rate_n = tgti)
      }
      
    },
    error = function(cond){
      data$gc_rate_n <- mean(Rx$tgti)
      return(data)
    }
  )
  
  if (gc_type == "uniform") {
    
    data$gc_rate_n <- gc_rate_n 
    
  } else {
    data <- data_temp
    
  }
  
  data$mrtn_rate <- mrtn 
  data$mrtn_min <- mrtn_min
  data$mrtn_max <- mrtn_max
  return(data)
}

make_data_for_eval <- function(data, est) {
  
  data_dt <- data.table(data)
  
  var_names_ls <- est$model %>% 
    data.table() %>% 
    dplyr::select(- any_of(c("n", "s", "yield"))) %>%
    names() 
  
    data_for_eval <- data_dt[, ..var_names_ls] %>% 
      .[, lapply(.SD, mean, na.rm=TRUE)]

  return(data_for_eval)
  
}

predict_yield_range <- function(data_for_eval, n_rate_seq, est) {
  
  eval_data <- cbind(data_for_eval,n_rate_seq)

  #--- predict yield ---#
  yield_prediction <- predict(est, newdata = eval_data, se = TRUE)
  
  eval_data <- eval_data %>% 
    .[, `:=`(
      yield_hat = yield_prediction$fit,
      yield_hat_se = yield_prediction$se.fit
    )] 
  
  return(eval_data)
  
}

predict_yield <- function(n_rate_seq, est){
  
  cov_matrix <- est$resvar[-c(1,2), -c(1,2)]
  
  pred_data <- expand.grid(n=n_rate_seq) %>% 
    data.table() %>% 
    .[,yield_hat:=est$coefficients[1]+est$coefficients[2]*n+est$coefficients[3]*n^2] %>%
    .[,yield_hat_se:=sqrt(cov_matrix[1,1] + n^2*cov_matrix[2,2] + n^4*cov_matrix[3,3] +
                            2*n*cov_matrix[2,1] + 2*n^2*cov_matrix[3,1] + 
                            2*n^3*cov_matrix[2,3])]
  
  return(pred_data)
  
}

predict_yield_qp <- function(n_rate_seq, est) {
  
  #--- predict yield ---#
  pred_data <- expand.grid(n=n_rate_seq) %>% 
    mutate(yield_hat = case_when(n <= -0.5 * est[2,]/est[3,] ~ est[1,] + est[2,] * n + est[3,] * n^2,  
                                 n >  -0.5 * est[2,]/est[3,]  ~ est[1,] -est[2,]^2/(4 * est[3,])
    ))%>%
    data.table()
  
  return(pred_data)
  
}
convert_N_unit <- function(form, unit, rate, reporting_unit, conversion_type = "to_n_equiv") {
  conv_table <-
    fromJSON(
      here("Data","nitrogen_conversion.json"),
      flatten = TRUE
    ) %>%
    data.table() %>%
    .[, conv_factor := as.numeric(conv_factor)] %>%
    .[, form_unit := paste(type, unit, sep = "_")] %>%
    as.data.frame()
  
  if (form == "N_equiv") {
    conv_factor_n <- 1
  } else {
    conv_factor_n <- which(conv_table[, "form_unit"] %in% paste(form, unit, sep = "_")) %>%
      conv_table[., "conv_factor"]
  }
  
  if (reporting_unit == "metric") {
    conv_factor_n <- conv_factor_n * conv_unit(1, "lbs", "kg") * conv_unit(1, "hectare", "acre")
  }
  
  if (conversion_type == "to_n_equiv") {
    converted_rate <- (conv_factor_n) * rate
  } else {
    converted_rate <- (1 / conv_factor_n) * rate
  }
  
  return(as.numeric(converted_rate))
}

get_dif_stat <- function(data, test_var, opt_var, gc_var, gam_res, crop_price, n_price){
  
  if ("scam" %in% class(gam_res)) {
    gam_coef <- gam_res$coefficients.t
    gam_V <- gam_res$Ve.t
  } else {
    gam_coef <- gam_res$coefficients
    gam_V <- gam_res$Ve
  }
  
  base_data <- data.table::copy(data) %>% 
    .[, (test_var) := get(gc_var)]
  
  comp_data <- data.table::copy(data) %>% 
    .[, (test_var) := get(opt_var)]
  
  #/*----------------------------------*/
  #' ## Profit (gc)
  #/*----------------------------------*/
  Xmat_base <- predict(gam_res, newdata = base_data, type = "lpmatrix") 
  # predict(gam_res, newdata = base_data) %>% mean
  
  #--- vector of 1s for summation divided by the number of observations for averaging ---#
  ones <- matrix(1 / dim(Xmat_base)[1], 1, dim(Xmat_base)[1])
  
  #--- average yield ---#
  yhat_base <- ones %*% Xmat_base %*% gam_coef
  
  #--- point estimate of profit differential ---#
  pi_gc <- crop_price * yhat_base - (n_price * ones %*% base_data$n)  
  
  big_mat_base <- ones %*% Xmat_base
  
  #--- se of the profit differential  ---# 
  pi_gc_se <- crop_price * sqrt(big_mat_base %*% gam_V %*% t(big_mat_base))
  
  #/*----------------------------------*/
  #' ## Profit (optimal)
  #/*----------------------------------*/
  Xmat_comp <- predict(gam_res, newdata = comp_data, type = "lpmatrix") 
  
  #--- vector of 1s for summation divided by the number of observations for averaging ---#
  ones <- matrix(1 / dim(Xmat_comp)[1], 1, dim(Xmat_comp)[1])
  
  #--- average yield ---#
  yhat_comp <- ones %*% Xmat_comp %*% gam_coef
  
  #--- point estimate of profit differential ---#
  pi_opt <- crop_price * yhat_comp - (n_price * ones %*% comp_data$n)  
  
  big_mat_comp <- ones %*% Xmat_comp
  
  #--- se of the profit differential  ---# 
  pi_opt_se <- crop_price * sqrt(big_mat_comp %*% gam_V %*% t(big_mat_comp))
  
  #/*----------------------------------*/
  #' ## Profit differential
  #/*----------------------------------*/
  #--- difference in X mat ---#
  X_dif_mat <- Xmat_comp - Xmat_base
  
  #--- vector of 1s for summation divided by the number of observations for averaging ---#
  ones <- matrix(1 / dim(X_dif_mat)[1], 1, dim(X_dif_mat)[1])
  
  #--- X_dif_mat summed ---#
  big_mat_dif <- ones %*% X_dif_mat
  
  #--- point estimate of profit differential ---#
  pi_dif <- ones %*% ((crop_price * X_dif_mat %*% gam_coef) - n_price * (comp_data$n - base_data$n))  
  
  #--- se of the profit differential  ---# 
  pi_dif_se <- crop_price * sqrt(big_mat_dif %*% gam_V %*% t(big_mat_dif))
  
  #--- t-stat ---#
  t_stat <- (pi_dif/pi_dif_se) %>% round(digits = 2)  
  
  return_data <- data.table(
    yhat_est_gc = yhat_base[1, 1],
    point_est_gc = pi_gc[1, 1],
    point_est_gc_se = pi_gc_se[1, 1],
    yhat_est_opt = yhat_comp[1, 1],
    point_est_opt = pi_opt[1, 1],
    point_est_opt_se = pi_opt_se[1, 1],
    point_est_dif = pi_dif[1, 1],
    point_est_dif_se = pi_dif_se[1, 1],
    t = t_stat[1, 1]
  )
  
  return(return_data)
  
}

get_dif_sem <- function(data, test_var, opt_var, gc_var, est, crop_price, n_price){
  
  cov_matrix <- est$resvar[-c(1,2), -c(1,2)]
  
  base_data <- data.table::copy(data) %>% 
    .[, (test_var) := get(gc_var)]
  
  comp_data <- data.table::copy(data) %>% 
    .[, (test_var) := get(opt_var)]
  
  #/*----------------------------------*/
  #' ## Profit (gc)
  #/*----------------------------------*/
  yhat_base <- expand.grid(n=base_data$n) %>% 
    data.table() %>% 
    .[,yield_hat:=est$coefficients[1]+est$coefficients[2]*n+est$coefficients[3]*n^2] 
  
  #--- point estimate of profit differential ---#
  pi_gc <- crop_price * yhat_base$yield_hat - n_price * yhat_base$n
  
  #--- se of the profit differential  ---# 
  pi_gc_se <- sqrt(
    crop_price^2*cov_matrix[1,1]                                     + 
      (crop_price-n_price)^2*yhat_base$n^2*cov_matrix[2,2]             + 
      crop_price^2*yhat_base$n^4*cov_matrix[3,3]                       +
      (crop_price-n_price)*2*yhat_base$n*cov_matrix[2,1]               + 
      crop_price*(crop_price-n_price)*2*yhat_base$n^2*cov_matrix[3,1]  + 
      crop_price^2*2*yhat_base$n^3*cov_matrix[2,3]
  )
  #/*----------------------------------*/
  #' ## Profit (optimal)
  #/*----------------------------------*/
  yhat_comp <- expand.grid(n=comp_data$n) %>% 
    data.table() %>% 
    .[,yield_hat:=est$coefficients[1]+est$coefficients[2]*n+est$coefficients[3]*n^2] 
  
  #--- point estimate of profit differential ---#
  pi_opt <- crop_price * yhat_comp$yield_hat - n_price * yhat_comp$n
  
  #--- se of the profit differential  ---# 
  pi_opt_se <- sqrt(
    crop_price^2*cov_matrix[1,1]                                     + 
      (crop_price-n_price)^2*yhat_comp$n^2*cov_matrix[2,2]             + 
      crop_price^2*yhat_comp$n^4*cov_matrix[3,3]                       +
      (crop_price-n_price)*2*yhat_comp$n*cov_matrix[2,1]               + 
      crop_price*(crop_price-n_price)*2*yhat_comp$n^2*cov_matrix[3,1]  + 
      crop_price^2*2*yhat_comp$n^3*cov_matrix[2,3]
  )
  
  pi_dif <- pi_opt - pi_gc
  
  #--- se of the profit differential  ---# 
  pi_dif_se <- sqrt(
    (crop_price-n_price)^2 * (yhat_comp$n-yhat_base$n)^2 * cov_matrix[2,2] +
      crop_price^2 * (yhat_comp$n^2-yhat_base$n^2)^2 * cov_matrix[3,3]            +
      2*(crop_price-n_price)*crop_price*(yhat_comp$n-yhat_base$n)*(yhat_comp$n^2-yhat_base$n^2) * cov_matrix[2,3]    
  )
  
  #--- t-stat ---#
  t_stat <- (pi_dif/pi_dif_se) %>% round(digits = 2)  
  
  return_data <- data.table(
    yhat_est_gc = mean(yhat_base$yield_hat),
    point_est_gc = mean(pi_gc),
    point_est_gc_se = mean(pi_gc_se),
    yhat_est_opt = mean(yhat_comp$yield_hat),
    point_est_opt = mean(pi_opt),
    point_est_opt_se = mean(pi_opt_se),
    point_est_dif = mean(pi_dif),
    point_est_dif_se = mean(pi_dif_se),
    t = mean(t_stat)
  )
  
  return(return_data)
  
}

get_dif_qp <- function(data, test_var, opt_var, gc_var, est, crop_price, n_price){
  
  base_data <- data.table::copy(data) %>% 
    .[, (test_var) := get(gc_var)]
  
  comp_data <- data.table::copy(data) %>% 
    .[, (test_var) := get(opt_var)]
  
  #/*----------------------------------*/
  #' ## Profit (gc)
  #/*----------------------------------*/
  yhat_base <- expand.grid(n=base_data$n) %>% 
    mutate(yield_hat = case_when(n <= -0.5 * est[2,]/est[3,] ~ est[1,] + est[2,] * n + est[3,] * n^2,  
                                 n >  -0.5 * est[2,]/est[3,]  ~ est[1,] -est[2,]^2/(4 * est[3,])
    )) %>%
    data.table() 
  
  #--- point estimate of profit differential ---#
  pi_gc <- crop_price * yhat_base$yield_hat - n_price * yhat_base$n
  
  #/*----------------------------------*/
  #' ## Profit (optimal)
  #/*----------------------------------*/
  yhat_comp <- expand.grid(n=comp_data$n) %>% 
    mutate(yield_hat = case_when(n <= -0.5 * est[2,]/est[3,] ~ est[1,] + est[2,] * n + est[3,] * n^2,  
                                 n >  -0.5 * est[2,]/est[3,]  ~ est[1,] -est[2,]^2/(4 * est[3,])
    )) %>%
    data.table() 
  
  #--- point estimate of profit differential ---#
  pi_opt <- crop_price * yhat_comp$yield_hat - n_price * yhat_comp$n
  
  #--- profit differential  ---# 
  
  pi_dif <- pi_opt - pi_gc
  
  return_data <- data.table(
    yhat_est_gc = mean(yhat_base$yield_hat),
    point_est_gc = mean(pi_gc),
    yhat_est_opt = mean(yhat_comp$yield_hat),
    point_est_opt = mean(pi_opt),
    point_est_dif = mean(pi_dif)
  )
  
  return(return_data)
  
}



