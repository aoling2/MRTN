#--- get field parameters for the field-year ---#
w_field_data <- field_data[field_year == ffy, ]

w_farm_field <- w_field_data$farm_field
w_year <- w_field_data$year

#--- get input data ---#
input_data <- dplyr::select(w_field_data, starts_with(
  "input")) %>%  map(1) %>% 
  rbindlist(fill = TRUE)

#/*----------------------------------*/
#' ## Crop information
#/*----------------------------------*/
crop <- w_field_data[, crop] 
crop_unit <- w_field_data[, crop_unit] 
crop_price <- as.numeric(w_field_data[, crop_price]) 
if(is.na(crop_price) == TRUE) {
  crop_price <- case_when(
    crop == "soy" ~ 14, # $/bu
    crop == "corn" ~ 5.5, # $/bu
    crop == "cotton" ~ 0.93 * 480, # $/bales
  )
}
land_unit <- w_field_data[, land_unit] 
reporting_unit <- w_field_data[, reporting_unit] 
harvester_width <- w_field_data[, h_width][[1]]

#/*----------------------------------*/
#' ## Input information
#/*----------------------------------*/
#/*~~~~~~~~~~~~~~~~~~~~~~*/
#' ### Seed
#/*~~~~~~~~~~~~~~~~~~~~~~*/
#--- should we process seed data? ---#
process_s <- "seed" %in% input_data[strategy == "trial", form]

if (process_s) {

  input_data_s <- input_data[strategy == "trial" & form == "seed", ]
  use_td_s <- input_data_s[, use_target_rate_instead]

  #++++++++++++++++
  # grower chosen rate and analysis type 
  #++++++++++++++++
  gc_rate_s <- input_data_s[, sq_rate]

  if(!is.numeric(gc_rate_s) | is.na(gc_rate_s)) {
      
    Rx_file_s <- file.path(
      here("Data/Growers", ffy, "Raw"), 
      paste0(gc_rate_s, ".shp")
    )

    if (file.exists(Rx_file_s)){
      #--- if the Rx file exists ---#
      gc_type_s <- "Rx"
      gc_rate_s <- Rx_file_s

    } else {
      #--- if the Rx file doe NOT exist ---#
      # default rate
      gc_rate_s <- case_when(
        crop == "corn" ~ 36,
        crop == "soy" ~ 120
      )
      gc_type_s <- "uniform"
    }
  } else {
    #--- seed rate conversion ---#
    if (gc_rate_s > 10000){
      #--- convert to K ---#
      gc_rate_s <- gc_rate_s / 1000
    }

    gc_type_s <- "uniform"
  }

  #++++++++++++++++
  # seed price
  #++++++++++++++++
  s_price <- as.numeric(input_data_s[, price])

  if(is.na(s_price) == TRUE) {
    #=== default price ===#
    if (crop == "corn") {
      s_price <- 0.00381 * 1000 # (thousand seed)
    } else if(crop == "soy") {
      s_price <- 0.000321 * 1000 # (thousand seed)
    }
  } else { # if price is numeric
    if (s_price < 0.01) {
      s_price <- s_price * 1000
    }
  }

  #++++++++++++++++
  # planter width
  #++++++++++++++++
  m_width_s <- input_data_s[, machine_width]
  
} else { # if not seed trial

  gc_type_s <- NA
  gc_rate_s <- NA
  s_price <- NA
  use_td_s <- NA
  planter_width <- NA
  m_width_s <- NA
}


#/*~~~~~~~~~~~~~~~~~~~~~~*/
#' ### Base nitrogen 
#/*~~~~~~~~~~~~~~~~~~~~~~*/
is_base_N <- "base" %in% input_data[, strategy]

if (is_base_N) {

  n_base_rate <- input_data[strategy == "base", ] %>% 
    rowwise() %>% 
    mutate(
      n_equiv_rate = convert_N_unit(
        form = form, 
        unit = unit, 
        rate = rate, 
        reporting_unit = reporting_unit
      ) 
    ) %>% 
    data.table() %>% 
    .[, sum(n_equiv_rate)]

} else {

  n_base_rate <- 0  

}

#/*~~~~~~~~~~~~~~~~~~~~~~*/
#' ### Nitrogen
#/*~~~~~~~~~~~~~~~~~~~~~~*/
#--- should we process N data? ---#
n_var_ls <- c("NH3", "urea", "uan32", "uan28", "1_2_1(36)", "LAN(26)", "MAP", "1_0_0", "1_0_1", "2_3_2(22)",
              "15_10_6", "3_0_1", "2_3_4(32)", "4_3_4(33)", "5_1_5", "Sp")
process_n_idv <- n_var_ls %in% input_data[strategy == "trial", form]
process_n <- any(process_n_idv) 

if (process_n) {
  n_var <- n_var_ls[process_n_idv]
  input_data_n <- input_data[strategy == "trial" & form == n_var, ]
  use_td_n <- input_data_n[, use_target_rate_instead]
  
  
  #++++++++++++++++
  # grower chosen rate and analysis type 
  #++++++++++++++++
  gc_rate_n <- input_data_n[, sq_rate]

  if (!is.numeric(gc_rate_n)) {
      
    Rx_file_n <- file.path(
      here("Data", ffy), 
      paste0(gc_rate_n, ".shp")
    )

    if (file.exists(Rx_file_n)){
      #--- if the Rx file exists ---#
      gc_type_n <- "Rx"
      gc_rate_n <- Rx_file_n

    } else {
      #--- if the Rx file doe NOT exist ---#
      # default rate
      gc_rate_n <- 180
      gc_type_n <- "uniform"
    }
  } else {

    gc_rate_n <- convert_N_unit(
      input_data_n[, form], 
      input_data_n[, unit], 
      gc_rate_n, 
      reporting_unit
    ) + 
    n_base_rate

    gc_type_n <- "uniform"

  }

  #++++++++++++++++
  # N price
  #++++++++++++++++
  if ("price" %in% names(input_data_n)) {
    n_price <- as.numeric(input_data_n[, price])
  } else {
    n_price <- NA
  }

  if(is.na(n_price) == TRUE) {
    n_price <- 0.4
  }

  #++++++++++++++++
  # applicator width
  #++++++++++++++++
  m_width_n <- input_data_n[, machine_width]

} else {

  n_var <- NA
  gc_type_n <- NA
  gc_rate_n <- NA
  n_price <- NA
  use_td_n <- NA
  m_width_n <- NA
 
}

#/*~~~~~~~~~~~~~~~~~~~~~~*/
#' ### K
#/*~~~~~~~~~~~~~~~~~~~~~~*/
#--- should we process N data? ---#
process_k <- "K" %in% input_data[strategy == "trial", form]

if (process_k) {

  input_data_k <- input_data[strategy == "trial" & form == "K", ]
  use_td_k <- input_data_k[, use_target_rate_instead]

  #++++++++++++++++
  # grower chosen rate and analysis type 
  #++++++++++++++++
  gc_rate_k <- input_data_k[, sq_rate]

  if (!is.numeric(gc_rate_k)) {
      
    Rx_file_k <- file.path(
      here("Data", ffy), 
      paste0(gc_rate_k, ".shp")
    )

    if (file.exists(Rx_file_k)){
      #--- if the Rx file exists ---#
      gc_type_k <- "Rx"
      gc_rate_k <- Rx_file_k

    } else {
      #--- if the Rx file doe NOT exist ---#
      # default rate
      gc_rate_k <- 40
      gc_type_k <- "uniform"
    }
  } else {

    gc_type_k <- "uniform"

  }

  #++++++++++++++++
  # k price
  #++++++++++++++++
  if ("price" %in% names(input_data_k)) {
    k_price <- as.numeric(input_data_k[, price])
  } else {
    k_price <- NA
  }

  if(is.na(k_price) == TRUE) {
    k_price <- 0.4
  }

  #++++++++++++++++
  # applicator width
  #++++++++++++++++
  m_width_k <- input_data_k[, machine_width]

} else {

  gc_type_k <- NA
  gc_rate_k <- NA
  k_price <- NA
  use_td_k <- NA
  m_width_k <- NA
}

#/*----------------------------------*/
#' ## Trial type
#/*----------------------------------*/
trial_type <- 
case_when(
  process_n & process_s ~ "SN",
  process_n & !process_s ~ "N",
  !process_n & process_s ~ "S",
  process_k ~ "K" 
)

trial_info <- 
tibble(
  crop = crop, 
  input_type = c("S", "N", "K"),
  process = c(process_s, process_n, process_k),
  use_td = c(use_td_s, use_td_n, use_td_k),
  price = c(s_price, n_price, k_price),
  gc_rate = list(gc_rate_s, gc_rate_n, gc_rate_k),
  gc_type = c(gc_type_s, gc_type_n, gc_type_k),
  machine_width = c(m_width_s, m_width_n, m_width_k)
) %>% 
filter(process == TRUE)
