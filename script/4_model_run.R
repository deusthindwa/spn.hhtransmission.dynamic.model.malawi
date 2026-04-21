#written by Deus
#10/04/2026
#WAIFW household pneumococcal carriage transmission modelling in Malawi

#====================================================================

#settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(1988L)

#age/HIV hroup labels
AGE_LABELS <- c(  "Y-Child", "O-Child", "Adult\nHIV\u2212", "Adult\nHIV+\nART-S", "Adult\nHIV+\nART-L")
AGE_LABELS_FLAT <- c("Y-Child", "O-Child", "Adult HIV-", "Adult HIV+ ART-S", "Adult HIV+ ART-L")
N_AGE <- 5L

for (model_case in 1:4) {
  
  scenario = model_case
  
  #adults with missing hiv status are grouped as hiv+ with longer ART duration in each household (scenario 1)
  if (scenario == 1){
    
    #load & encode data
    dat <- 
      spn_hhF %>% 
      dplyr::rename('agex' = 'agecat' ) %>%
      dplyr::rename('agecat' = 'agecatx' ) %>%
      dplyr::rename('agecatx' = 'agex' ) %>%
      dplyr::mutate(
        state = case_when(
          stg == 1L ~ 1L,
          stg == 2L ~ 2L,
          stg == 3L ~ 3L,
          TRUE      ~ NA_integer_)) %>%
      dplyr::mutate(
        age_grp = case_when(
          agecat == "ychild"                                    ~ 1L,
          agecat == "ochild"                                    ~ 2L,
          agecat == "adult" & hiv == "hiv-"                     ~ 3L,
          agecat == "adult" & hiv == "hiv+_artS"                ~ 4L,
          agecat == "adult" & (is.na(hiv) | hiv == "hiv+_artL") ~ 5L,
          TRUE ~ NA_integer_)) %>%
      dplyr::filter(!is.na(state), !is.na(age_grp))
    
    #adults with missing hiv status are grouped as hiv- in each household (scenario 2)
  } else if (scenario == 2){
    #load & encode data
    dat <- 
      spn_hhF %>% 
      dplyr::rename('agex' = 'agecat' ) %>%
      dplyr::rename('agecat' = 'agecatx' ) %>%
      dplyr::rename('agecatx' = 'agex' ) %>%
      dplyr::mutate(
        state = case_when(
          stg == 1L ~ 1L,
          stg == 2L ~ 2L,
          stg == 3L ~ 3L,
          TRUE      ~ NA_integer_)) %>%
      dplyr::mutate(
        age_grp = case_when(
        agecat == "ychild"                                ~ 1L,
        agecat == "ochild"                                ~ 2L,
        agecat == "adult" & (is.na(hiv) | hiv == "hiv-")  ~ 3L,
        agecat == "adult" & hiv == "hiv+_artS"            ~ 4L,
        agecat == "adult" & hiv == "hiv+_artL"            ~ 5L,
        TRUE ~ NA_integer_)) %>%
        dplyr::filter(!is.na(state), !is.na(age_grp))
      
      #adults with missing hiv status are grouped as hiv+ with short ART duration in each household (scenario 3)
  }  else if (scenario == 3) {
    #load & encode data
    dat <- 
      spn_hhF %>% 
      dplyr::rename('agex' = 'agecat' ) %>%
      dplyr::rename('agecat' = 'agecatx' ) %>%
      dplyr::rename('agecatx' = 'agex' ) %>%
      dplyr::mutate(
        state = case_when(
          stg == 1L ~ 1L,
          stg == 2L ~ 2L,
          stg == 3L ~ 3L,
          TRUE      ~ NA_integer_)) %>%
      dplyr::mutate(
      age_grp = case_when(
        agecat == "ychild"                                    ~ 1L,
        agecat == "ochild"                                    ~ 2L,
        agecat == "adult" & hiv == "hiv-"                     ~ 3L,
        agecat == "adult" & (is.na(hiv) | hiv == "hiv+_artS") ~ 4L,
        agecat == "adult" & hiv == "hiv+_artL"                ~ 5L,
        TRUE ~ NA_integer_)) %>%
        dplyr::filter(!is.na(state), !is.na(age_grp))
      
    #adults with missing hiv status have the same HIV status as HH index (baseline)
  } else { #(scenario == 4)
    #load & encode data
    dat <- 
      spn_hhF %>% 
      dplyr::rename('agex' = 'agecat' ) %>%
      dplyr::rename('agecat' = 'agecatx' ) %>%
      dplyr::rename('agecatx' = 'agex' ) %>%
      dplyr::mutate(
        state = case_when(
          stg == 1L ~ 1L,
          stg == 2L ~ 2L,
          stg == 3L ~ 3L,
          TRUE      ~ NA_integer_)) %>%
      dplyr::group_by(hhid) %>%
      dplyr::mutate(hh_hiv_adult = hiv[which(!is.na(hiv) & agecat %in% c("adult"))][1], #extract a single non-missing  HIV value per household (if any)
                    hiv = if_else(is.na(hiv) & agecat %in% c("adult"), hh_hiv_adult, hiv)) %>% #fill missing HIV for adults only
      dplyr::ungroup() %>%
      dplyr::select(-hh_hiv_adult) %>%
      dplyr::mutate(
        age_grp = case_when(
          agecat == "ychild"                     ~ 1L,
          agecat == "ochild"                     ~ 2L,
          agecat == "adult" & hiv == "hiv-"      ~ 3L,
          agecat == "adult" & hiv == "hiv+_artS" ~ 4L,
          agecat == "adult" & hiv == "hiv+_artL" ~ 5L,
          TRUE ~ NA_integer_)) %>%
      dplyr::filter(!is.na(state), !is.na(age_grp))
  }
  
  #ge/HIV group sizes across all visits
  dat %>%
    distinct(pid, age_grp) %>%
    count(age_grp) %>%
    mutate(label = AGE_LABELS_FLAT[age_grp]) %>%
    print()
  
  #arriage state distribution by group
  dat %>%
    dplyr::count(age_grp, state) %>%
    dplyr::mutate(
      group = AGE_LABELS_FLAT[age_grp],
      state_label = c("S", "VT", "NVT")[state]
    ) %>%
    dplyr::select(group, state_label, n) %>%
    tidyr::pivot_wider(names_from = state_label, values_from = n,
                       values_fill = 0L) %>%
    print()
  
  #household composition by age group at each visit
  long <- 
    dat %>% 
    dplyr::select(hhid, pid, vno, state, age_grp)
  
  hh_raw <- 
    long %>%
    dplyr::inner_join(
      long %>% dplyr::rename(pid_o = pid, state_o = state, age_o = age_grp),
      by = c("hhid", "vno"),
      relationship = "many-to-many") %>%
    dplyr::filter(pid != pid_o)
  
  hh_comp <- 
    hh_raw %>%
    dplyr::group_by(hhid, pid, vno) %>%
    dplyr::summarise(
      n_VT_1  = sum(state_o == 2L & age_o == 1L),
      n_VT_2  = sum(state_o == 2L & age_o == 2L),
      n_VT_3  = sum(state_o == 2L & age_o == 3L),
      n_VT_4  = sum(state_o == 2L & age_o == 4L),
      n_VT_5  = sum(state_o == 2L & age_o == 5L),
      n_NVT_1 = sum(state_o == 3L & age_o == 1L),
      n_NVT_2 = sum(state_o == 3L & age_o == 2L),
      n_NVT_3 = sum(state_o == 3L & age_o == 3L),
      n_NVT_4 = sum(state_o == 3L & age_o == 4L),
      n_NVT_5 = sum(state_o == 3L & age_o == 5L),
      n_HH_others = n(),
      .groups = "drop")
  
  hh_comp_full <- 
    long %>%
    dplyr::select(hhid, pid, vno) %>%
    dplyr::distinct() %>%
    dplyr::left_join(hh_comp, by = c("hhid", "pid", "vno")) %>%
    tidyr::replace_na(as.list(setNames(
      rep(0L, 11),
      c(paste0("n_VT_", 1:5), paste0("n_NVT_", 1:5), "n_HH_others"))))
  
  #community prevalence by age group & visit
  comm_prev <- 
    dat %>%
    dplyr::filter(vno %in% c(1L, 2L)) %>%
    dplyr::group_by(vno, age_grp) %>%
    dplyr::summarise(prev_VT  = mean(state == 2L),
                     prev_NVT = mean(state == 3L),
                     .groups  = "drop")
  
  build_comm_mat <- function(col) {
    mat <- matrix(0.0, N_AGE, 2L)
    for (b in seq_len(N_AGE)) {
      for (t in 1:2) {
        tmp <- comm_prev %>% dplyr::filter(age_grp == b, vno == t)
        if (nrow(tmp) > 0) mat[b, t] <- pull(tmp, all_of(col))
      }
    }
    mat
  }
  
  comm_VT_mat  <- build_comm_mat("prev_VT")
  comm_NVT_mat <- build_comm_mat("prev_NVT")
  
  #community VT prevalence [group x visit]
  print(round(`rownames<-`(`colnames<-`(comm_VT_mat, c("Visit1","Visit2")),
                           AGE_LABELS_FLAT), 3))
  #community NVT prevalence [group x visit]
  print(round(`rownames<-`(`colnames<-`(comm_NVT_mat, c("Visit1","Visit2")),
                           AGE_LABELS_FLAT), 3))
  
  rownames(comm_VT_mat) <- rownames(comm_NVT_mat) <- NULL
  colnames(comm_VT_mat) <- colnames(comm_NVT_mat) <- NULL
  
  #reshape & merge HH covariates
  dat_wide <- 
    dat %>%
    dplyr::select(hhid, pid, vno, state, age_grp) %>%
    tidyr::pivot_wider(names_from = vno, values_from = state, names_prefix = "s") %>%
    dplyr::filter(!is.na(s1), !is.na(s2), !is.na(s3))
  
  get_hh_at_visit <- function(v) {
    hh_comp_full %>%
      dplyr::filter(vno == v) %>%
      dplyr::select(-vno) %>%
      dplyr::rename_with(~ paste0(.x, "_v", v), starts_with("n_"))
  }
  
  dat_model <- 
    dat_wide %>%
    dplyr::left_join(get_hh_at_visit(1), by = c("hhid", "pid")) %>%
    dplyr::left_join(get_hh_at_visit(2), by = c("hhid", "pid"))
  
  N <- nrow(dat_model)
  cat(sprintf("\nmodel dataset: %d individuals, 2 transitions, %d observations\n", N, 2L * N))
  
  #build Stan 3-D Arrays [N x N_AGE x 2
  n_VT_arr  <- array(0L, dim = c(N, N_AGE, 2L))
  n_NVT_arr <- array(0L, dim = c(N, N_AGE, 2L))
  for (b in seq_len(N_AGE)) {
    for (t in 1:2) {
      n_VT_arr[, b, t]  <- as.integer(dat_model[[paste0("n_VT_",  b, "_v", t)]])
      n_NVT_arr[, b, t] <- as.integer(dat_model[[paste0("n_NVT_", b, "_v", t)]])
    }
  }
  
  n_HH_others_mat <- cbind(as.integer(dat_model$n_HH_others_v1), as.integer(dat_model$n_HH_others_v2))
  
  stan_data <- list(
    N           = N,
    N_AGE       = N_AGE,
    s1          = as.integer(dat_model$s1),
    s2          = as.integer(dat_model$s2),
    s3          = as.integer(dat_model$s3),
    age_grp     = as.integer(dat_model$age_grp),
    n_VT        = n_VT_arr,
    n_NVT       = n_NVT_arr,
    n_HH_others = n_HH_others_mat,
    comm_VT     = comm_VT_mat,
    comm_NVT    = comm_NVT_mat
  )
  
  # #call Stan model
  # source(here("script", "3_stan_model.R"))
  # 
  # #compile Stan Model
  # cat("compiling Stan model\n")
  # spn_model <- stan_model(
  #   model_code = stan_code,
  #   model_name = "spn_waifw_5grp_foicomp",
  #   verbose    = FALSE
  # )
  # 
  # #fit via HMC NUTS
  # cat("\nrunning HMC sampler (4 chains x 5000 iterations)\n")
  # fit <- sampling(
  #   object  = spn_model,
  #   data    = stan_data,
  #   chains  = 4L,
  #   iter    = 5000L,
  #   warmup  = 1000L,
  #   seed    = 2026L,
  #   control = list(adapt_delta = 0.95, max_treedepth = 12L),
  #   verbose = TRUE
  # )
  # 
  # if (scenario == 1){
  #   #save scenario 1 model
  #   rio::export(fit, here::here("results", "fit_hhbayes_longART.rds"))
  #   
  # } else if (scenario == 2){
  #   #save scenario 2 model
  #   rio::export(fit, here::here("results", "fit_hhbayes_negHIV.rds"))
  #   
  # } else if (scenario == 3){
  #   #save scenario 3 model
  #   rio::export(fit, here::here("results", "fit_hhbayes_shortART.rds"))
  #   
  # } else{
  #   #save baseline model
  #   rio::export(fit, here::here("results", "fit_hhbayes_main.rds"))
  #   
  # }
  
}
