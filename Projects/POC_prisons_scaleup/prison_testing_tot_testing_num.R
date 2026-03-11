
res_2024 <- HCVMSM(POC_AU, best_estimates, best_est_pop,
       disease_progress, pop_array,
       dfList,  
       param_cascade_sc = dfList_NP_2024, 
       fib = fib, 
       modelrun = "UN", proj = "POC_AU", end_Y = endY, 
       cost = NULL, costflow = NULL, 
       costflow_Neg = NULL, 
       fc = fs[["2024"]])

test_compartment <- popResults_MidYear(POC_AU, res_2024,
                   Population = POC_AU$popNames,
                   Disease_prog = POC_AU$diseaseprogress_Name, 
                   Cascade = POC_AU$cascade_name, param = NULL, 
                   endYear = 100)%>%ungroup()%>%
  mutate(diag_status = ifelse(cascade %in% c("s", "cured") & disease_prog!= "a", "S", 
                              ifelse(cascade %in% c("undiag") &disease_prog != "a", "U", "diag")),
         year = year + POC_AU$cabY - 1)%>%
  group_by(population, year, diag_status)%>%
  summarise(best = sum(best))%>%filter(year== 2025)
View(test_compartment)
fm[["2024"]] <- c(1.1, 1.1, 18.5, 18.5, 1)

12000/(72180*0.1457498+ 2218*1.1 + 336729*0.1457498 + 5866*1.1)
13000/(5755*2*0.1049667+ 133*2*18.5 + 11237*2*0.1049667 + 59*18.5 + 20287*4*1) 




#### for 2025 #### 
odd_num_test <- 1

Ccal[[2025]] <- list("C" = 11000/(72180*0.1457498+ 2218*1.1 + 336729*0.1457498 + 5866*1.1),
                     "P" = 14000/(6191*2*0.1049667+ 618*2*18.5 + 12346*2*0.1049667 + 236*18.5 + 20311*4*1)*0.98)
frac_test[[2025]] <- frac_test[[2024]]

frac_ab[["2025"]] <- c(unlist(as.numeric(frac_test[[2025]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2025]]$P$reflex)))

dfList_NP_2025 <- dfList_NP_2024
for(i in param_var){  
  dfList_NP_2025[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2025, index = i, 
                                   frac_testing = frac_test[[2025]],
                                   S_Yint = 2025, S_Yend = 2026, r_Yend = 2026, NPlst = NPlst, 
                                   fp = c(Ccal[[2025]]$C*fm[["2025"]][1], 
                                          Ccal[[2025]]$C*fm[["2025"]][2], 
                                          Ccal[[2025]]$P*fm[["2025"]][3], 
                                          Ccal[[2025]]$P*fm[["2025"]][4], 
                                          Ccal[[2025]]$P*fm[["2025"]][5]))
  
}

for(i in param_var){
  # begining of 2025
  b_pt <- (2026 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2025[[i]])[3]
  dfList_NP_2025[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2025"]] <- fs[["2024"]]
fs[["2025"]][1, (ini_dt + 1/POC_AU$timestep): (end_dt + 1/POC_AU$timestep)] <-   fs[["2024"]][1, ini_dt: end_dt]
fs[["2025"]][2, (ini_dt + 1/POC_AU$timestep): (end_dt + 1/POC_AU$timestep)] <-   fs[["2024"]][2, ini_dt: end_dt]
fs[["2025"]][3, (ini_dt + 1/POC_AU$timestep): (end_dt + 1/POC_AU$timestep)] <-   fs[["2024"]][3, ini_dt: end_dt]
fs[["2025"]][4, (ini_dt + 1/POC_AU$timestep): (end_dt + 1/POC_AU$timestep)] <-   fs[["2024"]][4, ini_dt: end_dt]
fs[["2025"]][5, (ini_dt + 1/POC_AU$timestep): (end_dt + 1/POC_AU$timestep)] <-   fs[["2024"]][5, ini_dt: end_dt]


test_t <- HCVMSM(POC_AU, best_estimates, best_est_pop,
               disease_progress, pop_array,
               dfList,  
               param_cascade_sc = dfList_NP_2025, 
               fib = fib, 
               modelrun = "UN", proj = "POC_AU", end_Y = endY, 
               cost = NULL, costflow = NULL, 
               costflow_Neg = NULL, 
               fc = fs[["2025"]])




test_compartment <- popResults_MidYear(POC_AU, test_t,
                                       Population = POC_AU$popNames,
                                       Disease_prog = POC_AU$diseaseprogress_Name, 
                                       Cascade = POC_AU$cascade_name, param = NULL, 
                                       endYear = 100)%>%ungroup()%>%
  mutate(diag_status = ifelse(cascade %in% c("s", "cured") & disease_prog!= "a", "S", 
                              ifelse(cascade %in% c("undiag") &disease_prog != "a", "U", "diag")),
         year = year + POC_AU$cabY - 1)%>%
  group_by(population, year, diag_status)%>%
  summarise(best = sum(best))%>%filter(year== 2025)

View(test_compartment)

Ccal[[2025]]
#### for 2026 #### 
odd_num_test <- 1
fm[["2026"]] <- fm[["2025"]]
Ccal[[2026]] <- list("C" = 11000/(73311*0.1457498+ 1849*1.1 + 340119*0.1457498 + 2999*1.1),
                     "P" = 14000/(5980*2*0.1049667+ 104*2*15 + 11393*2*0.1049667 + 44*15 + 20287*4*1))
frac_test[[2026]] <- frac_test[[2025]]

frac_ab[["2026"]] <- c(unlist(as.numeric(frac_test[[2026]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2026]]$P$reflex)))

dfList_NP_2026 <- dfList_NP_2025
for(i in param_var){  
  dfList_NP_2026[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2025, index = i, 
                                   frac_testing = frac_test[[2026]],
                                   S_Yint = 2026, S_Yend = 2027, r_Yend = 2027, NPlst = NPlst, 
                                   fp = c(Ccal[[2026]]$C*fm[["2026"]][1], 
                                          Ccal[[2026]]$C*fm[["2026"]][2], 
                                          Ccal[[2026]]$P*fm[["2026"]][3], 
                                          Ccal[[2026]]$P*fm[["2026"]][4], 
                                          Ccal[[2026]]$P*fm[["2026"]][5]))
  
}

for(i in param_var){
  # begining of 2026
  b_pt <- (2027 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2026[[i]])[3]
  dfList_NP_2026[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 
ini_dt <- (2025 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2025 + 1 ) - POC_AU$cabY)/POC_AU$timestep
fs[["2026"]] <- fs[["2025"]]
fs[["2026"]][1, (ini_dt + 1/POC_AU$timestep ): (end_dt + 1/POC_AU$timestep )] <-   fs[["2025"]][1, ini_dt: end_dt]
fs[["2026"]][2, (ini_dt + 1/POC_AU$timestep ): (end_dt + 1/POC_AU$timestep )] <-   fs[["2025"]][2, ini_dt: end_dt]
fs[["2026"]][3, (ini_dt + 1/POC_AU$timestep ): (end_dt + 1/POC_AU$timestep )] <-   fs[["2025"]][3, ini_dt: end_dt]
fs[["2026"]][4, (ini_dt + 1/POC_AU$timestep ): (end_dt + 1/POC_AU$timestep )] <-   fs[["2025"]][4, ini_dt: end_dt]
fs[["2026"]][5, (ini_dt + 1/POC_AU$timestep ): (end_dt + 1/POC_AU$timestep )] <-   fs[["2025"]][5, ini_dt: end_dt]


test_t <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                 disease_progress, pop_array,
                 dfList,  
                 param_cascade_sc = dfList_NP_2026, 
                 fib = fib, 
                 modelrun = "UN", proj = "POC_AU", end_Y = endY, 
                 cost = NULL, costflow = NULL, 
                 costflow_Neg = NULL, 
                 fc = fs[["2026"]])






test <- list()
cl_ext <- names(test_t)[c(10:22)]
for(i in cl_ext){
  
  test[[i]] <- modres.flow.t(POC_AU, test_t, endYear = 100, 
                             allp = i)%>%
    ungroup()%>%
    group_by(year, population)%>%
    summarise(best = sum(best))
  
  
}

test_sq <- list()
for(i in cl_ext){
  
  test_sq[[i]] <- modres.flow.t(POC_AU, Sce_sq, endYear = 100, 
                                allp = i)%>%
    ungroup()%>%
    group_by(year, population)%>%
    summarise(best = sum(best))
  
  
}


test <- dplyr::bind_rows(test, .id = 'index')%>%group_by(year, population)%>%spread(index, best)
test_sq <- dplyr::bind_rows(test_sq, .id = 'index')%>%group_by(year, population)%>%spread(index, best)


test_fscal <- test%>%mutate(Ab = newTestingAb_sc + newTestingAb_sc_neg, 
                            RNA = (newTestingAg_sc+ newTestingAg_sc_neg + newTestingPOCT_sc + 
                                     newTestingPOCT_sc_neg ))%>%select(year, population, Ab, 
                                                                       RNA)%>%
  mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), "C", "P"))%>%
  ungroup()%>%
  select(-c(population))%>%
  gather(index, value, -c(year, setting))%>%
  group_by(year, setting, index)%>%summarise(value = sum(value))%>%
  mutate(scenario = "prisons_testing_I")
dtp <- c(6667,6667, 11476,11476,20393, 20393, 12000, 13000, 11000, 14000, 10000, 15000,9000, 16000,8000, 17000,7000, 18000 )
# View(test_fscal%>%filter(year%in% c(7:15)))
View(test_fscal%>%filter(year%in% c(7:15))%>%group_by(year, setting)%>%
       summarise(tot_NP_test = sum(value))%>%
       mutate(year = year + 2015))

################################################################################
#                        test 
################################################################################

calibrate_NP_year <- function(
    cal_year, prev_dfList_NP, prev_fs, prev_model, dfList_NP_base,
    pj, Ccal, fm, frac_test, NPlst, n_ab_np, frac_ab, param_var,
    best_estimates, best_est_pop, disease_progress,
    pop_array, dfList, fib, endY,
    target_C     = NULL,
    target_P     = NULL,
    adj_factor_C = NULL,
    adj_factor_P = NULL
) {
  yr_chr    <- as.character(cal_year)
  steps     <- 1 / pj$timestep
  ini_dt    <- (cal_year - pj$cabY) / pj$timestep + 1
  end_dt    <- ((cal_year + 1) - pj$cabY) / pj$timestep
  target_yr <- cal_year - pj$cabY
  
  # ── Input guards ────────────────────────────────────────────────────────────
  if (is.null(n_ab_np[[yr_chr]]))     stop(paste("n_ab_np[[", yr_chr, "]] is NULL"))
  if (is.null(fm[[yr_chr]]))          stop(paste("fm[[", yr_chr, "]] is NULL"))
  if (is.null(frac_test[[cal_year]])) stop(paste("frac_test[[", cal_year, "]] is NULL"))
  
  # ── Helper: extract per-population S and U ──────────────────────────────────
  extract_pop_US <- function(model, at_year) {
    popResults_MidYear(
      pj, model,
      Population   = pj$popNames,
      Disease_prog = pj$diseaseprogress_Name,
      Cascade      = pj$cascade_name,
      param        = NULL, endYear = 100
    ) %>%
      ungroup() %>%
      mutate(
        cascade      = as.character(cascade),
        disease_prog = as.character(disease_prog),
        diag_status  = ifelse(cascade %in% c("s","cured") & disease_prog != "a", "S",
                              ifelse(cascade %in% c("undiag")    & disease_prog != "a", "U", "diag")),
        cal_yr       = year + pj$cabY - 1
      ) %>%
      filter(cal_yr == at_year & diag_status %in% c("U","S")) %>%
      group_by(population, diag_status) %>%
      summarise(best = sum(best), .groups = "drop") %>%
      pivot_wider(names_from = diag_status, values_from = best, values_fill = 0) %>%
      mutate(population = as.character(population),
             S = as.numeric(S), U = as.numeric(U))
  }
  
  get_SU <- function(pop_data, pop_name) {
    row <- pop_data %>% filter(population == pop_name)
    if (nrow(row) == 0) return(list(S = 0, U = 0))
    list(S = as.numeric(row$S), U = as.numeric(row$U))
  }
  
  # ── Helper: get total testing numbers from a model run ──────────────────────
  get_test_totals <- function(model) {
    cl_ext <- c("newTestingAb_sc","newTestingAg_sc","newTestingPOCT_sc",
                "newTestingAb_sc_neg","newTestingAg_sc_neg","newTestingPOCT_sc_neg")
    tflow <- list()
    for (i in cl_ext) {
      tflow[[i]] <- modres.flow.t(pj, model, endYear = 100, allp = i) %>%
        ungroup() %>% group_by(year, population) %>%
        summarise(best = sum(best), .groups = "drop")
    }
    dplyr::bind_rows(tflow, .id = "index") %>%
      group_by(year, population) %>%
      spread(index, best) %>%
      filter(year == target_yr) %>%
      mutate(
        Ab  = rowSums(cbind(newTestingAb_sc,  newTestingAb_sc_neg),  na.rm = TRUE),
        RNA = rowSums(cbind(newTestingAg_sc,  newTestingAg_sc_neg,
                            newTestingPOCT_sc, newTestingPOCT_sc_neg), na.rm = TRUE),
        setting = ifelse(population %in% c("C_PWID","C_fPWID"), "C", "P")
      ) %>%
      group_by(setting) %>%
      summarise(tot = sum(Ab + RNA, na.rm = TRUE), .groups = "drop")
  }
  
  # ── Core run function ────────────────────────────────────────────────────────
  run_model <- function(adj_P, adj_C = 1.0) {
    
    pop_prev <- extract_pop_US(prev_model, at_year = cal_year)
    su_prev  <- lapply(pj$popNames, function(p) get_SU(pop_prev, p))
    names(su_prev) <- pj$popNames
    
    # fs_prev_val from previous year
    fs_prev_val <- as.numeric(prev_fs[, ini_dt - 1])
    if (any(!is.finite(fs_prev_val)) | any(fs_prev_val == 0))
      fs_prev_val <- as.numeric(prev_fs[, ini_dt])
    
    # Compute denoms using fs_prev_val + adj factors
    denom_C <- (su_prev[["C_PWID"]]$S  * fs_prev_val[1] +
                  su_prev[["C_PWID"]]$U  * fm[[yr_chr]][1] +
                  su_prev[["C_fPWID"]]$S * fs_prev_val[2] +
                  su_prev[["C_fPWID"]]$U * fm[[yr_chr]][2]) * adj_C
    
    denom_P <- (su_prev[["P_PWID"]]$S  * fs_prev_val[3] +
                  su_prev[["P_PWID"]]$U  * fm[[yr_chr]][3] +
                  su_prev[["P_fPWID"]]$S * fs_prev_val[4] +
                  su_prev[["P_fPWID"]]$U * fm[[yr_chr]][4] +
                  (su_prev[["P_nPWID"]]$S + su_prev[["P_nPWID"]]$U) * fm[[yr_chr]][5]) *
      (1/pj$timestep * adj_P)
    
    if (!is.finite(denom_C) | denom_C == 0) stop(paste("denom_C is", denom_C))
    if (!is.finite(denom_P) | denom_P == 0) stop(paste("denom_P is", denom_P))
    
    Ccal_yr <- list(
      C = as.numeric(n_ab_np[[yr_chr]][1]) / denom_C,
      P = as.numeric(n_ab_np[[yr_chr]][2]) / denom_P
    )
    
    fp_vec <- c(Ccal_yr$C * fm[[yr_chr]][1],
                Ccal_yr$C * fm[[yr_chr]][2],
                Ccal_yr$P * fm[[yr_chr]][3],
                Ccal_yr$P * fm[[yr_chr]][4],
                Ccal_yr$P * fm[[yr_chr]][5])
    
    # Build dfList_NP
    dfList_NP_year <- prev_dfList_NP
    for (i in param_var) {
      dfList_NP_year[[i]] <- Param_cal(
        pj = pj, dlist = dfList_NP_year, index = i,
        frac_testing = frac_test[[cal_year]],
        S_Yint = cal_year, S_Yend = cal_year + 1, r_Yend = cal_year + 1,
        NPlst = NPlst, fp = fp_vec)
    }
    for (i in param_var) {
      b_pt       <- ((cal_year + 1) - pj$cabY) / pj$timestep + 1
      dim_length <- dim(dfList_NP_year[[i]])[3]
      dfList_NP_year[[i]][, , b_pt:dim_length] <- dfList_NP_base[[i]][, , b_pt:dim_length]
    }
    
    # Initialize fs using prev_fs for ALL pops
    fs_new <- prev_fs
    for (pop in 1:pj$npops) {
      fs_new[pop, ini_dt:end_dt] <- prev_fs[pop, (ini_dt - steps):(end_dt - steps)]
    }
    
    # Run model
    model_new <- HCVMSM(
      pj, best_estimates, best_est_pop, disease_progress, pop_array, dfList,
      param_cascade_sc = dfList_NP_year, fib = fib,
      modelrun = "UN", proj = "POC_AU", end_Y = endY,
      cost = NULL, costflow = NULL, costflow_Neg = NULL, fc = fs_new
    )
    
    list(model = model_new, dfList_NP = dfList_NP_year,
         fs = fs_new, Ccal = Ccal_yr, fp_vec = fp_vec,
         fs_prev_val = fs_prev_val, denom_C = denom_C, denom_P = denom_P)
  }
  
  # ── Joint optimization of adj_C and adj_P ───────────────────────────────────
  if (!is.null(target_C) & !is.null(target_P) &
      (is.null(adj_factor_C) | is.null(adj_factor_P))) {
    
    cat("=== Joint optimization of adj_C and adj_P for year", cal_year, "===\n")
    
    obj_fn_joint <- function(params) {
      adj_c <- params[1]
      adj_p <- params[2]
      
      if (adj_c <= 0 | adj_p <= 0) return(1e10)
      
      res  <- run_model(adj_P = adj_p, adj_C = adj_c)
      tots <- get_test_totals(res$model)
      
      tot_C <- tots %>% filter(setting == "C") %>% pull(tot)
      tot_P <- tots %>% filter(setting == "P") %>% pull(tot)
      
      err_C <- (tot_C - target_C)^2 / target_C^2
      err_P <- (tot_P - target_P)^2 / target_P^2
      
      cat("  adj_C:", round(adj_c, 4), " adj_P:", round(adj_p, 4),
          "→ C:", round(tot_C), "/", target_C,
          " P:", round(tot_P), "/", target_P, "\n")
      
      err_C + err_P
    }
    
    # Starting values: use previous year adj_factors if available
    start_C <- ifelse(is.null(adj_factor_C), 1.0,  adj_factor_C)
    start_P <- ifelse(is.null(adj_factor_P), 0.44, adj_factor_P)
    
    opt <- optim(
      par     = c(start_C, start_P),
      fn      = obj_fn_joint,
      method  = "Nelder-Mead",
      control = list(maxit = 500, reltol = 1e-5)
    )
    
    adj_factor_C <- opt$par[1]
    adj_factor_P <- opt$par[2]
    cat("=== Best adj_factor_C:", round(adj_factor_C, 4),
        " adj_factor_P:", round(adj_factor_P, 4), "===\n")
    
  } else {
    # Use provided values or defaults
    if (is.null(adj_factor_C)) adj_factor_C <- 1.0
    if (is.null(adj_factor_P)) adj_factor_P <- 0.44
  }
  
  # ── Final model run ──────────────────────────────────────────────────────────
  cat("=== Final run: adj_C:", round(adj_factor_C, 4),
      " adj_P:", round(adj_factor_P, 4), "===\n")
  final <- run_model(adj_P = adj_factor_P, adj_C = adj_factor_C)
  
  # ── Extract U and S from final model ────────────────────────────────────────
  pop_new <- extract_pop_US(final$model, at_year = cal_year)
  su_new  <- lapply(pj$popNames, function(p) get_SU(pop_new, p))
  names(su_new) <- pj$popNames
  
  undiag_C <- su_new[["C_PWID"]]$U  + su_new[["C_fPWID"]]$U
  undiag_P <- su_new[["P_PWID"]]$U  + su_new[["P_fPWID"]]$U + su_new[["P_nPWID"]]$U
  s_bar_C  <- su_new[["C_PWID"]]$S  + su_new[["C_fPWID"]]$S
  s_bar_P  <- su_new[["P_PWID"]]$S  + su_new[["P_fPWID"]]$S + su_new[["P_nPWID"]]$S
  
  if (!is.finite(s_bar_C) | s_bar_C == 0) stop(paste("s_bar_C is", s_bar_C))
  if (!is.finite(s_bar_P) | s_bar_P == 0) stop(paste("s_bar_P is", s_bar_P))
  
  # ── Recalculate fs (xfs/fm — matches original fs_estimate) ──────────────────
  fab   <- frac_ab[[yr_chr]]
  n_ab  <- n_ab_np[[yr_chr]]
  cov_C <- final$Ccal$C
  cov_P <- final$Ccal$P
  
  xfs_1 <- (as.numeric(n_ab[1])/(cov_C*fab[1]) - fm[[yr_chr]][1]*undiag_C) / s_bar_C
  xfs_2 <- (as.numeric(n_ab[1])/(cov_C*fab[1]) - fm[[yr_chr]][2]*undiag_C) / s_bar_C
  xfs_3 <- (as.numeric(n_ab[2])/(cov_P*fab[2]) - fm[[yr_chr]][3]*undiag_P) / s_bar_P
  xfs_4 <- (as.numeric(n_ab[2])/(cov_P*fab[2]) - fm[[yr_chr]][4]*undiag_P) / s_bar_P
  
  cat("xfs values:", xfs_1, xfs_2, xfs_3, xfs_4, "\n")
  
  fs_final <- final$fs
  fs_final[1, ini_dt:end_dt] <- rep(xfs_1 / fm[[yr_chr]][1], steps)
  fs_final[2, ini_dt:end_dt] <- rep(xfs_2 / fm[[yr_chr]][2], steps)
  fs_final[3, ini_dt:end_dt] <- rep(xfs_3 / fm[[yr_chr]][3], steps)
  fs_final[4, ini_dt:end_dt] <- rep(xfs_4 / fm[[yr_chr]][4], steps)
  fs_final[5, ini_dt:end_dt] <- rep(1, steps)
  
  cat("fs values:", fs_final[1,ini_dt], fs_final[2,ini_dt],
      fs_final[3,ini_dt], fs_final[4,ini_dt], "\n")
  
  # ── Final testing totals ─────────────────────────────────────────────────────
  final_tots <- get_test_totals(final$model)
  tot_C <- final_tots %>% filter(setting == "C") %>% pull(tot)
  tot_P <- final_tots %>% filter(setting == "P") %>% pull(tot)
  
  cat("=== Final testing numbers ===\n")
  cat("C:", round(tot_C), "/", ifelse(is.null(target_C), "no target", target_C),
      "(", round((tot_C - ifelse(is.null(target_C), tot_C, target_C))/
                   ifelse(is.null(target_C), tot_C, target_C)*100, 2), "%)\n")
  cat("P:", round(tot_P), "/", ifelse(is.null(target_P), "no target", target_P),
      "(", round((tot_P - ifelse(is.null(target_P), tot_P, target_P))/
                   ifelse(is.null(target_P), tot_P, target_P)*100, 2), "%)\n")
  
  return(list(
    dfList_NP    = final$dfList_NP,
    fs           = fs_final,
    model        = final$model,
    Ccal         = final$Ccal,
    adj_factor_C = adj_factor_C,
    adj_factor_P = adj_factor_P,
    diagnostics  = list(
      tot_C       = tot_C,       tot_P    = tot_P,
      target_C    = target_C,    target_P = target_P,
      undiag_C    = undiag_C,    undiag_P = undiag_P,
      s_bar_C     = s_bar_C,     s_bar_P  = s_bar_P,
      denom_C     = final$denom_C, denom_P = final$denom_P,
      fp_vec      = final$fp_vec,  Ccal    = final$Ccal,
      fs_prev_val = final$fs_prev_val
    )
  ))
}


# ── Define testing targets per year ──────────────────────────────────────────
target_tests <- list(
  "2025" = list(C = 11000, P = 14000),
  "2026" = list(C = 11000, P = 14000),
  "2027" = list(C = 11000, P = 14000),
  "2028" = list(C = 11000, P = 14000),
  "2029" = list(C = 11000, P = 14000),
  "2030" = list(C = 11000, P = 14000)
)

# ── Define n_ab_np for all years ─────────────────────────────────────────────
n_ab_np[["2025"]] <- c(11000, 14000)
n_ab_np[["2026"]] <- c(11000, 14000)
n_ab_np[["2027"]] <- c(11000, 14000)
n_ab_np[["2028"]] <- c(11000, 14000)
n_ab_np[["2029"]] <- c(11000, 14000)
n_ab_np[["2030"]] <- c(11000, 14000)

# ── Define fm for all years ───────────────────────────────────────────────────
fm[["2025"]] <- fm[["2024"]]
fm[["2026"]] <- fm[["2024"]]
fm[["2027"]] <- fm[["2024"]]
fm[["2028"]] <- fm[["2024"]]
fm[["2029"]] <- fm[["2024"]]
fm[["2030"]] <- fm[["2024"]]

# ── Define frac_test for all years ───────────────────────────────────────────
frac_test[[2025]] <- frac_test[[2024]]
frac_test[[2026]] <- frac_test[[2024]]
frac_test[[2027]] <- frac_test[[2024]]
frac_test[[2028]] <- frac_test[[2024]]
frac_test[[2029]] <- frac_test[[2024]]
frac_test[[2030]] <- frac_test[[2024]]

frac_ab[["2025"]] <- c(unlist(as.numeric(frac_test[[2025]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2025]]$P$reflex)))
frac_ab[["2026"]] <- c(unlist(as.numeric(frac_test[[2026]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2026]]$P$reflex)))
frac_ab[["2027"]] <- c(unlist(as.numeric(frac_test[[2027]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2027]]$P$reflex)))
frac_ab[["2028"]] <- c(unlist(as.numeric(frac_test[[2028]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2028]]$P$reflex)))
frac_ab[["2029"]] <- c(unlist(as.numeric(frac_test[[2029]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2029]]$P$reflex)))

frac_ab[["2030"]] <- c(unlist(as.numeric(frac_test[[2030]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2030]]$P$reflex)))

# ── Storage for results ───────────────────────────────────────────────────────
res_list     <- list()
adj_factors  <- list()  # store found adj_factors for inspection

# ── 2025 ─────────────────────────────────────────────────────────────────────
cat("\n========== YEAR 2025 ==========\n")
res_list[["2025"]] <- calibrate_NP_year(
  cal_year       = 2025,
  prev_dfList_NP = dfList_NP_2024,
  prev_fs        = fs[["2024"]],
  prev_model     = res_2024,          # your existing 2024 NP model
  dfList_NP_base = dfList_NP_2024,
  pj             = POC_AU,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  n_ab_np = n_ab_np, frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY,
  target_C     = target_tests[["2025"]]$C,
  target_P     = target_tests[["2025"]]$P,
  adj_factor_P = NULL)

# Extract and store
dfList_NP_2025      <- res_list[["2025"]]$dfList_NP
fs[["2025"]]        <- res_list[["2025"]]$fs
Ccal[[2025]]        <- res_list[["2025"]]$Ccal
adj_factors[["2025"]] <- res_list[["2025"]]$adj_factor_P
cat("2025 adj_factor_P:", adj_factors[["2025"]], "\n")
cat("2025 diagnostics:\n"); print(res_list[["2025"]]$diagnostics)

adj_factors_C[["2025"]] <- res_list[["2025"]]$adj_factor_C
adj_factors_P[["2025"]] <- res_list[["2025"]]$adj_factor_P
# ── 2026 ─────────────────────────────────────────────────────────────────────
cat("\n========== YEAR 2026 ==========\n")
res_list[["2026"]] <- calibrate_NP_year(
  cal_year       = 2026,
  prev_dfList_NP = dfList_NP_2025,
  prev_fs        = fs[["2025"]],
  prev_model     = res_list[["2025"]]$model,
  dfList_NP_base = dfList_NP_2025,
  pj             = POC_AU,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  n_ab_np = n_ab_np, frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY,
  target_C     = target_tests[["2026"]]$C,
  target_P     = target_tests[["2026"]]$P,
  adj_factor_C = NULL,  # warm start
  adj_factor_P = NULL)

dfList_NP_2026        <- res_list[["2026"]]$dfList_NP
fs[["2026"]]          <- res_list[["2026"]]$fs
Ccal[[2026]]          <- res_list[["2026"]]$Ccal
adj_factors[["2026"]] <- res_list[["2026"]]$adj_factor_P
cat("2026 adj_factor_P:", adj_factors[["2026"]], "\n")
cat("2026 diagnostics:\n"); print(res_list[["2026"]]$diagnostics)

# ── 2027 ─────────────────────────────────────────────────────────────────────
cat("\n========== YEAR 2027 ==========\n")
res_list[["2027"]] <- calibrate_NP_year(
  cal_year       = 2027,
  prev_dfList_NP = dfList_NP_2026,
  prev_fs        = fs[["2026"]],
  prev_model     = res_list[["2026"]]$model,
  dfList_NP_base = dfList_NP_2026,
  pj             = POC_AU,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  n_ab_np = n_ab_np, frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY,
  target_C     = target_tests[["2027"]]$C,
  target_P     = target_tests[["2027"]]$P,
  adj_factor_P = NULL
)

dfList_NP_2027        <- res_list[["2027"]]$dfList_NP
fs[["2027"]]          <- res_list[["2027"]]$fs
Ccal[[2027]]          <- res_list[["2027"]]$Ccal
adj_factors[["2027"]] <- res_list[["2027"]]$adj_factor_P
cat("2027 adj_factor_P:", adj_factors[["2027"]], "\n")
cat("2027 diagnostics:\n"); print(res_list[["2027"]]$diagnostics)

# ── 2028 ─────────────────────────────────────────────────────────────────────
cat("\n========== YEAR 2028 ==========\n")
res_list[["2028"]] <- calibrate_NP_year(
  cal_year       = 2028,
  prev_dfList_NP = dfList_NP_2027,
  prev_fs        = fs[["2027"]],
  prev_model     = res_list[["2027"]]$model,
  dfList_NP_base = dfList_NP_2027,
  pj             = POC_AU,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  n_ab_np = n_ab_np, frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY,
  target_C     = target_tests[["2028"]]$C,
  target_P     = target_tests[["2028"]]$P,
  adj_factor_P = NULL
)

dfList_NP_2028        <- res_list[["2028"]]$dfList_NP
fs[["2028"]]          <- res_list[["2028"]]$fs
Ccal[[2028]]          <- res_list[["2028"]]$Ccal
adj_factors[["2028"]] <- res_list[["2028"]]$adj_factor_P
cat("2028 adj_factor_P:", adj_factors[["2028"]], "\n")
cat("2028 diagnostics:\n"); print(res_list[["2028"]]$diagnostics)

# ── 2029 ─────────────────────────────────────────────────────────────────────
cat("\n========== YEAR 2029 ==========\n")
res_list[["2029"]] <- calibrate_NP_year(
  cal_year       = 2029,
  prev_dfList_NP = dfList_NP_2028,
  prev_fs        = fs[["2028"]],
  prev_model     = res_list[["2028"]]$model,
  dfList_NP_base = dfList_NP_2028,
  pj             = POC_AU,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  n_ab_np = n_ab_np, frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY,
  target_C     = target_tests[["2029"]]$C,
  target_P     = target_tests[["2029"]]$P,
  adj_factor_P = NULL
)

dfList_NP_2029        <- res_list[["2029"]]$dfList_NP
fs[["2029"]]          <- res_list[["2029"]]$fs
Ccal[[2029]]          <- res_list[["2029"]]$Ccal
adj_factors[["2029"]] <- res_list[["2029"]]$adj_factor_P
cat("2029 adj_factor_P:", adj_factors[["2029"]], "\n")
cat("2029 diagnostics:\n"); print(res_list[["2029"]]$diagnostics)

# ── 2030 ─────────────────────────────────────────────────────────────────────
cat("\n========== YEAR 2030 ==========\n")
res_list[["2030"]] <- calibrate_NP_year(
  cal_year       = 2030,
  prev_dfList_NP = dfList_NP_2029,
  prev_fs        = fs[["2029"]],
  prev_model     = res_list[["2029"]]$model,
  dfList_NP_base = dfList_NP_2029,
  pj             = POC_AU,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  n_ab_np = n_ab_np, frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY,
  target_C     = target_tests[["2030"]]$C,
  target_P     = target_tests[["2030"]]$P,
  adj_factor_P = NULL
)

dfList_NP_2030        <- res_list[["2030"]]$dfList_NP
fs[["2030"]]          <- res_list[["2030"]]$fs
Ccal[[2030]]          <- res_list[["2030"]]$Ccal
adj_factors[["2030"]] <- res_list[["2030"]]$adj_factor_P
cat("2030 adj_factor_P:", adj_factors[["2030"]], "\n")
cat("2030 diagnostics:\n"); print(res_list[["2030"]]$diagnostics)

# ── Summary of all adj_factors and testing numbers ───────────────────────────
cat("\n========== SUMMARY ==========\n")
summary_df <- data.frame(
  year       = 2025:2026,
  adj_factor = unlist(adj_factors),
  tot_C      = sapply(res_list, function(r) r$diagnostics$tot_C),
  target_C   = sapply(target_tests, function(t) t$C),
  tot_P      = sapply(res_list, function(r) r$diagnostics$tot_P),
  target_P   = sapply(target_tests, function(t) t$P),
  pct_diff_C = sapply(res_list, function(r) 
    round((r$diagnostics$tot_C - r$diagnostics$target_C)/r$diagnostics$target_C*100, 2)),
  pct_diff_P = sapply(res_list, function(r) 
    round((r$diagnostics$tot_P - r$diagnostics$target_P)/r$diagnostics$target_P*100, 2))
)
print(summary_df)
