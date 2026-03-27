endY <- 100
res_2024 <- HCVMSM_cpp(POC_AU, best_estimates, best_est_pop,
       disease_progress, pop_array,
       dfList,  
       param_cascade_sc = dfList_NP_2024, 
       fib = fib, 
       modelrun = "UN", proj = "POC_AU", end_Y = endY, 
       
       fc = fs[["2024"]])


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
    result <- popResults_MidYear(
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
      summarise(best = sum(best), .groups = "drop")
    
    # ── Diagnostic: check what diag_status values exist ─────────────────────
    cat("diag_status values at year", at_year, ":", 
        unique(result$diag_status), "\n")
    cat("n rows:", nrow(result), "\n")
    
    result <- result %>%
      pivot_wider(names_from = diag_status, values_from = best, values_fill = 0)
    
    # ── Add missing S/U columns if pivot_wider didn't create them ───────────
    if (!"S" %in% names(result)) {
      cat("WARNING: S column missing — adding zeros\n")
      result$S <- 0
    }
    if (!"U" %in% names(result)) {
      cat("WARNING: U column missing — adding zeros\n")
      result$U <- 0
    }
    
    result %>%
      mutate(population = as.character(population),
             S = as.numeric(S),
             U = as.numeric(U))
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
    
    fs_prev_val <- as.numeric(prev_fs[, ini_dt - 1])
    if (any(!is.finite(fs_prev_val)) | any(fs_prev_val == 0))
      fs_prev_val <- as.numeric(prev_fs[, ini_dt])
    
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
    
    # Build dfList_NP for current year only
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
    model_new <- HCVMSM_cpp(
      pj, best_estimates, best_est_pop, disease_progress, pop_array, dfList,
      param_cascade_sc = dfList_NP_year, fib = fib,
      modelrun = "UN", proj = "POC_AU", end_Y = endY,
      fc_sc = fs_new
    )
    
    list(model     = model_new,
         dfList_NP = dfList_NP_year,
         fs        = fs_new,        # ← pre-xfs fc actually passed to HCVMSM
         Ccal      = Ccal_yr,
         fp_vec    = fp_vec,
         fs_prev_val = fs_prev_val,
         denom_C   = denom_C,
         denom_P   = denom_P)
  }
  
  # ── Joint optimization ───────────────────────────────────────────────────────
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
  
  # ── Recalculate fs ───────────────────────────────────────────────────────────
  fab   <- frac_ab[[yr_chr]]
  n_ab  <- n_ab_np[[yr_chr]]
  cov_C <- final$Ccal$C
  cov_P <- final$Ccal$P
  
  xfs_1 <- (as.numeric(n_ab[1])/(cov_C*fab[1]) - fm[[yr_chr]][1]*undiag_C) / s_bar_C
  xfs_2 <- (as.numeric(n_ab[1])/(cov_C*fab[1]) - fm[[yr_chr]][2]*undiag_C) / s_bar_C
  xfs_3 <- (as.numeric(n_ab[2])/(cov_P*fab[2]) - fm[[yr_chr]][3]*undiag_P) / s_bar_P
  xfs_4 <- (as.numeric(n_ab[2])/(cov_P*fab[2]) - fm[[yr_chr]][4]*undiag_P) / s_bar_P
  
  cat("xfs values:", xfs_1, xfs_2, xfs_3, xfs_4, "\n")
  
  # fs_final = post-xfs, used as prev_fs for NEXT year's chain
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
    fs           = fs_final,     # post-xfs → feed into next year as prev_fs
    fc_used      = final$fs,     # ← pre-xfs fc actually passed to HCVMSM
    model        = final$model,
    Ccal         = final$Ccal,
    adj_factor_C = adj_factor_C,
    adj_factor_P = adj_factor_P,
    diagnostics  = list(
      tot_C       = tot_C,         tot_P    = tot_P,
      target_C    = target_C,      target_P = target_P,
      undiag_C    = undiag_C,      undiag_P = undiag_P,
      s_bar_C     = s_bar_C,       s_bar_P  = s_bar_P,
      denom_C     = final$denom_C, denom_P  = final$denom_P,
      fp_vec      = final$fp_vec,  Ccal     = final$Ccal,
      fs_prev_val = final$fs_prev_val
    )
  ))
}

# ── function to simulate the model and optimise adjusting factors in sequnce ──
run_NP_scenario <- function(
    scenario_name,
    target_tests,        # named list: list("2025"=list(C=,P=), ..., "2030"=list(C=,P=))
    pj, fs_base,         # fs[["2024"]]
    prev_model_base,     # res_2024
    dfList_NP_base,      # dfList_NP_2024
    Ccal, fm, frac_test, NPlst, frac_ab, param_var,
    best_estimates, best_est_pop, disease_progress,
    pop_array, dfList, fib, endY
) {
  cat("\n══════════════════════════════════════════\n")
  cat("SCENARIO:", scenario_name, "\n")
  cat("══════════════════════════════════════════\n")
  
  years <- as.character(2025:2030)
  
  # ── Build n_ab_np from targets ──────────────────────────────────────────────
  n_ab_np_scen <- list()
  for (yr in years) {
    n_ab_np_scen[[yr]] <- c(target_tests[[yr]]$C, target_tests[[yr]]$P)
  }
  
  # ── Storage ──────────────────────────────────────────────────────────────────
  res_list    <- list()
  adj_factors <- list()
  fs_scen     <- fs_base  # start fresh from 2024 base
  prev_dfList <- dfList_NP_base
  prev_fs     <- fs_base[["2024"]]
  prev_model  <- prev_model_base
  
  # ── Run each year sequentially ───────────────────────────────────────────────
  for (yr in years) {
    cat("\n========== SCENARIO:", scenario_name, "— YEAR", yr, "==========\n")
    
    res <- calibrate_NP_year(
      cal_year       = as.numeric(yr),
      prev_dfList_NP = prev_dfList,
      prev_fs        = prev_fs,
      prev_model     = prev_model,
      dfList_NP_base = prev_dfList,
      pj             = pj,
      Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
      n_ab_np  = n_ab_np_scen, frac_ab = frac_ab, param_var = param_var,
      best_estimates = best_estimates, best_est_pop = best_est_pop,
      disease_progress = disease_progress, pop_array = pop_array,
      dfList = dfList, fib = fib, endY = endY,
      target_C     = target_tests[[yr]]$C,
      target_P     = target_tests[[yr]]$P,
      adj_factor_C = NULL,
      adj_factor_P = NULL
    )
    
    res_list[[yr]]    <- res
    adj_factors[[yr]] <- list(C = res$adj_factor_C, P = res$adj_factor_P)
    fs_scen[[yr]]     <- res$fs
    
    cat("adj_C:", round(res$adj_factor_C, 4),
        " adj_P:", round(res$adj_factor_P, 4), "\n")
    
    # ── Chain to next year ──────────────────────────────────────────────────
    prev_dfList <- res$dfList_NP
    prev_fs     <- res$fs
    prev_model  <- res$model
  }
  
  # ── Build combined fc and dfList_NP ─────────────────────────────────────────
  cat("\n========== BUILDING COMBINED MODEL:", scenario_name, "==========\n")
  
  fc_combined       <- fs_base[["2024"]]
  dfList_NP_combined <- dfList_NP_base
  
  for (yr in years) {
    yr_num <- as.numeric(yr)
    ini_dt <- (yr_num - pj$cabY) / pj$timestep + 1
    end_dt <- ((yr_num + 1) - pj$cabY) / pj$timestep
    
    fc_combined[, ini_dt:end_dt] <- res_list[[yr]]$fc_used[, ini_dt:end_dt]
    
    for (i in param_var) {
      dfList_NP_combined[[i]][, , ini_dt:end_dt] <-
        res_list[[yr]]$dfList_NP[[i]][, , ini_dt:end_dt]
    }
  }
  # ── check np testing  ───────────────────────────────────────────────────────
  get_NP_testing_check <- function(model, label) {
    cl_ext <- c("newTestingAb_sc","newTestingAg_sc","newTestingPOCT_sc",
                "newTestingAb_sc_neg","newTestingAg_sc_neg","newTestingPOCT_sc_neg")
    
    tflow <- list()
    for (i in cl_ext) {
      tflow[[i]] <- modres.flow.t(POC_AU, model, endYear = 100, allp = i) %>%
        ungroup() %>%
        group_by(year, population) %>%
        summarise(best = sum(best), .groups = "drop")
    }
    
    dplyr::bind_rows(tflow, .id = "index") %>%
      group_by(year, population) %>%
      spread(index, best) %>%
      filter(year %in% c(10:15)) %>%   # 2025-2030
      mutate(
        Ab  = rowSums(cbind(newTestingAb_sc,  newTestingAb_sc_neg),  na.rm = TRUE),
        RNA = rowSums(cbind(newTestingAg_sc,  newTestingAg_sc_neg,
                            newTestingPOCT_sc, newTestingPOCT_sc_neg), na.rm = TRUE),
        setting = ifelse(population %in% c("C_PWID","C_fPWID"), "C", "P"),
        cal_yr  = year + POC_AU$cabY
      ) %>%
      group_by(cal_yr, setting) %>%
      summarise(tot = sum(Ab + RNA, na.rm = TRUE), .groups = "drop") %>%
      mutate(source = label)
  }
  # ── Run combined model ───────────────────────────────────────────────────────
  model_combined <- HCVMSM_cpp(
    pj, best_estimates, best_est_pop, disease_progress, pop_array, dfList,
    param_cascade_sc = dfList_NP_combined, fib = fib,
    modelrun = "UN", proj = "POC_AU", end_Y = endY,
    fc_sc = fc_combined
  )
  
  # ── Verify ───────────────────────────────────────────────────────────────────
  targets_df <- data.frame(
    cal_yr  = rep(2025:2030, each = 2),
    setting = rep(c("C","P"), 6),
    target  = unlist(lapply(years, function(yr)
      c(target_tests[[yr]]$C, target_tests[[yr]]$P)))
  )
  
  cat("\n=== Verification:", scenario_name, "===\n")
  chk <- get_NP_testing_check(model_combined, scenario_name) %>%
    left_join(targets_df, by = c("cal_yr","setting")) %>%
    mutate(pct_diff = round((tot - target)/target*100, 1))
  print(chk)
  
  # ── Return ───────────────────────────────────────────────────────────────────
  return(list(
    model          = model_combined,
    res_list       = res_list,
    adj_factors    = adj_factors,
    fs             = fs_scen,
    dfList_NP      = dfList_NP_combined,
    fc_combined    = fc_combined,
    verification   = chk
  ))
} 




# ── Define fm for all years ───────────────────────────────────────────────────
fm <- list()
fm[["2024"]] <- c(1.1, 1.1, 18.5, 18.5, 1)

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


# ── Scenario 1: sustained ─────────────────────────────────────────────────────
scen1 <- run_NP_scenario(
  scenario_name  = "sustained",
  target_tests   = list(
    "2025" = list(C = 11000, P = 14000),
    "2026" = list(C = 11000, P = 14000),
    "2027" = list(C = 11000, P = 14000),
    "2028" = list(C = 11000, P = 14000),
    "2029" = list(C = 11000, P = 14000),
    "2030" = list(C = 11000, P = 14000)
  ),
  pj = POC_AU, fs_base = fs, prev_model_base = res_2024,
  dfList_NP_base = dfList_NP_2024,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY
)

# ── Scenario 2: prison shift ──────────────────────────────────────────────────
scen2 <- run_NP_scenario(
  scenario_name  = "prison_shift",
  target_tests   = list(
    "2025" = list(C = 11000, P = 14000),
    "2026" = list(C = 10000, P = 15000),
    "2027" = list(C =  9000, P = 16000),
    "2028" = list(C =  8000, P = 17000),
    "2029" = list(C =  7000, P = 18000),
    "2030" = list(C =  6000, P = 19000)
  ),
  pj = POC_AU, fs_base = fs, prev_model_base = res_2024,
  dfList_NP_base = dfList_NP_2024,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY
)

# ── Scenario 3: prison shift ──────────────────────────────────────────────────
scen3 <- run_NP_scenario(
  scenario_name  = "scale_up",
  target_tests   = list(
    "2025" = list(C = 11000, P = 14000),
    "2026" = list(C = 11000, P = 14000),
    "2027" = list(C = 13750, P = 17500),
    "2028" = list(C = 16500, P = 21000),
    "2029" = list(C = 19250, P = 24500),
    "2030" = list(C = 22000, P = 28000)
  ),
  pj = POC_AU, fs_base = fs, prev_model_base = res_2024,
  dfList_NP_base = dfList_NP_2024,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY
)

# ── Scenario 4: prison shift scale-up ──────────────────────────────────────────────────
scen4 <- run_NP_scenario(
  scenario_name  = "prison_shift_scale_up",
  target_tests   = list(
    "2025" = list(C = 11000, P = 14000),
    "2026" = list(C = 11000, P = 14000),
    "2027" = list(C = 11000, P = 20250),
    "2028" = list(C = 11000, P = 26500),
    "2029" = list(C = 11000, P = 32750),
    "2030" = list(C = 11000, P = 39000)
  ),
  pj = POC_AU, fs_base = fs, prev_model_base = res_2024,
  dfList_NP_base = dfList_NP_2024,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY
)


# ── Scenario 5: community shift scale-up ──────────────────────────────────────────────────
scen5 <- run_NP_scenario(
  scenario_name  = "community_shift_scale_up",
  target_tests   = list(
    "2025" = list(C = 11000, P = 14000),
    "2026" = list(C = 11000, P = 14000),
    "2027" = list(C = 17250, P = 14000),
    "2028" = list(C = 23500, P = 14000),
    "2029" = list(C = 29750, P = 14000),
    "2030" = list(C = 36000, P = 14000)
  ),
  pj = POC_AU, fs_base = fs, prev_model_base = res_2024,
  dfList_NP_base = dfList_NP_2024,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY
)
# ── Add to all_scenarios ──────────────────────────────────────────────────────
all_scenarios <- list(
  "Status quo"   = Sce_sq,
  "sustained"    = scen1,
  "prison_shift" = scen2,
  "scale_up"     = scen3,
  "prison_shift_scale_up" = scen4,
  "community_shift_scale_up" = scen5
)

save(all_scenarios,
     file = file.path(OutputFolder,
                      paste0("testing_simulations",".rda")))

# get_NP_testing_check(model_NP_combined, "combined") %>%
#  left_join(targets, by = c("cal_yr","setting")) %>%
#  mutate(pct_diff = round((tot - target)/target*100, 1)) %>%
#  print()
