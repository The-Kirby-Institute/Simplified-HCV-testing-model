rm(list = ls())
library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
library(readxl)
library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")
project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Projects/Simplified-HCV-testing-model")

data_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                   project_name, sep = "")
# Load useful libraries

Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
OutputFolder <- file.path(data_path, "02. Output")

load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "cali.rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_timev.rda")))
load(file.path(OutputFolder, paste0(project_name, "_NPcal_2024.rda")))
source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotManuscript.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R")) 
source(file.path(codefun_path, "/Projects/", project_name, "/Param_np_fc_cal.R")) 


# =============================================================================
# HELPERS from calibration_pipeline.R
# =============================================================================

endY <- 100   # must match Script 1

diag_RNA_states <- c(
  "f0_diag_RNA","f1_diag_RNA","f2_diag_RNA",
  "f3_diag_RNA","f4_diag_RNA","dc_diag_RNA","hcc_diag_RNA",
  "lt_diag_RNA","plt_diag_RNA"
)
undiag_states  <- gsub("diag_RNA", "undiag", diag_RNA_states)
neg_states     <- c("s","f0_cured","f1_cured","f2_cured","f3_cured",
                    "f4_cured","dc_cured","hcc_cured","lt_cured","plt_cured")
diag_ab_states <- c("f0_diag_ab","f1_diag_ab","f2_diag_ab",
                    "f3_diag_ab","f4_diag_ab","dc_diag_ab","hcc_diag_ab",
                    "lt_diag_ab","plt_diag_ab")
sab_states     <- c("f0_cured","f1_cured","f2_cured","f3_cured","f4_cured",
                    "dc_cured","hcc_cured","lt_cured","plt_cured")
s_states       <- c("s")

year_time_index <- function(pj, year) {
  steps_per_year <- 1 / pj$timestep
  k  <- year - pj$cabY + 1
  t0 <- (k - 1) * steps_per_year + 1
  t1 <-  k      * steps_per_year
  list(t0 = t0, t1 = t1, steps_per_year = steps_per_year)
}

monthly_setting_stock <- function(Sce, pj, pops, states, year) {
  tt <- year_time_index(pj, year)
  sapply(tt$t0:tt$t1, function(t) {
    sum(vapply(
      states,
      function(st) sum(Sce$allPops[pops, st, t], na.rm = TRUE),
      numeric(1)
    ))
  })
}

safe_divide <- function(num, den, eps = 1e-12) num / pmax(den, eps)


# =============================================================================
# SCENARIO DEFINITIONS
# =============================================================================
# Each scenario specifies total tests (T_target) per year per setting for
# 2025-2030. 2022-2024 are always taken from the calibration.
#
# T_target = 0 means "no NP this year" (Cov_np = Cov_neg = 0).
# =============================================================================

scenario_years <- 2025:2030

# Helper: build a scenario from per-year totals + split fraction
make_test_scenario <- function(totals_by_year, comm_frac = 11000 / 25000) {
  tibble(
    year      = scenario_years,
    T_total   = totals_by_year,
    T_target_C = totals_by_year * comm_frac,
    T_target_P = totals_by_year * (1 - comm_frac)
  )
}

scenarios <- list(
  # Both = zero NP from 2025 onward
  no_np        = make_test_scenario(rep(0, 6)),
  foundational = make_test_scenario(rep(0, 6)),
  
  # Succession: 25k in 2025-2026, then nothing
  succession   = make_test_scenario(c(25000, 25000, 0, 0, 0, 0)),
  
  # Sustained: 25k every year
  sustained    = make_test_scenario(rep(25000, 6)),
  
  # Scale-up: 25k in 2025-2026, then linear ramp to 50k by 2030
  scaleup      = make_test_scenario(c(25000, 25000, 31250, 37500, 43750, 50000))
)
#
#
# =============================================================================
# AUTO-FIT WRAPPER: Find Cov_np, Cov_neg to hit T_target within 5%
# Per year, per setting. Chained across years 2025-2030.
# =============================================================================
# Strategy:
#   Phase 1: Single multiplier m applied to both Cov_np and Cov_neg (fc fixed).
#            Bisection in [m_min, m_max] = [1/3, 3].
#   Phase 2: If Phase 1 can't hit 5%, allow fc to drift by searching Cov_neg
#            independently while Cov_np stays at the Phase 1 best m.
#   Max 20 simulations per year per setting (combined).
#   Each simulation call uses run_single_year_manual under the hood.
# =============================================================================

fit_year_to_target <- function(year, prev_values,
                               T_target_C, T_target_P,
                               tolerance = 0.03,
                               max_iter_total = 20,
                               m_min = 1/3, m_max = 3,
                               verbose = TRUE) {
  
  # prev_values = list(Cov_np_C, Cov_neg_C, Cov_np_P_monthly, Cov_neg_P_monthly)
  
  # Helper: run one trial, return c(ratio_C, ratio_P, results_C, results_P)
  trial <- function(Cnp_C, Cneg_C, Cnp_Pm, Cneg_Pm, label = "") {
    if (verbose) {
      cat(sprintf("    [trial%s] Cov_np_C=%.5f Cov_neg_C=%.5f | Cov_np_P_m=%.5f Cov_neg_P_m=%.5f\n",
                  ifelse(label == "", "", paste0(" ", label)),
                  Cnp_C, Cneg_C, Cnp_Pm, Cneg_Pm))
    }
    res <- run_single_year_manual(
      year = year,
      Cov_np_C          = Cnp_C,
      Cov_neg_C         = Cneg_C,
      Cov_np_P_monthly  = Cnp_Pm,
      Cov_neg_P_monthly = Cneg_Pm,
      verbose = FALSE
    )
    T_C <- res$results_C$Total_tests
    T_P <- res$results_P$Total_tests
    ratio_C <- if (T_target_C > 0) T_C / T_target_C else NA_real_
    ratio_P <- if (T_target_P > 0) T_P / T_target_P else NA_real_
    if (verbose) {
      cat(sprintf("       → T_C=%.0f (ratio=%.3f), T_P=%.0f (ratio=%.3f)\n",
                  T_C, ratio_C %||% NA, T_P, ratio_P %||% NA))
    }
    list(Sce = res$Sce, T_C = T_C, T_P = T_P,
         ratio_C = ratio_C, ratio_P = ratio_P,
         results_C = res$results_C, results_P = res$results_P,
         Cnp_C = Cnp_C, Cneg_C = Cneg_C, Cnp_Pm = Cnp_Pm, Cneg_Pm = Cneg_Pm)
  }
  
  `%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
  
  iter_count <- 0
  trial_log  <- list()
  
  # Handle zero-target settings separately
  do_C <- T_target_C > 0
  do_P <- T_target_P > 0
  
  # Starting values: Phase 1 starts with multiplier = 1 (= prev values)
  m_C <- 1.0
  m_P <- 1.0
  m_C_lo <- m_min; m_C_hi <- m_max
  m_P_lo <- m_min; m_P_hi <- m_max
  
  # -------------------------------------------------------------------------
  # PHASE 1: SINGLE-MULTIPLIER BISECTION (alternating C and P)
  # -------------------------------------------------------------------------
  if (verbose) cat(sprintf("\n  PHASE 1: single-multiplier search for %d\n", year))
  
  best_C <- list(ratio = Inf, m = 1.0); best_P <- list(ratio = Inf, m = 1.0)
  converged_C <- !do_C
  converged_P <- !do_P
  
  phase1_max <- 12  # leave headroom for phase 2
  
  while (iter_count < phase1_max && (!converged_C || !converged_P)) {
    iter_count <- iter_count + 1
    
    Cnp_C  <- prev_values$Cov_np_C          * m_C
    Cneg_C <- prev_values$Cov_neg_C         * m_C
    Cnp_P  <- prev_values$Cov_np_P_monthly  * m_P
    Cneg_P <- prev_values$Cov_neg_P_monthly * m_P
    
    # Clamp to legal probability range
    Cnp_C  <- pmin(pmax(Cnp_C,  0), 0.999999)
    Cneg_C <- pmin(pmax(Cneg_C, 0), 0.999999)
    Cnp_P  <- pmin(pmax(Cnp_P,  0), 0.999999)
    Cneg_P <- pmin(pmax(Cneg_P, 0), 0.999999)
    
    t <- trial(Cnp_C, Cneg_C, Cnp_P, Cneg_P,
               label = sprintf("p1.%d m_C=%.3f m_P=%.3f", iter_count, m_C, m_P))
    trial_log[[length(trial_log) + 1]] <- t
    
    # Track best
    if (do_C && !is.na(t$ratio_C) && abs(t$ratio_C - 1) < abs(best_C$ratio - 1)) {
      best_C <- list(ratio = t$ratio_C, m = m_C, trial = t)
    }
    if (do_P && !is.na(t$ratio_P) && abs(t$ratio_P - 1) < abs(best_P$ratio - 1)) {
      best_P <- list(ratio = t$ratio_P, m = m_P, trial = t)
    }
    
    # Convergence
    if (do_C && !is.na(t$ratio_C) && abs(t$ratio_C - 1) <= tolerance) converged_C <- TRUE
    if (do_P && !is.na(t$ratio_P) && abs(t$ratio_P - 1) <= tolerance) converged_P <- TRUE
    
    if (converged_C && converged_P) break
    
    # Bisection update
    if (do_C && !converged_C) {
      if (t$ratio_C > 1) { m_C_hi <- m_C } else { m_C_lo <- m_C }
      m_C <- (m_C_lo + m_C_hi) / 2
    }
    if (do_P && !converged_P) {
      if (t$ratio_P > 1) { m_P_hi <- m_P } else { m_P_lo <- m_P }
      m_P <- (m_P_lo + m_P_hi) / 2
    }
  }
  
  if (verbose) {
    cat(sprintf("\n  PHASE 1 results: m_C=%.4f (ratio=%.3f, conv=%s)  m_P=%.4f (ratio=%.3f, conv=%s)\n",
                best_C$m, best_C$ratio, converged_C,
                best_P$m, best_P$ratio, converged_P))
  }
  
  # -------------------------------------------------------------------------
  # PHASE 2: INDEPENDENT Cov_neg SEARCH (if Phase 1 failed)
  # -------------------------------------------------------------------------
  # Hold Cov_np at Phase 1 best, search Cov_neg independently
  # Favour neg in same direction as Cov_np (already moved by m)
  if ((!converged_C || !converged_P) && iter_count < max_iter_total) {
    if (verbose) cat(sprintf("\n  PHASE 2: independent Cov_neg search for %d\n", year))
    
    # Setup Phase 2 starting state from best Phase 1 trial
    Cnp_C_locked <- prev_values$Cov_np_C  * best_C$m
    Cnp_P_locked <- prev_values$Cov_np_P_monthly * best_P$m
    Cneg_C_cur   <- prev_values$Cov_neg_C * best_C$m
    Cneg_P_cur   <- prev_values$Cov_neg_P_monthly * best_P$m
    
    # Bisection bounds for Cov_neg (allow ±3x from current)
    Cneg_C_lo <- Cneg_C_cur / 3; Cneg_C_hi <- Cneg_C_cur * 3
    Cneg_P_lo <- Cneg_P_cur / 3; Cneg_P_hi <- Cneg_P_cur * 3
    
    while (iter_count < max_iter_total && (!converged_C || !converged_P)) {
      iter_count <- iter_count + 1
      
      Cnp_C  <- pmin(pmax(Cnp_C_locked, 0), 0.999999)
      Cneg_C <- pmin(pmax(Cneg_C_cur,    0), 0.999999)
      Cnp_P  <- pmin(pmax(Cnp_P_locked, 0), 0.999999)
      Cneg_P <- pmin(pmax(Cneg_P_cur,    0), 0.999999)
      
      t <- trial(Cnp_C, Cneg_C, Cnp_P, Cneg_P,
                 label = sprintf("p2.%d Cneg_C=%.5f Cneg_P=%.5f", iter_count, Cneg_C, Cneg_P))
      trial_log[[length(trial_log) + 1]] <- t
      
      if (do_C && !is.na(t$ratio_C) && abs(t$ratio_C - 1) < abs(best_C$ratio - 1)) {
        best_C <- list(ratio = t$ratio_C, m = best_C$m, trial = t,
                       Cneg = Cneg_C, phase = 2)
      }
      if (do_P && !is.na(t$ratio_P) && abs(t$ratio_P - 1) < abs(best_P$ratio - 1)) {
        best_P <- list(ratio = t$ratio_P, m = best_P$m, trial = t,
                       Cneg = Cneg_P, phase = 2)
      }
      
      if (do_C && !is.na(t$ratio_C) && abs(t$ratio_C - 1) <= tolerance) converged_C <- TRUE
      if (do_P && !is.na(t$ratio_P) && abs(t$ratio_P - 1) <= tolerance) converged_P <- TRUE
      
      if (converged_C && converged_P) break
      
      # Bisect Cov_neg
      if (do_C && !converged_C) {
        if (t$ratio_C > 1) { Cneg_C_hi <- Cneg_C_cur } else { Cneg_C_lo <- Cneg_C_cur }
        Cneg_C_cur <- (Cneg_C_lo + Cneg_C_hi) / 2
      }
      if (do_P && !converged_P) {
        if (t$ratio_P > 1) { Cneg_P_hi <- Cneg_P_cur } else { Cneg_P_lo <- Cneg_P_cur }
        Cneg_P_cur <- (Cneg_P_lo + Cneg_P_hi) / 2
      }
    }
  }
  
  # -------------------------------------------------------------------------
  # Return best trial (final solved values to chain into next year)
  # -------------------------------------------------------------------------
  # Pick the best COMBINED trial: closest to ratio=1 on both settings
  best_combined <- trial_log[[ which.min(sapply(trial_log, function(t) {
    rc <- if (!do_C) 0 else abs(t$ratio_C - 1)
    rp <- if (!do_P) 0 else abs(t$ratio_P - 1)
    rc + rp
  })) ]]
  
  if (verbose) {
    cat(sprintf("\n  FINAL year %d: iter=%d, ratio_C=%.3f, ratio_P=%.3f, converged_C=%s converged_P=%s\n",
                year, iter_count,
                best_combined$ratio_C %||% NA,
                best_combined$ratio_P %||% NA,
                converged_C, converged_P))
  }
  
  list(
    year             = year,
    converged_C      = converged_C,
    converged_P      = converged_P,
    iterations       = iter_count,
    Cov_np_C         = best_combined$Cnp_C,
    Cov_neg_C        = best_combined$Cneg_C,
    Cov_np_P_monthly = best_combined$Cnp_Pm,
    Cov_neg_P_monthly = best_combined$Cneg_Pm,
    ratio_C          = best_combined$ratio_C,
    ratio_P          = best_combined$ratio_P,
    T_C              = best_combined$T_C,
    T_P              = best_combined$T_P,
    Sce              = best_combined$Sce,
    results_C        = best_combined$results_C,
    results_P        = best_combined$results_P,
    trial_log        = trial_log
  )
}

run_single_year_manual <- function(year = 2025,
                                   Cov_np_C, Cov_neg_C,
                                   Cov_np_P_monthly, Cov_neg_P_monthly,
                                   verbose = TRUE) {
  
  zero_cascade <- lapply(dfList, function(x) x * 0)
  steps_per_year <- 1 / POC_AU$timestep
  m2a <- function(x) {
    x <- pmin(pmax(x, 0), 1)
    1 - (1 - x) ^ steps_per_year
  }
  
  # ---- Build starting cascade (2024 calibrated, all future zeroed) ----
  dfList_NP <- dfList_NP_2024
  dfList_CT <- dfList_CT_NP_2024
  fc_full_local <- fc_full
  
  tt_year <- year_time_index(POC_AU, year)
  zero_range <- tt_year$t0 : dim(dfList$eta)[3]
  for (i in param_var) {
    dfList_NP[[i]][, , zero_range] <- 0
    dfList_CT[[i]][, , zero_range] <- dfList[[i]][, , zero_range]
  }
  fc_full_local[, zero_range] <- 0
  
  # ---- Convert monthly prison values to annual-equivalent ----
  Cov_np_P_annual  <- m2a(Cov_np_P_monthly)
  Cov_neg_P_annual <- m2a(Cov_neg_P_monthly)
  
  # ---- Build fp_yr following calibration convention ----
  # C_PWID, C_fPWID get community Cov_np (annual)
  # P_PWID, P_fPWID get prison Cov_np (annual-equiv of monthly)
  # P_nPWID gets Cov_x = Cov_neg (annual-equiv of monthly Cov_neg)
  fp_yr <- c(
    Cov_np_C,           # C_PWID
    Cov_np_C,           # C_fPWID
    Cov_np_P_annual,    # P_PWID
    Cov_np_P_annual,    # P_fPWID
    Cov_neg_P_annual    # P_nPWID  <-- this is Cov_x, NOT positive Cov_np
  )
  fp_yr <- pmin(pmax(fp_yr, 0), 0.999999)
  
  # ---- fc per pop following calibration convention ----
  fc_C         <- if (Cov_np_C > 1e-12) Cov_neg_C / Cov_np_C else 0
  fc_P_monthly <- if (Cov_np_P_monthly > 1e-12) Cov_neg_P_monthly / Cov_np_P_monthly else 0
  fc_per_pop <- c(
    fc_C,           # C_PWID
    fc_C,           # C_fPWID
    fc_P_monthly,   # P_PWID
    fc_P_monthly,   # P_fPWID
    1               # P_nPWID  <-- fc=1 because its Cov_np already IS Cov_neg
  )
  
  if (verbose) {
    cat("=== MANUAL INPUTS (year", year, ") ===\n")
    cat(sprintf("  Community:  Cov_np  = %.6f   Cov_neg = %.6f   fc = %.4f\n",
                Cov_np_C, Cov_neg_C, fc_C))
    cat(sprintf("  Prison:     Cov_np_monthly  = %.6f   Cov_neg_monthly = %.6f   fc = %.4f\n",
                Cov_np_P_monthly, Cov_neg_P_monthly, fc_P_monthly))
    cat(sprintf("  fp_yr passed to Param_cal (annual-equivalent for prison):\n"))
    cat(sprintf("    C_PWID=%.4f  C_fPWID=%.4f  P_PWID=%.4f  P_fPWID=%.4f  P_nPWID=%.4f\n",
                fp_yr[1], fp_yr[2], fp_yr[3], fp_yr[4], fp_yr[5]))
    cat(sprintf("  fc per pop applied via fc_full:\n"))
    cat(sprintf("    C_PWID=%.4f  C_fPWID=%.4f  P_PWID=%.4f  P_fPWID=%.4f  P_nPWID=%.4f\n",
                fc_per_pop[1], fc_per_pop[2], fc_per_pop[3], fc_per_pop[4], fc_per_pop[5]))
  }
  
  # ---- Inject NP cascade for this year only ----
  for (i in param_var) {
    dfList_NP[[i]] <- Param_cal(
      pj = POC_AU, dlist = dfList_NP, index = i,
      frac_testing = frac_test[[as.character(max(cal_years))]],
      S_Yint = year, S_Yend = year + 1, r_Yend = year + 1,
      NPlst = NPlst, fp = fp_yr
    )
  }
  
  # Carry prison values forward
  b_pt   <- ((year + 1) - POC_AU$cabY) / POC_AU$timestep + 1
  end_dt <- ((year + 1) - POC_AU$cabY) / POC_AU$timestep
  for (i in param_var) {
    dfList_NP[[i]] <- carry_forward_prison(
      dlist_np = dfList_NP, dlist_base = zero_cascade,
      index = i, b_pt = b_pt, end_dt = end_dt
    )
  }
  
  # Re-zero everything AFTER the target year for clean single-year test
  if (year < POC_AU$cabY + endY - 1) {
    tt_next <- year_time_index(POC_AU, year + 1)
    post_range <- tt_next$t0 : dim(dfList$eta)[3]
    for (i in param_var) {
      dfList_NP[[i]][, , post_range] <- 0
      dfList_CT[[i]][, , post_range] <- dfList[[i]][, , post_range]
    }
    fc_full_local[, post_range] <- 0
  }
  
  # CT displacement for target year
  
  dfList_CT <- scale_CT_eta(dfList_CT, dfList, dfList_NP, POC_AU, year, year+1)
  dfList_CT <- scale_CT_RNA(dfList_CT, dfList, dfList_NP, POC_AU, year, year+1)
  dfList_CT <- scale_CT_ab( dfList_CT, dfList, dfList_NP, POC_AU, year, year+1)
  
  # fc matrix for target year
  for (p in seq_len(5)) {
    fc_full_local[p, tt_year$t0:tt_year$t1] <- fc_per_pop[p]
  }
  
  # ---- Run scenario ----
  Sce <- HCV_np(
    POC_AU, best_estimates, best_est_pop, disease_progress, pop_array,
    param_cascade    = dfList_CT,
    param_cascade_sc = dfList_NP,
    fib              = fib, modelrun = "UN", proj = "POC_AU",
    end_Y            = endY,
    cost = NULL, costflow = NULL, costflow_Neg = NULL,
    fc_sc            = fc_full_local
  )
  
  # ---- Extract results ----
  comm <- 1:2
  pris <- 3:5
  
  extract_events <- function(pops, label) {
    ab_pos   <- sum(Sce$newTestingAb_sc[pops, tt_year$t0:tt_year$t1], na.rm = TRUE)
    ab_neg   <- sum(Sce$newTestingAb_sc_neg[pops, tt_year$t0:tt_year$t1], na.rm = TRUE)
    ag_pos   <- sum(Sce$newTestingAg_sc[pops, tt_year$t0:tt_year$t1], na.rm = TRUE)
    ag_neg   <- sum(Sce$newTestingAg_sc_neg[pops, tt_year$t0:tt_year$t1], na.rm = TRUE)
    poct_pos <- sum(Sce$newTestingPOCT_sc[pops, tt_year$t0:tt_year$t1], na.rm = TRUE)
    poct_neg <- sum(Sce$newTestingPOCT_sc_neg[pops, tt_year$t0:tt_year$t1], na.rm = TRUE)
    tr       <- sum(Sce$newTreatment_sc[pops, tt_year$t0:tt_year$t1], na.rm = TRUE)
    
    cat(sprintf("\n--- %s, year %d ---\n", label, year))
    cat(sprintf("  Ab_pos:    %8.1f      Ab_neg:    %8.1f\n",   ab_pos,   ab_neg))
    cat(sprintf("  RNA_pos:   %8.1f      RNA_neg:   %8.1f\n",   ag_pos,   ag_neg))
    cat(sprintf("  POCT_pos:  %8.1f      POCT_neg:  %8.1f\n",   poct_pos, poct_neg))
    cat(sprintf("  Total tests:        %10.1f\n",
                ab_pos + ab_neg + ag_pos + ag_neg + poct_pos + poct_neg))
    cat(sprintf("  Treatments:         %10.1f\n", tr))
    
    tibble(setting = label, year = year,
           Ab_pos = ab_pos, Ab_neg = ab_neg,
           RNA_pos = ag_pos, RNA_neg = ag_neg,
           POCT_pos = poct_pos, POCT_neg = poct_neg,
           Total_tests = ab_pos+ab_neg+ag_pos+ag_neg+poct_pos+poct_neg,
           Treatments = tr)
  }
  
  if (verbose) cat("\n=== MODEL OUTPUTS ===\n")
  results_C <- extract_events(comm, "Community")
  results_P <- extract_events(pris, "Prison")
  
  invisible(list(Sce = Sce, results_C = results_C, results_P = results_P))
}
# =============================================================================
# RUN ONE SCENARIO ACROSS 2025-2030, CHAINED
# =============================================================================


run_scenario_chained <- function(scn_df, scn_name,
                                 init_values = NULL,
                                 tolerance = 0.03,
                                 max_iter_per_year = 20) {
  
  # init_values: list(Cov_np_C, Cov_neg_C, Cov_np_P_monthly, Cov_neg_P_monthly)
  # If NULL, use 2024 calibrated values
  if (is.null(init_values)) {
    cs_2024 <- cal_list[["2024"]]$ccal_setting
    init_values <- list(
      Cov_np_C          = cs_2024$Cov_np[cs_2024$setting == "Community"],
      Cov_neg_C         = cs_2024$Cov_neg[cs_2024$setting == "Community"],
      Cov_np_P_monthly  = cs_2024$Cov_np_monthly[cs_2024$setting == "Prison"],
      Cov_neg_P_monthly = cs_2024$Cov_neg_monthly[cs_2024$setting == "Prison"]
    )
  }
  
  cat(sprintf("\n=============================================================================\n"))
  cat(sprintf("SCENARIO: %s\n", scn_name))
  cat(sprintf("=============================================================================\n"))
  cat("Initial values (from 2024):\n")
  cat(sprintf("  Cov_np_C=%.5f  Cov_neg_C=%.5f  Cov_np_P_m=%.5f  Cov_neg_P_m=%.5f\n",
              init_values$Cov_np_C, init_values$Cov_neg_C,
              init_values$Cov_np_P_monthly, init_values$Cov_neg_P_monthly))
  
  prev_values <- init_values
  year_results <- list()
  
  for (yr in scn_df$year) {
    yr_str <- as.character(yr)
    row <- scn_df |> dplyr::filter(year == yr)
    
    cat(sprintf("\n--- Fitting year %d (target_C=%.0f target_P=%.0f) ---\n",
                yr, row$T_target_C, row$T_target_P))
    
    if (row$T_target_C == 0 && row$T_target_P == 0) {
      # No NP this year - skip fitting, set Cov values to 0
      cat("  Both targets = 0, skipping fit (Cov values = 0)\n")
      year_results[[yr_str]] <- list(
        year = yr,
        converged_C = TRUE, converged_P = TRUE,
        Cov_np_C = 0, Cov_neg_C = 0,
        Cov_np_P_monthly = 0, Cov_neg_P_monthly = 0,
        T_C = 0, T_P = 0,
        ratio_C = NA, ratio_P = NA
      )
      prev_values <- list(
        Cov_np_C = 0, Cov_neg_C = 0,
        Cov_np_P_monthly = 0, Cov_neg_P_monthly = 0
      )
      next
    }
    
    fit <- fit_year_to_target(
      year         = yr,
      prev_values  = prev_values,
      T_target_C   = row$T_target_C,
      T_target_P   = row$T_target_P,
      tolerance    = tolerance,
      max_iter_total = max_iter_per_year,
      verbose      = TRUE
    )
    
    year_results[[yr_str]] <- fit
    
    # Chain into next year
    prev_values <- list(
      Cov_np_C          = fit$Cov_np_C,
      Cov_neg_C         = fit$Cov_neg_C,
      Cov_np_P_monthly  = fit$Cov_np_P_monthly,
      Cov_neg_P_monthly = fit$Cov_neg_P_monthly
    )
  }
  
  # Summary table
  summary_df <- bind_rows(lapply(year_results, function(r) {
    tibble(
      year = r$year,
      Cov_np_C  = r$Cov_np_C,  Cov_neg_C  = r$Cov_neg_C,
      Cov_np_P_monthly = r$Cov_np_P_monthly,
      Cov_neg_P_monthly = r$Cov_neg_P_monthly,
      T_target_C = scn_df$T_target_C[scn_df$year == r$year],
      T_mod_C    = r$T_C,
      ratio_C    = r$ratio_C,
      T_target_P = scn_df$T_target_P[scn_df$year == r$year],
      T_mod_P    = r$T_P,
      ratio_P    = r$ratio_P,
      converged  = isTRUE(r$converged_C) && isTRUE(r$converged_P)
    )
  }))
  
  cat(sprintf("\n--- Summary for scenario '%s' ---\n", scn_name))
  print(summary_df, n = Inf)
  
  list(
    scenario_name = scn_name,
    summary       = summary_df,
    year_results  = year_results
  )
}


# =============================================================================
# RUN ALL SCENARIOS
# =============================================================================

scenario_fits <- list()
for (scn_name in names(scenarios)) {
  scenario_fits[[scn_name]] <- run_scenario_chained(
    scn_df  = scenarios[[scn_name]],
    scn_name = scn_name,
    tolerance = 0.03,
    max_iter_per_year = 20
  )
}


# =============================================================================
# EXTRACT FITTED Cov VALUES INTO A FLAT TABLE FOR REPRODUCIBILITY
# =============================================================================

fitted_coverages <- bind_rows(lapply(names(scenario_fits), function(scn_name) {
  sf <- scenario_fits[[scn_name]]
  sf$summary |>
    dplyr::mutate(scenario = scn_name) |>
    dplyr::select(scenario, year,
                  Cov_np_C, Cov_neg_C,
                  Cov_np_P_monthly, Cov_neg_P_monthly,
                  T_target_C, T_mod_C, ratio_C,
                  T_target_P, T_mod_P, ratio_P,
                  converged)
}))

cat("Fitted coverages across all scenarios:\n")
print(fitted_coverages, n = Inf)


# =============================================================================
# SAVE
# =============================================================================

save(
  scenario_fits,        # full nested results (Sce objects, trial logs, etc.)
  scenarios,            # scenario definitions (T_targets)
  fitted_coverages,     # flat reproducibility table
  
  file = file.path(OutputFolder, paste0(project_name, "_NPscenarios_fitted.rda"))
)

cat(sprintf("\nSaved to: %s\n",
            file.path(OutputFolder, paste0(project_name, "_NPscenarios_fitted.rda"))))

# =============================================================================
# BUILD FINAL Sce_np LIST FROM FITTED SCENARIOS
# =============================================================================
# For each scenario:
#   - Start from 2024 calibrated cascade (dfList_NP_2024, dfList_CT_NP_2024, fc_full)
#   - Inject fitted Cov values for 2025-2030 sequentially
#   - Run full HCV_np simulation
#   - Store result keyed by scenario name
# =============================================================================

build_scenario_cascade <- function(fitted_df_scn, verbose = FALSE) {
  # fitted_df_scn: rows for one scenario, years 2025-2030, with fitted Cov values
  
  zero_cascade <- lapply(dfList, function(x) x * 0)
  steps_per_year <- 1 / POC_AU$timestep
  m2a <- function(x) {
    x <- pmin(pmax(x, 0), 1)
    1 - (1 - x) ^ steps_per_year
  }
  
  # Start from 2024 calibrated state
  dfList_NP <- dfList_NP_2024
  dfList_CT <- dfList_CT_NP_2024
  fc_full_scn <- fc_full
  
  # Zero out 2025+ (override carry_forward_prison from 2024)
  tt_2025 <- year_time_index(POC_AU, 2025)
  zero_range <- tt_2025$t0 : dim(dfList$eta)[3]
  for (i in param_var) {
    dfList_NP[[i]][, , zero_range] <- 0
    dfList_CT[[i]][, , zero_range] <- dfList[[i]][, , zero_range]
  }
  fc_full_scn[, zero_range] <- 0
  
  # Inject each year 2025-2030
  for (k in seq_len(nrow(fitted_df_scn))) {
    row <- fitted_df_scn[k, ]
    yr  <- row$year
    
    # Skip zero-coverage years
    if (row$Cov_np_C == 0 && row$Cov_np_P_monthly == 0) {
      if (verbose) cat(sprintf("  %d: all-zero NP, skipping injection\n", yr))
      next
    }
    
    # Build fp_yr (follows calibration convention: P_nPWID = Cov_neg)
    Cov_np_P_annual  <- m2a(row$Cov_np_P_monthly)
    Cov_neg_P_annual <- m2a(row$Cov_neg_P_monthly)
    fp_yr <- c(
      row$Cov_np_C, row$Cov_np_C,
      Cov_np_P_annual, Cov_np_P_annual,
      Cov_neg_P_annual    # P_nPWID
    )
    fp_yr <- pmin(pmax(fp_yr, 0), 0.999999)
    
    # fc per pop
    fc_C         <- if (row$Cov_np_C > 1e-12) row$Cov_neg_C / row$Cov_np_C else 0
    fc_P_monthly <- if (row$Cov_np_P_monthly > 1e-12) row$Cov_neg_P_monthly / row$Cov_np_P_monthly else 0
    fc_per_pop   <- c(fc_C, fc_C, fc_P_monthly, fc_P_monthly, 1)
    
    if (verbose) {
      cat(sprintf("  %d: fp_yr=[%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                  yr, fp_yr[1], fp_yr[2], fp_yr[3], fp_yr[4], fp_yr[5]))
    }
    
    # Inject NP cascade
    for (i in param_var) {
      dfList_NP[[i]] <- Param_cal(
        pj = POC_AU, dlist = dfList_NP, index = i,
        frac_testing = frac_test[[as.character(max(cal_years))]],
        S_Yint = yr, S_Yend = yr + 1, r_Yend = yr + 1,
        NPlst = NPlst, fp = fp_yr
      )
    }
    
    # Carry prison values forward
    b_pt   <- ((yr + 1) - POC_AU$cabY) / POC_AU$timestep + 1
    end_dt <- ((yr + 1) - POC_AU$cabY) / POC_AU$timestep
    for (i in param_var) {
      dfList_NP[[i]] <- carry_forward_prison(
        dlist_np = dfList_NP, dlist_base = zero_cascade,
        index = i, b_pt = b_pt, end_dt = end_dt
      )
    }
    
    # CT displacement
    
    dfList_CT <- scale_CT_eta(dfList_CT, dfList, dfList_NP, POC_AU, yr, yr+1)
    dfList_CT <- scale_CT_RNA(dfList_CT, dfList, dfList_NP, POC_AU, yr, yr+1)
    dfList_CT <- scale_CT_ab( dfList_CT, dfList, dfList_NP, POC_AU, yr, yr+1)
    
    # fc matrix for this year
    tt <- year_time_index(POC_AU, yr)
    for (p in seq_len(5)) {
      fc_full_scn[p, tt$t0:tt$t1] <- fc_per_pop[p]
    }
  }
  
  # Re-zero any years AFTER the last scenario year (clean tail)
  last_yr <- max(fitted_df_scn$year)
  if (last_yr < POC_AU$cabY + endY - 1) {
    tt_after <- year_time_index(POC_AU, last_yr + 1)
    tail_range <- tt_after$t0 : dim(dfList$eta)[3]
    for (i in param_var) {
      dfList_NP[[i]][, , tail_range] <- 0
      dfList_CT[[i]][, , tail_range] <- dfList[[i]][, , tail_range]
    }
    fc_full_scn[, tail_range] <- 0
  }
  
  list(
    scenario_cascade_NP = dfList_NP,
    scenario_cascade_CT = dfList_CT,
    scenario_fc         = fc_full_scn
  )
}


# =============================================================================
# RUN HCV_np FOR EACH SCENARIO USING ITS FITTED CASCADE
# =============================================================================
endY <- 100

# Initialize cascade containers (built once, reused across cost types)
scenario_cascade    <- list()   # NP cascade per scenario
scenario_cascade_CT <- list()   # CT cascade per scenario (CT-displaced)
scenario_fc         <- list()   # fc matrix per scenario

# --- Step 1: build all scenario cascades (cost-independent) ---
cat("\n=============================================================================\n")
cat("Building cascades for all scenarios\n")
cat("=============================================================================\n")
zero_cascade_base <- lapply(dfList, function(x) x * 0)
fc_zero           <- matrix(0, nrow = nrow(fc_full), ncol = ncol(fc_full))

for (scn_name in names(scenarios)) {
  cat(sprintf("\n  Building cascade for: %s\n", scn_name))
  
  if (scn_name == "no_np") {
    # Status quo (= Sce_sq): no national program at any time
    cat("    Status quo: no NP injection ever — using baseline dfList\n")
    scenario_cascade[[scn_name]]    <- zero_cascade_base   # no scenario testing
    scenario_cascade_CT[[scn_name]] <- dfList              # routine testing, un-displaced
    scenario_fc[[scn_name]]         <- fc_zero             # no fc adjustments
    next
  }
  
  fitted_df_scn <- fitted_coverages |>
    dplyr::filter(scenario == scn_name) |>
    dplyr::arrange(year)
  
  built <- build_scenario_cascade(fitted_df_scn, verbose = TRUE)
  
  scenario_cascade[[scn_name]]    <- built$scenario_cascade_NP
  scenario_cascade_CT[[scn_name]] <- built$scenario_cascade_CT
  scenario_fc[[scn_name]]         <- built$scenario_fc
}


cat("\nAll scenario cascades built.\n")


# --- Step 2: outer cost loop, inner scenario loop ---
cost_types <- c("fixednvariable", "total", "DAAcost_reduchalf")

for (cost_type in cost_types) {
  cat(sprintf("\n=============================================================================\n"))
  cat(sprintf("Running simulations for cost type: %s\n", cost_type))
  cat(sprintf("=============================================================================\n"))
  
  # Load cost data for this cost type
  cost_dir <- switch(
    cost_type,
    "fixednvariable"    = file.path(DataFolder, "cost"),
    "total"             = file.path(DataFolder, "cost", "sensitivity_total"),
    "DAAcost_reduchalf" = file.path(DataFolder, "cost", "sensitivity_DAAcost_reduchalf")
  )
  
  files <- list.files(path = cost_dir, pattern = "\\.csv$", full.names = TRUE)
  
  costdfList <- lapply(files, function(fp) {
    df <- read.csv(fp, header = TRUE, check.names = FALSE)
    df <- df[, -1, drop = FALSE]
    as.matrix(df)
  })
  names(costdfList) <- tools::file_path_sans_ext(basename(files))
  
  cost_state <- costdfList[["state"]]
  
  costflow <- list(
    costdfList[["costFlow"]],
    costdfList[["costFlow_POCRNA"]]
  )
  
  costflow_Neg <- list(
    costdfList[["costFlow_NEG"]],
    costdfList[["costFlow_POCRNA _NEG"]]
  )
  
  # Run scenarios with this cost configuration
  Sce_np <- list()
  
  tic <- proc.time()
  
  for (scn_name in names(scenarios)) {
    cat(sprintf("\n  Running scenario: %s (cost: %s)\n", scn_name, cost_type))
    
    Sce_np[[scn_name]] <- HCV_np(
      POC_AU, best_estimates, best_est_pop, disease_progress, pop_array,
      param_cascade    = scenario_cascade_CT[[scn_name]],
      param_cascade_sc = scenario_cascade[[scn_name]],
      fib              = fib,
      modelrun         = "UN",
      proj             = "POC_AU",
      end_Y            = endY,
      cost             = costdfList,
      costflow         = costflow,
      costflow_Neg     = costflow_Neg,
      fc_sc            = scenario_fc[[scn_name]]
    )
    
    cat(sprintf("    Scenario '%s' complete.\n", scn_name))
  }
  
  toc <- proc.time() - tic
  cat(sprintf("\nCost type '%s' total time: %.1f sec\n", cost_type, toc[3]))
  
  # Save per-cost-type file
  save(
    Sce_np,                  # final simulations per scenario for this cost type
    scenario_cascade,        # NP cascade per scenario (same across cost types)
    scenario_cascade_CT,     # CT cascade per scenario
    scenario_fc,             # fc matrix per scenario
    scenarios,               # scenario definitions
    fitted_coverages,        # per-year per-scenario Cov values
    
    file = file.path(OutputFolder,
                     paste0(project_name, "Simulations_", cost_type, ".rda"))
  )
  
  cat(sprintf("Saved: %sSimulations_%s.rda\n", project_name, cost_type))
  
  rm(Sce_np, costdfList, costflow, costflow_Neg, cost_state)
  gc()
}


# =============================================================================
# VALIDATE: check that fitted values reproduce the targets
# (run on the LAST cost-type's Sce_np, since dynamics don't depend on cost)
# =============================================================================

# Load back the last saved file for validation
load(file.path(OutputFolder,
               paste0(project_name, "Simulations_", tail("DAAcost_reduchalf", 1), ".rda")))

validation <- bind_rows(lapply(names(Sce_np), function(scn_name) {
  Sce <- Sce_np[[scn_name]]
  scn_df <- scenarios[[scn_name]]
  bind_rows(lapply(scn_df$year, function(yr) {
    tt <- year_time_index(POC_AU, yr)
    t_C <- sum(Sce$newTestingAb_sc[1:2, tt$t0:tt$t1],
               Sce$newTestingAb_sc_neg[1:2, tt$t0:tt$t1],
               Sce$newTestingAg_sc[1:2, tt$t0:tt$t1],
               Sce$newTestingAg_sc_neg[1:2, tt$t0:tt$t1],
               Sce$newTestingPOCT_sc[1:2, tt$t0:tt$t1],
               Sce$newTestingPOCT_sc_neg[1:2, tt$t0:tt$t1], na.rm = TRUE)
    t_P <- sum(Sce$newTestingAb_sc[3:5, tt$t0:tt$t1],
               Sce$newTestingAb_sc_neg[3:5, tt$t0:tt$t1],
               Sce$newTestingAg_sc[3:5, tt$t0:tt$t1],
               Sce$newTestingAg_sc_neg[3:5, tt$t0:tt$t1],
               Sce$newTestingPOCT_sc[3:5, tt$t0:tt$t1],
               Sce$newTestingPOCT_sc_neg[3:5, tt$t0:tt$t1], na.rm = TRUE)
    row <- scn_df |> dplyr::filter(year == yr)
    tibble(
      scenario = scn_name, year = yr,
      T_target_C = row$T_target_C, T_mod_C = t_C,
      ratio_C = if (row$T_target_C > 0) t_C / row$T_target_C else NA,
      T_target_P = row$T_target_P, T_mod_P = t_P,
      ratio_P = if (row$T_target_P > 0) t_P / row$T_target_P else NA
    )
  }))
}))

cat("\n=============================================================================\n")
cat("VALIDATION: chained scenario tests vs targets\n")
cat("=============================================================================\n")
print(validation, n = Inf)


