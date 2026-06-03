# ══════════════════════════════════════════════════════════════════════════════
# RL Option Coverage v3: Multiplier-on-Coverage Action
# ──────────────────────────────────────────────────────────────────────────────
# Action per year (4 dims, all multipliers ≥ 1.0):
#   m_np_C   — multiplier on Cov_np_C
#   m_neg_C  — multiplier on Cov_neg_C
#   m_np_P   — multiplier on Cov_np_P_monthly
#   m_neg_P  — multiplier on Cov_neg_P_monthly
# Full episode: 5 years × 4 = 20-dim action.
#
# Hard constraints (clipped before passing to fit):
#   • Cov_*_new ≥ Cov_*_prev      (multiplier ≥ 1.0)
#   • Cov_np ≤ 0.99               (community annual + prison monthly)
#   • Cov_neg ≤ 0.99
#   • Cov_neg ≤ Cov_np            (negatives never tested above positives)
#
# Soft constraint (penalty in reward):
#   • Prison total tests ≤ 60,000 / year
#
# Implicit "fc drift": because both Cov_np and Cov_neg are non-decreasing
# (and Cov_neg ≤ Cov_np), fc = Cov_neg / Cov_np can only drift in a bounded
# way; no separate fc enforcement.
#
# Episode = one HCVMSM_cpp call. Incidence = 100 × newInf / (2 × N_midyear)
# (existing P_PWID formula). Trajectory: linear from inc_2025 → 2.0 by 2030.
# ══════════════════════════════════════════════════════════════════════════════

load("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/POC_AU/02. Output/POC_AU.rda")
load("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/POC_AU/02. Output/POC_AUcali.rda")
load("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/POC_AU/02. Output/POC_AUcali_timev.rda")
load("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/POC_AU/02. Output/POC_AU_NPcal_2024.rda")
# Load fitted 5-scenario state (needed for the 2025 anchor)
load("/Users/jjwu/Projects/Simplified-HCV-testing-model/Projects/POC_prisons_scaleup/02. Output/POC_prisons_scaleupSimulations_epi.rda")

source("/Users/jjwu/Projects/Simplified-HCV-testing-model/03. Code/Functions/HCV_model.R")
source("/Users/jjwu/Projects/Simplified-HCV-testing-model/03. Code/Functions/plotFunctions.R")
source("~/Projects/Simplified-HCV-testing-model/Projects/POC_prisons_scaleup/Param_np_fc_cal.R")

library(Rcpp); library(RcppArmadillo); library(parallel); library(dplyr); library(tidyr)
Rcpp::sourceCpp("~/Projects/Simplified-HCV-testing-model/Projects/POC_prisons_scaleup/hcvmsm_core.cpp")
source(        "~/Projects/Simplified-HCV-testing-model/Projects/POC_prisons_scaleup/HCVMSM_cpp.R")

if (!exists("endY")) endY <- 100

required_objs <- c("POC_AU","best_estimates","best_est_pop","disease_progress",
                   "pop_array","fib","dfList","dfList_NP_2024","dfList_CT_NP_2024",
                   "fc_full","param_var","frac_test","NPlst","endY",
                   "scenario_cascade","scenario_cascade_CT","scenario_fc",
                   "fitted_coverages","Sce_np")
missing_objs <- required_objs[!sapply(required_objs, exists)]
if (length(missing_objs) > 0) stop("Missing: ", paste(missing_objs, collapse=", "))

required_fns <- c("Param_cal","carry_forward_prison","scale_CT_eta","scale_CT_RNA",
                  "scale_CT_ab","popResults_MidYear","modres.flow.t",
                  "HCVMSM_cpp","hcvmsm_loop_cpp")
missing_fns <- required_fns[!sapply(required_fns, exists)]
if (length(missing_fns) > 0) stop("Missing functions: ", paste(missing_fns, collapse=", "))
cat("All required objects and functions loaded ✓\n")

# ══════════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════
RL_YEARS         <- 2026:2030
ANCHOR_YEAR      <- 2025
TARGET_INC_P     <- 2.0
P_PWID_POP       <- "P_PWID"
PRISON_POPS_IDX  <- 3:5
PRISON_TESTS_CAP <- 60000

COV_MAX        <- 0.99      # cap for prison Cov_np_monthly and any Cov_neg
COV_MAX_C_NP   <- 0.60      # community Cov_np cap (new — was 0.99)
M_MIN          <- 1.0
M_THRESHOLD    <- 0.20      # Cov_np threshold for switching m cap
M_MAX_LOW      <- 2.0       # m cap when Cov_np_prev < threshold (can ~double)
M_MAX_HIGH     <- 1.5       # m cap when Cov_np_prev ≥ threshold (gentler)
M_MAX          <- M_MAX_LOW # used by policy sampling Gaussian clip

# Reward
GAMMA        <- 0.95
ALPHA_TRAJ   <- 1.0
ALPHA_PRISON <- 5.0    # soft penalty for prison >60k
TERMINAL_B   <- 10.0

# Training
N_CORES      <- 6
N_EPISODES   <- 500
BATCH_SIZE   <- 20
LR           <- 0.05
ENTROPY_COEF <- 0.1

# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════
year_time_index <- function(pj, year) {
  steps_per_year <- 1 / pj$timestep
  k  <- year - pj$cabY + 1
  t0 <- (k - 1) * steps_per_year + 1
  t1 <-  k      * steps_per_year
  list(t0 = t0, t1 = t1, steps_per_year = steps_per_year)
}

safe_num <- function(x, default = 0) {
  if (is.null(x) || length(x) == 0) return(as.numeric(default))
  if (is.list(x)) x <- x[[1]]
  if (is.null(x) || length(x) == 0 || is.na(x)) return(as.numeric(default))
  as.numeric(x)
}

safe_extract <- function(results_list, field, default = 0) {
  vapply(results_list, function(r) safe_num(r[[field]], default), numeric(1))
}

m2a <- function(x, dt) {
  x <- pmin(pmax(x, 0), 1)
  1 - (1 - x)^(1/dt)
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# ══════════════════════════════════════════════════════════════════════════════
# ACTION → COVERAGES
# ══════════════════════════════════════════════════════════════════════════════
# action_4: c(m_np_C, m_neg_C, m_np_P, m_neg_P), each ≥ 1
# prev_cov: list with Cov_np_C, Cov_neg_C, Cov_np_P_monthly, Cov_neg_P_monthly
#
# Hard clips applied (in order):
#   1. m ≥ 1.0                  → Cov_new ≥ Cov_prev
#   2. Cov_np_new ≤ COV_MAX
#   3. Cov_neg_new ≤ COV_MAX
#   4. Cov_neg_new ≤ Cov_np_new (negatives never above positives)
#   5. Cov_neg_new ≥ Cov_neg_prev (floor)
#   6. Cov_np_new  ≥ Cov_np_prev  (floor)
#   7. fc_new ≥ fc_prev (Option A: raise Cov_neg to fc_prev * Cov_np;
#      if that exceeds COV_MAX, cap Cov_np at COV_MAX / fc_prev so the
#      ratio can still be satisfied)
apply_action_to_coverages <- function(action_4, prev_cov) {
  m_np_C  <- max(M_MIN, min(M_MAX, action_4[1]))
  m_neg_C <- max(M_MIN, min(M_MAX, action_4[2]))
  m_np_P  <- max(M_MIN, min(M_MAX, action_4[3]))
  m_neg_P <- max(M_MIN, min(M_MAX, action_4[4]))
  
  # Apply multipliers, then clip to coverage cap
  Cov_np_C   <- min(prev_cov$Cov_np_C          * m_np_C,  COV_MAX)
  Cov_neg_C  <- min(prev_cov$Cov_neg_C         * m_neg_C, COV_MAX)
  Cov_np_P   <- min(prev_cov$Cov_np_P_monthly  * m_np_P,  COV_MAX)
  Cov_neg_P  <- min(prev_cov$Cov_neg_P_monthly * m_neg_P, COV_MAX)
  
  # Enforce Cov_neg ≤ Cov_np
  Cov_neg_C  <- min(Cov_neg_C, Cov_np_C)
  Cov_neg_P  <- min(Cov_neg_P, Cov_np_P)
  
  # Floor: never below previous year
  Cov_np_C   <- max(Cov_np_C,  prev_cov$Cov_np_C)
  Cov_neg_C  <- max(Cov_neg_C, prev_cov$Cov_neg_C)
  Cov_np_P   <- max(Cov_np_P,  prev_cov$Cov_np_P_monthly)
  Cov_neg_P  <- max(Cov_neg_P, prev_cov$Cov_neg_P_monthly)
  
  # If after floor enforcement Cov_neg > Cov_np (degenerate edge case), drag
  # Cov_np up
  if (Cov_neg_C > Cov_np_C) Cov_np_C <- min(Cov_neg_C, COV_MAX)
  if (Cov_neg_P > Cov_np_P) Cov_np_P <- min(Cov_neg_P, COV_MAX)
  
  # ── fc ≥ fc_prev (Option A: raise Cov_neg) ────────────────────────────────
  fc_prev_C <- if (prev_cov$Cov_np_C > 1e-12)
    prev_cov$Cov_neg_C / prev_cov$Cov_np_C else 0
  fc_prev_P <- if (prev_cov$Cov_np_P_monthly > 1e-12)
    prev_cov$Cov_neg_P_monthly / prev_cov$Cov_np_P_monthly else 0
  
  # Community: enforce fc_C_new ≥ fc_prev_C
  if (Cov_np_C > 1e-12) {
    Cov_neg_C_needed <- fc_prev_C * Cov_np_C
    if (Cov_neg_C_needed > COV_MAX) {
      # Can't raise Cov_neg high enough — cap Cov_np to make ratio achievable
      Cov_np_C  <- COV_MAX / max(fc_prev_C, 1e-12)
      Cov_neg_C <- COV_MAX
    } else if (Cov_neg_C < Cov_neg_C_needed) {
      Cov_neg_C <- Cov_neg_C_needed
    }
  }
  
  # Prison: enforce fc_P_new ≥ fc_prev_P
  if (Cov_np_P > 1e-12) {
    Cov_neg_P_needed <- fc_prev_P * Cov_np_P
    if (Cov_neg_P_needed > COV_MAX) {
      Cov_np_P  <- COV_MAX / max(fc_prev_P, 1e-12)
      Cov_neg_P <- COV_MAX
    } else if (Cov_neg_P < Cov_neg_P_needed) {
      Cov_neg_P <- Cov_neg_P_needed
    }
  }
  
  # Final safety: Cov_neg ≤ Cov_np still holds (raising Cov_neg to fc_prev × Cov_np
  # satisfies this since fc_prev ≤ 1)
  Cov_neg_C  <- min(Cov_neg_C, Cov_np_C)
  Cov_neg_P  <- min(Cov_neg_P, Cov_np_P)
  
  list(
    Cov_np_C          = Cov_np_C,
    Cov_neg_C         = Cov_neg_C,
    Cov_np_P_monthly  = Cov_np_P,
    Cov_neg_P_monthly = Cov_neg_P,
    m_np_C  = m_np_C,  m_neg_C = m_neg_C,
    m_np_P  = m_np_P,  m_neg_P = m_neg_P
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# BUILD EPISODE CASCADE
# ══════════════════════════════════════════════════════════════════════════════
# Mirrors build_scenario_cascade from the 5-scenario pipeline.
build_episode_cascade <- function(coverages_by_year, pj,
                                  dfList_in, frac_test_in,
                                  NPlst_in, param_var_in, endY_in) {
  
  dt <- pj$timestep
  zero_cascade <- lapply(dfList_in, function(x) x * 0)
  
  # Start from 2025-anchored scenario cascade
  dfList_NP   <- scenario_cascade[["sustained_equal_split"]]
  dfList_CT   <- scenario_cascade_CT[["sustained_equal_split"]]
  fc_full_scn <- scenario_fc[["sustained_equal_split"]]
  
  # Zero out 2026+
  tt_2026 <- year_time_index(pj, 2026)
  zero_range <- tt_2026$t0 : dim(dfList_in$eta)[3]
  for (i in param_var_in) {
    dfList_NP[[i]][, , zero_range] <- 0
    dfList_CT[[i]][, , zero_range] <- dfList_in[[i]][, , zero_range]
  }
  fc_full_scn[, zero_range] <- 0
  
  # Inject each year
  for (yr in RL_YEARS) {
    cov_yr <- coverages_by_year[[as.character(yr)]]
    
    Cov_np_P_annual  <- m2a(cov_yr$Cov_np_P_monthly,  dt)
    Cov_neg_P_annual <- m2a(cov_yr$Cov_neg_P_monthly, dt)
    
    fp_yr <- c(
      cov_yr$Cov_np_C, cov_yr$Cov_np_C,
      Cov_np_P_annual, Cov_np_P_annual,
      Cov_neg_P_annual)
    fp_yr <- pmin(pmax(fp_yr, 0), 0.999999)
    
    fc_C <- if (cov_yr$Cov_np_C > 1e-12)
      cov_yr$Cov_neg_C / cov_yr$Cov_np_C else 0
    fc_P <- if (cov_yr$Cov_np_P_monthly > 1e-12)
      cov_yr$Cov_neg_P_monthly / cov_yr$Cov_np_P_monthly else 0
    fc_per_pop <- c(fc_C, fc_C, fc_P, fc_P, 1.0)
    
    for (i in param_var_in) {
      dfList_NP[[i]] <- Param_cal(
        pj = pj, dlist = dfList_NP, index = i,
        frac_testing = frac_test_in[[as.character(yr)]],
        S_Yint = yr, S_Yend = yr + 1, r_Yend = yr + 1,
        NPlst = NPlst_in, fp = fp_yr)
    }
    b_pt   <- ((yr + 1) - pj$cabY) / dt + 1
    end_dt <- ((yr + 1) - pj$cabY) / dt
    for (i in param_var_in) {
      dfList_NP[[i]] <- carry_forward_prison(
        dlist_np = dfList_NP, dlist_base = zero_cascade,
        index = i, b_pt = b_pt, end_dt = end_dt)
    }
    dfList_CT <- scale_CT_eta(dfList_CT, dfList_in, dfList_NP, pj, yr, yr+1)
    dfList_CT <- scale_CT_RNA(dfList_CT, dfList_in, dfList_NP, pj, yr, yr+1)
    dfList_CT <- scale_CT_ab( dfList_CT, dfList_in, dfList_NP, pj, yr, yr+1)
    
    tt <- year_time_index(pj, yr)
    for (p in seq_len(5)) {
      fc_full_scn[p, tt$t0:tt$t1] <- fc_per_pop[p]
    }
  }
  
  # Tail cleanup
  if (max(RL_YEARS) < pj$cabY + endY_in - 1) {
    tt_after <- year_time_index(pj, max(RL_YEARS) + 1)
    tail_range <- tt_after$t0 : dim(dfList_in$eta)[3]
    for (i in param_var_in) {
      dfList_NP[[i]][, , tail_range] <- 0
      dfList_CT[[i]][, , tail_range] <- dfList_in[[i]][, , tail_range]
    }
    fc_full_scn[, tail_range] <- 0
  }
  
  list(dfList_NP = dfList_NP, dfList_CT = dfList_CT, fc_full = fc_full_scn)
}

# ══════════════════════════════════════════════════════════════════════════════
# INCIDENCE (P_PWID: 100 × newInf / (2 × N_midyear))
# ══════════════════════════════════════════════════════════════════════════════
inc_PPWID_year <- function(Sce, pj, year, endY_in) {
  yr_idx <- year - pj$cabY + 1
  tt     <- year_time_index(pj, year)
  new_inf_P <- sum(Sce$newInfections["P_PWID", tt$t0:tt$t1], na.rm = TRUE)
  N_mid <- popResults_MidYear(
    pj, Sce,
    Population = P_PWID_POP,
    Disease_prog = NULL, Cascade = NULL, param = NULL,
    endYear = endY_in) %>%
    ungroup() %>%
    dplyr::filter(year == yr_idx) %>%
    dplyr::pull(best) %>% sum()
  if (N_mid <= 0) return(NA_real_)
  100 * new_inf_P / (2 * N_mid)
}

prison_tests_year <- function(Sce, pj, year) {
  tt <- year_time_index(pj, year)
  flows <- c("newTestingAb_sc","newTestingAg_sc","newTestingPOCT_sc",
             "newTestingAb_sc_neg","newTestingAg_sc_neg","newTestingPOCT_sc_neg")
  sum(sapply(flows, function(f) sum(Sce[[f]][PRISON_POPS_IDX, tt$t0:tt$t1], na.rm=TRUE)))
}

# ══════════════════════════════════════════════════════════════════════════════
# REWARD
# ══════════════════════════════════════════════════════════════════════════════
make_trajectory <- function(inc_2025) {
  yrs <- RL_YEARS
  inc_2025 + (TARGET_INC_P - inc_2025) * (yrs - 2025) / 5
}

compute_episode_reward <- function(inc_trajectory, prison_vec, inc_2025) {
  inc_targets <- make_trajectory(inc_2025)
  
  # Trajectory penalty
  r_traj <- 0
  for (k in seq_along(RL_YEARS)) {
    overshoot <- max(0, inc_trajectory[k] - inc_targets[k]) / max(inc_targets[k], 1e-6)
    r_traj    <- r_traj - ALPHA_TRAJ * (GAMMA^(k-1)) * overshoot
  }
  
  # Prison soft penalty (over 60k)
  r_prison <- 0
  for (k in seq_along(RL_YEARS)) {
    over_cap <- max(0, prison_vec[k] - PRISON_TESTS_CAP) / PRISON_TESTS_CAP
    r_prison <- r_prison - ALPHA_PRISON * (GAMMA^(k-1)) * over_cap
  }
  
  # Terminal
  inc_2030    <- inc_trajectory[length(inc_trajectory)]
  eliminated  <- as.numeric(inc_2030 <= TARGET_INC_P)
  r_terminal  <- (GAMMA^4) * (TERMINAL_B * eliminated - TERMINAL_B * (1 - eliminated))
  
  list(r_total = r_traj + r_prison + r_terminal,
       r_traj = r_traj, r_prison = r_prison, r_terminal = r_terminal,
       eliminated = eliminated, inc_2030 = inc_2030)
}

# ══════════════════════════════════════════════════════════════════════════════
# ANCHOR
# ══════════════════════════════════════════════════════════════════════════════
get_initial_coverages <- function() {
  cov_2025 <- fitted_coverages |>
    dplyr::filter(scenario == "sustained_equal_split", year == 2025) |>
    dplyr::slice(1)
  list(
    Cov_np_C          = cov_2025$Cov_np_C,
    Cov_neg_C         = cov_2025$Cov_neg_C,
    Cov_np_P_monthly  = cov_2025$Cov_np_P_monthly,
    Cov_neg_P_monthly = cov_2025$Cov_neg_P_monthly
  )
}

compute_inc_2025 <- function() {
  inc_PPWID_year(Sce_np[["sustained_equal_split"]], POC_AU, ANCHOR_YEAR, endY)
}

# ══════════════════════════════════════════════════════════════════════════════
# RUN ONE EPISODE
# ══════════════════════════════════════════════════════════════════════════════
# action_vec: length 20, indexing [4(k-1)+1..4k] = (m_np_C, m_neg_C, m_np_P, m_neg_P)
run_episode <- function(action_vec, inc_2025) {
  pj      <- POC_AU
  endY_in <- endY
  
  # Chain coverages through years
  prev_cov <- get_initial_coverages()
  cov_by_year <- list()
  cov_traj    <- list()
  
  for (k in seq_along(RL_YEARS)) {
    yr  <- RL_YEARS[k]
    a4  <- action_vec[(4*(k-1)+1):(4*k)]
    new_cov <- apply_action_to_coverages(a4, prev_cov)
    cov_by_year[[as.character(yr)]] <- new_cov[c("Cov_np_C","Cov_neg_C",
                                                 "Cov_np_P_monthly",
                                                 "Cov_neg_P_monthly")]
    cov_traj[[as.character(yr)]] <- new_cov
    prev_cov <- cov_by_year[[as.character(yr)]]
  }
  
  # Build cascade and run one sim
  Sce <- tryCatch({
    built <- build_episode_cascade(
      cov_by_year, pj,
      dfList_in = dfList, frac_test_in = frac_test,
      NPlst_in = NPlst, param_var_in = param_var, endY_in = endY_in)
    HCVMSM_cpp(
      pj, best_estimates, best_est_pop, disease_progress, pop_array,
      param_cascade = built$dfList_CT, param_cascade_sc = built$dfList_NP,
      fib = fib, modelrun = "UN", proj = "POC_AU",
      end_Y = endY_in, fc_sc = built$fc_full)
  }, error = function(e) {
    cat("  [Episode error]:", conditionMessage(e), "\n"); NULL
  })
  
  if (is.null(Sce)) return(list(total_return = -50, eliminated = 0, error = TRUE))
  
  inc_traj   <- sapply(RL_YEARS, function(y) inc_PPWID_year(Sce, pj, y, endY_in))
  prison_vec <- sapply(RL_YEARS, function(y) prison_tests_year(Sce, pj, y))
  
  reward <- compute_episode_reward(inc_traj, prison_vec, inc_2025)
  
  list(
    total_return = reward$r_total,
    r_traj       = reward$r_traj,
    r_prison     = reward$r_prison,
    r_terminal   = reward$r_terminal,
    eliminated   = reward$eliminated,
    inc_2030     = reward$inc_2030,
    inc_traj     = inc_traj,
    prison_vec   = prison_vec,
    cov_traj     = cov_traj,
    action_vec   = action_vec,
    error        = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# PERSISTENT CLUSTER
# ══════════════════════════════════════════════════════════════════════════════
setup_cluster <- function(n_cores = N_CORES) {
  cl <- makeCluster(n_cores)
  
  # Export ONLY what the workers need (exclude giant unused objects)
  # We deliberately exclude Sce_np here — it's only used in the main process
  # for compute_inc_2025() before training starts.
  needed_objs <- c("POC_AU","best_estimates","best_est_pop","disease_progress",
                   "pop_array","fib","dfList","dfList_NP_2024","dfList_CT_NP_2024",
                   "fc_full","param_var","frac_test","NPlst","endY",
                   "scenario_cascade","scenario_cascade_CT","scenario_fc",
                   "fitted_coverages",
                   # Constants
                   "RL_YEARS","ANCHOR_YEAR","TARGET_INC_P","P_PWID_POP",
                   "PRISON_POPS_IDX","PRISON_TESTS_CAP","COV_MAX",
                   "M_MIN","M_MAX","GAMMA","ALPHA_TRAJ","ALPHA_PRISON","TERMINAL_B")
  needed_fns <- c("year_time_index","safe_num","m2a",
                  "apply_action_to_coverages","build_episode_cascade",
                  "inc_PPWID_year","prison_tests_year",
                  "make_trajectory","compute_episode_reward",
                  "get_initial_coverages","run_episode",
                  "Param_cal","carry_forward_prison",
                  "scale_CT_eta","scale_CT_RNA","scale_CT_ab",
                  "popResults_MidYear","modres.flow.t",
                  "HCVMSM_cpp")
  cat("Exporting", length(needed_fns), "fns +", length(needed_objs), "objs...\n")
  t0 <- Sys.time()
  clusterExport(cl, varlist = c(needed_fns, needed_objs), envir = .GlobalEnv)
  cat("  done in", round(as.numeric(Sys.time()-t0, units="secs"), 1), "sec\n")
  
  # Source C++ engine + project helpers on workers.
  # Use absolute paths — '~' expansion can be unreliable on parallel workers.
  # We source ALL the project R files that the main process sources, so
  # internal helpers (scale_CT_key, MidyearIndex, etc.) are available.
  clusterEvalQ(cl, {
    library(dplyr); library(tidyr); library(Rcpp); library(RcppArmadillo)
    Rcpp::sourceCpp("/Users/jjwu/Projects/Simplified-HCV-testing-model/Projects/POC_prisons_scaleup/hcvmsm_core.cpp")
    source(        "/Users/jjwu/Projects/Simplified-HCV-testing-model/Projects/POC_prisons_scaleup/HCVMSM_cpp.R")
    source(        "/Users/jjwu/Projects/Simplified-HCV-testing-model/Projects/POC_prisons_scaleup/Param_np_fc_cal.R")
    source(        "/Users/jjwu/Projects/Simplified-HCV-testing-model/03. Code/Functions/HCV_model.R")
    source(        "/Users/jjwu/Projects/Simplified-HCV-testing-model/03. Code/Functions/plotFunctions.R")
    `%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
    NULL
  })
  
  # Verify (include internal helpers — these are the ones we've discovered
  # going silently missing in past sessions)
  check <- clusterEvalQ(cl, {
    needed <- c("HCVMSM_cpp","hcvmsm_loop_cpp","Param_cal","run_episode",
                "scenario_cascade","POC_AU","dfList","pop_array",
                "scale_CT_key","scale_CT_eta","scale_CT_RNA","scale_CT_ab",
                "carry_forward_prison","popResults_MidYear","MidyearIndex",
                "modres.flow.t")
    missing <- needed[!sapply(needed, exists)]
    if (length(missing) > 0) paste("MISSING:", paste(missing, collapse=","))
    else "OK"
  })
  bad <- check[check != "OK"]
  if (length(bad) > 0) {
    stopCluster(cl)
    stop("Worker setup failed: ", paste(unique(unlist(bad)), collapse="; "))
  }
  cat("Worker check: all functions/objects present ✓\n")
  cl
}

run_episodes_on_cluster <- function(cl, action_vecs, inc_2025) {
  local_env <- environment()
  clusterExport(cl, varlist = "inc_2025", envir = local_env)
  parLapply(cl, action_vecs, function(a) run_episode(a, inc_2025))
}

# ══════════════════════════════════════════════════════════════════════════════
# SYNC TO WORKERS
# ══════════════════════════════════════════════════════════════════════════════
# Re-push the current values of constants and functions from .GlobalEnv to all
# workers. Call this AFTER you change anything in the main session
# (e.g., updated PRISON_TESTS_CAP, ALPHA_PRISON, compute_episode_reward, etc.)
# and BEFORE training, so workers don't run with stale copies.
sync_to_workers <- function(cl, also = character(0)) {
  # Constants the reward/episode logic reads
  consts <- c("RL_YEARS","ANCHOR_YEAR","TARGET_INC_P","P_PWID_POP",
              "PRISON_POPS_IDX","PRISON_TESTS_CAP","COV_MAX",
              "M_MIN","M_MAX","GAMMA","ALPHA_TRAJ","ALPHA_PRISON","TERMINAL_B")
  # Functions that contain logic worth re-syncing
  fns <- c("apply_action_to_coverages","build_episode_cascade",
           "inc_PPWID_year","prison_tests_year",
           "make_trajectory","compute_episode_reward",
           "get_initial_coverages","run_episode",
           "year_time_index","m2a","safe_num")
  varlist <- unique(c(consts, fns, also))
  
  # Filter to objects that actually exist in .GlobalEnv (avoids error if any missing)
  present <- varlist[sapply(varlist, exists, envir = .GlobalEnv)]
  absent  <- setdiff(varlist, present)
  if (length(absent) > 0) {
    cat("  (skipping non-existent:", paste(absent, collapse=", "), ")\n")
  }
  
  clusterExport(cl, varlist = present, envir = .GlobalEnv)
  cat("Synced", length(present), "objects to workers ✓\n")
  invisible(present)
}

# ══════════════════════════════════════════════════════════════════════════════
# SAC TRAINING
# ══════════════════════════════════════════════════════════════════════════════
# Policy: 20-dim Gaussian. Mean init at 1.0 (no change from prev = Phase-1 default).
train_sac_gaussian <- function(inc_2025,
                               n_episodes   = N_EPISODES,
                               n_cores      = N_CORES,
                               batch_size   = BATCH_SIZE,
                               lr           = LR,
                               entropy_coef = ENTROPY_COEF,
                               top_k        = 20,
                               cl           = NULL) {
  cat("\n=== SAC Training (Coverage v3: 4 multipliers/year) ===\n")
  cat("Anchor: P_PWID inc_2025 =", round(inc_2025, 3), "per 100 PY\n")
  cat("Target: P_PWID inc_2030 ≤", TARGET_INC_P, "per 100 PY\n")
  cat("Trajectory:", round(make_trajectory(inc_2025), 2), "\n")
  cat("Tracking top", top_k, "episodes by return.\n\n")
  
  own_cluster <- is.null(cl)
  if (own_cluster) {
    cl <- setup_cluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
  }
  
  policy_mean    <- rep(1.0, 20)              # all multipliers = 1.0 (no change)
  policy_log_std <- rep(log(0.3), 20)
  best_return    <- -Inf
  best_action    <- NULL
  history        <- data.frame()
  top_k_list     <- list()   # list of full episode results, sorted desc by return
  n_batches      <- ceiling(n_episodes / batch_size)
  
  for (batch in seq_len(n_batches)) {
    cat("── Batch", batch, "/", n_batches, "──\n")
    policy_std  <- exp(policy_log_std)
    action_vecs <- lapply(seq_len(batch_size), function(i) {
      vec <- rnorm(20, mean = policy_mean, sd = policy_std)
      pmax(M_MIN, pmin(M_MAX, vec))
    })
    
    t_batch <- system.time({
      results <- run_episodes_on_cluster(cl, action_vecs, inc_2025)
    })
    
    returns    <- safe_extract(results, "total_return", default = -50)
    eliminated <- safe_extract(results, "eliminated",   default = 0)
    inc_2030s  <- safe_extract(results, "inc_2030",     default = NA)
    
    cat(sprintf("  returns: mean=%.2f  max=%.2f  elim=%d/%d  time=%.1fs\n",
                mean(returns, na.rm=TRUE), max(returns, na.rm=TRUE),
                sum(eliminated), batch_size, t_batch["elapsed"]))
    
    # Policy gradient
    advantages   <- returns - mean(returns, na.rm = TRUE)
    grad_mean    <- rep(0, 20); grad_log_std <- rep(0, 20)
    for (i in seq_along(results)) {
      vec <- action_vecs[[i]]; adv <- advantages[i]
      g_mean    <- (vec - policy_mean) / (policy_std^2)
      g_log_std <- ((vec - policy_mean)^2 / policy_std^2) - 1
      grad_mean    <- grad_mean    + adv * g_mean
      grad_log_std <- grad_log_std + adv * g_log_std + entropy_coef
    }
    grad_mean    <- grad_mean    / batch_size
    grad_log_std <- grad_log_std / batch_size
    
    policy_mean    <- pmax(M_MIN, pmin(M_MAX, policy_mean + lr * grad_mean))
    policy_log_std <- pmax(-3, pmin(0, policy_log_std + lr * grad_log_std))
    
    # Track single best (for backward-compat)
    best_idx <- which.max(returns)
    if (returns[best_idx] > best_return) {
      best_return <- returns[best_idx]
      best_action <- action_vecs[[best_idx]]
      cat("  ** new best return:", round(best_return,3),
          " | inc_2030 =", round(inc_2030s[best_idx], 2), "\n")
    }
    
    # ── Update top-K list ─────────────────────────────────────────────────
    # Each entry stores the full episode result (cov_traj, inc_traj, prison_vec,
    # action_vec, return, etc.) plus the batch and episode-in-batch index.
    for (i in seq_along(results)) {
      r <- results[[i]]
      if (isTRUE(r$error) || is.null(r$total_return)) next
      entry <- list(
        return       = r$total_return,
        eliminated   = r$eliminated,
        inc_2030     = r$inc_2030,
        inc_traj     = r$inc_traj,
        prison_vec   = r$prison_vec,
        cov_traj     = r$cov_traj,
        action_vec   = r$action_vec,
        episode_idx  = (batch - 1) * batch_size + i,
        batch_idx    = batch
      )
      top_k_list[[length(top_k_list) + 1L]] <- entry
    }
    # Sort by return desc, keep top_k
    if (length(top_k_list) > top_k) {
      rets <- sapply(top_k_list, function(e) e$return)
      keep <- order(rets, decreasing = TRUE)[seq_len(top_k)]
      top_k_list <- top_k_list[keep]
    }
    
    history <- rbind(history, data.frame(
      episode    = (batch-1)*batch_size + seq_along(results),
      batch      = batch, ret = returns,
      eliminated = eliminated, inc_2030 = inc_2030s))
  }
  
  # Final sort of top-K
  if (length(top_k_list) > 0) {
    rets <- sapply(top_k_list, function(e) e$return)
    top_k_list <- top_k_list[order(rets, decreasing = TRUE)]
  }
  
  list(policy_mean = policy_mean, policy_log_std = policy_log_std,
       best_return = best_return, best_action = best_action,
       top_k       = top_k_list,
       history     = history)
}

# ══════════════════════════════════════════════════════════════════════════════
# INSPECT TOP-K POLICIES
# ══════════════════════════════════════════════════════════════════════════════
# Returns a tidy data.frame: one row per (rank × year) for the top-K episodes.
# Columns: rank, episode_idx, return, year, Cov_np_C, Cov_neg_C,
#          Cov_np_P_monthly, Cov_neg_P_monthly, fc_C, fc_P, inc_P_PWID,
#          prison_tests.
summarise_top_k <- function(sac_result, k = NULL) {
  top <- sac_result$top_k
  if (is.null(top) || length(top) == 0) stop("sac_result$top_k is empty.")
  if (is.null(k)) k <- length(top)
  k <- min(k, length(top))
  
  rows <- list()
  for (rank_i in seq_len(k)) {
    ep <- top[[rank_i]]
    for (j in seq_along(RL_YEARS)) {
      yr <- RL_YEARS[j]
      cv <- ep$cov_traj[[as.character(yr)]]
      rows[[length(rows) + 1L]] <- data.frame(
        rank              = rank_i,
        episode_idx       = ep$episode_idx,
        return            = ep$return,
        eliminated_2030   = ep$eliminated,
        year              = yr,
        Cov_np_C          = cv$Cov_np_C,
        Cov_neg_C         = cv$Cov_neg_C,
        Cov_np_P_monthly  = cv$Cov_np_P_monthly,
        Cov_neg_P_monthly = cv$Cov_neg_P_monthly,
        fc_C              = if (cv$Cov_np_C > 0) cv$Cov_neg_C / cv$Cov_np_C else NA,
        fc_P              = if (cv$Cov_np_P_monthly > 0)
          cv$Cov_neg_P_monthly / cv$Cov_np_P_monthly else NA,
        inc_P_PWID        = ep$inc_traj[j],
        prison_tests      = ep$prison_vec[j]
      )
    }
  }
  do.call(rbind, rows)
}

# Print a compact one-line-per-episode summary
print_top_k_summary <- function(sac_result, k = NULL) {
  top <- sac_result$top_k
  if (is.null(k)) k <- length(top)
  k <- min(k, length(top))
  cat(sprintf("\n── Top %d policies by return ──\n", k))
  cat(sprintf("%-6s %-12s %-9s %-9s %-30s\n",
              "rank", "return", "elim_2030", "inc_2030", "prison_tests (2026..2030)"))
  for (i in seq_len(k)) {
    ep <- top[[i]]
    cat(sprintf("%-6d %-12.3f %-9d %-9.2f %s\n",
                i, ep$return, ep$eliminated, ep$inc_2030,
                paste(round(ep$prison_vec), collapse=" ")))
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# EVALUATE
# ══════════════════════════════════════════════════════════════════════════════
evaluate_policy <- function(policy_mean, inc_2025, n_eval = 10, cl = NULL) {
  cat("\n=== Evaluating learned policy ===\n")
  action_vecs <- rep(list(policy_mean), n_eval)
  
  own_cluster <- is.null(cl)
  if (own_cluster) {
    cl <- setup_cluster(min(n_eval, N_CORES))
    on.exit(stopCluster(cl), add = TRUE)
  }
  
  results <- run_episodes_on_cluster(cl, action_vecs, inc_2025)
  
  returns    <- safe_extract(results, "total_return", default = -50)
  eliminated <- safe_extract(results, "eliminated",   default = 0)
  inc_2030s  <- safe_extract(results, "inc_2030",     default = NA)
  
  cat(sprintf("Mean return: %.3f\n", mean(returns, na.rm=TRUE)))
  cat(sprintf("Eliminated:  %d/%d\n", sum(eliminated), n_eval))
  cat(sprintf("inc_2030:    mean=%.3f  range=[%.3f, %.3f]\n",
              mean(inc_2030s, na.rm=TRUE),
              min(inc_2030s, na.rm=TRUE), max(inc_2030s, na.rm=TRUE)))
  
  rep_r <- results[[1]]
  if (!isTRUE(rep_r$error)) {
    cat("\n── Learned per-year trajectory (representative) ──\n")
    for (k in seq_along(RL_YEARS)) {
      yr <- RL_YEARS[k]
      cv <- rep_r$cov_traj[[as.character(yr)]]
      fc_C <- if (cv$Cov_np_C > 1e-12) cv$Cov_neg_C/cv$Cov_np_C else NA
      fc_P <- if (cv$Cov_np_P_monthly > 1e-12) cv$Cov_neg_P_monthly/cv$Cov_np_P_monthly else NA
      cat(sprintf(
        "  %d: Cov_np_C=%.4f Cov_neg_C=%.4f (fc=%.3f) | Cov_np_P_m=%.5f Cov_neg_P_m=%.5f (fc=%.3f) | inc=%.2f prison=%d\n",
        yr, cv$Cov_np_C, cv$Cov_neg_C, fc_C %||% NA_real_,
        cv$Cov_np_P_monthly, cv$Cov_neg_P_monthly, fc_P %||% NA_real_,
        rep_r$inc_traj[k], round(rep_r$prison_vec[k])))
    }
  }
  
  list(returns = returns, eliminated = eliminated, inc_2030 = inc_2030s,
       results = results, policy_mean = policy_mean)
}

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
try(stopCluster(master_cl), silent = TRUE)

# Build a fresh one (this exports current constants automatically)
master_cl <- setup_cluster(N_CORES)
sync_to_workers(master_cl)

# 3. Verify new caps on workers
clusterEvalQ(master_cl, c(
  COV_MAX_C_NP = COV_MAX_C_NP,
  M_THRESHOLD  = M_THRESHOLD,
  M_MAX_LOW    = M_MAX_LOW,
  M_MAX_HIGH   = M_MAX_HIGH
))

# 4. Smoke test — try aggressive multipliers, see them get capped
clusterEvalQ(master_cl, {
  prev <- get_initial_coverages()
  cat("prev Cov_np_C       =", prev$Cov_np_C, "\n")
  cat("prev Cov_np_P_m     =", prev$Cov_np_P_monthly, "\n")
  res1 <- apply_action_to_coverages(c(2.0, 1.0, 2.0, 1.0), prev)
  cat("after m=2 attempt:\n")
  cat("  Cov_np_C =", res1$Cov_np_C, " (m_used =", res1$m_np_C, ", cap was", res1$m_np_C_cap, ")\n")
  cat("  Cov_np_P_m =", res1$Cov_np_P_monthly, " (m_used =", res1$m_np_P, ", cap was", res1$m_np_P_cap, ")\n")
})[[1]]




cat("\n=== Step 1: P_PWID 2025 anchor incidence ===\n")
inc_2025 <- compute_inc_2025()
cat("inc_P_PWID_2025 =", round(inc_2025, 3), "per 100 PY\n")
cat("Trajectory:", round(make_trajectory(inc_2025), 2), "\n")

cat("\n=== Step 2: Single-episode smoke test ===\n")
set.seed(42)
test_action <- rep(1.2, 20)  # uniform 20% scale-up per year
t_test <- system.time({ test_ep <- run_episode(test_action, inc_2025) })
cat("  time:", round(t_test["elapsed"], 1), "sec\n")
cat("  return:", round(test_ep$total_return, 3), "\n")
cat("  inc trajectory:", round(test_ep$inc_traj, 2), "\n")
cat("  prison tests:", round(test_ep$prison_vec), "\n")
cat("  per-year coverages:\n")
for (k in seq_along(RL_YEARS)) {
  yr <- RL_YEARS[k]
  cv <- test_ep$cov_traj[[as.character(yr)]]
  cat(sprintf("    %d: Cov_np_C=%.4f Cov_neg_C=%.4f | Cov_np_P_m=%.5f Cov_neg_P_m=%.5f\n",
              yr, cv$Cov_np_C, cv$Cov_neg_C, cv$Cov_np_P_monthly, cv$Cov_neg_P_monthly))
}

cat("\n=== Step 3: SAC training ===\n")
master_cl <- setup_cluster(N_CORES)
on.exit(try(stopCluster(master_cl), silent = TRUE), add = TRUE)

t_train <- system.time({
  sac_result <- train_sac_gaussian(inc_2025 = inc_2025, cl = master_cl)
})
cat("\nTotal training time:", round(t_train["elapsed"]/60, 1), "min\n")

eval_result <- evaluate_policy(sac_result$policy_mean, inc_2025, cl = master_cl)
try(stopCluster(master_cl), silent = TRUE)

save(sac_result, eval_result, inc_2025,
     file = "/Users/jjwu/Projects/Simplified-HCV-testing-model/Projects/POC_prisons_scaleup/02. Output/rl_coverage_v3_results.RData")
cat("\nSaved to rl_coverage_v3_results.RData\n")






#### extract the learned parameters #### 
# Decode the policy mean into per-year coverages

decode_policy <- function(action_vec) {
  prev_cov <- get_initial_coverages()
  rows <- list()
  for (k in seq_along(RL_YEARS)) {
    yr <- RL_YEARS[k]
    a4 <- action_vec[(4*(k-1)+1):(4*k)]
    new_cov <- apply_action_to_coverages(a4, prev_cov)
    rows[[k]] <- data.frame(
      year              = yr,
      m_np_C            = new_cov$m_np_C,
      m_neg_C           = new_cov$m_neg_C,
      m_np_P            = new_cov$m_np_P,
      m_neg_P           = new_cov$m_neg_P,
      Cov_np_C          = new_cov$Cov_np_C,
      Cov_neg_C         = new_cov$Cov_neg_C,
      Cov_np_P_monthly  = new_cov$Cov_np_P_monthly,
      Cov_neg_P_monthly = new_cov$Cov_neg_P_monthly,
      fc_C              = new_cov$Cov_neg_C / new_cov$Cov_np_C,
      fc_P              = new_cov$Cov_neg_P_monthly / new_cov$Cov_np_P_monthly
    )
    prev_cov <- new_cov[c("Cov_np_C","Cov_neg_C","Cov_np_P_monthly","Cov_neg_P_monthly")]
  }
  do.call(rbind, rows)
}

# Mean policy (deterministic, what RL settled on)
mean_policy <- decode_policy(sac_result$policy_mean)
cat("── Learned mean policy ──\n")
print(mean_policy, digits = 4)

# Best single episode (highest-return rollout seen during training)
best_policy <- decode_policy(sac_result$best_action)
cat("\n── Best single episode action ──\n")
print(best_policy, digits = 4)
