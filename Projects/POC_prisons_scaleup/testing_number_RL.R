# ══════════════════════════════════════════════════════════════════════════════
# RL ENVIRONMENT 
# Action: relative change δ ∈ [-0.5, +1.0]
# Episode: year-by-year, calibrate_NP_year per step
# ══════════════════════════════════════════════════════════════════════════════
# ── Run ONCE before training — not inside each episode ───────────────────────
cat("=== Pre-computing res_2025 anchor ===\n")

n_ab_np_2025 <- list()
n_ab_np_2025[["2025"]] <- c(ANCHOR_C, ANCHOR_P)

res_2025_anchor <- calibrate_NP_year(
  cal_year       = 2025,
  prev_dfList_NP = dfList_NP_2024,
  prev_fs        = fs[["2024"]],
  prev_model     = res_2024,
  dfList_NP_base = dfList_NP_2024,
  pj             = POC_AU,
  Ccal = Ccal, fm = fm, frac_test = frac_test, NPlst = NPlst,
  n_ab_np  = n_ab_np_2025, frac_ab = frac_ab, param_var = param_var,
  best_estimates = best_estimates, best_est_pop = best_est_pop,
  disease_progress = disease_progress, pop_array = pop_array,
  dfList = dfList, fib = fib, endY = endY,
  target_C     = ANCHOR_C,
  target_P     = ANCHOR_P,
  adj_factor_C = NULL,
  adj_factor_P = NULL)

cat("2025 anchor C:", round(res_2025_anchor$diagnostics$tot_C),
    "/ ", ANCHOR_C, "\n")
cat("2025 anchor P:", round(res_2025_anchor$diagnostics$tot_P),
    "/ ", ANCHOR_P, "\n")
# ── Constants ─────────────────────────────────────────────────────────────────
RL_YEARS     <- 2026:2030
ANCHOR_C     <- 11000
ANCHOR_P     <- 14000
GAMMA        <- 0.95
ALPHA_SLOPE  <- 0.5
TERMINAL_B   <- 10.0
MIN_DELTA    <- -0.5
MAX_DELTA    <-  1.0
N_CORES      <- 6
N_EPISODES   <- 500
FM_P_MIN   <- 5.0
FM_P_MAX   <- 18.5
N_AB_P_MAX <- 75000 

apply_prison_action <- function(delta_fm_P, delta_Ccal_P,
                                prev_fm_P, prev_Ccal_P) {
  new_fm_P   <- max(FM_P_MIN, min(FM_P_MAX,
                                  prev_fm_P * (1 + delta_fm_P)))
  new_Ccal_P <- max(CCAL_P_BASE * 0.5, min(1.0,
                                           prev_Ccal_P * (1 + delta_Ccal_P)))
  
  list(fm_P   = new_fm_P,
       Ccal_P = new_Ccal_P,
       fp_vec = new_Ccal_P * new_fm_P)
}
# ── Safe extractor ────────────────────────────────────────────────────────────
safe_num <- function(x, default = 0) {
  if (is.null(x) || length(x) == 0) return(as.numeric(default))
  if (is.list(x)) x <- x[[1]]
  if (is.null(x) || length(x) == 0 || is.na(x)) return(as.numeric(default))
  as.numeric(x)
}
safe_extract <- function(results_list, field, default = 0) {
  vapply(results_list, function(r) safe_num(r[[field]], default), numeric(1))
}

# ── Action helpers ────────────────────────────────────────────────────────────
MIN_DELTA   <- -0.5
MAX_DELTA_C <-  1.0   # community: up to 2× per year
MAX_DELTA_P <-  1.0   # prison: up to 2× per year

clip_delta <- function(delta_C, delta_P) {
  list(
    delta_C = max(MIN_DELTA, min(MAX_DELTA_C, delta_C)),
    delta_P = max(MIN_DELTA, min(MAX_DELTA_P, delta_P))
  )
}

get_effective_fp <- function(fp_vec, NPlst) {
  # Param_cal clips: scVal_dt = NPlst × fp_vec, then min(scVal_dt, 1)
  # Effective fp = min(NPlst × fp_vec, 1) / NPlst = min(fp_vec, 1/NPlst)
  c(
    min(fp_vec[1], 1/NPlst[["C"]][["tau_ab"]]),
    min(fp_vec[2], 1/NPlst[["C"]][["tau_ab"]]),
    min(fp_vec[3], 1/NPlst[["P"]][["tau_ab"]]),
    min(fp_vec[4], 1/NPlst[["P"]][["tau_ab"]]),
    min(fp_vec[5], 1/NPlst[["P"]][["tau_ab"]])
  )
}

delta_to_target <- function(delta_C, delta_P, prev_C, prev_P) {
  list(
    C = max(3000,  min(30000, round(prev_C * (1 + delta_C)))),
    P = max(5000,  min(70000, round(prev_P * (1 + delta_P))))
  )
}

actions_to_vec <- function(actions) {
  unlist(lapply(actions, function(a) c(a$delta_C, a$delta_P)))
}

vec_to_actions <- function(vec) {
  lapply(seq_along(RL_YEARS), function(i)
    list(delta_C = vec[2*i - 1], delta_P = vec[2*i]))
}

random_actions <- function() {
  lapply(seq_along(RL_YEARS), function(i)
    list(delta_C = runif(1, MIN_DELTA, MAX_DELTA),
         delta_P = runif(1, MIN_DELTA, MAX_DELTA)))
}

# ── Baseline incidence ────────────────────────────────────────────────────────
compute_baseline_incidence <- function(model, pj, endY) {
  yr_idx <- 1  # year 1 = 2015
  
  inf_all <- modres.flow.t(pj, model, endYear = endY,
                           allp = "newInfections") %>%
    ungroup() %>%
    group_by(year, population) %>%
    summarise(best = sum(best), .groups = "drop") %>%
    filter(year == yr_idx)
  
  pop_n <- popResults_MidYear(
    pj, model,
    Population   = pj$popNames,
    Disease_prog = NULL, Cascade = NULL,
    param = NULL, endYear = endY) %>%
    ungroup() %>%
    filter(year == yr_idx) %>%
    group_by(population) %>%
    summarise(N = sum(best), .groups = "drop")
  
  inf_C <- inf_all %>%
    filter(population %in% c("C_PWID","C_fPWID")) %>%
    pull(best) %>% sum()
  inf_P <- inf_all %>%
    filter(population %in% c("P_PWID","P_fPWID","P_nPWID")) %>%
    pull(best) %>% sum()
  N_C <- pop_n %>%
    filter(population %in% c("C_PWID","C_fPWID")) %>%
    pull(N) %>% sum()
  N_P <- pop_n %>%
    filter(population %in% c("P_PWID","P_fPWID","P_nPWID")) %>%
    pull(N) %>% sum()
  
  list(inc_C_2015 = 1000 * inf_C / N_C,
       inc_P_2015 = 1000 * inf_P / (2 * N_P))
}

# ── Incidence from model at a given year ──────────────────────────────────────
compute_incidence_year <- function(model, yr_idx, pj, endY) {
  
  inf_all <- modres.flow.t(pj, model, endYear = endY,
                           allp = "newInfections") %>%
    ungroup() %>%
    group_by(year, population) %>%
    summarise(best = sum(best), .groups = "drop") %>%
    filter(year == yr_idx)
  
  pop_n <- popResults_MidYear(
    pj, model,
    Population   = pj$popNames,
    Disease_prog = NULL, Cascade = NULL,
    param = NULL, endYear = endY) %>%
    ungroup() %>%
    filter(year == yr_idx) %>%
    group_by(population) %>%
    summarise(N = sum(best), .groups = "drop")
  
  inf_C <- inf_all %>%
    filter(population %in% c("C_PWID","C_fPWID")) %>%
    pull(best) %>% sum()
  inf_P <- inf_all %>%
    filter(population %in% c("P_PWID","P_fPWID","P_nPWID")) %>%
    pull(best) %>% sum()
  N_C <- pop_n %>%
    filter(population %in% c("C_PWID","C_fPWID")) %>%
    pull(N) %>% sum()
  N_P <- pop_n %>%
    filter(population %in% c("P_PWID","P_fPWID","P_nPWID")) %>%
    pull(N) %>% sum()
  
  list(inc_C = 1000 * inf_C / N_C,
       inc_P = 1000 * inf_P / (2 * N_P),
       N_C = N_C, N_P = N_P)
}

# ── State vector ──────────────────────────────────────────────────────────────
compute_state <- function(model, cal_year, prev_n_C, prev_n_P, pj, endY) {
  
  yr_idx <- cal_year - pj$cabY
  inc    <- compute_incidence_year(model, yr_idx, pj, endY)
  
  undiag <- popResults_MidYear(
    pj, model,
    Population   = pj$popNames,
    Disease_prog = pj$diseaseprogress_Name,
    Cascade      = pj$cascade_name,
    param = NULL, endYear = endY) %>%
    ungroup() %>%
    mutate(cascade      = as.character(cascade),
           disease_prog = as.character(disease_prog)) %>%
    filter(year == yr_idx,
           cascade == "undiag",
           disease_prog != "a") %>%
    group_by(population) %>%
    summarise(U = sum(best), .groups = "drop")
  
  U_C <- undiag %>%
    filter(population %in% c("C_PWID","C_fPWID")) %>%
    pull(U) %>% sum()
  U_P <- undiag %>%
    filter(population %in% c("P_PWID","P_fPWID","P_nPWID")) %>%
    pull(U) %>% sum()
  
  c(inc_C       = inc$inc_C,
    inc_P       = inc$inc_P,
    N_C         = inc$N_C,
    N_P         = inc$N_P,
    undiag_C    = U_C / max(inc$N_C, 1),
    undiag_P    = U_P / max(inc$N_P, 1),
    prev_n_C    = prev_n_C,
    prev_n_P    = prev_n_P,
    t_remaining = 2030 - cal_year)
}

# ── Reward ────────────────────────────────────────────────────────────────────
compute_reward <- function(s_t, s_prev, cal_year, baseline,
                           alpha = ALPHA_SLOPE, B = TERMINAL_B) {
  
  inc_C_target <- 0.20 * baseline$inc_C_2015
  inc_P_target <- 0.20 * baseline$inc_P_2015
  
  gap_C   <- max(0, s_t["inc_C"] - inc_C_target) / baseline$inc_C_2015
  gap_P   <- max(0, s_t["inc_P"] - inc_P_target) / baseline$inc_P_2015
  delta_C <- (s_t["inc_C"] - s_prev["inc_C"]) / baseline$inc_C_2015
  delta_P <- (s_t["inc_P"] - s_prev["inc_P"]) / baseline$inc_P_2015
  
  # Slope reward only when above target
  slope_r_C <- ifelse(s_t["inc_C"] > inc_C_target, -alpha * delta_C, 0)
  slope_r_P <- ifelse(s_t["inc_P"] > inc_P_target, -alpha * delta_P, 0)
  r_step    <- -gap_C - gap_P + slope_r_C + slope_r_P
  
  r_terminal <- 0; elim_C <- NA; elim_P <- NA
  if (cal_year == 2030) {
    elim_C     <- as.numeric(s_t["inc_C"] <= inc_C_target)
    elim_P     <- as.numeric(s_t["inc_P"] <= inc_P_target)
    r_terminal <- B*elim_C + B*elim_P - B*(1-elim_C) - B*(1-elim_P)
  }
  
  list(r_total = r_step + r_terminal, r_step = r_step,
       r_terminal = r_terminal, gap_C = gap_C, gap_P = gap_P,
       delta_C = delta_C, delta_P = delta_P,
       elim_C = elim_C, elim_P = elim_P)
}

# ── Single episode ────────────────────────────────────────────────────────────
MIN_DELTA   <- -0.5
MAX_DELTA_C <-  1.0
MAX_DELTA_P <-  1.0

clip_delta <- function(delta_C, delta_P) {
  list(
    delta_C = max(MIN_DELTA, min(MAX_DELTA_C, delta_C)),
    delta_P = max(MIN_DELTA, min(MAX_DELTA_P, delta_P))
  )
}

delta_to_target <- function(delta_C, delta_P, prev_C, prev_P) {
  list(
    C = max(3000,  min(40000, round(prev_C * (1 + delta_C)))),
    P = max(5000,  min(60000, round(prev_P * (1 + delta_P))))
  )
}

# ── run_episode — direct mode ─────────────────────────────────────────────────
run_episode <- function(actions, baseline) {
  
  pj          <- POC_AU
  endY_val    <- endY
  prev_model  <- res_2025_anchor$model
  prev_dfList <- res_2025_anchor$dfList_NP
  prev_fs     <- res_2025_anchor$fs
  prev_n_C    <- ANCHOR_C
  prev_fm_P   <- FM_P_BASE
  prev_Ccal_P <- CCAL_P_BASE
  total_return <- 0
  trajectory   <- list()
  
  s_prev <- compute_state(res_2025_anchor$model, 2025,
                          ANCHOR_C, ANCHOR_C, pj, endY_val)
  
  for (step in seq_along(RL_YEARS)) {
    yr  <- RL_YEARS[step]
    yrc <- as.character(yr)
    
    a <- clip_action(actions[[step]]$delta_C,
                     actions[[step]]$delta_fm_P,
                     actions[[step]]$delta_Ccal_P)
    
    # ── Community: delta → target_C ──────────────────────────────────────────
    new_C <- delta_to_C(a$delta_C, prev_n_C)
    
    # ── Prison: delta → fm_P and Ccal_P ──────────────────────────────────────
    pa <- apply_prison_action(a$delta_fm_P, a$delta_Ccal_P,
                              prev_fm_P, prev_Ccal_P)
    
    res_yr <- tryCatch(
      calibrate_NP_year(
        cal_year       = yr,
        prev_dfList_NP = prev_dfList,
        prev_fs        = prev_fs,
        prev_model     = prev_model,
        dfList_NP_base = prev_dfList,
        pj=pj, Ccal=Ccal, fm=fm, frac_test=frac_test, NPlst=NPlst,
        n_ab_np=NULL, frac_ab=frac_ab, param_var=param_var,
        best_estimates=best_estimates, best_est_pop=best_est_pop,
        disease_progress=disease_progress, pop_array=pop_array,
        dfList=dfList, fib=fib, endY=endY_val,
        mode          = "fm_ccal",
        target_C      = new_C,
        target_fm_P   = pa$fm_P,
        target_Ccal_P = pa$Ccal_P),
      error = function(e) {
        cat("  [Error yr", yr, "]:", conditionMessage(e), "\n"); NULL
      })
    
    if (is.null(res_yr)) {
      total_return <- total_return - 20; break
    }
    
    s_t  <- compute_state(res_yr$model, yr, new_C,
                          pa$Ccal_P, pj, endY_val)
    r    <- compute_reward(s_t, s_prev, yr, baseline)
    disc <- GAMMA^(step - 1)
    total_return <- total_return + disc * r$r_total
    
    trajectory[[yrc]] <- list(
      year      = yr,
      target_C  = new_C,
      fm_P      = pa$fm_P,
      Ccal_P    = pa$Ccal_P,
      fp_vec_P  = pa$fp_vec,
      n_ab_P    = res_yr$diagnostics$n_ab_P,
      tot_C     = res_yr$diagnostics$tot_C,
      tot_P     = res_yr$diagnostics$tot_P,
      delta_C   = a$delta_C,
      delta_fm_P= a$delta_fm_P,
      delta_Ccal_P = a$delta_Ccal_P,
      state     = s_t,
      reward    = r)
    
    s_prev      <- s_t
    prev_n_C    <- new_C
    prev_fm_P   <- pa$fm_P
    prev_Ccal_P <- pa$Ccal_P
    prev_model  <- res_yr$model
    prev_dfList <- res_yr$dfList_NP
    prev_fs     <- res_yr$fs
  }
  
  list(total_return=total_return, trajectory=trajectory,
       elim_C=safe_num(trajectory[["2030"]]$reward$elim_C),
       elim_P=safe_num(trajectory[["2030"]]$reward$elim_P),
       actions=actions)
}







# ── Parallel runner ───────────────────────────────────────────────────────────
run_episodes_parallel <- function(all_actions_list, n_cores = N_CORES,
                                  baseline) {
  cl <- makeCluster(n_cores)
  
  all_fns  <- Filter(function(x) is.function(get(x, envir=.GlobalEnv)),
                     ls(envir=.GlobalEnv))
  all_objs <- Filter(function(x) !is.function(get(x, envir=.GlobalEnv)),
                     ls(envir=.GlobalEnv))
  clusterExport(cl, varlist = c(all_fns, all_objs), envir = .GlobalEnv)
  
  local_env <- environment()
  clusterExport(cl, varlist = "baseline", envir = local_env)
  
  clusterEvalQ(cl, {
    library(dplyr); library(tidyr)
    library(Rcpp); library(RcppArmadillo)
    sourceCpp("~/Projects/Simplified-HCV-testing-model/Projects/POC_prisons_scaleup/hcvmsm_core.cpp")
    source("~/Projects/Simplified-HCV-testing-model/Projects/POC_prisons_scaleup/HCVMSM_cpp.R")
  })
  
  check <- clusterEvalQ(cl, {
    needed  <- c("calibrate_NP_year","HCVMSM_cpp","hcvmsm_loop_cpp",
                 "modres.flow.t","popResults_MidYear","Param_cal",
                 "compute_state","compute_reward","run_episode", "res_2025_anchor")
    missing <- needed[!sapply(needed, exists)]
    if (length(missing) > 0) paste("MISSING:", paste(missing, collapse=","))
    else "OK"
  })
  missing_report <- unique(unlist(check))
  if (any(missing_report != "OK")) {
    stopCluster(cl)
    stop("Workers missing: ", paste(missing_report[missing_report!="OK"],
                                    collapse="\n"))
  }
  cat("Worker check: all functions found ✓\n")
  
  results <- parLapply(cl, all_actions_list, function(actions)
    run_episode(actions, baseline = baseline))
  
  stopCluster(cl)
  results
}

# ── SAC training ──────────────────────────────────────────────────────────────
train_sac_gaussian <- function(n_episodes   = N_EPISODES,
                               n_cores      = N_CORES,
                               baseline,
                               batch_size   = 20,
                               lr           = 0.05,
                               entropy_coef = 0.1) {
  
  cat("=== SAC Training ===\n")
  cat("Target C ≤", round(0.20*baseline$inc_C_2015,3),
      "| P ≤", round(0.20*baseline$inc_P_2015,3), "\n\n")
  
  policy_mean    <- rep(0.0, 10)
  policy_log_std <- rep(log(0.3), 10)
  best_return    <- -Inf
  best_actions   <- NULL
  history        <- data.frame()
  n_batches      <- ceiling(n_episodes / batch_size)
  
  for (batch in 1:n_batches) {
    cat("\n── Batch", batch, "/", n_batches, "──\n")
    
    policy_std  <- exp(policy_log_std)
    action_vecs <- lapply(1:batch_size, function(i) {
      vec <- rnorm(10, mean=policy_mean, sd=policy_std)
      pmax(MIN_DELTA, pmin(MAX_DELTA, vec))
    })
    all_actions <- lapply(action_vecs, vec_to_actions)
    
    t_batch <- system.time({
      results <- run_episodes_parallel(all_actions, n_cores, baseline)
    })
    
    returns <- safe_extract(results, "total_return", default=-20)
    elim_C  <- safe_extract(results, "elim_C",       default=0)
    elim_P  <- safe_extract(results, "elim_P",       default=0)
    
    cat("Returns: mean=", round(mean(returns),3),
        "| max=",  round(max(returns),3),
        "| elim_C:", sum(elim_C), "/", batch_size,
        "| elim_P:", sum(elim_P), "/", batch_size,
        "| time:", round(t_batch["elapsed"],1), "sec\n")
    
    # ── Policy gradient update ────────────────────────────────────────────────
    advantages     <- returns - mean(returns)
    grad_mean      <- rep(0, 10)
    grad_log_std   <- rep(0, 10)
    
    for (i in seq_along(results)) {
      vec <- action_vecs[[i]]; adv <- advantages[i]
      g_mean    <- (vec - policy_mean) / (policy_std^2)
      g_log_std <- ((vec - policy_mean)^2 / policy_std^2) - 1
      grad_mean    <- grad_mean    + adv * g_mean
      grad_log_std <- grad_log_std + adv * g_log_std + entropy_coef
    }
    grad_mean    <- grad_mean    / batch_size
    grad_log_std <- grad_log_std / batch_size
    
    policy_mean    <- pmax(MIN_DELTA, pmin(MAX_DELTA,
                                           policy_mean + lr * grad_mean))
    policy_log_std <- pmax(-3, pmin(0,
                                    policy_log_std + lr * grad_log_std))
    
    best_idx <- which.max(returns)
    if (returns[best_idx] > best_return) {
      best_return  <- returns[best_idx]
      best_actions <- all_actions[[best_idx]]
      cat("New best return:", round(best_return,3), "\n")
    }
    
    history <- rbind(history, data.frame(
      episode = (batch-1)*batch_size + seq_along(results),
      batch=batch, ret=returns, elim_C=elim_C, elim_P=elim_P))
    
    # ── Print policy ──────────────────────────────────────────────────────────
    cat("Policy [δC, δP] → [target_C, target_P] per year:\n")
    n_C <- ANCHOR_C; n_P <- ANCHOR_P
    for (i in seq_along(RL_YEARS)) {
      dC <- policy_mean[2*i-1]; dP <- policy_mean[2*i]
      tg <- delta_to_target(dC, dP, n_C, n_P)
      cat("  ", RL_YEARS[i], ": δC=", round(dC,3), "δP=", round(dP,3),
          "→ C=", tg$C, "P=", tg$P, "\n")
      n_C <- tg$C; n_P <- tg$P
    }
  }
  
  list(policy_mean=policy_mean, policy_log_std=policy_log_std,
       best_return=best_return, best_actions=best_actions,
       history=history)
}

# ── Evaluate ──────────────────────────────────────────────────────────────────
evaluate_policy <- function(policy_mean, baseline, n_eval=10) {
  cat("\n=== Evaluating policy ===\n")
  best_actions <- vec_to_actions(pmax(MIN_DELTA, pmin(MAX_DELTA, policy_mean)))
  results <- run_episodes_parallel(rep(list(best_actions), n_eval),
                                   n_cores=min(n_eval, N_CORES),
                                   baseline=baseline)
  returns <- safe_extract(results, "total_return", default=-20)
  elim_C  <- safe_extract(results, "elim_C", default=0)
  elim_P  <- safe_extract(results, "elim_P", default=0)
  
  cat("Mean return:", round(mean(returns),3), "\n")
  cat("Elim C:", sum(elim_C), "/", n_eval, "\n")
  cat("Elim P:", sum(elim_P), "/", n_eval, "\n")
  cat("Both:  ", sum(elim_C & elim_P), "/", n_eval, "\n")
  
  cat("\nOptimal target tests per year:\n")
  n_C <- ANCHOR_C; n_P <- ANCHOR_P
  for (i in seq_along(RL_YEARS)) {
    tg <- delta_to_target(policy_mean[2*i-1], policy_mean[2*i], n_C, n_P)
    cat("  ", RL_YEARS[i], ": target_C =", tg$C, "| target_P =", tg$P, "\n")
    n_C <- tg$C; n_P <- tg$P
  }
  
  list(returns=returns, elim_C=elim_C, elim_P=elim_P,
       best_actions=best_actions)
}

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
cat("=== Step 1: Baseline incidence ===\n")
baseline <- compute_baseline_incidence(res_2025, POC_AU, endY)
cat("inc_C_2015:", round(baseline$inc_C_2015,4), "per 1000 PY\n")
cat("inc_P_2015:", round(baseline$inc_P_2015,4), "per 1000 PY\n")
cat("Target C ≤", round(0.20*baseline$inc_C_2015,4), "\n")
cat("Target P ≤", round(0.20*baseline$inc_P_2015,4), "\n")

cat("\n=== Step 2: Single episode test ===\n")
t_test <- system.time({
  test_ep <- run_episode(random_actions(), baseline)
})
cat("Time:", round(t_test["elapsed"],1), "sec\n")
cat("Return:", round(test_ep$total_return,3), "\n")
cat("Trajectory years:", names(test_ep$trajectory), "\n")
cat("Elim C:", test_ep$elim_C, "| Elim P:", test_ep$elim_P, "\n")

# ── Start SAC training ────────────────────────────────────────────────────────
cat("=== Starting SAC training ===\n")
cat("Episodes:", N_EPISODES, "| Cores:", N_CORES, "| Batch:", 20, "\n")
cat("Estimated time: ~", round(N_EPISODES * 136 / N_CORES / 60, 0), "min\n\n")

t_train <- system.time({
  sac_result <- train_sac_gaussian(
    n_episodes   = N_EPISODES,
    n_cores      = N_CORES,
    baseline     = baseline,
    batch_size   = 20,
    lr           = 0.1,
    entropy_coef = 0.3)
})
cat("\nTotal training time:", round(t_train["elapsed"]/60, 1), "min\n")

# ── Evaluate ──────────────────────────────────────────────────────────────────
eval_result <- evaluate_policy(sac_result$policy_mean, baseline)

# ── Save everything ───────────────────────────────────────────────────────────
save(sac_result, eval_result, baseline, res_2025_anchor,
     file = "rl_sac_results.RData")
cat("Saved to rl_sac_results.RData\n")



