# =============================================================================
# calibration_pipeline.R
# =============================================================================
# Calibrates NP coverage parameters across cal_years through iterative
# stock-based passes. Pass 1 chains across years; passes 2..max_passes use
# the post-NP scenario from the previous pass as a universal stock source
# until convergence (or max_passes is reached).
# =============================================================================
# Requires in workspace before sourcing:
#   POC_AU, NPlst, dfList, best_estimates, best_est_pop, disease_progress,
#   pop_array, fib, fm, param_var, endY, max_passes, tolerance,
#   Param_cal(), scale_CT_eta(), scale_CT_RNA(), scale_CT_ab(),
#   carry_forward_prison(), HCV_test()
# =============================================================================

library(dplyr); library(tidyr); library(readxl)
load("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/POC_AU/02. Output/POC_AU.rda")
load("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/POC_AU/02. Output/POC_AUcali.rda")
load("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/POC_AU/02. Output/POC_AUcali_timev.rda")


project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Documents/Simplified-HCV-testing-model")

data_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                   project_name, sep = "")
# =============================================================================
# 1. INPUTS
# =============================================================================

obs_data <- tribble(
  ~year, ~setting,    ~T_ab, ~T_rna, ~M,    ~Tr,
  2022,  "Community",   926,  1217, 1701,    97,
  2022,  "Prison",        0,  4524, 4082,   543,
  2023,  "Community",  2407,  3049, 4613,   233,
  2023,  "Prison",     1030,  4990, 5104,   588,
  2024,  "Community",  5254,  4241, 7881,   406,
  2024,  "Prison",     4697,  6201, 9600,   754
)
Num_test_person_NP <- read_excel(paste0(data_path, "/01. DATA/Num_test_person_NP.xlsx"))%>%
  select(year, settings, num_ab, num_RNA, num_tests, num_person)%>%
  mutate(num_RNA = as.numeric(num_RNA), 
         num_tests = as.numeric(num_tests),
         num_person = as.numeric(num_person))%>%
  filter(year%in% c(2022:2024))

Num_test_person_NP <- Num_test_person_NP%>%
  mutate(reflex_frac = num_ab/num_person, 
         immedRNA_frac = 1 - reflex_frac) 


# current achievement data 
current_data <- function(dt, y, s, index){ 
  x <- dt%>%filter(year == y & settings == s)
  x <- x[[index]]
  return(x)
}

# extracted the coverage of national program and fraction of testing pathways 

frac_test <- list() 
reflex_frac_C <- list()
reflex_frac_P <- list()


for(i in c(2022:2024)){ 
  
  reflex_frac_C[[i]] <- current_data(Num_test_person_NP, y = i, 
                                     s = "community", index = "reflex_frac")
  
  reflex_frac_P[[i]] <- current_data(Num_test_person_NP, y = i, 
                                     s = "prison", index = "reflex_frac")
  
  frac_test[[i]] <- list("C" = list("reflex" = as.numeric(reflex_frac_C[[i]]), 
                                    "immeRNA" = 1- as.numeric(reflex_frac_C[[i]])),
                         "P" = list("reflex" = as.numeric(reflex_frac_P[[i]]), 
                                    "immeRNA" = 1- as.numeric(reflex_frac_P[[i]])))
   
}
# for year of 2022, using RNA testing in prisons for n_ab to estimate the testing coverage 
n_ab <- Num_test_person_NP%>%select(year, settings, num_ab)%>%
  spread(settings, -c(year))

n_ab[1, "prison"] <- Num_test_person_NP%>%
  filter(year == 2022 & settings == "prison")%>%select(num_RNA)%>%unlist()%>%c()

# efficacy of national program 

np_effect <- read_excel(paste0(data_path, "/01. DATA/effect_np.xlsx"))

NP_tauab_C <- unlist(as.numeric(np_effect[1,2]))

NP_tauRNA_C <- unlist(as.numeric(np_effect[2,2]))

NP_tauRNAonly_C <- unlist(as.numeric(np_effect[3,2]))

# the treatment initiation is the % of people tested RNA+ initiated DAA within 120 days 
# old 
# NP_eta_C <-  1- (1- unlist(as.numeric(np_effect[4,2])))^(1/POC_AU$timestep/4)

# test on new data informed by national program [2025/07/30]
# NP_eta_C <- 0.6
NP_eta_C <- 1- (1- unlist(as.numeric(np_effect[4,2])))

NP_tauab_P <- unlist(as.numeric(np_effect[1,3]))

NP_tauRNA_P <- unlist(as.numeric(np_effect[2,3]))


NP_tauRNAonly_P <- unlist(as.numeric(np_effect[3,3]))

# NP_eta_P <- 0.89 for 1 month 
# turn it back to annual probability
NP_eta_P <- 1- (1-unlist(as.numeric(np_effect[4,3])))^(1/POC_AU$timestep/1)
NPlst <- list("C" = list("tau_ab" = NP_tauab_C,
                         "tau_RNA" = NP_tauRNA_C, 
                         "tau_poct" = NP_tauRNAonly_C,
                         "eta" = NP_eta_C),
              "P" = list("tau_ab" = NP_tauab_P,
                         "tau_RNA" = NP_tauRNA_P, 
                         "tau_poct" = NP_tauRNAonly_P,
                         "eta" = NP_eta_P))
NP_eta_C <- NPlst$C$eta
NP_eta_P <- NPlst$P$eta

cat("Inputs:\n")
cat(sprintf("  NP_eta_C = %.3f\n", NP_eta_C))
cat(sprintf("  NP_eta_P = %.15f\n", NP_eta_P))


# =============================================================================
# 2. STATE-SET CONSTANTS
# =============================================================================

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


# =============================================================================
# 3. TIME HELPERS
# =============================================================================

year_time_index <- function(pj, year) {
  steps_per_year <- 1 / pj$timestep
  k  <- year - pj$cabY + 1
  t0 <- (k - 1) * steps_per_year + 1
  t1 <-  k      * steps_per_year
  list(t0 = t0, t1 = t1, steps_per_year = steps_per_year)
}

yr_sum_stock <- function(Sce, pj, pop_idx, states, year) {
  tt <- year_time_index(pj, year)
  sum(vapply(
    states,
    function(s) sum(Sce$allPops[pop_idx, s, tt$t0:tt$t1], na.rm = TRUE),
    numeric(1)
  ))
}

yr_avg_stock <- function(Sce, pj, pop_idx, states, year) {
  tt <- year_time_index(pj, year)
  yr_sum_stock(Sce, pj, pop_idx, states, year) / tt$steps_per_year
}

yr_sum_event <- function(mat, pj, pops, year) {
  tt <- year_time_index(pj, year)
  sum(mat[pops, tt$t0:tt$t1], na.rm = TRUE)
}

monthly_stock_sum <- function(Sce, pj, pop_idx, states, year) {
  tt <- year_time_index(pj, year)
  sapply(tt$t0:tt$t1, function(t) {
    sum(vapply(
      states,
      function(st) sum(Sce$allPops[pop_idx, st, t], na.rm = TRUE),
      numeric(1)
    ))
  })
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


# =============================================================================
# 4. STOCK EXTRACTION
# =============================================================================

get_subpop_stocks <- function(Sce, pj, year) {
  bind_rows(lapply(seq_len(5), function(p) {
    tibble(
      year     = year,
      pop      = pj$popNames[p],
      undiag   = yr_avg_stock(Sce, pj, p, undiag_states,  year),
      diag_ab  = yr_avg_stock(Sce, pj, p, diag_ab_states, year),
      diag_RNA = yr_avg_stock(Sce, pj, p, diag_RNA_states, year),
      neg      = yr_avg_stock(Sce, pj, p, neg_states,     year)
    )
  }))
}

add_setting <- function(dat) {
  dat |>
    mutate(setting = ifelse(pop %in% c("C_PWID", "C_fPWID"), "Community", "Prison"))
}

get_setting_stocks <- function(subpop_stocks) {
  subpop_stocks |>
    add_setting() |>
    group_by(year, setting) |>
    summarise(
      undiag   = sum(undiag,   na.rm = TRUE),
      diag_ab  = sum(diag_ab,  na.rm = TRUE),
      diag_RNA = sum(diag_RNA, na.rm = TRUE),
      neg      = sum(neg,      na.rm = TRUE),
      .groups  = "drop"
    )
}

# Replaces 6 near-identical fallback blocks: if "X_avg" exists and "X" doesn't,
# rename it to "X". Applied to setting_stocks and subpop_stocks before use.
rename_avg_cols <- function(df) {
  for (col in c("undiag", "diag_RNA", "neg")) {
    avg_col <- paste0(col, "_avg")
    if (!col %in% names(df) && avg_col %in% names(df)) {
      df <- df |> dplyr::mutate("{col}" := .data[[avg_col]])
    }
  }
  df
}


# =============================================================================
# 5. PARAMETER BUILDERS
# =============================================================================

build_fp <- function(yr_str, Ccal_list, fm) {
  c(
    Ccal_list[[yr_str]]$C_PWID  * fm[[yr_str]][1],
    Ccal_list[[yr_str]]$C_fPWID * fm[[yr_str]][2],
    Ccal_list[[yr_str]]$P_PWID  * fm[[yr_str]][3],
    Ccal_list[[yr_str]]$P_fPWID * fm[[yr_str]][4],
    Ccal_list[[yr_str]]$P_nPWID * fm[[yr_str]][5]
  )
}

make_fc_matrix_for_year <- function(pj, year, subpop_table, ntime) {
  tt <- year_time_index(pj, year)
  fc_mat <- matrix(0, nrow = pj$npops, ncol = ntime)
  for (p in seq_len(5)) {
    pop_name <- pj$popNames[p]
    fc_val   <- subpop_table$fc[subpop_table$pop == pop_name]
    fc_mat[p, tt$t0:tt$t1] <- fc_val
  }
  fc_mat
}

build_frac_test <- function(obs_data, years) {
  frac_test <- list()
  for (yr in years) {
    o_C <- obs_data |> filter(year == yr, setting == "Community")
    o_P <- obs_data |> filter(year == yr, setting == "Prison")
    
    reflex_C <- o_C$T_ab / max(o_C$M, 1)
    reflex_P <- o_P$T_ab / max(o_P$M, 1)
    
    frac_test[[as.character(yr)]] <- list(
      C = list(reflex = reflex_C, immeRNA = 1 - reflex_C),
      P = list(reflex = reflex_P, immeRNA = 1 - reflex_P)
    )
  }
  frac_test
}

cal_years <- c(2022, 2023, 2024)
frac_test <- build_frac_test(obs_data, cal_years)


# =============================================================================
# 6. DIAGNOSTICS HELPERS
# =============================================================================

check_param <- function(dlist, index, year, pj) {
  end_dt <- ((year + 1) - pj$cabY) / pj$timestep
  vals <- round(dlist[[index]][, 3, end_dt], 6)
  names(vals) <- pj$popNames
  cat(sprintf("\n-- %s end of %d --\n", index, year))
  print(vals)
}

print_calibration_summary <- function(cal) {
  cat("\nSetting-level stocks:\n")
  print(cal$setting_stocks)
  
  cat("\nSetting-level coverages:\n")
  print(cal$ccal_setting |>
          select(year, setting, Cov_np, M_pos, M_neg, Cov_neg, fc))
  
  cat("\nVerification:\n")
  print(cal$ccal_setting |>
          select(year, setting, M, M_recovered, Tr, Tr_check))
  
  cat("\nSanity flags:\n")
  print(cal$ccal_setting |>
          mutate(
            flag_Mpos = case_when(
              M_pos > M ~ "M_pos > M_total: undiag too large or diag_RNA too small",
              M_pos < 0 ~ "negative",
              TRUE ~ "ok"
            ),
            flag_Mneg = case_when(
              M_neg <= 0 ~ "M_neg <= 0: system cannot fit M",
              TRUE ~ "ok"
            ),
            flag_Cov = case_when(
              Cov_np > 1 ~ "Cov_np > 1",
              Cov_np < 1e-4 ~ "Cov_np very low",
              TRUE ~ "ok"
            ),
            flag_fc = case_when(
              fc > 1 ~ "fc > 1: neg-arm coverage > pos-arm coverage",
              fc < 0 ~ "fc < 0",
              TRUE ~ "ok"
            )
          ) |>
          select(year, setting, M_pos, flag_Mpos, M_neg, flag_Mneg, Cov_np, flag_Cov, fc, flag_fc))
}

print_model_vs_obs <- function(Sce, pj, obs_data, years) {
  comm <- 1:2
  pris <- 3:5
  
  cat(sprintf(
    "%-5s  %-10s  %-7s  %-7s  %-5s   %-7s  %-7s  %-5s   %-6s  %-6s  %-5s\n",
    "Year", "Setting", "Ab_mod", "Ab_obs", "%", "RNA_mod", "RNA_obs", "%", "Tr_mod", "Tr_obs", "%"
  ))
  cat(strrep("-", 95), "\n")
  
  for (yr in years) {
    for (sett in c("Community", "Prison")) {
      pops <- if (sett == "Community") comm else pris
      
      ab_mod <- round(
        yr_sum_event(Sce$newTestingAb_sc,     pj, pops, yr) +
          yr_sum_event(Sce$newTestingAb_sc_neg, pj, pops, yr)
      )
      rna_mod <- round(
        yr_sum_event(Sce$newTestingAg_sc,       pj, pops, yr) +
          yr_sum_event(Sce$newTestingAg_sc_neg, pj, pops, yr) +
          yr_sum_event(Sce$newTestingPOCT_sc,     pj, pops, yr) +
          yr_sum_event(Sce$newTestingPOCT_sc_neg, pj, pops, yr)
      )
      tr_mod <- round(yr_sum_event(Sce$newTreatment_sc, pj, pops, yr))
      
      o <- obs_data |> filter(year == yr, setting == sett)
      
      cat(sprintf(
        "%-5d  %-10s  %-7d  %-7d  %+4.0f%%  %-7d  %-7d  %+4.0f%%  %-6d  %-6d  %+4.0f%%\n",
        yr, sett,
        ab_mod,  o$T_ab,  100 * (ab_mod  - o$T_ab)  / max(o$T_ab,  1),
        rna_mod, o$T_rna, 100 * (rna_mod - o$T_rna) / max(o$T_rna, 1),
        tr_mod,  o$Tr,    100 * (tr_mod  - o$Tr)    / max(o$Tr,    1)
      ))
    }
  }
}

audit_common_prison_negative_coverage <- function(cal, year) {
  cal$ccal_setting |>
    dplyr::filter(year == !!year, setting == "Prison") |>
    dplyr::select(
      year, setting,
      M, Tr,
      M_pos, M_residual,
      M_neg_PWID_fPWID, M_x_P_nPWID, M_neg, M_recovered,
      Cov_np_monthly, Cov_neg_monthly, Cov_x_monthly,
      Cov_np, Cov_neg, Cov_x,
      fc, Tr_check
    )
}

audit_subpop_pass_through <- function(cal, year) {
  cal$subpop_table |>
    dplyr::filter(year == !!year) |>
    dplyr::select(
      year, setting, pop,
      Cov_np, Cov_neg, fc,
      Tr_predicted, M_pos_predicted, M_neg_predicted
    )
}


# =============================================================================
# 7. CALIBRATION FUNCTION
# =============================================================================

safe_divide <- function(num, den, eps = 1e-12) {
  num / pmax(den, eps)
}

calibrate_np_coverage <- function(Sce, pj, year, obs_data, NP_eta_C, NP_eta_P) {
  
  subpop_stocks  <- get_subpop_stocks(Sce, pj, year)
  setting_stocks <- get_setting_stocks(subpop_stocks)
  
  steps_per_year <- 1 / pj$timestep
  
  m2a <- function(x) {
    x <- pmin(pmax(x, 0), 1)
    1 - (1 - x)^steps_per_year
  }
  
  # Robust column names (handles legacy *_avg names)
  setting_stocks <- rename_avg_cols(setting_stocks)
  subpop_stocks  <- rename_avg_cols(subpop_stocks)
  
  # Population indices
  prison_pwid_pops <- which(pj$popNames %in% c("P_PWID", "P_fPWID"))
  prison_npwid_pop <- which(pj$popNames == "P_nPWID")
  
  if (length(prison_pwid_pops) != 2) {
    stop("Could not identify P_PWID and P_fPWID in pj$popNames.")
  }
  if (length(prison_npwid_pop) != 1) {
    stop("Could not identify P_nPWID in pj$popNames.")
  }
  
  # Prison monthly stocks
  prison_pwid_diag_RNA_monthly <- monthly_setting_stock(
    Sce = Sce, pj = pj, pops = prison_pwid_pops,
    states = diag_RNA_states, year = year
  )
  prison_pwid_undiag_monthly <- monthly_setting_stock(
    Sce = Sce, pj = pj, pops = prison_pwid_pops,
    states = undiag_states, year = year
  )
  prison_pwid_neg_monthly <- monthly_setting_stock(
    Sce = Sce, pj = pj, pops = prison_pwid_pops,
    states = neg_states, year = year
  )
  prison_npwid_neg_monthly <- monthly_setting_stock(
    Sce = Sce, pj = pj, pops = prison_npwid_pop,
    states = neg_states, year = year
  )
  
  prison_pwid_diag_RNA_year_sum <- sum(prison_pwid_diag_RNA_monthly, na.rm = TRUE)
  prison_pwid_undiag_year_sum   <- sum(prison_pwid_undiag_monthly,   na.rm = TRUE)
  prison_pwid_neg_year_sum      <- sum(prison_pwid_neg_monthly,      na.rm = TRUE)
  prison_npwid_neg_year_sum     <- sum(prison_npwid_neg_monthly,     na.rm = TRUE)
  prison_total_neg_year_sum     <- prison_pwid_neg_year_sum + prison_npwid_neg_year_sum
  
  # Debug
  obs_prison_debug <- obs_data |>
    dplyr::filter(year == !!year, setting == "Prison")
  
  cat(sprintf("\n[DEBUG %d prison denominator]\n", year))
  cat(sprintf("  P_PWID + P_fPWID diag_RNA_year_sum = %.15f\n",
              prison_pwid_diag_RNA_year_sum))
  cat(sprintf("  P_PWID + P_fPWID undiag_year_sum   = %.15f\n",
              prison_pwid_undiag_year_sum))
  cat(sprintf("  P_PWID + P_fPWID neg_year_sum      = %.15f\n",
              prison_pwid_neg_year_sum))
  cat(sprintf("  P_nPWID neg_year_sum               = %.15f\n",
              prison_npwid_neg_year_sum))
  cat(sprintf("  Prison total neg_year_sum          = %.15f\n",
              prison_total_neg_year_sum))
  cat(sprintf("  Tr_obs = %.15f\n", obs_prison_debug$Tr))
  cat(sprintf("  Cov_np_monthly expected = %.15f\n",
              safe_divide(obs_prison_debug$Tr,
                          NP_eta_P * prison_pwid_diag_RNA_year_sum)))
  
  # Setting-level calibration
  ccal_setting <- obs_data |>
    dplyr::filter(year == !!year) |>
    dplyr::left_join(setting_stocks, by = c("year", "setting")) |>
    dplyr::mutate(
      eta_NP = ifelse(setting == "Community", NP_eta_C, NP_eta_P),
      
      # Positive coverage
      Cov_np_monthly = dplyr::case_when(
        setting == "Community" ~ NA_real_,
        setting == "Prison"    ~ safe_divide(Tr, NP_eta_P * prison_pwid_diag_RNA_year_sum),
        TRUE ~ NA_real_
      ),
      Cov_np = dplyr::case_when(
        setting == "Community" ~ safe_divide(Tr, eta_NP * diag_RNA),
        setting == "Prison"    ~ m2a(Cov_np_monthly),
        TRUE ~ NA_real_
      ),
      
      # Positive tested people
      M_pos = dplyr::case_when(
        setting == "Community" ~ Cov_np * undiag,
        setting == "Prison"    ~ Cov_np_monthly * prison_pwid_undiag_year_sum,
        TRUE ~ NA_real_
      ),
      M_residual = pmax(M - M_pos, 0),
      
      # Negative/general coverage
      # Assumption: P_PWID/P_fPWID HCV-negative and P_nPWID share monthly coverage.
      Cov_neg_monthly = dplyr::case_when(
        setting == "Community" ~ NA_real_,
        setting == "Prison"    ~ safe_divide(M_residual, prison_total_neg_year_sum),
        TRUE ~ NA_real_
      ),
      Cov_x_monthly = dplyr::case_when(
        setting == "Community" ~ NA_real_,
        setting == "Prison"    ~ Cov_neg_monthly,
        TRUE ~ NA_real_
      ),
      Cov_neg = dplyr::case_when(
        setting == "Community" ~ safe_divide(M_residual, neg),
        setting == "Prison"    ~ m2a(Cov_neg_monthly),
        TRUE ~ NA_real_
      ),
      Cov_x = dplyr::case_when(
        setting == "Community" ~ NA_real_,
        setting == "Prison"    ~ Cov_neg,
        TRUE ~ NA_real_
      ),
      
      # fc is applied to monthly rates -> prison fc uses monthly coverages.
      fc = dplyr::case_when(
        setting == "Community" ~ safe_divide(Cov_neg, Cov_np),
        setting == "Prison"    ~ safe_divide(Cov_neg_monthly, Cov_np_monthly),
        TRUE ~ NA_real_
      ),
      
      # Verification
      M_neg_PWID_fPWID = dplyr::case_when(
        setting == "Community" ~ NA_real_,
        setting == "Prison"    ~ Cov_neg_monthly * prison_pwid_neg_year_sum,
        TRUE ~ NA_real_
      ),
      M_x_P_nPWID = dplyr::case_when(
        setting == "Community" ~ NA_real_,
        setting == "Prison"    ~ Cov_x_monthly * prison_npwid_neg_year_sum,
        TRUE ~ NA_real_
      ),
      M_neg = dplyr::case_when(
        setting == "Community" ~ M_residual,
        setting == "Prison"    ~ M_neg_PWID_fPWID + M_x_P_nPWID,
        TRUE ~ NA_real_
      ),
      M_recovered = M_pos + M_neg,
      Tr_check = dplyr::case_when(
        setting == "Community" ~ eta_NP * Cov_np * diag_RNA,
        setting == "Prison"    ~ NP_eta_P * Cov_np_monthly * prison_pwid_diag_RNA_year_sum,
        TRUE ~ NA_real_
      )
    )
  
  # Add yearly-sum columns to subpopulation stocks
  if (!"diag_RNA_year_sum" %in% names(subpop_stocks)) {
    subpop_stocks <- subpop_stocks |>
      dplyr::mutate(diag_RNA_year_sum = diag_RNA * steps_per_year)
  }
  if (!"undiag_year_sum" %in% names(subpop_stocks)) {
    subpop_stocks <- subpop_stocks |>
      dplyr::mutate(undiag_year_sum = undiag * steps_per_year)
  }
  if (!"neg_year_sum" %in% names(subpop_stocks)) {
    subpop_stocks <- subpop_stocks |>
      dplyr::mutate(neg_year_sum = neg * steps_per_year)
  }
  
  # Subpopulation pass-through table
  subpop_table <- subpop_stocks |>
    add_setting() |>
    dplyr::left_join(
      ccal_setting |>
        dplyr::select(
          year, setting,
          Cov_np_set          = Cov_np,
          Cov_np_monthly_set  = Cov_np_monthly,
          Cov_neg_set         = Cov_neg,
          Cov_neg_monthly_set = Cov_neg_monthly,
          Cov_x_set           = Cov_x,
          Cov_x_monthly_set   = Cov_x_monthly,
          fc_setting          = fc
        ),
      by = c("year", "setting")
    ) |>
    dplyr::mutate(
      # Coverage passed to Param_cal via Ccal_list / fp_yr
      Cov_np = dplyr::case_when(
        setting == "Community"          ~ Cov_np_set,
        pop %in% c("P_PWID", "P_fPWID") ~ Cov_np_set,
        pop == "P_nPWID"                ~ Cov_x_set,
        TRUE ~ NA_real_
      ),
      Cov_neg = dplyr::case_when(
        setting == "Community"          ~ Cov_neg_set,
        pop %in% c("P_PWID", "P_fPWID") ~ Cov_neg_set,
        pop == "P_nPWID"                ~ Cov_x_set,
        TRUE ~ NA_real_
      ),
      # fc passed to HCV_test through fc_full
      fc = dplyr::case_when(
        setting == "Community"          ~ fc_setting,
        pop %in% c("P_PWID", "P_fPWID") ~ safe_divide(Cov_neg_monthly_set,
                                                      Cov_np_monthly_set),
        pop == "P_nPWID"                ~ 1,
        TRUE ~ NA_real_
      ),
      # Display-only predicted treatment
      Tr_predicted = dplyr::case_when(
        setting == "Community"          ~ NP_eta_C * Cov_np * diag_RNA,
        pop %in% c("P_PWID", "P_fPWID") ~ NP_eta_P * Cov_np_monthly_set * diag_RNA_year_sum,
        pop == "P_nPWID"                ~ 0,
        TRUE ~ NA_real_
      ),
      # Display-only people-tested predictions
      M_pos_predicted = dplyr::case_when(
        setting == "Community"          ~ Cov_np * undiag,
        pop %in% c("P_PWID", "P_fPWID") ~ Cov_np_monthly_set * undiag_year_sum,
        pop == "P_nPWID"                ~ 0,
        TRUE ~ NA_real_
      ),
      M_neg_predicted = dplyr::case_when(
        setting == "Community"          ~ Cov_neg * neg,
        pop %in% c("P_PWID", "P_fPWID") ~ Cov_neg_monthly_set * neg_year_sum,
        pop == "P_nPWID"                ~ Cov_x_monthly_set * neg_year_sum,
        TRUE ~ NA_real_
      )
    )
  
  prison_row <- ccal_setting |> dplyr::filter(setting == "Prison")
  prison_tr_monthly_check <- NP_eta_P *
    prison_row$Cov_np_monthly * prison_pwid_diag_RNA_monthly
  
  list(
    subpop_stocks  = subpop_stocks,
    setting_stocks = setting_stocks,
    ccal_setting   = ccal_setting,
    subpop_table   = subpop_table,
    
    prison_pwid_diag_RNA_monthly  = prison_pwid_diag_RNA_monthly,
    prison_pwid_undiag_monthly    = prison_pwid_undiag_monthly,
    prison_pwid_neg_monthly       = prison_pwid_neg_monthly,
    prison_npwid_neg_monthly      = prison_npwid_neg_monthly,
    
    prison_pwid_diag_RNA_year_sum = prison_pwid_diag_RNA_year_sum,
    prison_pwid_undiag_year_sum   = prison_pwid_undiag_year_sum,
    prison_pwid_neg_year_sum      = prison_pwid_neg_year_sum,
    prison_npwid_neg_year_sum     = prison_npwid_neg_year_sum,
    prison_total_neg_year_sum     = prison_total_neg_year_sum,
    
    prison_tr_monthly_check = prison_tr_monthly_check
  )
}


# =============================================================================
# 8. STEP 0: STATUS QUO BASELINE
# =============================================================================

zero_cascade <- lapply(dfList, function(x) x * 0)
zero_fc      <- matrix(0, nrow = POC_AU$npops, ncol = dim(dfList$eta)[3])

Sce_sq <- HCV_np(
  POC_AU,
  best_estimates,
  best_est_pop,
  disease_progress,
  pop_array,
  param_cascade    = dfList,
  param_cascade_sc = zero_cascade,
  fib              = fib,
  modelrun         = "UN",
  proj             = "POC_AU",
  end_Y            = endY,
  cost             = NULL,
  costflow         = NULL,
  costflow_Neg     = NULL,
  fc_sc            = zero_fc
)

cat("Status quo run complete.\n")


# =============================================================================
# 9. STEP 1: CHAINED CALIBRATION (PASS 1)
# =============================================================================
fm <-list("2022" = c(1,1,1,1,1),
          "2023" = c(1,1,1,1,1),
          "2024" = c(1,1,1,1,1))
param_var <- c("tau_ab", "tau_RNA", "tau_poct", "eta")

ntime <- dim(dfList$eta)[3]

Ccal_list <- list()
fc_list   <- list()
cal_list  <- list()



dfList_NP_current <- zero_cascade
dfList_CT_current <- dfList
fc_full <- matrix(0, nrow = POC_AU$npops, ncol = ntime)

Sce_for_stock <- Sce_sq
Sce_chain     <- NULL

alpha_eta <- 1.0
alpha_rna <- 1.0
alpha_ab  <- 1.0

Sce_stock_list      <- list()
Sce_stock_name_list <- list()

for (yr in cal_years) {
  yr_str <- as.character(yr)
  cat(sprintf("\n=============================================================================\n"))
  cat(sprintf("Calibrating NP parameters for %d using stocks from current chained simulation\n", yr))
  cat(sprintf("=============================================================================\n"))
  
  Sce_stock_list[[yr_str]] <- Sce_for_stock
  cat(sprintf("\n[STOCK SOURCE CHECK %s]\n", yr_str))
  
  # A. Calibrate Cov_np and fc using stocks from the current simulation.
  cal <- calibrate_np_coverage(
    Sce = Sce_for_stock, pj = POC_AU, year = yr,
    obs_data = obs_data,
    NP_eta_C = NP_eta_C, NP_eta_P = NP_eta_P
  )
  cal_list[[yr_str]] <- cal
  print_calibration_summary(cal)
  
  Ccal_list[[yr_str]] <- setNames(as.list(cal$subpop_table$Cov_np),
                                  cal$subpop_table$pop)
  fc_list[[yr_str]]   <- make_fc_matrix_for_year(POC_AU, yr,
                                                 cal$subpop_table, ntime)
  
  tt <- year_time_index(POC_AU, yr)
  fc_full[, tt$t0:tt$t1] <- fc_list[[yr_str]][, tt$t0:tt$t1]
  
  # B. Build NP parameter vector for this year.
  fp_yr <- build_fp(yr_str, Ccal_list, fm)
  cat(sprintf("\nfp_yr used by Param_cal in %s:\n", yr_str))
  print(setNames(fp_yr, POC_AU$popNames[1:5]))
  
  cat(sprintf("\nfc used by HCV_test in %s, first month of year:\n", yr_str))
  print(setNames(fc_full[1:5, tt$t0], POC_AU$popNames[1:5]))
  
  # C. Add this year's NP parameters to the cumulative NP cascade.
  for (i in param_var) {
    dfList_NP_current[[i]] <- Param_cal(
      pj           = POC_AU,
      dlist        = dfList_NP_current,
      index        = i,
      frac_testing = frac_test[[yr_str]],
      S_Yint       = yr,
      S_Yend       = yr + 1,
      r_Yend       = yr + 1,
      NPlst        = NPlst,
      fp           = fp_yr
    )
  }
  
  # Preserve prison values forward where needed.
  b_pt   <- ((yr + 1) - POC_AU$cabY) / POC_AU$timestep + 1
  end_dt <- ((yr + 1) - POC_AU$cabY) / POC_AU$timestep
  
  for (i in param_var) {
    dfList_NP_current[[i]] <- carry_forward_prison(
      dlist_np   = dfList_NP_current,
      dlist_base = zero_cascade,
      index      = i,
      b_pt       = b_pt,
      end_dt     = end_dt
    )
  }
  
  cat(sprintf("\nNP parameters after adding %d:\n", yr))
  for (i in param_var) check_param(dfList_NP_current, i, yr, POC_AU)
  
  # D. Apply CT displacement for this year only.
  C_np <- pmin(fp_yr, 1)
  dfList_CT_current <- scale_CT_eta(dfList_CT_current, dfList, dfList_NP_current, 
                                    POC_AU,
                                    yr, yr + 1)
  dfList_CT_current <- scale_CT_RNA(dfList_CT_current, dfList,dfList_NP_current, 
                                    POC_AU,
                                    yr, yr + 1)
  dfList_CT_current <- scale_CT_ab( dfList_CT_current, dfList, dfList_NP_current, POC_AU,
                                    yr, yr + 1)
  
  # E. Run the chained scenario after this year's calibration.
  Sce_chain <- HCV_np(
    POC_AU,
    best_estimates,
    best_est_pop,
    disease_progress,
    pop_array,
    param_cascade    = dfList_CT_current,
    param_cascade_sc = dfList_NP_current,
    fib              = fib,
    modelrun         = "UN",
    proj             = "POC_AU",
    end_Y            = endY,
    cost             = NULL,
    costflow         = NULL,
    costflow_Neg     = NULL,
    fc_sc            = fc_full
  )
  
  cat(sprintf("\nChained scenario run complete after %d calibration.\n", yr))
  Sce_for_stock <- Sce_chain
}

x <- Sce_chain
dfList_NP_2024    <- dfList_NP_current
dfList_CT_NP_2024 <- dfList_CT_current

cat("\nFinal chained NP test scenario complete.\n\n")
print_model_vs_obs(x, POC_AU, obs_data, cal_years)


# =============================================================================
# 10. STEP 2: ITERATIVE REFINEMENT (PASSES 2..max_passes)
# =============================================================================

all_converged <- FALSE
max_passes <- 4
tolerance <- 0.05
for (pass in 2:max_passes) {
  cat(sprintf("\n=== Pass %d ===\n", pass))
  
  # Reset working objects for this pass
  dfList_NP_current <- zero_cascade
  dfList_CT_current <- dfList
  fc_full <- matrix(0, nrow = POC_AU$npops, ncol = ntime)
  Sce_chain_pass <- NULL
  
  # Universal stock source = final scenario from previous pass
  Sce_universal <- x
  
  cal_list_pass       <- list()
  Ccal_list_pass      <- list()
  fc_list_pass        <- list()
  Sce_stock_list_pass <- list()
  
  for (yr in cal_years) {
    yr_str <- as.character(yr)
    cat(sprintf("\n--- Pass %d, calibrating %s ---\n", pass, yr_str))
    
    Sce_stock_list_pass[[yr_str]] <- Sce_universal
    
    # A. Re-calibrate coverage using post-NP stocks from previous pass
    cal <- calibrate_np_coverage(
      Sce = Sce_universal, pj = POC_AU, year = yr,
      obs_data = obs_data,
      NP_eta_C = NP_eta_C, NP_eta_P = NP_eta_P
    )
    cal_list_pass[[yr_str]] <- cal
    
    cat("\nCalibration summary:\n")
    print_calibration_summary(cal)
    
    cat("\nCommon prison negative-coverage audit:\n")
    print(audit_common_prison_negative_coverage(cal, yr), n = Inf)
    
    cat("\nSubpopulation pass-through audit:\n")
    print(audit_subpop_pass_through(cal, yr), n = Inf)
    
    # B. Build year-specific Ccal and fc matrices
    Ccal_list_pass[[yr_str]] <- setNames(
      as.list(cal$subpop_table$Cov_np),
      cal$subpop_table$pop
    )
    fc_list_pass[[yr_str]] <- make_fc_matrix_for_year(
      pj = POC_AU, year = yr,
      subpop_table = cal$subpop_table, ntime = ntime
    )
    
    tt <- year_time_index(POC_AU, yr)
    fc_full[, tt$t0:tt$t1] <- fc_list_pass[[yr_str]][, tt$t0:tt$t1]
    
    # C. Build fp vector for Param_cal
    fp_yr <- build_fp(yr_str, Ccal_list_pass, fm)
    fp_yr <- pmin(pmax(fp_yr, 0), 0.999999)
    
    cat(sprintf("\nfp_yr used by Param_cal in %s:\n", yr_str))
    print(setNames(fp_yr, POC_AU$popNames[1:5]))
    
    cat(sprintf("\nfc used by HCV_test in %s, first month of year:\n", yr_str))
    print(setNames(fc_full[1:5, tt$t0], POC_AU$popNames[1:5]))
    
    # D. Add this year's NP parameters to cumulative NP cascade
    for (i in param_var) {
      dfList_NP_current[[i]] <- Param_cal(
        pj           = POC_AU,
        dlist        = dfList_NP_current,
        index        = i,
        frac_testing = frac_test[[yr_str]],
        S_Yint       = yr,
        S_Yend       = yr + 1,
        r_Yend       = yr + 1,
        NPlst        = NPlst,
        fp           = fp_yr
      )
    }
    
    cat(sprintf("\nNP parameters after Param_cal for %d:\n", yr))
    for (i in param_var) check_param(dfList_NP_current, i, yr, POC_AU)
    
    # E. Carry prison values forward where needed
    b_pt   <- ((yr + 1) - POC_AU$cabY) / POC_AU$timestep + 1
    end_dt <- ((yr + 1) - POC_AU$cabY) / POC_AU$timestep
    
    for (i in param_var) {
      dfList_NP_current[[i]] <- carry_forward_prison(
        dlist_np   = dfList_NP_current,
        dlist_base = zero_cascade,
        index      = i,
        b_pt       = b_pt,
        end_dt     = end_dt
      )
    }
    
    cat(sprintf("\nNP parameters after carry_forward_prison for %d:\n", yr))
    for (i in param_var) check_param(dfList_NP_current, i, yr, POC_AU)
    
    # F. Apply CT displacement for this year only
    C_np <- pmin(pmax(fp_yr, 0), 1)
    dfList_CT_current <- scale_CT_eta(dfList_CT_current, dfList,dfList_NP_current, POC_AU,
                                      yr, yr + 1)
    dfList_CT_current <- scale_CT_RNA(dfList_CT_current, dfList, dfList_NP_current,POC_AU,
                                      yr, yr + 1)
    dfList_CT_current <- scale_CT_ab( dfList_CT_current, dfList, dfList_NP_current, POC_AU,
                                      yr, yr + 1)
  }
  
  # G. Run full scenario with this pass's calibration
  Sce_chain_pass <- HCV_np(
    POC_AU,
    best_estimates,
    best_est_pop,
    disease_progress,
    pop_array,
    param_cascade    = dfList_CT_current,
    param_cascade_sc = dfList_NP_current,
    fib              = fib,
    modelrun         = "UN",
    proj             = "POC_AU",
    end_Y            = endY,
    cost             = NULL,
    costflow         = NULL,
    costflow_Neg     = NULL,
    fc_sc            = fc_full
  )
  
  cat(sprintf("\nScenario run complete for pass %d.\n", pass))
  
  # H. Diagnostics against observed data
  print_model_vs_obs(Sce_chain_pass, POC_AU, obs_data, cal_years)
  
  # I. Convergence check
  cat("\nConvergence check - prison PWID/fPWID diag_RNA stock and treatment:\n")
  all_converged <- TRUE
  prison_pwid_pops <- which(POC_AU$popNames %in% c("P_PWID", "P_fPWID"))
  
  for (yr in cal_years) {
    yr_str <- as.character(yr)
    tt <- year_time_index(POC_AU, yr)
    ts_yr <- tt$t0:tt$t1
    
    sum_cal <- cal_list_pass[[yr_str]]$prison_pwid_diag_RNA_year_sum
    sum_final <- sum(vapply(ts_yr, function(t) {
      sum(vapply(diag_RNA_states, function(s) {
        sum(Sce_chain_pass$allPops[prison_pwid_pops, s, t], na.rm = TRUE)
      }, numeric(1)))
    }, numeric(1)), na.rm = TRUE)
    
    stock_ratio <- sum_final / pmax(sum_cal, 1e-12)
    tr_actual <- sum(Sce_chain_pass$newTreatment_sc[3:5, ts_yr], na.rm = TRUE)
    tr_obs <- obs_data |>
      dplyr::filter(year == yr, setting == "Prison") |>
      dplyr::pull(Tr)
    tr_ratio <- tr_actual / pmax(tr_obs, 1e-12)
    
    cat(sprintf(
      "  %d: stock_ratio=%.3f  Tr_mod=%.0f  Tr_obs=%.0f  Tr_ratio=%.3f  Tr_diff=%+.1f%%\n",
      yr, stock_ratio, tr_actual, tr_obs, tr_ratio,
      100 * (tr_actual - tr_obs) / pmax(tr_obs, 1e-12)
    ))
    
    if (abs(stock_ratio - 1) > tolerance) all_converged <- FALSE
    if (abs(tr_ratio   - 1) > tolerance) all_converged <- FALSE
  }
  
  # J. Update objects for next pass
  x <- Sce_chain_pass
  cal_list       <- cal_list_pass
  Ccal_list      <- Ccal_list_pass
  fc_list        <- fc_list_pass
  Sce_stock_list <- Sce_stock_list_pass
  Sce_np_test    <- Sce_chain_pass
  
  dfList_NP_2024    <- dfList_NP_current
  dfList_CT_NP_2024 <- dfList_CT_current
  
  if (all_converged) {
    cat(sprintf("\nCONVERGED after pass %d\n", pass))
    break
  }
}

if (!all_converged) {
  cat(sprintf(
    "\nDid not fully converge after %d passes within %.0f%% tolerance.\n",
    max_passes, 100 * tolerance
  ))
}



