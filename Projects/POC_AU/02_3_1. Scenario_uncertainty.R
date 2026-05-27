# =============================================================================
# Parameter & scenario uncertainty sampling — MERGED SCRIPT
#
# Outputs kept unchanged for downstream scripts:
#   - <project_name>lhs_sampling.rda
#   - <project_name>param.rda
#   - <project_name>paramDflist.rda
#   - <project_name>param_scenario_<scenario>.rda
#   - <project_name>param_cost_<cost_type>.rda
#
# Main update:\n#   - SOC cascade parameters are sampled jointly as competing probabilities.\n#   - Prison-specific annual probabilities are converted back to monthly\n#     probabilities before applying the competing-probability constraint.\n#   - Competing probabilities are only scaled down if their sum exceeds the\n#     allowed total. They are not forced to sum to 1.
# =============================================================================

rm(list = ls())
gc()

project_name <- "POC_AU"
options(digits = 15)

codefun_path <- "/Users/jjwu/Projects/Simplified-HCV-testing-model"
data_path    <- paste0(
  "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/",
  "05. PhD Project/Simplified HCV testing model_/Projects/",
  project_name
)

library("lhs")
library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")
library("ggplot2")
library("viridis")
library("openxlsx")
library("gt")

Rcode        <- file.path(codefun_path, "03. Code")
DataFolder   <- file.path(data_path, "01. DATA/model input")
OutputFolder <- file.path(data_path, "02. Output")
OutputFig    <- file.path(OutputFolder, "Figs")

load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "cali.rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_timev.rda")))
load(file.path(OutputFolder, paste0(project_name, "Simulations_DAAcost_reduchalf.rda")))

source(file.path(Rcode, "Functions/HCV_model.R"))
source(file.path(Rcode, "Functions/plotManuscript.R"))
source(file.path(Rcode, "Functions/plotFunctions.R"))
source(file.path(Rcode, "Functions/check_steady.R"))

runSamples     <- TRUE
saveAsBase     <- TRUE
number_samples <- 1000
start_lower    <- 0.75
start_upper    <- 1.25
end_lower      <- 0.75
end_upper      <- 1.25

POC_AU$numberSamples <- number_samples


# =============================================================================
# Helper functions
# =============================================================================

clamp01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

sample_from_bounds <- function(u, LL, UU) {
  clamp01(u * (UU - LL) + LL)
}

annual_to_monthly_prob <- function(x) {  # Converts annual probability to equivalent monthly probability.\n  # For example, annual p = 1 becomes monthly p = 1, and annual p = 0 remains 0.\n  clamp01(1 - (1 - clamp01(x))^(1 / 12))\n}\n\nget_prison_rows <- function(template) {\n  if (!is.null(rownames(template))) {\n    return(grepl("^P_", rownames(template)))\n  }\n\n  if (exists("POC_AU") && !is.null(POC_AU$popNames)) {\n    return(grepl("^P_", POC_AU$popNames))\n  }\n\n  rep(FALSE, nrow(template))\n}\n\nconvert_prison_annual_to_monthly <- function(p, keys) {\n  # Prison cascade values are stored as annual probabilities but the model\n  # competing-probability constraint is applied at the monthly time step.\n  prison_rows <- get_prison_rows(p[[keys[1]]])\n\n  if (!any(prison_rows)) {\n    return(p)\n  }\n\n  for (key in keys) {\n    p[[key]][prison_rows, ] <- annual_to_monthly_prob(p[[key]][prison_rows, ])\n  }\n\n  p\n}\n\nscale_competing_if_needed <- function(p, keys, target = 1, eps = 1e-12) {\n  # Enforces sum(keys) <= target.\n  # It does NOT force the probabilities to add to target.\n  template <- p[[keys[1]]]\n\n  target_mat <- if (length(target) == 1) {\n    template * 0 + target\n  } else {\n    target\n  }\n\n  target_mat <- pmax(target_mat, 0)\n\n  denom <- Reduce(`+`, lapply(keys, function(k) pmax(p[[k]], 0)))\n  over  <- denom > target_mat & denom > eps\n\n  for (key in keys) {\n    raw <- clamp01(p[[key]])\n\n    raw[over] <- raw[over] / denom[over] * target_mat[over]\n\n    p[[key]] <- clamp01(raw)\n  }\n\n  p\n}\n\nresidual_after_fixed <- function(fixed_list, keys, template) {
  if (is.null(fixed_list)) {
    return(template * 0 + 1)
  }
  
  fixed_sum <- template * 0
  
  for (key in keys) {
    if (!is.null(fixed_list[[key]])) {
      fixed_sum <- fixed_sum + fixed_list[[key]]
    }
  }
  
  pmax(1 - fixed_sum, 0)
}

safe_ratio <- function(num, den, eps = 1e-12) {
  out <- num / den
  out[!is.finite(out) | abs(den) < eps] <- 1
  out
}

check_competing_sum <- function(p, keys, fixed_list = NULL) {
  total <- Reduce(`+`, lapply(keys, function(k) p[[k]]))
  
  if (!is.null(fixed_list)) {
    fixed_total <- total * 0
    
    for (key in keys) {
      if (!is.null(fixed_list[[key]])) {
        fixed_total <- fixed_total + fixed_list[[key]]
      }
    }
    
    total <- total + fixed_total
  }
  
  range(total, na.rm = TRUE)
}


# =============================================================================
# Cascade sampling objects
# =============================================================================
# These are the mutually exclusive cascade components that should add to 1.
# Keep the object names eta, rho, tau_ab, tau_poct and tau_RNA unchanged because
# downstream scripts use these names inside dfList, paramDflist and scenario_p.

cascade_keys <- c("eta", "rho", "tau_ab", "tau_poct", "tau_RNA")

soc_lhs_cols <- setNames(
  paste0("SOC_", cascade_keys),
  cascade_keys
)

np_effect_lhs_cols <- setNames(
  paste0("NP_effect_", cascade_keys),
  cascade_keys
)

extra_lhs_cols <- c(unname(soc_lhs_cols), unname(np_effect_lhs_cols))


# =============================================================================
# STEP 1 — LHS matrix, one call used by all downstream sampling
# =============================================================================

param_constant <- read.csv(
  file.path(DataFolder, "parameters_constants.csv"),
  header = TRUE
) %>%
  dplyr::select(-"X") %>%
  as.data.frame()

# Disease-progress UL/LL bounds
disease_progUL <- read.csv(
  file.path(DataFolder, "diseaseProgress_UU.csv"),
  header = TRUE
)

disease_progressUL <- c(disease_progUL[, -1]) %>%
  as_tibble() %>%
  mutate(
    a_f0   = disease_progress$a_f0   * start_upper,
    f3_hcc = disease_progress$f3_hcc * start_upper
  )

disease_progLL <- read.csv(
  file.path(DataFolder, "diseaseProgress_LL.csv"),
  header = TRUE
)

disease_progressLL <- c(disease_progLL[, -1]) %>%
  as_tibble() %>%
  mutate(
    a_f0   = disease_progress$a_f0   * start_lower,
    f3_hcc = disease_progress$f3_hcc * start_lower
  )

# Cured-progress UL/LL bounds
cured_progUL <- read.csv(
  file.path(DataFolder, "curedProgress_UU.csv"),
  header = TRUE
)

cured_progressUL <- c(cured_progUL[, -1]) %>%
  as_tibble() %>%
  mutate(
    lt_plt             = fib$lt_plt             * start_upper,
    lt_cured_plt_cured = fib$lt_cured_plt_cured * start_upper
  )

cured_progLL <- read.csv(
  file.path(DataFolder, "curedProgress_LL.csv"),
  header = TRUE
)

cured_progressLL <- c(cured_progLL[, -1]) %>%
  as_tibble() %>%
  mutate(
    lt_plt             = fib$lt_plt             * start_lower,
    lt_cured_plt_cured = fib$lt_cured_plt_cured * start_lower
  )

set.seed(123456)
lhs_samples <- randomLHS(
  n = POC_AU$numberSamples,
  k = nrow(param_constant) +
    dim(disease_progress)[2] +
    dim(fib)[2] +
    1 +
    POC_AU$npops +
    length(extra_lhs_cols)
)

colnames(lhs_samples) <- c(
  param_constant$parameter,
  names(disease_progress),
  names(fib),
  "poparray",
  "C_PWID", "C_fPWID", "P_PWID", "P_fPWID", "P_nPWID",
  extra_lhs_cols
)

save(
  lhs_samples,
  file = file.path(OutputFolder, paste0(project_name, "lhs_sampling.rda"))
)


# =============================================================================
# STEP 2 — Non-cascade params: constants, disease, fib, poparray and Pops
# =============================================================================

param_constant[is.na(param_constant$lower), "lower"] <- start_lower
param_constant[is.na(param_constant$upper), "upper"] <- start_upper

param_constant <- param_constant %>%
  mutate(
    lower = ifelse(lower == start_lower, start_lower * value, lower),
    upper = ifelse(upper == start_upper, start_upper * value, upper)
  )

init_pop <- filter(initialPops, parameter == "init_pop")$value

pop_prop <- initialPops %>%
  filter(parameter %in% c("pop_prop1", "pop_prop2", "pop_prop3", "pop_prop4", "pop_prop5")) %>%
  dplyr::select(value) %>%
  unlist() %>%
  as.vector()

popProp <- as.numeric(init_pop) * pop_prop

init_prop_I <- c(
  best_estimates$HCVP1[1],
  best_estimates$HCVP2[1],
  best_estimates$HCVP3[1],
  best_estimates$HCVP4[1],
  best_estimates$HCVP5[1]
)

init_prop_S <- 1 - init_prop_I

estPops <- read.csv(
  file.path(DataFolder, "Estimate_initial_pop.csv"),
  header = TRUE
)

x            <- list()
randompara   <- lhs_samples[, param_constant$parameter]
parampop     <- lhs_samples[, POC_AU$popNames]

y <- matrix(
  NA,
  nrow = POC_AU$npops,
  ncol = length(names(disease_progress))
)
colnames(y) <- names(disease_progress)
Param_disease_progress <- rep(list(y), POC_AU$numberSamples)

z <- matrix(
  NA,
  nrow = POC_AU$npops,
  ncol = length(names(fib))
)
colnames(z) <- names(fib)
Param_fib <- rep(list(z), POC_AU$numberSamples)

disease_progparam <- list()
cured_progparam   <- list()
parampopsize      <- list()

for (i in 1:POC_AU$numberSamples) {
  x[[i]] <- randompara[i, ] *
    (param_constant$upper - param_constant$lower) +
    param_constant$lower
  
  parampopsize[[i]] <- parampop[i, ] *
    (popProp * start_upper - popProp * start_lower) +
    popProp * start_lower
  
  disease_progparam[[i]] <- do.call(
    rbind,
    replicate(
      POC_AU$npops,
      lhs_samples[i, names(disease_progress)],
      simplify = FALSE
    )
  )
  
  cured_progparam[[i]] <- do.call(
    rbind,
    replicate(
      POC_AU$npops,
      lhs_samples[i, names(fib)],
      simplify = FALSE
    )
  )
  
  Param_disease_progress[[i]] <- disease_progparam[[i]] *
    (disease_progressUL - disease_progressLL) +
    disease_progressLL
  
  Param_fib[[i]] <- cured_progparam[[i]] *
    (cured_progressUL - cured_progressLL) +
    cured_progressLL
}

Param_Cparam <- lapply(x, function(m) as_tibble(t(m)))
Param_estimates <- lapply(
  Param_Cparam,
  function(x) x[rep(1, POC_AU$npts), ]
)

param_poparray <- rep(list(pop_array), POC_AU$numberSamples)
pop_array_LL   <- pop_array * start_lower
pop_array_UU   <- pop_array * start_upper

Param_estPops <- list()
Param_Pops    <- list()

for (i in 1:POC_AU$numberSamples) {
  Param_estPops[[i]] <- estPops %>%
    mutate(
      pop_group = rep(
        c(parampopsize[[i]]),
        dim(estPops)[1] / POC_AU$npops
      ),
      SIprop = ifelse(
        estPops$SI == "S",
        rep(init_prop_S, POC_AU$diseaseprogress_n * POC_AU$npops),
        rep(
          init_prop_I,
          POC_AU$ncomponent * POC_AU$npops -
            POC_AU$diseaseprogress_n * POC_AU$npops
        )
      ),
      est_pop = value * pop_group * SIprop
    )
  
  Param_Pops[[i]] <- as.matrix(
    as.data.frame(
      matrix(
        Param_estPops[[i]]$est_pop,
        ncol = POC_AU$ncomponent,
        nrow = POC_AU$npops
      )
    )
  )
}

Param_Pops <- lapply(Param_Pops, "colnames<-", POC_AU$component_name)

for (i in 1:POC_AU$numberSamples) {
  param_poparray[[i]] <- lhs_samples[i, "poparray"] *
    (pop_array_UU - pop_array_LL) +
    pop_array_LL
  
  param_poparray[[i]][param_poparray[[i]] > 1] <- 1
}

save(
  Param_estimates,
  Param_disease_progress,
  param_poparray,
  parampopsize,
  Param_Pops,
  Param_fib,
  lhs_samples,
  file = file.path(OutputFolder, paste0(project_name, "param.rda"))
)

rm(
  Param_estimates,
  Param_disease_progress,
  param_poparray,
  parampopsize,
  Param_Pops,
  Param_fib
)
gc()

# =============================================================================
# STEP 3 — scenarios uncertainty
# =============================================================================
annual_to_monthly_prob <- function(p_annual) {
  p_annual <- pmin(pmax(p_annual, 0), 1)
  1 - (1 - p_annual)^(1 / 12)
}

monthly_to_annual_prob <- function(p_monthly) {
  p_monthly <- pmin(pmax(p_monthly, 0), 1)
  1 - (1 - p_monthly)^12
}


make_default_bounds <- function(param_list, start_lower = 0.95, start_upper = 1.05) {
  
  LL <- lapply(param_list, function(x) {
    pmin(x * start_lower, 1)
  })
  
  UU <- lapply(param_list, function(x) {
    pmin(x * start_upper, 1)
  })
  
  return(list(LL = LL, UU = UU))
}


apply_prison_monthly_bounds <- function(bounds, param_list,
                                        params_monthly_prison = c("eta", "rho", "tau_ab", "tau_RNA", "tau_poct"),
                                        start_lower = 0.95,
                                        start_upper = 1.05) {
  
  for (p in params_monthly_prison) {
    
    if (!p %in% names(param_list)) {
      warning(paste0("Parameter '", p, "' not found in param_list. Skipping."))
      next
    }
    
    x <- param_list[[p]]
    
    if (is.null(rownames(x))) {
      warning(paste0("Parameter '", p, "' has no rownames. Cannot identify prison rows. Skipping."))
      next
    }
    
    prison_rows <- grepl("^P_", rownames(x))
    
    if (!any(prison_rows)) {
      warning(paste0("No prison rows found for parameter '", p, "'. Skipping."))
      next
    }
    
    # Convert prison annual probabilities to monthly probabilities
    x_monthly <- annual_to_monthly_prob(x[prison_rows, , , drop = FALSE])
    
    # Apply uncertainty on monthly scale
    LL_monthly <- pmin(x_monthly * start_lower, 1)
    UU_monthly <- pmin(x_monthly * start_upper, 1)
    
    # Convert back to annual scale
    bounds$LL[[p]][prison_rows, , ] <- monthly_to_annual_prob(LL_monthly)
    bounds$UU[[p]][prison_rows, , ] <- monthly_to_annual_prob(UU_monthly)
  }
  
  return(bounds)
}


make_scenario_bounds <- function(param_list,
                                 start_lower = 0.95,
                                 start_upper = 1.05,
                                 params_monthly_prison = c("eta", "rho", "tau_ab", "tau_RNA", "tau_poct")) {
  
  bounds <- make_default_bounds(
    param_list = param_list,
    start_lower = start_lower,
    start_upper = start_upper
  )
  
  bounds <- apply_prison_monthly_bounds(
    bounds = bounds,
    param_list = param_list,
    params_monthly_prison = params_monthly_prison,
    start_lower = start_lower,
    start_upper = start_upper
  )
  
  return(bounds)
}


paramset_scenario_memsafe <- function(param_template, lhs, dflist_LL, dflist_UU,
                                      n_samples = POC_AU$numberSamples,
                                      lhs_col = "poparray") {
  
  midpoint_params <- c("cured", "lota")
  sampled_params <- c("eta", "rho", "tau_ab", "tau_poct", "tau_RNA")
  
  if (!lhs_col %in% colnames(lhs)) {
    stop(paste0("Column '", lhs_col, "' not found in lhs."))
  }
  
  n_samples <- min(n_samples, nrow(lhs))
  
  out <- vector("list", n_samples)
  
  for (i in seq_len(n_samples)) {
    
    # Start from one template only, not rep(list(...), n)
    tmp <- param_template
    
    for (p in midpoint_params) {
      if (p %in% names(tmp) && p %in% names(dflist_LL) && p %in% names(dflist_UU)) {
        tmp[[p]] <- (dflist_LL[[p]] + dflist_UU[[p]]) / 2
        tmp[[p]] <- pmin(pmax(tmp[[p]], 0), 1)
      }
    }
    
    for (p in sampled_params) {
      if (p %in% names(tmp) && p %in% names(dflist_LL) && p %in% names(dflist_UU)) {
        tmp[[p]] <- lhs[i, lhs_col] * 
          (dflist_UU[[p]] - dflist_LL[[p]]) + 
          dflist_LL[[p]]
        
        tmp[[p]] <- pmin(pmax(tmp[[p]], 0), 1)
      }
    }
    
    out[[i]] <- tmp
    
    rm(tmp)
    
    if (i %% 10 == 0) gc()
  }
  
  return(out)
}



sce_name <- names(scenario_cascade)

start_lower <- 0.95
start_upper <- 1.05

params_monthly_prison <- c("eta", "rho", "tau_ab", "tau_RNA", "tau_poct")

for (i in sce_name) {
  
  message("Processing scenario: ", i)
  
  # -----------------------------
  # Scenario cascade bounds
  # -----------------------------
  scenario_bounds <- make_scenario_bounds(
    param_list = scenario_cascade[[i]],
    start_lower = start_lower,
    start_upper = start_upper,
    params_monthly_prison = params_monthly_prison
  )
  
  scenario_dfList_LL <- scenario_bounds$LL
  scenario_dfList_UU <- scenario_bounds$UU
  
  rm(scenario_bounds)
  gc()
  
  # -----------------------------
  # Generate scenario_p
  # -----------------------------
  scenario_p <- paramset_scenario_memsafe(
    param_template = scenario_cascade[[i]],
    lhs = lhs_samples,
    dflist_LL = scenario_dfList_LL,
    dflist_UU = scenario_dfList_UU,
    n_samples = POC_AU$numberSamples,
    lhs_col = "poparray"
  )
  
  save(
    scenario_p,
    file = file.path(
      OutputFolder,
      paste0(project_name, "param_scenario_", i, ".rda")
    ),
    compress = "xz"
  )
  
  rm(scenario_p, scenario_dfList_LL, scenario_dfList_UU)
  gc()
  
  # -----------------------------
  # CT cascade bounds
  # -----------------------------
  CT_bounds <- make_scenario_bounds(
    param_list = scenario_cascade_CT[[i]],
    start_lower = start_lower,
    start_upper = start_upper,
    params_monthly_prison = params_monthly_prison
  )
  
  param_CT_dfList_LL <- CT_bounds$LL
  param_CT_dfList_UU <- CT_bounds$UU
  
  rm(CT_bounds)
  gc()
  
  # -----------------------------
  # Generate scenario_ct
  # -----------------------------
  scenario_ct <- paramset_scenario_memsafe(
    param_template = scenario_cascade_CT[[i]],
    lhs = lhs_samples,
    dflist_LL = param_CT_dfList_LL,
    dflist_UU = param_CT_dfList_UU,
    n_samples = POC_AU$numberSamples,
    lhs_col = "poparray"
  )
  
  save(
    scenario_ct,
    file = file.path(
      OutputFolder,
      paste0(project_name, "param_scenario_CT_", i, ".rda")
    ),
    compress = "xz"
  )
  
  rm(scenario_ct, param_CT_dfList_LL, param_CT_dfList_UU)
  gc()
}

# =============================================================================
# STEP 5 — Cost uncertainty per cost_type
# =============================================================================

cost_types <- c("fixednvariable", "total", "DAAcost_reduchalf")

for (cost_type in cost_types) {
  if (cost_type == "fixednvariable") {
    cost_dir <- file.path(DataFolder, "cost")
  } else if (cost_type == "total") {
    cost_dir <- file.path(DataFolder, "cost/sensitivity_total")
  } else if (cost_type == "DAAcost_reduchalf") {
    cost_dir <- file.path(DataFolder, "cost/sensitivity_DAAcost_reduchalf")
  }
  
  files <- list.files(path = cost_dir, pattern = "*.csv")
  
  costdfList <- lapply(files, function(f) {
    df <- read.csv(file.path(cost_dir, f), header = TRUE)
    df <- df[, -1]
    df <- df %>% as_tibble()
    as.matrix(df, nrow = npops, ncol = length(.) + 1)
  })
  
  names(costdfList) <- gsub(".csv", "", files, fixed = TRUE)
  
  cost_state   <- costdfList$state
  costflow     <- list(costdfList$costFlow, costdfList$costFlow_POCRNA)
  costflow_Neg <- list(costdfList$costFlow_NEG, costdfList$`costFlow_POCRNA _NEG`)
  
  set.seed(123456)
  rand_multiply <- runif(number_samples, 0.9, 1.1)
  
  param_cost <- lapply(
    rand_multiply,
    function(x) lapply(costdfList, function(y) y * x)
  )
  
  for (i in 1:length(rand_multiply)) {
    names(param_cost[[i]]) <- names(costdfList)
    param_cost[[i]]$QALY <- lhs_samples[i, "poparray"] *
      (costdfList$QALYPops_UU - costdfList$QALYPops_LL) +
      costdfList$QALYPops_LL
  }
  
  param_cost_flow    <- list()
  param_costflow_Neg <- list()
  param_QALY         <- list()
  
  tic <- proc.time()
  
  for (i in 1:number_samples) {
    param_cost_flow[[i]] <- list(
      param_cost[[i]]$costFlow,
      param_cost[[i]]$costFlow_POCRNA
    )
    
    param_costflow_Neg[[i]] <- list(
      param_cost[[i]]$costFlow_NEG,
      param_cost[[i]]$`costFlow_POCRNA _NEG`
    )
    
    param_QALY[[i]] <- param_cost[[i]]$QALY
  }
  
  toc <- proc.time() - tic
  cat(sprintf("Cost loop %s done in %.1f sec\n", cost_type, toc[3]))
  
  save(
    param_cost,
    param_cost_flow,
    param_costflow_Neg,
    param_QALY,
    rand_multiply,
    file = file.path(OutputFolder, paste0(project_name, "param_cost_", cost_type, ".rda"))
  )
  
  cat(sprintf("Saved param_cost_%s.rda\n", cost_type))
}

