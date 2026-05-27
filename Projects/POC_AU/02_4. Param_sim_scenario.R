# =============================================================================
# Memory-safe simulation loop
# Three changes vs previous version:
#   1. paramDflist dropped — param_cascade now uses scenario_CT_p (per scenario)
#   2. Inner 1000-sample loop is chunked — limits peak HCV_np output in RAM
#   3. Scenario cascade files loaded once per (scn × cost) then immediately freed
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
RDAFolder <- file.path(data_path, "02. Output")
OutputFolder <- file.path(codefun_path, "Projects/POC_AU/Output")
OutputFig    <- file.path(OutputFolder, "Figs")

source(file.path(Rcode, "Functions/HCV_model.R"))
source(file.path(Rcode, "Functions/plotManuscript.R"))
source(file.path(Rcode, "Functions/plotFunctions.R"))
source(file.path(Rcode, "Functions/check_steady.R"))

load(file.path(RDAFolder, paste0(project_name, ".rda")))
load(file.path(RDAFolder, paste0(project_name, "param.rda")))

number_samples       <- 1000
POC_AU$numberSamples <- number_samples
endY                 <- 100
trim_pt              <- 100 * (1 / POC_AU$timestep)

# Trim Param_* lists to simulation horizon
Param_estimates <- lapply(Param_estimates, function(x) x[1:trim_pt, ] %>% as.data.frame)
param_poparray  <- lapply(param_poparray,  function(x) x[, , 1:trim_pt])
gc()

# Scenario metadata: scenario names + scenario_fc (cost-independent so loaded once)
load(file.path(RDAFolder, paste0(project_name, "Simulations_DAAcost_reduchalf.rda")))
sce_name    <- names(scenario_fc)
scenario_fc <- lapply(scenario_fc, function(m) m[, 1:trim_pt])
rm(scenario_cascade, scenario_cascade_CT, Sce_np,
   scenarios, fitted_coverages); gc()


# =============================================================================
# Pre-run checks — fail fast if any required input is missing
# =============================================================================
cat("\n=== Pre-run checks ===\n")

stopifnot("no_np" %in% sce_name)
cat("[ok] 'no_np' in sce_name\n")

for (scn in sce_name) {
  fp_ct <- file.path(RDAFolder, paste0(project_name, "param_scenario_CT_", scn, ".rda"))
  fp_np <- file.path(RDAFolder, paste0(project_name, "param_scenario_",    scn, ".rda"))
  if (!file.exists(fp_ct)) stop("Missing: ", fp_ct)
  if (!file.exists(fp_np)) stop("Missing: ", fp_np)
}
cat("[ok] all param_scenario_CT_<scn>.rda and param_scenario_<scn>.rda exist\n")

cost_types <- c("fixednvariable", "total", "DAAcost_reduchalf")
for (ct in cost_types) {
  fp <- file.path(RDAFolder, paste0(project_name, "param_cost_", ct, ".rda"))
  if (!file.exists(fp)) stop("Missing: ", fp)
}
cat("[ok] all param_cost_<cost>.rda exist\n")

for (scn in sce_name) {
  stopifnot(scn %in% names(scenario_fc))
}
cat("[ok] scenario_fc has all scenario keys\n\n")


# =============================================================================
# Main loop: chunked, memory-safe
# =============================================================================
chunk_size   <- 100
chunk_starts <- seq(1, number_samples, by = chunk_size)

for (cost_type in cost_types) {
  cat(sprintf("\n=============================================================================\n"))
  cat(sprintf("COST TYPE: %s\n", cost_type))
  cat(sprintf("=============================================================================\n"))
  
  load(file.path(RDAFolder, paste0(project_name, "param_cost_", cost_type, ".rda")))
  # -> param_cost, param_cost_flow, param_costflow_Neg, param_QALY, rand_multiply
  
  for (scn in sce_name) {
    cat(sprintf("\n  --- Scenario: %s ---\n", scn))
    
    # Load scenario-specific cascade pairs
    load(file.path(RDAFolder,
                   paste0(project_name, "param_scenario_CT_", scn, ".rda")))
    scenario_ct <- lapply(scenario_ct,
                          function(x) lapply(x, function(y) y[, , 1:trim_pt]))
    gc()
    
    load(file.path(RDAFolder,
                   paste0(project_name, "param_scenario_", scn, ".rda")))
    scenario_p <- lapply(scenario_p,
                         function(x) lapply(x, function(y) y[, , 1:trim_pt]))
    gc()
    
    tic <- proc.time()
    tmp_files <- character(length(chunk_starts))
    
    # Chunked inner loop
    for (ci in seq_along(chunk_starts)) {
      idx <- chunk_starts[ci]:min(chunk_starts[ci] + chunk_size - 1, number_samples)
      
      chunk_result <- vector("list", length(idx))
      for (j in seq_along(idx)) {
        x <- idx[j]
        chunk_result[[j]] <- HCV_np(
          POC_AU, Param_estimates[[x]], Param_Pops[[x]],
          Param_disease_progress[[x]], param_poparray[[x]],
          scenario_ct[[x]],
          param_cascade_sc = scenario_p[[x]],
          fib              = Param_fib[[x]],
          modelrun         = "UN",
          proj             = "POC_AU",
          end_Y            = endY,
          cost             = param_cost[[x]],
          costflow         = param_cost_flow[[x]],
          costflow_Neg     = param_costflow_Neg[[x]],
          fc_sc            = scenario_fc[[scn]],
          fp               = NULL
        )
      }
      
      tmp_files[ci] <- tempfile(
        pattern = sprintf("chunk%02d_%s_%s_", ci, scn, cost_type),
        fileext = ".rda",
        tmpdir  = OutputFolder
      )
      save(chunk_result, file = tmp_files[ci])
      cat(sprintf("    chunk %d/%d (%d samples) written\n",
                  ci, length(chunk_starts), length(idx)))
      rm(chunk_result); gc()
    }
    
    # Reassemble
    cat(sprintf("  Reassembling %d chunks...\n", length(chunk_starts)))
    param_scenario <- vector("list", number_samples)
    for (ci in seq_along(chunk_starts)) {
      idx <- chunk_starts[ci]:min(chunk_starts[ci] + chunk_size - 1, number_samples)
      load(tmp_files[ci])
      for (j in seq_along(idx)) {
        param_scenario[[idx[j]]] <- chunk_result[[j]]
      }
      file.remove(tmp_files[ci])
      rm(chunk_result); gc()
    }
    
    toc <- proc.time() - tic
    cat(sprintf("    %s × %s: %d samples in %.1f sec\n",
                scn, cost_type, number_samples, toc[3]))
    
    save(param_scenario,
         file = file.path(OutputFolder,
                          paste0(project_name, "param_sc_", scn, "_", cost_type, ".rda")))
    cat(sprintf("    Saved: param_sc_%s_%s.rda\n", scn, cost_type))
    
    rm(param_scenario, scenario_ct, scenario_p); gc()
  }
  
  rm(param_cost, param_cost_flow, param_costflow_Neg, param_QALY, rand_multiply); gc()
}

cat("\n=== ALL DONE ===\n")