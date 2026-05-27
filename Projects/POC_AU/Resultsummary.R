# =============================================================================
# Results summary across 5 scenarios — memory-safe + scenario-chunk saving
#
# Main change vs previous version:
#   STEP 1 now treats each scenario as one chunk:
#     1. Run one scenario
#     2. Save one scenario-level epi result file
#     3. Remove that scenario from memory
#     4. After all scenarios finish, reload scenario files and combine
#
# This protects completed scenarios if the script crashes later.
# =============================================================================


gc(); rm(list = ls()); gc()

project_name <- "POC_AU"

codefun_path <- "/Users/jjwu/Projects/Simplified-HCV-testing-model"
data_path    <- paste0("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/",
                       "05. PhD Project/Simplified HCV testing model_/Projects/",
                       project_name)

library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")

Rcode        <- file.path(codefun_path, "03. Code")
DataFolder   <- file.path(data_path, "01. DATA/model input")
RDAFolder    <- file.path(data_path, "02. Output")
OutputFolder <- file.path(codefun_path, "Projects/POC_AU/Output")
OutputFig    <- file.path(codefun_path, "Projects/POC_AU/Figs")
Proj_code    <- file.path(codefun_path, paste0("projects/", project_name))

dir.create(OutputFolder, recursive = TRUE, showWarnings = FALSE)

load(file.path(RDAFolder, paste0(project_name, ".rda")))
source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R"))
source(file.path(Proj_code, "/model_timestep.R"))

endY    <- 100
par_col <- c("best", paste0("set", seq(1, 1000, 1)))


# =============================================================================
# TUNING — adjust to your machine's RAM
# -----------------------------------------------------------------------------
#   16 GB : n_cores = 3,  par_chunk = 100
#   32 GB : n_cores = 5,  par_chunk = 200
#   64 GB+: n_cores = 7,  par_chunk = 300
# If you see "cannot allocate vector" or the machine swaps, lower n_cores first.
# =============================================================================
n_cores   <- 5L
par_chunk <- 200L

n_cores  <- min(n_cores, max(1L, parallel::detectCores() - 1L))
mc_opts  <- list(mc.cores = n_cores, mc.preschedule = TRUE, mc.cleanup = TRUE)
cat(sprintf("Parallel workers: %d (detected %d cores), par_chunk = %d\n",
            n_cores, parallel::detectCores(), par_chunk))


# =============================================================================
# modres.t — LOCAL OVERRIDE
# =============================================================================
modres.t <- function(pg, Best, endYear, allp = NULL) {
  allpop <- as.data.frame.table(Best[[if (!is.null(allp)) allp else "allPops"]])
  
  allpop <- allpop %>%
    mutate(time = rep(seq(1.0, (endYear - pg$timestep), pg$timestep),
                      each = pg$ncomponent * pg$npops),
           cascade = sub("^[^_]*_", "", Var2),
           disease_prog = sub("\\_.*", "", Var2)) %>%
    dplyr::select(-Var3) %>%
    ungroup()
  
  allpop <- allpop %>%
    mutate(time = c(rep(seq(pg$startYear, endYear - 1 * pg$timestep, pg$timestep),
                        each = pg$npops * pg$ncomponent)))
  
  names(allpop) <- c("population", "state", "best", "timestep",
                     "cascade", "disease_prog")
  
  timelong <- seq(pg$startYear, endYear, pg$timestep)
  
  allpop <- allpop %>%
    filter(timestep %in% timelong) %>%
    mutate(year = timestep %/% 1)
}


# =============================================================================
# Timing helpers
# =============================================================================
fmt_time <- function(secs) {
  if (secs < 60)   return(sprintf("%.1fs", secs))
  if (secs < 3600) return(sprintf("%.1fm", secs / 60))
  return(sprintf("%.2fh", secs / 3600))
}

print_checkpoint <- function(label, t_start, n_done, n_total, t_block_start = NULL) {
  elapsed <- as.numeric(proc.time()["elapsed"] - t_start)
  block_elapsed <- if (!is.null(t_block_start)) {
    as.numeric(proc.time()["elapsed"] - t_block_start)
  } else NA
  pct <- 100 * n_done / n_total
  if (n_done > 0 && n_done < n_total) {
    est_remaining <- (elapsed / n_done) * (n_total - n_done)
    cat(sprintf("    [%s] step:%s  total:%s  %d/%d (%.0f%%)  ~%s remaining\n",
                label,
                if (!is.na(block_elapsed)) fmt_time(block_elapsed) else "n/a",
                fmt_time(elapsed), n_done, n_total, pct, fmt_time(est_remaining)))
  } else {
    cat(sprintf("    [%s] step:%s  total:%s  %d/%d (%.0f%%)\n",
                label,
                if (!is.na(block_elapsed)) fmt_time(block_elapsed) else "n/a",
                fmt_time(elapsed), n_done, n_total, pct))
  }
}


# =============================================================================
# Parallel helpers — chunked mclapply with worker-error detection
# =============================================================================
chunked_mclapply <- function(X, FUN, label = "mclapply") {
  n      <- length(X)
  starts <- seq(1L, n, by = par_chunk)
  out    <- vector("list", n)
  
  for (s in starts) {
    idx <- s:min(s + par_chunk - 1L, n)
    res <- do.call(mclapply, c(list(X = X[idx], FUN = FUN), mc_opts))
    
    errs <- which(vapply(res, inherits, logical(1), "try-error"))
    if (length(errs) > 0) {
      stop(sprintf("%s: worker errors at samples %s",
                   label, paste(idx[errs], collapse = ",")))
    }
    out[idx] <- res
    rm(res)
    gc()
  }
  out
}

extract_par_modres <- function(param_scenario, endY, allp = NULL) {
  chunked_mclapply(
    param_scenario,
    FUN = function(x) {
      modres.t(POC_AU, x, endYear = endY, allp = allp) %>%
        tibble::as_tibble() %>%
        select(best)
    },
    label = paste0("extract_par_modres(", if (is.null(allp)) "Num_box" else allp, ")")
  )
}

extract_par_modres_flow <- function(param_scenario, flow_names) {
  chunked_mclapply(
    param_scenario,
    FUN = function(x) {
      out <- lapply(flow_names, function(y) {
        tibble::tibble(best = as.vector(x[[y]]))
      })
      names(out) <- flow_names
      out
    },
    label = "extract_par_modres_flow"
  )
}


# =============================================================================
# STEP 1 — EPI OUTCOMES, saved by scenario, then combined
# =============================================================================
epi_cost_type <- "DAAcost_reduchalf"

cat(sprintf("\n=============================================================\n"))
cat(sprintf("STEP 1: Epi outcomes (cost type: %s)\n", epi_cost_type))
cat(sprintf("=============================================================\n"))

t_global_start <- proc.time()["elapsed"]
t_epi_start    <- proc.time()["elapsed"]

load(file.path(RDAFolder, paste0(project_name, "Simulations_", epi_cost_type, ".rda")))
sce_name <- names(Sce_np)
cat("Scenarios:", paste(sce_name, collapse = ", "), "\n")

# ---- Indicator routing ----
all_flow_indicators <- names(Sce_np[[1]])

exclude_structural  <- c("allPops", "HCVdeathState", "newDeathState",
                         "costPops", "QALYPops", "death_hcv")

cost_flow_names     <- setdiff(grep("^cost", all_flow_indicators, value = TRUE),
                               exclude_structural)
non_cost_flow_names <- setdiff(all_flow_indicators,
                               c(exclude_structural, cost_flow_names))

cat(sprintf("Non-cost flow indicators: %d\n", length(non_cost_flow_names)))
cat(sprintf("Cost flow indicators:     %d\n", length(cost_flow_names)))

.shape_ok <- sapply(c(non_cost_flow_names, cost_flow_names), function(n) {
  d <- dim(Sce_np[[1]][[n]])
  length(d) == 2 && d[2] == 1188
})
if (!all(.shape_ok)) {
  stop("Mis-routed indicator(s): ",
       paste(names(.shape_ok)[!.shape_ok], collapse = ", "))
}
cat("[ok] all flow indicators are 2D [pop x time]\n\n")

# ---- Result object names used for scenario save/combine ----
epi_result_names <- c(
  "Num_box", "pop_N", "commu_N", "prison_N", "prisonPWID_N",
  "Num_diag", "Num_diag_ab", "Num_diag_Treated",
  "Num_chronic_cured", "Num_curInf",
  "Num_dc", "Num_hcc", "Num_lt", "Num_plt",
  "Sce_flow", "tempNOTInfected_subpop",
  "tempChronic_subpop", "tempPrev_subpop",
  "tempNOTInfectedRNA_subpop", "tempPrevRNA_subpop", "HCVInc_subpop",
  "tempNOTInfected_commu", "tempNOTInfected_prison", "tempNOTInfected_prisonPWID",
  "tempPrev_setting", "tempNOTInfectedRNA_commu", "tempNOTInfectedRNA_prison",
  "tempNOTInfectedRNA_prisonPWID", "tempPrevRNA_setting",
  "newInf_commu", "newInf_prison", "newInf_prisonPWID", "HCVInc_setting"
)

epi_scn_file <- function(scn) {
  file.path(
    OutputFolder,
    paste0(project_name, "epiRes_timestep_", scn, "_", epi_cost_type, ".rda")
  )
}

save_epi_scenario <- function(scn) {
  epi_one_scenario <- list(scn = scn)
  
  for (nm in epi_result_names) {
    epi_one_scenario[[nm]] <- get(nm, envir = .GlobalEnv)[[scn]]
  }
  
  f <- epi_scn_file(scn)
  save(epi_one_scenario, file = f)
  cat(sprintf("  Saved scenario result: %s\n", f))
  
  rm(epi_one_scenario)
  invisible(f)
}

drop_epi_scenario_from_memory <- function(scn) {
  for (nm in epi_result_names) {
    tmp <- get(nm, envir = .GlobalEnv)
    tmp[[scn]] <- NULL
    assign(nm, tmp, envir = .GlobalEnv)
    rm(tmp)
  }
  gc()
}

# ---- Containers ----
Num_box <- list(); pop_N <- list(); commu_N <- list(); prison_N <- list(); prisonPWID_N <- list()
Num_diag <- list(); Num_diag_ab <- list(); Num_diag_Treated <- list()
Num_chronic_cured <- list(); Num_curInf <- list()
Num_dc <- list(); Num_hcc <- list(); Num_lt <- list(); Num_plt <- list()
Sce_flow <- list()
tempNOTInfected_subpop <- list(); tempChronic_subpop <- list(); tempPrev_subpop <- list()
tempNOTInfectedRNA_subpop <- list(); tempPrevRNA_subpop <- list(); HCVInc_subpop <- list()
tempNOTInfected_commu <- list(); tempNOTInfected_prison <- list(); tempNOTInfected_prisonPWID <- list()
tempPrev_setting <- list()
tempNOTInfectedRNA_commu <- list(); tempNOTInfectedRNA_prison <- list(); tempNOTInfectedRNA_prisonPWID <- list()
tempPrevRNA_setting <- list()
newInf_commu <- list(); newInf_prison <- list(); newInf_prisonPWID <- list()
HCVInc_setting <- list()

# ---- Run one scenario at a time ----
epi_scn_count <- 0
for (scn in sce_name) {
  epi_scn_count <- epi_scn_count + 1
  
  if (file.exists(epi_scn_file(scn))) {
    cat(sprintf("\n--- [Step 1: %d/%d] Epi: %s already saved; skipping ---\n",
                epi_scn_count, length(sce_name), scn))
    next
  }
  
  cat(sprintf("\n--- [Step 1: %d/%d] Epi: %s ---\n",
              epi_scn_count, length(sce_name), scn))
  t_scn_start <- proc.time()["elapsed"]
  
  load(file.path(OutputFolder,
                 paste0(project_name, "param_sc_", scn, "_", epi_cost_type, ".rda")))
  
  # ---- Num_box ----
  t_block <- proc.time()["elapsed"]
  Num_box[[scn]] <- modres.t(POC_AU, Sce_np[[scn]], endYear = endY) %>%
    tibble::as_tibble()
  
  par_Num_box <- extract_par_modres(param_scenario, endY)
  for (i in 1:length(par_Num_box)) {
    Num_box[[scn]][, paste0("set", i)] <- par_Num_box[[i]]$best
  }
  rm(par_Num_box)
  gc()
  
  Num_box[[scn]] <- Num_box[[scn]] %>%
    select(year, population, state, timestep, cascade, disease_prog,
           best, paste0("set", seq(1, 1000, 1)))
  cat(sprintf("    Num_box built in %s\n", fmt_time(proc.time()["elapsed"] - t_block)))
  
  # ---- Population totals ----
  pop_N[[scn]]        <- N_pop_sum(Num_box[[scn]], pop = NULL,                              param = "y", name_parset = par_col)
  commu_N[[scn]]      <- N_pop_sum(Num_box[[scn]], pop = c("C_PWID", "C_fPWID"),            param = "y", name_parset = par_col)
  prison_N[[scn]]     <- N_pop_sum(Num_box[[scn]], pop = c("P_PWID", "P_fPWID", "P_nPWID"), param = "y", name_parset = par_col)
  prisonPWID_N[[scn]] <- N_pop_sum(Num_box[[scn]], pop = c("P_PWID", "P_fPWID"),            param = "y", name_parset = par_col)
  
  # ---- Cascade counts ----
  Num_diag[[scn]]          <- N_pop_casdisprog(Num_box[[scn]], pop = NULL,
                                               cas = c("diag_RNA", "treat", "treat_f"),
                                               disprog = POC_AU$progress_name[-1],
                                               param = "y", name_parset = par_col)
  Num_diag_ab[[scn]]       <- N_pop_casdisprog(Num_box[[scn]], pop = NULL,
                                               cas = c("diag_ab", "diag_RNA", "treat", "treat_f", "cured"),
                                               disprog = POC_AU$progress_name[-1],
                                               param = "y", name_parset = par_col)
  Num_diag_Treated[[scn]]  <- N_pop_casdisprog(Num_box[[scn]], pop = NULL,
                                               cas = c("treat", "treat_f", "cured"),
                                               disprog = POC_AU$progress_name[-1],
                                               param = "y", name_parset = par_col)
  Num_chronic_cured[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL,
                                               cas = c("cured"),
                                               disprog = POC_AU$progress_name[-1],
                                               param = "y", name_parset = par_col)
  Num_curInf[[scn]]        <- N_pop_casdisprog(Num_box[[scn]], pop = NULL,
                                               cas = c("undiag", "diag_ab", "diag_RNA", "treat", "treat_f"),
                                               disprog = NULL,
                                               param = "y", name_parset = par_col)
  Num_dc[[scn]]  <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, cas = NULL, disprog = c("dc"),  param = "y", name_parset = par_col)
  Num_hcc[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, cas = NULL, disprog = c("hcc"), param = "y", name_parset = par_col)
  Num_lt[[scn]]  <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, cas = NULL, disprog = c("lt"),  param = "y", name_parset = par_col)
  Num_plt[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, cas = NULL, disprog = c("plt"), param = "y", name_parset = par_col)
  
  # ---- Sce_flow non-cost ----
  t_block <- proc.time()["elapsed"]
  Sce_flow[[scn]] <- list()
  for (x in non_cost_flow_names) {
    Sce_flow[[scn]][[x]] <- modres.flow.t(POC_AU, Sce_np[[scn]], endYear = endY, allp = x)
  }
  
  par_Sce_flow <- extract_par_modres_flow(param_scenario, non_cost_flow_names)
  for (i in 1:length(par_Sce_flow)) {
    for (x in non_cost_flow_names) {
      Sce_flow[[scn]][[x]][, paste0("set", i)] <- par_Sce_flow[[i]][[x]]$best
    }
  }
  rm(par_Sce_flow)
  gc()
  cat(sprintf("    Sce_flow (non-cost) built in %s\n",
              fmt_time(proc.time()["elapsed"] - t_block)))
  
  rm(param_scenario)
  gc()
  
  # ---- Prev & inc by subpop ----
  tempNOTInfected_subpop[[scn]] <- Num_box[[scn]] %>%
    filter(disease_prog != "a") %>%
    filter(state == "s") %>%
    group_by(timestep, population) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep, population)
  
  tempChronic_subpop[[scn]] <- Num_box[[scn]] %>%
    filter(disease_prog != "a") %>%
    group_by(timestep, population) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep, population)
  
  tempPrev_subpop[[scn]] <- cbind(
    timestep = pop_N[[scn]]$timestep,
    population = POC_AU$popNames,
    as.data.frame(100 * (pop_N[[scn]][, par_col] - tempNOTInfected_subpop[[scn]][, par_col]) /
                    pop_N[[scn]][, par_col])
  ) %>% tibble::as_tibble()
  
  tempNOTInfectedRNA_subpop[[scn]] <- Num_box[[scn]] %>%
    filter(cascade %in% c("s", "cured")) %>%
    group_by(timestep, population) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep, population)
  
  tempPrevRNA_subpop[[scn]] <- cbind(
    timestep = pop_N[[scn]]$timestep,
    population = POC_AU$popNames,
    as.data.frame(100 * (pop_N[[scn]][, par_col] - tempNOTInfectedRNA_subpop[[scn]][, par_col]) /
                    pop_N[[scn]][, par_col])
  ) %>% tibble::as_tibble()
  
  HCVInc_subpop[[scn]] <- cbind(
    timestep = pop_N[[scn]]$timestep,
    population = POC_AU$popNames,
    as.data.frame(100 * Sce_flow[[scn]]$newInfections[, par_col] / pop_N[[scn]][, par_col])
  ) %>% tibble::as_tibble()
  
  # ---- Prev by setting ----
  tempNOTInfected_commu[[scn]] <- Num_box[[scn]] %>%
    filter(population %in% c("C_PWID", "C_fPWID") & disease_prog == "s") %>%
    group_by(timestep) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep)
  
  tempPrev_setting[[scn]] <- list()
  tempPrev_setting[[scn]][["commu"]] <- cbind(
    timestep = commu_N[[scn]]$timestep,
    as.data.frame(100 * (commu_N[[scn]][, par_col] - tempNOTInfected_commu[[scn]][, par_col]) /
                    commu_N[[scn]][, par_col])
  ) %>% tibble::as_tibble()
  
  tempNOTInfected_prison[[scn]] <- Num_box[[scn]] %>%
    filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID") & disease_prog == "s") %>%
    group_by(timestep) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep)
  
  tempPrev_setting[[scn]][["prisons"]] <- cbind(
    timestep = prison_N[[scn]]$timestep,
    as.data.frame(100 * (prison_N[[scn]][, par_col] - tempNOTInfected_prison[[scn]][, par_col]) /
                    prison_N[[scn]][, par_col])
  ) %>% tibble::as_tibble()
  
  tempNOTInfected_prisonPWID[[scn]] <- Num_box[[scn]] %>%
    filter(population %in% c("P_PWID", "P_fPWID") & disease_prog == "s") %>%
    group_by(timestep) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep)
  
  tempPrev_setting[[scn]][["prisonsPWID"]] <- cbind(
    timestep = prisonPWID_N[[scn]]$timestep,
    as.data.frame(100 * (prisonPWID_N[[scn]][, par_col] - tempNOTInfected_prisonPWID[[scn]][, par_col]) /
                    prisonPWID_N[[scn]][, par_col])
  ) %>% tibble::as_tibble()
  
  tempNOTInfectedRNA_commu[[scn]] <- Num_box[[scn]] %>%
    filter(cascade %in% c("s", "cured") & population %in% c("C_PWID", "C_fPWID")) %>%
    group_by(timestep) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep)
  
  tempPrevRNA_setting[[scn]] <- list()
  tempPrevRNA_setting[[scn]][["commu"]] <- cbind(
    timestep = commu_N[[scn]]$timestep,
    as.data.frame(100 * (commu_N[[scn]][, par_col] - tempNOTInfectedRNA_commu[[scn]][, par_col]) /
                    commu_N[[scn]][, par_col])
  ) %>% tibble::as_tibble()
  
  tempNOTInfectedRNA_prison[[scn]] <- Num_box[[scn]] %>%
    filter(cascade %in% c("s", "cured") & population %in% c("P_PWID", "P_fPWID", "P_nPWID")) %>%
    group_by(timestep) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep)
  
  tempPrevRNA_setting[[scn]][["prisons"]] <- cbind(
    timestep = prison_N[[scn]]$timestep,
    as.data.frame(100 * (prison_N[[scn]][, par_col] - tempNOTInfectedRNA_prison[[scn]][, par_col]) /
                    prison_N[[scn]][, par_col])
  ) %>% tibble::as_tibble()
  
  tempNOTInfectedRNA_prisonPWID[[scn]] <- Num_box[[scn]] %>%
    filter(cascade %in% c("s", "cured") & population %in% c("P_PWID", "P_fPWID")) %>%
    group_by(timestep) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep)
  
  tempPrevRNA_setting[[scn]][["prisonsPWID"]] <- cbind(
    timestep = prisonPWID_N[[scn]]$timestep,
    as.data.frame(100 * (prisonPWID_N[[scn]][, par_col] - tempNOTInfectedRNA_prisonPWID[[scn]][, par_col]) /
                    prisonPWID_N[[scn]][, par_col])
  ) %>% tibble::as_tibble()
  
  # ---- Incidence by setting ----
  newInf_commu[[scn]] <- Sce_flow[[scn]]$newInfections %>%
    filter(population %in% c("C_PWID", "C_fPWID")) %>%
    group_by(timestep) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep)
  
  newInf_prison[[scn]] <- Sce_flow[[scn]]$newInfections %>%
    filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID")) %>%
    group_by(timestep) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep)
  
  newInf_prisonPWID[[scn]] <- Sce_flow[[scn]]$newInfections %>%
    filter(population %in% c("P_PWID", "P_fPWID")) %>%
    group_by(timestep) %>%
    summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
    arrange(timestep)
  
  HCVInc_setting[[scn]] <- list()
  HCVInc_setting[[scn]][["commu"]] <- cbind(
    timestep = commu_N[[scn]]$timestep,
    as.data.frame(100 * (newInf_commu[[scn]][, par_col] / commu_N[[scn]][, par_col]))
  ) %>% tibble::as_tibble()
  HCVInc_setting[[scn]][["prison"]] <- cbind(
    timestep = prison_N[[scn]]$timestep,
    as.data.frame(100 * (newInf_prison[[scn]][, par_col] / prison_N[[scn]][, par_col]))
  ) %>% tibble::as_tibble()
  HCVInc_setting[[scn]][["prisonPWID"]] <- cbind(
    timestep = prisonPWID_N[[scn]]$timestep,
    as.data.frame(100 * (newInf_prisonPWID[[scn]][, par_col] / prisonPWID_N[[scn]][, par_col]))
  ) %>% tibble::as_tibble()
  
  scn_elapsed <- proc.time()["elapsed"] - t_scn_start
  cat(sprintf("  Scenario '%s' epi done in %s\n", scn, fmt_time(scn_elapsed)))
  print_checkpoint(sprintf("STEP 1: epi %s", scn),
                   t_global_start, epi_scn_count, length(sce_name), t_epi_start)
  
  # ---- Save this completed scenario and clear it from memory ----
  save_epi_scenario(scn)
  drop_epi_scenario_from_memory(scn)
  gc()
}

cat(sprintf("\n[STEP 1 SCENARIO RUNS COMPLETE] %s for %d scenarios\n",
            fmt_time(proc.time()["elapsed"] - t_epi_start), length(sce_name)))


# =============================================================================
# STEP 1B — COMBINE SAVED SCENARIO FILES
# =============================================================================
cat("\nCombining saved scenario epi result files...\n")

# Re-create empty combined containers
for (nm in epi_result_names) {
  assign(nm, list())
}

for (scn in sce_name) {
  f <- epi_scn_file(scn)
  if (!file.exists(f)) {
    stop("Missing scenario epi result file: ", f)
  }
  
  load(f)  # loads epi_one_scenario
  
  for (nm in epi_result_names) {
    tmp <- get(nm)
    tmp[[scn]] <- epi_one_scenario[[nm]]
    assign(nm, tmp)
    rm(tmp)
  }
  
  rm(epi_one_scenario)
  gc()
  cat(sprintf("  Combined: %s\n", scn))
}

epi_combined_file <- file.path(
  OutputFolder,
  paste0(project_name, "epiRes_timestep_", epi_cost_type, ".rda")
)

save(Num_box, pop_N, commu_N, prison_N, prisonPWID_N,
     Num_diag, Num_diag_ab, Num_diag_Treated,
     Num_chronic_cured, Num_curInf,
     Num_dc, Num_hcc, Num_lt, Num_plt,
     Sce_flow, tempNOTInfected_subpop,
     tempChronic_subpop, tempPrev_subpop,
     tempNOTInfectedRNA_subpop, tempPrevRNA_subpop, HCVInc_subpop,
     tempNOTInfected_commu, tempNOTInfected_prison, tempNOTInfected_prisonPWID,
     tempPrev_setting, tempNOTInfectedRNA_commu, tempNOTInfectedRNA_prison,
     tempNOTInfectedRNA_prisonPWID, tempPrevRNA_setting,
     newInf_commu, newInf_prison, newInf_prisonPWID, HCVInc_setting,
     file = epi_combined_file)

cat(sprintf("Saved combined epi result: %s\n", epi_combined_file))


# =============================================================================
# Build Res_dt_scenarios from combined epi objects
# =============================================================================
Resflow_dt    <- list()
Resflow_sc_dt <- list()
for (scn in sce_name) {
  Resflow_dt[[scn]] <- list(
    newInfections       = Sce_flow[[scn]]$newInfections,
    HCVdeath            = Sce_flow[[scn]]$newHCVdeaths,
    Treatment           = Sce_flow[[scn]]$newTreatment,
    Retreat             = Sce_flow[[scn]]$newRetreat,
    Testing_ab          = Sce_flow[[scn]]$newTestingAb,
    Testing_RNA         = Sce_flow[[scn]]$newTestingAg,
    Testing_POCT        = Sce_flow[[scn]]$newTestingPOCT,
    Testing_ab_neg      = Sce_flow[[scn]]$newTestingAb_neg,
    Testing_RNA_neg     = Sce_flow[[scn]]$newTestingAg_neg,
    Testing_POCT_neg    = Sce_flow[[scn]]$newTestingPOCT_neg,
    Cured               = Sce_flow[[scn]]$newCured,
    Treatment_sc        = Sce_flow[[scn]]$newTreatment_sc,
    Testing_ab_sc       = Sce_flow[[scn]]$newTestingAb_sc,
    Testing_RNA_sc      = Sce_flow[[scn]]$newTestingAg_sc,
    Testing_POCT_sc     = Sce_flow[[scn]]$newTestingPOCT_sc,
    Testing_ab_sc_neg   = Sce_flow[[scn]]$newTestingAb_sc_neg,
    Testing_RNA_sc_neg  = Sce_flow[[scn]]$newTestingAg_sc_neg,
    Testing_POCT_sc_neg = Sce_flow[[scn]]$newTestingPOCT_sc_neg
  )
  Resflow_sc_dt[[scn]] <- list(
    Treatment_sc        = Sce_flow[[scn]]$newTreatment_sc,
    Testing_ab_sc       = Sce_flow[[scn]]$newTestingAb_sc,
    Testing_RNA_sc      = Sce_flow[[scn]]$newTestingAg_sc,
    Testing_POCT_sc     = Sce_flow[[scn]]$newTestingPOCT_sc,
    Testing_ab_sc_neg   = Sce_flow[[scn]]$newTestingAb_sc_neg,
    Testing_RNA_sc_neg  = Sce_flow[[scn]]$newTestingAg_sc_neg,
    Testing_POCT_sc_neg = Sce_flow[[scn]]$newTestingPOCT_sc_neg
  )
}

res_dt_file <- file.path(OutputFolder, paste0(project_name, "Res_dt_scenarios.rda"))
save(Num_box, Resflow_dt, Resflow_sc_dt, file = res_dt_file)
cat(sprintf("Saved: %s\n", res_dt_file))

rm(Num_box, pop_N, commu_N, prison_N, prisonPWID_N,
   Num_diag, Num_diag_ab, Num_diag_Treated, Num_chronic_cured, Num_curInf,
   Num_dc, Num_hcc, Num_lt, Num_plt,
   Sce_flow, tempNOTInfected_subpop, tempChronic_subpop, tempPrev_subpop,
   tempNOTInfectedRNA_subpop, tempPrevRNA_subpop, HCVInc_subpop,
   tempNOTInfected_commu, tempNOTInfected_prison, tempNOTInfected_prisonPWID,
   tempPrev_setting, tempNOTInfectedRNA_commu, tempNOTInfectedRNA_prison,
   tempNOTInfectedRNA_prisonPWID, tempPrevRNA_setting,
   newInf_commu, newInf_prison, newInf_prisonPWID, HCVInc_setting,
   Resflow_dt, Resflow_sc_dt, Sce_np)
gc()


# =============================================================================
# STEP 2 — COST & QALY
# Already saved by scn x cost_type. Adds skip logic for reruns.
# =============================================================================
t_cost_start    <- proc.time()["elapsed"]
cost_types      <- c("fixednvariable", "total", "DAAcost_reduchalf")
n_cost_iter     <- length(cost_types) * length(sce_name)
cost_iter_count <- 0

for (cost_type in cost_types) {
  cat(sprintf("\n=============================================================\n"))
  cat(sprintf("STEP 2: Cost & QALY for cost type: %s\n", cost_type))
  cat(sprintf("=============================================================\n"))
  
  load(file.path(RDAFolder, paste0(project_name, "Simulations_", cost_type, ".rda")))
  
  for (scn in sce_name) {
    cost_iter_count <- cost_iter_count + 1
    
    rescost_file <- file.path(
      OutputFolder,
      paste0(project_name, "Rescost_dt_", scn, "_", cost_type, ".rda")
    )
    
    if (file.exists(rescost_file)) {
      cat(sprintf("\n--- [Step 2: %d/%d] Cost: %s x %s already saved; skipping ---\n",
                  cost_iter_count, n_cost_iter, scn, cost_type))
      next
    }
    
    cat(sprintf("\n--- [Step 2: %d/%d] Cost: %s x %s ---\n",
                cost_iter_count, n_cost_iter, scn, cost_type))
    t_iter_start <- proc.time()["elapsed"]
    
    load(file.path(OutputFolder,
                   paste0(project_name, "param_sc_", scn, "_", cost_type, ".rda")))
    
    # ---- cost_box ----
    t_block <- proc.time()["elapsed"]
    cost_box <- modres.t(POC_AU, Sce_np[[scn]], endYear = endY, allp = "costPops") %>%
      as.data.frame()
    par_cost_box <- extract_par_modres(param_scenario, endY, allp = "costPops")
    for (i in 1:length(par_cost_box)) {
      cost_box[, paste0("set", i)] <- par_cost_box[[i]]$best
    }
    rm(par_cost_box)
    gc()
    cost_box <- cost_box %>%
      select(year, population, state, timestep, cascade, disease_prog,
             best, paste0("set", seq(1, 1000, 1)))
    cost_box_sum <- N_pop_sum(cost_box, pop = NULL, param = "y", name_parset = par_col)
    cat(sprintf("    cost_box built in %s\n", fmt_time(proc.time()["elapsed"] - t_block)))
    
    # ---- QALY_box ----
    t_block <- proc.time()["elapsed"]
    QALY_box <- modres.t(POC_AU, Sce_np[[scn]], endYear = endY, allp = "QALYPops") %>%
      tibble::as_tibble()
    par_QALY_box <- extract_par_modres(param_scenario, endY, allp = "QALYPops")
    for (i in 1:length(par_QALY_box)) {
      QALY_box[, paste0("set", i)] <- par_QALY_box[[i]]$best
    }
    rm(par_QALY_box)
    gc()
    QALY_box <- QALY_box %>%
      select(year, population, state, timestep, cascade, disease_prog,
             best, paste0("set", seq(1, 1000, 1)))
    QALY_box_sum <- N_pop_sum(QALY_box, pop = NULL, param = "y", name_parset = par_col)
    cat(sprintf("    QALY_box built in %s\n", fmt_time(proc.time()["elapsed"] - t_block)))
    
    # ---- Cost flows ----
    t_block <- proc.time()["elapsed"]
    Sce_flow_cost <- list()
    for (x in cost_flow_names) {
      Sce_flow_cost[[x]] <- modres.flow.t(POC_AU, Sce_np[[scn]], endYear = endY, allp = x)
    }
    par_Sce_flow_cost <- extract_par_modres_flow(param_scenario, cost_flow_names)
    for (i in 1:length(par_Sce_flow_cost)) {
      for (x in cost_flow_names) {
        Sce_flow_cost[[x]][, paste0("set", i)] <- par_Sce_flow_cost[[i]][[x]]$best
      }
    }
    rm(par_Sce_flow_cost)
    gc()
    cat(sprintf("    cost flows built in %s\n", fmt_time(proc.time()["elapsed"] - t_block)))
    
    rm(param_scenario)
    gc()
    
    # ---- Assemble Rescost_dt ----
    Rescost_dt <- list()
    Rescost_dt[[scn]] <- list(
      cost_compartment  = cost_box_sum,
      cost_ab           = Sce_flow_cost$costTestingAb,
      cost_RNA          = Sce_flow_cost$costTestingAg,
      cost_POCT         = Sce_flow_cost$costTestingPOCT,
      cost_Treatment    = Sce_flow_cost$costTreatment,
      cost_Retreat      = Sce_flow_cost$costRetreat,
      cost_Cured        = Sce_flow_cost$costCured,
      cost_ab_sc        = Sce_flow_cost$costTestingAb_sc,
      cost_RNA_sc       = Sce_flow_cost$costTestingAg_sc,
      cost_POCT_sc      = Sce_flow_cost$costTestingPOCT_sc,
      cost_Treatment_sc = Sce_flow_cost$costTreatment_sc,
      QALY_compartment  = QALY_box_sum
    )
    Rescost_dt[[scn]][["cost_total"]] <- Rescost_dt[[scn]]$cost_ab
    for (i in par_col) {
      Rescost_dt[[scn]][["cost_total"]][, i] <-
        Rescost_dt[[scn]]$cost_compartment[, i] +
        Rescost_dt[[scn]]$cost_ab[, i] +
        Rescost_dt[[scn]]$cost_RNA[, i] +
        Rescost_dt[[scn]]$cost_POCT[, i] +
        Rescost_dt[[scn]]$cost_Treatment[, i] +
        Rescost_dt[[scn]]$cost_Retreat[, i] +
        Rescost_dt[[scn]]$cost_Cured[, i]
    }
    
    save(Rescost_dt, file = rescost_file)
    
    iter_elapsed <- proc.time()["elapsed"] - t_iter_start
    cat(sprintf("  '%s x %s' done in %s\n", scn, cost_type, fmt_time(iter_elapsed)))
    cat(sprintf("  Saved: %s\n", rescost_file))
    print_checkpoint(sprintf("STEP 2: %s x %s", scn, cost_type),
                     t_global_start, cost_iter_count, n_cost_iter, t_cost_start)
    
    rm(cost_box, cost_box_sum, QALY_box, QALY_box_sum,
       Sce_flow_cost, Rescost_dt)
    gc()
  }
  
  rm(Sce_np)
  gc()
}

cat(sprintf("\n=============================================================\n"))
cat(sprintf("=== ALL DONE in %s ===\n", fmt_time(proc.time()["elapsed"] - t_global_start)))
cat(sprintf("  Step 1 (epi):  %s\n", fmt_time(t_cost_start - t_epi_start)))
cat(sprintf("  Step 2 (cost): %s\n", fmt_time(proc.time()["elapsed"] - t_cost_start)))
cat(sprintf("  Workers: %d  |  par_chunk: %d\n", n_cores, par_chunk))
cat(sprintf("=============================================================\n"))
