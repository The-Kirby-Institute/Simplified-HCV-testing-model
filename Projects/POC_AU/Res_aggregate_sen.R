# =============================================================================
# Res_aggregate — v2, with DAA unit-price sensitivity sweep
# =============================================================================
# WHY THE STRUCTURE CHANGED
# ---------------------------------------------------------------------------
# In v1, DAA price was tied to `cost_type`, so each price needed its own
# simulation run. But this script deletes cost_Treatment / cost_Retreat and
# rebuilds DAA cost from treatment VOLUMES:
#
#     cost_TreatDAA = Resflow_year_pop$Treatment[, par_col] * unitC_DAA
#
# Volumes come from the epi output, which is cost-independent, and cost never
# feeds back into model dynamics. DAA price is therefore a pure post-processing
# multiplier. A full 100%->10% sweep costs ~10 extra multiplications, not 10
# extra simulation runs.
#
# NEW LOOP SHAPE
#   load epi once
#     build Resflow_* aggregations once        (basis- and price-independent)
#     for (basis in c("fixednvariable","total"))
#       load + annualise Rescost_dt once       (expensive, price-independent)
#       for (mult in DAA_multipliers)
#         build fibroscan + TreatOther + RetreatOther   (Block C1)
#         build DAA cost at this price                  (Block C2)
#         -> cap -> discount -> cumulate -> save
#
# Only the file read and the timestep->year collapse are hoisted. Every cost
# component that ends up in a saved object is rebuilt inside the price loop
# from an untouched raw copy, so no price level can inherit mutated state from
# another. This matters for fibroscan in particular: that step SUBTRACTS the
# old bundled rate from cost_RNA / cost_POCT, so it must never run twice on
# the same object.
#
# ---------------------------------------------------------------------------
# ASSUMPTION FLAGGED FOR REVIEW
# ---------------------------------------------------------------------------
# In v1 the "total" basis used DAA = 17978.18, exactly half of the
# "fixednvariable" figure of 35956.37. That meant "total at 100%" and
# "fixednvariable at 50%" were the same DAA price, so the DAA axis and the
# diagnosis-cost axis were entangled and a sweep would double-count.
#
# This script assumes 17978.18 was itself an already-discounted figure, and
# sweeps BOTH bases from the same list price. The two bases then differ only
# in diagnosis / compartment cost, which is what you confirmed.
#
# If that is wrong, set USE_SHARED_LIST_PRICE <- FALSE below and each basis
# reverts to its own 100% reference.
#
# ---------------------------------------------------------------------------
# ON "63% off" vs "67% off"
# ---------------------------------------------------------------------------
# Both are included as named points on the grid (0.37 and 0.33), along with
# 0.50 for continuity with the old DAAcost_reduchalf run. Pick the output file
# that matches the discount you meant; no rerun needed either way.
# =============================================================================

gc(); rm(list = ls()); gc()
tic <- proc.time()

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
library("openxlsx")

Rcode        <- file.path(codefun_path, "03. Code")
DataFolder   <- file.path(data_path, "01. DATA/model input")
RDAFolder    <- file.path(data_path, "02. Output")
OutputFolder <- file.path(codefun_path, "Projects/POC_AU/Output")
OutputFig    <- file.path(codefun_path, "Projects/POC_AU/Figs")
Proj_code    <- file.path(codefun_path, paste0("projects/", project_name))

load(file.path(RDAFolder, paste0(project_name, ".rda")))
source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R"))
source(file.path(Proj_code, "/model_timestep.R"))


# =============================================================================
# CONFIGURATION
# =============================================================================
sce_name   <- c("no_np", "foundational", "succession", "sustained", "scaleup")

# Cost bases that genuinely differ in NON-DAA (diagnosis / compartment) cost.
# "DAAcost_reduchalf" is deliberately absent — it is now mult = 0.50.
cost_bases <- c("fixednvariable", "total")

USE_SHARED_LIST_PRICE <- TRUE

# 100% reference (list) price
unitC_DAA_list     <- 35956.37
unitC_secline_list <- 44613.66

# Per-basis fallback if USE_SHARED_LIST_PRICE is FALSE
unitC_by_basis <- list(
  "fixednvariable" = c(DAA = 35956.37*0.37, secline = 44613.66*0.37),
  "total"          = c(DAA = 35956.37*0.37, secline = 44613.66*0.37)
)

# ---- DAA price sweep -------------------------------------------------------
# 100% down to 10% in 10-point steps, plus the named policy points.
DAA_multipliers <- sort(unique(c(seq(1.00, 0.10, by = -0.10),
                                 0.50,    # old DAAcost_reduchalf
                                 0.37,    # "63% off"
                                 0.33)),  # "67% off"
                        decreasing = TRUE)
DAA_multipliers <- round(DAA_multipliers, 4)

# Label used in filenames: 1.00 -> "DAA100", 0.37 -> "DAA037"
mult_label <- function(m) sprintf("DAA%03d", round(m * 100))

cat("DAA price grid:\n")
print(data.frame(label = vapply(DAA_multipliers, mult_label, ""),
                 multiplier = DAA_multipliers,
                 DAA_price = round(unitC_DAA_list * DAA_multipliers, 2)))

AUdiscount <- 0.05
cap        <- 200000000
endY       <- 100

# ---------------------------------------------------------------------------
# FORCE_REBUILD
#
# The resume guard below skips any price level whose output file already
# exists. That is right after a crash, and WRONG after a code change -- the
# stale file survives and the new logic never runs.
#
# Set TRUE whenever the cost logic has changed since the files were written.
# TRUE is the safe default; flip to FALSE only to resume an interrupted run.
# ---------------------------------------------------------------------------
FORCE_REBUILD <- TRUE

unitC_fibroscan_old <- 62.52                        # bundled cost previously in ctau_ag/ctau_poct
unitC_fibroscan_SOC <- c(C = 73.44, P = 82.30)      # SOC: community / prison
unitC_fibroscan_POC <- c(C = 75.59, P = 83.77)      # POC: community / prison
unitC_eta_other_SOC <- c(C = 1341.14, P = 1063.77)  # community / prison
unitC_eta_other_POC <- c(C = 1342.61, P = 1059.82)


# =============================================================================
# HELPERS
# =============================================================================

# NaN-safe element-wise: replace NaN with 0 so NaN + x does not poison the sum
z <- function(d) { m <- as.matrix(d); m[is.nan(m)] <- 0; m }

rda2list <- function(file) {
  e <- new.env(); load(file, envir = e); as.list(e)
}

fmt_time <- function(secs) {
  if (secs < 60)   return(sprintf("%.1fs", secs))
  if (secs < 3600) return(sprintf("%.1fm", secs / 60))
  sprintf("%.2fh", secs / 3600)
}

# Collapse year x population -> year. na.rm is passed explicitly because v1
# used TRUE for the main path and FALSE for the _sc path; that asymmetry is
# preserved rather than silently harmonised.
agg_year <- function(df, na_rm) {
  df %>% as.data.frame() %>% ungroup() %>%
    group_by(year) %>%
    summarise(across(all_of(par_col), ~ sum(.x, na.rm = na_rm))) %>%
    ungroup()
}

add_discount <- function(df) {
  df %>% ungroup() %>%
    mutate(id = year - POC_AU$simY,
           discount = ifelse(id >= 0, (1 + AUdiscount)^id, NA))
}

cumsum_by_pop <- function(df) {
  df %>% as_tibble() %>% arrange(population) %>% group_by(population) %>%
    mutate(across(all_of(par_col), cumsum, .names = "{col}"))
}

# ---------------------------------------------------------------------------
# build_unit_cost_components()
#
# Called once per (basis x scenario x price level), from Block C1.
#
# Takes the RAW annualised cost list for a scenario and returns:
#   $rc                 -- same list, with the old bundled fibroscan rate
#                          unbundled out of cost_RNA / cost_RNA_sc /
#                          cost_POCT / cost_POCT_sc
#   $cost_fibroscan     -- rebuilt at population-specific SOC/POC rates
#   $cost_fibroscan_sc
#   $cost_TreatOther    -- non-DAA treatment cost (eta_other)
#   $cost_RetreatOther
#   $cost_TreatOther_sc
#
# rc_raw is never modified in place -- the subtraction happens on a local copy,
# so this is safe to call repeatedly for every price level.
# ---------------------------------------------------------------------------
build_unit_cost_components <- function(scn, rc_raw, ic) {
  
  rc         <- rc_raw                       # local copy; copy-on-modify
  pop_prefix <- substr(as.character(ic$population), 1, 1)
  
  # ---- non-DAA treatment cost (eta_other) ----
  SOC_eta <- unitC_eta_other_SOC[pop_prefix]
  POC_eta <- unitC_eta_other_POC[pop_prefix]
  
  n_treat    <- as.matrix(Resflow_year_pop[[scn]]$Treatment[,       par_col])
  n_treat_sc <- as.matrix(Resflow_sc_year_pop[[scn]]$Treatment_sc[, par_col])
  n_retreat  <- as.matrix(Resflow_year_pop[[scn]]$Retreat[,         par_col])
  
  # SOC unit x Treatment + POC unit x Treatment_sc
  # (sq has Treatment_sc = 0, so this collapses cleanly to SOC x Treatment)
  cost_TreatOther    <- cbind(ic, as.data.frame(SOC_eta * n_treat + POC_eta * n_treat_sc))
  cost_RetreatOther  <- cbind(ic, as.data.frame(SOC_eta * n_retreat))
  cost_TreatOther_sc <- cbind(ic, as.data.frame(POC_eta * n_treat_sc))
  
  # ---- fibroscan: unbundle old flat rate, rebuild pop-specific ----
  SOC_fib <- unitC_fibroscan_SOC[pop_prefix]
  POC_fib <- unitC_fibroscan_POC[pop_prefix]
  
  n_RNA_total  <- as.matrix(Resflow_year_pop[[scn]]$Testing_RNA[,     par_col])
  n_RNA_sc     <- as.matrix(Resflow_year_pop[[scn]]$Testing_RNA_sc[,  par_col])
  n_POCT_total <- as.matrix(Resflow_year_pop[[scn]]$Testing_POCT[,    par_col])
  n_POCT_sc    <- as.matrix(Resflow_year_pop[[scn]]$Testing_POCT_sc[, par_col])
  
  rc$cost_RNA[,     par_col] <- rc$cost_RNA[,     par_col] - (n_RNA_total  + n_RNA_sc)  * unitC_fibroscan_old
  rc$cost_RNA_sc[,  par_col] <- rc$cost_RNA_sc[,  par_col] -  n_RNA_sc                  * unitC_fibroscan_old
  rc$cost_POCT[,    par_col] <- rc$cost_POCT[,    par_col] - (n_POCT_total + n_POCT_sc) * unitC_fibroscan_old
  rc$cost_POCT_sc[, par_col] <- rc$cost_POCT_sc[, par_col] -  n_POCT_sc                 * unitC_fibroscan_old
  
  cost_fibroscan <- cbind(
    ic,
    as.data.frame((n_RNA_total + n_POCT_total) * SOC_fib +
                    (n_RNA_sc  + n_POCT_sc)    * POC_fib))
  
  cost_fibroscan_sc <- cbind(
    ic,
    as.data.frame((n_RNA_sc + n_POCT_sc) * POC_fib))
  
  list(rc                 = rc,
       cost_fibroscan     = cost_fibroscan,
       cost_fibroscan_sc  = cost_fibroscan_sc,
       cost_TreatOther    = cost_TreatOther,
       cost_RetreatOther  = cost_RetreatOther,
       cost_TreatOther_sc = cost_TreatOther_sc)
}


# =============================================================================
# LOAD EPI ONCE — basis- and price-independent
# =============================================================================
cat("\n=== Loading epi flows (Res_dt_scenarios.rda) ===\n")
epi_env <- new.env()
load(file.path(OutputFolder, paste0(project_name, "Res_dt_scenarios.rda")),
     envir = epi_env)

# Derive sample count from the data rather than trusting POC_AU$numberSamples,
# which is set at runtime in param_sim_scenario and not always re-saved.
.set_cols  <- grep("^set[0-9]+$", names(epi_env$Num_box[[1]]), value = TRUE)
n_samples  <- length(.set_cols)
par_col    <- c("best", paste0("set", seq_len(n_samples)))
cat(sprintf("Detected %d parameter sets\n", n_samples))

missing_scn <- setdiff(sce_name, names(epi_env$Num_box))
if (length(missing_scn)) stop("Scenario(s) missing from epi file: ",
                              paste(missing_scn, collapse = ", "))

# ---- Pre-flight: every Rescost_dt file must exist before we start ----
for (basis in cost_bases) for (scn in sce_name) {
  f <- file.path(OutputFolder,
                 paste0(project_name, "Rescost_dt_", scn, "_", basis, ".rda"))
  if (!file.exists(f)) stop("Missing cost file: ", f)
}
cat("[ok] all Rescost_dt files present\n")


# =============================================================================
# BLOCK A — epi flow aggregation to year (done ONCE for everything)
# =============================================================================
cat("\n=== Aggregating epi flows to annual (once) ===\n")
t_a <- proc.time()["elapsed"]

Resflow_year_pop    <- list()
Resflow_year_all    <- list()
Resflow_sc_year_pop <- list()
Resflow_sc_year_all <- list()

for (scn in sce_name) {
  # --- main flows ---
  src <- epi_env$Resflow_dt[[scn]]
  Resflow_year_pop[[scn]] <- lapply(src, function(d) {
    d %>% as_tibble() %>% ungroup() %>% arrange(population, year) %>%
      group_by(year, population) %>%
      summarise(across(all_of(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
      arrange(year)
  })
  Resflow_year_all[[scn]] <- lapply(Resflow_year_pop[[scn]], agg_year, na_rm = FALSE)
  
  # --- scenario-specific (_sc) flows ---
  src_sc <- epi_env$Resflow_sc_dt[[scn]]
  Resflow_sc_year_pop[[scn]] <- lapply(src_sc, function(d) {
    d %>% as_tibble() %>% ungroup() %>% arrange(population, year) %>%
      group_by(year, population) %>%
      summarise(across(all_of(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop") %>%
      arrange(year)
  })
  Resflow_sc_year_all[[scn]] <- lapply(Resflow_sc_year_pop[[scn]], agg_year, na_rm = FALSE)
  
  cat(sprintf("  %s aggregated\n", scn))
}

# Num_box carried through unchanged for downstream plotting
Res_numbox <- epi_env$Num_box[sce_name]
save(Res_numbox, file = file.path(OutputFolder,
                                  paste0(project_name, "Res_numbox.rda")))
cat(sprintf("  Res_numbox saved (no longer duplicated per cost type)\n"))
rm(Res_numbox); gc()

rm(epi_env); gc()
cat(sprintf("Block A done in %s\n", fmt_time(proc.time()["elapsed"] - t_a)))


# =============================================================================
# MAIN LOOP
# =============================================================================
manifest <- data.frame()
t_global <- proc.time()["elapsed"]

for (basis in cost_bases) {
  
  cat(sprintf("\n=============================================================\n"))
  cat(sprintf("COST BASIS: %s\n", basis))
  cat(sprintf("=============================================================\n"))
  
  if (USE_SHARED_LIST_PRICE) {
    base_DAA     <- unitC_DAA_list
    base_secline <- unitC_secline_list
  } else {
    base_DAA     <- unname(unitC_by_basis[[basis]]["DAA"])
    base_secline <- unname(unitC_by_basis[[basis]]["secline"])
  }
  message(sprintf("100%% reference for '%s': DAA %.2f | secline %.2f",
                  basis, base_DAA, base_secline))
  
  # ---------------------------------------------------------------------------
  # BLOCK B — load and annualise, ONCE per basis
  #
  # Only the two genuinely expensive, genuinely price-independent operations
  # live here: reading the .rda and collapsing timestep -> year. The raw copy
  # is kept pristine; all derived components are built in Block C1.
  # ---------------------------------------------------------------------------
  t_b <- proc.time()["elapsed"]
  
  Rescost_annual_raw <- list()   # annualised, fibroscan NOT yet unbundled
  index_col          <- list()
  
  for (scn in sce_name) {
    cost_env <- new.env()
    load(file.path(OutputFolder,
                   paste0(project_name, "Rescost_dt_", scn, "_", basis, ".rda")),
         envir = cost_env)
    rc <- cost_env$Rescost_dt[[scn]]
    rm(cost_env)
    
    # ---- annualise: timestep -> year, by population ----
    rc <- lapply(rc, function(d) {
      d %>% as_tibble() %>% ungroup() %>% arrange(timestep, population) %>%
        mutate(year = c(rep(seq(1, endY - 2, 1),
                            each = POC_AU$npops * 1 / POC_AU$timestep),
                        rep(endY - 1, POC_AU$npops * (1 / POC_AU$timestep)))) %>%
        group_by(year, population) %>%
        summarise(across(all_of(par_col), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
        arrange(year)
    })
    
    index_col[[scn]] <- cbind.data.frame(
      year       = rc$cost_Treatment$year,
      population = rc$cost_Treatment$population
    )
    
    Rescost_annual_raw[[scn]] <- rc
    rm(rc); gc()
    cat(sprintf("  [%s] %s loaded and annualised\n", basis, scn))
  }
  
  cat(sprintf("Block B (%s) done in %s\n", basis,
              fmt_time(proc.time()["elapsed"] - t_b)))
  
  # ---------------------------------------------------------------------------
  # BLOCK C — price-DEPENDENT loop
  # ---------------------------------------------------------------------------
  for (mult in DAA_multipliers) {
    
    lab       <- mult_label(mult)
    unitC_DAA         <- base_DAA     * mult
    unitC_secline_DAA <- base_secline * mult
    
    cat(sprintf("\n--- %s x %s | DAA %.2f | secline %.2f ---\n",
                basis, lab, unitC_DAA, unitC_secline_DAA))
    t_c <- proc.time()["elapsed"]
    
    out_file <- file.path(
      OutputFolder,
      paste0(project_name, "Res_flowcost_", basis, "_", lab, ".rda"))
    
    # Manifest is recorded BEFORE the skip guard. It is an index of what the
    # run is meant to produce, not a log of what this invocation rewrote --
    # otherwise a fully-resumed run leaves downstream scripts with an empty
    # index and nothing to read.
    manifest <- rbind(manifest, data.frame(
      cost_basis = basis, label = lab, multiplier = mult,
      unitC_DAA = unitC_DAA, unitC_secline = unitC_secline_DAA,
      file = basename(out_file), stringsAsFactors = FALSE))
    
    if (file.exists(out_file) && !FORCE_REBUILD) {
      cat("    already saved; skipping (FORCE_REBUILD is FALSE)\n")
      next
    }
    
    RescostDAA        <- list()
    cost_TreatDAA     <- list()
    cost_RetreatDAA   <- list()
    cost_TreatDAA_sc  <- list()
    cost_totalDAA     <- list()
    cost_TreatOther   <- list()
    cost_RetreatOther <- list()
    cost_TreatOther_sc<- list()
    cost_fibroscan    <- list()
    cost_fibroscan_sc <- list()
    
    for (scn in sce_name) {
      ic <- index_col[[scn]]
      
      # ---- C1: fibroscan + non-DAA treatment cost, rebuilt for this scenario
      #      from the pristine raw annual copy ----
      comp <- build_unit_cost_components(scn, Rescost_annual_raw[[scn]], ic)
      
      rc                      <- comp$rc      # fibroscan already unbundled
      cost_fibroscan[[scn]]    <- comp$cost_fibroscan
      cost_fibroscan_sc[[scn]] <- comp$cost_fibroscan_sc
      cost_TreatOther[[scn]]   <- comp$cost_TreatOther
      cost_RetreatOther[[scn]] <- comp$cost_RetreatOther
      cost_TreatOther_sc[[scn]]<- comp$cost_TreatOther_sc
      rm(comp)
      
      # ---- C2: DAA cost at this price level ----
      cost_TreatDAA[[scn]] <- cbind(
        ic, Resflow_year_pop[[scn]]$Treatment[, par_col] * unitC_DAA)
      cost_RetreatDAA[[scn]] <- cbind(
        ic, Resflow_year_pop[[scn]]$Retreat[, par_col] * unitC_secline_DAA)
      cost_TreatDAA_sc[[scn]] <- cbind(
        ic, Resflow_sc_year_pop[[scn]]$Treatment_sc[, par_col] * unitC_DAA)
      
      cost_totalDAA[[scn]] <- cbind(
        ic,
        cost_TreatDAA[[scn]][, par_col] +
          cost_RetreatDAA[[scn]][, par_col] +
          cost_TreatDAA_sc[[scn]][, par_col])
      
      RescostDAA[[scn]] <- append(
        rc,
        list("cost_TreatDAA"      = cost_TreatDAA[[scn]],
             "cost_RetreatDAA"    = cost_RetreatDAA[[scn]],
             "cost_TreatDAA_sc"   = cost_TreatDAA_sc[[scn]],
             "cost_TreatOther"    = cost_TreatOther[[scn]],
             "cost_RetreatOther"  = cost_RetreatOther[[scn]],
             "cost_TreatOther_sc" = cost_TreatOther_sc[[scn]],
             "cost_totalDAA"      = cost_totalDAA[[scn]],
             "cost_fibroscan"     = cost_fibroscan[[scn]],
             "cost_fibroscan_sc"  = cost_fibroscan_sc[[scn]]))
      rm(rc)
    }
    
    # ---- totals across populations ----
    RescostDAA_totalpop <- lapply(RescostDAA, function(sc_list)
      lapply(sc_list, agg_year, na_rm = TRUE))
    
    # ---- scenario-only (_sc) cost view ----
    # NOTE: v1 assigned cost_ab (not cost_ab_sc) to the "cost_ab_sc" slot.
    # Preserved verbatim to keep results comparable. Set FIX_COST_AB_SC <- TRUE
    # if that was in fact a typo.
    FIX_COST_AB_SC <- FALSE
    
    RescostDAA_sc <- list()
    Rescost_DAA_sc_dt <- list()
    for (scn in sce_name) {
      RescostDAA_sc[[scn]] <- list(
        "cost_TreatDAA"     = cost_TreatDAA_sc[[scn]],
        "cost_TreatOther"   = cost_TreatOther_sc[[scn]],
        "cost_fibroscan_sc" = cost_fibroscan_sc[[scn]],
        "cost_totalDAA"     = cost_TreatDAA_sc[[scn]])
      
      Rescost_DAA_sc_dt[[scn]] <- list(
        "cost_ab_sc"        = if (FIX_COST_AB_SC) RescostDAA[[scn]]$cost_ab_sc
        else               RescostDAA[[scn]]$cost_ab,
        "cost_RNA_sc"       = RescostDAA[[scn]]$cost_RNA_sc,
        "cost_POCT_sc"      = RescostDAA[[scn]]$cost_POCT_sc,
        "cost_TreatDAA"     = cost_TreatDAA_sc[[scn]],
        "cost_TreatOther"   = cost_TreatOther_sc[[scn]],
        "cost_fibroscan_sc" = cost_fibroscan_sc[[scn]],
        "cost_totalDAA"     = cost_TreatDAA_sc[[scn]])
    }
    
    RescostDAA_sc_totalpop <- lapply(RescostDAA_sc, function(sc_list)
      lapply(sc_list, agg_year, na_rm = TRUE))
    
    # ---- drop pre-split treatment costs; DAA + Other now replace them ----
    Rescost_DAA_dt <- RescostDAA
    for (scn in sce_name) {
      Rescost_DAA_dt[[scn]]$cost_Treatment <- NULL
      Rescost_DAA_dt[[scn]]$cost_Retreat   <- NULL
    }
    
    # ---- shift to calendar year, attach discount factor ----
    shift_year <- function(d) {
      d %>% as_tibble() %>%
        mutate(year     = year + POC_AU$cabY - 1,
               id       = year - POC_AU$simY,
               discount = ifelse(id >= 0, (1 + AUdiscount)^id, NA))
    }
    
    Rescost_year    <- lapply(Rescost_DAA_dt,    function(l) lapply(l, shift_year))
    Rescost_sc_year <- lapply(Rescost_DAA_sc_dt, function(l) lapply(l, shift_year))
    
    # ---- cumulative by population ----
    Rescost_yearcum_pop    <- lapply(Rescost_year,    function(l) lapply(l, cumsum_by_pop))
    Rescost_yearcum_sc_pop <- lapply(Rescost_sc_year, function(l) lapply(l, cumsum_by_pop))
    
    # ---- collapse to year, all populations ----
    Rescost_year_all    <- lapply(Rescost_year,    function(l) lapply(l, agg_year, na_rm = TRUE))
    Rescost_year_sc_all <- lapply(Rescost_sc_year, function(l) lapply(l, agg_year, na_rm = FALSE))
    
    # ---- budget cap on DAA spend, then NaN-safe grand totals ----
    # NOTE: cost_total inherited from Rescost_dt still carries the DAA price
    # baked into param_cost and is NOT price-responsive. Both totals below are
    # rebuilt from components so they track unitC_DAA. Downstream CEA must use
    # cost_total_uncapped / cost_total_Cap, never cost_total.
    grand_total <- function(scn, daa_component) {
      cbind(
        year = Rescost_year_all[[scn]]$cost_totalDAA$year,
        as.data.frame(
          z(Rescost_year_all[[scn]][["cost_compartment"]][,  par_col]) +
            z(Rescost_year_all[[scn]][["cost_ab"]][,           par_col]) +
            z(Rescost_year_all[[scn]][["cost_RNA"]][,          par_col]) +
            z(Rescost_year_all[[scn]][["cost_POCT"]][,         par_col]) +
            z(Rescost_year_all[[scn]][["cost_fibroscan"]][,    par_col]) +
            z(Rescost_year_all[[scn]][[daa_component]][,       par_col]) +
            z(Rescost_year_all[[scn]][["cost_TreatOther"]][,   par_col]) +
            z(Rescost_year_all[[scn]][["cost_RetreatOther"]][, par_col]) +
            z(Rescost_year_all[[scn]][["cost_Cured"]][,        par_col])))
    }
    
    for (scn in sce_name) {
      Rescost_year_all[[scn]][["cost_totalDAA_Cap"]] <-
        Rescost_year_all[[scn]][["cost_totalDAA"]] %>%
        mutate(across(all_of(par_col), ~ ifelse(. >= cap, cap, .), .names = "{col}"))
      
      Rescost_year_all[[scn]][["cost_total_Cap"]] <-
        grand_total(scn, "cost_totalDAA_Cap")
      
      Rescost_year_all[[scn]][["cost_total_uncapped"]] <-
        grand_total(scn, "cost_totalDAA")
    }
    
    # ---- discounted annual ----
    discount_cols <- function(d) {
      add_discount(as.data.frame(d)) %>%
        mutate(across(all_of(par_col), ~ . / discount, .names = "{col}"))
    }
    Rescost_disyear_all    <- lapply(Rescost_year_all,    function(l) lapply(l, discount_cols))
    Rescost_disyear_sc_all <- lapply(Rescost_year_sc_all, function(l) lapply(l, discount_cols))
    
    # ---- cumulative, and discounted cumulative ----
    cum_cols <- function(d) {
      d %>% ungroup() %>%
        mutate(across(all_of(par_col), cumsum, .names = "{col}")) %>%
        mutate(id = year - POC_AU$simY,
               discount = ifelse(id >= 0, (1 + AUdiscount)^id, NA))
    }
    Rescost_yearcum_all    <- lapply(Rescost_year_all,    function(l) lapply(l, cum_cols))
    Rescost_yearcum_sc_all <- lapply(Rescost_year_sc_all, function(l) lapply(l, cum_cols))
    
    discum <- function(d) {
      d %>% ungroup() %>% mutate(across(all_of(par_col), ~ . / discount, .names = "{col}"))
    }
    Rescost_discum_all    <- lapply(Rescost_yearcum_all,    function(l) lapply(l, discum))
    Rescost_discum_sc_all <- lapply(Rescost_yearcum_sc_all, function(l) lapply(l, discum))
    
    # ---- save ----
    DAA_multiplier <- mult
    cost_basis     <- basis
    
    save(Resflow_year_pop, Resflow_year_all,
         Resflow_sc_year_pop, Resflow_sc_year_all,
         RescostDAA, RescostDAA_totalpop,
         RescostDAA_sc, RescostDAA_sc_totalpop,
         Rescost_year, Rescost_sc_year,
         Rescost_yearcum_pop, Rescost_yearcum_sc_pop,
         Rescost_year_all, Rescost_year_sc_all,
         Rescost_disyear_all, Rescost_disyear_sc_all,
         Rescost_yearcum_all, Rescost_yearcum_sc_all,
         Rescost_discum_all, Rescost_discum_sc_all,
         AUdiscount, cap,
         cost_basis, DAA_multiplier,
         unitC_DAA, unitC_secline_DAA,
         unitC_fibroscan_old, unitC_fibroscan_SOC, unitC_fibroscan_POC,
         unitC_eta_other_SOC, unitC_eta_other_POC,
         file = out_file)
    
    cat(sprintf("    saved %s  (%s)\n", basename(out_file),
                fmt_time(proc.time()["elapsed"] - t_c)))
    
    rm(RescostDAA, RescostDAA_totalpop, RescostDAA_sc, RescostDAA_sc_totalpop,
       Rescost_DAA_dt, Rescost_DAA_sc_dt,
       Rescost_year, Rescost_sc_year,
       Rescost_yearcum_pop, Rescost_yearcum_sc_pop,
       Rescost_year_all, Rescost_year_sc_all,
       Rescost_disyear_all, Rescost_disyear_sc_all,
       Rescost_yearcum_all, Rescost_yearcum_sc_all,
       Rescost_discum_all, Rescost_discum_sc_all,
       cost_TreatDAA, cost_RetreatDAA, cost_TreatDAA_sc, cost_totalDAA,
       cost_TreatOther, cost_RetreatOther, cost_TreatOther_sc,
       cost_fibroscan, cost_fibroscan_sc)
    gc()
  }
  
  rm(Rescost_annual_raw, index_col)
  gc()
}

# =============================================================================
# MANIFEST — index of everything written
# =============================================================================
manifest_file <- file.path(OutputFolder,
                           paste0(project_name, "DAA_sweep_manifest.csv"))
write.csv(manifest, manifest_file, row.names = FALSE)

cat(sprintf("\n=============================================================\n"))
cat(sprintf("ALL DONE in %s\n", fmt_time(proc.time()["elapsed"] - t_global)))
cat(sprintf("%d output files (%d bases x %d price levels)\n",
            nrow(manifest), length(cost_bases), length(DAA_multipliers)))
cat(sprintf("Manifest: %s\n", manifest_file))
cat(sprintf("=============================================================\n"))
print(manifest)