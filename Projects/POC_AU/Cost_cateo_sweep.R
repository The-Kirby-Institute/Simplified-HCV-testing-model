# =============================================================================
# Cost categories + CEA datasets across the DAA sweep — memory-safe
# =============================================================================
# Reproduces the doc 7 / doc 8 logic, but loads ONE Res_flowcost_*.rda at a
# time, reduces it to the small (summarised) per-scenario datasets, tags them
# with a sensitivity code + label, drops the big object, and moves on. The
# 1001-column raw objects never accumulate across files.
#
# SENSITIVITY SET (each becomes one 'sensitivity' level)
#   fixednvariable x DAA010 ... DAA100   (10 levels; DAA050 = base case)
#   total          x DAA050              (1 level)
#
# NOTE ON popResults_range
#   Your popResults_range appends best/min/max/Med/Mu/q5/q25/q75/q95 while
#   KEEPING the raw par_col draws. Every downstream step here relies on that
#   (e.g. y_cost_disyear_categories sums across par_col after the range call).
#   If a step errors with "column set1 doesn't exist", that assumption is the
#   thing to check first.
# =============================================================================

gc(); rm(list = ls()); gc()

project_name <- "POC_AU"
codefun_path <- "/Users/jjwu/Projects/Simplified-HCV-testing-model"
data_path    <- paste0("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/",
                       "05. PhD Project/Simplified HCV testing model_/Projects/",
                       project_name)

library("readr"); library("dplyr"); library("tidyr"); library("purrr")
library("ggplot2"); library("openxlsx")

Rcode        <- file.path(codefun_path, "03. Code")
RDAFolder    <- file.path(data_path, "02. Output")
OutputFolder <- file.path(codefun_path, "Projects/POC_AU/Output")
OutputFig    <- file.path(codefun_path, "Projects/POC_AU/Figs")
OutputFig_y_cum_avert <- file.path(OutputFig, "y_cum_avert")
Proj_code    <- file.path(codefun_path, paste0("projects/", project_name))

dir.create(OutputFig,             recursive = TRUE, showWarnings = FALSE)
dir.create(OutputFig_y_cum_avert, recursive = TRUE, showWarnings = FALSE)

load(file.path(RDAFolder, paste0(project_name, ".rda")))
source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R"))
source(file.path(Proj_code, "/model_timestep.R"))

endY     <- 100
year_obs <- c(POC_AU$simY + 5 - 1, POC_AU$simY + 10 - 1, POC_AU$simY + 20 - 1)
par_col  <- c("best", paste0("set", seq_len(POC_AU$numberSamples)))

sce_level <- c("no_np", "foundational", "succession", "sustained", "scaleup")
sce_label <- c("(1) No national program", "(2) Foundational implementation",
               "(3) Program succession", "(4) Program sustained",
               "(5) Program scale-up")


# =============================================================================
# SENSITIVITY PLAN — the 11 files to process, with code + pretty label
# =============================================================================
make_plan <- function() {
  sweep_one <- function(ct) {
    tibble::tibble(cost_type = ct, pct = seq(10, 100, 10)) %>%
      mutate(
        daa_label  = sprintf("DAA%03d", pct),
        sens_code  = paste0(cost_type, "_", daa_label),
        sens_label = dplyr::case_when(
          ct == "fixednvariable" & pct == 30  ~ "Base case: 30% PBS-listed DAA price",
          ct == "fixednvariable" & pct == 100 ~ "PBS-listed full DAA price",
          ct == "fixednvariable"              ~ sprintf("PBS-listed DAA price at %d%%", pct),
          ct == "total"          & pct == 30  ~ "NP total program cost (30% DAA)",
          ct == "total"          & pct == 100 ~ "NP total program cost (full DAA)",
          ct == "total"                       ~ sprintf("NP total program cost (%d%% DAA)", pct)
        )
      )
  }
  bind_rows(sweep_one("fixednvariable"), sweep_one("total"))
}

plan <- make_plan()
cat("Sensitivity plan:\n"); print(plan[, c("sens_code", "sens_label")], n = Inf)

# ordering used for factor levels in combined outputs (sweep ascending, base
# case and total last so they read as distinct references)
sens_levels <- plan$sens_code
sens_labels <- plan$sens_label
names(sens_labels) <- sens_levels


# =============================================================================
# PER-FILE BUILDER
# =============================================================================
# Returns, for one (cost_type, daa_label) file, the small tagged datasets.
# Everything returned is already collapsed by popResults_range, so it is safe
# to accumulate across all 11 files.
# =============================================================================
build_one <- function(cost_type, daa_label, sens_code) {
  
  f <- file.path(OutputFolder,
                 paste0(project_name, "Res_flowcost_", cost_type, "_", daa_label, ".rda"))
  if (!file.exists(f)) stop("Missing file: ", f)
  
  e <- new.env()
  load(f, envir = e)                       # Rescost_year_all, Rescost_disyear_all, ...
  Rescost_year_all    <- e$Rescost_year_all
  Rescost_disyear_all <- e$Rescost_disyear_all
  rm(e)
  
  cost_y_categories       <- list()
  cost_disyear_categories <- list()
  
  cost_cols_to_clean <- c("cost_ab", "cost_RNA", "cost_POCT",
                          "cost_compartment", "cost_Cured",
                          "cost_TreatOther", "cost_RetreatOther",
                          "cost_fibroscan", "cost_totalDAA",
                          "cost_totalDAA_Cap")
  
  for (n in names(Rescost_year_all)) {
    
    yr_obj  <- Rescost_year_all[[n]]
    dis_obj <- Rescost_disyear_all[[n]]
    
    # ---- NA / negative cleanup ----
    for (col_name in cost_cols_to_clean) {
      yr_obj [[col_name]][is.na(yr_obj [[col_name]])] <- 0
      dis_obj[[col_name]][is.na(dis_obj[[col_name]])] <- 0
    }
    for (cc in c("cost_Cured", "cost_TreatOther", "cost_RetreatOther")) {
      yr_obj [[cc]][yr_obj [[cc]] < 0] <- 0
      dis_obj[[cc]][dis_obj[[cc]] < 0] <- 0
    }
    
    # ---- helper: build one category, blank zeros to NA, summarise ----
    mk_cat <- function(obj, components) {
      m <- Reduce(`+`, lapply(components, function(cc) obj[[cc]][, par_col]))
      d <- cbind(year = obj[[components[1]]]$year, as.data.frame(m))
      d[d == 0] <- NA
      popResults_range(POC_AU, d, Population = NULL, end_Y = endY - 1)
    }
    mk_single <- function(obj, comp) {
      d <- obj[[comp]]
      d[d == 0] <- NA
      popResults_range(POC_AU, d, Population = NULL, end_Y = endY - 1)
    }
    
    # ---- undiscounted categories ----
    cost_y_categories[[n]][["Diagnosis"]]     <- mk_cat(yr_obj, c("cost_ab","cost_RNA","cost_POCT"))
    cost_y_categories[[n]][["Treatment_cap"]] <- mk_single(yr_obj, "cost_totalDAA_Cap")
    cost_y_categories[[n]][["Treatment"]]     <- mk_single(yr_obj, "cost_totalDAA")
    cost_y_categories[[n]][["Management"]]    <- mk_cat(yr_obj,
                                                        c("cost_compartment","cost_Cured","cost_TreatOther","cost_RetreatOther","cost_fibroscan"))
    cost_y_categories[[n]] <- dplyr::bind_rows(cost_y_categories[[n]], .id = "Categories")
    
    # ---- discounted categories ----
    cost_disyear_categories[[n]][["Diagnosis"]]     <- mk_cat(dis_obj, c("cost_ab","cost_RNA","cost_POCT"))
    cost_disyear_categories[[n]][["Treatment_cap"]] <- mk_single(dis_obj, "cost_totalDAA_Cap")
    cost_disyear_categories[[n]][["Treatment"]]     <- mk_single(dis_obj, "cost_totalDAA")
    cost_disyear_categories[[n]][["Management"]]    <- mk_cat(dis_obj,
                                                              c("cost_compartment","cost_Cured","cost_TreatOther","cost_RetreatOther","cost_fibroscan"))
    cost_disyear_categories[[n]] <- dplyr::bind_rows(cost_disyear_categories[[n]], .id = "Categories")
  }
  
  # ---- cross-scenario discounted turning-point series (per file) ----
  y_disy <- cost_disyear_categories %>%
    dplyr::bind_rows(.id = "scenario") %>%
    ungroup() %>%
    group_by(scenario, year) %>%
    dplyr::summarise(across(all_of(par_col), ~ sum(.x, na.rm = FALSE)), .groups = "drop")
  
  y_disy_range <- y_disy %>%
    tidyr::gather("simulation", "estimate", -c(year, scenario)) %>%
    group_by(year, scenario) %>%
    summarise(
      best = estimate[simulation == "best"][1],
      min  = min(estimate,    na.rm = TRUE),
      max  = max(estimate,    na.rm = TRUE),
      Med  = median(estimate, na.rm = TRUE),
      Mu   = mean(estimate,   na.rm = TRUE),
      q5   = quantile(estimate, 0.025, na.rm = TRUE),
      q25  = quantile(estimate, 0.25,  na.rm = TRUE),
      q75  = quantile(estimate, 0.75,  na.rm = TRUE),
      q95  = quantile(estimate, 0.975, na.rm = TRUE),
      .groups = "drop")
  
  ref <- y_disy %>% filter(scenario == sce_level[1])
  y_turning <- y_disy %>%
    mutate(best_turning = best - ref$best[match(year, ref$year)]) %>%
    select(scenario, year, best_turning, all_of(par_col)) %>%
    mutate(scenario = factor(scenario, levels = sce_level, labels = sce_label))
  
  # ---- tag scenario + sensitivity on the category frames ----
  tag <- function(x) {
    dplyr::bind_rows(x, .id = "scenario") %>%
      mutate(scenario   = factor(scenario, levels = sce_level, labels = sce_label),
             sens_code  = sens_code)
  }
  cost_y   <- tag(cost_y_categories)
  cost_dis <- tag(cost_disyear_categories)
  
  # ---- CEA cost/QALY tables (doc 8 tab_costqaly_lst), per file ----
  # Uses the raw Rescost objects (kept only for the duration of this call).
  cost_qaly        <- list(); cost_qaly_disy    <- list()
  cost_qaly_ycum   <- list(); cost_qaly_disycum <- list()
  
  for (n in names(Rescost_year_all)) {
    for (m in names(Rescost_year_all[[n]])) {
      yd <- Rescost_year_all[[n]][[m]]
      dd <- Rescost_disyear_all[[n]][[m]]
      if (!"year" %in% names(yd) || !"year" %in% names(dd)) next
      yd[yd == 0] <- NA; dd[dd == 0] <- NA
      
      cost_qaly[[n]][[m]]        <- popResults_range(POC_AU, yd, end_Y = endY - 1) %>% as_tibble()
      cost_qaly_disy[[n]][[m]]   <- popResults_range(POC_AU, dd, end_Y = endY - 1) %>% as_tibble()
      cost_qaly_ycum[[n]][[m]]   <- yd %>% filter(year >= 2022) %>%
        mutate(across(all_of(par_col), cumsum)) %>%
        popResults_range(POC_AU, ., end_Y = endY - 1) %>% as_tibble()
      cost_qaly_disycum[[n]][[m]] <- dd %>% filter(year >= 2022) %>%
        mutate(across(all_of(par_col), cumsum)) %>%
        popResults_range(POC_AU, ., end_Y = endY - 1) %>% as_tibble()
    }
  }
  
  flatten_cq <- function(cq) {
    cq %>% purrr::transpose() %>%
      lapply(function(y) bind_rows(y, .id = "scenario") %>%
               mutate(scenario  = factor(scenario, levels = sce_level, labels = sce_label),
                      sens_code = sens_code))
  }
  
  rm(Rescost_year_all, Rescost_disyear_all); gc()
  
  list(
    cost_y            = cost_y,
    cost_dis          = cost_dis,
    y_turning         = y_turning,
    y_disy_range      = y_disy_range %>% mutate(sens_code = sens_code),
    cq_year           = flatten_cq(cost_qaly),
    cq_disy           = flatten_cq(cost_qaly_disy),
    cq_ycum           = flatten_cq(cost_qaly_ycum),
    cq_disycum        = flatten_cq(cost_qaly_disycum)
  )
}


# =============================================================================
# LOOP over the plan, accumulate small results
# =============================================================================
acc <- list(cost_y = list(), cost_dis = list(), y_turning = list(),
            y_range = list(),
            cq_year = list(), cq_disy = list(), cq_ycum = list(), cq_disycum = list())

for (r in seq_len(nrow(plan))) {
  pr <- plan[r, ]
  cat(sprintf("\n[%2d/%2d] %s  (%s)\n", r, nrow(plan), pr$sens_code, pr$sens_label))
  t0 <- proc.time()["elapsed"]
  
  res <- build_one(pr$cost_type, pr$daa_label, pr$sens_code)
  
  acc$cost_y[[pr$sens_code]]     <- res$cost_y
  acc$cost_dis[[pr$sens_code]]   <- res$cost_dis
  acc$y_turning[[pr$sens_code]]  <- res$y_turning %>% mutate(sens_code = pr$sens_code)
  acc$y_range[[pr$sens_code]]    <- res$y_disy_range
  acc$cq_year[[pr$sens_code]]    <- res$cq_year
  acc$cq_disy[[pr$sens_code]]    <- res$cq_disy
  acc$cq_ycum[[pr$sens_code]]    <- res$cq_ycum
  acc$cq_disycum[[pr$sens_code]] <- res$cq_disycum
  
  rm(res); gc()
  cat(sprintf("       done in %.1fs\n", proc.time()["elapsed"] - t0))
}


# =============================================================================
# COMBINE into sensitivity-keyed frames
# =============================================================================
add_label <- function(df) {
  df %>% mutate(sensitivity = factor(sens_code, levels = sens_levels,
                                     labels = sens_labels[sens_levels]))
}

cost_y_categories       <- bind_rows(acc$cost_y)   %>% add_label()
cost_disyear_categories <- bind_rows(acc$cost_dis) %>% add_label()
y_turning_all           <- bind_rows(acc$y_turning) %>% add_label()
y_range_all             <- bind_rows(acc$y_range)   %>% add_label()

# cap / nocap splits (doc 7 tail)
cost_ydaanocap_categories    <- cost_y_categories       %>% filter(Categories != "Treatment_cap")
cost_ydaacap_categories      <- cost_y_categories       %>% filter(Categories != "Treatment")
cost_disydaanocap_categories <- cost_disyear_categories %>% filter(Categories != "Treatment_cap")
cost_disydaacap_categories   <- cost_disyear_categories %>% filter(Categories != "Treatment")


# =============================================================================
# WRITE category xlsx (one workbook each, sensitivity as a column)
# =============================================================================
sel <- function(df) df %>% select(sensitivity, sens_code, scenario, Categories, year,
                                  best, min, max, Med, Mu, q5, q25, q75, q95)

write.xlsx(sel(cost_ydaacap_categories),
           file.path(OutputFig, "cost_y_daacap_sweep.xlsx"))
write.xlsx(sel(cost_ydaanocap_categories),
           file.path(OutputFig, "cost_y_daanocap_sweep.xlsx"))
write.xlsx(sel(cost_disydaacap_categories),
           file.path(OutputFig, "cost_disy_daacap_sweep.xlsx"))
write.xlsx(sel(cost_disydaanocap_categories),
           file.path(OutputFig, "cost_disy_daanocap_sweep.xlsx"))

cat("\nCategory xlsx written.\n")


# =============================================================================
# PER-LEVEL turning plot (safe for any number of levels)
# =============================================================================
col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")

cost_turning_plot <- function(dt, ttl) {
  ggplot(dt, aes(x = year, colour = scenario)) +
    geom_line(aes(y = best_turning)) +
    scale_x_continuous(expand = c(0, 0), limits = c(2022, 2045),
                       breaks = seq(2022, 2045, 1)) +
    scale_colour_manual(name = "Scenarios", values = col_pal) +
    labs(x = "Year", y = "Annual discounted cost vs no program, billions",
         title = ttl) +
    theme_Publication()
}

for (sc in sens_levels) {
  dt <- y_turning_all %>% filter(sens_code == sc)
  if (nrow(dt) == 0) next
  p <- cost_turning_plot(dt, sens_labels[[sc]])
  ggsave(file.path(OutputFig_y_cum_avert, paste0("p_cost_y_turning_", sc, ".png")),
         p, width = 8, height = 8, bg = "white", dpi = 300)
}
cat("Per-level turning plots written.\n")


# =============================================================================
# SAVE combined objects for the downstream figure/table step
# =============================================================================
save(cost_y_categories, cost_disyear_categories,
     cost_ydaacap_categories, cost_ydaanocap_categories,
     cost_disydaacap_categories, cost_disydaanocap_categories,
     y_turning_all, y_range_all,
     acc,                       # cq_* live here, keyed by sens_code then indicator
     plan, sens_levels, sens_labels,
     file = file.path(OutputFolder, paste0(project_name, "cost_categories_sweep.rda")))

cat(sprintf("\nSaved combined datasets: %s\n",
            file.path(OutputFolder, paste0(project_name, "cost_categories_sweep.rda"))))
cat("=== per-file build complete ===\n")