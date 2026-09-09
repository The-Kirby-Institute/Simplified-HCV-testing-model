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
  fixed <- tibble::tibble(
    cost_type  = "fixednvariable",
    pct        = seq(10, 100, 10),
    daa_label  = sprintf("DAA%03d", pct)
  ) %>%
    mutate(
      sens_code  = paste0(cost_type, "_", daa_label),
      sens_label = dplyr::case_when(
        pct == 30  ~ "Base case: 30% PBS-listed DAA price",
        pct == 100 ~ "PBS-listed full DAA price",
        TRUE       ~ sprintf("PBS-listed DAA price at %d%%", pct)
      )
    )
  
  total <- tibble::tibble(
    cost_type  = "total",
    pct        = 30,
    daa_label  = "DAA030",
    sens_code  = "total_DAA030",
    sens_label = "NP total program cost"
  )
  
  bind_rows(fixed, total)
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

acc$cq_disycum$fixednvariable_DAA030$cost_total_Cap
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
  p <- ggplot(dt, aes(x = year, colour = scenario)) +
    geom_line(aes(y = best), size = 1.05) +
    scale_x_continuous(expand = c(0, 0), limits = c(2022, 2041),
                       breaks = seq(2022, 2041, 1)) +
    scale_colour_manual(name = "Scenarios", values = col_pal) +
    labs(x = "Year", y = "Incremental costs (in millions)",
         title = ttl) +
    scale_y_continuous(expand = c(0,0), limits = c(-200000000, 100000000), breaks = seq(-200000000, 100000000,10000000), 
                       labels = seq(-200000000, 100000000,10000000)/1000000) + 
    theme_Publication()
  return(p)  
}

base_cumy_turning <- lapply(sens_levels, function(x) acc$cq_disycum[[x]]$cost_total_Cap%>%
  filter(scenario == "(1) No national program")%>%select(scenario, year, best))
names(base_cumy_turning) <- sens_levels
cum_y_turning <- lapply(sens_levels, function(x) acc$cq_disycum[[x]]$cost_total_Cap%>%
  select(scenario, year, best)%>%
  mutate(best = best - base_cumy_turning[[x]]$best))
names(cum_y_turning) <- sens_levels

sens_labels
cost_turning_plot(cum_y_turning$fixednvariable_DAA010, sens_labels[1])
turning_plot <- lapply(1: length(sens_levels), function(x) 
  p <- cost_turning_plot(cum_y_turning[[x]], sens_labels[x])
)

lapply((1:length(sens_levels)), function(x)  
  ggsave(file.path(OutputFig_y_cum_avert, paste0("p_cost_y_turning_", sens_levels[x], ".png")),
                   turning_plot[[x]], width = 8, height = 8, bg = "white", dpi = 300)
)
cat("Per-level turning plots written.\n")

print(dt[[1]])
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



# =============================================================================
# Cost-category figures + CEA/ICER tables across the DAA sweep
# =============================================================================
# Consumes  POC_AUcost_categories_sweep.rda  (from cost_categories_sweep.R).
#
# Produces, per sensitivity level (11 of them):
#   - stacked cost-category bar at 20y, capped and uncapped        (2 PNG each)
#   - cost-saving-by-category bar vs no program                    (1 PNG each)
# and, across all levels:
#   - CEA/ICER table (cost from capped category sums, QALY from the
#     cumulative discounted QALY), grouped by sensitivity           (gt .docx)
#   - ICER long table                                              (CSV)
#
# CONVENTIONS (locked in earlier)
#   Cost  = capped category sum: Diagnosis + Treatment_cap + Management,
#           cumulated from 2022 on the DISCOUNTED series. Not cost_total_Cap.
#   QALY  = cumulative discounted QALY_compartment (CEAanalysis$disycum).
#   Headline point estimate = median across parameter sets, matching the
#           published tables (their "best" column was formatC(Med, ...)).
# =============================================================================

# =============================================================================
# Cost-category figures + CEA/ICER tables across the DAA sweep
# =============================================================================
# Consumes  POC_AUcost_categories_sweep.rda  (from cost_categories_sweep.R).
#
# Produces, per sensitivity level (11 of them):
#   - stacked cost-category bar at 20y, capped and uncapped        (2 PNG each)
#   - cost-saving-by-category bar vs no program                    (1 PNG each)
# and, across all levels:
#   - CEA/ICER table (cost from capped category sums, QALY from the
#     cumulative discounted QALY), grouped by sensitivity           (gt .docx)
#   - ICER long table                                              (CSV)
#
# CONVENTIONS (locked in earlier)
#   Cost  = capped category sum: Diagnosis + Treatment_cap + Management,
#           cumulated from 2022 on the DISCOUNTED series. Not cost_total_Cap.
#   QALY  = cumulative discounted QALY_compartment (CEAanalysis$disycum).
#   Headline point estimate = median across parameter sets, matching the
#           published tables (their "best" column was formatC(Med, ...)).
# =============================================================================

gc(); rm(list = ls()); gc()

project_name <- "POC_AU"
codefun_path <- "/Users/jjwu/Projects/Simplified-HCV-testing-model"

library("dplyr"); library("tidyr"); library("purrr")
library("ggplot2"); library("gt")

OutputFolder <- file.path(codefun_path, "Projects/POC_AU/Output")
OutputFig    <- file.path(codefun_path, "Projects/POC_AU/Figs")
RDAFolder    <- file.path(paste0("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/",
                                 "05. PhD Project/Simplified HCV testing model_/Projects/",
                                 project_name), "02. Output")

load(file.path(RDAFolder, paste0(project_name, ".rda")))
source(file.path(codefun_path, "03. Code/Functions/plotManuscript.R"))
source(file.path(codefun_path, "03. Code/Functions/plotFunctions.R"))

load(file.path(OutputFolder, paste0(project_name, "cost_categories_sweep.rda")))
# -> cost_*categories, y_turning_all, y_range_all, acc, plan, sens_levels, sens_labels

par_col   <- c("best", paste0("set", seq_len(POC_AU$numberSamples)))
set_cols  <- paste0("set", seq_len(POC_AU$numberSamples))
WTP       <- 50000
bench_scn <- "(1) No national program"

sce_label <- c("(1) No national program", "(2) Foundational implementation",
               "(3) Program succession", "(4) Program sustained",
               "(5) Program scale-up")

timeframes      <- c(5, 10, 20)
timeframe_names <- c("5y", "10y", "20y")
year_obs20      <- POC_AU$simY + 20 - 1


# =============================================================================
# CUMULATION HELPERS
# =============================================================================
# Summary cumulation (best/q5/q95) per category — drives the stacked bars and
# the cost-saving bars. Mirrors document-5 x_catcost.
cat_cum_summary <- function(df) {
  df %>%
    filter(year >= 2022) %>%
    group_by(sens_code, scenario, Categories) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(across(c(best, q5, q95), cumsum)) %>%
    ungroup()
}

# Per-set cumulation, then sum across categories — drives the ICER.
cost_cum_sets <- function(df) {
  df %>%
    filter(year >= 2022) %>%
    group_by(sens_code, scenario, Categories) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(across(all_of(par_col), cumsum)) %>%
    ungroup() %>%
    group_by(sens_code, scenario, year) %>%
    summarise(across(all_of(par_col), ~ sum(.x)), .groups = "drop")   # na.rm=FALSE, matches doc 5
}

cap_summary  <- cat_cum_summary(cost_disydaacap_categories)     # Diagnosis, Treatment_cap, Management
nocap_summary<- cat_cum_summary(cost_disydaanocap_categories)   # Diagnosis, Treatment,     Management
cap_sets     <- cost_cum_sets(cost_disydaanocap_categories)       # per-set capped total


# =============================================================================
# 1. PER-LEVEL STACKED COST-CATEGORY BARS (20y)  — capped and uncapped
# =============================================================================
grey3 <- c("grey10", "grey40", "grey80")

stacked_bar <- function(dat, ttl, top_cat = NULL, thin_frac = 0.05) {
  
  # total per scenario, to know where the bar top is and what's "thin"
  dat <- dat %>%
    group_by(scenario) %>%
    mutate(total = sum(best),
           frac  = best / total) %>%
    ungroup()
  
  # which category is the thin top one: use the name you pass, else auto-detect
  is_thin <- if (!is.null(top_cat)) dat$Categories == top_cat else dat$frac < thin_frac
  
  dat_big  <- dat[!is_thin, ]
  dat_thin <- dat[ is_thin, ]
  
  ggplot(dat, aes(fill = Categories, y = best, x = scenario)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = grey3) +
    theme_Publication(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_y_continuous(limit = c(0, 5e9),
                       breaks = seq(0, 5e9, 5e8),
                       labels = seq(0, 5e9, 5e8) / 1e6) +
    labs(y = "Cost (discounted, millions)", x = NULL) +
    # big segments: centred as before
    geom_text(data = dat_big,
              aes(y = best,
                  label = paste0(format(round(best / 1e6, 1), nsmall = 1), "m"),
                  group = Categories),
              position = position_stack(vjust = 0.5), size = 5) +
    # thin top segment: pulled above the bar top, no overlap
    geom_text(data = dat_thin,
              aes(x = scenario, y = total,
                  label = paste0(format(round(best / 1e6, 1), nsmall = 1), "m")),
              vjust = -0.4, size = 4.5, inherit.aes = FALSE) +
    ggtitle(ttl)
}
sen_maintext <- sens_labels
sen_maintext[[3]] <- "Main analysis"
for (sc in sens_levels) {
  lab <- sen_maintext[[sc]]
  
  d_cap <- cap_summary %>% filter(sens_code == sc, year == year_obs20) %>% arrange(Categories)
  d_noc <- nocap_summary %>% filter(sens_code == sc, year == year_obs20) %>% arrange(Categories)
  if (nrow(d_cap) == 0) next
  
  ggsave(file.path(OutputFig, paste0("cost_catego_20y_", sc, "_cap.png")),
         stacked_bar(d_cap, paste0(lab, " (DAA capped)")),
         width = 10, height = 8, bg = "white", dpi = 300)
  ggsave(file.path(OutputFig, paste0("cost_catego_20y_", sc, "_nocap.png")),
         stacked_bar(d_noc, lab),
         width = 10, height = 8, bg = "white", dpi = 300)
}

# =============================================================================
# 2. PER-LEVEL COST-SAVING-BY-CATEGORY BARS (vs no program, 20y)
# =============================================================================
# incre_best = -(scenario - no_np), per category, at 2041.
incre_by_cat <- function(summary_df, treat_name) {
  ref <- summary_df %>%
    filter(year == year_obs20, scenario == bench_scn) %>%
    select(sens_code, Categories, ref_best = best, ref_q5 = q5, ref_q95 = q95)
  
  summary_df %>%
    filter(year == year_obs20) %>%
    left_join(ref, by = c("sens_code", "Categories")) %>%
    mutate(incre_best = -(best - ref_best),
           incre_q5   = -(q5  - ref_q5),
           incre_q95  = -(q95 - ref_q95),
           Categories = if_else(Categories == treat_name, "Treatment", as.character(Categories))) %>%
    filter(scenario != bench_scn) %>%
    ungroup()
}

incre_cap <- incre_by_cat(cap_summary,   "Treatment_cap")

saving_bar <- function(dat, ttl) {
  net <- dat %>% group_by(scenario) %>%
    summarise(net = sum(incre_best),
              y   = sum(incre_best[incre_best > 0]) + 8e6, .groups = "drop")
  ggplot(dat, aes(x = scenario, y = incre_best, fill = Categories)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = c("grey30", "grey40", "grey80")) +
    theme_Publication(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "right", legend.direction = "vertical") +
    scale_y_continuous(labels = function(v) v / 1e6) +
    labs(y = "Cost saving (discounted, millions)", x = "Scenarios") +
    geom_text(aes(label = paste0(format(round(incre_best / 1e6, 1), nsmall = 1), "m")),
              position = position_stack(vjust = 0.5), size = 5, check_overlap = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(data = net, aes(x = scenario, y = y,
                              label = paste0("Net: ", format(round(net / 1e6, 1), nsmall = 1), "m")),
              inherit.aes = FALSE, fontface = "bold", size = 4) +
    ggtitle(ttl)
}

for (sc in sens_levels) {
  d <- incre_cap %>% filter(sens_code == sc)
  if (nrow(d) == 0) next
  ggsave(file.path(OutputFig, paste0("cost_saving_20y_", sc, ".png")),
         saving_bar(d, sens_labels[[sc]]),
         width = 12, height = 8, bg = "white", dpi = 300)
}
cat("Per-level cost-saving bars written.\n")


# =============================================================================
# 3. CEA / ICER
# =============================================================================
# Cost = capped category total (cap_sets, cumulative discounted).
# QALY = cumulative discounted QALY_compartment.
# Reporting:
#   Costs, QALYs, incremental cost, QALYs gained -> BEST (q2.5, q97.5)
#   ICER                                          -> MEDIAN (q2.5, q97.5)
# The reference scenario (no program) keeps absolute Cost/QALY; its
# incremental columns and ICER are NA.
# =============================================================================


get_qaly <- function(sc) acc$cq_disycum[[sc]]$QALY_compartment   # scenario,year,par_col,...

#### PSA ####
psa_rows <- list()
sens_labels <- c(
  "fixednvariable_DAA010" = "10% of PBS-listed DAA price",
  "fixednvariable_DAA020" = "20% of PBS-listed DAA price",
  "fixednvariable_DAA030" = "Base case: 30% of PBS-listed DAA price",
  "fixednvariable_DAA040" = "40% of PBS-listed DAA price",
  "fixednvariable_DAA050" = "50% of PBS-listed DAA price",
  "fixednvariable_DAA060" = "60% of PBS-listed DAA price",
  "fixednvariable_DAA070" = "70% of PBS-listed DAA price",
  "fixednvariable_DAA080" = "80% of PBS-listed DAA price",
  "fixednvariable_DAA090" = "90% of PBS-listed DAA price",
  "fixednvariable_DAA100" = "PBS-listed full DAA price",
  "total_DAA030"          = "Australian National Program total costs and 30% of PBS-listed DAA price"
)
for (sc in sens_levels) {
  qaly_all <- get_qaly(sc)
  if (is.null(qaly_all)) { warning("No QALY for ", sc); next }
  
  for (ti in seq_along(timeframes)) {
    yr <- POC_AU$simY + timeframes[ti] - 1
    
    cost_at <- cap_sets %>% filter(sens_code == sc, year == yr)
    qaly_at <- qaly_all %>% filter(year == yr)
    if (nrow(cost_at) == 0 || nrow(qaly_at) == 0) next
    
    ref_cost <- cost_at %>% filter(scenario == bench_scn)
    ref_qaly <- qaly_at %>% filter(scenario == bench_scn)
    if (nrow(ref_cost) == 0 || nrow(ref_qaly) == 0) next
    rc <- as.numeric(ref_cost[1, par_col])
    rq <- as.numeric(ref_qaly[1, par_col])
    
    for (slab in sce_label) {
      cc <- cost_at %>% filter(scenario == slab)
      qq <- qaly_at %>% filter(scenario == slab)
      if (nrow(cc) == 0 || nrow(qq) == 0) next
      
      vc <- as.numeric(cc[1, par_col])   # element 1 = best; 2..n = sets
      vq <- as.numeric(qq[1, par_col])
      
      is_bench <- identical(slab, bench_scn)
      
      # pairing guard: draws must align 1:1 with the benchmark
      if (length(vc) != length(rc) || length(vq) != length(rq)) {
        warning("Length mismatch: ", sc, " / ", slab, " @ ", yr); next
      }
      
      n <- length(vc)
      psa_rows[[length(psa_rows) + 1L]] <- data.frame(
        sens_code = sc,
        timeframe = timeframes[ti],
        year      = yr,
        scenario  = slab,
        is_bench  = is_bench,
        draw      = seq_len(n),               # 1 = best/deterministic
        is_best   = c(TRUE, rep(FALSE, n - 1)),
        cost      = vc,
        qaly      = vq,
        inc_cost  = vc - rc,
        inc_qaly  = vq - rq,
        icer      = (vc - rc) / (vq - rq),
        stringsAsFactors = FALSE
      )
    }
  }
}


psa_df <- dplyr::bind_rows(psa_rows)

psa_plot_dt <- psa_df %>%
  filter(timeframe == 20) %>%
  filter(!is_bench) %>%
  mutate(
    Cost_cap   = inc_cost,
    QALY       = inc_qaly,
    outline    = ifelse(is_best, 1L, 0L),
    sens_label = factor(recode(sens_code, !!!sens_labels),
                        levels = unname(sens_labels))
  )



PSA <- ggplot(psa_plot_dt, aes(y = Cost_cap/1000000, x = QALY)) +
  geom_point(aes(colour = scenario)) +
  facet_grid(sens_label ~ scenario,
             labeller = labeller(sens_label = label_wrap_gen(width = 18))) +   # rows = sens, cols = scenario
  scale_color_manual(name = "Scenarios",
                     values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
  geom_point(data = psa_plot_dt %>% filter(outline == 1),
             color = "gray50", size = 1.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  geom_vline(xintercept = 0, color = "black", linewidth = 1) +
  labs(colour = "Scenarios", x = "QALY", y = "Costs (millions)") +
  theme_bw() +
  scale_y_continuous(breaks = seq(-300, 300, 100),
                     labels = seq(-300, 300, 100)) +
  scale_x_continuous(breaks = seq(-1000, 3000, 1000)) +
  coord_cartesian(xlim = c(-1000, 3000),
                  ylim = c(-300, 300)) + 
  theme(
    panel.background = element_rect(colour = "white"),
    plot.background  = element_rect(colour = "white"),
    panel.border     = element_rect(fill = NA, colour = "black"),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.title       = element_text(face = "bold", size = 12),
    axis.title.y     = element_text(angle = 90, vjust = 1),
    axis.title.x     = element_text(vjust = -0.2),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                               size = 9, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black"),
    strip.background   = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    strip.text.y.right = element_text(size = 9, colour = "black", face = "bold", angle = 0),
    strip.text.x       = element_text(size = 11, colour = "black", face = "bold"),
    strip.clip         = "off",
    plot.margin        = unit(c(10, 12, 5, 5), "mm"),  # horizontal row labels, readable
    legend.text      = element_text(size = 14, face = "bold"),
    legend.key       = element_rect(colour = NA),
    legend.position  = "",
    legend.title     = element_text(face = "bold", size = 14),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  #guides(color = guide_legend(direction = "vertical", override.aes = list(size = 2)),
  #       fill = "none") +
  geom_abline(aes(slope = 50000/1000000, intercept = 0, linetype = "WTP: A$50,000"),
              colour = "black") +
  scale_linetype_manual(name = "", values = c(2),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  stat_ellipse(color = "gray50", alpha = 0.7, linewidth = 0.6,
               show.legend = FALSE, level = 0.95)

ggsave(file.path(OutputFig, paste0("PSA.psa.pdf")), plot = PSA,
       width = 13, height = 22, units = "in",
       device = cairo_pdf, limitsize = FALSE)

#### CEA ####
cea_rows <- list()

for (sc in sens_levels) {
  qaly_all <- get_qaly(sc)
  if (is.null(qaly_all)) { warning("No QALY for ", sc); next }
  
  for (ti in seq_along(timeframes)) {
    yr <- POC_AU$simY + timeframes[ti] - 1
    
    cost_at <- cap_sets %>% filter(sens_code == sc, year == yr)
    qaly_at <- qaly_all %>% filter(year == yr)
    if (nrow(cost_at) == 0 || nrow(qaly_at) == 0) next
    
    ref_cost <- cost_at %>% filter(scenario == bench_scn)
    ref_qaly <- qaly_at %>% filter(scenario == bench_scn)
    if (nrow(ref_cost) == 0 || nrow(ref_qaly) == 0) next
    rc <- as.numeric(ref_cost[1, par_col])
    rq <- as.numeric(ref_qaly[1, par_col])
    
    for (slab in sce_label) {
      cc <- cost_at %>% filter(scenario == slab)
      qq <- qaly_at %>% filter(scenario == slab)
      if (nrow(cc) == 0 || nrow(qq) == 0) next
      
      vc <- as.numeric(cc[1, par_col])   # element 1 = best; 2..n = sets
      vq <- as.numeric(qq[1, par_col])
      
      is_bench <- identical(slab, bench_scn)
      
      if (is_bench) {
        inc_cost_best <- NA; inc_cost_q2.5 <- NA; inc_cost_q97.5 <- NA
        inc_qaly_best <- NA; inc_qaly_q2.5 <- NA; inc_qaly_q97.5 <- NA
        icer_med <- NA; icer_q2.5 <- NA; icer_q97.5 <- NA
        prob_ce <- NA; dominant <- NA
      } else {
        d_cost <- vc - rc
        d_qaly <- vq - rq
        icer   <- d_cost / d_qaly
        s_ic <- d_cost[-1]; s_iq <- d_qaly[-1]; s_ir <- icer[-1]
        
        inc_cost_best  <- d_cost[1]
        inc_cost_q2.5  <- quantile(s_ic, 0.025, na.rm = TRUE)
        inc_cost_q97.5 <- quantile(s_ic, 0.975, na.rm = TRUE)
        inc_qaly_best  <- d_qaly[1]
        inc_qaly_q2.5  <- quantile(s_iq, 0.025, na.rm = TRUE)
        inc_qaly_q97.5 <- quantile(s_iq, 0.975, na.rm = TRUE)
        icer_med   <- median(s_ir, na.rm = TRUE)
        icer_q2.5  <- quantile(s_ir, 0.025, na.rm = TRUE)
        icer_q97.5 <- quantile(s_ir, 0.975, na.rm = TRUE)
        prob_ce  <- mean(s_ir < WTP & s_iq > 0, na.rm = TRUE)
        dominant <- mean(s_ic < 0   & s_iq > 0, na.rm = TRUE)
      }
      
      cea_rows[[length(cea_rows) + 1]] <- data.frame(
        sens_code = sc, sensitivity = sens_labels[[sc]],
        scenario = slab, timeframe = timeframe_names[ti],
        
        # absolute, BEST (q2.5, q97.5)
        cost_best = vc[1],
        cost_q2.5 = quantile(vc[-1], 0.025, na.rm = TRUE),
        cost_q97.5 = quantile(vc[-1], 0.975, na.rm = TRUE),
        qaly_best = vq[1],
        qaly_q2.5 = quantile(vq[-1], 0.025, na.rm = TRUE),
        qaly_q97.5 = quantile(vq[-1], 0.975, na.rm = TRUE),
        
        # incremental, BEST (q2.5, q97.5)
        inc_cost_best, inc_cost_q2.5, inc_cost_q97.5,
        inc_qaly_best, inc_qaly_q2.5, inc_qaly_q97.5,
        
        # ICER, MEDIAN (q2.5, q97.5)
        icer_med, icer_q2.5, icer_q97.5,
        prob_ce, dominant,
        stringsAsFactors = FALSE)
    }
  }
}

cea_full <- bind_rows(cea_rows) %>%
  mutate(sensitivity = factor(sensitivity, levels = sens_labels[sens_levels]),
         scenario    = factor(scenario,    levels = sce_label))
View(cea_full)
write.csv(cea_full,
          file.path(OutputFolder, paste0(project_name, "DAA_sweep_CEA.csv")),
          row.names = FALSE)


# =============================================================================
# ICER vs DAA price figure (20y) — fixednvariable sweep, non-reference scenarios
# =============================================================================
col_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")

sweep_only <- cea_full %>%
  filter(timeframe == "20y",
         grepl("^fixednvariable_", sens_code),
         scenario != bench_scn) %>%
  mutate(pct = as.integer(sub(".*DAA0*", "", sens_code)))

if (nrow(sweep_only) > 0) {
  p_icer <- ggplot(sweep_only, aes(pct, icer_med, colour = scenario)) +
    geom_hline(yintercept = WTP, linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    scale_colour_manual(name = "Scenarios", values = col_pal) +
    scale_x_reverse(breaks = seq(100, 10, -10)) +
    labs(x = "DAA price (% of PBS list)", y = "ICER (A$ / QALY, median)",
         caption = sprintf("Dashed = WTP A$%s/QALY. fixednvariable basis, 20-year horizon.",
                           format(WTP, big.mark = ","))) +
    theme_Publication()
  ggsave(file.path(OutputFig, "ICER_vs_DAA_price_20y.png"),
         p_icer, width = 11, height = 7, bg = "white", dpi = 300)
  cat("ICER-vs-price figure written.\n")
}


# =============================================================================
# CEA gt table — matches the published column layout
# =============================================================================
# ---- formatters ----
# absolute cost: A$3,550 (A$2,913 - A$5,335) Million   [best, rounded to $M]
f_cost_abs <- function(v, lo, hi) {
  m <- function(x) paste0("A$", format(round(x / 1e6), big.mark = ",", scientific = FALSE))
  sprintf("%s\n(%s - %s) Million", m(v), m(lo), m(hi))
}
# absolute QALY: 5,318 (4,315 - 6,447) thousand   [best, rounded to thousand]
f_qaly_abs <- function(v, lo, hi) {
  k <- function(x) format(round(x / 1e3), big.mark = ",", scientific = FALSE)
  sprintf("%s\n(%s - %s) thousand", k(v), k(lo), k(hi))
}
# incremental cost: -A$12.5 (-A$62.2 - A$1.0) Million  [best, 1 dp, signed]
f_cost_inc <- function(v, lo, hi) {
  m <- function(x) sprintf("%sA$%s", ifelse(x < 0, "-", ""),
                           format(round(abs(x) / 1e3, 1), nsmall = 1, big.mark = ","))
  ifelse(is.na(v), "NA", sprintf("%s\n(%s - %s) thousand", m(v), m(lo), m(hi)))
}
# QALYs gained: 375 (28 - 595)   [best, raw QALYs]
f_qaly_inc <- function(v, lo, hi) {
  n <- function(x) format(round(x), big.mark = ",", scientific = FALSE)
  ifelse(is.na(v), "NA", sprintf("%s\n(%s - %s)", n(v), n(lo), n(hi)))
}
# ICER: dominant (dominant - 13,351)  [median; negative -> dominant]
r_icer <- function(x) ifelse(x < 0, "dominant",
                             format(round(x), big.mark = ",", scientific = FALSE))
f_icer <- function(v, lo, hi) ifelse(is.na(v), "NA",
                                     sprintf("%s\n(%s - %s)", r_icer(v), r_icer(lo), r_icer(hi)))

cea_gt <- cea_full %>%
  filter(timeframe == "20y") %>%
  arrange(sensitivity, scenario) %>%
  transmute(
    sensitivity, scenario,
    `Costs, million\n(discounted)`             = f_cost_abs(cost_best, cost_q2.5, cost_q97.5),
    `QALYs, thousand\n(discounted)`            = f_qaly_abs(qaly_best, qaly_q2.5, qaly_q97.5),
    `Incremental costs, million\n(discounted)` = f_cost_inc(inc_cost_best, inc_cost_q2.5, inc_cost_q97.5),
    `QALYs gained\n(discounted)`               = f_qaly_inc(inc_qaly_best, inc_qaly_q2.5, inc_qaly_q97.5),
    ICER                                       = f_icer(inc_cost_best/inc_qaly_best, icer_q2.5, icer_q97.5)
  )

cea_gt %>%
  gt(groupname_col = "sensitivity", rowname_col = "scenario") %>%
  tab_header(title = "Cost-effectiveness across DAA price levels (20-year horizon)") %>%
  cols_align("center", columns = -scenario) %>%
  tab_footnote(footnote = paste0(
    "Costs, QALYs and incremental values are best estimate (95% credible ",
    "interval); ICER is median (95% credible interval). Incremental values ",
    "and ICER are relative to no national program. 'dominant' = the scenario ",
    "both saves cost and gains QALYs.")) %>%
  tab_options(table.font.size = px(11),
              heading.title.font.weight = "bold",
              column_labels.font.weight = "bold",
              row_group.font.weight = "bold") %>%
  gtsave(file = file.path(OutputFig, "HCV_ICER_DAA_sweep.docx"))

cat("CEA gt table written.\n")
cat("\n=== downstream figures + CEA complete ===\n")

####################################################################
# =============================================================================
# ROI across the DAA sweep — updated for the current (per-file) structure
# =============================================================================
# Run after cost_figures_CEA_sweep.R (or standalone — it reloads what it needs).
#
#   ROI = total cost savings / program cost
#   total cost savings = -incremental cost + program cost
#
# WHERE THE NUMBERS COME FROM
#   incremental cost : cap_sets (per-set cumulative discounted CAPPED category
#                      total) at 2041, scenario - no program. Same object the
#                      CEA table uses, so ROI and CEA reconcile.
#   program cost     : NP-lane spend = cost_ab_sc + cost_RNA_sc + cost_POCT_sc
#                      + cost_TreatOther_sc, cumulative discounted to 2030.
#
# WHY PROGRAM COST IS RECOMPUTED (not read from acc$cq_disycum)
#   The generic cq_disycum in cost_categories_sweep.R blanks zeros to NA BEFORE
#   cumsum. For sparse _sc series (zero in early years) that NA propagates
#   through the whole cumulative column. The doc-9 program-cost path na0's
#   BEFORE cumsum, so it must be rebuilt from the raw yearly _sc costs. We
#   reload each Res_flowcost_*.rda once, extract only the four _sc objects, and
#   drop the file — cheap and memory-safe.
# =============================================================================
# =============================================================================
# ROI across the DAA sweep — updated for the current (per-file) structure
# =============================================================================
# Run after cost_figures_CEA_sweep.R (or standalone — it reloads what it needs).
#
#   ROI = total cost savings / program cost
#   total cost savings = -incremental cost + program cost
#
# WHERE THE NUMBERS COME FROM
#   incremental cost : cap_sets (per-set cumulative discounted CAPPED category
#                      total) at 2041, scenario - no program. Same object the
#                      CEA table uses, so ROI and CEA reconcile.
#   program cost     : NP-lane spend = cost_ab_sc + cost_RNA_sc + cost_POCT_sc
#                      + cost_TreatOther_sc, cumulative discounted to 2030.
#
# WHY PROGRAM COST IS RECOMPUTED (not read from acc$cq_disycum)
#   The generic cq_disycum in cost_categories_sweep.R blanks zeros to NA BEFORE
#   cumsum. For sparse _sc series (zero in early years) that NA propagates
#   through the whole cumulative column. The doc-9 program-cost path na0's
#   BEFORE cumsum, so it must be rebuilt from the raw yearly _sc costs. We
#   reload each Res_flowcost_*.rda once, extract only the four _sc objects, and
#   drop the file — cheap and memory-safe.
# =============================================================================


load(file.path(OutputFolder, paste0(project_name, "cost_categories_sweep.rda")))
# -> cost_disydaacap_categories, plan, sens_levels, sens_labels, ...

par_col   <- c("best", paste0("set", seq_len(POC_AU$numberSamples)))
bench_scn <- "(1) No national program"
AUdiscount <- 0.05

sce_level <- c("no_np", "foundational", "succession", "sustained", "scaleup")
sce_label <- c("(1) No national program", "(2) Foundational implementation",
               "(3) Program succession", "(4) Program sustained",
               "(5) Program scale-up")

prog_ind <- c("cost_ab_sc", "cost_RNA_sc", "cost_POCT_sc", "cost_TreatOther_sc")
pc_year  <- 2030                        # program cost measured to end of rollout
inc_year <- POC_AU$simY + 20 - 1        # savings measured over 20 years (2041)

na0 <- function(v) { v[is.na(v)] <- 0; v }


# =============================================================================
# INCREMENTAL COST — per-set capped category total at inc_year (as in the CEA)
# =============================================================================
cost_disydaanocap_categories
cap_sets <- cost_disydaanocap_categories %>%
  filter(year >= 2022) %>%
  group_by(sens_code, scenario, Categories) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(across(all_of(par_col), cumsum)) %>%
  ungroup() %>%
  group_by(sens_code, scenario, year) %>%
  summarise(across(all_of(par_col), ~ sum(.x)), .groups = "drop")   # na.rm=FALSE, matches CEA


# =============================================================================
# PROGRAM COST — recomputed the doc-9 way, one file at a time
# =============================================================================
# program_cost_2030[[sens_code]][[scenario_label]] = per-set numeric vector
#   (element 1 = best, 2..n = sets), cumulative discounted to pc_year.
# =============================================================================
program_cost_2030 <- list()

for (r in seq_len(nrow(plan))) {
  pr <- plan[r, ]; sc <- pr$sens_code
  f  <- file.path(OutputFolder,
                  paste0(project_name, "Res_flowcost_", pr$cost_type, "_", pr$daa_label, ".rda"))
  if (!file.exists(f)) stop("Missing file: ", f)
  
  e <- new.env(); load(f, envir = e)
  RY <- e$Rescost_year_all; rm(e)
  
  for (m in names(RY)) {                       # m = scenario CODE (no_np, ...)
    mlab <- sce_label[match(m, sce_level)]
    if (is.na(mlab)) next
    
    yr  <- RY[[m]][[prog_ind[1]]]$year
    mat <- Reduce(`+`, lapply(prog_ind, function(ind)
      na0(as.matrix(RY[[m]][[ind]][, par_col]))))
    
    id   <- yr - POC_AU$simY
    disc <- ifelse(id >= 0, (1 + AUdiscount)^id, NA)
    mat  <- sweep(mat, 1, disc, "/")           # discount each year
    
    d <- data.frame(year = yr, mat); names(d) <- c("year", par_col)
    d <- d %>% filter(year >= 2022) %>% arrange(year) %>%
      mutate(across(all_of(par_col), ~ cumsum(na0(.x))))   # na0 BEFORE cumsum
    
    row <- d %>% filter(year == pc_year)
    if (nrow(row)) {
      v <- as.numeric(row[1, par_col])
      # A draw with no valid NP-lane spend (missing/all-NA component that na0
      # turned into a propagated 0) is not a real zero-cost program. Convert
      # those back to NA at source so they never enter cumsum-derived summaries
      # or the quantiles. The "best" element (index 1) is left as-is.
      bad <- c(FALSE, !is.finite(v[-1]) | v[-1] == 0)
      v[bad] <- NA_real_
      program_cost_2030[[sc]][[mlab]] <- v
      n_bad <- sum(bad[-1])
      if (n_bad > 0)
        cat(sprintf("      [%s / %s] %d/%d sets NA/zero program cost at %d -> set NA\n",
                    sc, mlab, n_bad, length(v) - 1, pc_year))
    }
  }
  
  rm(RY); gc()
  cat(sprintf("  program cost: %s\n", sc))
}


# =============================================================================
# ROI  =  (-incremental cost + program cost) / program cost, per set
# =============================================================================
# q() summarises the valid draws. NA/zero handling now happens at source
# (program cost is NA for broken draws), so here we only need na.rm.
q <- function(v, p) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(NA_real_)
  quantile(v, p, na.rm = TRUE)
}

roi_rows <- list()

for (sc in sens_levels) {
  cost_at <- cap_sets %>% filter(sens_code == sc, year == inc_year)
  ref <- cost_at %>% filter(scenario == bench_scn)
  if (nrow(ref) == 0) next
  rc <- as.numeric(ref[1, par_col])
  
  for (slab in sce_label[-1]) {
    cc <- cost_at %>% filter(scenario == slab)
    pc <- program_cost_2030[[sc]][[slab]]
    if (nrow(cc) == 0 || is.null(pc)) next
    
    vc  <- as.numeric(cc[1, par_col])
    inc <- vc - rc                    # incremental cost (negative = cheaper)
    sav <- -inc + pc                  # total cost savings
    roi <- sav / pc
    
    # per-set draws (drop the leading "best" element). Program cost is the
    # gating quantity: where it is NA (broken draw), the set is excluded from
    # all three CrIs so savings/ROI/pc stay on the same valid set.
    s_pc  <- pc[-1]
    valid <- is.finite(s_pc)
    s_sav <- (sav[-1])[valid]
    s_roi <- (roi[-1])[valid]
    s_pc  <- s_pc[valid]
    
    roi_rows[[length(roi_rows) + 1]] <- data.frame(
      sens_code = sc, sensitivity = sens_labels[[sc]], scenario = slab,
      n_valid   = length(s_pc),
      sav_best  = sav[1],
      sav_q2.5  = q(s_sav, 0.025),
      sav_q97.5 = q(s_sav, 0.975),
      pc_best   = pc[1],
      pc_q2.5   = q(s_pc, 0.025),
      pc_q97.5  = q(s_pc, 0.975),
      roi_best  = roi[1],
      roi_q2.5  = q(s_roi, 0.025),
      roi_q97.5 = q(s_roi, 0.975),
      stringsAsFactors = FALSE)
  }
}

roi_full <- bind_rows(roi_rows) %>%
  mutate(sensitivity = factor(sensitivity, levels = sens_labels[sens_levels]),
         scenario    = factor(scenario,    levels = sce_label[-1]))

write.csv(roi_full,
          file.path(OutputFolder, paste0(project_name, "DAA_sweep_ROI.csv")),
          row.names = FALSE)


# =============================================================================
# ROI gt table — grouped by sensitivity
# =============================================================================
f_cur <- function(v, lo, hi) {
  m <- function(x) paste0("A$", format(round(x / 1e6), big.mark = ",", scientific = FALSE))
  sprintf("%s\n(%s - %s)", m(v), m(lo), m(hi))
}
f_roi <- function(v, lo, hi) sprintf("%.2f\n(%.2f - %.2f)", v, lo, hi)

roi_gt_df <- roi_full %>%
  arrange(sensitivity, scenario) %>%
  transmute(
    sensitivity, scenario,
    `Total cost savings, million\n(discounted)` = f_cur(sav_best, sav_q2.5, sav_q97.5),
    `Program cost, million\n(discounted)`       = f_cur(pc_best,  pc_q2.5,  pc_q97.5),
    `ROI\n(savings / program cost)`             = f_roi(roi_best, roi_q2.5, roi_q97.5)
  )

roi_gt_df %>%
  gt(groupname_col = "sensitivity", rowname_col = "scenario") %>%
  tab_header(title = "Return on investment across DAA price levels") %>%
  cols_align("center", columns = -scenario) %>%
  tab_footnote(footnote = paste0(
    "Estimate (95% credible interval). ROI = total cost savings / program cost; ",
    "ROI > 1 indicates a positive return. Savings measured over 20 years (to ",
    inc_year, "); program cost cumulative to ", pc_year, ".")) %>%
  tab_options(table.font.size = px(11),
              heading.title.font.weight = "bold",
              column_labels.font.weight = "bold",
              row_group.font.weight = "bold") %>%
  gtsave(file = file.path(OutputFig, "HCV_ROI_DAA_sweep.docx"))

cat("ROI table written.\n=== ROI complete ===\n")



# =============================================================================
# Per-simulation tallies: how many draws are dominant (CEA) and ROI > 1
# =============================================================================
# Standalone. Reproduces the per-set incremental cost / QALY / program cost the
# CEA and ROI tables are built from, then COUNTS draws rather than summarising
# them. Output: one row per (sensitivity x scenario) with counts, valid N, and
# percentages.
#
#   dominant : inc_cost < 0 AND inc_qaly > 0   (strict dominance; the scenario
#              both saves cost and gains QALYs). Counted over all 1000 draws.
#   ROI > 1  : total cost savings / program cost > 1, over the draws with a
#              valid (non-NA) program cost.
#
# Both use 20-year horizon (inc_year) for cost/QALY and program cost to 2030,
# matching the CEA and ROI scripts.
# =============================================================================


AUdiscount <- 0.05
WTP       <- 50000                      # A$ per QALY, for the CE threshold



# =============================================================================
# Per-set capped category total (cost) at inc_year  — same as CEA/ROI
# =============================================================================
cap_sets <- cost_disydaanocap_categories %>%
  filter(year >= 2022) %>%
  group_by(sens_code, scenario, Categories) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(across(all_of(par_col), cumsum)) %>%
  ungroup() %>%
  group_by(sens_code, scenario, year) %>%
  summarise(across(all_of(par_col), ~ sum(.x)), .groups = "drop")

get_qaly <- function(sc) acc$cq_disycum[[sc]]$QALY_compartment


# =============================================================================
# Program cost per set at pc_year — recomputed the ROI way (NA at source)
# =============================================================================
program_cost_2030 <- list()

for (r in seq_len(nrow(plan))) {
  pr <- plan[r, ]; sc <- pr$sens_code
  f  <- file.path(OutputFolder,
                  paste0(project_name, "Res_flowcost_", pr$cost_type, "_", pr$daa_label, ".rda"))
  if (!file.exists(f)) stop("Missing file: ", f)
  
  e <- new.env(); load(f, envir = e); RY <- e$Rescost_year_all; rm(e)
  
  for (m in names(RY)) {
    mlab <- sce_label[match(m, sce_level)]
    if (is.na(mlab)) next
    
    yr  <- RY[[m]][[prog_ind[1]]]$year
    mat <- Reduce(`+`, lapply(prog_ind, function(ind)
      na0(as.matrix(RY[[m]][[ind]][, par_col]))))
    id   <- yr - POC_AU$simY
    disc <- ifelse(id >= 0, (1 + AUdiscount)^id, NA)
    mat  <- sweep(mat, 1, disc, "/")
    
    d <- data.frame(year = yr, mat); names(d) <- c("year", par_col)
    d <- d %>% filter(year >= 2022) %>% arrange(year) %>%
      mutate(across(all_of(par_col), ~ cumsum(na0(.x))))
    
    row <- d %>% filter(year == pc_year)
    if (nrow(row)) {
      v <- as.numeric(row[1, par_col])
      bad <- c(FALSE, !is.finite(v[-1]) | v[-1] == 0)
      v[bad] <- NA_real_
      program_cost_2030[[sc]][[mlab]] <- v
    }
  }
  rm(RY); gc()
}


# =============================================================================
# TALLY per (sensitivity x scenario)
# =============================================================================
tally_rows <- list()

for (sc in sens_levels) {
  cost_at  <- cap_sets %>% filter(sens_code == sc, year == inc_year)
  qaly_all <- get_qaly(sc)
  if (nrow(cost_at) == 0 || is.null(qaly_all)) next
  qaly_at  <- qaly_all %>% filter(year == inc_year)
  
  ref_c <- cost_at %>% filter(scenario == bench_scn)
  ref_q <- qaly_at %>% filter(scenario == bench_scn)
  if (nrow(ref_c) == 0 || nrow(ref_q) == 0) next
  rc <- as.numeric(ref_c[1, set_cols])   # sets only (exclude "best")
  rq <- as.numeric(ref_q[1, set_cols])
  
  for (slab in sce_label[-1]) {
    cc <- cost_at %>% filter(scenario == slab)
    qq <- qaly_at %>% filter(scenario == slab)
    pc <- program_cost_2030[[sc]][[slab]]
    if (nrow(cc) == 0 || nrow(qq) == 0) next
    
    vc <- as.numeric(cc[1, set_cols])
    vq <- as.numeric(qq[1, set_cols])
    
    inc_cost <- vc - rc
    inc_qaly <- vq - rq
    
    # ---- dominant: cheaper AND more effective ----
    dom_ok  <- is.finite(inc_cost) & is.finite(inc_qaly)
    n_dom_denom <- sum(dom_ok)
    n_dominant  <- sum(inc_cost < 0 & inc_qaly > 0, na.rm = TRUE)
    
    # ---- cost-effective at WTP: net monetary benefit > 0 ----
    # NMB = WTP * inc_qaly - inc_cost. Linear, valid in all four quadrants
    # (no ICER ratio, no sign-flip). Dominant draws are a subset of these.
    nmb <- WTP * inc_qaly - inc_cost
    n_ce <- sum(nmb > 0, na.rm = TRUE)
    
    # ---- also useful: cost-saving regardless of QALY ----
    n_cost_saving <- sum(inc_cost < 0, na.rm = TRUE)
    
    # ---- ROI > 1, over draws with valid program cost ----
    if (!is.null(pc)) {
      pc_sets <- pc[-1]                       # sets only
      sav <- -inc_cost + pc_sets
      roi <- sav / pc_sets
      roi_ok    <- is.finite(pc_sets) & is.finite(roi)
      n_roi_denom <- sum(roi_ok)
      n_roi_gt1   <- sum(roi > 1 & roi_ok, na.rm = TRUE)
    } else {
      n_roi_denom <- NA_integer_; n_roi_gt1 <- NA_integer_
    }
    
    tally_rows[[length(tally_rows) + 1]] <- data.frame(
      sens_code   = sc,
      sensitivity = sens_labels[[sc]],
      scenario    = slab,
      # --- probability of cost-effectiveness (CEAC value at WTP = 50,000) ---
      n_ce         = n_ce,
      n_ce_denom   = n_dom_denom,
      prob_ce      = n_ce / n_dom_denom,
      # --- probability ROI > 1 ---
      n_roi_gt1    = n_roi_gt1,
      n_roi_denom  = n_roi_denom,
      prob_roi_gt1 = n_roi_gt1 / n_roi_denom,
      stringsAsFactors = FALSE)
  }
}

library(tidyr)

tally <- bind_rows(tally_rows) %>%
  mutate(sensitivity = factor(sensitivity, levels = sens_labels[sens_levels]),
         scenario    = factor(scenario,    levels = sce_label[-1])) %>%
  arrange(sensitivity, scenario)

# --- tidy long frame: one row per (sens x scenario x metric) -> heatmap-ready
heat_dt <- tally %>%
  select(sens_code, sensitivity, scenario, prob_ce, prob_roi_gt1) %>%
  pivot_longer(c(prob_ce, prob_roi_gt1),
               names_to = "metric", values_to = "probability") %>%
  mutate(metric = recode(metric,
                         prob_ce      = "P(cost-effective)",
                         prob_roi_gt1 = "P(ROI > 1)"))

# --- CSVs ---
write.csv(tally,   file.path(OutputFolder, paste0(project_name, "DAA_sweep_prob.csv")),
          row.names = FALSE)                                   # wide, with counts + denoms
write.csv(heat_dt, file.path(OutputFolder, paste0(project_name, "DAA_sweep_prob_long.csv")),
          row.names = FALSE)                                   # tidy long, for the heatmap

cat("\n=== Probability of CE and ROI>1 (20-year horizon, WTP = 50,000) ===\n")
print(as.data.frame(tally), row.names = FALSE)

cat("\n=== Dominant / ROI>1 tally (20-year horizon) ===\n")
print(as.data.frame(tally), row.names = FALSE)
cat(sprintf("\nSaved: %s\n",
            file.path(OutputFolder, paste0(project_name, "DAA_sweep_tally.csv"))))

library(scales)
# heatmap
Here's the complete block — from the tally loop through the heatmap. I've folded in the best-estimate ROI computation and set the cell text to show probability on top with ROI: x.xx beneath.

r
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# =============================================================================
# TALLY per (sensitivity x scenario) — probabilities + best-estimate ROI
# =============================================================================
tally_rows <- list()

for (sc in sens_levels) {
  cost_at  <- cap_sets %>% filter(sens_code == sc, year == inc_year)
  qaly_all <- get_qaly(sc)
  if (nrow(cost_at) == 0 || is.null(qaly_all)) next
  qaly_at  <- qaly_all %>% filter(year == inc_year)
  
  ref_c <- cost_at %>% filter(scenario == bench_scn)
  ref_q <- qaly_at %>% filter(scenario == bench_scn)
  if (nrow(ref_c) == 0 || nrow(ref_q) == 0) next
  
  rc <- as.numeric(ref_c[1, set_cols])     # sets only (draws 2..n)
  rq <- as.numeric(ref_q[1, set_cols])
  rc_best <- as.numeric(ref_c[1, par_col])[1]   # element 1 = best estimate
  
  for (slab in sce_label[-1]) {
    cc <- cost_at %>% filter(scenario == slab)
    qq <- qaly_at %>% filter(scenario == slab)
    pc <- program_cost_2030[[sc]][[slab]]
    if (nrow(cc) == 0 || nrow(qq) == 0) next
    
    vc <- as.numeric(cc[1, set_cols])
    vq <- as.numeric(qq[1, set_cols])
    
    inc_cost <- vc - rc
    inc_qaly <- vq - rq
    
    # ---- probability of cost-effectiveness (NMB > 0 at WTP) ----
    dom_ok <- is.finite(inc_cost) & is.finite(inc_qaly)
    n_ce_denom <- sum(dom_ok)
    nmb  <- WTP * inc_qaly - inc_cost
    n_ce <- sum(nmb > 0, na.rm = TRUE)
    
    # ---- probability ROI > 1 (over draws with valid program cost) ----
    if (!is.null(pc)) {
      pc_sets <- pc[-1]                       # sets only
      sav     <- -inc_cost + pc_sets
      roi     <- sav / pc_sets
      roi_ok  <- is.finite(pc_sets) & is.finite(roi)
      n_roi_denom <- sum(roi_ok)
      n_roi_gt1   <- sum(roi > 1 & roi_ok, na.rm = TRUE)
    } else {
      n_roi_denom <- NA_integer_; n_roi_gt1 <- NA_integer_
    }
    
    # ---- best-estimate ROI (deterministic, element 1) ----
    vc_best       <- as.numeric(cc[1, par_col])[1]
    inc_cost_best <- vc_best - rc_best
    roi_best <- if (!is.null(pc)) {
      pcb <- pc[1]
      if (is.finite(pcb) && pcb != 0) (-inc_cost_best + pcb) / pcb else NA_real_
    } else NA_real_
    
    tally_rows[[length(tally_rows) + 1]] <- data.frame(
      sens_code    = sc,
      sensitivity  = sens_labels[[sc]],
      scenario     = slab,
      n_ce         = n_ce,
      n_ce_denom   = n_ce_denom,
      prob_ce      = n_ce / n_ce_denom,
      n_roi_gt1    = n_roi_gt1,
      n_roi_denom  = n_roi_denom,
      prob_roi_gt1 = n_roi_gt1 / n_roi_denom,
      roi_best     = roi_best,
      stringsAsFactors = FALSE)
  }
}

tally <- bind_rows(tally_rows) %>%
  mutate(sensitivity = factor(sensitivity, levels = sens_labels[sens_levels]),
         scenario    = factor(scenario,    levels = sce_label[-1])) %>%
  arrange(sensitivity, scenario)

# --- CSVs ---
write.csv(tally,
          file.path(OutputFolder, paste0(project_name, "DAA_sweep_prob.csv")),
          row.names = FALSE)

heat_dt <- tally %>%
  select(sens_code, sensitivity, scenario, prob_ce, prob_roi_gt1, roi_best) %>%
  pivot_longer(c(prob_ce, prob_roi_gt1),
               names_to = "metric", values_to = "probability") %>%
  mutate(metric = recode(metric,
                         prob_ce      = "P(cost-effective)",
                         prob_roi_gt1 = "P(ROI > 1)"))
write.csv(heat_dt,
          file.path(OutputFolder, paste0(project_name, "DAA_sweep_prob_long.csv")),
          row.names = FALSE)

cat("\n=== Probability of CE / ROI>1 + best-estimate ROI (20-year horizon) ===\n")
print(as.data.frame(tally), row.names = FALSE)


# =============================================================================
# HEATMAP: fill = P(ROI > 1), high = good (reversed ramp); text = ROI: x.xx
# =============================================================================
library(ggplot2); library(scales)

heat_roi <- ggplot(tally, aes(x = scenario, y = sensitivity, fill = prob_roi_gt1)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%\nROI: %.2f", 100 * prob_roi_gt1, roi_best),
                colour = prob_roi_gt1 > 0.55),
            size = 3.6, fontface = "bold", lineheight = 1.0, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
  scale_fill_viridis_c(
    option    = "mako",
    direction = 1,
    limits    = c(0, 1),
    breaks    = c(0, 0.25, 0.5, 0.75, 1),
    labels    = percent,
    name      = "P(ROI > 1)"
  ) +
  scale_x_discrete(expand = c(0, 0), position = "top") +   # scenario labels on top, above the grid
  scale_y_discrete(expand = c(0, 0), limits = rev(levels(tally$sensitivity))) +
  labs(x = NULL, y = NULL,
       title = "Probability ROI > 1 (fill) and best-estimate ROI (value)") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x.top = element_text(angle = 30, hjust = 0, face = "bold", size = 11),
    axis.text.y     = element_text(face = "bold", size = 11),
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 13),
    legend.key.height = unit(1.6, "cm"),
    plot.margin     = margin(10, 40, 10, 10)   # right margin for angled top labels
  )

heat_roi

check <- tally %>% filter(sens_code %in% c("fixednvariable_DAA030","fixednvariable_DAA090"))
