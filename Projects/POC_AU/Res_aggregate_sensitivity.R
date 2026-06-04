####Res_aggregate update#####
#### this script work with the sensitivty analysis on the unit cost of DAA 

# this script imports the datasets tidy up in Resultssummary.R 
# {Res_dt_scenarios.rda} and {Rescost_dt_[scn]_[cost_type].rda}
# we aggregate the number to annual numbers and apply cost discount in this script. 
# three datasets we work on in this script: {Num_box}, {Resflow_dt}, {Rescost_dt}
gc()
rm(list = ls())
tic <- proc.time()

project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Projects/Simplified-HCV-testing-model")

data_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                   project_name, sep = "")
# Load useful libraries

library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")
library("openxlsx")

Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
RDAFolder <- file.path(data_path, "02. Output" )
OutputFolder <- file.path(codefun_path, "Projects/POC_AU/Output")
OutputFig <- file.path(codefun_path, "Projects/POC_AU/Figs")
# project specific code path 
Proj_code <- file.path(codefun_path, paste0("projects/", project_name))

load(file.path(RDAFolder, paste0(project_name, ".rda")))


# load rda into list of list 
rda2list <- function(file) {
  e <- new.env()
  load(file, envir = e)
  as.list(e)
}

sce_name <- c("no_np","foundational", "succession",
              "sustained", "scaleup")

cost_types <- c("fixednvariable", "total",  "DAAcost_reduchalf")
source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Proj_code, "/model_timestep.R")) 
AUdiscount <- 0.05

cap <- 200000000
unitC_fibroscan_old <- 62.52                       # bundled cost previously in ctau_ag/ctau_poct
unitC_fibroscan_SOC <- c(C = 73.44, P = 82.30)     # SOC: community / prison
unitC_fibroscan_POC <- c(C = 75.59, P = 83.77)     # POC: community / prison
unitC_eta_other_SOC <- c(C = 1341.14, P = 1063.77)   # community / prison
unitC_eta_other_POC <- c(C = 1342.61, P = 1059.82)
endY <- 100
par_col <- c("best", paste0("set", seq(1, POC_AU$numberSamples,1)))

# NaN-safe element-wise: replace NaN with 0 so NaN + x does not poison the sum
z <- function(d) { m <- as.matrix(d); m[is.nan(m)] <- 0; m }

# ---------------------------------------------------------------------------
# Load epi flows ONCE: Res_dt_scenarios.rda has no cost_type dependence.
# It contains Num_box, Resflow_dt, Resflow_sc_dt -- each keyed by scenario.
# ---------------------------------------------------------------------------
epi_env <- new.env()
load(file.path(OutputFolder, paste0(project_name, "Res_dt_scenarios.rda")),
     envir = epi_env)

for(cost_type in cost_types){
  
  # --- Reassemble Res_dt into the nested shape the rest of the script expects:
  #     Res_dt[[scn]]$Num_box[[1]], $Resflow_dt[[1]], $Resflow_sc_dt[[1]],
  #     $Rescost_dt[[1]]
  Res_dt <- list()
  for(scn in sce_name){
    cost_env <- new.env()
    load(file.path(OutputFolder,
                   paste0(project_name, "Rescost_dt_", scn, "_", cost_type, ".rda")),
         envir = cost_env)
    # cost_env$Rescost_dt is a one-element list keyed by scn
    
    Res_dt[[scn]] <- list(
      Num_box       = list(epi_env$Num_box[[scn]]),
      Resflow_dt    = list(epi_env$Resflow_dt[[scn]]),
      Resflow_sc_dt = list(epi_env$Resflow_sc_dt[[scn]]),
      Rescost_dt    = list(cost_env$Rescost_dt[[scn]])
    )
    rm(cost_env)
  }
  
  Res_numbox <- list()
  Resflow_year_pop <- list()
  
  for(i in names(Res_dt)){ 
    Res_numbox[[i]] <- Res_dt[[i]]$Num_box[[1]]
  }
  
  save(Res_numbox,
       file = file.path(OutputFolder,
                        paste0(project_name,"Res_numbox_",cost_type,".rda")))
  
  for(i in names(Res_dt)){
    Resflow_year_pop[[i]] <- Res_dt[[i]]$Resflow_dt[[1]]
    for(indic in names(Resflow_year_pop[[1]])){
      Resflow_year_pop[[i]][[indic]] <- Resflow_year_pop[[i]][[indic]]%>%
        as_tibble()%>%ungroup()%>%arrange(population, year)%>%
        group_by(year, population)%>%
        summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
        arrange(year)
    }
  }
  Resflow_year_all <- list()
  for(i in names(Resflow_year_pop)){
    for(indic in names(Resflow_year_pop[[1]])){
      
      Resflow_year_all[[i]][[indic]] <- Resflow_year_pop[[i]][[indic]]%>%
        ungroup()%>%group_by(year)%>%
        summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
        arrange(year)
      
    }
  }
  
  # scenarios 
  Resflow_sc_year_pop <- list()
  
  for(n in names(Res_dt)){ 
    Resflow_sc_year_pop[[n]] <- Res_dt[[n]]$Resflow_sc_dt[[1]]
  }
  for(n in names(Res_dt)){
    for(indic in names(Resflow_sc_year_pop[[1]])){
      Resflow_sc_year_pop[[n]][[indic]] <- Resflow_sc_year_pop[[n]][[indic]]%>%
        as_tibble()%>%ungroup()%>%arrange(population, year)%>%
        group_by(year, population)%>%
        summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
        arrange(year)
    }
  }
  
  
  
  Resflow_sc_year_all <- list()
  for(i in names(Resflow_sc_year_pop)){
    for(indic in names(Resflow_sc_year_pop[[1]])){
      
      Resflow_sc_year_all[[i]][[indic]] <- Resflow_sc_year_pop[[i]][[indic]]%>%
        ungroup()%>%group_by(year)%>%
        summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
        arrange(year)
      
    }
  }
  
  RescostDAA <- list()
  RescostDAA <- lapply(names(Res_dt), function(x) Res_dt[[x]]$Rescost_dt[[1]])
  names(RescostDAA) <- names(Res_dt)
  
  for(i in names(Res_dt)){
    for(indic in names(RescostDAA[[i]])){
      RescostDAA[[i]][[indic]] <- RescostDAA[[i]][[indic]]%>%
        as_tibble()%>%ungroup()%>%arrange(timestep,population)%>%
        mutate(year = c(rep(seq(1, endY - 2, 1), each = POC_AU$npops*1/POC_AU$timestep),
                        c(rep(endY - 1, POC_AU$npops*(1/POC_AU$timestep)))))%>%
        group_by(year, population)%>%
        summarise(across(c(par_col),~ sum(.x, na.rm = TRUE)))%>%
        arrange(year)
    }
  }
  
  cost_TreatDAA <- list()
  cost_RetreatDAA <- list()
  cost_TreatDAA_sc <- list()
  cost_TreatOther <- list()
  cost_RetreatOther <- list()
  cost_TreatOther_sc <- list()
  cost_totalDAA <- list()
  index_col <- list()
  #### applied different unit cost of DAA to different cost_type sensitivity analysis 
  unit_costs <- list(
    "fixednvariable" = c(DAA = 35956.37, secline = 44613.66),
    "total" = c(DAA = 17978.18, secline = 22306.83),
    "DAAcost_reduchalf" = c(DAA = 17978.18, secline = 22306.83)
  )
  idx <- match(cost_type, cost_types)
  unitC_DAA <- unit_costs[[idx]]["DAA"]
  unitC_secline_DAA <- unit_costs[[idx]]["secline"]
  
  # Verify (optional)
  message(paste0("Cost type: ", cost_type, " | DAA: ", unitC_DAA, " | Secline: ", unitC_secline_DAA))
  
  # Apply costs
  for(i in names(RescostDAA)){
    index_col[[i]] <- cbind.data.frame(
      year = RescostDAA[[i]]$cost_Treatment$year,
      population = RescostDAA[[i]]$cost_Treatment$population
    )
    
    cost_TreatDAA[[i]] <- cbind(
      as.data.frame(index_col[[i]]), 
      Resflow_year_pop[[i]]$Treatment[, c(par_col)] * unitC_DAA
    )
    
    cost_RetreatDAA[[i]] <- cbind(
      as.data.frame(index_col[[i]]), 
      Resflow_year_pop[[i]]$Retreat[, c(par_col)] * unitC_secline_DAA
    )
    
    cost_TreatDAA_sc[[i]] <- cbind(
      as.data.frame(index_col[[i]]), 
      Resflow_sc_year_pop[[i]]$Treatment_sc[, c(par_col)] * unitC_DAA
    )
  }
  index_col <- list()
  for(i in names(RescostDAA)){
    
    index_col[[i]] <- cbind.data.frame(year = RescostDAA[[i]]$cost_Treatment$year, 
                                       population = RescostDAA[[i]]$cost_Treatment$population)
    
    pop_prefix_eta <- substr(as.character(index_col[[i]]$population), 1, 1)
    SOC_eta_other  <- unitC_eta_other_SOC[pop_prefix_eta]
    POC_eta_other  <- unitC_eta_other_POC[pop_prefix_eta]
    
    # Volumes (year x population, already aggregated)
    n_treat    <- as.matrix(Resflow_year_pop[[i]]$Treatment[,    par_col])
    n_treat_sc <- as.matrix(Resflow_sc_year_pop[[i]]$Treatment_sc[, par_col])
    n_retreat  <- as.matrix(Resflow_year_pop[[i]]$Retreat[,     par_col])
    
    # cost_TreatOther = SOC_unit x Treatment + POC_unit x Treatment_sc
    # (sq has Treatment_sc = 0, so collapses to SOC_unit x Treatment cleanly)
    cost_TreatOther[[i]] <- cbind(
      as.data.frame(index_col[[i]]),
      as.data.frame(SOC_eta_other * n_treat + POC_eta_other * n_treat_sc)
    )
    
    # cost_RetreatOther = SOC_unit x Retreat (no _sc lane, SOC unit only per spec)
    cost_RetreatOther[[i]] <- cbind(
      as.data.frame(index_col[[i]]),
      as.data.frame(SOC_eta_other * n_retreat)
    )
    
    # cost_TreatOther_sc = POC_unit x Treatment_sc only
    cost_TreatOther_sc[[i]] <- cbind(
      as.data.frame(index_col[[i]]),
      as.data.frame(POC_eta_other * n_treat_sc)
    )
    
    cost_totalDAA[[i]] <- cbind(as.data.frame(index_col[[i]]), 
                                (cost_TreatDAA[[i]][, c(par_col)] + 
                                   cost_RetreatDAA[[i]][, c(par_col)] + 
                                   cost_TreatDAA_sc[[i]][, c(par_col)]))
  }  
  
  cost_fibroscan    <- list()
  cost_fibroscan_sc <- list()
  
  for(i in names(RescostDAA)){
    # Population-specific unit costs (C* = community, P* = prison)
    pop_vec    <- as.character(RescostDAA[[i]]$cost_ab$population)
    pop_prefix <- substr(pop_vec, 1, 1)
    SOC_unit   <- unitC_fibroscan_SOC[pop_prefix]
    POC_unit   <- unitC_fibroscan_POC[pop_prefix]
    
    # Testing volumes (already aggregated to year x population)
    n_RNA_total  <- as.matrix(Resflow_year_pop[[i]]$Testing_RNA[,     par_col])
    n_RNA_sc     <- as.matrix(Resflow_year_pop[[i]]$Testing_RNA_sc[,  par_col])
    n_POCT_total <- as.matrix(Resflow_year_pop[[i]]$Testing_POCT[,    par_col])
    n_POCT_sc    <- as.matrix(Resflow_year_pop[[i]]$Testing_POCT_sc[, par_col])
    
    # 1. Strip bundled $62.52 fibroscan from existing cost columns
    RescostDAA[[i]]$cost_RNA[,    par_col] <-
      RescostDAA[[i]]$cost_RNA[,    par_col] - (n_RNA_total  + n_RNA_sc)  * unitC_fibroscan_old
    RescostDAA[[i]]$cost_RNA_sc[, par_col] <-
      RescostDAA[[i]]$cost_RNA_sc[, par_col] -  n_RNA_sc                  * unitC_fibroscan_old
    RescostDAA[[i]]$cost_POCT[,    par_col] <-
      RescostDAA[[i]]$cost_POCT[,    par_col] - (n_POCT_total + n_POCT_sc) * unitC_fibroscan_old
    RescostDAA[[i]]$cost_POCT_sc[, par_col] <-
      RescostDAA[[i]]$cost_POCT_sc[, par_col] -  n_POCT_sc                 * unitC_fibroscan_old
    
    # 2. Build cost_fibroscan: SOC_unit x total + POC_unit x sc  (RNA + POCT combined)
    cost_fibroscan[[i]] <- cbind(
      as.data.frame(index_col[[i]]),
      as.data.frame(
        (n_RNA_total + n_POCT_total) * SOC_unit +
          (n_RNA_sc    + n_POCT_sc)    * POC_unit
      )
    )
    
    # 3. Build cost_fibroscan_sc: POC_unit x sc only
    cost_fibroscan_sc[[i]] <- cbind(
      as.data.frame(index_col[[i]]),
      as.data.frame((n_RNA_sc + n_POCT_sc) * POC_unit)
    )
  }
  
  for(i in names(RescostDAA)){
    RescostDAA[[i]] <- append(RescostDAA[[i]],
                              list("cost_TreatDAA"      = cost_TreatDAA[[i]],
                                   "cost_RetreatDAA"    = cost_RetreatDAA[[i]],
                                   "cost_TreatDAA_sc"   = cost_TreatDAA_sc[[i]],
                                   "cost_TreatOther"    = cost_TreatOther[[i]],
                                   "cost_RetreatOther"  = cost_RetreatOther[[i]],
                                   "cost_TreatOther_sc" = cost_TreatOther_sc[[i]],
                                   "cost_totalDAA"      = cost_totalDAA[[i]],
                                   "cost_fibroscan"     = cost_fibroscan[[i]],       # NEW
                                   "cost_fibroscan_sc"  = cost_fibroscan_sc[[i]])    # NEW
    )
  }
  
  RescostDAA_totalpop <- list()
  
  for(i in names(RescostDAA)){
    for(indic in names(RescostDAA[[1]])){
      RescostDAA_totalpop[[i]][[indic]] <- RescostDAA[[i]][[indic]]%>%as_tibble()%>%
        ungroup()%>%arrange(year)%>%group_by(year)%>%
        summarise(across(c(par_col),~ sum(.x, na.rm = TRUE)))%>%ungroup()
    }
    
  }
  
  RescostDAA_sc <- list()
  for(i in names(RescostDAA)){
    
    RescostDAA_sc[[i]] <- list("cost_TreatDAA"     = cost_TreatDAA_sc[[i]],
                               "cost_TreatOther"   = cost_TreatOther_sc[[i]],
                               "cost_fibroscan_sc" = cost_fibroscan_sc[[i]],   # NEW
                               "cost_totalDAA"     = cost_TreatDAA_sc[[i]])
  }
  
  RescostDAA_sc_totalpop <- list()
  
  for(i in names(RescostDAA_sc)){
    for(indic in names(RescostDAA_sc[[1]])){ 
      RescostDAA_sc_totalpop[[i]][[indic]] <- RescostDAA_sc[[i]][[indic]]%>%ungroup()%>%
        arrange(year)%>%group_by(year)%>%
        summarise(across(c(par_col),~ sum(.x, na.rm = TRUE)))%>%ungroup()
    }
  }
  
  
  # remove the original cost_Treatment and cost_Retreat since we breakdown the cost of treatment and retreat 
  # into DAA,and nonDAA cost.
  
  Rescost_DAA_dt <- RescostDAA
  
  for(i in names(RescostDAA)){ 
    Rescost_DAA_dt[[i]]$cost_Treatment <- NULL
    Rescost_DAA_dt[[i]]$cost_Retreat <- NULL
    
  }
  
  
  Rescost_DAA_sc_dt <- list()
  
  for(i in names(RescostDAA)){ 
    Rescost_DAA_sc_dt[[i]] <- list("cost_ab_sc"        = RescostDAA[[i]]$cost_ab,
                                   "cost_RNA_sc"       = RescostDAA[[i]]$cost_RNA_sc,
                                   "cost_POCT_sc"      = RescostDAA[[i]]$cost_POCT_sc,
                                   "cost_TreatDAA"     = cost_TreatDAA_sc[[i]],
                                   "cost_TreatOther"   = cost_TreatOther_sc[[i]],
                                   "cost_fibroscan_sc" = cost_fibroscan_sc[[i]],   # NEW
                                   "cost_totalDAA"     = cost_TreatDAA_sc[[i]])
    
    
  }
  
  # cost_total DAA, exculdeDAA replaced cost_Treatment & cost_Retreat 
  # valide with cost_total, 
  # calculating total_cost_cap  
  # extract the columns to replace in cost_treatment& cost_retreat 
  # "cost_totalDAA_y": "cost_totalDAAcap_discum"
  Rescost_year <- list()
  
  for(i in names(Rescost_DAA_dt)){ 
    for(indic in names(Rescost_DAA_dt[[1]])){ 
      Rescost_year[[i]][[indic]] <- Rescost_DAA_dt[[i]][[indic]]%>%
        as_tibble()%>%mutate(year = year + POC_AU$cabY - 1,
                             id = year - POC_AU$simY, 
                             discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))
      
      
    }
    
  }
  
  
  Rescost_sc_year <- list() 
  
  for(i in names(Rescost_DAA_sc_dt)){ 
    
    for(indic in names(Rescost_DAA_sc_dt[[1]])){ 
      Rescost_sc_year[[i]][[indic]] <- Rescost_DAA_sc_dt[[i]][[indic]]%>%
        as_tibble()%>%mutate(year = year + POC_AU$cabY - 1,
                             id = year - POC_AU$simY, 
                             discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))
    }
  }
  
  
  # cumulative yearly cost by population  & discount value 
  Rescost_yearcum_pop <- list() 
  
  for(i in names(Rescost_year)){ 
    
    for(indic in names(Rescost_year[[1]])){ 
      Rescost_yearcum_pop[[i]][[indic]] <- Rescost_year[[i]][[indic]]%>%
        as_tibble()%>%arrange(population)%>%group_by(population)%>%
        mutate(across(c(par_col), cumsum, .names = "{col}"))
      
    }
    
  }
  
  Rescost_yearcum_sc_pop <- list()
  
  for(i in names(Rescost_sc_year)){ 
    for(indic in names(Rescost_sc_year[[1]])){
      # FIX: was Rescost_yearcum_sc_pop[[i]] <- ... which overwrote every
      # indicator, leaving only the last one. Now indexed [[i]][[indic]].
      Rescost_yearcum_sc_pop[[i]][[indic]] <- Rescost_sc_year[[i]][[indic]]%>%
        as_tibble()%>%arrange(population)%>%group_by(population)%>%
        mutate(across(c(par_col), cumsum, .names = "{col}"))
      
    }
  }
  
  Rescost_year_all <- list()
  
  for(i in names(Rescost_year)){ 
    
    for(indic in names(Rescost_year[[1]])){ 
      Rescost_year_all[[i]][[indic]] <- Rescost_year[[i]][[indic]]%>%
        as.data.frame()%>%
        ungroup()%>%
        group_by(year)%>%
        summarise(across(c(par_col),~ sum(.x, na.rm = TRUE)))%>%ungroup()
    }
  }
  
  
  for(i in names(Rescost_year)){ 
    Rescost_year_all[[i]][["cost_totalDAA_Cap"]] <- 
      Rescost_year_all[[i]][["cost_totalDAA"]]%>%
      mutate(across(c(par_col), ~ ifelse(.>=cap, cap, . ),.names = "{col}"))
    
    # NaN-safe: z() turns NaN into 0 in each component before summing, so a
    # NaN in any one component no longer poisons the whole cell.
    Rescost_year_all[[i]][["cost_total_Cap"]] <-
      cbind(year = Rescost_year_all[[i]]$cost_totalDAA$year,
            as.data.frame(z(Rescost_year_all[[i]][["cost_compartment"]][,    c(par_col)]) +
                            z(Rescost_year_all[[i]][["cost_ab"]][,           c(par_col)]) +
                            z(Rescost_year_all[[i]][["cost_RNA"]][,          c(par_col)]) +
                            z(Rescost_year_all[[i]][["cost_POCT"]][,         c(par_col)]) +
                            z(Rescost_year_all[[i]][["cost_fibroscan"]][,    c(par_col)]) +   # NEW
                            z(Rescost_year_all[[i]][["cost_totalDAA_Cap"]][, c(par_col)]) +
                            z(Rescost_year_all[[i]][["cost_TreatOther"]][,   c(par_col)]) +
                            z(Rescost_year_all[[i]][["cost_RetreatOther"]][, c(par_col)]) +
                            z(Rescost_year_all[[i]][["cost_Cured"]][,        c(par_col)])))
  }
  
  Rescost_year_sc_all <- list()
  
  for(i in names(Rescost_sc_year)){ 
    for(indic in names(Rescost_sc_year[[1]])){ 
      
      Rescost_year_sc_all[[i]][[indic]] <- Rescost_sc_year[[i]][[indic]]%>%
        as.data.frame()%>%
        ungroup()%>%
        group_by(year)%>%
        summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%ungroup()
      
    }
  }
  
  Rescost_disyear_all <- list()
  
  for(i in names(Rescost_year_all)){ 
    for(indic in names(Rescost_year_all[[1]])){ 
      
      Rescost_disyear_all[[i]][[indic]] <- Rescost_year_all[[i]][[indic]]%>%
        as.data.frame()%>%
        ungroup()%>%
        mutate(id = year - POC_AU$simY, 
               discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))%>%
        mutate(across(c(par_col), ~./discount,
                      .names = "{col}"))
      
    }
  }
  
  Rescost_disyear_sc_all <- list()
  
  for(i in names(Rescost_sc_year)){ 
    for(indic in names(Rescost_sc_year[[1]])){ 
      
      Rescost_disyear_sc_all[[i]][[indic]] <- Rescost_year_sc_all[[i]][[indic]]%>%
        as.data.frame()%>%
        ungroup()%>%
        mutate(id = year - POC_AU$simY, 
               discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))%>%
        mutate(across(c(par_col), ~./discount,
                      .names = "{col}"))
      
    }
  }
  
  
  # cumulative cost each year for overall pops
  Rescost_yearcum_all <- list()
  
  for(i in names(Rescost_year_all)){ 
    for(indic in names(Rescost_year_all[[1]])){
      Rescost_yearcum_all[[i]][[indic]] <- Rescost_year_all[[i]][[indic]]%>%
        ungroup()%>%
        mutate(across(c(par_col), cumsum,
                      .names = "{col}"))%>%
        mutate(id = year - POC_AU$simY, 
               discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))
    }
  }
  
  # cumulative cost each year for overall pops
  Rescost_yearcum_sc_all <- list()
  
  for(i in names(Rescost_year_sc_all)){ 
    for(indic in names(Rescost_year_sc_all[[1]])){
      Rescost_yearcum_sc_all[[i]][[indic]] <- Rescost_year_sc_all[[i]][[indic]]%>%
        ungroup()%>%
        mutate(across(c(par_col), cumsum,
                      .names = "{col}"))%>%
        mutate(id = year - POC_AU$simY, 
               discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))
      
    }
  }
  
  # discount cumulative cost each year for overall pops
  
  Rescost_discum_all <- list()
  
  for(i in names(Rescost_year_all)){ 
    for(indic in names(Rescost_year_all[[1]])){
      Rescost_discum_all[[i]][[indic]] <- Rescost_yearcum_all[[i]][[indic]]%>%
        ungroup()%>%
        mutate(across(c(par_col), ~./discount,
                      .names = "{col}"))
    }
  }
  
  # discount cumulative cost each year for overall pops
  Rescost_discum_sc_all <- list()
  
  for(i in names(Rescost_year_sc_all)){ 
    for(indic in names(Rescost_year_sc_all[[1]])){
      Rescost_discum_sc_all[[i]][[indic]] <- Rescost_yearcum_sc_all[[i]][[indic]]%>%
        ungroup()%>%
        mutate(across(c(par_col), ~ ./discount,
                      .names = "{col}"))
      
    }
  }
  
  # save rda files 
  save(Resflow_year_pop, Resflow_year_all, 
       Resflow_sc_year_pop, Resflow_sc_year_all,
       RescostDAA, RescostDAA_totalpop, 
       Rescost_year, Rescost_sc_year,
       Rescost_yearcum_pop, Rescost_yearcum_sc_pop, 
       Rescost_year_all, Rescost_year_sc_all, 
       Rescost_disyear_all, Rescost_disyear_sc_all,
       Rescost_yearcum_all, Rescost_yearcum_sc_all,
       Rescost_discum_all, Rescost_discum_sc_all,
       AUdiscount,
       unitC_DAA,
       unitC_secline_DAA,
       unitC_fibroscan_old,
       unitC_fibroscan_SOC,
       unitC_fibroscan_POC,
       cap,
       file = file.path(OutputFolder,
                        paste0(project_name,"Res_flowcost_", cost_type ,".rda"))) 
  
  gc()
}
