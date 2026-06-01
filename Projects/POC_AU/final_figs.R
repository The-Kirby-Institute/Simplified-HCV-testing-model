#This script using the Rda file {Res_flowcost.rda} from {Res_aggregate.R} to 
# generate the Tables and Figures 
# Tables and figures are generated based on overall pops, all subpops, and settings(community and prisons)
# the indicators include epi indicators and cost indicators 
# epi indicators include 
# 1. number of infections, number of infections averted 
# 2. number of advanced liver disease and averted 
#    (a) number of dc, (b) number of hcc, (c) number of lt and its averted cases 
# 3. number of HCV deaths and its averted cases 
# 4. Number of people living with HCV: excluding "s" and "cured" 
# 5. Number of people living with HCV get tested
# 6. Number of people living with diagnosed HCV
# 7. Number of reinfection 
# cost 
gc()
rm(list = ls())
gc()
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
library("ggplot2")
library("viridis") 
library("openxlsx")
library("writexl")
urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
Rcode <- file.path(codefun_path, "03. Code")
DataFolder <- file.path(data_path, "01. DATA/model input" )
RDAFolder <- file.path(data_path, "02. Output" )
OutputFolder <- file.path(codefun_path, "Projects/POC_AU/Output")
OutputFig <- file.path(codefun_path, "Projects/POC_AU/Figs")
OutputFig_y_cum_avert <- file.path(OutputFig, "y_cum_avert")
# Create directory if it doesn't exist
if (!dir.exists(OutputFig)) {
  dir.create(OutputFig, recursive = TRUE)
}
if (!dir.exists(OutputFig_y_cum_avert)) {
  dir.create(OutputFig_y_cum_avert, recursive = TRUE)
}
load(file.path(RDAFolder, paste0(project_name, ".rda")))
source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(codefun_path, "/Projects/POC_AU/figtable_code.R"))

endY <- 100
year_obs <- c(POC_AU$simY  +5 - 1  , POC_AU$simY + 10 - 1, POC_AU$simY + 20 - 1)
par_col <- c("best", paste0("set", seq(1,1000,1)))

cost_types <- c("fixednvariable", "total",  "DAAcost_reduchalf")
sce_level <- c("no_np","foundational", "succession",
               "sustained", "scaleup")
sce_label <- c("(1) No national program", "(2) Foundational implementation", 
               "(3) Program succession", "(4) Program sustained", 
               "(5) Program scale-up")
cost_types <- c("fixednvariable", "total", "DAAcost_reduchalf")

load(file.path(OutputFolder, paste0(project_name, "Res_flowcost_",cost_types[1],".rda")))
load(file.path(OutputFolder, paste0(project_name, "Res_numbox_",cost_types[1],".rda")))

for(i in names(Resflow_year_all)){ 
  for(n in names(Resflow_year_all[[1]])){ 
    Resflow_year_all[[i]][[n]] <- Resflow_year_all[[i]][[n]]%>%ungroup()
    Resflow_year_all[[i]][[n]][is.na(Resflow_year_all[[i]][[n]])] <- 0
    Resflow_year_pop[[i]][[n]] <-  Resflow_year_pop[[i]][[n]]%>%ungroup()
    Resflow_year_pop[[i]][[n]][is.na(Resflow_year_pop[[i]][[n]])]  <- 0 
    
    Resflow_year_all[[i]][[n]] <- Resflow_year_all[[i]][[n]]%>%
      mutate(year = ifelse(is.na(year), POC_AU$cabY, year +POC_AU$cabY ))
    
    Resflow_year_pop[[i]][[n]] <- Resflow_year_pop[[i]][[n]]%>%
      mutate(year = ifelse(is.na(year), POC_AU$cabY, year +POC_AU$cabY ))
    
  }
}

for(i in names(Resflow_year_all)){ 
  
  Resflow_year_all[[i]][["Tot_Treatment"]] <- 
    cbind(year = Resflow_year_all[[i]][["Treatment"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Treatment"]][, par_col] + 
                             Resflow_year_all[[i]][["Retreat"]][, par_col] + 
                             Resflow_year_all[[i]][["Treatment_sc"]][, par_col])
    )
  
  Resflow_year_all[[i]][["Tot_Testing_ab"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_ab"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Testing_ab"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_ab_neg"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_ab_sc"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_ab_sc_neg"]][, par_col])
    )
  
  Resflow_year_all[[i]][["Tot_Testing_RNA"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_RNA"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Testing_RNA"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_RNA_neg"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_RNA_sc"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_RNA_sc_neg"]][, par_col])
    )
  
  Resflow_year_all[[i]][["Tot_Testing_POCT"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_POCT"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Testing_POCT"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_POCT_neg"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_POCT_sc"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_POCT_sc_neg"]][, par_col])
    )
}

for(i in names(Resflow_year_all)){
  Resflow_year_pop[[i]][["Tot_Treatment"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Treatment"]]$year, 
          population = Resflow_year_pop[[i]][["Treatment"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Treatment"]][, par_col] + 
                             Resflow_year_pop[[i]][["Retreat"]][, par_col] + 
                             Resflow_year_pop[[i]][["Treatment_sc"]][, par_col])
    )
  
  Resflow_year_pop[[i]][["Tot_Testing_ab"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_ab"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_ab"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Testing_ab"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_ab_neg"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_ab_sc"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_ab_sc_neg"]][, par_col])
    )
  
  
  Resflow_year_pop[[i]][["Tot_Testing_RNA"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_RNA"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_RNA"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Testing_RNA"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_RNA_neg"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_RNA_sc"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_RNA_sc_neg"]][, par_col])
    )
  
  Resflow_year_pop[[i]][["Tot_Testing_POCT"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_POCT"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_POCT"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Testing_POCT"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_POCT_neg"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_POCT_sc"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_POCT_sc_neg"]][, par_col])
    )
  
}
Resflow_year_setting <- list()

for(i in names(Resflow_year_pop)){ 
  for(indic in names(Resflow_year_pop[[1]])){ 
    Resflow_year_setting[[i]][[indic]] <- 
      Resflow_year_pop[[i]][[indic]]%>%
      mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), 
                              "commu", "prisons"))%>%
      group_by(year, setting)%>%
      summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
      ungroup()%>%
      mutate(population = setting)%>%
      select(-setting)%>%
      mutate(population = factor(population, 
                                 levels = c("commu", "prisons"), 
                                 labels = c("Community", "Prisons")))%>%
      select(year, population, par_col)
    
    
    
  }
}
for(i in names(Resflow_year_all)){ 
  for(n in names(Resflow_year_all[[1]])){ 
    Resflow_year_pop[[i]][[n]][Resflow_year_pop[[i]][[n]] ==0]  <- NA 
    Resflow_year_setting[[i]][[n]][Resflow_year_setting[[i]][[n]] == 0]  <- NA  
    Resflow_year_all[[i]][[n]][Resflow_year_all[[i]][[n]] == 0]  <- NA  
  }
}


Resflow_year_setting_range <- list()
Resflow_year_pop_range <- list()
Resflow_year_all_range <- list()
for(i in names(Resflow_year_pop)){ 
  for(indic in names(Resflow_year_pop[[1]])){ 
    Resflow_year_setting_range[[i]][[indic]] <- 
      Resflow_year_setting[[i]][[indic]]%>%ungroup()%>%
      popResults_range(POC_AU, ., Population = c("Community", "Prisons"),
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
    
    Resflow_year_pop_range[[i]][[indic]] <- 
      Resflow_year_pop[[i]][[indic]]%>%ungroup()%>%
      popResults_range(POC_AU, ., Population = POC_AU$popNames,
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
    
    Resflow_year_all_range[[i]][[indic]] <- 
      Resflow_year_all[[i]][[indic]]%>%ungroup()%>%
      popResults_range(POC_AU, ., Population = NULL,
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
  }
}


#### cumulative: number of infections: overall in 5, 10, 20y ####
# starting at year of 2022
# tidy up data in dt to yearly and cumulative 
Resflow_cum_all <- list()
par_col <- c("best", paste0("set", seq(1, POC_AU$numberSamples,1)))

for(i in names(Resflow_year_all)){ 
  for(indic in names(Resflow_year_all[[1]])){ 
    Resflow_cum_all[[i]][[indic]] <- Resflow_year_all[[i]][[indic]]%>%
      ungroup()%>%
      filter(year>= POC_AU$simY )%>%
      ungroup()%>%
      mutate(across(all_of(par_col), ~{
        i <- !is.na(.x)
        out <- rep(NA_real_, length(.x))
        out[i] <- cumsum(.x[i])
        out
      }, .names = "{col}"))
  }
}



# avert number of transitions 
Resflow_cum_all_avert <- list()

for(i in names(Resflow_cum_all)){ 
  for(indic in names(Resflow_cum_all[[1]])){
    Resflow_cum_all[[i]][[indic]][is.na(Resflow_cum_all[[i]][[indic]])] <- 0
    Resflow_cum_all_avert[[i]][[indic]] <- cbind(
      year = Resflow_cum_all[[i]][[indic]]$year, 
      data.frame(Resflow_cum_all[["no_np"]][[indic]][, c(par_col)] - 
                   Resflow_cum_all[[i]][[indic]][, c(par_col)]))
    
    
    
    
  }
}

for(i in names(Resflow_cum_all)){ 
  for(indic in names(Resflow_cum_all[[1]])){
    
    # Set those columns to NA
    Resflow_cum_all_avert[[i]][[indic]][Resflow_cum_all_avert[[i]][[indic]] == 0] <- NA
    
    
  }
}



Resflow_sc_cum_all <- list()
for(i in names(Resflow_sc_year_all)){ 
  for(indic in names(Resflow_sc_year_all[[1]])){
    
    Resflow_sc_year_all[[i]][[indic]] <- Resflow_sc_year_all[[i]][[indic]]%>%mutate(year = year + POC_AU$cabY)
  }
}

for(i in names(Resflow_sc_year_all)){ 
  for(indic in names(Resflow_sc_year_all[[1]])){
    Resflow_sc_cum_all[[i]][[indic]] <- Resflow_sc_year_all[[i]][[indic]]%>%
      ungroup()%>%
      mutate(across(all_of(par_col), ~{
        i <- !is.na(.x)
        out <- rep(NA_real_, length(.x))
        out[i] <- cumsum(.x[i])
        out
      }, .names = "{col}"))
    
    # Set those columns to NA
    Resflow_sc_cum_all[[i]][[indic]][Resflow_sc_cum_all[[i]][[indic]] == 0] <- NA
    
  }
}


# avert number of transitions 
Resflow_sc_cum_all_avert <- list()
for(i in names(Resflow_sc_cum_all)){ 
  for(indic in names(Resflow_sc_cum_all[[1]])){
    Resflow_sc_cum_all[[i]][[indic]][is.na(Resflow_sc_cum_all[[i]][[indic]])] <- 0
    Resflow_sc_cum_all_avert[[i]][[indic]] <- cbind(
      year = Resflow_sc_cum_all[[i]][[indic]]$year, 
      dplyr::bind_cols(Resflow_sc_cum_all[["no_np"]][[indic]][, c(par_col)] - 
                         Resflow_sc_cum_all[[i]][[indic]][, c(par_col)]))
  }
}

for(i in names(Resflow_sc_cum_all)){ 
  for(indic in names(Resflow_sc_cum_all[[1]])){
    Resflow_sc_cum_all_avert[[i]][[indic]][Resflow_sc_cum_all_avert[[i]][[indic]] == 0] <- NA
  }}


# aggregate q2.5 - q97.5
# Resflow_year_all, Resflow_sc_year_all
# Resflow_cum_all, Resflow_cum_all_avert, 
# Resflow_sc_cum_all, Resflow_sc_cum_all_avert 


# output as xlsx 
Resflow_year_all_range <- list()
for(i in names(Resflow_year_all)){ 
  for(indic in names(Resflow_year_all[[1]])){
    Resflow_year_all_range[[i]][[indic]] <- Resflow_year_all[[i]][[indic]]%>%
      popResults_range(POC_AU, . , Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
  }
}

for(i in names(Resflow_year_all)){ 
  write.xlsx(Resflow_year_all_range[[i]], file = file.path(OutputFig, paste0("Resflow_year_all_range_", i,cost_types[1],".xlsx")), 
             append=TRUE) 
  
}

Resflow_sc_year_all_range <- list() 

for(i in names(Resflow_sc_year_all)){ 
  for(indic in names(Resflow_sc_year_all[[1]])){
    Resflow_sc_year_all_range[[i]][[indic]] <- Resflow_sc_year_all[[i]][[indic]]%>%
      popResults_range(POC_AU, . , Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
    
  }
}

Resflow_cum_all_range <- list() 

for(i in names(Resflow_cum_all)){ 
  for(indic in names(Resflow_cum_all[[1]])){
    
    # Set those columns to NA
    Resflow_cum_all[[i]][[indic]][Resflow_cum_all[[i]][[indic]] == 0] <- NA
    
    
    Resflow_cum_all_range[[i]][[indic]] <- 
      popResults_range(POC_AU, Resflow_cum_all[[i]][[indic]], Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
  }
}

Resflow_cum_all_avert_range <- list() 

for(i in names(Resflow_cum_all_avert)){ 
  for(indic in names(Resflow_cum_all_avert[[1]])){
    Resflow_cum_all_avert_range[[i]][[indic]] <- 
      popResults_range(POC_AU, Resflow_cum_all_avert[[i]][[indic]], Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
  }
}  

Resflow_sc_cum_all_range <- list()
for(i in names(Resflow_sc_cum_all)){ 
  for(indic in names(Resflow_sc_cum_all[[1]])){
    
    # Set those columns to NA
    Resflow_sc_cum_all[[i]][[indic]][Resflow_sc_cum_all[[i]][[indic]] ==0] <- NA
    
    Resflow_sc_cum_all_range[[i]][[indic]] <- 
      popResults_range(POC_AU, Resflow_sc_cum_all[[i]][[indic]], Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
  }
}  

Resflow_sc_cum_all_avert_range <- list()
for(i in names(Resflow_sc_cum_all)){ 
  for(indic in names(Resflow_sc_cum_all[[1]])){
    Resflow_sc_cum_all_avert_range[[i]][[indic]] <- 
      popResults_range(POC_AU, Resflow_sc_cum_all_avert[[i]][[indic]], Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
  }
}

Resflow_all <- list()
Resflow_all <- list("Resflow_year" = Resflow_year_all_range, 
                    "Resflow_cum" = Resflow_cum_all_range, 
                    "Resflow_cum_avert" = Resflow_cum_all_avert_range)


Resflow_all <- lapply(Resflow_all, function(x) {
  names(x) <- sce_level 
  return(x)
})

Resflow_all_lst <- list()

for(i in names(Resflow_all)){ 
  Resflow_all_lst[[i]] <- Resflow_all[[i]]%>%purrr::transpose()%>%
    lapply(., function(x) dplyr::bind_rows(x, .id = 'scenario'))
  
  
  for(n in names(Resflow_all_lst[[i]])){ 
    Resflow_all_lst[[i]][[n]] <- Resflow_all_lst[[i]][[n]]%>%
      mutate(scenario = factor(scenario, 
                               levels = sce_level, 
                               labels = sce_label),
             sensitivity = cost_types[1])%>%
      arrange(scenario, year)
  }
}

# extracting the scenario wanna present 
Resflow_all_lst_subsce <- list()
for(i in names(Resflow_all_lst)){
  for(n in names(Resflow_all_lst[[1]])){ 
    Resflow_all_lst_subsce[[i]][[n]] <- Resflow_all_lst[[i]][[n]]
    
  }
}

for(i in names(Resflow_all_lst_subsce$Resflow_cum_avert)){ 
  Resflow_all_lst_subsce$Resflow_cum_avert[[i]] <- 
    Resflow_all_lst_subsce$Resflow_cum_avert[[i]]%>%
    filter(scenario != sce_label[1])
}
p_pocau_y <- list()
p_pocau_cum <- list()
p_pocau_avert <- list()
lim_dt_y <- list()
lim_dt_cum <- list()
lim_dt_avert <- list()
ylab_name <- list()
ylab_name <- list("New HCV Infections", 
                  "Number of HCV releated deaths",
                  "Number of Treatment",
                  "Number of Retreat",
                  "Number of antibody testing\n(HCV positive, out of National Program)",
                  "Number of RNA testing\n(two-step, HCV positive, out of National Program)",
                  "Number of point-of-care RNA testing\n(HCV positive, out of National Program)",
                  "Number of antibody testing\n(HCV negative, out of National Program)",
                  "Number of RNA testing\n(two-step, HCV negative, out of National Program)",
                  "Number of point-of-care RNA testing\n(HCV negative, out of National Program)",
                  "Number of SVR achieved",
                  "Number of treatment initiated\nvia National Program",
                  "Number of antibody testing via National Program\n(HCV positive)",
                  "Number of RNA testing via\nNational Program(HCV positive)",
                  "Number of point-of-care RNA testing\nvia National Program(HCV positive)",
                  "Number of antibody testing via\nNational Program(HCV negative)",
                  "Number of RNA testing via\nNational Program(HCV negative)",
                  "Number of point-of-care RNA testing via\nNational Program(HCV negative)",
                  "Number of total treatment initiated",
                  "Number of total antibody testing",
                  "Number of total two-step RNA testing",
                  "Number of total point-of-care RNA testing")

yavert_lab_name <- paste0("Averted ", tolower(ylab_name))
names(ylab_name) <- names(Resflow_all_lst_subsce[[1]])
names(yavert_lab_name) <- names(Resflow_all_lst_subsce[[1]])

for(indic in names(Resflow_all_lst_subsce[[1]])){
  
  lim_dt_y[[indic]] <- lim_ident(Resflow_all_lst_subsce$Resflow_year[[indic]], seq(2021, 2050, 1))
  lim_dt_cum[[indic]] <- lim_ident(Resflow_all_lst_subsce$Resflow_cum[[indic]], seq(2022, year_obs[3], 1))
  lim_dt_avert[[indic]] <- lim_ident(Resflow_all_lst_subsce$Resflow_cum_avert[[indic]], seq(2022, year_obs[3], 1))
  
  p_pocau_y[[indic]] <- 
    plot_pocau(POC_AU, Resflow_all_lst_subsce$Resflow_year, indicator = indic,
               type = "new", 5) + 
    labs( y = ylab_name[[indic]], x = "Year") + 
    theme(axis.title = element_text()) + 
    #scale_y_continuous(limits = c(0, as.numeric(lim_dt_y[[indic]][, "lim"])),
    #                   breaks = seq(0, as.numeric(lim_dt_y[[indic]][, "lim"]),
    #                                (as.numeric(lim_dt_y[[indic]][, "lim"] - 0))/10)) + 
    theme(legend.position = "right", legend.direction="vertical")
  
  p_pocau_cum[[indic]] <- 
    plot_pocau(POC_AU, Resflow_all_lst_subsce$Resflow_cum, type = "cum", 
               indicator = indic, year_obs = year_obs) + 
    ylab(ylab_name[[indic]]) + 
    theme(axis.title = element_text()) + 
    # scale_y_continuous(limits = c(0, as.numeric(lim_dt_cum[[indic]][, "lim"])),
    #                   breaks = seq(0, as.numeric(lim_dt_cum[[indic]][, "lim"]),
    #                                (as.numeric(lim_dt_cum[[indic]][, "lim"] - 0))/10)) + 
    theme(legend.position = "right", legend.direction="vertical")
  p_pocau_avert[[indic]] <- 
    plot_pocau(POC_AU, 
               Resflow_all_lst_subsce$Resflow_cum_avert, type = "avert", 
               indicator = indic, year_obs = year_obs) + 
    ylab(yavert_lab_name[[indic]]) + 
    #scale_y_continuous(limits = c(0, as.numeric(lim_dt_avert[[indic]][, "lim"])),
    #                   breaks = seq(0, as.numeric(lim_dt_avert[[indic]][, "lim"]),
    #                                (as.numeric(lim_dt_avert[[indic]][, "lim"] - 0))/10)) + 
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) + 
    
    theme(axis.title = element_text()) + 
    theme(legend.position = "right", legend.direction="vertical")
  
}

p_pocau_y$newInfections <- p_pocau_y$newInfections + 
  scale_y_continuous(limits = c(0, 6000), 
                     breaks = seq(0, 6000, 1000)) 
p_pocau_y$HCVdeath <- p_pocau_y$HCVdeath + scale_y_continuous(limits = c(0, 400), 
                                                              breaks = seq(0, 400, 100))

p_pocau_y$Treatment <- p_pocau_y$Treatment +
  scale_y_continuous(limits = c(0, 8000),breaks = seq(0, 8000, 1000)) 

p_pocau_y$Retreat <- p_pocau_y$Retreat + 
  scale_y_continuous(limits = c(0, 800),breaks = seq(0, 800, 100)) 

p_pocau_y$Testing_ab <- p_pocau_y$Testing_ab + 
  scale_y_continuous(limits = c(0, 12000),breaks = seq(0, 12000, 1000)) 

p_pocau_y$Testing_RNA <- p_pocau_y$Testing_RNA + 
  scale_y_continuous(limits = c(0, 12000),breaks = seq(0, 12000, 1000)) 

p_pocau_y$Testing_ab_neg <- p_pocau_y$Testing_ab_neg + 
  scale_y_continuous(limits = c(0, 30000),breaks = seq(0, 30000, 5000))

p_pocau_y$Testing_RNA_neg <- p_pocau_y$Testing_RNA_neg + 
  scale_y_continuous(limits = c(0, 20000),breaks = seq(0, 20000, 5000))

p_pocau_y$Cured <- p_pocau_y$Cured + 
  scale_y_continuous(limits = c(0, 8000),breaks = seq(0, 8000, 1000))

p_pocau_y$Treatment_sc <- p_pocau_y$Treatment_sc + 
  scale_y_continuous(limits = c(0, 900),breaks = seq(0, 900, 100)) 

p_pocau_y$Testing_ab_sc <- p_pocau_y$Testing_ab_sc + 
  scale_y_continuous(limits = c(0, 800),breaks = seq(0, 800, 100))

p_pocau_y$Testing_RNA_sc <- p_pocau_y$Testing_RNA_sc + 
  scale_y_continuous(limits = c(0, 100),breaks = seq(0, 100, 10))

p_pocau_y$Testing_POCT_sc <- p_pocau_y$Testing_POCT_sc + 
  scale_y_continuous(limits = c(0, 2000),breaks = seq(0, 2000, 100))

p_pocau_y$Testing_ab_sc_neg <- p_pocau_y$Testing_ab_sc_neg + 
  scale_y_continuous(limits = c(0, 25000),breaks = seq(0, 25000, 5000))

p_pocau_y$Testing_RNA_sc_neg <- p_pocau_y$Testing_RNA_sc_neg + 
  scale_y_continuous(limits = c(0, 3000),breaks = seq(0, 3000, 1000)) 

p_pocau_y$Testing_POCT_sc_neg <- p_pocau_y$Testing_POCT_sc_neg + 
  scale_y_continuous(limits = c(0, 20000),breaks = seq(0, 20000, 5000))

p_pocau_y$Tot_Treatment <- p_pocau_y$Tot_Treatment + 
  scale_y_continuous(limits = c(0, 9000),breaks = seq(0, 9000, 1000))

p_pocau_y$Tot_Testing_ab <- p_pocau_y$Tot_Testing_ab + 
  scale_y_continuous(limits = c(0, 50000),breaks = seq(0, 50000, 5000))

p_pocau_y$Tot_Testing_RNA <- p_pocau_y$Tot_Testing_RNA + 
  scale_y_continuous(limits = c(0, 30000),breaks = seq(0, 30000, 5000))

p_pocau_y$Tot_Testing_POCT <- p_pocau_y$Tot_Testing_POCT + 
  scale_y_continuous(limits = c(0, 20000),breaks = seq(0, 20000, 5000))

p_pocau_avert$HCVdeath <- p_pocau_avert$HCVdeath + 
  scale_y_continuous(limits = c(0, 30),breaks = seq(0, 30, 10))

p_pocau_avert$Tot_Treatment <- p_pocau_avert$Tot_Treatment + 
  scale_y_continuous(limits = c(-1500, 1000),breaks = seq(-1500, 1000, 100))

p_pocau_avert$Tot_Testing_ab <- p_pocau_avert$Tot_Testing_ab + 
  scale_y_continuous(limits = c(-25000, 2000),breaks = seq(-25000, 2000, 1000))

p_pocau_avert$Tot_Testing_RNA <- p_pocau_avert$Tot_Testing_RNA + 
  scale_y_continuous(limits = c(-4000, 1000),breaks = seq(-4000, 1000, 500))

p_pocau_cum$newInfections <- p_pocau_cum$newInfections + 
  scale_y_continuous(limits = c(0, 250000), breaks = seq(0, 250000, 50000))

p_pocau_cum$HCVdeath <- p_pocau_cum$HCVdeath + 
  scale_y_continuous(limits = c(0, 15000), breaks = seq(0, 15000, 1000))

p_pocau_cum$Tot_Treatment <- p_pocau_cum$Tot_Treatment + 
  scale_y_continuous(limits = c(0, 250000), breaks = seq(0, 250000, 50000))

p_pocau_cum$Tot_Testing_ab <- p_pocau_cum$Tot_Testing_ab + 
  scale_y_continuous(limits = c(0, 2000000), breaks = seq(0, 2000000, 500000))
p_pocau_cum$Tot_Testing_RNA <- p_pocau_cum$Tot_Testing_RNA + 
  scale_y_continuous(limits = c(0, 2000000), breaks = seq(0, 2000000, 500000))
p_pocau_cum$Testing_POCT_sc <- p_pocau_cum$Testing_POCT_sc +
  scale_y_continuous(limits = c(0, 2000000), breaks = seq(0, 2000000, 500000))
p_pocau_avert$HCVdeath <- p_pocau_avert$HCVdeath + 
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 100))
p_pocau_avert$Tot_Treatment <- p_pocau_avert$Tot_Treatment + 
  scale_y_continuous(limits = c(-10000, 10000), breaks = seq(-10000, 10000, 1000))
for(i in names(p_pocau_avert)){ 
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0(i,"_avert" ,".png")), 
         p_pocau_avert[[i]], 
         width = 9, height = 6, bg = "white", dpi = 300)
  
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0(i,"_cum" ,".png")), 
         p_pocau_cum[[i]], 
         width = 9, height = 6, bg = "white", dpi = 300)
  
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0(i,"_y" ,".png")), 
         p_pocau_y[[i]], 
         width = 9, height = 6, bg = "white", dpi = 300)
}


#### number of treatment, testing plot with calibration points####



names(Resflow_year_setting_range) <- sce_level
Resflow_year_setting_range <- Resflow_year_setting_range%>%
  purrr::transpose()%>%
  lapply(., function(x) dplyr::bind_rows(x, .id = 'scenario'))
lst_name <- names(Resflow_year_setting_range)
Resflow_year_setting_range <- lapply(Resflow_year_setting_range, function(x){ 
  
  x <- x%>%
    mutate(scenario = factor(scenario, 
                             levels = sce_level, 
                             labels = sce_label),
           sensitivity = cost_types[1])%>%
    arrange(scenario, year)
  return(x)
  } )


  
  
HCVtreatinitN_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatinitN_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year - POC_AU$cabY + 1, 
                           realPop = realpop,
                           up = upper,
                           low = lower,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))
HCVtreatinitN_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatinitN_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year - POC_AU$cabY + 1, 
                           realPop = realpop,
                           up = upper,
                           low = lower)

HCVtreatinitN_NP_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatinitN_NP_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year - POC_AU$cabY + 1, 
                           realPop = realpop,
                           up = upper,
                           low = lower, 
                           population = factor(population, 
                                               levels = c("Total", "Natioanl Program"), 
                                               labels = c("Total", "National Program")))




total_treatm <- Resflow_year_setting_range$Tot_Treatment%>%
  group_by(year, scenario)%>%summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))
total_treatm_lst <- list()
for(i in sce_label){ 
  
  total_treatm_lst[[i]] <- total_treatm%>%filter(scenario == i)%>%
    popResults_range(POC_AU, .)
}


total_treatm_lst <- total_treatm_lst%>%bind_rows(., .id = 'scenario')%>%
  mutate(scenario = factor(scenario, levels = sce_label, label = sce_label))%>%
  mutate(year = year )

p_T_num <- Cas_num_plot(POC_AU,Resflow_year_setting_range$Tot_Treatment%>%
                          filter(scenario %in% sce_label[1:2]) , 
                        obdt = HCVtreatinitN_setting_fit, 
                        xlimits = c(2015, 2030, 5), UI = sce_label[1]) + 
  labs(x = "Year", y = "Number of total treatment initiated") + 
  scale_x_continuous(limits = c(2015,2030), breaks = seq(2015,2030,1))
p_T_num <- p_T_num + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 40000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 8000)))
                  )) +
  theme(legend.position = "right", legend.direction="vertical")

# total treatment number 
p_T_num_total <- Cas_num_plot(POC_AU,total_treatm_lst , 
                              obdt =HCVtreatinitN_fit, 
                              xlimits = c(2015, 2030, 5), UI = sce_label[1], 
                              population = "n") + 
  labs(x = "Year", y = "Number of total treatment initiated") +
  theme(legend.position = "",
        legend.direction = "vertical")

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("NP_numTreat_all", ".png")), 
       p_T_num_total , 
       width = 8, height = 6, bg = "white", dpi = 300)

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("NP_numTreat_setting", ".png")), 
       p_T_num , 
       width = 12, height = 6, bg = "white", dpi = 300)



# new inf and treatment number for main text 
p_T_num_total_maintext <- Cas_num_plot(POC_AU,total_treatm_lst , 
                                       obdt = NULL, 
                                       xlimits = c(2021, 2030, 1), UI = NULL, 
                                       population = "n") + 
  labs(x = "Year", y = "Treatment initiations", tag ="B") +
  theme(legend.position = "right",
        legend.direction = "vertical") + 
  scale_y_continuous(limits = c(0, 9000), breaks = seq(0, 9000,1000)) + 
  scale_x_continuous(expand = c(0.005, 0),limits = c(2021, 2030), breaks = seq(2021,2030, 1))

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("T_num_total_maintext", ".png")), 
       p_T_num_total_maintext, 
       width = 8, height = 6, bg = "white", dpi = 300)

p_newinf_num_total_maintext <- p_pocau_y$newInfections + 
  scale_x_continuous(expand = c(0,0), limits = c(2021, 2030), breaks = seq(2021, 2030, 1)) + 
  scale_y_continuous(limits = c(0, 9000), breaks = seq(0, 9000, 1000)) + 
  labs(x = "Year", y = "New HCV infections", tag ="A") + 
  theme(legend.position = "",
        legend.direction = "vertical") 


ggsave(file=file.path(OutputFig_y_cum_avert, paste0("newinf_num_total_maintext", ".png")), 
       p_newinf_num_total_maintext , 
       width = 8, height = 6, bg = "white", dpi = 300)

p_newinfnT_maintext  <- ggarrange(plotlist = list(p_newinf_num_total_maintext, 
                                                  p_T_num_total_maintext), ncol = 2, nrow = 1, 
                                  common.legend = TRUE, legend="bottom")

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("newinfnT_maintext ", ".png")), 
       p_newinfnT_maintext , 
       width = 10, height = 6, bg = "white", dpi = 300)


p_HCVdeath_supp <- p_pocau_y$HCVdeath + 
  scale_x_continuous(expand = c(0,0), limits = c(2021, 2041), breaks = seq(2021, 2041, 1)) + 
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 50)) + 
  labs(x = "Year", y = "Number of HCV-related deaths") + 
  theme(legend.position = "bottom",
        legend.direction = "vertical") 


ggsave(file=file.path(OutputFig_y_cum_avert, paste0("HCVdeath_supp", ".png")), 
       p_HCVdeath_supp , 
       width = 8, height = 6, bg = "white", dpi = 300)

#### testing numbers #### 

# testing 
xt_np_ab <- cbind(year = Resflow_year_setting_range$Testing_ab_sc$year, 
                  population = Resflow_year_setting_range$Testing_ab_sc$population, 
                  scenario = Resflow_year_setting_range$Testing_ab_sc$scenario, 
                  as.data.frame(Resflow_year_setting_range$Testing_ab_sc[, par_col] + 
                                  Resflow_year_setting_range$Testing_ab_sc_neg[, par_col]))%>%
  as.data.frame()%>%
  split(., .$scenario)%>%
  lapply(., function(x) popResults_range(POC_AU, x, Population = c("Community", "Prisons"),
                                         Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)%>%
           mutate(NP = "NP"))
  

xt_pnp_ab <- cbind(year = Resflow_year_setting_range$Testing_ab$year, 
                   population = Resflow_year_setting_range$Testing_ab_sc$population, 
                   scenario = Resflow_year_setting_range$Testing_ab$scenario, 
                   as.data.frame(Resflow_year_setting_range$Testing_ab[, par_col] + 
                                   Resflow_year_setting_range$Testing_ab_neg[, par_col]))%>%
  as.data.frame()%>%
  split(., .$scenario)%>%
  lapply(., function(x) popResults_range(POC_AU, x, Population = c("Community", "Prisons"),
                                         Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)%>%
  mutate(NP = "out of NP"))

# number of treatment, testing plot with calibration points 



Resflow_year_setting_range_trajectory <- lapply(Resflow_year_setting_range, function(x)
  x%>%mutate(scenario = factor(scenario, levels = sce_level, 
                               labels = sce_label)))

xt_np_rna <- cbind(year = Resflow_year_setting_range$Testing_RNA_sc$year, 
                   population = Resflow_year_setting_range$Testing_RNA_sc$population, 
                   scenario = Resflow_year_setting_range$Testing_RNA_sc$scenario, 
                   as.data.frame(
                     replace(Resflow_year_setting_range$Testing_RNA_sc[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_RNA_sc[, par_col]), 0) + 
                       replace(Resflow_year_setting_range$Testing_RNA_sc_neg[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_RNA_sc_neg[, par_col]), 0) + 
                       replace(Resflow_year_setting_range$Testing_POCT_sc[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_POCT_sc[, par_col]), 0) + 
                       replace(Resflow_year_setting_range$Testing_POCT_sc_neg[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_POCT_sc_neg[, par_col]), 0)
                   ))%>%
  as.data.frame()%>%mutate(across(all_of(par_col), ~na_if(., 0)))%>%
  split(., .$scenario)%>%
  lapply(., function(x) popResults_range(POC_AU, x, Population = c("Community", "Prisons"),
                                         Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)%>%
           mutate(NP = "NP"))

xt_pnp_rna <- cbind(year = Resflow_year_setting_range$Testing_ab$year, 
                    population = Resflow_year_setting_range$Testing_ab_sc$population, 
                    scenario = Resflow_year_setting_range$Testing_ab$scenario, 
                    as.data.frame(
                      replace(Resflow_year_setting_range$Testing_RNA[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_RNA[, par_col]), 0) + 
                        replace(Resflow_year_setting_range$Testing_RNA_neg[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_RNA_neg[, par_col]), 0) + 
                        replace(Resflow_year_setting_range$Testing_POCT[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_POCT[, par_col]), 0) + 
                        replace(Resflow_year_setting_range$Testing_POCT_neg[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_POCT_neg[, par_col]), 0)
                    ))%>%
  as.data.frame()%>%mutate(across(all_of(par_col), ~na_if(., 0)))%>%
  split(., .$scenario)%>%
  lapply(., function(x) popResults_range(POC_AU, x, Population = c("Community", "Prisons"),
                                         Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)%>%
           mutate(NP = "out of NP"))

# screened (ab + poct)
xt_np_screened <- cbind(year = Resflow_year_setting_range$Testing_ab_sc$year, 
                        population = Resflow_year_setting_range$Testing_ab_sc$population, 
                        scenario = Resflow_year_setting_range$Testing_ab_sc$scenario, 
                        as.data.frame(
                          replace(Resflow_year_setting_range$Testing_ab_sc[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_ab_sc[, par_col]), 0) +
                            replace(Resflow_year_setting_range$Testing_ab_sc_neg[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_ab_sc_neg[, par_col]), 0)+ 
                            replace(Resflow_year_setting_range$Testing_POCT_sc[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_POCT_sc[, par_col]), 0) + 
                            replace(Resflow_year_setting_range$Testing_POCT_sc_neg[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_POCT_sc_neg[, par_col]), 0)
                          
                       
                   ))%>%
  as.data.frame()%>%mutate(across(all_of(par_col), ~na_if(., 0)))%>%
  split(., .$scenario)%>%
  lapply(., function(x) popResults_range(POC_AU, x, Population = c("Community", "Prisons"),
                                         Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)%>%
           mutate(NP = "NP"))

xt_pnp_screened  <- cbind(year = Resflow_year_setting_range$Testing_ab$year, 
                    population = Resflow_year_setting_range$Testing_ab_sc$population, 
                    scenario = Resflow_year_setting_range$Testing_ab$scenario, 
                    as.data.frame(
                      Resflow_year_setting_range$Testing_ab[, par_col] + 
                        Resflow_year_setting_range$Testing_ab_neg[, par_col] + 
                        replace(Resflow_year_setting_range$Testing_POCT[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_POCT[, par_col]), 0) + 
                        replace(Resflow_year_setting_range$Testing_POCT_neg[, par_col], is.na(Resflow_year_setting_range_trajectory$Testing_POCT_neg[, par_col]), 0)
                    ))%>%
  as.data.frame()%>%mutate(across(all_of(par_col), ~na_if(., 0)))%>%
  split(., .$scenario)%>%
  lapply(., function(x) popResults_range(POC_AU, x, Population = c("Community", "Prisons"),
                                         Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)%>%
           mutate(NP = "out of NP"))

 


xt_np_ab <- dplyr::bind_rows(xt_np_ab, .id = 'scenario') 
xt_pnp_ab <- dplyr::bind_rows(xt_pnp_ab, .id = 'scenario') 
xt_np_rna <- dplyr::bind_rows(xt_np_rna , .id = 'scenario') 
xt_pnp_rna <- dplyr::bind_rows(xt_pnp_rna , .id = 'scenario') 
xt_np_screened <- dplyr::bind_rows(xt_np_screened , .id = 'scenario') 
xt_pnp_screened <- dplyr::bind_rows(xt_pnp_screened , .id = 'scenario') 

xt_ab <- rbind(xt_np_ab, xt_pnp_ab)
xt_rna <- rbind(xt_np_rna, xt_pnp_rna)
xt_screened <- rbind(xt_np_screened, xt_pnp_screened)
xt_ab <- xt_ab%>%mutate(NP = factor(NP, levels = c("out of NP", "NP"),
                                    labels = c("Out of the National Program", 
                                               "National Program")))

xt_rna <- xt_rna%>%mutate(NP = factor(NP, levels = c("out of NP", "NP"),
                                      labels = c("Out of the National Program", 
                                                 "National Program")))

xt_screened <- xt_screened%>%mutate(NP = factor(NP, levels = c("out of NP", "NP"),
                                    labels = c("Out of the National Program", 
                                               "National Program")))


HCVNP_ab_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVNPab_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year - POC_AU$cabY , 
                           best = ab,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))%>%
  mutate(scenario = ifelse(year + POC_AU$cabY == 2024, "dfList_NP_2023",  "foundational"),
         NP = "NP")%>%
  
  mutate( NP = factor(NP, 
                      labels = c(
                        "National program")))

HCVNP_rna_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVNPab_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year - POC_AU$cabY , 
                           best = RNA,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))%>%
  mutate(scenario = ifelse(year + POC_AU$cabY == 2024, "dfList_NP_2023" ,"foundational"),
         NP = "NP")%>%
  
  mutate(NP = factor(NP, 
                     labels = c(
                       "National Program")))
HCVNP_ab_setting_fit_all <- HCVNP_ab_setting_fit%>%ungroup()%>%
  group_by(year, scenario)%>%summarise(best= sum(best))%>%ungroup()
View(HCVNP_ab_setting_fit_all)
HCVNP_ab_setting_fit_lst <- list()
HCVNP_ab_setting_fit_lst_all <- list()
HCVNP_rna_setting_fit_lst <- list()
HCVNP_rna_setting_fit_lst_all <- list()
View(HCVNP_ab_setting_fit)
for(i in sce_label){ 
  HCVNP_ab_setting_fit_lst[[i]] <- as.data.frame(HCVNP_ab_setting_fit)%>%
    mutate(scenario = i)
  
  HCVNP_rna_setting_fit_lst[[i]] <- as.data.frame(HCVNP_rna_setting_fit)%>%
    mutate(scenario = i)
  
}
HCVNP_ab_setting_fit_all <- list()
HCVNP_rna_setting_fit_all <- list()
for(i in names(HCVNP_ab_setting_fit_lst)){ 
  HCVNP_ab_setting_fit_all[[i]] <- HCVNP_ab_setting_fit_lst[[i]]%>%
    group_by(year, scenario)%>%summarise(best= sum(best)) 
  
  HCVNP_rna_setting_fit_all[[i]] <- HCVNP_rna_setting_fit_lst[[i]]%>%
    group_by(year, scenario)%>%summarise(best= sum(best))
  
}

View(xt_ab)


xt_ab[is.na(xt_ab)] <- 0
xt_rna[is.na(xt_rna)] <- 0
parea_ab <- list()
parea_ab[["No National Program"]] <- 
  ggplot(data = xt_ab%>%filter(scenario == sce_label[1]), aes(x = year, y = best)) +
  geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
            alpha=0.6 , size=0.2, colour="black") +
  scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
  facet_wrap(~population, scales='free_y') + 
  theme_Publication_facet() + 
  scale_x_continuous(limits = c(2015,2030), breaks = seq(2015,2030,5), 
                     labels = seq(2015,2030,5) )  + 
  labs(x = "Year", y = "Number of antibody testing", title = "No National Program") + 
  facet_custom (~population,
                scales = "free", ncol = 1,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 30000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 30000)))))
parea_rna <- list()
parea_rna[["No National Program"]] <- 
  ggplot(data = xt_rna%>%filter(scenario == sce_label[1]), aes(x = year, y = best)) +
  geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
            alpha=0.6 , size=0.2, colour="black") +
  scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
  facet_wrap(~population, scales='free_y') + 
  theme_Publication_facet() + 
  scale_x_continuous(limits = c(2015,2030), breaks = seq(2015,2030,5), 
                     labels = seq(2015,2030,5))  + 
  labs(x = "Year", y = "Number of RNA testing", title = "No National Program") + 
  facet_custom (~population,
                scales = "free", ncol = 1,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 30000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 25000)))))


for(i in names(HCVNP_ab_setting_fit_lst)[-1]){ 
  parea_ab[[i]] <- ggplot(data = xt_ab %>% 
                            filter(scenario == i) %>%
                            mutate(best = replace(best, is.na(best), 0)), 
                          aes(x = year, y = best)) +
    geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
              alpha=0.6 , size=0.2, colour="black") +
    scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
    geom_point(data=HCVNP_ab_setting_fit_lst[[i]], aes(y=best, x = time), 
               colour = "black", size = 1) + 
    facet_wrap(~population, scales='free_y') + 
    theme_Publication_facet() + 
    scale_x_continuous(limits = c(1,16), breaks = seq(1,16,5), 
                       labels = seq(1,16,5) + POC_AU$cabY - 1)  + 
    labs(x = "Year", y = "Number of antibody testing", title = i) + 
    facet_custom (~population,
                  scales = "free", ncol = 1,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 50000))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 50000)))))
  
  parea_rna[[i]] <- ggplot(data = xt_rna %>% 
                             filter(scenario == i) %>%
                             mutate(best = replace(best, is.na(best), 0)), 
                           aes(x = year, y = best)) +
    geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
              alpha=0.6 , size=0.2, colour="black") +
    scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
    geom_point(data=HCVNP_rna_setting_fit_lst[[i]], aes(y=best, x = time), 
               colour = "black", size = 1) + 
    facet_wrap(~population, scales='free_y') + 
    theme_Publication_facet() + 
    scale_x_continuous(limits = c(1,16), breaks = seq(1,16,5), 
                       labels = seq(1,16,5) + POC_AU$cabY - 1)  + 
    labs(x = "Year", y = "Number of RNA testing", title = i) + 
    facet_custom (~population,
                  scales = "free", ncol = 1,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 50000))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 50000)))))
  
}



# total tests done in np 

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("T_num" ,".png")), 
       p_T_num, 
       width = 12, height = 8, bg = "white", dpi = 300)

# total test 
xt_ab <- xt_ab%>%arrange(year, population, scenario, NP)%>%ungroup()
xt_rna <- xt_rna%>%arrange(year, population, scenario, NP)%>%ungroup()
xt_toltest <- cbind(year = xt_ab$year, 
                    population = xt_ab$population,
                    scenario = xt_ab$scenario, 
                    NP = xt_ab$NP,
                    as.data.frame(xt_ab[, c(par_col)] + 
                                    xt_rna[, c(par_col)]))%>%as.data.frame()%>%
  group_by(year, scenario, NP)%>%
  summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))


xt_toltest_lst <- list()
for(i in unique(xt_toltest$scenario)){ 
  xt_toltest_lst[[i]] <- xt_toltest%>%filter(scenario == i)
  
}

xt_toltest_lst[[sce_label[1]]] <- 
  xt_toltest_lst[[sce_label[1]]]%>% ungroup() %>%
  mutate(dt = NA, dt_exp = NA)

xt_toltest_lst[[sce_label[2]]]<- 
  xt_toltest_lst[[sce_label[2]]]%>% ungroup() %>%
  mutate(dt = ifelse(year == 2022 & NP != "National Program", 6667, 
                     ifelse(year == 2023 & NP != "National Program", 11476, 
                            ifelse(year == 2024 & NP != "National Program", 20393, NA))), 
         dt_exp = NA)

xt_toltest_lst[[sce_label[3]]]<- 
  xt_toltest_lst[[sce_label[3]]]%>% ungroup() %>%
  mutate(dt = ifelse(year == 2022 & NP != "National Program", 6667, 
                     ifelse(year == 2023 & NP != "National Program", 11476, 
                            ifelse(year == 2024 & NP != "National Program", 20393, NA))), 
         dt_exp = ifelse(year %in% c(2025) & NP != "National Program", 25000,
                         ifelse(year %in% c(2026) & NP != "National Program", 25000, NA)))

xt_toltest_lst[[sce_label[4]]]<- 
  xt_toltest_lst[[sce_label[4]]]%>% ungroup() %>%
  mutate(dt = ifelse(year == 2022 & NP != "National Program", 6667, 
                     ifelse(year == 2023 & NP != "National Program", 11476, 
                            ifelse(year == 2024 & NP != "National Program", 20393, NA))),  
         dt_exp = ifelse(year %in% c(2025:2030) & NP != "National Program", 25000, NA))

xt_toltest_lst[[sce_label[5]]]<- 
  xt_toltest_lst[[sce_label[5]]]%>% ungroup() %>%
  mutate(dt = ifelse(year == 2022 & NP != "National Program", 6667, 
                     ifelse(year == 2023 & NP != "National Program", 11476, 
                            ifelse(year == 2024 & NP != "National Program", 20393, NA))),  
         dt_exp = ifelse(year %in% c(2025:2026) & NP != "National Program", 25000, 
                         ifelse(year %in% c(2027) & NP != "National Program", 31250,
                                ifelse(year %in% c(2028) & NP != "National Program", 37500,
                                       ifelse(year %in% c(2029) & NP != "National Program", 43750,
                                              ifelse(year %in% c(2030) & NP != "National Program", 50000, NA))))))

parea_tol <- list()

for(i in names(HCVNP_ab_setting_fit_lst[-1])){ 
  if(i == sce_label[2]){ 
    parea_tol[[i]] <- ggplot(xt_toltest_lst[[i]]) +
      geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
                alpha=0.6 , size=0.2, colour="black") +
      scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
      geom_point(data = xt_toltest_lst[[i]] %>% 
                   filter(!is.na(dt)) %>% 
                   distinct(year, dt),
                 aes(y = dt, x = year), colour = "black", size = 2)  +
      theme_Publication() + 
      scale_x_continuous(limits = c(2015,2030), breaks = seq(2015,2030,5), 
                         labels = seq(2015,2030,5) )  + 
      labs(x = "Year", y = "Number of total tests (thousands)", title = i) + 
      scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, 10000),
                         labels = seq(0, 100000, 10000)/1000)
  } 
  else{ 
    parea_tol[[i]] <- ggplot(xt_toltest_lst[[i]]) +
      geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
                alpha=0.6 , size=0.2, colour="black") +
      scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
      geom_point(aes(y=dt, x = year), colour = "black", size = 2
      ) + 
      geom_point(aes(y=dt_exp, x = year) , shape = 1, colour = "black", size = 2
      ) +
      theme_Publication() + 
      scale_x_continuous(limits = c(2015,2030), breaks = seq(2015,2030,5), 
                         labels = seq(2015,2030,5))  + 
      labs(x = "Year", y = "Number of total tests (thousands)", title = i) + 
      scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, 10000),
                         labels = seq(0, 100000, 10000)/1000)
  }
  
}


parea_tol[[names(HCVNP_ab_setting_fit_lst)[1]]] <- 
  ggplot(xt_toltest_lst[[names(HCVNP_ab_setting_fit_lst)[1]]]) +
  geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
            alpha=0.6 , size=0.2, colour="black") +
  scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
  
  theme_Publication() + 
  scale_x_continuous(limits = c(2015,2030), breaks = seq(2015,2030,5), 
                     labels = seq(2015,2030,5))  + 
  labs(x = "Year", y = "Number of total tests (thousands)", title = names(HCVNP_ab_setting_fit_lst)[1]) + 
  scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, 10000),
                     labels = seq(0, 100000, 10000)/1000)


parea_tol <- lapply(parea_tol, function(x) x + rremove("ylab") + rremove("xlab"))

parea_tol
parea_tol_ggarrange <- ggarrange(plotlist = list(parea_tol[[sce_label[1]]],
                                              parea_tol[[sce_label[2]]],
                                              parea_tol[[sce_label[3]]],
                                              parea_tol[[sce_label[4]]],
                                              parea_tol[[sce_label[5]]]), ncol = 3, nrow = 2, 
                                 common.legend = TRUE)
# adding common x, y label 
parea_tol_ggarrange <- annotate_figure(parea_tol_ggarrange, 
                                       left = textGrob("Number of total tests (thousands)", 
                                                       rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                       bottom = textGrob("Year", gp = gpar(cex = 1.3)))


ggsave(file=file.path(OutputFig_y_cum_avert, paste0("NP_numtotal_arrange", ".png")), 
       parea_tol_ggarrange, 
       width = 12, height = 8, bg = "white", dpi = 300)


################################################################################
# number of advanced liver diseases 
Res_Numbox_y <- list()
Res_Numbox_cum <- list()
Res_Numbox_avert <- list()
# output in timestep finding mid-year 

for(i in names(Res_numbox)){ 
  Res_Numbox_y[[i]][["DC"]] <- Res_numbox[[i]]%>%filter(disease_prog == "dc")%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
    filter(timestep%%1 == 0.5)%>%mutate(year = timestep%/%1)%>%
    select(year, par_col)
  
  Res_Numbox_y[[i]][["HCC"]] <- Res_numbox[[i]]%>%filter(disease_prog == "hcc")%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
    filter(timestep%%1 == 0.5)%>%mutate(year = timestep%/%1)%>%
    select(year, par_col)
  
  Res_Numbox_y[[i]][["LT"]] <- Res_numbox[[i]]%>%filter(disease_prog == "lt")%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
    filter(timestep%%1 == 0.5)%>%mutate(year = timestep%/%1)%>%
    select(year, par_col)
  
  Res_Numbox_y[[i]][["PLT"]] <- Res_numbox[[i]]%>%filter(disease_prog == "plt")%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ mean(.x, na.rm = FALSE)))%>%
    filter(timestep%%1 == 0.5)%>%mutate(year = timestep%/%1)%>%
    select(year, par_col)
  
}



for(i in names(Res_Numbox_y)){
  for(indic in names(Res_Numbox_y[[1]])){ 
    Res_Numbox_cum[[i]][[indic]] <- Res_Numbox_y[[i]][[indic]]%>%
      mutate(year = year + POC_AU$cabY - 1)%>%
      filter(year >= POC_AU$simY )%>%ungroup()%>%
      mutate(across(c(par_col), cumsum, .names = "{col}"))
    
  }
}

for(i in names(Res_Numbox_y)){ 
  for(indic in names(Res_Numbox_y[[1]])){
    Res_Numbox_avert[[i]][[indic]] <- cbind(
      year = Res_Numbox_cum[[i]][[indic]]$year, 
      data.frame(Res_Numbox_cum[["no_np"]][[indic]][, c(par_col)] - 
                   Res_Numbox_cum[[i]][[indic]][, c(par_col)]))
  }
}
Res_Numbox_y_range <- list()
Res_Numbox_cum_range <- list()
Res_Numbox_avert_range <- list()
for(i in names(Res_Numbox_y)){ 
  for(indic in names(Res_Numbox_y[[1]])){
    Res_Numbox_y_range[[i]][[indic]] <- Res_Numbox_y[[i]][[indic]]%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    Res_Numbox_cum_range[[i]][[indic]] <- Res_Numbox_cum[[i]][[indic]]%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    Res_Numbox_avert_range[[i]][[indic]] <- Res_Numbox_avert[[i]][[indic]]%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
  }
  
}
names(Res_Numbox_y_range) <- sce_level
names(Res_Numbox_cum_range) <- sce_level
names(Res_Numbox_avert_range) <- sce_level
Res_Numbox_y_range <- Res_Numbox_y_range%>%purrr::transpose()
Res_Numbox_cum_range <- Res_Numbox_cum_range%>%purrr::transpose()
Res_Numbox_avert_range <- Res_Numbox_avert_range%>%purrr::transpose()

for(i in names(Res_Numbox_y_range)){ 
  Res_Numbox_y_range[[i]] <- dplyr::bind_rows(Res_Numbox_y_range[[i]], .id = "scenario")%>%
    mutate(year = year + POC_AU$cabY - 1)%>%
    mutate(scenario = factor(scenario, 
                             levels = sce_level, 
                             labels = sce_label))
  
  Res_Numbox_cum_range[[i]] <- dplyr::bind_rows(Res_Numbox_cum_range[[i]], .id = "scenario")%>%
    mutate(scenario = factor(scenario, 
                             levels = sce_level, 
                             labels = sce_label))
  Res_Numbox_avert_range[[i]] <- dplyr::bind_rows(Res_Numbox_avert_range[[i]], .id = "scenario")%>%
    mutate(scenario = factor(scenario, 
                             levels = sce_level, 
                             labels = sce_label))
  
}

for(i in names(Res_Numbox_avert_range)){ 
  Res_Numbox_avert_range[[i]] <- Res_Numbox_avert_range[[i]]%>%
    filter(scenario != sce_label[1])
}



lim_adliver_y <- list()
lim_adliver_cum <- list()
lim_adliver_avert <- list()
for(i in names(Res_Numbox_y_range)){ 
  
  lim_adliver_y[[i]] <- lim_ident(Res_Numbox_y_range[[i]], seq(2021, 2050, 1))
  lim_adliver_cum[[i]] <- lim_ident(Res_Numbox_cum_range[[i]], seq(2021, 2050, 1))
  lim_adliver_avert[[i]] <- lim_ident(Res_Numbox_avert_range[[i]], seq(2021, 2050, 1))
}


p_num_adliver_y <- list()
p_num_adliver_cum <- list()
p_num_adliver_avert <- list()
y_adlab_name <- list("Decompensated cirrhosis", 
                     "Hepatocellular carcinoma", 
                     "Liver transplant", 
                     "Post-liver transplant")
names(y_adlab_name) <- names(Res_Numbox_y_range)

for(i in names(Res_Numbox_y_range)){ 
  
  p_num_adliver_y[[i]] <- plot_pocau(POC_AU, 
                                     Res_Numbox_y_range, type = "new", 
                                     indicator = i) +
    scale_x_continuous(limits = c(2022, 2041), breaks = c(seq(2022, 2037, 5), 2041)) + 
    scale_y_continuous(limits = c(0, as.numeric(lim_adliver_y[[i]][, "lim"]))) + 
    theme(legend.position = "right", legend.direction="vertical") + 
    labs(y = y_adlab_name[[i]], x = "Year")  + 
    theme(theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  p_num_adliver_cum[[i]] <- plot_pocau(POC_AU, Res_Numbox_cum_range, type = "cum", 
                                       indicator = i, year_obs = year_obs) + 
    
    scale_y_continuous(limits = c(0, as.numeric(lim_adliver_cum[[i]][, "lim"]))) +
    theme(legend.position = "right", legend.direction="vertical") + 
    labs(y = y_adlab_name[[i]], x = "Year") 
  
  p_num_adliver_avert[[i]] <- plot_pocau(POC_AU, Res_Numbox_avert_range, type = "avert", 
                                         indicator = i, year_obs = year_obs) + 
    
    scale_y_continuous(limits = c(0, as.numeric(lim_adliver_avert[[i]][, "lim"]))) +
    theme(legend.position = "right", legend.direction="vertical") + 
    labs(y = y_adlab_name[[i]] , x = "Year") 
}

p_num_adliver_y$DC <- p_num_adliver_y$DC + scale_y_continuous(limits = c(0, 2000))
p_num_adliver_y$HCC <- p_num_adliver_y$HCC + scale_y_continuous(limits = c(0, 500))
p_num_adliver_y$LT <- p_num_adliver_y$LT + scale_y_continuous(limits = c(0, 100))
p_num_adliver_y$PLT <- p_num_adliver_y$PLT + scale_y_continuous(limits = c(0, 100))

p_num_adliver_avert$DC <- p_num_adliver_avert$DC + 
  scale_y_continuous(limits = c(0, 1000)) + 
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") + 
  labs(y = "Number averted of Decompensated cirrhosis")

p_num_adliver_avert$HCC <- p_num_adliver_avert$HCC + 
  scale_y_continuous(limits = c(-5, 150)) + 
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") +
  labs(y = "Number averted of Hepatocellular carcinoma")

p_num_adliver_avert$LT <- p_num_adliver_avert$LT + 
  scale_y_continuous(limits = c(-5, 60)) + 
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") +
  labs(y = "Number averted of liver transplant")

p_num_adliver_avert$PLT <- p_num_adliver_avert$PLT + 
  scale_y_continuous(limits = c(-1, 10)) + 
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") +
  labs(y = "Number averted of post-liver transplant") 

p_y_adliver_arrange <- ggarrange(plotlist = p_num_adliver_y, ncol = 2, nrow = 2, 
                                 common.legend = TRUE, legend = "bottom")
p_cum_adliver_arrange <- ggarrange(plotlist = p_num_adliver_cum, ncol = 2, nrow = 2, 
                                   common.legend = TRUE, legend = "bottom")
p_avert_adliver_arrange <- ggarrange(plotlist = p_num_adliver_avert, ncol = 2, nrow = 2, 
                                     common.legend = TRUE, legend = "bottom")

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("num_adliver_y"  ,".png")), 
       p_y_adliver_arrange, 
       width = 8, height = 10, bg = "white", dpi = 300)

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("cumnum_adliver",".png")), 
       p_cum_adliver_arrange, 
       width = 8, height = 10, bg = "white", dpi = 300)  

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("avert_adliver_", ".png")), 
       p_avert_adliver_arrange , 
       width = 8, height = 11, bg = "white", dpi = 300)


# save individual figures 
for(i in names(p_num_adliver_y)){ 
  ggsave(file=file.path(OutputFig_y_cum_avert, 
                        paste0("num_adliver_",i  ,".png")), 
         p_num_adliver_y[[i]], 
         width = 10, height = 6, bg = "white", dpi = 300)
  
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0("cumnum_adliver_",i  ,".png")), 
         p_num_adliver_cum[[i]], 
         width = 10, height = 6, bg = "white", dpi = 300)
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0("avert_adliver_",i  ,".png")), 
         p_num_adliver_avert[[i]], 
         width = 10, height = 6, bg = "white", dpi = 300)
  
}

################################################################################

####                Simple and easy to understand plots                     #### 


################################################################################
# plots: x: 5y, 10y, 20y; y: the averted number of {a,b,c,d} in each scenario
# (a) newinfections; (b) HCV-related Deathes; (a) DC; (b) HCC; (c) LT; (d) PLT 
# (a) & (b) flows 
# (c) to (f) Num_box

# (a) to (b) 
# first select the list of dtset we would like to generate the plots, which is aligned with the sce_lab below 
# chose the indicator want to generate plots: new infections and HCV releated death 
Resflow_cumavert <- Resflow_cum_all_avert_range
nindic <- list("newInfections" = "New HCV Infections", "HCVdeath" = "HCV-related Deathes")
plot_num_avert <- list()
nlim_avert <- list("newInfections" = 20000, "HCVdeath" = 500)
nnumtag <- list("a", "b")
names(nnumtag) <- names(nindic)
for(indic in names(nindic)){ 
  for(s in names(Resflow_cumavert)){ 
    plot_num_avert[[s]][[indic]] <- 
      ggplot(Resflow_cumavert[[s]][[indic]]%>%
               filter(year %in% (year_obs)), 
             aes(x = as.character(year), y = best)) + 
      geom_bar(stat="identity") +
      geom_text(aes(label=round(best, digits = 0)), position = position_dodge(width = 0.55),
                hjust = 0.5, vjust = -0.5) + 
      scale_x_discrete(labels = c(paste0((year_obs - POC_AU$simY + 1), "-Year", sep = ""))) + 
      labs(y = "Number averted", x = "Time frame", tag = nnumtag[[indic]]) + 
      ggtitle(nindic[[indic]]) + 
      scale_y_continuous(limits = c(0,nlim_avert [[indic]])) +
      theme_Publication() 
  }
  
}



names(plot_num_avert) <- sce_level
ntitle_avert <- list("Decompensated cirrhosis", "Hepatocellular carcinoma", 
                     "Liver transplant", "Post-liver transplant") 
ntag <- list("a", "b", "c", "d")
names(ntag) <- names(Res_Numbox_avert_range)
names(ntitle_avert) <- names(Res_Numbox_avert_range)
lim_avert <- list(30, 10, 5, 5 )
names(lim_avert) <- names(Res_Numbox_avert_range)
# (c) to (f)

for(indic in names(Res_Numbox_avert_range)){
  for(s in 1: (length(sce_level)-1)){ 
    plot_num_avert[[s]][[indic]] <- 
      ggplot(Res_Numbox_avert_range[[indic]]%>%
               filter(scenario == sce_label[s+1] & year %in% (year_obs)), 
             aes(x = as.character(year), y = best)) + 
      geom_bar(stat="identity") + 
      geom_text(aes(label=round(best, digits = 0)),position = position_dodge(width = 0.55),
                hjust = 0.5, vjust = -0.5) + 
      scale_x_discrete(labels = c(paste0((year_obs - POC_AU$simY + 1), "-Year", sep = ""))) + 
      labs(y = "Number averted", x = "Time frame", tags = ntag[[indic]]) + 
      ggtitle(ntitle_avert[[indic]]) +
      theme_Publication() + 
      scale_y_continuous(limits = c(0,lim_avert[[indic]]))
    
  }
} 


names(plot_num_avert) <- c(sce_label[2:5], sce_label[1])
# adjust achievement 2024

for(i in names(plot_num_avert)){ 
  plot_num_avert[[i]]$newInfections <- 
    plot_num_avert[[i]]$newInfections + scale_y_continuous(limits = c(0,15000))
  plot_num_avert[[i]]$HCVdeath <- 
    plot_num_avert[[i]]$HCVdeath + scale_y_continuous(limits = c(0,100))
  plot_num_avert[[i]]$DC <- plot_num_avert[[i]]$DC + scale_y_continuous(limits = c(0,200))
  
  plot_num_avert[[i]]$HCC <- plot_num_avert[[i]]$HCC + scale_y_continuous(limits = c(0,100))
  
  plot_num_avert[[i]]$LT <- plot_num_avert[[i]]$LT + scale_y_continuous(limits = c(0,10))
  
  plot_num_avert[[i]]$PLT <- plot_num_avert[[i]]$PLT + scale_y_continuous(limits = c(0,5))
}

## arrange the plot (a) to (f) adding common title and tag 

plot_numavert <- list() 
plot_flowavert <- list()
for(i in names(plot_num_avert)[1:4]){ 
  plot_flowavert[[i]] <- ggarrange(plotlist = 
                                     list(plot_num_avert[[i]]$newInfections, 
                                          plot_num_avert[[i]]$HCVdeath),
                                   ncol = 2, nrow = 1, common.legend = TRUE,
                                   legend="bottom")
}
for(i in names(plot_num_avert)[1:4]){ 
  plot_flowavert[[i]] <- annotate_figure(plot_flowavert[[i]], top = text_grob(i, 
                                                                              color = "black", 
                                                                              face = "bold", size = 14))
}
plot_numavert_arrange <- list()
for(i in names(plot_num_avert)[1:4]){
  plot_numavert_arrange[[i]] <- ggarrange(plotlist = 
                                    list(plot_num_avert[[i]]$DC, 
                                         plot_num_avert[[i]]$HCC, 
                                         plot_num_avert[[i]]$LT, 
                                         plot_num_avert[[i]]$PLT),
                                  ncol = 2, nrow = 2, common.legend = TRUE,
                                  legend="bottom")

  plot_numavert_arrange[[i]] <- annotate_figure(plot_numavert_arrange[[i]], 
                                        top = text_grob(i, 
                                                        color = "black", 
                                                        face = "bold", size = 14))
  
}


#dir.create(file.path(paste0(OutputFig, "/Reports")))
for(i in names(plot_flowavert)){ 
  
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0("flow_cumavert_",i,".png")), 
         plot_flowavert[[i]], 
         width = 8 , height = 6, bg = "white")
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0("box_cumavert_",i,".png")), 
         plot_numavert_arrange[[i]], 
         width =12, height = 8, bg = "white")
  
}


##### 
#####
#####################################
# cost data plot 
# extract yearly value 
cap <- 200000000
cost_year_all <- list()
# load(file.path(OutputFolder, paste0(project_name, "Res_flowcost_",cost_types[1],".rda")))
cost_y_categories <- list()
cost_disyear_categories <- list()
rda2list <- function(file) {
  e <- new.env()
  load(file, envir = e)
  as.list(e)
}

files <- list.files(OutputFolder, pattern = paste0(project_name,"Res_flowcost_"))

Resflowcost_dt <- Map(rda2list, file.path(OutputFolder, files))

name_file <- sub("POC_AURes_flowcost_", "", files)

names(Resflowcost_dt) <- tools::file_path_sans_ext(name_file)



Rescost_year_all <- list()
Rescost_disyear_all <- list()
# Rescost_year_all$cost_type[1:4]$sce_level

for(i in c("DAAcost_reduchalf" ,"fixednvariable" ,"total")){
  
  for(n in sce_level){
    Rescost_year_all[[i]][[n]] <- Resflowcost_dt[[i]]$Rescost_year_all[[n]]
    Rescost_disyear_all[[i]][[n]] <- Resflowcost_dt[[i]]$Rescost_disyear_all[[n]]
  }
  names(Rescost_year_all[[i]]) <- c(sce_label)
  names(Rescost_disyear_all[[i]]) <- c(sce_label)
}


cost_y_categories <- list()
cost_disyear_categories <- list()
for (i in names(Rescost_year_all)) {
  for (n in names(Rescost_year_all[[i]])) {
    
    yr_obj  <- Rescost_year_all[[i]][[n]]
    dis_obj <- Rescost_disyear_all[[i]][[n]]
    
    # ---- NA / negative cleanup, once up front ----
    cost_cols_to_clean <- c("cost_ab", "cost_RNA", "cost_POCT",
                            "cost_compartment", "cost_Cured",
                            "cost_TreatOther", "cost_RetreatOther",
                            "cost_fibroscan", "cost_totalDAA",
                            "cost_totalDAA_Cap")
    
    for (col_name in cost_cols_to_clean) {
      yr_obj [[col_name]][is.na(yr_obj [[col_name]])] <- 0
      dis_obj[[col_name]][is.na(dis_obj[[col_name]])] <- 0
    }
    
    yr_obj $cost_Cured       [yr_obj $cost_Cured        < 0] <- 0
    dis_obj$cost_Cured       [dis_obj$cost_Cured        < 0] <- 0
    yr_obj $cost_TreatOther  [yr_obj $cost_TreatOther   < 0] <- 0
    dis_obj$cost_TreatOther  [dis_obj$cost_TreatOther   < 0] <- 0
    yr_obj $cost_RetreatOther[yr_obj $cost_RetreatOther < 0] <- 0
    dis_obj$cost_RetreatOther[dis_obj$cost_RetreatOther < 0] <- 0
    
    # ============================================================
    # Undiscounted yearly categories
    # ============================================================
    # Diagnosis = ab + RNA + POCT  (was previously only cost_ab — bug fixed)
    cost_y_categories[[i]][[n]][["Diagnosis"]] <- cbind(
      year = yr_obj$cost_ab$year,
      as.data.frame(
        yr_obj$cost_ab  [, par_col] +
          yr_obj$cost_RNA [, par_col] +
          yr_obj$cost_POCT[, par_col]
      )
    )
    cost_y_categories[[i]][[n]][["Diagnosis"]][cost_y_categories[[i]][[n]][["Diagnosis"]] == 0] <- NA
    cost_y_categories[[i]][[n]][["Diagnosis"]] <- cost_y_categories[[i]][[n]][["Diagnosis"]] %>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY - 1)
    
    # Treatment (capped)
    cost_y_categories[[i]][[n]][["Treatment_cap"]] <- yr_obj$cost_totalDAA_Cap
    cost_y_categories[[i]][[n]][["Treatment_cap"]][cost_y_categories[[i]][[n]][["Treatment_cap"]] == 0] <- NA
    cost_y_categories[[i]][[n]][["Treatment_cap"]] <- cost_y_categories[[i]][[n]][["Treatment_cap"]] %>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY - 1)
    
    # Treatment (uncapped)
    cost_y_categories[[i]][[n]][["Treatment"]] <- yr_obj$cost_totalDAA
    cost_y_categories[[i]][[n]][["Treatment"]][cost_y_categories[[i]][[n]][["Treatment"]] == 0] <- NA
    cost_y_categories[[i]][[n]][["Treatment"]] <- cost_y_categories[[i]][[n]][["Treatment"]] %>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY - 1)
    
    # Management = compartment + Cured + TreatOther + RetreatOther + fibroscan
    cost_y_categories[[i]][[n]][["Management"]] <- cbind(
      year = yr_obj$cost_compartment$year,
      as.data.frame(
        yr_obj$cost_compartment  [, par_col] +
          yr_obj$cost_Cured        [, par_col] +
          yr_obj$cost_TreatOther   [, par_col] +
          yr_obj$cost_RetreatOther [, par_col] +
          yr_obj$cost_fibroscan    [, par_col]
      )
    )
    cost_y_categories[[i]][[n]][["Management"]][cost_y_categories[[i]][[n]][["Management"]] == 0] <- NA
    cost_y_categories[[i]][[n]][["Management"]] <- cost_y_categories[[i]][[n]][["Management"]] %>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY - 1)
    
    cost_y_categories[[i]][[n]] <- cost_y_categories[[i]][[n]] %>%
      dplyr::bind_rows(.id = "Categories")
    
    # ============================================================
    # Discounted yearly categories
    # ============================================================
    cost_disyear_categories[[i]][[n]][["Diagnosis"]] <- cbind(
      year = dis_obj$cost_ab$year,
      as.data.frame(
        dis_obj$cost_ab  [, par_col] +
          dis_obj$cost_RNA [, par_col] +
          dis_obj$cost_POCT[, par_col]
      )
    )
    cost_disyear_categories[[i]][[n]][["Diagnosis"]][cost_disyear_categories[[i]][[n]][["Diagnosis"]] == 0] <- NA
    cost_disyear_categories[[i]][[n]][["Diagnosis"]] <- cost_disyear_categories[[i]][[n]][["Diagnosis"]] %>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY - 1)
    
    cost_disyear_categories[[i]][[n]][["Treatment_cap"]] <- dis_obj$cost_totalDAA_Cap
    cost_disyear_categories[[i]][[n]][["Treatment_cap"]][cost_disyear_categories[[i]][[n]][["Treatment_cap"]] == 0] <- NA
    cost_disyear_categories[[i]][[n]][["Treatment_cap"]] <- cost_disyear_categories[[i]][[n]][["Treatment_cap"]] %>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY - 1)
    
    cost_disyear_categories[[i]][[n]][["Treatment"]] <- dis_obj$cost_totalDAA
    cost_disyear_categories[[i]][[n]][["Treatment"]][cost_disyear_categories[[i]][[n]][["Treatment"]] == 0] <- NA
    cost_disyear_categories[[i]][[n]][["Treatment"]] <- cost_disyear_categories[[i]][[n]][["Treatment"]] %>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY - 1)
    
    cost_disyear_categories[[i]][[n]][["Management"]] <- cbind(
      year = dis_obj$cost_compartment$year,
      as.data.frame(
        dis_obj$cost_compartment  [, par_col] +
          dis_obj$cost_Cured        [, par_col] +
          dis_obj$cost_TreatOther   [, par_col] +
          dis_obj$cost_RetreatOther [, par_col] +
          dis_obj$cost_fibroscan    [, par_col]
      )
    )
    cost_disyear_categories[[i]][[n]][["Management"]][cost_disyear_categories[[i]][[n]][["Management"]] == 0] <- NA
    cost_disyear_categories[[i]][[n]][["Management"]] <- cost_disyear_categories[[i]][[n]][["Management"]] %>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY - 1)
    
    cost_disyear_categories[[i]][[n]] <- cost_disyear_categories[[i]][[n]] %>%
      dplyr::bind_rows(.id = "Categories")
  }
}

# ============================================================
# Cross-scenario aggregation (discounted)
# ============================================================
y_cost_disyear_categories <- list()
for (i in names(cost_disyear_categories)) {
  y_cost_disyear_categories[[i]] <- cost_disyear_categories[[i]] %>%
    dplyr::bind_rows(.id = "scenario") %>%
    ungroup() %>%
    mutate(sensitivity = i) %>%
    group_by(scenario, year) %>%
    dplyr::summarise(across(c(par_col), ~ sum(.x, na.rm = FALSE)))
}

y_cost_disyear_categories_range <- list()
ref_sce <- list()
for(i in names(cost_disyear_categories)){
  y_cost_disyear_categories_range[[i]] <- y_cost_disyear_categories[[i]]%>%
    gather("simulation", "estimate", -c(year, scenario))%>%
    group_by(year, scenario)%>%
    summarise(min = min(estimate, na.rm = TRUE),
              max = max(estimate, na.rm = TRUE),
              Med = median(estimate, na.rm = TRUE),
              Mu = mean(estimate, na.rm = TRUE),
              q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
              q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
              q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
              q95 = quantile(estimate, prob = 0.975, na.rm = TRUE)
    )%>%ungroup
  
  
  ref_sce[[i]] <- y_cost_disyear_categories[[i]]%>%filter(scenario == sce_label[1]) 
  
  y_cost_disyear_categories[[i]] <- y_cost_disyear_categories[[i]]%>%
    mutate(best_turning = best - ref_sce[[i]]$best)
  y_cost_disyear_categories[[i]] <- y_cost_disyear_categories[[i]]%>%
    select(scenario, year, best_turning,par_col)%>%
    mutate(scenario = factor(scenario, 
                             levels =sce_label , 
                             labels = sce_label))


}


col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
cost_turning_plot <- function(dt){ 
  col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
  fig <- ggplot(dt, aes(x = year, colour = scenario)) + 
    geom_line(aes(y = best)) + 
    scale_x_continuous(expand = c(0,0),limits = c(2022, 2045), breaks = seq(2022, 2045, 1)) + 
    
    scale_y_continuous(expand = c(0,0), limits = c(0, 700000000), breaks = seq(0, 700000000,100000000), 
                       labels = seq(0, 700000000,100000000)/1000000000) + 
    theme_Publication() + 
    
    scale_colour_manual(name = "Scenarios", values = col_pal) + 
    labs( x = "Year", y = "Annual discounted costs, billions")
  return (fig)
}

p_cost_y_turning <- lapply(y_cost_disyear_categories, function(x) 
  cost_turning_plot(x)) 

names(p_cost_y_turning) <- names(y_cost_disyear_categories)

for(i in names(p_cost_y_turning)){ 
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0("p_cost_y_turning_",i,".png")), 
         p_cost_y_turning[[i]], 
         width = 8, height = 8, bg = "white", dpi = 300)
  }

for(i in names(cost_y_categories)){ 
  cost_y_categories[[i]] <- cost_y_categories[[i]]%>%
    dplyr::bind_rows(., .id = "scenario")%>%
    mutate(scenario = factor(scenario, 
                             levels = sce_label,
                             labels = sce_label))
  
  cost_disyear_categories[[i]] <- cost_disyear_categories[[i]]%>%
    dplyr::bind_rows(., .id = "scenario")%>%
    mutate(scenario = factor(scenario, 
                             levels = sce_label,
                             labels = sce_label))
  }

cost_ydaanocap_categories <- list()
cost_ydaacap_categories <- list()
cost_disydaanocap_categories <- list()
cost_disydaacap_categories <- list()
for(i in names(cost_y_categories)){ 
  cost_ydaanocap_categories[[i]] <- cost_y_categories[[i]]%>%filter(Categories != "Treatment_cap")
  
  # undiscount cap 
  
  cost_ydaacap_categories[[i]] <- cost_y_categories[[i]]%>%filter(Categories != "Treatment")
  
  cost_disydaanocap_categories[[i]] <- cost_disyear_categories[[i]]%>%filter(Categories != "Treatment_cap")
  cost_disydaacap_categories[[i]] <- cost_disyear_categories[[i]]%>%filter(Categories != "Treatment")

  }
names(cost_ydaacap_categories)
cost_ydaacap_categories_bind <- dplyr::bind_rows(cost_ydaacap_categories, .id = "sensitivity")
cost_ydaanocap_categories_bind <- dplyr::bind_rows(cost_ydaanocap_categories, .id = "sensitivity")
cost_disydaacap_categories_bind <- dplyr::bind_rows(cost_disydaacap_categories, .id = "sensitivity")
cost_disydaanocap_categories_bind <- dplyr::bind_rows(cost_disydaanocap_categories, .id = "sensitivity")




write.xlsx(cost_ydaacap_categories_bind%>%
             select(scenario,sensitivity, Categories, year, best, min, max, 
                    Med, Mu, q5, q25, q75, q95), file = file.path(OutputFig, paste0("cost_y_daacap.xlsx")), 
           append=TRUE) 
write.xlsx(cost_ydaanocap_categories_bind%>%
             select(scenario, Categories, year, best, min, max, 
                    Med, Mu, q5, q25, q75, q95), file = file.path(OutputFig, paste0("cost_y_daanocap.xlsx")), 
           append=TRUE) 

write.xlsx(cost_disydaacap_categories_bind%>%
             select(scenario, Categories, year, best, min, max, 
                    Med, Mu, q5, q25, q75, q95), file = file.path(OutputFig, paste0("cost_disy_daacap.xlsx")), 
           append=TRUE) 

write.xlsx(cost_disydaanocap_categories_bind%>%
             select(scenario, Categories, year, best, min, max, 
                    Med, Mu, q5, q25, q75, q95), file = file.path(OutputFig, paste0("cost_disy_daanocap.xlsx")), 
           append=TRUE) 

View(cost_disydaanocap_categories_bind)
# gt_table: 4 tables by categories
# categories yearly cost and discount yearly cost to 2022- 2080 
# columns: scenarios 

# output excel files 

# benefit: Lifetime cost averted = total lifetime cost_ref -  total lifetime cost_program 
# cost: program cost: 5y

# diagnosis cost 
x_catcost <- list()
x_catcost <- lapply(list(cost_disydaacap_categories_bind,cost_disydaanocap_categories_bind), function(x) x%>%
                      filter(year>= 2022)%>%
                      group_by(scenario, sensitivity, Categories)%>%
                      mutate(across(c(par_col, "min", "max", "Med", "Mu", "q5", 
                                      "q25", "q75", "q95"), cumsum, .names = "{col}"))%>%ungroup()%>%
                      arrange(scenario, sensitivity)%>%
                      select(scenario,sensitivity, Categories, year, best, q5, q95))

names(x_catcost) <- c("discount_cap", "discount_nocap")

col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")

pcatcost <- list()
pcatcost_nocap <- list()
title_name <- c("5-Year: 2022-2026", 
                "10-Year: 2022-2031",
                "20-Year: 2022-2041")
unique(x_catcost$discount_cap$sensitivity)

x_catcost <- lapply(x_catcost, function(x){ 
  
  a <- x%>%
    mutate(sensitivity = factor(sensitivity, 
                                levels = c("fixednvariable", 
                                           "total", 
                                           "DAAcost_reduchalf"), 
                                labels = c("PBS-listed full price of DAA", 
                                           "NP total program cost",
                                           "Main analysis")))
  return(a)
  })


for(i in unique(x_catcost$discount_cap$sensitivity)) {
  pcatcost[[i]] <- list()

}
#### 20 years #### 
pcatcost$`Main analysis`$discount_nocap
for(i in unique(x_catcost$discount_cap$sensitivity)){ 
  pcatcost[[i]][[names(x_catcost)[1]]] <- x_catcost[[1]]%>%
    filter(year == year_obs[3] & sensitivity == i)%>%arrange(Categories)%>%
    ggplot(., aes(fill = Categories, y = best, x = scenario, 
                  label = round(best/1000000, digits = 1))) + 
    geom_bar(position="stack", stat="identity") + 
    theme(panel.spacing = unit(0, 'lines')) +
    scale_fill_manual(values = c( "grey10", "grey40","grey80")) + 
    theme_Publication(base_size = 16) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
    scale_y_continuous(limit = c(0, 5000000000), 
                       breaks = seq(0, 5000000000, 500000000),
                       labels = seq(0, 5000000000, 500000000)/1000000) + 
    labs(y = "Cost (discounted, millions)") + 
    geom_text(aes(x = scenario, y = best + 50000000, 
                  label = paste0(format(round(best/1000000, digits = 1), nsmall = 1), "m"), 
                  group = Categories),
              position = position_stack(vjust = 0.5), size = 6) + 
    ggtitle(paste0(i,"( DAA capped)" ))


  pcatcost[[i]][[names(x_catcost)[2]]] <- x_catcost[[2]]%>%
    filter(year == year_obs[3] & sensitivity == i)%>%arrange(Categories)%>%
    ggplot(., aes(fill = Categories, y = best, x = scenario, 
                  label = round(best/1000000, digits = 1))) + 
    geom_bar(position="stack", stat="identity") + 
    theme(panel.spacing = unit(0, 'lines')) +
    scale_fill_manual(values = c( "grey10", "grey40","grey80")) + 
    theme_Publication(base_size = 16) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
    scale_y_continuous(limit = c(0, 5000000000), 
                       breaks = seq(0, 5000000000, 500000000),
                       labels = seq(0, 5000000000, 500000000)/1000000) + 
    labs(y = "Cost (discounted, millions)") + 
    geom_text(aes(x = scenario, y = best + 110000000, 
                  label = paste0(format(round(best/1000000, digits = 1), nsmall = 1), "m"), 
                  group = Categories),
              position = position_stack(vjust = 0.5), size = 6) + 
    ggtitle(paste0(i ))
  

}
names(pcatcost)
pcatcost$`Main analysis`$discount_nocap
file_name <- c( cost_types[3], cost_types[1], cost_types[2])
for(i in 1:length(pcatcost)){ 
  for(n in names(pcatcost[[1]])){
    
    ggsave(file=file.path(OutputFig, paste0("cost_catego_20y_",file_name[i], "_", n,".png")), 
           pcatcost[[i]][[n]],  width = 10, height = 8, bg = "white", dpi = 300)  
  }
}

#### categories increase #### 
ref_catcost_cap <- x_catcost$discount_cap%>%group_by(sensitivity)%>%filter(year == 2041 & scenario == sce_label[1] & sensitivity == "PBS-listed full price of DAA")%>%
  select(sensitivity, Categories, ref_best = best, ref_q5 = q5, ref_q95 = q95)
incre_catcost_cap <- x_catcost$discount_cap%>%group_by(sensitivity)%>%filter(year == 2041& sensitivity == "PBS-listed full price of DAA" )%>%
  left_join(ref_catcost_cap, by = c("sensitivity", "Categories"))%>%
  mutate(incre_best = -(best- ref_best), 
         incre_q5 = -(q5 - ref_q5), 
         incre_q95 = -(q95 - ref_q95))%>%
  ungroup()
incre_catcost_cap <- incre_catcost_cap%>%mutate(sensitivity = "DAA capped")
incre_catcost_cap <- incre_catcost_cap%>%mutate(Categories = if_else(Categories == "Treatment_cap", "Treatment", Categories))

ref_catcost <- x_catcost$discount_nocap%>%group_by(sensitivity)%>%filter(year == 2041 & scenario == sce_label[1])%>%
  select(sensitivity, Categories, ref_best = best, ref_q5 = q5, ref_q95 = q95)
incre_catcost <- x_catcost$discount_nocap%>%group_by(sensitivity)%>%filter(year == 2041 )%>%
  left_join(ref_catcost, by = c("sensitivity", "Categories"))%>%
  mutate(incre_best = -(best- ref_best), 
         incre_q5 = -(q5 - ref_q5), 
         incre_q95 = -(q95 - ref_q95))%>%
  ungroup()

incre_catcost <- rbind(incre_catcost, incre_catcost_cap)%>%filter(scenario != sce_label[1])
View(incre_catcost)

incre_plot_sens <- list()
for(i in unique(incre_catcost$sensitivity)){ 
  if(i == "Main analysis"){ 
    gtitle <- paste0("Base case: 50% PBS-listed DAA price" )

  }else if(i == "NP total program cost"){
    gtitle <- (paste0("National Program total costs" ))
  }else if(i == "DAA capped"){
    gtitle <- "PBS-listed DAA price with annual cap"
  }else{gtitle <- i }
  incre_plot_sens[[i]] <- 
    ggplot(incre_catcost%>%filter(sensitivity == i), 
         aes(x = scenario, y = incre_best, fill = Categories)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = c( "grey30", "grey40","grey80")) + 
    theme_Publication(base_size = 14) + 
    theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
          legend.direction = "vertical") +
    theme(axis.text.y = element_text(size = 14)) +
    theme(
      legend.key.size = unit(1.5, "cm"),
      legend.key.width = unit(1.5, "cm"),    # width of legend keys
      legend.key.height = unit(1, "cm"),     # height of legend keys
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.spacing.y = unit(0.5, "cm")     # spacing between legend items
    ) + theme(legend.position = "right") +
    scale_y_continuous(limit = c(-40000000, 200000000), 
                       breaks = seq(-40000000, 200000000, 20000000),
                       labels = seq(-40000000, 200000000, 20000000)/1000000) + 
    labs(y = "Cost savings (discounted, millions)", x = "Scenarios") + 
    
    geom_text( aes(label = paste0(format(round(incre_best/1000000, 1), nsmall = 1), "m")),
               position = position_stack(vjust = 0.5), 
               size = 5, check_overlap = TRUE)  + geom_hline(yintercept = 0, linetype ="dashed") + 
    ggtitle(gtitle )
  
  
}

length(incre_plot_sens)
library(ggpubr)
library(cowplot)
print(incre_plot_sens[[1]])
incre_plot_sens[[1]] <- incre_plot_sens[[1]] + 
  theme(
    legend.key.size = unit(1.5, "cm"),
    legend.key.width = unit(1.5, "cm"),    # width of legend keys
    legend.key.height = unit(1, "cm"),     # height of legend keys
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.spacing.y = unit(0.5, "cm")     # spacing between legend items
  ) + theme(legend.position = "right")

get_legend_manual <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}
incre_plot_sens[[4]]
legend <- get_legend_manual(incre_plot_sens[[1]])

# Remove legends from all plots
plots_no_legend <- lapply(incre_plot_sens, function(p) {
  p + theme(legend.position = "none")
})

# Arrange: 3 plots + legend in 4th cell
incre_cost_sen <- ggarrange(
  plots_no_legend[[1]], 
  plots_no_legend[[3]], 
  plots_no_legend[[4]],
  legend,
  nrow = 2, ncol = 2
)

ggsave(file=file.path(OutputFig, paste0("incre_cost_sen.png")), 
       incre_cost_sen,  width = 20, height = 16, bg = "white", dpi = 300) 

for(i in 1: length(names(cost_disydaanocap_categories))){ 
  
  ggsave(file=file.path(OutputFig, paste0("incre_cost_sen_", names(cost_disydaanocap_categories)[i],".png")), 
         incre_plot_sens[[i]],  width = 12.5, height = 12, bg = "white", dpi = 300) 
  }

ggsave(file=file.path(OutputFig, paste0("incre_cost_maintext.png")), 
       incre_plot_sens[[1]],  width = 10, height = 12, bg = "white", dpi = 300) 

ggsave(file=file.path(OutputFig, paste0("incre_cost_legend.png")), 
       legend,  width = 10, height = 12, bg = "white", dpi = 300) 



tab_epi <- Resflow_all_lst
cost_qaly_range <- list() 
cost_qaly_range_disy <- list()
cost_qaly_range_ycum <- list()
cost_qaly_range_disycum <- list()
temp_df <- list()
temp_disdf <- list()

for(i in names(Rescost_year_all)){ 
  for(n in names(Rescost_year_all[[1]])){ 
    for(m in names(Rescost_year_all[[1]][[1]])){
      Rescost_year_all[[i]][[n]][[m]][Rescost_year_all[[i]][[n]][[m]] ==0 ]  <- NA
      
      Rescost_disyear_all[[i]][[n]][[m]][Rescost_disyear_all[[i]][[n]][[m]] ==0 ]  <- NA
      temp_df <- Rescost_year_all[[i]][[n]][[m]]
      temp_disdf <- Rescost_disyear_all[[i]][[n]][[m]]
      
      if (!"year" %in% names(temp_df)) {
        stop(paste("No 'year' column in temp_df for", i, n, m))
      }
      if (!"year" %in% names(temp_disdf)) {
        stop(paste("No 'year' column in temp_disdf for", i, n, m))
      }
      
      cost_qaly_range[[i]][[n]][[m]] <- 
        popResults_range(POC_AU, temp_df, end_Y = 100-1)%>%
        as_tibble()
   
      cost_qaly_range_disy[[i]][[n]][[m]] <- 
        popResults_range(POC_AU, temp_disdf, end_Y = 100-1)%>%
        as_tibble()
      
      cost_qaly_range_ycum[[i]][[n]][[m]] <- temp_df%>%
        filter(year>=2022)%>%
        mutate(across(par_col, list(cum=cumsum), .names = "{col}"))%>%
        popResults_range(POC_AU, ., end_Y = 100-1)%>%as_tibble()
      
      cost_qaly_range_disycum[[i]][[n]][[m]] <- temp_disdf%>%
        filter(year>=2022)%>%
        mutate(across(par_col, list(cum=cumsum), .names = "{col}"))%>%
        popResults_range(POC_AU, ., end_Y = 100-1)%>%
        as_tibble()
      
    }
    
    
    
  }
}


tab_costqaly <- lapply(cost_qaly_range, function(x) x%>%purrr::transpose())
tab_costqaly_disy <- lapply(cost_qaly_range_disy, function(x) x%>%purrr::transpose())
tab_costqaly_ycum <- lapply(cost_qaly_range_ycum, function(x) x%>%purrr::transpose())
tab_costqaly_disycum <- lapply(cost_qaly_range_disycum, function(x) x%>%purrr::transpose())

tab_costqaly <- lapply(tab_costqaly, function(x) lapply(x, function(y)bind_rows(y, .id = "scenario")%>%
                                                          mutate(scenario = factor(scenario, levels = sce_label,
                                                                                   labels = sce_label))))

tab_costqaly_disy <- lapply(tab_costqaly_disy, function(x) lapply(x, function(y)bind_rows(y, .id = "scenario")%>%
                                                          mutate(scenario = factor(scenario, levels = sce_label,
                                                                                   labels = sce_label))))

tab_costqaly_ycum <- lapply(tab_costqaly_ycum, function(x) lapply(x, function(y)bind_rows(y, .id = "scenario")%>%
                                                                    mutate(scenario = factor(scenario, levels = sce_label,
                                                                                             labels = sce_label))))

tab_costqaly_disycum  <- lapply(tab_costqaly_disycum , function(x) lapply(x, function(y)bind_rows(y, .id = "scenario")%>%
                                                                    mutate(scenario = factor(scenario, levels = sce_label,
                                                                                             labels = sce_label))))




tab_costqaly_lst <- list()
  
for(i in names(tab_costqaly_disy)){ 
  tab_costqaly_lst[[i]] <- list("year" = tab_costqaly[[i]], 
                                "disy" = tab_costqaly_disy[[i]], 
                                "ycum" = tab_costqaly_ycum[[i]], 
                                "disycum" = tab_costqaly_disycum[[i]])
  } 

#### table generation ####  

library("data.table")
library("formattable")
library("gt")
library("writexl")
library("gtsummary")
# name required packages
list.of.packages <- c("gapminder", "gt", "tidyverse")

# install required packages, if necessary, and load them ----
{
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages, require, character.only = TRUE)
}

# table default settings 
n = 0
c_col = c("#1e3048", "#274060", "#2f5375", "#4073a0", "#5088b9")
c_col_light_blue = c("#edf2fb", "#e2eafc", "#d7e3fc", "#ccdbfd", "#c1d3fe")
c_container_width = px(800)
c_table_width = px(650)
c_rn = 30
c_save = TRUE
c_format = "html"



# grouping rows
tab_x <- list()
for(i in names(tab_epi)){ 
  tab_x[[i]] <- tab_epi[[i]]%>%dplyr::bind_rows(., .id ="Indicators")%>%
    filter(year %in% seq(2022,2041,1))%>%
    select(Indicators, scenario, year, best, q5, q95)%>%
    mutate(best = formatC(best,  format = "fg", big.mark = ","),
           q5 = formatC(q5,  format = "fg", big.mark = ","),
           q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
    mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
    select(-c(best, q5, q95))%>%ungroup()%>%spread(scenario, vv)
  
  tab_x[[i]]%>%group_by(Indicators)%>%
    gt(groupname_col = "Indicators",
       rowname_col = "year")%>%
    gtsave(., file = file.path(OutputFig, paste0(i, ".docx")))
}

tab_cost <- list()
for(n in names(tab_costqaly_lst)){ 
  for(i in names(tab_costqaly_lst[[1]])){ 
    if(i %in% c("ycum", "disycum")){ 
      tab_cost[[n]][[i]] <- tab_costqaly_lst[[n]][[i]]%>%dplyr::bind_rows(., .id ="Indicators")%>%
        filter(year %in% seq(2022,2041,1))%>%
        select(Indicators, scenario, year, best, q5, q95)%>%
        mutate(best = formatC(best,  format = "fg", big.mark = ","),
               q5 = formatC(q5,  format = "fg", big.mark = ","),
               q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
        mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
        select(-c(best, q5, q95))%>%ungroup()%>%spread(scenario, vv)
    }else if(i %in% c("year", "disy")){
      tab_cost[[n]][[i]] <- tab_costqaly_lst[[n]][[i]]%>%dplyr::bind_rows(., .id ="Indicators")%>%
        filter(year %in% seq(2015,2041,1))%>%
        select(Indicators, scenario, year, best, q5, q95)%>%
        mutate(best = formatC(best,  format = "fg", big.mark = ","),
               q5 = formatC(q5,  format = "fg", big.mark = ","),
               q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
        mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
        select(-c(best, q5, q95))%>%ungroup()%>%spread(scenario, vv)
      
      }
    
    
    tab_cost[[n]][[i]]%>%group_by(Indicators)%>%
      gt(groupname_col = "Indicators",
         rowname_col = "year")%>%
      gtsave(., file = file.path(OutputFig, paste0("costqaly_", i,"_",n, ".docx")))
    
  }
  }


# numbox
tab_numbox <- list()

tab_numbox <- list("year" = Res_Numbox_y_range, 
                  "cumyear" = Res_Numbox_cum_range, 
                  "cumavert" = Res_Numbox_avert_range)

for(i in names(tab_numbox)){
  if(i == "year"){ 
    tab_numbox[[i]] <- tab_numbox[[i]]%>%dplyr::bind_rows(., .id ="Indicators")%>%
      filter(year %in% seq(2015,2041,1))%>%
      select(Indicators, scenario, year, best, q5, q95)%>%
      mutate(best = formatC(best,  format = "fg", big.mark = ","),
             q5 = formatC(q5,  format = "fg", big.mark = ","),
             q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
      mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
      select(-c(best, q5, q95))%>%ungroup()%>%spread(scenario, vv)
  } else{ 
    tab_numbox[[i]] <- tab_numbox[[i]]%>%dplyr::bind_rows(., .id ="Indicators")%>%
      filter(year %in% seq(2022,2041,1))%>%
      select(Indicators, scenario, year, best, q5, q95)%>%
      mutate(best = formatC(best,  format = "fg", big.mark = ","),
             q5 = formatC(q5,  format = "fg", big.mark = ","),
             q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
      mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
      select(-c(best, q5, q95))%>%ungroup()%>%spread(scenario, vv)
    
    }
  
  
  tab_numbox[[i]]%>%group_by(Indicators)%>%
    gt(groupname_col = "Indicators",
       rowname_col = "year")%>%
    gtsave(., file = file.path(OutputFig, paste0("numbox_", i, ".docx")))
}


#### CEA ####

CEAanalysis <- list()
timeframe <- c(5, 10, 20, 30, 40, 50, 60)
timeframe_name <- c("5y", "10y", "20y", "30y", "40y", "50y", "60y")
for(n in names(tab_costqaly_lst)){ 
  for(i in 1: length(timeframe)){ 
    CEAanalysis[[n]][[timeframe_name[i]]][["QALY"]] <- tab_costqaly_lst[[n]]$disycum$QALY_compartment%>%
      filter(year == POC_AU$simY + timeframe[i] - 1)%>%
      split(.$scenario)%>%map(~.x %>% select(-scenario))%>%
      replace(is.na(.), 0)
    
    CEAanalysis[[n]][[timeframe_name[i]]][["Cost"]] <- tab_costqaly_lst[[n]]$disycum$cost_total%>%
      filter(year == POC_AU$simY + timeframe[i] - 1)%>%
      split(.$scenario)%>%map(~.x %>% select(-scenario))%>%
      replace(is.na(.), 0)
    
    CEAanalysis[[n]][[timeframe_name[i]]][["Cost_cap"]] <- tab_costqaly_lst[[n]]$disycum$cost_total_Cap%>%
      filter(year == POC_AU$simY + timeframe[i] - 1)%>%
      split(.$scenario)%>%map(~.x %>% select(-scenario))%>%
      replace(is.na(.), 0)
    }
  
  }

Incre <- list()
for(m in names(CEAanalysis)){
  for(i in names(CEAanalysis[[1]])){ 
    for(n in names(CEAanalysis[[1]][[1]])){ 
      for(q in names(CEAanalysis[[1]][[1]][[1]])){ 
        Incre[[m]][[i]][[n]][[q]] <- cbind(year = CEAanalysis[[m]][[i]][[n]][[q]]$year, 
                                      dplyr::bind_cols(CEAanalysis[[m]][[i]][[n]][[q]][, par_col] - 
                                                         CEAanalysis[[m]][[i]][[n]][[sce_label[1]]][, par_col]))%>%
          as_tibble()%>%
          popResults_range(POC_AU, .)
        }
      
    }
  }
  
}


Increx <- lapply(Incre, function(x) lapply(x, function(y) y%>%purrr::transpose())) 

Increx$DAAcost_reduchalf$`20y`[[sce_label[2]]]$Cost%>%select(year, Med, q5, q95, Mu)%>%
  mutate(Mu = round(Med/1000000, digits = 1))
CEA <- list()
CEA_cap <- list()
for(m in names(Increx)){
  for(i in names(Increx[[1]])){ 
    for(n in names(Increx[[1]][[1]])){
      CEA[[m]][[i]][[n]] <- 
        cbind(year = Increx[[m]][[i]][[n]][["Cost"]]$year,
              as.data.frame(Increx[[m]][[i]][[n]][["Cost"]][, par_col]/Increx[[m]][[i]][[n]][["QALY"]][, par_col]))
      
      
      CEA_cap[[m]][[i]][[n]] <- cbind(year = Increx[[m]][[i]][[n]][["Cost_cap"]]$year,
                                as.data.frame(Increx[[m]][[i]][[n]][["Cost_cap"]][, par_col]/Increx[[m]][[i]][[n]][["QALY"]][, par_col]))
    }
  }
}


   

for(m in names(CEA)){ 
  for(i in names(CEA[[1]])){ 
    for(n in names(CEA[[1]][[1]])){ 
      CEA[[m]][[i]][[n]] <- CEA[[m]][[i]][[n]]%>%as_tibble()%>%popResults_range(POC_AU, .)
      CEA_cap[[m]][[i]][[n]] <- CEA_cap[[m]][[i]][[n]]%>%as_tibble()%>%popResults_range(POC_AU, .)
      }
    
  }
}
CEA$DAAcost_reduchalf$`20y`[[sce_label[5]]]%>%select(year, Mu, min, q5, q95)
te <- lapply(Increx, function(x) lapply(x, function(y) lapply(y, function(m) dplyr::bind_rows(m, .id = "indicator"))))


te <- lapply(te, function(x) lapply(x, function(y) dplyr::bind_rows(y, .id = "scenario")))

x_cap <- list()
x_nocap <- list()
for(i in names(te)){ 
  for(n in names(te[[1]])){ 
    x_cap[[i]][[n]] <- te[[i]][[n]]%>%filter(indicator != "Cost")%>%
      gather(sim, val, -c(indicator, scenario, year))%>%
      spread(indicator, val)
    
    x_nocap[[i]][[n]] <- te[[i]][[n]]%>%filter(indicator != "Cost_cap")%>%
      gather(sim, val, -c(indicator, scenario, year))%>%
      spread(indicator, val)
    }
  
}


PSA_dt <- list()
PSA_nocap_dt <- list()
for(i in names(x_cap)){ 
  PSA_dt[[i]] <- x_cap[[i]][["20y"]]%>%mutate(outline = ifelse(sim == "Mu", 1,0))%>%
    mutate(category = paste0(scenario, outline)) 
  PSA_nocap_dt[[i]] <- x_nocap[[i]][["20y"]]%>%mutate(outline = ifelse(sim == "Mu", 1,0))%>%
    mutate(category = paste0(scenario, outline))
  
  }

PSA_dt <- lapply(PSA_dt, function(x) x%>%filter(scenario != sce_label[1]))


PSA <- lapply(PSA_dt, function(x)
  ggplot(x, 
              aes(y = Cost_cap, x = QALY)) + 
  geom_point(aes(colour = scenario))  +
  facet_wrap(~scenario) +
  scale_color_manual(name = "Scenarios", 
                     values = c("#E69F00", "#56B4E9","#009E73", "#F0E442")) +
  geom_point(data = x%>%filter(outline == 1), color = "gray50", size = 1.5)  +
  facet_wrap(~scenario) +
  
  geom_hline(yintercept=0, color = "black", linewidth = 1) +
  geom_vline(xintercept=0, color = "black", linewidth = 1)   +
  labs(colour = "Scenarios",  x = "QALY", 
       y = "Costs (millions)") + 
  theme_bw() + 
  scale_y_continuous(limits = c(-200000000, 200000000),
                     breaks = seq(-200000000, 200000000, 100000000), 
                     labels = seq(-200000000, 200000000, 100000000)/1000000) + 
  scale_x_continuous(limits = c(-1000, 3000), breaks = seq(-1000, 3000, 1000)) + 
  theme(plot.title = element_text(hjust = 0.5)) +  
  theme(panel.background = element_rect(colour = "white"),
        plot.background = element_rect(colour = "white"),
        panel.border     = element_rect(fill = NA, colour = "black", 
                                        size = NA),
        plot.title = element_text(face = "bold",
                                  size = 12, hjust = 0.5),
        text = element_text(),
        axis.title = element_text(face = "bold",size = 12),
        axis.title.y = element_text(angle=90,vjust =1),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1,
                                   face = "bold",size = 12, colour = "black"), 
        axis.text.y = element_text(
          face = "bold",size = 14,colour = "black"),
        strip.text.x = element_text(size=14, color="black",
                                    face="bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.title = element_text(face="bold", size= 14),
        plot.margin = unit(c(10,5,5,5),"mm")) + 
  guides(color =guide_legend(direction='vertical',
                             override.aes = list(size=2)),
         fill = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_abline(aes(slope = 50000, intercept = 0,linetype = "WTP: A50,000"), colour = "black") + 
  scale_linetype_manual(name = "", values = c(2), 
                        guide = guide_legend(override.aes = list(color = c("black")))) +
  stat_ellipse(color = "gray50", 
               alpha = 0.7,
               linewidth = 0.6,
               show.legend = FALSE, 
               level = 0.95) 
)

PSA_nocap_dt <- lapply(PSA_nocap_dt, function(x) x%>%filter(scenario != sce_label[1]))
PSA_nocap <- lapply(PSA_nocap_dt, function(x)
  ggplot(x, 
         aes(y = Cost, x = QALY)) + 
    geom_point(aes(colour = scenario))  +
    facet_wrap(~scenario) +
    scale_color_manual(name = "Scenarios", 
                       values = c("#E69F00", "#56B4E9","#009E73", "#F0E442")) +
    geom_point(data = x%>%filter(outline == 1), color = "gray50")  +
    facet_wrap(~scenario) +
    
    geom_hline(yintercept=0, color = "black", linewidth = 1) +
    geom_vline(xintercept=0, color = "black", linewidth = 1)   +
    labs(colour = "Scenarios",  x = "QALY", 
         y = "Costs (millions)") + 
    theme_bw() + 
    scale_y_continuous(limits = c(-500000000, 200000000),
                       breaks = seq(-500000000, 200000000, 100000000), 
                       labels = seq(-500000000, 200000000, 100000000)/1000000) + 
    scale_x_continuous(limits = c(-1000, 3000), breaks = seq(-1000, 3000, 1000)) + 
    theme(plot.title = element_text(hjust = 0.5)) +  
    theme(panel.background = element_rect(colour = "white"),
          plot.background = element_rect(colour = "white"),
          panel.border     = element_rect(fill = NA, colour = "black", 
                                          size = NA),
          plot.title = element_text(face = "bold",
                                    size = 12, hjust = 0.5),
          text = element_text(),
          axis.title = element_text(face = "bold",size = 12),
          axis.title.y = element_text(angle=90,vjust =1),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1,
                                     face = "bold",size = 12, colour = "black"), 
          axis.text.y = element_text(
            face = "bold",size = 14,colour = "black"),
          strip.text.x = element_text(size=14, color="black",
                                      face="bold"),
          legend.text = element_text(size = 14, face = "bold"),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "vertical",
          legend.title = element_text(face="bold", size= 14),
          plot.margin = unit(c(10,5,5,5),"mm")) + 
    guides(color =guide_legend(direction='vertical',
                               override.aes = list(size=2)),
           fill = "none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
    geom_abline(aes(slope = 50000, intercept = 0,linetype = "WTP: A50,000"), colour = "black") + 
    scale_linetype_manual(name = "", values = c(2), 
                          guide = guide_legend(override.aes = list(color = c("black")))) +
    stat_ellipse(color = "gray50", 
                 alpha = 0.7,
                 linewidth = 0.6,
                 show.legend = FALSE, 
                 level = 0.95) 
)

for(i in names(PSA)){
  ggsave(file=file.path(OutputFig, paste0("PSA_cap_",i,".png")), 
         PSA[[i]], 
         width = 12, height = 8, bg = "white", dpi = 300)  
  
  ggsave(file=file.path(OutputFig, paste0("PSA_nocap_",i,".png")), 
         PSA_nocap[[i]], 
         width = 12, height = 8, bg = "white", dpi = 300) 
  
}

PSA_nocap$DAAcost_reduchalf
CEA_capbind <- list()
CEA_capbind <- lapply(CEA_cap, function(x) lapply(x, function(y) dplyr::bind_rows(y, .id = "Scenario")))
CEA_capbind <- lapply(CEA_capbind, function(x) bind_rows(x, .id = "Timeframe"))


CEA_capbind <- bind_rows(CEA_capbind, .id = "sensitivity")

CEA_capbind <- CEA_capbind%>%mutate(Timframe = factor(Timeframe, 
                                              levels = c("5y", "10y", "20y", "30y", "40y", "50y", "60y"), 
                                              labels = c("5y", "10y", "20y", "30y", "40y", "50y", "60y")))



CEA_capbind_gt <- list()
for(i in unique(CEA_capbind$sensitivity)){ 
  CEA_capbind_gt[[i]] <- CEA_capbind%>%select(Timeframe,sensitivity ,Scenario,Med, q5, q95)%>%
    filter(sensitivity == i)%>%
    mutate(best = formatC(Med,  format = "fg", big.mark = ","),
           q5 = formatC(q5,  format = "fg", big.mark = ","),
           q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
    mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
    select(-c(best, Med, q5, q95, sensitivity))%>%ungroup()%>%
    spread(Scenario, vv)
  
    CEA_capbind_gt[[i]]%>%
    gt(
      rowname_col = "Timeframe")%>%
    tab_header(
      title = md(i))%>%
    gtsave(., file = file.path(OutputFig, paste0("costqaly_CEA_",i,".docx")))
  
  }


# one-way sensitivity analysis 

CEA_bind <- lapply(CEA, function(x) lapply(x, function(y) dplyr::bind_rows(y, .id = "Scenario")))
CEA_bind <- lapply(CEA_bind, function(x) bind_rows(x, .id = "Timeframe"))
CEA_bind <- bind_rows(CEA_bind, .id = "sensitivity")

CEA_bind_gt <- list()
for(i in unique(CEA_bind$sensitivity)){
  CEA_bind_gt[[i]] <- CEA_bind%>%select(Timeframe,sensitivity ,Scenario,Med, q5, q95)%>%
    filter(sensitivity == i)%>%
    mutate(best = formatC(Med,  format = "fg", big.mark = ","),
           q5 = formatC(q5,  format = "fg", big.mark = ","),
           q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
    mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
    select(-c(best, Med, q5, q95, sensitivity))%>%ungroup()%>%
    spread(Scenario, vv)

  CEA_bind_gt[[i]]%>%
    gt(
      rowname_col = "Timeframe")%>%
    tab_header(
      title = md(i))%>%
    gtsave(., file = file.path(OutputFig, paste0("costqaly_CEA_nocap",i,".docx")))
  }





#### incremental line #### 

x_total_ref_cap <- lapply(cost_disydaacap_categories, 
                      function(x) x%>%
                        filter(year>= 2022)%>%
                        group_by(scenario, Categories)%>%
                        mutate(across(c(par_col, "min", "max", "Med", "Mu", "q5", 
                                        "q25", "q75", "q95"), cumsum, .names = "{col}"))%>%ungroup()%>%
                        arrange(scenario)%>%group_by(year, scenario)%>%
                        summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%filter(scenario == sce_label[1]))


x_catcost_total_cap <- lapply(cost_disydaacap_categories, 
                          function(x) x%>%
                            filter(year>= 2022)%>%
                            group_by(scenario, Categories)%>%
                            mutate(across(c(par_col, "min", "max", "Med", "Mu", "q5", 
                                            "q25", "q75", "q95"), cumsum, .names = "{col}"))%>%ungroup()%>%
                            arrange(scenario)%>%group_by(year, scenario)%>%
                            summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))
)

x_total_ref <- lapply(cost_disydaanocap_categories, 
                          function(x) x%>%
                            filter(year>= 2022)%>%
                            group_by(scenario, Categories)%>%
                            mutate(across(c(par_col, "min", "max", "Med", "Mu", "q5", 
                                            "q25", "q75", "q95"), cumsum, .names = "{col}"))%>%ungroup()%>%
                            arrange(scenario)%>%group_by(year, scenario)%>%
                            summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%filter(scenario == sce_label[1]))


x_catcost_total <- lapply(cost_disydaanocap_categories, 
                              function(x) x%>%
                                filter(year>= 2022)%>%
                                group_by(scenario, Categories)%>%
                                mutate(across(c(par_col, "min", "max", "Med", "Mu", "q5", 
                                                "q25", "q75", "q95"), cumsum, .names = "{col}"))%>%ungroup()%>%
                                arrange(scenario)%>%group_by(year, scenario)%>%
                                summarise(across(c(par_col, "min", "max", "Med", "Mu", "q5", 
                                                   "q25", "q75", "q95"),~ sum(.x, na.rm = FALSE)))
)

x_catcost_total$DAAcost_reduchalf%>%filter(year == 2041)%>%select(year, scenario, best, min, q5, q95)
x_total_ref_cap 
x_total_ref
x_catcost_total_cap
x_catcost_total

x_catcost_total_incre <- list()
x_catcost_total_cap_incre <- list()
x <- list()
x_cap <- list()

for(i in unique(x_catcost_total$scenario)) {
  x[[i]] <- list()
  x_catcost_total_incre[[i]] <- list()
  
  x_cap[[i]] <- list()
  x_catcost_total_cap_incre[[i]] <- list()
}


for(i in names(x_catcost_total)){
  for(n in unique(x_catcost_total[[1]]$scenario)){ 
    x[[i]][[n]] <- x_catcost_total[[i]]%>%filter(scenario == n)
    x_catcost_total_incre[[i]][[n]] <- 
      cbind(year = x[[i]][[n]]$year, scenario = x[[i]][[n]]$scenario, 
            as.data.frame(x[[i]][[n]][, c(par_col)] - x_total_ref[[i]][, c(par_col)]))%>%
      as.data.frame()
    
    x_catcost_total_incre[[i]][[n]][x_catcost_total_incre[[i]][[n]] == 0] <- NA
    x_catcost_total_incre[[i]][[n]] <- x_catcost_total_incre[[i]][[n]]%>%
      popResults_range(POC_AU, ., end_Y = 100)
    
    
    # cap 
    x_cap[[i]][[n]] <- x_catcost_total_cap[[i]]%>%filter(scenario == n)
    x_catcost_total_cap_incre[[i]][[n]] <- 
      cbind(year = x_cap[[i]][[n]]$year, scenario = x_cap[[i]][[n]]$scenario, 
            as.data.frame(x_cap[[i]][[n]][, c(par_col)] - x_total_ref_cap[[i]][, c(par_col)]))%>%
      as.data.frame()
    
    x_catcost_total_cap_incre[[i]][[n]][x_catcost_total_cap_incre[[i]][[n]] == 0] <- NA
    x_catcost_total_cap_incre[[i]][[n]] <- x_catcost_total_cap_incre[[i]][[n]]%>%
      popResults_range(POC_AU, ., end_Y = 100)
    }
  }



x_catcost_total_cap_incre <- lapply(x_catcost_total_cap_incre, function(x) bind_rows(x, .id = "scenario"))
x_catcost_total_incre <- lapply(x_catcost_total_incre, function(x) bind_rows(x, .id = "scenario"))

incremental_cost <- list()
incremental_cost_cap <- list()
sensi_name <- c("Main analysis: 50% discounted DAA cost", "25% discounted DAA cost", "PBS-listed full price of DAA", "Total Program cost")
for(i in 1: length(names(x_catcost_total_incre))){ 
  incremental_cost[[i]] <- ggplot(x_catcost_total_incre[[i]]%>%
                                 mutate(scenario = factor(scenario, levels = c(sce_label))), 
                               aes(x = year, colour = scenario) ) + 
    geom_line(aes(x = year, y = best, colour = scenario,
                  linetype = scenario), size = 1
    ) + 
    theme_Publication() + 
    scale_x_continuous(expand = c(0,0), limits = c(2022,2041), breaks = seq(2022, 2041, 1)) + 
    scale_y_continuous(limits = c(-200000000, 100000000), 
                       breaks = seq(-200000000, 100000000, 10000000), 
                       labels = seq(-200000000, 100000000, 10000000)/1000000) + 
    scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) + 
    scale_linetype_manual(values = c("dashed", "solid", "solid", "solid", "solid")) + 
    labs( y = "Incremental cost (in millions)", x = "Year") + 
    geom_hline(linetype = "dashed", yintercept = 0, size = 1) +
    ggtitle(sensi_name[i]) + 
    theme(
      legend.position = c(0.02, 0.02),
      legend.justification = c("left", "bottom"),
      legend.direction = "vertical",
      legend.title = element_text(face = "bold"),
      legend.background = element_blank(),  # no box
      legend.key = element_blank()  # no key background
    ) +
    labs(linetype = "Scenarios",
         color = "Scenarios")
  


  ggsave(file=file.path(OutputFig, paste0("incremental_cost",names(x_catcost_total_incre)[i],".png")), 
         incremental_cost[[i]], 
         width = 12, height = 8, bg = "white", dpi = 300)   
  
  
  View(x_catcost_total_incre$DAAcost_reduchalf)
  
  incremental_cost_cap[[i]] <- ggplot(x_catcost_total_cap_incre[[i]]%>%
                                    mutate(scenario = factor(scenario, levels = c(sce_label))), 
                                  aes(x = year, colour = scenario) ) + 
    geom_line(aes(x = year, y = best, colour = scenario,
                  linetype = scenario), size = 1
    ) + 
    theme_Publication() + 
    scale_x_continuous(expand = c(0,0), limits = c(2022,2041), breaks = seq(2022, 2041, 1)) + 
    scale_y_continuous(limits = c(-200000000, 100000000), 
                       breaks = seq(-200000000, 100000000, 10000000), 
                       labels = seq(-200000000, 100000000, 10000000)/1000000) + 
    scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) + 
    scale_linetype_manual(values = c("dashed", "solid", "solid", "solid", "solid")) + 
    labs( y = "Incremental cost (in millions)", x = "Year") + 
    geom_hline(linetype = "dashed", yintercept = 0, size = 1) +
    ggtitle(sensi_name[i]) +
    theme(
      legend.position = c(0.02, 0.02),
      legend.justification = c("left", "bottom"),
      legend.direction = "vertical",
      legend.title = element_text(face = "bold"),
      legend.background = element_blank(),  # no box
      legend.key = element_blank()  # no key background
    ) +
    labs(linetype = "Scenarios",
         color = "Scenarios")
  
  
  ggsave(file=file.path(OutputFig, paste0("incremental_cost_cap",names(x_catcost_total_incre)[i],".png")), 
         incremental_cost_cap[[i]], 
         width = 12, height = 8, bg = "white", dpi = 300)   
  
  }
incremental_cost_cap
#### prevalence plots for manuscript  ####
PrevInc_plot <- function(pj, dt, obdt =NULL, xlimits, UI = NULL){ 
  if(length(unique(dt$scenario)) == 2){ 
    col_pal <- c("#000000", "#E69F00")
    
  } 
  else{col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
  }
  
  
  
  if(is.null(obdt) & is.null(UI)){ 
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) +
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
  }
  else if(is.null(obdt) & !is.null(UI)){ 
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      geom_ribbon(aes(ymin = q5, ymax = q95, fill = scenario), alpha = 0.2) +
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_fill_manual(name = "Scenarios", values = col_pal ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
    
  }
  
  else if(!is.null(obdt) & is.null(UI)){
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      geom_point(data=obdt, aes(y=realPop, x = time), 
                 colour = "black", size = 1) +
      geom_segment(data = obdt, 
                   aes ( y = low, yend = up, x = time, xend = time)) +
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
    
    
  }
  else if(!is.null(obdt) & !is.null(UI)){
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      geom_ribbon(aes(ymin = q5, ymax = q95, fill = scenario), alpha = 0.2) +
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_fill_manual(name = "Scenarios", values = col_pal ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      geom_point(data=obdt, aes(y=realPop, x = time), 
                 colour = "black", size = 1) +
      geom_segment(data = obdt, 
                   aes ( y = low, yend = up, x = time, xend = time)) +
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
  }
  
  return(traj_plot)
}


# function for getting ceiling number for plots 
lim_ident <- function(dt, group_index, year_range, summar_col){ 
  if(summar_col == "max"){
    if(group_index == "pop"){ 
      lim <- dt%>%group_by(population)%>%
        filter(year %in% year_range)%>%
        summarise(x = max(max))%>%
        mutate(lim = case_when( 
          x <1 ~5, 
          x>=1 & x <10 ~ 10, 
          x>=10 & x< 30 ~ 30, 
          x>=30 & x<60 ~ 60, 
          x>=60 & x<80 ~ 80, 
          x>80 & x<=100 ~100,
          x>100 & x<=1000 ~ (x%/%100 + 1)*100,
          x >1000 & x<10000 ~ (x%/%1000 + 1)*1000
        ))
    }
    else if(group_index == "setting"){ 
      lim <- dt%>%group_by(setting)%>%
        filter(year %in% year_range)%>%
        summarise(x = max(max))%>%
        mutate(lim = case_when( 
          x <1 ~5, 
          x>=1 & x <10 ~ 10, 
          x>=10 & x< 30 ~ 30, 
          x>=30 & x<60 ~ 60, 
          x>=60 & x<80 ~ 80, 
          x>80 & x<=100 ~100,
          x>100 & x<=1000 ~ (x%/%100 + 1)*100,
          x >1000 & x<10000 ~ (x%/%1000 + 1)*1000
        ))
    }
    return(lim)  
  }
  else if(summar_col == "best"){ 
    if(group_index == "pop"){ 
      lim <- dt%>%group_by(population)%>%
        filter(year %in% year_range)%>%
        summarise(x = max(best))%>%
        mutate(lim = case_when( 
          x <1 ~5, 
          x>=1 & x <10 ~ 10, 
          x>=10 & x< 30 ~ 30, 
          x>=30 & x<60 ~ 60, 
          x>=60 & x<80 ~ 80, 
          x>80 & x<=100 ~100,
          x>100 & x<=1000 ~ (x%/%100 + 1)*100,
          x >1000 & x<10000 ~ (x%/%1000 + 1)*1000
        ))
    }
    else if(group_index == "setting"){ 
      lim <- dt%>%group_by(setting)%>%
        filter(year %in% year_range)%>%
        summarise(x = max(best))%>%
        mutate(lim = case_when( 
          x <1 ~5, 
          x>=1 & x <10 ~ 10, 
          x>=10 & x< 30 ~ 30, 
          x>=30 & x<60 ~ 60, 
          x>=60 & x<80 ~ 80, 
          x>80 & x<=100 ~100,
          x>100 & x<=1000 ~ (x%/%100 + 1)*100,
          x >1000 & x<10000 ~ (x%/%1000 + 1)*1000
        ))
    }
    
  }
  
  return(lim)
} 

library(readxl)
output_path <- "/Users/jjwu/Projects/Simplified-HCV-testing-model/Projects/POC_AU/Output"

Prev_dt <- read_excel(file.path(paste0(output_path, "/PrevInc_epi.xlsx")), sheet = "tempPrevRNA_setting")
Prev_dt <- Prev_dt%>%mutate(scenario = factor(scenario, levels = sce_level, labels = sce_label))
Prev_dt <- Prev_dt%>%mutate(setting = factor(setting, levels = c("commu", "prisons"), 
                                             labels = c("Community", "Prison")))

View(Prev_dt)

RNA_prev <- ggplot(Prev_dt%>%mutate(year = year + POC_AU$cabY - 1)%>%filter(setting%in%c("Community", "Prison") ), 
       aes(x = year, y = best)) + 
  geom_line(aes(colour = scenario, linetype = scenario), size = 1) + 
  facet_wrap(~ setting, scale ="free", ncol = 2 ) + 
  scale_color_manual(name = "Scenarios", values = col_pal ) + 
  scale_fill_manual(name = "Scenarios", values = col_pal ) + 
  scale_linetype_manual(name = "Scenarios", 
                        values = c("dashed", rep("solid", length(unique(Prev_dt$scenario)) - 1))) + 
  coord_cartesian(xlim = c(2021,2030)) +
  scale_x_continuous(expand = c(0, 0), limits =c(2021,2030) ,
                     breaks = seq(2021,2030, 
                                  by = 1)) + 
  theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
  theme(legend.key.size = unit(1,"line"),
        legend.direction = "vertical") +
  theme(legend.position = c(0.2, 0.15)) + 
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 1)) + 
  labs(x = "Year", y = "HCV RNA prevalence") 
  
ggsave(file=file.path(OutputFig, paste0("RNAprev_setting_maintext.png")), 
       RNA_prev, 
       width = 12, height = 8, bg = "white", dpi = 300)   



#### cost_saving plot #### 
# net savings
net_labels <- incre_catcost %>%
  group_by(scenario, sensitivity) %>%
  summarise(
    net_saving = sum(incre_best),
    label_y = sum(incre_best[incre_best > 0]) + 8000000,
    .groups = "drop"
  )


incre_plot_sens <- list()
for(i in unique(incre_catcost$sensitivity)){ 
  if(i == "Main analysis"){ 
    gtitle <- paste0("Main analysis: 50% discounted DAA cost" )
    
  }else if(i == "NP total program cost"){
    gtitle <- (paste0("National Program total costs" ))
  }else {
    gtitle <- i
  }
  incre_plot_sens[[i]] <- 
    ggplot(incre_catcost%>%filter(sensitivity == i), 
           aes(x = scenario, y = incre_best, fill = Categories)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = c( "grey30", "grey40","grey80")) + 
    theme_Publication(base_size = 16) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          legend.direction = "vertical") + 
    scale_y_continuous(limit = c(-40000000, 100000000), 
                       breaks = seq(-40000000, 100000000, 20000000),
                       labels = seq(-40000000, 100000000,  20000000)/1000000) + 
    labs(y = "Cost saving (discounted, millions)", x = "Scenarios") + 
    geom_text( aes(label = paste0(format(round(incre_best/1000000, 1), nsmall = 1), "m")),
               position = position_stack(vjust = 0.5), 
               size = 6)  + geom_hline(yintercept = 0, linetype ="dashed") + 
    geom_text(
      data = net_labels%>%filter(sensitivity == i),
      aes(
        x = scenario,
        y = label_y,
        label = paste0("Net saving: ", format(round(net_saving/1000000, 1),, nsmall = 1), "m")
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4
    )
    ggtitle(gtitle )
  
  
}






library(ggpubr)

library(cowplot)

incre_plot_sens[[1]] <- incre_plot_sens[[1]] + 
  theme(
    legend.key.size = unit(1.5, "cm"),
    legend.key.width = unit(1.5, "cm"),    # width of legend keys
    legend.key.height = unit(1, "cm"),     # height of legend keys
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.spacing.y = unit(0.5, "cm")     # spacing between legend items
  ) + theme(legend.position = "right") + 
  scale_y_continuous(limit = c(-40000000, 80000000), 
                     breaks = seq(-40000000, 80000000, 10000000),
                     labels = seq(-40000000, 80000000,  10000000)/1000000)

get_legend_manual <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

legend <- get_legend_manual(incre_plot_sens[[1]])

# Remove legends from all plots
plots_no_legend <- lapply(incre_plot_sens, function(p) {
  p + theme(legend.position = "none")
})

# Arrange: 3 plots + legend in 4th cell
incre_cost_sen <- ggarrange(
  plots_no_legend[[1]], plots_no_legend[[3]],
  plots_no_legend[[2]], legend,
  nrow = 2, ncol = 2
)

ggsave(file=file.path(OutputFig, paste0("incre_cost_sen.png")), 
       incre_cost_sen,  width = 20, height = 16, bg = "white", dpi = 600) 

for(i in 1: length(names(cost_disydaanocap_categories))){ 
  
  ggsave(file=file.path(OutputFig, paste0("incre_cost_sen_", names(cost_disydaanocap_categories)[i],".png")), 
         incre_plot_sens[[i]],  width = 14, height = 12, bg = "white", dpi = 300) 
}






ggsave(file=file.path(OutputFig, paste0("incre_cost_maintext.png")), 
       incre_plot_sens[[1]],  width = 14, height = 12, bg = "white", dpi = 300) 

incre_plot_sens[[1]] <- incre_plot_sens[[1]] + 
  theme(
    legend.key.size = unit(1.5, "cm"),
    legend.key.width = unit(1.5, "cm"),    # width of legend keys
    legend.key.height = unit(1, "cm"),     # height of legend keys
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.spacing.y = unit(0.5, "cm")     # spacing between legend items
  ) + theme(legend.position = "right")

get_legend_manual <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

legend <- get_legend_manual(incre_plot_sens[[1]])

# Remove legends from all plots
plots_no_legend <- lapply(incre_plot_sens, function(p) {
  p + theme(legend.position = "none")
})

incre_plot_sens[[1]] <- incre_plot_sens[[1]] + ggtitle("")
ggsave(file=file.path(OutputFig%>%dirname(), paste0("costsaving_category_maintext.png")), 
       incre_plot_sens[[1]],  width = 13, height = 8, bg = "white", dpi = 300) 




