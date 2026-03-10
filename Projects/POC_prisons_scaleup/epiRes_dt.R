
gc()
rm(list = ls())

project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Projects/Simplified-HCV-testing-model")

data_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                   project_name, sep = "")
output_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                   "POC_prisons scale-up", sep = "")
# Load useful libraries

library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")

Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
RdaFolder <- file.path(data_path, "02. Output")
OutputFolder <- file.path(output_path, "02. Output")
OutputFig <- file.path(paste0(OutputFolder, "/Figs/PrevInc"))
if (!dir.exists(OutputFig)) {
  dir.create(OutputFig, recursive = TRUE, showWarnings = FALSE)
}
# project specific code path 
Proj_code <- file.path(codefun_path, paste0("Projects/", project_name))



load(file.path(RdaFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0("Simulations", ".rda")))
source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Proj_code, "/model_timestep.R")) 
# simulation outcomes 
# Sce_sq: pre-national program scenario 
# Sce_np: national program scenarios 

#### epi outcomes ####  

all_scenarios <- c(list("Status quo" = Sce_sq), Sce_np)

# the flow here is only for epi 
epi_only <- TRUE

if(isTRUE(epi_only)){ 
  indicator_flow <- Sce_sq[!names(Sce_sq) %in% c("allPops", "newpop_tran", 
                                                 "newpop_tranState", "HCVdeathState",
                                                 "newDeathState", "death_hcv", 
                                                 "costPops", "QALYPops",
                                                 "costTestingAb", "costTestingAg", 
                                                 "costTestingPOCT", "costTreatment",
                                                 "costCured", "costRetreat",
                                                 "costTestingAb_sc", "costTestingAg_sc",
                                                 "costTestingPOCT_sc", "costTreatment_sc")]
} else { 
  indicator_flow <- Sce_sq[!names(Sce_sq) %in% c("allPops", "newpop_tran", 
                                                 "newpop_tranState", "HCVdeathState",
                                                 "newDeathState", "death_hcv")]
}

endY <- 100

# par_col <- c("best", paste0("set", seq(1,1000,1)))
par_col <- c("best")

Num_box <- list()
# par_Num_box <- list()
# get number in each component in each timestep 

pop_N <- list()

commu_N <- list()

prison_N <- list()

prisonPWID_N <- list()

# number in the each cascade box 
Num_diag <- list()

Num_diag_ab <- list()

Num_diag_Treated <- list() 

Num_chronic_cured <- list() 

Num_curInf <- list()

Num_dc <- list()

Num_hcc <- list()

Num_lt <- list()

Num_plt <- list()

# number in flows in each timestep
Sce_flow <- list() 

#### prev & inc in each population ####
# the prevalence and incidence only include chronic stage, excluding acute stage (a)

tempNOTInfected_subpop <- list()

tempChronic_subpop <- list()

tempPrev_subpop <- list()

tempNOTInfectedRNA_subpop <- list()

tempPrevRNA_subpop <- list()

HCVInc_subpop <- list()


for (name in names(all_scenarios)) { 
  Num_box[[name]] <- modres.t(POC_AU, all_scenarios[[name]], endYear = 100)%>%
    tibble::as_tibble() 
  #par_Num_box <- lapply(param_sq, function(x) modres.t(POC_AU, x, endYear = 100)%>%
  #                        tibble::as_tibble()%>%select(best))
  
  # for(i in 1:length(par_Num_box)){ 
  #  Num_box[[name]][, paste0("set",i)] <- par_Num_box[[i]]$best
  
  # }
  # Num_box[[name]] <- Num_box[[name]]%>%
  #  select(year,population, state, timestep, cascade, disease_prog,
  #         best, paste0("set", seq(1,1000,1)))
  
  # total N of all compartments in each timestep 
  
  
  pop_N[[name]] <- N_pop_sum(Num_box[[name]], 
                                     pop = NULL, param = NULL, name_parset = NULL)
  
  commu_N[[name]] <- N_pop_sum(Num_box[[name]], 
                                       pop = c("C_PWID", "C_fPWID"),param = NULL, name_parset = NULL)
  
  prison_N[[name]] <- N_pop_sum(Num_box[[name]], 
                                        pop = c("P_PWID", "P_fPWID", "P_nPWID"),param = NULL, name_parset = NULL)
  
  prisonPWID_N[[name]] <- N_pop_sum(Num_box[[name]], 
                                            pop = c("P_PWID", "P_fPWID"),param = NULL, name_parset = NULL)
  
  # diagnosis: compartments includes those cured from treatment: excluding those achieved cured at acute stage
  
  Num_diag[[name]] <- N_pop_casdisprog(Num_box[[name]], pop = NULL, 
                                               cas = c("diag_RNA", "treat", "treat_f"),
                                               disprog = c(POC_AU$progress_name)[-1], 
                                               param = NULL, name_parset = NULL)
  
  Num_diag_ab[[name]] <- N_pop_casdisprog(Num_box[[name]], pop = NULL, 
                                                  cas = c("diag_ab","diag_RNA", "treat", 
                                                          "treat_f", "cured"), 
                                                  disprog = c(POC_AU$progress_name)[-1],
                                                  param = NULL, name_parset = NULL)
  
  Num_diag_Treated[[name]] <- N_pop_casdisprog(Num_box[[name]], pop = NULL, 
                                                       cas = c("treat", "treat_f", "cured"),
                                                       disprog = c(POC_AU$progress_name)[-1],
                                                       param = NULL, name_parset = NULL)
  
  Num_chronic_cured[[name]] <- N_pop_casdisprog(Num_box[[name]], pop = NULL, 
                                                        cas = c("cured"),
                                                        disprog = c(POC_AU$progress_name)[-1],
                                                        param = NULL, name_parset = NULL)
  
  Num_curInf[[name]] <- N_pop_casdisprog(Num_box[[name]], pop = NULL, 
                                                 cas = c("undiag", "diag_ab", "diag_RNA", 
                                                         "treat", "treat_f"),
                                                 disprog = NULL,
                                                 param = NULL, name_parset = NULL)
  
  Num_dc[[name]] <- N_pop_casdisprog(Num_box[[name]], pop = NULL, cas = NULL,
                                             disprog = c("dc"), 
                                             param = NULL, name_parset = NULL)
  
  Num_hcc[[name]] <- N_pop_casdisprog(Num_box[[name]], pop = NULL, cas = NULL,
                                              disprog = c("hcc"),
                                              param = NULL, name_parset = NULL)
  
  Num_lt[[name]] <- N_pop_casdisprog(Num_box[[name]], pop = NULL, cas = NULL,
                                             disprog = c("lt"),
                                             param = NULL, name_parset = NULL)
  
  Num_plt[[name]] <- N_pop_casdisprog(Num_box[[name]], pop = NULL, cas = NULL,
                                              disprog = c("plt"),
                                              param = NULL, name_parset = NULL)
  
  
  
 
  }

for(name in names(all_scenarios)){ 
  Sce_flow[[name]] <- list()
  for(x in names(indicator_flow)){ 
    
    Sce_flow[[name]][[x]] <- modres.flow.t(POC_AU, all_scenarios[[name]], endYear = endY, 
                                           allp = x)
    # par_Sce_flow <- lapply(param_sq, 
    #                       function(x) lapply(names(indicator_flow), 
    #                                          function(y) modres.flow.t(POC_AU, x, 
    #                                                                    endYear = endY, 
    #                                                                    allp = y)%>%
    #                                            tibble::as_tibble()%>%select(best)))
    
    
    #for(i in 1:length(par_Sce_flow)){ 
    #  names(par_Sce_flow[[i]]) <- names(indicator_flow)
    #  for(x in names(indicator_flow)){
    #    Sce_flow[["Status quo"]][[x]][, paste0("set",i)] <- par_Sce_flow[[i]][[x]]$best
    #    
    #  }
    #}
    
    
  }
}



for(i in names(Num_box)){ 
  tempNOTInfected_subpop[[i]] <- Num_box[[i]]%>%
    filter(disease_prog!= "a")%>%
    filter(state == "s")%>%group_by(timestep, population)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep, population)
  
  tempChronic_subpop[[i]] <- Num_box[[i]]%>%filter(disease_prog!= "a")%>%
    group_by(timestep, population)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep, population) 
  
  
  # arrange order to align with other dts 
  tempPrev_subpop[[i]] <- cbind(timestep = pop_N[[i]]$timestep,
                                population = POC_AU$popNames,
                                as.data.frame(100*(pop_N[[i]][, par_col] - 
                                                     tempNOTInfected_subpop[[i]][ ,par_col])/ 
                                                pop_N[[i]][ ,par_col]))%>%
    tibble::as_tibble() 
  
  
  # RNA prevalence 
  tempNOTInfectedRNA_subpop[[i]] <- Num_box[[i]]%>%
    filter(cascade%in% c("s", "cured"))%>%group_by(timestep, population)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep, population)
  
  tempPrevRNA_subpop[[i]] <- 
    cbind(timestep = pop_N[[i]]$timestep,
          population = POC_AU$popNames,
          
          as.data.frame(100*(pop_N[[i]][, par_col] - 
                               tempNOTInfectedRNA_subpop[[i]][ ,par_col])/ 
                          pop_N[[i]][ ,par_col]))%>%
    tibble::as_tibble() 
  
  
  # incidence 
  HCVInc_subpop[[i]] <- cbind(timestep = pop_N[[i]]$timestep,
                              population = POC_AU$popNames,
                              as.data.frame(100*Sce_flow[[i]]$newInfections[, par_col] / 
                                              pop_N[[i]][ ,par_col]))%>%
    tibble::as_tibble() 
  
}  


#### HCV incidence by setting #### 



#### prev & inc by settings ####
tempNOTInfected_commu <- list()

tempNOTInfected_prison <- list()

tempNOTInfected_prisonPWID <- list()

tempPrev_setting <- list()

tempNOTInfectedRNA_commu <- list()

tempNOTInfectedRNA_prison <- list()

tempNOTInfectedRNA_prisonPWID <- list()

tempPrevRNA_setting <- list()

for(n in names(Num_box)){ 
  
  tempNOTInfected_commu[[n]] <- Num_box[[n]]%>%
    filter(population %in% c("C_PWID", "C_fPWID") & disease_prog == "s")%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep)
  
  tempPrev_setting[[n]][["commu"]] <- cbind(timestep = commu_N[[n]]$timestep,
                                            as.data.frame(100*(commu_N[[n]][, par_col] - 
                                                                 tempNOTInfected_commu[[n]][ ,par_col])/ 
                                                            commu_N[[n]][ ,par_col]))%>%tibble::as_tibble()
  
  
  tempNOTInfected_prison[[n]] <- Num_box[[n]]%>%
    filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID") & disease_prog == "s")%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep)
  
  tempPrev_setting[[n]][["prisons"]] <- cbind(timestep = prison_N[[n]]$timestep,
                                              as.data.frame(100*(prison_N[[n]][, par_col] - 
                                                                   tempNOTInfected_prison[[n]][ ,par_col])/ 
                                                              prison_N[[n]][ ,par_col]))%>%tibble::as_tibble()
  
  
  tempNOTInfected_prisonPWID[[n]] <- Num_box[[n]]%>%
    filter(population %in% c("P_PWID", "P_fPWID") & disease_prog == "s")%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep)
  
  
  tempPrev_setting[[n]][["prisonsPWID"]] <- 
    cbind(timestep = prisonPWID_N[[n]]$timestep,
          as.data.frame(100*(prisonPWID_N[[n]][, par_col] - 
                               tempNOTInfected_prisonPWID[[n]][ ,par_col])/ 
                          prisonPWID_N[[n]][, par_col]))%>%tibble::as_tibble() 
  
  
  # RNA prevalence 
  
  
  tempNOTInfectedRNA_commu[[n]] <- Num_box[[n]]%>%
    filter(cascade%in% c("s", "cured") & population %in% c("C_PWID", "C_fPWID"))%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep)
  
  tempPrevRNA_setting[[n]][["commu"]] <- cbind(timestep = commu_N[[n]]$timestep,
                                               as.data.frame(100*(commu_N[[n]][, par_col] - 
                                                                    tempNOTInfectedRNA_commu[[n]][ ,par_col])/ 
                                                               commu_N[[n]][ ,par_col]))%>%tibble::as_tibble()
  
  tempNOTInfectedRNA_prison[[n]] <- Num_box[[n]]%>%
    filter(cascade%in% c("s", "cured") & population %in% c("P_PWID", "P_fPWID", "P_nPWID"))%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep)
  
  tempPrevRNA_setting[[n]][["prisons"]] <- 
    cbind(timestep = prison_N[[n]]$timestep,
          as.data.frame(100*(prison_N[[n]][, par_col] - 
                               tempNOTInfectedRNA_prison[[n]][ ,par_col])/ 
                          prison_N[[n]][ ,par_col]))%>%tibble::as_tibble()
  
  tempNOTInfectedRNA_prisonPWID[[n]] <- Num_box[[n]]%>%
    filter(cascade%in% c("s", "cured") & population %in% c("P_PWID", "P_fPWID"))%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep)
  
  tempPrevRNA_setting[[n]][["prisonsPWID"]] <- 
    cbind(timestep = prisonPWID_N[[n]]$timestep,
          as.data.frame(100*(prisonPWID_N[[n]][, par_col] - 
                               tempNOTInfectedRNA_prisonPWID[[n]][ ,par_col])/ 
                          prisonPWID_N[[n]][ ,par_col]))%>%tibble::as_tibble()
  
}


#### incidence in settings ####
newInf_commu <- list()
newInf_prison <- list()
newInf_prisonPWID <- list()
HCVInc_setting <- list()
for(n in names(Sce_flow)){ 
  newInf_commu[[n]] <- Sce_flow[[n]]$newInfections%>%filter(population %in% c("C_PWID", "C_fPWID"))%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep)
  
  newInf_prison[[n]] <- Sce_flow[[n]]$newInfections%>%
    filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))%>%group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep)
  
  newInf_prisonPWID[[n]] <- Sce_flow[[n]]$newInfections%>%
    filter(population %in% c("P_PWID", "P_fPWID"))%>%group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%arrange(timestep)
  
  HCVInc_setting[[n]][["commu"]] <- cbind(timestep = commu_N[[n]]$timestep,
                                          as.data.frame(100*(newInf_commu[[n]][ , par_col]/ 
                                                               commu_N[[n]][ ,par_col])))%>%
    tibble::as_tibble()
  
  HCVInc_setting[[n]][["prison"]] <- cbind(timestep = prison_N[[n]]$timestep,
                                           as.data.frame(100*(newInf_prison[[n]][ , par_col]/ 
                                                                prison_N[[n]][ ,par_col])))%>%
    tibble::as_tibble()
  
  HCVInc_setting[[n]][["prisonPWID"]] <- cbind(timestep = prisonPWID_N[[n]]$timestep,
                                               as.data.frame(100*(newInf_prisonPWID[[n]][ , par_col]/ 
                                                                    prisonPWID_N[[n]][ ,par_col])))%>%
    tibble::as_tibble()
  
}

save(Num_box, pop_N, commu_N, prison_N, prisonPWID_N, 
     Num_diag, Num_diag_ab, Num_diag_Treated, 
     Num_chronic_cured, Num_curInf,
     Num_dc, Num_hcc, Num_lt, Num_plt, 
     Sce_flow, tempNOTInfected_subpop, 
     tempChronic_subpop, tempPrev_subpop,
     tempNOTInfectedRNA_subpop, tempPrevRNA_subpop, 
     HCVInc_subpop, 
     tempNOTInfected_commu, tempNOTInfected_prison, 
     tempNOTInfected_prisonPWID,
     tempPrev_setting, tempNOTInfectedRNA_commu, 
     tempNOTInfectedRNA_prison,
     tempNOTInfectedRNA_prisonPWID, tempPrevRNA_setting, 
     newInf_commu, newInf_prison, 
     newInf_prisonPWID, HCVInc_setting,
     # file = file.path(OutputFolder,paste0(project_name,"epiRes_timestep_sq" ,".rda"))
     file = file.path(OutputFolder,paste0("epiRes_timestep",".rda"))
     
)

Resflow_dt <- list()
Resflow_sc_dt <- list()
for(scn in names(Sce_flow)){ 
  Resflow_dt[[scn]] <- list(newInfections = Sce_flow[[scn]]$newInfections,
                            HCVdeath = Sce_flow[[scn]]$newHCVdeaths,
                            Treatment = Sce_flow[[scn]]$newTreatment,
                            Retreat = Sce_flow[[scn]]$newRetreat,
                            Testing_ab = Sce_flow[[scn]]$newTestingAb,
                            Testing_RNA = Sce_flow[[scn]]$newTestingAg,
                            Testing_POCT = Sce_flow[[scn]]$newTestingPOCT, 
                            Testing_ab_neg = Sce_flow[[scn]]$newTestingAb_neg,
                            Testing_RNA_neg = Sce_flow[[scn]]$newTestingAg_neg,
                            Testing_POCT_neg = Sce_flow[[scn]]$newTestingPOCT_neg,
                            Cured = Sce_flow[[scn]]$newCured, 
                            Treatment_sc = Sce_flow[[scn]]$newTreatment_sc,
                            Testing_ab_sc = Sce_flow[[scn]]$newTestingAb_sc,
                            Testing_RNA_sc = Sce_flow[[scn]]$newTestingAg_sc,
                            Testing_POCT_sc = Sce_flow[[scn]]$newTestingPOCT_sc, 
                            Testing_ab_sc_neg = Sce_flow[[scn]]$newTestingAb_sc_neg,
                            Testing_RNA_sc_neg = Sce_flow[[scn]]$newTestingAg_sc_neg,
                            Testing_POCT_sc_neg = Sce_flow[[scn]]$newTestingPOCT_sc_neg)
  
  
  Resflow_sc_dt[[scn]] <- list(Treatment_sc = Sce_flow[[scn]]$newTreatment_sc,
                               Testing_ab_sc = Sce_flow[[scn]]$newTestingAb_sc,
                               Testing_RNA_sc = Sce_flow[[scn]]$newTestingAg_sc,
                               Testing_POCT_sc = Sce_flow[[scn]]$newTestingPOCT_sc,
                               Testing_ab_sc_neg = Sce_flow[[scn]]$newTestingAb_sc_neg,
                               Testing_RNA_sc_neg = Sce_flow[[scn]]$newTestingAg_sc_neg,
                               Testing_POCT_sc_neg = Sce_flow[[scn]]$newTestingPOCT_sc_neg
  )
  }

save(Num_box, Resflow_dt, Resflow_sc_dt,
     # file = file.path(OutputFolder,paste0(project_name,"epiRes_timestep_sq" ,".rda"))
     file = file.path(OutputFolder,paste0("Res_dt",".rda")))





#===============================================================================
#                              function to generate plots
#===============================================================================

library(paletteer)

scale_colour_paletteer_d("nationalparkcolors::Acadia")
PrevInc_plot <- function(pj, dt, obdt =NULL, xlimits, UI = NULL){ 
  if(length(unique(dt$scenario)) == 2){ 
    col_pal <- c("#000000", "#E69F00")
    
  } 
  else{col_pal <- c(paletteer_d("nationalparkcolors::Acadia"))
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

tempPrevRNA_setting_bind <- lapply(tempPrevRNA_setting, function(x) bind_rows(x,
                                                                               .id = 'setting'))

tempPrevRNA_setting_bind <- bind_rows(tempPrevRNA_setting_bind, .id = 'scenario')

tempPrevRNA_setting_bind_mid <- tempPrevRNA_setting_bind%>%filter(round(timestep*12)%%12 ==6)%>%
  mutate(year = timestep%/%1)


tempPrevRNA_setting_bind_mid <- tempPrevRNA_setting_bind_mid%>%filter(setting%in% c("commu", "prisons"))%>%
  mutate(setting = factor(setting, levels = c("commu", "prisons"),
                          labels = c("Community", "Prison")))%>%
  mutate(scenario = factor(scenario, 
                           levels = unique(tempPrevRNA_setting_bind_mid$scenario),
                           labels = c("Status Quo", "Prison_testing_I",
                                      "Prison_testing_II", 
                                      "Prison_testing_III",
                                      "Program sustained",
                                      "Program scale-up")))


col_pal <- c(paletteer_d("nationalparkcolors::Acadia"))
#### RNA prevalence ####
RNA_prev <- list()
lab_name <- c("Status Quo", "Prison_testing_I",
              "Prison_testing_II", 
              "Prison_testing_III",
              "Program sustained",
              "Program scale-up")

for(i in lab_name){
  RNA_prev[[i]] <- ggplot(tempPrevRNA_setting_bind_mid%>%
                            mutate(year = year + POC_AU$cabY - 1)%>%
           filter(setting%in%c("Community", "Prison"))%>%
                    filter(scenario == i), aes(x = year, y = best)) + 
    geom_line(aes(colour = scenario, linetype = scenario), size = 0.5) + 
    facet_wrap(~ setting, scale ="free", ncol = 2 ) + 
    scale_color_manual(name = "Scenarios", values = col_pal ) + 
    scale_fill_manual(name = "Scenarios", values = col_pal ) + 
    scale_linetype_manual(name = "Scenarios", 
                          values = c("solid")) + 
    coord_cartesian(xlim = c(2021,2030)) +
    scale_x_continuous(expand = c(0, 0), limits =c(2021,2030) ,
                       breaks = seq(2021,2030, 
                                    by = 1)) + 
    theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
    theme(legend.key.size = unit(1,"line")) + 
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 1)) + 
    labs(x = "Year", y = "HCV RNA prevalence") 
  
  ggsave(file=file.path(OutputFig, paste0("PrevRNA_setting_",i ,".png")), 
         RNA_prev[[i]], 
         width = 12, height = 8, bg = "white", dpi = 300) 
}








#### result aggregate ####
gc()
rm(list = ls())

project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Projects/Simplified-HCV-testing-model")

data_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                   project_name, sep = "")
output_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                     "POC_prisons scale-up", sep = "")
# Load useful libraries

library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")

Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
RdaFolder <- file.path(data_path, "02. Output")
OutputFolder <- file.path(output_path, "02. Output")
OutputFig <- file.path(paste0(OutputFolder, "/Figs/PrevInc"))
if (!dir.exists(OutputFig)) {
  dir.create(OutputFig, recursive = TRUE, showWarnings = FALSE)
}
# project specific code path 
Proj_code <- file.path(codefun_path, paste0("Projects/", project_name))



load(file.path(RdaFolder, paste0(project_name, ".rda")))

source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Proj_code, "/model_timestep.R")) 
load(file.path(OutputFolder, paste0("Res_dt", ".rda")))



View(Resflow_dt$Prison_testing_I$newInfections)
Resflow_year_pop <- list()

for(x in names(Resflow_dt)){ 
  Resflow_year_pop[[x]] <- Resflow_dt[[x]]
  
}

names(Resflow_year_pop) <- names(Resflow_dt)
par_col <- c("best")
for(i in names(Resflow_year_pop)){
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

for(x in names(Resflow_dt)){ 
  Resflow_sc_year_pop[[x]] <- Resflow_sc_dt[[x]]
  
}

names(Resflow_sc_year_pop) <- names(Resflow_dt)



for(i in names(Resflow_sc_year_pop)){
  for(indic in names(Resflow_sc_year_pop[[1]])){
    Resflow_sc_year_pop[[i]][[indic]] <- Resflow_sc_year_pop[[i]][[indic]]%>%
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
    )%>%as.data.frame()
  
  Resflow_year_all[[i]][["Tot_Treatment_np"]] <- 
    cbind(year = Resflow_year_all[[i]][["Treatment"]]$year, 
          dplyr::bind_cols(
                             Resflow_year_all[[i]][["Treatment_sc"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_all[[i]][["Tot_Testing_ab"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_ab"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Testing_ab"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_ab_neg"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_ab_sc"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_ab_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_all[[i]][["Tot_Testing_ab_np"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_ab"]]$year, 
          dplyr::bind_cols(
                             Resflow_year_all[[i]][["Testing_ab_sc"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_ab_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_all[[i]][["Tot_Testing_RNA"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_RNA"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Testing_RNA"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_RNA_neg"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_RNA_sc"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_RNA_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_all[[i]][["Tot_Testing_RNA_np"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_RNA"]]$year, 
          dplyr::bind_cols(
                             Resflow_year_all[[i]][["Testing_RNA_sc"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_RNA_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_all[[i]][["Tot_Testing_POCT"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_POCT"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Testing_POCT"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_POCT_neg"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_POCT_sc"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_POCT_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_all[[i]][["Tot_Testing_POCT_np"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_POCT"]]$year, 
          dplyr::bind_cols( 
                             Resflow_year_all[[i]][["Testing_POCT_sc"]][, par_col] + 
                             Resflow_year_all[[i]][["Testing_POCT_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
 
}

for(i in names(Resflow_year_all)){
  Resflow_year_pop[[i]][["Tot_Treatment"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Treatment"]]$year, 
          population = Resflow_year_pop[[i]][["Treatment"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Treatment"]][, par_col] + 
                             Resflow_year_pop[[i]][["Retreat"]][, par_col] + 
                             Resflow_year_pop[[i]][["Treatment_sc"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_pop[[i]][["Tot_Treatment_np"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Treatment"]]$year, 
          population = Resflow_year_pop[[i]][["Treatment"]]$population,
          dplyr::bind_cols( 
                             Resflow_year_pop[[i]][["Treatment_sc"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_pop[[i]][["Tot_Testing_ab"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_ab"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_ab"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Testing_ab"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_ab_neg"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_ab_sc"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_ab_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_pop[[i]][["Tot_Testing_ab_np"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_ab"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_ab"]]$population,
          dplyr::bind_cols(
                             Resflow_year_pop[[i]][["Testing_ab_sc"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_ab_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
  
  Resflow_year_pop[[i]][["Tot_Testing_RNA"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_RNA"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_RNA"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Testing_RNA"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_RNA_neg"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_RNA_sc"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_RNA_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_pop[[i]][["Tot_Testing_RNA_np"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_RNA"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_RNA"]]$population,
          dplyr::bind_cols(
                             Resflow_year_pop[[i]][["Testing_RNA_sc"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_RNA_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_pop[[i]][["Tot_Testing_POCT"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_POCT"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_POCT"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Testing_POCT"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_POCT_neg"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_POCT_sc"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_POCT_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_pop[[i]][["Tot_Testing_POCT_np"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_POCT"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_POCT"]]$population,
          dplyr::bind_cols(
                             Resflow_year_pop[[i]][["Testing_POCT_sc"]][, par_col] + 
                             Resflow_year_pop[[i]][["Testing_POCT_sc_neg"]][, par_col])
    )%>%as.data.frame()
  
 
  
}



for(i in  names(Resflow_year_all)){ 
  
  Resflow_year_pop[[i]][["Tot_Testing"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_POCT"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_POCT"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Tot_Testing_ab"]][, par_col] + 
                                    Resflow_year_pop[[i]][["Tot_Testing_RNA"]][, par_col] + 
                                    Resflow_year_pop[[i]][["Tot_Testing_POCT"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_pop[[i]][["Tot_Testing_np"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_POCT"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_POCT"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Tot_Testing_ab_np"]][, par_col] + 
                             Resflow_year_pop[[i]][["Tot_Testing_RNA_np"]][, par_col] + 
                             Resflow_year_pop[[i]][["Tot_Testing_POCT_np"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_pop[[i]][["Tot_screened"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_POCT"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_POCT"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Tot_Testing_ab"]][, par_col] +
                                    Resflow_year_pop[[i]][["Tot_Testing_POCT"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_pop[[i]][["Tot_screened_np"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_POCT"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_POCT"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Tot_Testing_ab_np"]][, par_col] +
                             Resflow_year_pop[[i]][["Tot_Testing_POCT_np"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_all[[i]][["Tot_Testing"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_POCT"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Tot_Testing_ab"]][, par_col] + 
                             Resflow_year_all[[i]][["Tot_Testing_RNA"]][, par_col] + 
                             Resflow_year_all[[i]][["Tot_Testing_POCT"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_all[[i]][["Tot_Testing_np"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_POCT"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Tot_Testing_ab_np"]][, par_col] + 
                             Resflow_year_all[[i]][["Tot_Testing_RNA_np"]][, par_col] + 
                             Resflow_year_all[[i]][["Tot_Testing_POCT_np"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_all[[i]][["Tot_screened"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_POCT"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Tot_Testing_ab"]][, par_col] + 
                             Resflow_year_all[[i]][["Tot_Testing_POCT"]][, par_col])
    )%>%as.data.frame()
  
  Resflow_year_all[[i]][["Tot_screened_np"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_POCT"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Tot_Testing_ab_np"]][, par_col] + 
                             Resflow_year_all[[i]][["Tot_Testing_POCT_np"]][, par_col])
    )%>%as.data.frame()
  
  
}

for(i in  names(Resflow_year_all)){ 
  colnames(Resflow_year_all[[i]][["Tot_Testing"]]) <- c("year", "best")
  colnames(Resflow_year_all[[i]][["Tot_screened"]]) <- c("year", "best")
  colnames(Resflow_year_pop[[i]][["Tot_Testing"]]) <- c("year", "population", "best")
  colnames(Resflow_year_pop[[i]][["Tot_screened"]]) <- c("year", "population", "best")
  colnames(Resflow_year_all[[i]][["Tot_Testing_np"]]) <- c("year", "best")
  colnames(Resflow_year_all[[i]][["Tot_screened_np"]]) <- c("year", "best")
  colnames(Resflow_year_pop[[i]][["Tot_Testing_np"]]) <- c("year", "population", "best")
  colnames(Resflow_year_pop[[i]][["Tot_screened_np"]]) <- c("year", "population", "best")
  
  
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

#### generating plots

##### new infection by setting ####
lab_name <- c("Status Quo", "Prison_testing_I",
              "Prison_testing_II", 
              "Prison_testing_III",
              "Program sustained",
              "Program scale-up")
newinf_setting <- list()
names(Resflow_year_setting) <- lab_name 
col_pal <- c(paletteer_d("nationalparkcolors::Acadia"))
for(i in lab_name){
  newinf_setting[[i]] <- ggplot(Resflow_year_setting[[i]]$newInfections,
                                aes(x = year, y = best)) + 
    geom_bar(stat = "identity", size = 0.5) + 
    geom_text(data = Resflow_year_setting[[i]]$newInfections,
              aes(x = year, y = best, label = round(best, 0)),
              size = 3, vjust = -0.5) + 
    facet_wrap(~ population, scale ="free", ncol = 2 ) + 
    scale_color_manual(name = "Scenarios", values = col_pal ) + 
    scale_fill_manual(name = "Scenarios", values = col_pal ) + 
    scale_linetype_manual(name = "Scenarios", 
                          values = c("solid")) + 
    coord_cartesian(xlim = c(2021,2031)) +
    scale_x_continuous(expand = c(0, 0), limits = c(2021, 2031),
                       breaks = seq(2022, 2030, by = 1)) + 
    theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
    theme(legend.key.size = unit(1,"line")) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 5000), breaks = seq(0, 5000, 500)) + 
    labs(x = "Year", y = "Number of HCV new infections", title = i) 

  ggsave(file=file.path(OutputFig, paste0("newinf_setting_",i ,".png")), 
         newinf_setting[[i]], 
         width = 12, height = 8, bg = "white", dpi = 300) 
}

treatment_setting <- list()
for(i in lab_name){
  treatment_setting[[i]] <- ggplot(Resflow_year_setting[[i]]$Tot_Treatment,
                                aes(x = year, y = best)) + 
    geom_bar(stat = "identity", size = 0.5) + 
    geom_text(data = Resflow_year_setting[[i]]$Tot_Treatment,
              aes(x = year, y = best, label = round(best, 0)),
              size = 3, vjust = -0.5) + 
    
    facet_wrap(~ population, scale ="free", ncol = 2 ) + 
    scale_color_manual(name = "Scenarios", values = col_pal ) + 
    scale_fill_manual(name = "Scenarios", values = col_pal ) + 
    scale_linetype_manual(name = "Scenarios", 
                          values = c("solid")) + 
    coord_cartesian(xlim = c(2021,2031)) +
    scale_x_continuous(expand = c(0, 0), limits = c(2021, 2031),
                       breaks = seq(2022, 2030, by = 1)) + 
    theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
    theme(legend.key.size = unit(1,"line")) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 6000), breaks = seq(0, 6000, 1000)) + 
    labs(x = "Year", y = "Number of treatment initiation", title = i) 
    

  ggsave(file=file.path(OutputFig, paste0("treatment_setting_",i ,".png")), 
         treatment_setting[[i]], 
         width = 12, height = 8, bg = "white", dpi = 300) 
}

#### Tot_testing 
tot_testing_setting <- list() 
for(i in 1:length(names(Resflow_year_setting))){ 
  
  p <- ggplot() +
    geom_area(data = Resflow_year_setting[[i]]$Tot_Testing,
              aes(x = year, y = best, fill = "Total"), 
              alpha = 0.6, size = 0.2, colour = "black")
  
  # only add National Program layers if not Status Quo
  if(lab_name[i] != "Status Quo"){
    p <- p +
      geom_area(data = Resflow_year_setting[[i]]$Tot_Testing_np,
                aes(x = year, y = best, fill = "National Program"),
                alpha = 0.6, size = 0.2, colour = "black") +
      geom_text(data = Resflow_year_setting[[i]]$Tot_Testing_np,
                aes(x = year, y = best, label = round(best, 0)),
                size = 3, vjust = -0.5)
  }
  
  p <- p +
    facet_wrap(~ population) +
    theme_Publication_facet() + 
    scale_fill_manual(name = "Testing pathway",
                      values = c("National Program" = "#F8766D",
                                 "Total"            = "#00BFC4")) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.6, colour = NA))) +
    scale_x_continuous(limits = c(2022, 2030), breaks = seq(2022, 2030, 1)) +
    scale_y_continuous(limits = c(0, 80000), breaks = seq(0, 80000, 10000),
                       labels = seq(0, 80000, 10000) / 1000) +
    labs(x = "Year", y = "Number of total tests (thousands)", title = lab_name[i])
  
  tot_testing_setting[[i]] <- p
  
  ggsave(file = file.path(OutputFig, paste0("tot_testing_setting_", lab_name[i], ".png")), 
         tot_testing_setting[[i]], 
         width = 12, height = 8, bg = "white", dpi = 300) 
}





#### combine figures ####
list.files("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/POC_prisons scale-up/02. Output/Figs/PrevInc")

library(patchwork)
library(png)


library(grid)
library(htmltools)

# define order
scenarios <- c("Status Quo", "Prison_testing_I", "Prison_testing_II", 
               "Prison_testing_III", "Program sustained", "Program scale-up")

indicators <- c("tot_testing_setting", "PrevRNA_setting", 
                "newinf_setting", "treatment_setting")

# helper to load png as ggplot grob
load_png <- function(path){
  img <- readPNG(path)
  ggplot() + 
    annotation_custom(rasterGrob(img, width = unit(1,"npc"), height = unit(1,"npc"))) + 
    theme_void()
}

# build list of plots: row = indicator, col = scenario
plot_list <- list()
for(ind in indicators){
  for(sc in scenarios){
    fname <- file.path(OutputFig, paste0(ind, "_", sc, ".png"))
    plot_list[[paste(ind, sc, sep = "_")]] <- load_png(fname)
  }
}

# assemble with patchwork: 4 rows x 6 cols
panel <- wrap_plots(plot_list, nrow = 4, ncol = 6) +
  plot_annotation(
    title = "HCV testing model outcomes by scenario",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(file.path(OutputFig%>%dirname(), "panel_all.png"), panel,
       width = 30, height = 20, dpi = 300, bg = "white") 




# build HTML table
indicators <- c("tot_testing_setting", "PrevRNA_setting", 
                "newinf_setting", "treatment_setting")

indicator_labels <- c("Total Testing", "RNA Prevalence", "New Infections", "Treatment")

header_row <- tags$tr(
  tags$th(""),  # empty corner
  lapply(scenarios, function(sc) tags$th(sc, style = "text-align:center; padding:8px; font-size:13px"))
)

body_rows <- lapply(seq_along(indicators), function(i){
  tags$tr(
    tags$td(indicator_labels[i], 
            style = "font-weight:bold; writing-mode:vertical-rl; text-align:center; padding:8px; font-size:13px"),
    lapply(scenarios, function(sc){
      fname <- file.path(OutputFig, paste0(indicators[i], "_", sc, ".png"))
      tags$td(
        tags$img(src = fname, width = "100%"),
        style = "padding:4px"
      )
    })
  )
})

html_out <- tags$html(
  tags$head(
    tags$title("HCV Model Panel"),
    tags$style("body { font-family: Arial; } table { width: 100%; border-collapse: collapse; } th, td { border: 1px solid #ddd; }")
  ),
  tags$body(
    tags$h2("HCV testing model outcomes by scenario", style = "text-align:center"),
    tags$table(header_row, tags$tbody(body_rows))
  )
)

save_html(html_out, file = file.path(OutputFig%>%dirname(), "panel_all.html"))
