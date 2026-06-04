#### Results summary for each NP scenarios ####

# with sensitivity analysis 

# simulation outcomes 
# Sce_sq: pre-national program scenario 
# Sce_np: national program scenarios 
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

Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
RDAFolder <- file.path(data_path, "02. Output" )
OutputFolder <- file.path(codefun_path, "Projects/POC_AU/Output")
OutputFig <- file.path(codefun_path, "Projects/POC_AU/Figs")
# project specific code path 
Proj_code <- file.path(codefun_path, paste0("projects/", project_name))
load(file.path(RDAFolder, paste0(project_name, ".rda")))
source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Proj_code, "/model_timestep.R")) 

sce_name <- c("dfList_NP_2024", "dfList_NPPhaseII",
              "dfList_NPPhaseIII_A", "dfList_NPPhaseIII_B")

cost_types <- c("fixednvariable", "total",  "DAAcost_reduchalf")
for(scn in sce_name){
  
  for(cost_type in cost_types){ 
    
    load(file.path(OutputFolder, paste0(project_name, "param_sc_",scn ,"_", cost_type, ".rda")))
    
    load(file.path(RDAFolder, paste0(project_name, "Simulations_", cost_type,".rda")))
    
    indicator_flow <- Sce_sq[!names(Sce_sq)%in% c("allPops", "newpop_tran", 
                                                  "newpop_tranState", "HCVdeathState",
                                                  "newDeathState", "death_hcv", 
                                                  "costPops", "QALYPops")]
    endY <- 100
    
    par_col <- c("best", paste0("set", seq(1,1000,1)))
    rm(Sce_sq)
    gc()
    
    Num_box <- list()
    
    # get number in each component in each timestep 
    
    Num_box[[scn]] <- modres.t(POC_AU, Sce_np[[scn]], endYear = 100)%>%
      tibble::as_tibble() 
    par_Num_box <- list() 
    par_Num_box <- lapply(param_scenario, function(x) modres.t(POC_AU, x, endYear = 100)%>%
                            tibble::as_tibble()%>%select(best))
    
    for(i in 1:length(par_Num_box)){ 
      Num_box[[scn]][, paste0("set",i)] <- par_Num_box[[i]]$best
      
    }
    Num_box[[scn]] <- Num_box[[scn]]%>%
      select(year,population, state, timestep, cascade, disease_prog,
             best, paste0("set", seq(1,1000,1)))
    
    
    
    
    pop_N <- list()
    
    commu_N <- list()
    
    prison_N <- list()
    
    prisonPWID_N <- list()
    
    # total N of all compartments in each timestep 
    
    pop_N[[scn]] <- N_pop_sum(Num_box[[scn]], 
                                     pop = NULL, param = "y", name_parset = par_col)
    
    commu_N[[scn]] <- N_pop_sum(Num_box[[scn]], 
                                       pop = c("C_PWID", "C_fPWID"),param = "y", name_parset = par_col)
    
    prison_N[[scn]] <- N_pop_sum(Num_box[[scn]], 
                                        pop = c("P_PWID", "P_fPWID", "P_nPWID"),param = "y", name_parset = par_col)
    
    prisonPWID_N[[scn]] <- N_pop_sum(Num_box[[scn]], 
                                            pop = c("P_PWID", "P_fPWID"),param = "y", name_parset = par_col)
    
    
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
    # diagnosis: compartments includes those cured from treatment: excluding those achieved cured at acute stage
    
    Num_diag[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, 
                                               cas = c("diag_RNA", "treat", "treat_f"),
                                               disprog = c(POC_AU$progress_name)[-1], 
                                               param = "y", name_parset = par_col)
    
    Num_diag_ab[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, 
                                                  cas = c("diag_ab","diag_RNA", "treat", 
                                                          "treat_f", "cured"), 
                                                  disprog = c(POC_AU$progress_name)[-1],
                                                  param = "y", name_parset = par_col)
    
    Num_diag_Treated[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, 
                                                       cas = c("treat", "treat_f", "cured"),
                                                       disprog = c(POC_AU$progress_name)[-1],
                                                       param = "y", name_parset = par_col)
    
    Num_chronic_cured[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, 
                                                        cas = c("cured"),
                                                        disprog = c(POC_AU$progress_name)[-1],
                                                        param = "y", name_parset = par_col)
    
    Num_curInf[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, 
                                                 cas = c("undiag", "diag_ab", "diag_RNA", 
                                                         "treat", "treat_f"),
                                                 disprog = NULL,
                                                 param = "y", name_parset = par_col)
    
    Num_dc[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, cas = NULL,
                                             disprog = c("dc"), 
                                             param = "y", name_parset = par_col)
    
    Num_hcc[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, cas = NULL,
                                              disprog = c("hcc"),
                                              param = "y", name_parset = par_col)
    
    Num_lt[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, cas = NULL,
                                             disprog = c("lt"),
                                             param = "y", name_parset = par_col)
    
    Num_plt[[scn]] <- N_pop_casdisprog(Num_box[[scn]], pop = NULL, cas = NULL,
                                              disprog = c("plt"),
                                              param = "y", name_parset = par_col)
    
    
    
    
    
    
    # number in flows in each timestep
    Sce_flow <- list() 
    
    par_Sce_flow <- list()
    for(x in names(indicator_flow)){ 
      
      Sce_flow[[scn]][[x]] <- modres.flow.t(POC_AU, Sce_np[[scn]], endYear = endY, 
                                                   allp = x)
      
    }
    
    par_Sce_flow <- lapply(param_scenario, 
                           function(x) lapply(names(indicator_flow), 
                                              function(y) modres.flow.t(POC_AU, x, 
                                                                        endYear = endY, 
                                                                        allp = y)%>%
                                                tibble::as_tibble()%>%select(best)))
    
    
    for(i in 1:length(par_Sce_flow)){ 
      names(par_Sce_flow[[i]]) <- names(indicator_flow)
      for(x in names(indicator_flow)){
        Sce_flow[[scn]][[x]][, paste0("set",i)] <- par_Sce_flow[[i]][[x]]$best
        
      }
    }
    
    
    #### prev & inc in each population ####
    # the prevalence and incidence only include chronic stage, excluding acute stage (a)
    
    tempNOTInfected_subpop <- list()
    
    tempChronic_subpop <- list()
    
    tempPrev_subpop <- list()
    
    tempNOTInfectedRNA_subpop <- list()
    
    tempPrevRNA_subpop <- list()
    
    HCVInc_subpop <- list()
    
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
      newInf_commu[[n]] <- Sce_flow[[n]]$newInfections%>%
        filter(population %in% c("C_PWID", "C_fPWID"))%>%
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
    
    
    save(Num_box, pop_N, commu_N, 
         prison_N, prisonPWID_N, 
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
         file = file.path(OutputFolder,
                          paste0(project_name,"epiRes_timestep", scn,"_", cost_type,".rda")))
    #                 paste0(project_name,"epiRes_timestep", sc_name,"totalcost" ,".rda")))
    
    #### tidy up cost #### 
    ##### total cost in each time step ##### 
    
    # cost attached to each compartment
    cost_box <- list()
    cost_box_sum <- list()
    par_cost_box <- list()
    cost_box[[scn]] <- modres.t(POC_AU, Sce_np[[scn]], endYear = 100, 
                                       allp = "costPops")%>%as.data.frame() 
    
    par_cost_box <- lapply(param_scenario, function(x) modres.t(POC_AU, x, endYear = 100, allp = "costPops")%>%
                             tibble::as_tibble()%>%select(best))
    
    for(i in 1:length(par_cost_box)){ 
      cost_box[[scn]][, paste0("set",i)] <- par_cost_box[[i]]$best
      
    }
    
    cost_box[[scn]] <- cost_box[[scn]]%>%
      select(year,population, state, timestep, cascade, disease_prog,
             best, paste0("set", seq(1,1000,1)))
    
    
    cost_box_sum[[scn]] <- N_pop_sum(cost_box[[scn]], 
                                            pop = NULL, param = "y", name_parset = par_col)
    
    # QALY
    QALY_box <- list()
    QALY_box_sum <- list()
    par_QALY_box <- list()
    QALY_box[[scn]] <- modres.t(POC_AU, Sce_np[[scn]], endYear = 100, allp = "QALYPops")%>%
      tibble::as_tibble()
    
    par_QALY_box <- lapply(param_scenario, function(x) modres.t(POC_AU, x, endYear = 100, allp = "QALYPops")%>%
                             tibble::as_tibble()%>%select(best))
    
    for(i in 1:length(par_cost_box)){ 
      QALY_box[[scn]][, paste0("set",i)] <- par_QALY_box[[i]]$best
      
    }
    
    QALY_box[[scn]] <- QALY_box[[scn]]%>%
      select(year,population, state, timestep, cascade, disease_prog,
             best, paste0("set", seq(1,1000,1)))
    
    
    QALY_box_sum[[scn]] <- N_pop_sum(QALY_box[[scn]], 
                                            pop = NULL, param = "y", name_parset = par_col)
    
    
    # cost dataset: costs attached to compartments and costs attaced to each flow
    # flow numbers: numbers transit between compartments at each timestep
    Rescost_dt <- list()
    Resflow_dt <- list()
    costqaly_name <- c("cost_box_sum" ,names(Sce_flow[[1]])[27:36])
    
    Rescost_dt[[scn]] <- list(cost_compartment = cost_box_sum[[scn]],
                                     cost_ab = Sce_flow[[scn]]$costTestingAb, 
                                     cost_RNA = Sce_flow[[scn]]$costTestingAg,
                                     cost_POCT = Sce_flow[[scn]]$costTestingPOCT,
                                     cost_Treatment = Sce_flow[[scn]]$costTreatment,
                                     cost_Retreat = Sce_flow[[scn]]$costRetreat,
                                     cost_Cured = Sce_flow[[scn]]$costCured,
                                     cost_ab_sc = Sce_flow[[scn]]$costTestingAb_sc, 
                                     cost_RNA_sc = Sce_flow[[scn]]$costTestingAg_sc,
                                     cost_POCT_sc = Sce_flow[[scn]]$costTestingPOCT_sc,
                                     cost_Treatment_sc = Sce_flow[[scn]]$costTreatment_sc,
                                     QALY_compartment = QALY_box_sum[[scn]])
    
    Rescost_dt[[scn]][["cost_total"]] <- Rescost_dt[[scn]]$cost_ab
    for(i in par_col){ 
      Rescost_dt[[scn]][["cost_total"]][, i] <- 
        Rescost_dt[[scn]]$cost_compartment[,i] + Rescost_dt[[scn]]$cost_ab[, i] + 
        Rescost_dt[[scn]]$cost_RNA[, i] + Rescost_dt[[scn]]$cost_POCT[, i] + 
        Rescost_dt[[scn]]$cost_Treatment[, i] + Rescost_dt[[scn]]$cost_Retreat[, i] + 
        Rescost_dt[[scn]]$cost_Cured[, i]
    }
    
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
    
    Resflow_sc_dt <- list()
    Resflow_sc_dt[[scn]] <- list(Treatment_sc = Sce_flow[[scn]]$newTreatment_sc,
                                        Testing_ab_sc = Sce_flow[[scn]]$newTestingAb_sc,
                                        Testing_RNA_sc = Sce_flow[[scn]]$newTestingAg_sc,
                                        Testing_POCT_sc = Sce_flow[[scn]]$newTestingPOCT_sc,
                                        Testing_ab_sc_neg = Sce_flow[[scn]]$newTestingAb_sc_neg,
                                        Testing_RNA_sc_neg = Sce_flow[[scn]]$newTestingAg_sc_neg,
                                        Testing_POCT_sc_neg = Sce_flow[[scn]]$newTestingPOCT_sc_neg
    )
    
    
    
    save(Num_box, Resflow_dt, Resflow_sc_dt, 
         Rescost_dt, 
         file = file.path(OutputFolder,
                          #               paste0(project_name,"Res_dt_",scn,".rda")))
                          paste0(project_name,"Res_dt_",scn,"_", cost_type ,".rda")))
    
    
    rm(Num_box, par_Num_box, param_scenario, Sce_np, 
       pop_N, commu_N, prison_N, prisonPWID_N, 
       Num_diag, Num_diag_ab, Num_diag_Treated, Num_chronic_cured, Num_curInf,
       Num_dc, Num_hcc, Num_lt, Num_plt, 
       Sce_flow, tempNOTInfected_subpop, tempChronic_subpop, tempPrev_subpop,
       tempNOTInfectedRNA_subpop, tempPrevRNA_subpop, HCVInc_subpop,
       tempNOTInfected_commu, tempNOTInfected_prison, tempNOTInfected_prisonPWID,
       tempPrev_setting, tempNOTInfectedRNA_commu,tempNOTInfectedRNA_prison,
       tempNOTInfectedRNA_prisonPWID,tempPrevRNA_setting,
       newInf_commu, newInf_prison, newInf_prisonPWID, HCVInc_setting,
       cost_box, cost_box_sum,par_cost_box, QALY_box, QALY_box_sum, par_QALY_box,
       Resflow_dt, Resflow_sc_dt, Rescost_dt)
    gc()
    
    
  }
  
}

