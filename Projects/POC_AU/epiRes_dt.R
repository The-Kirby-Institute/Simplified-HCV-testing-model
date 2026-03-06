
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
  Num_box[[name]] <- modres.t(POC_AU, Sce_sq, endYear = 100)%>%
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




