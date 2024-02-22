# result summary for each scenario

rm(list = ls())

project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Documents/Simplified-HCV-testing-model")

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
OutputFolder <- file.path(data_path, "02. Output")
# project specific code path 
Proj_code <- file.path(codefun_path, paste0("projects/", project_name))



load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "Simulations.rda")))

source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Proj_code, "/model_timestep.R")) 
# simulation outcomes 
# Sce_sq: pre-national program scenario 
# Sce_np: national program scenarios 

#### epi outcomes ####  
#### HCV diagnosis and treatment #### 
# Number of people living with diagnosed RNA 
# Number of people living with diagnosed ab
# Number of people who are treated: on treatment & treatment failure 
# Number of people living with SVR: chronic_cured  
# number of people current HCV infection 


indicator_flow <- Sce_sq[!names(Sce_sq)%in% c("allPops", "newpop_tran", 
                                              "newpop_tranState", "HCVdeathState",
                                              "newDeathState", "death_hcv", 
                                              "costPops", "QALYPops")]
endY <- 100


Num_box <- list()

# get number in each component in each timestep 

Num_box[["Status quo"]] <- modres.t(POC_AU, Sce_sq, endYear = 100)%>%as.data.frame() 

for(i in names(Sce_np)){ 
  
  Num_box[[i]] <- modres.t(POC_AU, Sce_np[[i]], endYear = 100)%>%as.data.frame() 
  
  }

pop_N <- list()

commu_N <- list()

prison_N <- list()

prisonPWID_N <- list()

# total N of all compartments in each timestep 

for(i in names(Num_box)){ 
  pop_N[[i]] <- N_pop_sum(Num_box[[i]], pop = NULL)
  
  commu_N[[i]] <- N_pop_sum(Num_box[[i]], 
                            pop = c("C_PWID", "C_fPWID"))
  
  prison_N[[i]] <- N_pop_sum(Num_box[[i]], 
                             pop = c("P_PWID", "P_fPWID", "P_nPWID"))
  
  prisonPWID_N[[i]] <- N_pop_sum(Num_box[[i]], 
                                 pop = c("P_PWID", "P_fPWID"))
}
  


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
for(i in names(Num_box)){
  Num_diag[[i]] <- N_pop_casdisprog(Num_box[[i]], pop = NULL, 
                                    cas = c("diag_RNA", "treat", "treat_f"),
                                    disprog = c(POC_AU$progress_name)[-1])
  
  Num_diag_ab[[i]] <- N_pop_casdisprog(Num_box[[i]], pop = NULL, 
                                       cas = c("diag_ab","diag_RNA", "treat", 
                                               "treat_f", "cured"), 
                                       disprog = c(POC_AU$progress_name)[-1])
  
  Num_diag_Treated[[i]] <- N_pop_casdisprog(Num_box[[i]], pop = NULL, 
                                            cas = c("treat", "treat_f", "cured"),
                                            disprog = c(POC_AU$progress_name)[-1])
  
  Num_chronic_cured[[i]] <- N_pop_casdisprog(Num_box[[i]], pop = NULL, 
                                             cas = c("cured"),
                                             disprog = c(POC_AU$progress_name)[-1])
  
  Num_curInf[[i]] <- N_pop_casdisprog(Num_box[[i]], pop = NULL, 
                                      cas = c("undiag", "diag_ab", "diag_RNA", 
                                              "treat", "treat_f"),
                                      disprog = NULL)
  
  Num_dc[[i]] <- N_pop_casdisprog(Num_box[[i]], pop = NULL, cas = NULL,
                                  disprog = c("dc"))
  
  Num_hcc[[i]] <- N_pop_casdisprog(Num_box[[i]], pop = NULL, cas = NULL,
                                   disprog = c("hcc"))
  
  Num_lt[[i]] <- N_pop_casdisprog(Num_box[[i]], pop = NULL, cas = NULL,
                                  disprog = c("lt"))
  
  Num_plt[[i]] <- N_pop_casdisprog(Num_box[[i]], pop = NULL, cas = NULL,
                                   disprog = c("plt"))
  
}

# number in flows in each timestep
Sce_flow <- list() 


for(x in names(indicator_flow)){ 
  
  Sce_flow[["Status quo"]][[x]] <- modres.flow.t(POC_AU, Sce_sq, endYear = endY, 
                                                 allp = x)
}

for(i in names(Sce_np)){
  for(x in names(indicator_flow)){ 
    Sce_flow[[i]][[x]] <- modres.flow.t(POC_AU, Sce_np[[i]], endYear = endY, 
                                        allp = x)
  
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
    summarise(best = sum(best))%>%arrange(timestep, population)
  
  tempChronic_subpop[[i]] <- Num_box[[i]]%>%filter(disease_prog!= "a")%>%
    group_by(timestep, population)%>%
    summarise(best = sum(best))%>%arrange(timestep, population) 
  
  
  # arrange order to align with other dts 
  tempPrev_subpop[[i]] <- cbind(timestep = pop_N[[i]]$timestep,
                                population = POC_AU$popNames,
                                as.data.frame(100*(pop_N[[i]][, -c(1,2)] - 
                                                     tempNOTInfected_subpop[[i]][ ,-c(1,2)])/ 
                                                pop_N[[i]][ ,-c(1,2)]))%>%
    tibble::as_tibble() 
  
  
  # RNA prevalence 
  tempNOTInfectedRNA_subpop[[i]] <- Num_box[[i]]%>%
    filter(cascade%in% c("s", "cured"))%>%group_by(timestep, population)%>%
    summarise(best = sum(best))%>%arrange(timestep, population)
  
  tempPrevRNA_subpop[[i]] <- 
    cbind(timestep = pop_N[[i]]$timestep,
          population = POC_AU$popNames,
          
          as.data.frame(100*(pop_N[[i]][, -c(1,2)] - 
                               tempNOTInfectedRNA_subpop[[i]][ ,-c(1,2)])/ 
                          pop_N[[i]][ ,-c(1,2)]))%>%
    tibble::as_tibble() 
  
  
  # incidence 
  HCVInc_subpop[[i]] <- cbind(timestep = pop_N[[i]]$timestep,
                              population = POC_AU$popNames,
                              as.data.frame(100*Sce_flow[[i]]$newInfections[, "best"] / 
                                              pop_N[[i]][ ,-c(1,2)]))%>%
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
    summarise(best = sum(best))%>%arrange(timestep)
  
  tempPrev_setting[[n]][["commu"]] <- cbind(timestep = commu_N[[n]]$timestep,
                                            as.data.frame(100*(commu_N[[n]][, -c(1)] - 
                                                                 tempNOTInfected_commu[[n]][ ,-c(1)])/ 
                                                            commu_N[[n]][ ,-c(1)]))%>%tibble::as_tibble()
  
  
  tempNOTInfected_prison[[n]] <- Num_box[[n]]%>%
    filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID") & disease_prog == "s")%>%
    group_by(timestep)%>%
    summarise(best = sum(best))%>%arrange(timestep)
  
  tempPrev_setting[[n]][["prisons"]] <- cbind(timestep = prison_N[[n]]$timestep,
                                              as.data.frame(100*(prison_N[[n]][, -c(1)] - 
                                                                   tempNOTInfected_prison[[n]][ ,-c(1)])/ 
                                                              prison_N[[n]][ ,-c(1)]))%>%tibble::as_tibble()
  
  
  tempNOTInfected_prisonPWID[[n]] <- Num_box[[n]]%>%
    filter(population %in% c("P_PWID", "P_fPWID") & disease_prog == "s")%>%
    group_by(timestep)%>%
    summarise(best = sum(best))%>%arrange(timestep)
  
  
  tempPrev_setting[[n]][["prisonsPWID"]] <- 
    cbind(timestep = prisonPWID_N[[n]]$timestep,
          as.data.frame(100*(prisonPWID_N[[n]][, -c(1)] - 
                               tempNOTInfected_prisonPWID[[n]][ ,-c(1)])/ 
                          prisonPWID_N[[n]][, -c(1)]))%>%tibble::as_tibble() 
  
  
  # RNA prevalence 
  
  
  tempNOTInfectedRNA_commu[[n]] <- Num_box[[n]]%>%
    filter(cascade%in% c("s", "cured") & population %in% c("C_PWID", "C_fPWID"))%>%
    group_by(timestep)%>%
    summarise(best = sum(best))%>%arrange(timestep)
  
  tempPrevRNA_setting[[n]][["commu"]] <- cbind(timestep = commu_N[[n]]$timestep,
                                               as.data.frame(100*(commu_N[[n]][, -c(1)] - 
                                                                    tempNOTInfectedRNA_commu[[n]][ ,-c(1)])/ 
                                                               commu_N[[n]][ ,-c(1)]))%>%tibble::as_tibble()
  
  tempNOTInfectedRNA_prison[[n]] <- Num_box[[n]]%>%
    filter(cascade%in% c("s", "cured") & population %in% c("P_PWID", "P_fPWID", "P_nPWID"))%>%
    group_by(timestep)%>%
    summarise(best = sum(best))%>%arrange(timestep)
  
  tempPrevRNA_setting[[n]][["prisons"]] <- 
    cbind(timestep = prison_N[[n]]$timestep,
          as.data.frame(100*(prison_N[[n]][, -c(1)] - 
                               tempNOTInfectedRNA_prison[[n]][ ,-c(1)])/ 
                          prison_N[[n]][ ,-c(1)]))%>%tibble::as_tibble()
  
  tempNOTInfectedRNA_prisonPWID[[n]] <- Num_box[[n]]%>%
    filter(cascade%in% c("s", "cured") & population %in% c("P_PWID", "P_fPWID"))%>%
    group_by(timestep)%>%
    summarise(best = sum(best))%>%arrange(timestep)
  
  tempPrevRNA_setting[[n]][["prisonsPWID"]] <- 
    cbind(timestep = prisonPWID_N[[n]]$timestep,
          as.data.frame(100*(prisonPWID_N[[n]][, -c(1)] - 
                               tempNOTInfectedRNA_prisonPWID[[n]][ ,-c(1)])/ 
                          prisonPWID_N[[n]][ ,-c(1)]))%>%tibble::as_tibble()
  
  }

#### incidence in settings ####
newInf_commu <- list()
newInf_prison <- list()
newInf_prisonPWID <- list()
HCVInc_setting <- list()
for(n in names(Sce_flow)){ 
  newInf_commu[[n]] <- Sce_flow[[n]]$newInfections%>%filter(population %in% c("C_PWID", "C_fPWID"))%>%
    group_by(timestep)%>%
    summarise_at("best", sum)%>%arrange(timestep)
  
  newInf_prison[[n]] <- Sce_flow[[n]]$newInfections%>%
    filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))%>%group_by(timestep)%>%
    summarise_at("best", sum)%>%arrange(timestep)
  
  newInf_prisonPWID[[n]] <- Sce_flow[[n]]$newInfections%>%
    filter(population %in% c("P_PWID", "P_fPWID"))%>%group_by(timestep)%>%
    summarise_at("best", sum)%>%arrange(timestep)
  
  HCVInc_setting[[n]][["commu"]] <- cbind(timestep = commu_N[[n]]$timestep,
                                          as.data.frame(100*(newInf_commu[[n]][ , "best"]/ 
                                                               commu_N[[n]][ ,-c(1)])))%>%
    tibble::as_tibble()
  
  HCVInc_setting[[n]][["prison"]] <- cbind(timestep = prison_N[[n]]$timestep,
                                          as.data.frame(100*(newInf_prison[[n]][ , "best"]/ 
                                                               prison_N[[n]][ ,-c(1)])))%>%
    tibble::as_tibble()
  
  HCVInc_setting[[n]][["prisonPWID"]] <- cbind(timestep = prisonPWID_N[[n]]$timestep,
                                           as.data.frame(100*(newInf_prisonPWID[[n]][ , "best"]/ 
                                                                prisonPWID_N[[n]][ ,-c(1)])))%>%
    tibble::as_tibble()
  
  }


save(Num_box, pop_N, commu_N, prison_N, prisonPWID_N, 
     Num_diag, Num_diag_ab, Num_diag_Treated, Num_chronic_cured, Num_curInf,
     Num_dc, Num_hcc, Num_lt, Num_plt, 
     Sce_flow, tempNOTInfected_subpop, tempChronic_subpop, tempPrev_subpop,
     tempNOTInfectedRNA_subpop, tempPrevRNA_subpop, HCVInc_subpop, 
     tempNOTInfected_commu, tempNOTInfected_prison, tempNOTInfected_prisonPWID,
     tempPrev_setting, tempNOTInfectedRNA_commu, tempNOTInfectedRNA_prison,
     tempNOTInfectedRNA_prisonPWID, tempPrevRNA_setting, 
     newInf_commu, newInf_prison, newInf_prisonPWID, HCVInc_setting, 
     file = file.path(OutputFolder,
                      paste0(project_name,"epiRes_timestep" ,".rda")))


#### tidy up cost #### 
##### total cost in each time step ##### 

# cost attached to each compartment
cost_box <- list()
cost_box_sum <- list()
cost_box[["Status quo"]] <- modres.t(POC_AU, Sce_sq, endYear = 100, 
                                     allp = "costPops")%>%as.data.frame() 
cost_box_sum[["Status quo"]] <- N_pop_sum(cost_box[["Status quo"]], pop = NULL)

for(i in names(Sce_np)){ 
  
  cost_box[[i]] <- modres.t(POC_AU, Sce_np[[i]], endYear = 100,
                           allp = "costPops")%>%as.data.frame()
  
  # sum all compartment at each timestep 
  cost_box_sum[[i]] <- N_pop_sum(cost_box[[i]], pop = NULL)
  
}

# cost dataset: costs attached to compartments and costs attaced to each flow
# flow numbers: numbers transit between compartments at each timestep
Rescost_dt <- list()
Resflow_dt <- list()
for(i in names(cost_box)){ 
  
  Rescost_dt[[i]] <- cbind(cost_box_sum[[i]], 
                           cost_ab = Sce_flow[[i]]$costTestingAb$best, 
                           cost_RNA = Sce_flow[[i]]$costTestingAg$best,
                           cost_POCT = Sce_flow[[i]]$costTestingPOCT$best,
                           cost_Treatment = Sce_flow[[i]]$costTreatment$best,
                           cost_Retreat = Sce_flow[[i]]$costRetreat$best,
                           cost_Cured = Sce_flow[[i]]$costCured$best)
  
  Rescost_dt[[i]] <- Rescost_dt[[i]]%>%
    mutate(cost_compartment = best,
           cost_total = cost_compartment + cost_ab + cost_RNA + cost_POCT +
                            cost_Treatment + cost_Retreat + cost_Cured)%>%
    select(timestep, population, cost_total, cost_compartment, 
           cost_ab, cost_RNA, cost_POCT, 
           cost_Treatment, cost_Retreat, cost_Cured)
  
  
  
  Resflow_dt[[i]] <- cbind(year = Sce_flow[[i]]$newInfections$year,
                           timestep = Sce_flow[[i]]$newInfections$timestep, 
                           population = Sce_flow[[i]]$newInfections$population,
                           newInfections = Sce_flow[[i]]$newInfections$best,
                           HCVdeath = Sce_flow[[i]]$newHCVdeaths$best,
                           Treatment = Sce_flow[[i]]$newTreatment$best,
                           Retreat = Sce_flow[[i]]$newRetreat$best,
                           Testing_ab = Sce_flow[[i]]$newTestingAb$best,
                           Testing_RNA = Sce_flow[[i]]$newTestingAg$best,
                           Testing_POCT = Sce_flow[[i]]$newTestingPOCT$best, 
                           Cured = Sce_flow[[i]]$newCured$best)

}

save(Num_box, Resflow_dt, Rescost_dt, 
     file = file.path(OutputFolder,
                      paste0(project_name,"Res_dt" ,".rda")))

