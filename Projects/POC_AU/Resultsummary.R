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

load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "cali.rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_timev.rda")))
load(file.path(OutputFolder, paste0(project_name, "Sce_sq.rda")))
load(file.path(OutputFolder, paste0(project_name, "Sce_np.rda")))
load(file.path(OutputFolder, paste0(project_name, "S_NP_test.rda")))
load(file.path(OutputFolder, paste0(project_name, "S_NPscale_test.rda")))

load(file.path(OutputFolder, paste0(project_name, "coverage_sq.rda")))
load(file.path(OutputFolder, paste0(project_name, "coverage_np.rda")))

load(file.path(OutputFolder, paste0(project_name, "cost_np.rda")))

source(file.path(Rcode, "/Functions/HCV_model.R"))
source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R"))
source(file.path(Rcode, "/Functions/check_steady.R"))
source(file.path(codefun_path,paste0("Projects/", project_name, "/cost_model.R")))

Sce_epi <- list("status quo" = Sce_sq, 
                 "National program" = Sce_np)
#### HCV diagnosis and treatment #### 
# Number of people living with diagnosed RNA 
# Number of people living with diagnosed ab
# Number of people who are treated: on treatment & treatment failure 
# Number of people living with SVR: chronic_cured  
# number of people current HCV infection 


indicator_flow <- Sce_epi$`status quo`[!names(Sce_epi$`status quo`)%in%
                                         c("allPops", "newpop_tran", "newpop_tranState",
                                           "HCVdeathState",
                                           "newDeathState", "death_hcv", 
                                           "costPops", "QALYPops")]
endY <- 100



Num_component <- list()

# get number in each component
for(n in names(Sce_epi)){ 
  Num_component[[n]] <- popResults_MidYear(POC_AU, Sce_epi[[n]], 
                     Population = POC_AU$popNames,
                     Disease_prog = POC_AU$progress_name, 
                     Cascade = c("s",POC_AU$cascade_name), 
                     param = NULL,
                     endYear = 100)%>%
    as.data.frame() 
}

pop_N <- list()
commu_N <- list()
prison_N <- list()
prisonPWID_N <- list()
for(n in names(Sce_epi)){ 
  pop_N[[n]] <- Num_component[[n]]%>%group_by(year, population)%>%
    summarise_at("best", sum)%>%arrange(year)
    
  commu_N[[n]] <- pop_N[[n]]%>%
    filter(population %in% c("C_PWID", "C_fPWID"))%>%
    summarise_at("best", sum)%>%arrange(year)
  
  prison_N[[n]] <- pop_N[[n]]%>%
    filter(population %in% c("P_PWID", "P_fPWID","P_nPWID"))%>%
    summarise_at("best", sum)%>%arrange(year)
  
  prisonPWID_N[[n]] <- pop_N[[n]]%>%
    filter(population %in% c("P_PWID", "P_fPWID"))%>%
    summarise_at("best", sum)%>%arrange(year)
  
}

Sce_flow_sub  <- list()
Sce_flow_all <- list()
Sce_flow_setting <- list()

Num_diag <- list()
Num_diag_ab <- list()
Num_Treated <- list()
Num_chronic_cured <- list()
Num_curInf <- list()

Num_diag_all <- list()
Num_diag_ab_all <- list()
Num_Treated_all <- list()
Num_chronic_cured_all <- list()
Num_curInf_all <- list()

Num_diag_setting <- list()
Num_diag_ab_setting <- list()
Num_Treated_setting <- list()
Num_chronic_cured_setting <- list()
Num_curInf_setting <- list()

# advanced liver stage
Num_dc <- list()
Num_hcc <- list()
Num_lt <- list()

Num_dc_all <- list()
Num_hcc_all <- list()
Num_lt_all <- list()

Num_dc_setting <- list()
Num_hcc_setting <- list()
Num_lt_setting <- list()

for(n in names(Sce_epi)){ 
  for(x in names(indicator_flow)){ 
    
    Sce_flow_sub[[n]][[x]] <- indicatorResults(POC_AU, Sce_epi[[n]], x, 
                                               pop=POC_AU$popNames,
                                               paramR = NULL, range = NULL,
                                               endY = endY)%>%
      mutate(scenario = n)
 
    Sce_flow_setting[["commu"]][[n]][[x]] <- 
      Sce_flow_sub[[n]][[x]]%>%filter(population %in% c("C_PWID", "C_fPWID"))%>%
      group_by(year)%>%summarise_at("best", sum)%>%
      mutate(scenario = n)%>%arrange(year)
    
    Sce_flow_setting[["prisons"]][[n]][[x]] <- 
      Sce_flow_sub[[n]][[x]]%>%
      filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))%>%
      group_by(year)%>%summarise_at("best", sum)%>%
      mutate(scenario = n)%>%arrange(year)
    
    Sce_flow_setting[["prisonsPWID"]][[n]][[x]] <- 
      Sce_flow_sub[[n]][[x]]%>%
      filter(population %in% c("P_PWID", "P_fPWID"))%>%
      group_by(year)%>%summarise_at("best", sum)%>%
      mutate(scenario = n)%>%arrange(year)
    
    Sce_flow_all[[n]][[x]] <- Sce_flow_sub[[n]][[x]]%>%
      group_by(year)%>%summarise_at("best", sum)%>%
      mutate(scenario = n)
  }
  
  Num_diag[[n]] <- Num_component[[n]]%>%
    filter(cascade%in%c("diag_RNA", "treat", "treat_f"))%>%
    group_by(year, population)%>%summarize(best = sum(best))%>%
    mutate(scenario = n)
  
  Num_diag_all[[n]] <- Num_diag[[n]]%>%
    group_by(year)%>%summarise_at("best", sum)%>%
    mutate(scenario = n)
  
  Num_diag_ab[[n]] <- Num_component[[n]]%>%
    filter(cascade%in%c("diag_ab","diag_RNA", "treat", "treat_f"))%>%
    group_by(year, population)%>%summarize(best = sum(best))%>%
    mutate(scenario = n)
  
  Num_diag_ab_all[[n]] <- Num_diag_ab[[n]]%>%
    group_by(year)%>%summarise_at("best", sum)%>%
    mutate(scenario = n)
  
  Num_Treated[[n]] <- Num_component[[n]]%>%
    filter(cascade%in%c("treat", "treat_f"))%>%
    group_by(year, population)%>%summarize(best = sum(best))%>%
    mutate(scenario = n)
  
  Num_Treated_all[[n]] <- Num_Treated[[n]]%>%
    group_by(year)%>%summarise_at("best", sum)%>%
    mutate(scenario = n)
  
  Num_chronic_cured[[n]] <- Num_component[[n]]%>%
    filter(cascade%in%c("cured"))%>%
    group_by(year, population)%>%summarize(best = sum(best))%>%
    mutate(scenario = n)
  
  Num_chronic_cured_all[[n]] <- Num_chronic_cured[[n]]%>%
    group_by(year)%>%summarise_at("best", sum)%>%
    mutate(scenario = n)
  
  Num_curInf[[n]] <- Num_component[[n]]%>%
    filter(cascade%in%c("undiag", "diag_ab",  "diag_RNA", "treat",  
                        "treat_f"))%>%
    group_by(year, population)%>%summarize(best = sum(best))%>%
    mutate(scenario = n)
  
  Num_curInf_all[[n]] <- Num_curInf[[n]]%>%
    group_by(year)%>%summarise_at("best", sum)%>%
    mutate(scenario = n)
  
  Num_dc[[n]] <- Num_component[[n]]%>%
    filter(disease_prog%in%c("dc"))%>%
    group_by(year, population)%>%summarize(best = sum(best))%>%
    mutate(scenario = n)
  
  Num_dc_all[[n]] <- Num_dc[[n]]%>%
    group_by(year)%>%summarise_at("best", sum)%>%
    mutate(scenario = n)
  
  Num_hcc[[n]] <- Num_component[[n]]%>%
    filter(disease_prog%in%c("hcc"))%>%
    group_by(year, population)%>%summarize(best = sum(best))%>%
    mutate(scenario = n)
  
  Num_hcc_all[[n]] <- Num_hcc[[n]]%>%
    group_by(year)%>%summarise_at("best", sum)%>%
    mutate(scenario = n)
  
  Num_lt[[n]] <- Num_component[[n]]%>%
    filter(disease_prog%in%c("lt"))%>%
    group_by(year, population)%>%summarize(best = sum(best))%>%
    mutate(scenario = n)
  
  Num_lt_all[[n]] <- Num_lt[[n]]%>%
    group_by(year)%>%summarise_at("best", sum)%>%
    mutate(scenario = n) 
} 
setting_lab <- c("commu", "prisons", "prisonsPWID")
popfil_lab <- list("commu" = c("C_PWID", "C_fPWID"),
                   "prisons" = c("P_PWID", "P_fPWID", "P_nPWID"),
                   "prisonsPWID" = c("P_PWID", "P_fPWID"))

for(s in setting_lab){ 
  for(n in names(Sce_epi)){ 
    Num_diag_setting[[s]][[n]] <- Num_diag[[n]]%>%
      filter(population%in%popfil_lab[[s]])%>%
      group_by(year)%>%summarize(best = sum(best))%>%
      mutate(scenario = n)
    
    Num_diag_ab_setting[[s]][[n]] <- Num_diag_ab[[n]]%>%
      filter(population%in%popfil_lab[[s]])%>%
      group_by(year)%>%summarize(best = sum(best))%>%
      mutate(scenario = n)
    
    Num_Treated_setting[[s]][[n]] <- Num_Treated[[n]]%>%
      filter(population%in%popfil_lab[[s]])%>%
      group_by(year)%>%summarize(best = sum(best))%>%
      mutate(scenario = n)
    
    Num_chronic_cured_setting[[s]][[n]] <- Num_chronic_cured[[n]]%>%
      filter(population%in%popfil_lab[[s]])%>%
      group_by(year)%>%summarize(best = sum(best))%>%
      mutate(scenario = n)
    
    Num_curInf_setting[[s]][[n]] <- Num_curInf[[n]]%>%
      filter(population%in%popfil_lab[[s]])%>%
      group_by(year)%>%summarize(best = sum(best))%>%
      mutate(scenario = n)
    
    Num_dc_setting[[s]][[n]] <- Num_dc[[n]]%>%
      filter(population%in%popfil_lab[[s]])%>%
      group_by(year)%>%summarize(best = sum(best))%>%
      mutate(scenario = n)
    
    Num_hcc_setting[[s]][[n]] <- Num_hcc[[n]]%>%
      filter(population%in%popfil_lab[[s]])%>%
      group_by(year)%>%summarize(best = sum(best))%>%
      mutate(scenario = n)
    
    Num_lt_setting[[s]][[n]] <- Num_lt[[n]]%>%
      filter(population%in%popfil_lab[[s]])%>%
      group_by(year)%>%summarize(best = sum(best))%>%
      mutate(scenario = n)
    
    
  }
  
  }


# prev & inc 
tempNOTInfected_subpop <- list()
tempChronic_subpop <- list()
tempPrev_subpop <- list()
tempPrev_setting <- list()

tempNOTInfected_commu <- list()
tempNOTInfected_prison <- list()
tempNOTInfected_prisonPWID <- list()

tempNOTInfectedRNA_subpop <- list()
tempNOTInfectedRNA_commu <- list()
tempNOTInfectedRNA_prison <- list()
tempNOTInfectedRNA_prisonPWID <- list()

tempPrevRNA_subpop <- list()
tempPrevRNA_setting <- list()

HCVInc_subpop <- list()
HCVInc_setting <- list()

for(n in names(Sce_epi)){   
  # Prevalence 
  tempNOTInfected_subpop[[n]] <- Num_component[[n]]%>%filter(disease_prog!= "a")%>%
    filter(state == "s")%>%group_by(year, population)%>%
    summarise(best = sum(best))%>%arrange(year, population)
  
  tempChronic_subpop[[n]] <- Num_component[[n]]%>%filter(disease_prog!= "a")%>%
    group_by(year, population)%>%
    summarise(best = sum(best))%>%arrange(year, population)
  
  # arrange order to align with other dts 
  tempPrev_subpop[[n]] <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                           population = POC_AU$popNames,
                           
                           as.data.frame(100*(pop_N[[n]][, -c(1,2)] - 
                                                tempNOTInfected_subpop[[n]][ ,-c(1,2)])/ 
                                           pop_N[[n]][ ,-c(1,2)]))%>%
    tibble::as_tibble() 
  

  tempNOTInfected_commu[[n]] <- Num_component[[n]]%>%
    filter(population %in% c("C_PWID", "C_fPWID") & disease_prog == "s")%>%
    group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  tempPrev_setting[[n]][["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                       as.data.frame(100*(commu_N[[n]][, -c(1)] - 
                                                            tempNOTInfected_commu[[n]][ ,-c(1)])/ 
                                                       commu_N[[n]][ ,-c(1)]))%>%tibble::as_tibble()
  
  # prison 
  ## PWID + former PWID + nonPWID 

  
  tempNOTInfected_prison[[n]] <- Num_component[[n]]%>%
    filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID") & disease_prog == "s")%>%
    group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  
  tempPrev_setting[[n]][["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                         as.data.frame(100*(prison_N[[n]][, -c(1)] - 
                                                              tempNOTInfected_prison[[n]][ ,-c(1)])/ 
                                                         prison_N[[n]][ ,-c(1)]))%>%tibble::as_tibble()
  
  
  # prison_PWID experienced  
  ## PWID + former PWID + nonPWID 
  
  
  tempNOTInfected_prisonPWID[[n]] <- Num_component[[n]]%>%
    filter(population %in% c("P_PWID", "P_fPWID") & disease_prog == "s")%>%
    group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  
  tempPrev_setting[[n]][["prisonsPWID"]] <- 
    cbind(year = seq(POC_AU$startYear , endY-1 ,1),
          as.data.frame(100*(prisonPWID_N[[n]][, -c(1)] - 
                               tempNOTInfected_prisonPWID[[n]][ ,-c(1)])/ 
                          prisonPWID_N[[n]][, -c(1)]))%>%tibble::as_tibble()%>%
    mutate()
  
  # RNA prevalence 
  tempNOTInfectedRNA_subpop[[n]] <- Num_component[[n]]%>%
    filter(cascade%in% c("s", "cured"))%>%group_by(year, population)%>%
    summarise(best = sum(best))%>%arrange(year, population)
  
  tempPrevRNA_subpop[[n]] <- 
    cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
          population = POC_AU$popNames,
          
          as.data.frame(100*(pop_N[[n]][, -c(1,2)] - 
                               tempNOTInfectedRNA_subpop[[n]][ ,-c(1,2)])/ 
                          pop_N[[n]][ ,-c(1,2)]))%>%
    tibble::as_tibble()
  
  tempNOTInfectedRNA_commu[[n]] <- Num_component[[n]]%>%
    filter(cascade%in% c("s", "cured") & population %in% c("C_PWID", "C_fPWID"))%>%
    group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  tempPrevRNA_setting[[n]][["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                          as.data.frame(100*(commu_N[[n]][, -c(1)] - 
                                                               tempNOTInfectedRNA_commu[[n]][ ,-c(1)])/ 
                                                          commu_N[[n]][ ,-c(1)]))%>%tibble::as_tibble()
  
  tempNOTInfectedRNA_prison[[n]] <- Num_component[[n]]%>%
    filter(cascade%in% c("s", "cured") & population %in% c("P_PWID", "P_fPWID", "P_nPWID"))%>%
    group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  tempPrevRNA_setting[[n]][["prisons"]] <- 
    cbind(year = seq(POC_AU$startYear , endY-1 ,1),
          as.data.frame(100*(prison_N[[n]][, -c(1)] - 
                               tempNOTInfectedRNA_prison[[n]][ ,-c(1)])/ 
                          prison_N[[n]][ ,-c(1)]))%>%tibble::as_tibble()
  
  tempNOTInfectedRNA_prisonPWID[[n]] <- Num_component[[n]]%>%
    filter(cascade%in% c("s", "cured") & population %in% c("P_PWID", "P_fPWID"))%>%
    group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  tempPrevRNA_setting[[n]][["prisonsPWID"]] <- 
    cbind(year = seq(POC_AU$startYear , endY-1 ,1),
          as.data.frame(100*(prisonPWID_N[[n]][, -c(1)] - 
                               tempNOTInfectedRNA_prisonPWID[[n]][ ,-c(1)])/ 
                          prisonPWID_N[[n]][ ,-c(1)]))%>%tibble::as_tibble()
  
  # incidence 
  HCVInc_subpop[[n]] <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                    each = POC_AU$npops),
                         population = POC_AU$popNames,
                         as.data.frame(100*Sce_flow_sub[[n]]$newInfections[, "best"] / 
                                         pop_N[[n]][ ,-c(1,2)]))%>%
    tibble::as_tibble() 
  
  HCVInc_setting[[n]][["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                     as.data.frame(100*(Sce_flow_setting[["commu"]][[n]]$newInfections[ , "best"]/ 
                                                          commu_N[[n]][ ,-c(1)])))%>%
    tibble::as_tibble()
  
  HCVInc_setting[[n]][["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                          as.data.frame(100*(Sce_flow_setting[["prisons"]][[n]]$newInfections[ , "best"]/ 
                                                               prison_N[[n]][ ,-c(1)])))%>%
    tibble::as_tibble()
  
  HCVInc_setting[[n]][["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                            as.data.frame(100*(Sce_flow_setting[["prisonsPWID"]][[n]]$newInfections[ , "best"]/ 
                                                                 prisonPWID_N[[n]][ ,-c(1)])))%>%
    tibble::as_tibble()
  
    
  
  }

save(Sce_flow_sub, Sce_flow_all, Sce_flow_setting,
     Num_component, 
     Num_diag, Num_diag_ab, Num_Treated, Num_chronic_cured, Num_curInf, 
     Num_diag_all, Num_diag_ab_all, Num_Treated_all, Num_chronic_cured_all, Num_curInf_all, 
     Num_diag_setting, Num_diag_ab_setting, Num_Treated_setting, Num_chronic_cured_setting, Num_curInf_setting,
     Num_dc, Num_hcc, Num_lt, 
     Num_dc_all, Num_hcc_all, Num_lt_all, 
     Num_dc_setting, Num_hcc_setting, Num_lt_setting,
     tempPrev_subpop,
     tempPrev_setting,
     tempPrevRNA_subpop,
     tempPrevRNA_setting,
     HCVInc_subpop,
     HCVInc_setting,
     
     file = file.path(OutputFolder ,
                                  paste0(project_name,"ResSum" ,".rda")))

endY_plot <- 35
x_plot <- list()
for(n in names(Sce_flow_sub)){ 
  
  x_plot[[n]] <- indicatorPlot(POC_AU,Sce_flow_all[[n]]$newTreatment, 
                               ylabel = "N",
                               xlimits = c(POC_AU$startYear, 
                                           (POC_AU$startYear+endY_plot), 5),
                               calibration_Y = POC_AU$cabY,
                               rangeun = NULL, 
                               groupPlot = NULL, 
                               facetPlot = NULL,
                               observationData = NULL, 
                               simulateYear = POC_AU$simY) + 
    ggtitle("number of diagnosed HCV") +
    scale_y_continuous(limits = c(0,60000))
  }

ggplot(data= Num_component$`status quo`) + 
  geom_line(aes(x = year, y = best, group = population)) + 
  facet_wrap(.~ state, scale = "free")
