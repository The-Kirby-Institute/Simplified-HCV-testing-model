# this is the script for generating the results outputs of simulations
# since rda file is around 10GB, each main scenario would have 



#=========================# note 05/30/2023 ====================================
# old codes, once doulbe checked the new codes, remove this script to bin. 
#=========================# note 05/30/2023 ====================================

rm(list = ls()) 

library(here)
here()

# Load useful libraries
library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("gridExtra")
library(doParallel)
library("ggthemes")
registerDoParallel(cores = detectCores() - 1)
# Setup directories after setting working directory to source file 
# directory 

#file path of "TWHCV-model" project
codepath <- file.path(here() %>% dirname(), 'TWHCV-model/03. Code/Functions')

DataFolder <- file.path(here(), "01. DATA/model input")

# output file path 
# dir.create("02. Output") # create subdircetory 
# dir.create("02. Output/RDA") 
# dir.create("02. Output/Results")
# dir.create("02. Output/Figs")
ResultsFolder <- file.path(here(), "02. Output/Results")

outputdt <- here("02. Output/RDA")

outputfig <- here("02. Output/Figs")

# source 
source(file.path(codepath, "HCV_model.R"))

source(file.path(codepath, "plotFunctions.R"))

source(file.path(codepath, "plotOptions.R"))

source(file.path(codepath, "plotManuscript.R")) 

#load rda files 

files <- list.files(path = paste0(outputdt, "/", sep = ""),
                    pattern = '^simBase.*\\.rda')
for( f in files){ 
  
  load(file.path(paste0(outputdt, "/" ,f, sep ="")))
  
}

# Rda file path 
# load the .rda file of base estimate 
rdapath <- file.path(here()%>%dirname(), "Taiwan-MSM-HCV-model")

projectFile <- file.path(rdapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile)

# simulation timeframe 
simY <- 100

# scenario names 
PrEPS <- c("PrEP", "HIVD", "PrEPnHIVD")

#### combine base scenarios in a list 
bestResults_PrEPHIV[["Base"]] <- bestResults

paramResults_PrEPHIV <- list(paramResults_PrEP, paramResults_HIVD, 
                             paramResults_PrEPnHIVD, paramResults)

names(paramResults_PrEPHIV) <- names(bestResults_PrEPHIV)

# setting names for scenarios
SceLab <- list("Increasing PrEP coverage",
               "Increasing HIV diagnosis rate",
               "Increasing PrEP and HIV diagnosis",
               "PrEP coverage and HIV diangosis rate no change")

names(SceLab) <- names(bestResults_PrEPHIV)

#### Section 1. population calibration: HIV PrEP% and HIVD% ####
options(warn=-1)

popS <- popResults_MidYear(HCV, bestResults_PrEPHIV$Base ,Population = NULL,
                           Disease_prog = NULL , 
                           Cascade = NULL, param = paramResults_PrEPHIV$Base, 
                           endYear = 100)%>%
  as.data.frame()

HIVdiag <- read.csv(file.path(paste0(rdapath, "/01. DATA", sep="/"), 
                              "HIVdiag.csv"), header = TRUE)%>%
  as.data.frame()%>%
  mutate(time = year - 2003, realPop = HIV.diagnosis.rate*100,
         low = lower*100,
         up = upper*100)

# subtract the number of diagnosed in the model
Allpop <- list()

subpop_percent <- list()

subpop_percent_range <- list()

HIVInfectedD <- list()

HIVInfectedUND <- list()

HIVDiagRate <- list()

HIVDiagRate_range <- list()

HIVInfected <- list() 

HIVInf <- list()

for(i in names(bestResults_PrEPHIV)){
  
  Allpop[[i]] <- popResults_MidYear(HCV, bestResults_PrEPHIV[[i]],
                                    Population = NULL,
                                    Disease_prog = NULL , 
                                    Cascade = NULL, 
                                    param = paramResults_PrEPHIV[[i]], 
                                    endYear = 100)%>%
    as.data.frame()
  
  HIVInfected[[i]] <- popResults_MidYear(HCV, bestResults_PrEPHIV[[i]],
                     Population = c("HIV+", "HIV+d"),
                     Disease_prog = NULL, 
                     Cascade = NULL, param = paramResults_PrEPHIV[[i]], 
                     endYear = 100) 
  
  HIVInf[[i]] <- HIVInfected[[i]]%>%ungroup()%>% dplyr::group_by(year)%>%
    summarise_at(.vars = c(colnames(popS)[-1]),sum)
  
  HIVInfectedD[[i]] <-  HIVInfected[[i]]%>%filter(population =="HIV+d")
    
  HIVInfectedUND[[i]] <-  HIVInfected[[i]]%>%filter(population =="HIV+")
  
  HIVDiagRate[[i]] <- cbind(year = seq(HCV$startYear , 100-1 ,1),
                       as.data.frame(HIVInfectedD[[i]][, -c(1,2)]/ HIVInf[[i]][ ,-1])*100)%>%
    tibble::as_tibble() 
  
  
  HIVDiagRate_range[[i]] <- popResults_range(HCV, HIVDiagRate[[i]], Population = NULL,
                                        Disease_prog = NULL , 
                                        Cascade = NULL, end_Y = 100) 
  
  subpop_percent[["HIV+d"]][[i]] <- cbind(year = seq(HCV$startYear , 100-1 ,1),
                                         as.data.frame(HIVInfectedD[[i]][, -c(1,2)]/ Allpop[[i]][ ,-1])*100)%>%
    tibble::as_tibble()
  
  subpop_percent_range[["HIV+d"]][[i]] <- popResults_range(HCV, 
                                                          subpop_percent[["HIV+d"]][[i]], 
                                                          Population = NULL,
                                                          Disease_prog = NULL , 
                                                          Cascade = NULL, end_Y = 100) 
  
  subpop_percent[["HIV+"]][[i]] <- cbind(year = seq(HCV$startYear , 100-1 ,1),
                                         as.data.frame(HIVInfectedUND[[i]][, -c(1,2)]/ Allpop[[i]][ ,-1])*100)%>%
    tibble::as_tibble()
  
  subpop_percent_range[["HIV+"]][[i]] <- popResults_range(HCV, 
                                                          subpop_percent[["HIV+"]][[i]], 
                                                          Population = NULL,
                                                          Disease_prog = NULL , 
                                                          Cascade = NULL, end_Y = 100) 
  
  
}


# get the quantile of HIV diagnosed number  
HIVDiagPlot <- list()
for(i in names(bestResults_PrEPHIV)){
  HIVDiagPlot[[i]] <- indicatorPlot(HIVDiagRate_range[[i]] , 
                             ylabel = "HIV Diagnosis among HIV positive MSM (%)",
                             xlimits = c(HCV$startYear, 
                                         HCV$startYear+26, 1),
                             calibration_Y = 2004,
                             rangeun = "y", 
                             groupPlot = NULL, 
                             facetPlot = NULL,
                             observationData = HIVdiag, 
                             simulateYear = HCV$simY) + 
  ggtitle(SceLab[[i]]) + 
    geom_vline(xintercept  = (HCV$simY - HCV$cabY +1), linetype = "dotted", size = 1.2)
  
  HIVDiagPlot[[i]] <- HIVDiagPlot[[i]] + scale_y_continuous(limits = c(0, 100)) + 
    theme_Publication(base_size = 22) + 
    scale_x_continuous(limits = c(1, 27), breaks = c(2,7,12,17,19,22,27),
                       labels = c(2005, 2010, 2015, 2020, 2022, 2025, 2030))
    
    

}

# HIV PrEP %
HIVN <- list()

HIVPrEPRate <- list()

HIVPrEPRate_range <- list()

HIVPrEP <- list() 

HIVnotPrEP <- list()

HIVNeg <- list()

for(i in names(bestResults_PrEPHIV)){
  HIVN[[i]] <- popResults_MidYear(HCV, bestResults_PrEPHIV[[i]],
                                         Population = c("HIV-", "HIV-PrEP"),
                                         Disease_prog = NULL, 
                                         Cascade = NULL, param = paramResults_PrEPHIV[[i]], 
                                         endYear = 100) 
  
  HIVNeg[[i]] <- HIVN[[i]]%>%ungroup()%>% dplyr::group_by(year)%>%
    summarise_at(.vars = c(colnames(popS)[-1]),sum)
  
  HIVPrEP[[i]] <- HIVN[[i]]%>%filter(population =="HIV-PrEP")
  
  HIVnotPrEP[[i]] <- HIVN[[i]]%>%filter(population =="HIV-")

  HIVPrEPRate[[i]] <- cbind(year = seq(HCV$startYear , 100-1 ,1),
                            as.data.frame(HIVPrEP[[i]][, -c(1,2)]/HIVNeg[[i]][ ,-1])*100)%>%
    tibble::as_tibble() 
  
  
  HIVPrEPRate_range[[i]] <- popResults_range(HCV, HIVPrEPRate[[i]], Population = NULL,
                                             Disease_prog = NULL , 
                                             Cascade = NULL, end_Y = 30) 
  
  subpop_percent[["HIV-"]][[i]] <- cbind(year = seq(HCV$startYear , 100-1 ,1),
                                          as.data.frame(HIVnotPrEP[[i]][, -c(1,2)]/ Allpop[[i]][ ,-1])*100)%>%
    tibble::as_tibble()
  
  subpop_percent_range[["HIV-"]][[i]] <- popResults_range(HCV, 
                                                           subpop_percent[["HIV-"]][[i]], 
                                                           Population = NULL,
                                                           Disease_prog = NULL , 
                                                           Cascade = NULL, end_Y = 100) 
  
  subpop_percent[["HIV-PrEP"]][[i]] <- cbind(year = seq(HCV$startYear , 100-1 ,1),
                                         as.data.frame(HIVPrEP[[i]][, -c(1,2)]/ Allpop[[i]][ ,-1])*100)%>%
    tibble::as_tibble()
  
  subpop_percent_range[["HIV-PrEP"]][[i]] <- popResults_range(HCV, 
                                                          subpop_percent[["HIV-PrEP"]][[i]], 
                                                          Population = NULL,
                                                          Disease_prog = NULL , 
                                                          Cascade = NULL, end_Y = 100)
  
  
}

save(HIVDiagRate_range, HIVPrEPRate_range, subpop_percent_range,
     file = file.path(outputdt, "pop_rate.rda"))
# get the quantile of HIV diagnosed number  
HIVPrEPPlot <- list()

for(i in names(bestResults_PrEPHIV)){
  HIVPrEPPlot[[i]] <- indicatorPlot(HIVPrEPRate_range[[i]] , 
                                    ylabel = "Number of PrEP users among HIV negative MSM (%)",
                                    xlimits = c(HCV$startYear, 
                                                HCV$startYear+26, 1),
                                    calibration_Y = 2004,
                                    rangeun = "y", 
                                    groupPlot = NULL, 
                                    facetPlot = NULL,
                                    observationData = NULL, 
                                    simulateYear = HCV$simY) + 
    ggtitle(SceLab[[i]]) + 
    geom_vline(xintercept  = (HCV$simY - HCV$cabY +1), linetype = "dotted", size = 1.2)
  
  HIVPrEPPlot[[i]] <- HIVPrEPPlot[[i]] + scale_y_continuous(limits = c(0, 25)) + 
    theme_Publication(base_size = 22) + 
    scale_x_continuous(limits = c(1, 27), breaks = c(2,7,12,17,19,22,27),
                       labels = c(2005, 2010, 2015, 2020, 2022, 2025, 2030))
  
}

simpop <- list()
for (i in names(bestResults_PrEPHIV)){ 
  simpop[[i]] <- ggpubr::ggarrange(HIVDiagPlot[[i]], HIVPrEPPlot[[i]], ncol = 2, 
                                   common.legend = TRUE,  legend="bottom") 
  
  ggsave(simpop[[i]], filename = here(outputfig, paste0( "simpop_PrEPHIV_", i,".png")),
         width = 15 , height = 9)
  
}
#### Epi outcomes ####
# including HCV new infections and HCV incidences 

# base 
bR <- list()

bPar <- list()

btt <- list()

bttnumpop <- list()

thres <- list()

bttpop <- list()

threspop <- list()

simY <- 100

for (i in names(bestResults_PrEPHIV)){ 
  
  bR[["Base"]][[i]] <- list(bestResults_PrEPHIV[[i]])
  
  bPar[["Base"]][[i]] <- list(paramResults_PrEPHIV[[i]])
  
  btt[[i]] <- scenario_Incidence_allset(bR[["Base"]][[i]],bPar[["Base"]][[i]],
                                 pop = NULL, statusQ = "y", endY = simY)
  
  thres[[i]] <- btt[[i]]%>%mutate(year = 2003 + year)%>%filter(year == 2015)%>%
    select(best)%>%mutate(thres = best*0.2)
  
  bttpop[[i]] <- scenario_Incidence_allset(bR[["Base"]][[i]],bPar[["Base"]][[i]],
                                pop = "pop", statusQ = "y", endY = simY)
  
  threspop[[i]] <- bttpop[[i]]%>%mutate(year = 2003 + year)%>%
    filter(year == 2015)%>%select(best)%>%mutate(thres = best*0.2)
  
  threspop[[i]][2,] <- threspop[[i]][1,]
  
  bttnumpop[[i]] <- scenario_Incidence_allset(bR[["Base"]][[i]],bPar[["Base"]][[i]],
                                      pop = "pop", statusQ = "y", 
                                      indicator = "newInfections", endY = simY)
  
}

save(btt, bttpop, thres, threspop, bttnumpop,
     file = file.path(outputdt, "Outcome_Base_epi.rda"))

#### cost outcomes #### 
# base 
indicators <- c("costTestingAb", "costTestingAg", "costnewTestingPOCT", 
                "costTreatment", "costCured", "costRetreat") 

indicators_allp <- c("QALYPops", "costPops")

Base_cost_flow <- list()

Base_cost_box <- list()

Base_cost_box_Range <-list() 

for(i in indicators){
  for(m in names(bestResults_PrEPHIV)){
    
    Base_cost_flow[[i]][[m]] <- indicatorResults(HCV, bestResults_PrEPHIV[[m]],
                                                      i, 
                                                      paramR = paramResults_PrEPHIV[[m]],
                                                      pop = NULL, range = "y", 
                                                      simY, scenario = NULL)
  }
}

save(Base_cost_flow,
     file = file.path(outputdt, "Outcome_Base_cost_flow.rda"))


for(i in indicators_allp){
  for(m in names(bestResults_PrEPHIV)){
    
    Base_cost_box[[i]][[m]] <- 
      popResults_MidYear(HCV, bestResults_PrEPHIV[[m]], 
                         Population = HCV$popNames,
                         Disease_prog = HCV$progress_name, 
                         Cascade = c("s",HCV$cascade_name), 
                         param = paramResults_PrEPHIV[[m]],
                         endYear = simY, 
                         allp = i)%>%as.data.frame() 
    
    
    Base_cost_box_Range[[i]][[m]] <- 
      popResults_range(HCV, Base_cost_box[[i]][[m]], 
                       Population = NULL,
                       Disease_prog = NULL, 
                       Cascade = NULL, 
                       end_Y = simY)
  }
}

save(Base_cost_box, Base_cost_box_Range,
     file = file.path(outputdt, "Outcome_Base_cost_box.rda"))

rm(bestResults_PrEPHIV, paramResults_PrEPHIV, btt, bttpop, thres, threspop, 
   bttnumpop,Base_cost_flow, Base_cost_box, Base_cost_box_Range)

#### outcomes for scenarios #### 
load(file.path(paste0(outputdt, "/simScen.rda")))

tic <- proc.time()
# base scenarios 
Sttpop_base <- scenario_Incidence_allset(ScenResults, 
                                         ScenparamResults, 
                                         pop = "pop", endY = simY) 

Stt_base <- scenario_Incidence_allset(ScenResults, 
                                         ScenparamResults, 
                                         pop = NULL, endY = simY) 

Sttnumpop_base <- scenario_Incidence_allset(ScenResults, 
                                            ScenparamResults,
                                            pop = "pop", 
                                            indicator = "newInfections", 
                                            endY = simY)

save(Sttpop_base,Stt_base, Sttnumpop_base, 
     file = file.path(outputdt, "Outcome_Scen_epi_base.rda"))

rm(Sttpop_base,Stt_base, Sttnumpop_base)
gc()
#### cost: flow 
Scen_cost_flow_base <- list()

Scen_cost_box_base <- list()

Scen_cost_box_base_Range <-list()



for(i in indicators){
  for(m in names(ScenResults)){
    
    Scen_cost_flow_base[[i]][[m]] <- indicatorResults(HCV, ScenResults[[m]],
                                                      i, 
                                                      paramR = ScenparamResults[[m]],
                                                      pop = NULL, range = "y", 
                                                      simY, scenario = NULL)
  }
}

save(Scen_cost_flow_base,
     file = file.path(outputdt, "Outcome_Scen_cost_flow_base.rda"))

rm(Scen_cost_flow_base)



for(i in indicators_allp){
  for(m in names(ScenResults)){
    
    Scen_cost_box_base[[i]][[m]] <- 
      popResults_MidYear(HCV, ScenResults[[m]], 
                         Population = HCV$popNames,
                         Disease_prog = HCV$progress_name, 
                         Cascade = c("s",HCV$cascade_name), 
                         param = ScenparamResults[[m]],
                         endYear = simY, 
                         allp = i)%>%as.data.frame() 
    
    
    Scen_cost_box_base_Range[[i]][[m]] <- 
      popResults_range(HCV, Scen_cost_box_base[[i]][[m]], 
                       Population = NULL,
                       Disease_prog = NULL, 
                       Cascade = NULL, 
                       end_Y = simY)
  }
}

save(Scen_cost_box_base, Scen_cost_box_base_Range,
     file = file.path(outputdt, "Outcome_Scen_cost_box_base.rda"))

rm(Scen_cost_box_base, Scen_cost_box_base_Range, ScenResults, ScenparamResults)

toc <- proc.time() - tic
toc
gc()
# epi outcomes for scenarios
load(file.path(paste0(outputdt, "/simScen_PrEPHIV.rda")))

# PrEP
tic <- proc.time()

Outcome_Scen_epi_PrEP <- list()

a <- get(load(file.path(paste0(outputdt, "/simScen_param_PrEP.rda"))))

Outcome_Scen_epi_PrEP[["Sttpop"]] <- 
  scenario_Incidence_allset(ScenResults_PrEPHIV[["PrEP"]], a, 
                            pop = "pop", endY = simY)

Outcome_Scen_epi_PrEP[["Stt"]] <- 
  scenario_Incidence_allset(ScenResults_PrEPHIV[["PrEP"]], a, 
                            pop = NULL, endY = simY)

Outcome_Scen_epi_PrEP[["Sttnumpop"]] <- 
  scenario_Incidence_allset(ScenResults_PrEPHIV[["PrEP"]],
                            a, pop = "pop", indicator = "newInfections", 
                            endY = simY)

Scen_cost_flow_PrEP <- list()

Scen_cost_box_PrEP <- list()

Scen_cost_box_range_PrEP <- list()

for(f in indicators){
  for(m in names(ScenResults_PrEPHIV[[1]])){
    
    Scen_cost_flow_PrEP[[f]][[m]] <- 
      indicatorResults(HCV, ScenResults_PrEPHIV[["PrEP"]][[m]], 
                       f, 
                       paramR = a[[m]],
                       pop = NULL, range = "y", simY, scenario = NULL)
    
    
    
  }
}


for(f in indicators_allp){
  for(m in names(ScenResults_PrEPHIV[[1]])){
    Scen_cost_box_PrEP[[f]][[m]] <- 
      popResults_MidYear(HCV, ScenResults_PrEPHIV[["PrEP"]][[m]], 
                         Population = HCV$popNames,
                         Disease_prog = HCV$progress_name, 
                         Cascade = c("s",HCV$cascade_name), 
                         param = a[[m]],
                         endYear = simY, 
                         allp = f)%>%as.data.frame() 
    
    Scen_cost_box_range_PrEP[[f]][[m]] <- 
      popResults_range(HCV, Scen_cost_box_PrEP[[f]][[m]], 
                       Population = NULL,
                       Disease_prog = NULL, 
                       Cascade = NULL, 
                       end_Y = simY)
  }
}

save(Outcome_Scen_epi_PrEP, 
     file = file.path(outputdt, paste0("Outcome_Scen_epi_PrEP.rda")))

save(Scen_cost_flow_PrEP, 
     file = file.path(outputdt, paste0("Outcome_Scen_cost_flow_PrEP.rda")))

save(Scen_cost_box_PrEP, Scen_cost_box_range_PrEP,
     file = file.path(outputdt, paste0("Outcome_Scen_cost_box_PrEP.rda")))

rm(Outcome_Scen_epi_PrEP, Scen_cost_flow_PrEP, 
   Scen_cost_box_PrEP, Scen_cost_box_range_PrEP, 
   a, ScenparamResults_PrEP)

gc()
# HIVD
tic <- proc.time()

Outcome_Scen_epi_HIVD <- list()

a <- get(load(file.path(paste0(outputdt, "/simScen_param_HIVD.rda"))))

Outcome_Scen_epi_HIVD[["Sttpop"]] <- 
  scenario_Incidence_allset(ScenResults_PrEPHIV[["HIVD"]], a, 
                            pop = "pop", endY = simY)

Outcome_Scen_epi_HIVD[["Stt"]] <- 
  scenario_Incidence_allset(ScenResults_PrEPHIV[["HIVD"]], a, 
                            pop = NULL, endY = simY)

Outcome_Scen_epi_HIVD[["Sttnumpop"]] <- 
  scenario_Incidence_allset(ScenResults_PrEPHIV[["HIVD"]],
                            a, pop = "pop", indicator = "newInfections", 
                            endY = simY)

save(Outcome_Scen_epi_HIVD, 
     file = file.path(outputdt, paste0("Outcome_Scen_epi_HIVD.rda")))

rm(Outcome_Scen_epi_HIVD)

Scen_cost_flow_HIVD <- list()

Scen_cost_box_HIVD <- list()

Scen_cost_box_range_HIVD <- list()

for(f in indicators){
  for(m in names(ScenResults_PrEPHIV[[2]])){
    
    Scen_cost_flow_HIVD[[f]][[m]] <- 
      indicatorResults(HCV, ScenResults_PrEPHIV[["HIVD"]][[m]], 
                       f, 
                       paramR = a[[m]],
                       pop = NULL, range = "y", simY, scenario = NULL)
    
    
    
  }
}

save(Scen_cost_flow_HIVD, 
     file = file.path(outputdt, paste0("Outcome_Scen_cost_flow_HIVD.rda")))

rm(Scen_cost_flow_HIVD)
gc()
for(f in indicators_allp){
  for(m in names(ScenResults_PrEPHIV[[2]])){
    Scen_cost_box_HIVD[[f]][[m]] <- 
      popResults_MidYear(HCV, ScenResults_PrEPHIV[["HIVD"]][[m]], 
                         Population = HCV$popNames,
                         Disease_prog = HCV$progress_name, 
                         Cascade = c("s",HCV$cascade_name), 
                         param = a[[m]],
                         endYear = simY, 
                         allp = f)%>%as.data.frame() 
    
    Scen_cost_box_range_HIVD[[f]][[m]] <- 
      popResults_range(HCV, Scen_cost_box_HIVD[[f]][[m]], 
                       Population = NULL,
                       Disease_prog = NULL, 
                       Cascade = NULL, 
                       end_Y = simY)
  }
}


save(Scen_cost_box_HIVD, Scen_cost_box_range_HIVD,
     file = file.path(outputdt, paste0("Outcome_Scen_cost_box_HIVD.rda")))

rm(Scen_cost_box_HIVD, Scen_cost_box_range_HIVD, 
   a, ScenparamResults_HIVD)


gc()

# PrEPnHIVD
tic <- proc.time()

Outcome_Scen_epi_PrEPnHIVD <- list()

a <- get(load(file.path(paste0(outputdt, "/simScen_param_PrEPnHIVD.rda"))))

Outcome_Scen_epi_PrEPnHIVD[["Sttpop"]] <- 
  scenario_Incidence_allset(ScenResults_PrEPHIV[["PrEPnHIVD"]], a, 
                            pop = "pop", endY = simY)

Outcome_Scen_epi_PrEPnHIVD[["Stt"]] <- 
  scenario_Incidence_allset(ScenResults_PrEPHIV[["PrEPnHIVD"]], a, 
                            pop = NULL, endY = simY)

Outcome_Scen_epi_PrEPnHIVD[["Sttnumpop"]] <- 
  scenario_Incidence_allset(ScenResults_PrEPHIV[["PrEPnHIVD"]],
                            a, pop = "pop", indicator = "newInfections", 
                            endY = simY)

save(Outcome_Scen_epi_PrEPnHIVD, 
     file = file.path(outputdt, paste0("Outcome_Scen_epi_PrEPnHIVD.rda")))

rm(Outcome_Scen_epi_PrEPnHIVD)

Scen_cost_flow_PrEPnHIVD <- list()

Scen_cost_box_PrEPnHIVD <- list()

Scen_cost_box_range_PrEPnHIVD <- list()

for(f in indicators){
  for(m in names(ScenResults_PrEPHIV[[3]])){
    
    Scen_cost_flow_PrEPnHIVD[[f]][[m]] <- 
      indicatorResults(HCV, ScenResults_PrEPHIV[["PrEPnHIVD"]][[m]], 
                       f, 
                       paramR = a[[m]],
                       pop = NULL, range = "y", simY, scenario = NULL)
    
    
    
  }
}

save(Scen_cost_flow_PrEPnHIVD, 
     file = file.path(outputdt, paste0("Outcome_Scen_cost_flow_PrEPnHIVD.rda")))

rm(Scen_cost_flow_PrEPnHIVD)

for(f in indicators_allp){
  for(m in names(ScenResults_PrEPHIV[[3]])){
    Scen_cost_box_PrEPnHIVD[[f]][[m]] <- 
      popResults_MidYear(HCV, ScenResults_PrEPHIV[["PrEPnHIVD"]][[m]], 
                         Population = HCV$popNames,
                         Disease_prog = HCV$progress_name, 
                         Cascade = c("s",HCV$cascade_name), 
                         param = a[[m]],
                         endYear = simY, 
                         allp = f)%>%as.data.frame() 
    
    Scen_cost_box_range_PrEPnHIVD[[f]][[m]] <- 
      popResults_range(HCV, Scen_cost_box_PrEPnHIVD[[f]][[m]], 
                       Population = NULL,
                       Disease_prog = NULL, 
                       Cascade = NULL, 
                       end_Y = simY)
  }
}

toc <- proc.time() - tic

toc

save(Scen_cost_box_PrEPnHIVD, Scen_cost_box_range_PrEPnHIVD,
     file = file.path(outputdt, paste0("Outcome_Scen_cost_box_PrEPnHIVD.rda")))

rm(Scen_cost_box_PrEPnHIVD, Scen_cost_box_range_PrEPnHIVD, 
   a, ScenparamResults_PrEPnHIVD)
################################################################################




