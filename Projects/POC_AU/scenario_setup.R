# coverage estimation: number of people tested via national program/popultion size 
# This script replaced PAF_coverage_cal.R & scenarioEstimate_11202023.R (old estimation method) since 2024/02/19
# We used the data (number of people tested in 2022, 2023 in communities and prisons) provided by Simon to estimate the coverage
# further we used the coverage to estimate the parameters on tau_ab, tau_RNA, tau_poct, eta for the nationa program 
# set up the scenarios 

rm(list = ls())
library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
library(readxl)
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
source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotManuscript.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R"))

# import the number of test 

Num_test_person_NP <- read_excel(paste0(data_path, "/01. DATA/Num_test_person_NP.xlsx"))

Num_test_person_NP <- Num_test_person_NP%>%
  mutate(num_pop = c(480000, 480000, 80000, 80000))

Num_test_person_NP <- Num_test_person_NP%>%
  mutate(coverage = num_person/num_pop, 
         reflex_frac = num_ab/num_person, 
         immedRNA_frac = 1 - reflex_frac) 

# current achievement data 
current_data <- function(dt, y, s, index){ 
  x <- dt%>%filter(year == y & settings == s)
  x <- x[[index]]
  return(x)
  }

# extracted the coverage of national program and fraction of testing pathways 

frac_test <- list() 
reflex_frac_C <- list()
reflex_frac_P <- list()
Ccal <- list()

for(i in c(2022, 2023)){ 
  
  reflex_frac_C[[i]] <- current_data(Num_test_person_NP, y = i, 
                                s = "community", index = "reflex_frac")
  
  reflex_frac_P[[i]] <- current_data(Num_test_person_NP, y = i, 
                                s = "prison", index = "reflex_frac")
  
  frac_test[[i]] <- list("C" = list("reflex" = as.numeric(reflex_frac_C[[i]]), 
                                    "immeRNA" = 1- as.numeric(reflex_frac_C[[i]])),
                         "P" = list("reflex" = as.numeric(reflex_frac_P[[i]]), 
                                    "immeRNA" = 1- as.numeric(reflex_frac_P[[i]])))
  
  Ccal[[i]] <- list("C" =  current_data(Num_test_person_NP, y = i, 
                                           s = "community", index = "coverage"),
                       "P" = current_data(Num_test_person_NP, y = i, 
                                          s = "prison", index = "coverage")) 
  }

frac_test[[2022]][["C"]][["reflex"]]



NP_tauRNA_C <- 0.939

NP_tauRNAonly_C <- 0.709

NP_eta_C <- 0.62

NP_tauRNA_P <- 0.939



NP_tauRNAonly_P <- 0.709

NP_eta_P <- 0.89

NPlst <- list("C" = list("tau_RNA" = NP_tauRNA_C, 
                         "tau_poct" = NP_tauRNAonly_C,
                         "eta" = NP_eta_C),
              "P" = list("tau_RNA" = NP_tauRNA_P, 
                         "tau_poct" = NP_tauRNAonly_P,
                         "eta" = NP_eta_P))

Param_cal <- function(pj, dlist, index ,S_Yint, S_Yend, r_Yend, NPlst, Ccal, frac_testing = NULL){ 
  #proj: project name 
  # dlist: list of cascade parameters 
  # index: the parameter to calculate 
  # SYint: starting year of scenario 
  # S_Yend: end year of scenario 
  # NPlst: the odd ratio of the parameters in scenario (save in a list)
  # Ccal: the list of national program coverage 
  # frac_testing: te list of fraction of testing pathways in national program 
  # convert to monthly probability then convert back to yearly probability 
  SYpoint_int <- (S_Yint - pj$cabY)/pj$timestep + 1
  SYpoint_end <- (S_Yend - pj$cabY)/pj$timestep
  
  SY_leng <- SYpoint_end - SYpoint_int + 1 
  # remain the same level to end of which year
  rYpoint_end <- (r_Yend - pj$cabY)/POC_AU$timestep
  
  rY_leng <- rYpoint_end - SYpoint_end
  
  intVal <- dlist[[index]][, 3, (SYpoint_int - 1)]
  
  endVal <- c()  
  if(is.null(frac_testing) & isTRUE(index%in% c("tau_ab"))){ 
    endVal[1] <- 1-(1- intVal[1])^pj$timestep + (1-(1- intVal[1]*Ccal[["C"]])^pj$timestep)
    endVal[2] <- 1-(1- intVal[2])^pj$timestep + (1-(1- intVal[2]*Ccal[["C"]])^pj$timestep)
    endVal[3] <- 1-(1- intVal[3])^pj$timestep + (1-(1- intVal[3]*Ccal[["P"]])^pj$timestep)
    endVal[4] <- 1-(1- intVal[4])^pj$timestep + (1-(1- intVal[4]*Ccal[["P"]])^pj$timestep)
    endVal[5] <- 1-(1- intVal[5])^pj$timestep + (1-(1- intVal[5]*Ccal[["P"]])^pj$timestep)
    
  }
  else if(is.null(frac_testing)){ 
    endVal[1] <- 1-(1- intVal[1])^pj$timestep + (1-(1- (NPlst[["C"]][[index]])*Ccal[["C"]])^pj$timestep)
    endVal[2] <- 1-(1- intVal[2])^pj$timestep + (1-(1- (NPlst[["C"]][[index]])*Ccal[["C"]])^pj$timestep)
    endVal[3] <- 1-(1- intVal[3])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]])^pj$timestep)
    endVal[4] <- 1-(1- intVal[4])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]])^pj$timestep)
    endVal[5] <- 1-(1- intVal[5])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]])^pj$timestep)
    
  }
  else if(!is.null(frac_testing) & isTRUE(index%in% c("tau_ab"))){ 
    endVal[1] <- 1-(1- intVal[1])^pj$timestep + (1-(1- intVal[1]*Ccal[["C"]]*frac_testing[["C"]][["reflex"]])^pj$timestep)
    endVal[2] <- 1-(1- intVal[2])^pj$timestep + (1-(1- intVal[2]*Ccal[["C"]]*frac_testing[["C"]][["reflex"]])^pj$timestep)
    endVal[3] <- 1-(1- intVal[3])^pj$timestep + (1-(1- intVal[3]*Ccal[["P"]]*frac_testing[["P"]][["reflex"]])^pj$timestep)
    endVal[4] <- 1-(1- intVal[4])^pj$timestep + (1-(1- intVal[4]*Ccal[["P"]]*frac_testing[["P"]][["reflex"]])^pj$timestep)
    endVal[5] <- 1-(1- intVal[5])^pj$timestep + (1-(1- intVal[5]*Ccal[["P"]]*frac_testing[["P"]][["reflex"]])^pj$timestep)
    
  }
  else if(!is.null(frac_testing) & isTRUE(index%in% c("tau_RNA"))){ 
    endVal[1] <- 1-(1- intVal[1])^pj$timestep + (1-(1- (NPlst[["C"]][[index]])*Ccal[["C"]]*frac_testing[["C"]][["reflex"]])^pj$timestep)
    endVal[2] <- 1-(1- intVal[2])^pj$timestep + (1-(1- (NPlst[["C"]][[index]])*Ccal[["C"]]*frac_testing[["C"]][["reflex"]])^pj$timestep)
    endVal[3] <- 1-(1- intVal[3])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]]*frac_testing[["P"]][["reflex"]])^pj$timestep)
    endVal[4] <- 1-(1- intVal[4])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]]*frac_testing[["P"]][["reflex"]])^pj$timestep)
    endVal[5] <- 1-(1- intVal[5])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]]*frac_testing[["P"]][["reflex"]])^pj$timestep)
    
  }
  else if(!is.null(frac_testing) & isTRUE(index == "tau_poct")){ 
    endVal[1] <- 1-(1- intVal[1])^pj$timestep + (1-(1- (NPlst[["C"]][[index]])*Ccal[["C"]]*frac_testing[["C"]][["immeRNA"]])^pj$timestep)
    endVal[2] <- 1-(1- intVal[2])^pj$timestep + (1-(1- (NPlst[["C"]][[index]])*Ccal[["C"]]*frac_testing[["C"]][["immeRNA"]])^pj$timestep)
    endVal[3] <- 1-(1- intVal[3])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]]*frac_testing[["P"]][["immeRNA"]])^pj$timestep)
    endVal[4] <- 1-(1- intVal[4])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]]*frac_testing[["P"]][["immeRNA"]])^pj$timestep)
    endVal[5] <- 1-(1- intVal[5])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]]*frac_testing[["P"]][["immeRNA"]])^pj$timestep)
    
    
    
  }
  else if(!is.null(frac_testing) & isTRUE(index == "eta")){ 
    endVal[1] <- 1-(1- intVal[1])^pj$timestep + (1-(1- (NPlst[["C"]][[index]])*Ccal[["C"]])^pj$timestep)
    endVal[2] <- 1-(1- intVal[2])^pj$timestep + (1-(1- (NPlst[["C"]][[index]])*Ccal[["C"]])^pj$timestep)
    endVal[3] <- 1-(1- intVal[3])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]])^pj$timestep)
    endVal[4] <- 1-(1- intVal[4])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]])^pj$timestep)
    endVal[5] <- 1-(1- intVal[5])^pj$timestep + (1-(1- (NPlst[["P"]][[index]])*Ccal[["P"]])^pj$timestep)
    
    
    
  }
  
  
  endVal[1] <- 1- (1 - endVal[1])^(1/pj$timestep)
  endVal[2] <- 1- (1 - endVal[2])^(1/pj$timestep)
  endVal[3] <- 1- (1 - endVal[3])^(1/pj$timestep)
  endVal[4] <- 1- (1 - endVal[4])^(1/pj$timestep)
  endVal[5] <- 1- (1 - endVal[5])^(1/pj$timestep)
  
  
  if(index== "tau_poct"){ 
    avcov <- dlist[["tau_ab"]][, 3, SYpoint_int]
    endVal[1] <- 
      ifelse((endVal[1] + avcov[1])>= 1, (1 - avcov[1]), endVal[1])
    endVal[2] <- 
      ifelse((endVal[2] + avcov[2])>= 1, (1 - avcov[2]), endVal[2])
    endVal[3] <- 
      ifelse((endVal[3] + avcov[3])>= 1, (1 - avcov[3]), endVal[3])
    endVal[4] <- 
      ifelse((endVal[4] + avcov[4])>= 1, (1 - avcov[4]), endVal[4])
    endVal[5] <- 
      ifelse((endVal[5] + avcov[5])>= 1, (1 - avcov[5]), endVal[5])
  }
  
  
  if(isTRUE(rYpoint_end > SYpoint_end)){ 
    for ( i in 2:dim(dlist[[index]])[[2]]){
      dlist[[index]][1, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[1]), as.numeric(endVal[1]), length = (SY_leng + 1)), 
          rep(as.numeric(endVal[1]), length = (rY_leng)),
          rep(intVal[1], pj$npts - rYpoint_end))
      
      dlist[[index]][2, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[2]), as.numeric(endVal[2]), length = (SY_leng + 1)), 
          rep(as.numeric(endVal[2]), length = (rY_leng)),
          rep(intVal[2], pj$npts - rYpoint_end))
      
      dlist[[index]][3, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[3]), as.numeric(endVal[3]), length = (SY_leng + 1)), 
          rep(as.numeric(endVal[3]), length = (rY_leng)),
          rep(intVal[3], pj$npts - rYpoint_end))
      
      dlist[[index]][4, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[4]), as.numeric(endVal[4]), length = (SY_leng + 1)), 
          rep(as.numeric(endVal[4]), length = (rY_leng)),
          rep(intVal[4], pj$npts - rYpoint_end))
      
      dlist[[index]][5, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[5]), as.numeric(endVal[5]), length = (SY_leng + 1)), 
          rep(as.numeric(endVal[5]), length = (rY_leng)),
          rep(intVal[5], pj$npts - rYpoint_end)) 
    }
    
  } else{ 
    for ( i in 2:dim(dlist[[index]])[[2]]){
      dlist[[index]][1, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[1]), as.numeric(endVal[1]), length = (SY_leng + 1)), 
          rep(intVal[1],  pj$npts - SYpoint_end))
      
      dlist[[index]][2, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[2]), as.numeric(endVal[2]), length = (SY_leng + 1)), 
          rep(intVal[2], pj$npts - SYpoint_end))
      
      dlist[[index]][3, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[3]), as.numeric(endVal[3]), length = (SY_leng + 1)), 
          rep(intVal[3], pj$npts - SYpoint_end))
      
      dlist[[index]][4, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[4]), as.numeric(endVal[4]), length = (SY_leng + 1)), 
          rep(intVal[4], pj$npts - SYpoint_end))
      
      dlist[[index]][5, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[5]), as.numeric(endVal[5]), length = (SY_leng + 1)), 
          rep(intVal[5], pj$npts - SYpoint_end))
    }
    
  }

    
    
    
  
  return(dlist[[index]])
}

param_var <- c("tau_ab","tau_RNA", "tau_poct", "eta")

#### Scenario 1: current achievement: starting of 2022 to end of 2023, back to pre-national program level since starting of 2024 #### 
# year of 2022
dfList_NP <- dfList
for(i in param_var){ 
  dfList_NP[[i]] <- Param_cal(pj = POC_AU, dlist = dfList, index = i, frac_testing = frac_test[[2022]],
                              S_Yint = 2022, S_Yend = 2023, r_Yend = 2023, NPlst = NPlst, 
                              Ccal = Ccal[[2022]])
}
# no change in prison setting in reflex RNA pathway 
dfList_NP$tau_ab[c(3:5), , ] <- dfList$tau_ab[c(3:5) , ,]
dfList_NP$tau_RNA[c(3:5), , ] <- dfList$tau_RNA[c(3:5) , ,]

# year of 2023 
dfList_NP_2023 <- dfList_NP 

for(i in param_var){ 
  dfList_NP_2023[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP, index = i, frac_testing = frac_test[[2023]],
                              S_Yint = 2023, S_Yend = 2024, r_Yend = 2024, NPlst = NPlst, 
                              Ccal = Ccal[[2023]])
}

# back to pre-national program’s level from starting of 2024
for(i in param_var){ 
  # begining of 2024 
  b_pt <- (2024 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2023[[i]])[3]
  dfList_NP_2023[[i]][, , b_pt: dim_length] <- dfList[[i]][, , b_pt: dim_length]
} 

Ccal_2023 <- Ccal[[2023]]

frac_test_2023 <- frac_test[[2023]]

save(dfList_NP_2023, Ccal_2023,frac_test_2023,
     file = file.path(OutputFolder ,
                      paste0(project_name,"NP_2023" ,".rda")))



#### Scenario 2: current achievement + 2024 prediction: back to pre-national program level since starting of 2025  ####  
# odds of coverage: 1.08 (total test_predicited/ number test in 2023, 19480/18035) 
odd_num_test <- list()
odd_num_test[[2024]] <- 19480/18035
Ccal[[2024]] <- lapply(Ccal[[2023]],function(x) x*odd_num_test[[2024]])

dfList_NP_2024 <- dfList_NP_2023

for(i in param_var){ 
  dfList_NP_2024[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2023, index = i, 
                                   frac_testing = frac_test[[2023]], S_Yint = 2024, 
                                   S_Yend = 2025, r_Yend = 2025, NPlst = NPlst, 
                                   Ccal = Ccal[[2024]])
}

# back to pre-national program’s level from starting of 2025
for(i in param_var){ 
  # begining of 2025
  b_pt <- (2025 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2024[[i]])[3]
  dfList_NP_2024[[i]][, , b_pt: dim_length] <- dfList[[i]][, , b_pt: dim_length]
}

Ccal_2024 <- Ccal[[2024]]

frac_test_2024 <- frac_test[[2023]] # same distribution as 2023

save(dfList_NP_2024, Ccal_2024,frac_test_2024,
     file = file.path(OutputFolder ,
                      paste0(project_name,"NP_2024" ,".rda")))


####National program expansion ####
#### scenario 3A: Scenario 2A and remaining at level of 2024 to end of 2027, back to pre-naitonal program level since starting of 2028 ####

dfList_NPexp_A <- dfList_NP_2023
for(i in param_var){ 
  dfList_NPexp_A[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2023, index = i, 
                                   frac_testing = frac_test[[2023]], S_Yint = 2024, 
                                   S_Yend = 2025, r_Yend = 2028, NPlst = NPlst, 
                                   Ccal = Ccal[[2024]])

}
# back to pre-national program’s level from starting of 2028
for(i in param_var){ 
  # begining of 2028
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPexp_A[[i]])[3]
  dfList_NPexp_A[[i]][, , b_pt: dim_length] <- dfList[[i]][, , b_pt: dim_length]
}

save(dfList_NPexp_A, Ccal_2024,
     file = file.path(OutputFolder ,
                      paste0(project_name,"NP_expand_A" ,".rda")))



#### scenario 3B: further increase coverage in 2025 and remaining at level of 2025 to end of 2027, back to pre-naitonal program level since starting of 2028 #### 
# 2025 coverage estimation: number test_2025/number test_2024 = cov_2024* (30000/20000)
odd_num_test[[2025]] <- 30000/20000
Ccal_2025 <- lapply(Ccal_2024,function(x) x*odd_num_test[[2025]])

dfList_NPexp_B <- dfList_NP_2024
for(i in param_var){ 
  dfList_NPexp_B[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2024, index = i, 
                                     frac_testing = frac_test[[2023]],
                                     S_Yint = 2025, S_Yend = 2026, r_Yend = 2028,
                                     NPlst = NPlst, 
                                     Ccal = Ccal_2025)
}

# back to pre-national program’s level from starting of 2028
for(i in param_var){ 
  # begining of 2028
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPexp_B[[i]])[3]
  dfList_NPexp_B[[i]][, , b_pt: dim_length] <- dfList[[i]][, , b_pt: dim_length]
}

save(dfList_NPexp_B, Ccal_2025,
     file = file.path(OutputFolder ,
                      paste0(project_name,"NP_expand_A" ,".rda")))


#### scenario 3C: further increase coverage in 2026 and remaining at level of 2026 to end of 2027, back to pre-naitonal program level since starting of 2028 #### 
# 2026 coverage estimation: number test_2026/number test_2025 = cov_2025* (40000/30000)
odd_num_test[[2026]] <- 40000/30000
Ccal_2026 <- lapply(Ccal_2025,function(x) x*odd_num_test[[2026]])

dfList_NPexp_C <- dfList_NPexp_B
for(i in param_var){ 
  dfList_NPexp_C[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NPexp_B, index = i, 
                                   frac_testing = frac_test[[2023]],
                                   S_Yint = 2026, S_Yend = 2027, r_Yend = 2028,
                                   NPlst = NPlst, 
                                   Ccal = Ccal_2026)
}

# back to pre-national program’s level from starting of 2028
for(i in param_var){ 
  # begining of 2028
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPexp_C[[i]])[3]
  dfList_NPexp_C[[i]][, , b_pt: dim_length] <- dfList[[i]][, , b_pt: dim_length]
}

save(dfList_NPexp_C, Ccal_2026,
     file = file.path(OutputFolder ,
                      paste0(project_name,"NP_expand_C" ,".rda")))


#### scenario 3D: further increase coverage in 2027, back to pre-naitonal program level since starting of 2028 #### 
# 2027 coverage estimation: number test_2027/number test_2026 = cov_2026* (50000/40000)
odd_num_test[[2027]] <- 50000/40000
Ccal_2027 <- lapply(Ccal_2026,function(x) x*odd_num_test[[2027]])

dfList_NPexp_D <- dfList_NPexp_C
for(i in param_var){ 
  dfList_NPexp_D[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NPexp_C, index = i, 
                                   frac_testing = frac_test[[2023]],
                                   S_Yint = 2027, S_Yend = 2028, r_Yend = 2028,
                                   NPlst = NPlst, 
                                   Ccal = Ccal_2027)
}

# back to pre-national program’s level from starting of 2028
for(i in param_var){ 
  # begining of 2028
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPexp_D[[i]])[3]
  dfList_NPexp_D[[i]][, , b_pt: dim_length] <- dfList[[i]][, , b_pt: dim_length]
}

save(dfList_NPexp_D, Ccal_2027,
     file = file.path(OutputFolder ,
                      paste0(project_name,"NP_expand_D" ,".rda")))



# extract the additional probability contributed by national program 
extra_cascade <- function(dList, dList_sc){ 
  param_sc <- lapply(names(dList), function(x) a <- dList_sc[[x]] - dList[[x]])
  names(param_sc) <- names(dList)
  
  return(param_sc)
}


param_sc <- list()
 
scenario_cascade <- list("dfList_NP_2023" = dfList_NP_2023, 
                         "dfList_NP_2024" = dfList_NP_2024, 
                         "dfList_NPexp_A" = dfList_NPexp_A, 
                         "dfList_NPexp_B" = dfList_NPexp_B, 
                         "dfList_NPexp_C" = dfList_NPexp_C, 
                         "dfList_NPexp_D" = dfList_NPexp_D)

for(i in names(scenario_cascade)){
  param_sc[[i]] <- extra_cascade(dfList, scenario_cascade[[i]])
  
  
}

# save param_sc for HCV_model variable: param_cascade_sc
save(param_sc,
     file = file.path(OutputFolder ,
                      paste0(project_name,"param_sc" ,".rda"))) 


# import cost data
files <- list.files(path = paste0(DataFolder, 
                                  "/cost/", sep =  ""), pattern = '*.csv')


costdfList <- lapply(files, function(f) {
  
  df <- read.csv(file.path(paste0(DataFolder, "/cost/", f, sep = "")), header = TRUE)
  
  df <- df[, -1]
  
  df <- df%>%as_tibble()
  
  df <- as.matrix(df, nrow = npops, ncol = length(.) + 1)
  
})

names(costdfList) <- c(gsub("^|.csv", "", files)) # ^: from beginning, \ end before .csv


cost_state <- costdfList$state
costflow <- list()
costflow[[1]] <- costdfList$costFlow
costflow[[2]] <- costdfList$costFlow_POCRNA

costflow_Neg <- list()
costflow_Neg[[1]] <- costdfList$costFlow_NEG
costflow_Neg[[2]] <- costdfList$`costFlow_POCRNA _NEG`


endY <- 100

param_dfList <- lapply(dfList, function(x) x*0)
names(param_dfList) <- names(dfList)


Sce_sq <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                 disease_progress,pop_array,
                 dfList, param_cascade_sc = param_dfList , fib = fib, 
                 modelrun="UN", proj = "POC_AU", end_Y = endY, 
                 cost = costdfList, costflow = costflow, 
                 costflow_Neg = costflow_Neg)

tic <- proc.time()
Sce_np <- list()

for(scenario in names(param_sc)){ 
  Sce_np[[scenario]] <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                   disease_progress,pop_array,
                   scenario_cascade[[scenario]], 
                   param_cascade_sc = param_sc[[scenario]] , fib = fib, 
                   modelrun="UN", proj = "POC_AU", end_Y = endY, 
                   cost = costdfList, costflow = costflow, 
                   costflow_Neg = costflow_Neg)
  
  }

toc <- proc.time() - tic
toc
save(Sce_sq,Sce_np,
     file = file.path(OutputFolder ,
                      paste0(project_name,"Simulations" ,".rda"))) 
