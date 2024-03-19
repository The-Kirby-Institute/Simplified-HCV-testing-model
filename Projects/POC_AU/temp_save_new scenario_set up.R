# new methods 
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
# assuming actively C_PWID and C_fPWID roughly equals to the pop size of CPWID
# We accounted the transition in prison setting regarding its high dynamic nature 
Num_test_person_NP <- Num_test_person_NP%>%
  mutate(num_pop = c(80000, 80000, 80000, 80000))

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

# efficacy of national program 
NP_tauab_C <- 1

NP_tauRNA_C <- 0.939

NP_tauRNAonly_C <- 1

# the treatment initiation is the % of people tested RNA+ initiated DAA within 120 days 
NP_eta_C <-  1- (1- 0.62)^(1/POC_AU$timestep/4)

NP_tauab_P <- 1

NP_tauRNA_P <- 0.939


NP_tauRNAonly_P <- 1

# NP_eta_P <- 0.89 for 1 month 
# turn it back to annual probability
NP_eta_P <- 1- (1-0.89)^(1/POC_AU$timestep)
NPlst <- list("C" = list("tau_ab" = NP_tauab_C,
                         "tau_RNA" = NP_tauRNA_C, 
                         "tau_poct" = NP_tauRNAonly_C,
                         "eta" = NP_eta_C),
              "P" = list("tau_ab" = NP_tauab_P,
                         "tau_RNA" = NP_tauRNA_P, 
                         "tau_poct" = NP_tauRNAonly_P,
                         "eta" = NP_eta_P))

NPlst
Param_cal <- function(pj, dlist, index ,S_Yint, S_Yend, r_Yend, NPlst, 
                      frac_testing = NULL, fp ){ 
  #proj: project name 
  # dlist: list of cascade parameters 
  # index: the parameter to calculate 
  # SYint: starting year of scenario 
  # S_Yend: end year of scenario 
  # NPlst: the odd ratio of the parameters in scenario (save in a list)
 
  # frac_testing: the list of fraction of testing pathways in national program 
  # convert to monthly probability then convert back to yearly probability 
  # fp: np coverage in the people living with HCV RNA+ 
  SYpoint_int <- (S_Yint - pj$cabY)/pj$timestep + 1
  SYpoint_end <- (S_Yend - pj$cabY)/pj$timestep
  
  SY_leng <- SYpoint_end - SYpoint_int + 1 
  # remain the same level to end of which year
  rYpoint_end <- (r_Yend - pj$cabY)/POC_AU$timestep
  
  rY_leng <- rYpoint_end - SYpoint_end
  
  intVal <- dlist[[index]][, 3, (SYpoint_int -1)]
  intVal_dt <- c()
  intVal_dt[1] <- 1-(1- intVal[1])^pj$timestep
  intVal_dt[2] <- 1-(1- intVal[2])^pj$timestep
  intVal_dt[3] <- 1-(1- intVal[3])^pj$timestep
  intVal_dt[4] <- 1-(1- intVal[4])^pj$timestep
  intVal_dt[5] <- 1-(1- intVal[5])^pj$timestep
  
  
  # function to estimate coverage among positive and negative 
 
  scVal_dt <- c()
  if(is.null(frac_testing)& isTRUE(index%in% c("tau_ab", "tau_poct"))){
    scVal_dt[1] <- NPlst[["C"]][[index]]*fp[1]
    scVal_dt[2] <- NPlst[["C"]][[index]]*fp[2]
    scVal_dt[3] <- NPlst[["P"]][[index]]*fp[3]
    scVal_dt[4] <- NPlst[["P"]][[index]]*fp[4]
    scVal_dt[5] <- NPlst[["P"]][[index]]*fp[5]
  }
  else if(is.null(frac_testing)& isTRUE(index%in% c("tau_RNA", "eta"))){ 
    scVal_dt[1] <- NPlst[["C"]][[index]]
    scVal_dt[2] <- NPlst[["C"]][[index]]
    scVal_dt[3] <- NPlst[["P"]][[index]]
    scVal_dt[4] <- NPlst[["P"]][[index]]
    scVal_dt[5] <- NPlst[["P"]][[index]]
    
    }
  
  else if(!is.null(frac_testing) & isTRUE(index%in% c("tau_ab"))){ 
    scVal_dt[1] <- NPlst[["C"]][[index]]*frac_testing[["C"]][["reflex"]]*fp[1]
    scVal_dt[2] <- NPlst[["C"]][[index]]*frac_testing[["C"]][["reflex"]]*fp[2]
    scVal_dt[3] <- NPlst[["P"]][[index]]*frac_testing[["P"]][["reflex"]]*fp[3]
    scVal_dt[4] <- NPlst[["P"]][[index]]*frac_testing[["P"]][["reflex"]]*fp[4]
    scVal_dt[5] <- NPlst[["P"]][[index]]*frac_testing[["P"]][["reflex"]]*fp[5]
    
  }
  else if(!is.null(frac_testing) & isTRUE(index == "tau_poct")){
    
    scVal_dt[1] <- NPlst[["C"]][[index]]*frac_testing[["C"]][["immeRNA"]]*fp[1]
    scVal_dt[2] <- NPlst[["C"]][[index]]*frac_testing[["C"]][["immeRNA"]]*fp[2]
    scVal_dt[3] <- NPlst[["P"]][[index]]*frac_testing[["P"]][["immeRNA"]]*fp[3]
    scVal_dt[4] <- NPlst[["P"]][[index]]*frac_testing[["P"]][["immeRNA"]]*fp[4]
    scVal_dt[5] <- NPlst[["P"]][[index]]*frac_testing[["P"]][["immeRNA"]]*fp[5]
    
  }
  else if(!is.null(frac_testing) & isTRUE(index %in% c("tau_RNA"))){
    
    scVal_dt[1] <- NPlst[["C"]][[index]]
    scVal_dt[2] <- NPlst[["C"]][[index]]
    scVal_dt[3] <- NPlst[["P"]][[index]]
    scVal_dt[4] <- NPlst[["P"]][[index]]
    scVal_dt[5] <- NPlst[["P"]][[index]]
    
  } else if(!is.null(frac_testing) & isTRUE(index %in% c("eta"))){
    
    scVal_dt[1] <- NPlst[["C"]][[index]]
    scVal_dt[2] <- NPlst[["C"]][[index]]
    scVal_dt[3] <- NPlst[["P"]][[index]]
    scVal_dt[4] <- NPlst[["P"]][[index]]
    scVal_dt[5] <- NPlst[["P"]][[index]]
    
  }

  
      
  
  
  
  if(isTRUE(rYpoint_end > SYpoint_end)){ 
    for ( i in 2:dim(dlist[[index]])[[2]]){
      dlist[[index]][1, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[1]), as.numeric(scVal_dt[1]), length = (SY_leng + 1)), 
          rep(as.numeric(scVal_dt[1]), length = (rY_leng)),
          rep(intVal[1], pj$npts - rYpoint_end))
      
      dlist[[index]][2, i, c((SYpoint_int - 1):pj$npts)] <- 
          c(seq(as.numeric(intVal[2]), as.numeric(scVal_dt[2]), length = (SY_leng + 1)), 
            rep(as.numeric(scVal_dt[2]), length = (rY_leng)),
            rep(intVal[2], pj$npts - rYpoint_end))
      
      dlist[[index]][3, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[3]), as.numeric(scVal_dt[3]), length = (SY_leng + 1)), 
          rep(as.numeric(scVal_dt[3]), length = (rY_leng)),
          rep(intVal[3], pj$npts - rYpoint_end))
      
      dlist[[index]][4, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[4]), as.numeric(scVal_dt[4]), length = (SY_leng + 1)), 
          rep(as.numeric(scVal_dt[4]), length = (rY_leng)),
          rep(intVal[4], pj$npts - rYpoint_end))
      
      dlist[[index]][5, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[5]), as.numeric(scVal_dt[5]), length = (SY_leng + 1)), 
          rep(as.numeric(scVal_dt[5]), length = (rY_leng)),
          rep(intVal[5], pj$npts - rYpoint_end))
    }
    
  } else{ 
    for ( i in 2:dim(dlist[[index]])[[2]]){
      dlist[[index]][1, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[1]), as.numeric(scVal_dt[1]), length = (SY_leng + 1)), 
          rep(intVal[1],  pj$npts - SYpoint_end))
      
        dlist[[index]][2, i, c((SYpoint_int - 1):pj$npts)] <- 
            c(seq(as.numeric(intVal[2]), as.numeric(scVal_dt[2]), length = (SY_leng + 1)), 
              rep(intVal[2], pj$npts - SYpoint_end))
      
      dlist[[index]][3, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[3]), as.numeric(scVal_dt[3]), length = (SY_leng + 1)), 
          rep(intVal[3], pj$npts - SYpoint_end))
      
      dlist[[index]][4, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[4]), as.numeric(scVal_dt[4]), length = (SY_leng + 1)), 
          rep(intVal[4], pj$npts - SYpoint_end))
      
      dlist[[index]][5, i, c((SYpoint_int - 1):pj$npts)] <- 
        c(seq(as.numeric(intVal[5]), as.numeric(scVal_dt[5]), length = (SY_leng + 1)), 
          rep(intVal[5], pj$npts - SYpoint_end))
    }
    
  }
  
  
  
  
  return(dlist[[index]])
}


dfList_NP <- lapply(dfList, function(x) x*0)

param_var <- c("tau_ab","tau_RNA", "tau_poct", "eta") 
for(i in param_var){  
    dfList_NP[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP, index = i, 
                                frac_testing = frac_test[[2022]],
                                S_Yint = 2022, S_Yend = 2023, r_Yend = 2023, 
                                NPlst = NPlst,
                                fp = c(Ccal[[2023]]$C*1.2, 
                                       Ccal[[2023]]$C*1.2, 
                                       Ccal[[2022]]$P, 
                                       Ccal[[2022]]$P, 
                                       Ccal[[2022]]$P))
    
    dfList_NP[[i]][, , 85:96] <- dfList_NP[[i]][, ,96]
  }


dfList_NP_2023 <- dfList_NP

Ccal_2023 <- Ccal[[2023]]
frac_test_2023 <- frac_test[[2023]]
for(i in param_var){
    dfList_NP_2023[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP, index = i, frac_testing = frac_test[[2023]],
                                          S_Yint = 2023, S_Yend = 2024, r_Yend = 2024, NPlst = NPlst, 
                                          Ccal = Ccal[[2023]], 
                                     fp = c(Ccal[[2023]]$C*1.2*2.611835, 
                                            Ccal[[2023]]$C*1.2*2.611835, 
                                            1.2*Ccal[[2022]]$P, 
                                            1.2*Ccal[[2022]]$P, 
                                            Ccal[[2023]]$P))
    

    dfList_NP_2023[[i]][, , 97:108] <- dfList_NP_2023[[i]][, , 108]
    
}




for(i in param_var){
    # begining of 2024 
    b_pt <- (2024 - POC_AU$cabY)/POC_AU$timestep + 1 
    
    # length of the time points
    dim_length <- dim(dfList_NP_2023[[i]])[3]
    dfList_NP_2023[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
  } 

save(dfList_NP_2023, Ccal_2023,frac_test_2023,
     file = file.path(OutputFolder ,
                      paste0(project_name,"NP_2023" ,".rda")))

#### Scenario 2: current achievement + 2024 prediction: back to pre-national program level since starting of 2025  ####  
# odds of coverage: 1.08 (total test_predicited/ number test in 2023, 19480/18035) 
odd_num_test <- list()


odd_num_test[[2024]] <- 20000/11276
Ccal[[2024]] <- lapply(Ccal[[2023]],function(x) x*odd_num_test[[2024]])
frac_test[[2024]] <- frac_test[[2023]]

dfList_NP_2024 <- dfList_NP_2023

for(i in param_var){ 
  dfList_NP_2024[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2023, index = i, 
                                   frac_testing = frac_test[[2024]], S_Yint = 2024, 
                                   S_Yend = 2025, r_Yend = 2025, NPlst = NPlst, 
                                   Ccal = Ccal[[2024]], 
                                   fp = c(0.0601375*1.2*2.611835*odd_num_test[[2024]]*1.35, 
                                          0.0601375*1.2*2.611835*odd_num_test[[2024]]*1.35, 
                                          1.2*0.05195*odd_num_test[[2024]]*1.35, 
                                          1.2*0.05195*odd_num_test[[2024]]*1.35, 
                                          Ccal[[2023]]$P*odd_num_test[[2024]]))
}


# back to pre-national program’s level from starting of 2025
for(i in param_var){ 
  # begining of 2025
  b_pt <- (2025 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2024[[i]])[3]
  dfList_NP_2024[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
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
                                   frac_testing = frac_test[[2024]], S_Yint = 2024, 
                                   S_Yend = 2025, r_Yend = 2028, NPlst = NPlst, 
                                   Ccal = Ccal[[2024]],
                                   fp = c(0.0601375*1.2*2.611835*odd_num_test[[2024]]*1.35, 
                                          0.0601375*1.2*2.611835*odd_num_test[[2024]]*1.35, 
                                          1.2*0.05195*odd_num_test[[2024]]*1.35, 
                                          1.2*0.05195*odd_num_test[[2024]]*1.35, 
                                          Ccal[[2023]]$P*odd_num_test[[2024]]))
  
}


# back to pre-national program’s level from starting of 2028
for(i in param_var){ 
  # begining of 2028
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPexp_A[[i]])[3]
  dfList_NPexp_A[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
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
                                   frac_testing = frac_test[[2024]],
                                   S_Yint = 2025, S_Yend = 2026, r_Yend = 2028,
                                   NPlst = NPlst, 
                                   Ccal = Ccal_2025, 
                                   fp = c(0.0601375*1.2*2.611835*odd_num_test[[2024]]*1.35*odd_num_test[[2025]], 
                                          0.0601375*1.2*2.611835*odd_num_test[[2024]]*1.35*odd_num_test[[2025]], 
                                          1.2*0.05195*odd_num_test[[2024]]*1.35*odd_num_test[[2025]], 
                                          1.2*0.05195*odd_num_test[[2024]]*1.35*odd_num_test[[2025]], 
                                          Ccal[[2023]]$P*odd_num_test[[2024]]*odd_num_test[[2025]]))
}



# back to pre-national program’s level from starting of 2028
for(i in param_var){ 
  # begining of 2028
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPexp_B[[i]])[3]
  dfList_NPexp_B[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
}

save(dfList_NPexp_B, Ccal_2025,
     file = file.path(OutputFolder ,
                      paste0(project_name,"NP_expand_B" ,".rda")))


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
                                   Ccal = Ccal_2026, Prev = Prev2, fc = fc_sc)
}

# back to pre-national program’s level from starting of 2028
for(i in param_var){ 
  # begining of 2028
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPexp_C[[i]])[3]
  dfList_NPexp_C[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
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
                                   Ccal = Ccal_2027, Prev = Prev2, fc = fc_sc)
}

# back to pre-national program’s level from starting of 2028
for(i in param_var){ 
  # begining of 2028
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPexp_D[[i]])[3]
  dfList_NPexp_D[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
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
                         "dfList_NP_2024" = dfList_NP_2024)

for(i in names(scenario_cascade)){
  param_sc[[i]] <- scenario_cascade[[i]]
  
  
}

param_sc$dfList_NP_2023$tau_poct[, , 108]
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
fc <- matrix(0, ncol = dim(dfList$eta)[3], nrow = POC_AU$npops)

Sce_sq <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                 disease_progress,pop_array,
                 dfList, param_cascade_sc = param_dfList , fib = fib, 
                 modelrun="UN", proj = "POC_AU", end_Y = endY, 
                 cost = costdfList, costflow = costflow, 
                 costflow_Neg = costflow_Neg, fc_sc = fc,
                 fp = NULL)

tic <- proc.time()
Sce_np <- list()
ncol(fc)
# numb_ab ==0 , raplaced by immRNA
# frac_ab ==0, then frac_ab = 1 
fc <- matrix(0, ncol = dim(dfList$eta)[3], nrow = POC_AU$npops)
fp <- c()
num_ab <- c()
cov_np <- c()
frac_ab <- c()
fp <- c()
coverage <- c()
fs_estimate <- function(num_ab, cov_np, frac_ab,fp, year, coverage, modsim, endY){ 
  # num_ab in year <- c(num_ab_c, num_ab_c, num_ab_p, num_ab_p, num_ab_p)
  # fp is the factor multiplied the coverage of national program 
  if(is.null(endY)){ 
    
    endY <- 100}
  t_init <- (year - POC_AU$cabY)/POC_AU$timestep + 1
  t_end <- (year - POC_AU$cabY + 1)/POC_AU$timestep
  undiag <- modres.t(POC_AU, modsim, endYear = endY, allp = NULL)%>%
    filter(cascade == "undiag" & disease_prog != "a")%>%
    mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), "C", "P"))%>%
   group_by(timestep, setting, cascade)%>%
    summarize(best = sum(best))%>%ungroup()%>%group_by(setting)%>%
    filter(timestep >= year - POC_AU$cabY & timestep<year - POC_AU$cabY + 1)%>%
    head(., n = 2)
  
  s_bar <- modres.t(POC_AU, modsim, endYear = endY, allp = NULL)%>%
    filter(cascade %in% c("s", "cured") & disease_prog != "a")%>%
    mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), "C", "P"))%>%
    group_by(timestep, setting)%>%
    summarize(best = sum(best))%>%ungroup()%>%group_by(setting)%>%
    filter(timestep >= year - POC_AU$cabY & timestep<year - POC_AU$cabY + 1)%>%ungroup()%>%
    head(., n = 2)

  fs <- matrix(0, ncol = dim(dfList$eta)[3], nrow = POC_AU$npops)
  
  fs[1, t_init: t_end] <- 
      rep((num_ab[1]/(cov_np[1]*frac_ab[1]) - fp[1]*unlist(undiag[1, "best"]))/unlist(s_bar[1, "best"]), 
          t_end - t_init + 1)
  
  fs[2, t_init: t_end] <- 
    rep((num_ab[1]/(cov_np[1]*frac_ab[1]) - fp[2]*unlist(undiag[1, "best"]))/unlist(s_bar[1, "best"]), 
        t_end - t_init + 1)
  
  fs[3, t_init: t_end] <- 
    rep((num_ab[2]/(cov_np[2]*frac_ab[2]) - fp[3]*unlist(undiag[2, "best"]))/unlist(s_bar[2, "best"]), 
        t_end - t_init + 1)
  
  fs[4, t_init: t_end] <- 
    rep((num_ab[2]/(cov_np[2]*frac_ab[2]) - fp[4]*unlist(undiag[2, "best"]))/unlist(s_bar[2, "best"]), 
        t_end - t_init + 1)
  
  fs[5, t_init: t_end] <- rep((num_ab[2]/(cov_np[2]*frac_ab[2]) - fp[5]*unlist(undiag[2, "best"]))/unlist(s_bar[2, "best"]), 
                              t_end - t_init + 1)
    
    

  return(list(fs, undiag, s_bar))

  }
coverage
n_ab <- c(873, 4515)

frac_ab <- c(frac_test[[2022]]$C$reflex,  
             frac_test[[2022]]$P$immeRNA )

fm <- c(1.3,1.3,8.5,8.5,1.5)
coverage_np <- c(0.023, 0.05195)
xfs <- fs_estimate(num_ab = n_ab, cov_np = coverage_np, frac_ab = frac_ab, 
            fp = fm, year = 2022, modsim = Sce_sq)
fs <- xfs[[1]] 
fp
xfs[[3]]
fs[1,] <-   xfs[[1]][1, ]/fm[1]
fs[2,] <-   xfs[[1]][2, ]/fm[2]
fs[3,] <-   xfs[[1]][3, ]/fm[3]
fs[4,] <-   xfs[[1]][4, ]/fm[4]
fs[5,] <-   xfs[[1]][5, ]/fm[5]
dfList_NP <- lapply(dfList, function(x) x*0)
fs[,96]
xfs[[1]][, 96]
param_var <- c("tau_ab","tau_RNA", "tau_poct", "eta") 
for(i in param_var){  
  dfList_NP[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP, index = i, 
                              frac_testing = frac_test[[2022]],
                              S_Yint = 2022, S_Yend = 2023, r_Yend = 2023, NPlst = NPlst, 
                              Ccal = Ccal[[2022]], fp = c(0.023*fm[1], 
                                                          0.023*fm[2], 
                                                          0.05195*fm[3], 
                                                          0.05195*fm[4], 
                                                          Ccal[[2022]]$P))
  
  dfList_NP[[i]][, , 85:96] <- dfList_NP[[i]][, ,96]
}
dfList_NP_2023 <- dfList_NP
fm2 <- c(0.7,0.6,10,10,1.5)
for(i in param_var){  
  dfList_NP_2023[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP, index = i, 
                              frac_testing = frac_test[[2023]],
                              S_Yint = 2023, S_Yend = 2024, r_Yend = 2024, NPlst = NPlst, 
                              Ccal = Ccal[[2023]], fp = c(Ccal[[2023]]$C*fm2[1], 
                                                          Ccal[[2023]]$C*fm2[2], 
                                                          Ccal[[2023]]$P*fm2[3], 
                                                          Ccal[[2023]]$P*fm2[4], 
                                                          Ccal[[2023]]$P))
  
  dfList_NP_2023[[i]][, , 97:108] <- dfList_NP_2023[[i]][, ,108]
}

for(i in param_var){
  # begining of 2024 
  b_pt <- (2024 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2023[[i]])[3]
  dfList_NP_2023[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

frac_test_2023 <- c(0.5,0.1423649)
n_ab <- c(2406,771)
Ccal_2023
coverage_np <- c(0.0601375, 0.0580375)
xfs2 <- fs_estimate(num_ab = n_ab, cov_np = coverage_np, frac_ab = frac_test_2023, 
                   fp = fm2, year = 2023, modsim = Sce_sq)
xfs2[[1]][, 108]
fp
xfs[[3]]
fs[1, c(97:108)] <-   xfs2[[1]][1, c(97:108)]/fm2[1]
fs[2,c(97:108)] <-   xfs2[[1]][2, c(97:108)]/fm2[2]
fs[3,c(97:108)] <-   xfs2[[1]][3, c(97:108)]/fm2[3]
fs[4,c(97:108)] <-   xfs2[[1]][4, c(97:108)]/fm2[4]
fs[5,c(97:108)] <-   xfs2[[1]][5, c(97:108)]/fm2[5]
dfList_NP_2024 <- dfList_NP_2023
fract_test_2024 <- c(0.6, 0.2)
frac_test[[2024]] <- frac_test[[2023]]
frac_test[[2024]]$C$reflex <- 0.6 
frac_test[[2024]]$C$immeRNA <- 1 - 0.6
frac_test[[2024]]$P$reflex <-0.2 
frac_test[[2024]]$P$immeRNA <- 1 - 0.2
n_ab <- c(6660,1771)
fm3 <- c(0.6,0.6,10,10,5.5)
fs[, 120]
Ccal[[2024]] <- lapply(Ccal[[2023]], function(x) x*2)
for(i in param_var){  
  dfList_NP_2024[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2024, index = i, 
                                   frac_testing = frac_test[[2024]],
                                   S_Yint = 2024, S_Yend = 2025, r_Yend = 2025, NPlst = NPlst, 
                                   Ccal = Ccal[[2024]], fp = c(Ccal[[2024]]$C*fm3[1], 
                                                               Ccal[[2024]]$C*fm3[2], 
                                                               Ccal[[2024]]$P*fm3[3], 
                                                               Ccal[[2024]]$P*fm3[4], 
                                                               Ccal[[2024]]$P*fm3[5]))
  
}

for(i in param_var){
  # begining of 2024 
  b_pt <- (2025 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2024[[i]])[3]
  dfList_NP_2024[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

fract_test_2024
coverage_np <- c(Ccal[[2024]]$C, Ccal[[2024]]$P)
xfs3 <- fs_estimate(num_ab = n_ab, cov_np = coverage_np, frac_ab = fract_test_2024, 
                    fp = fm3, year = 2024, modsim = xxx, endY = 20)
xfs3[[1]][, 120]
fp

fs[1, c(109:120)] <-   xfs3[[1]][1, c(109:120)]/fm3[1]


fs[2,c(109:120)] <-   xfs3[[1]][2, c(109:120)]/fm3[2]


fs[3,c(109:120)] <-   xfs3[[1]][3, c(109:120)]/fm3[3]
fs[4,c(109:120)] <-   xfs3[[1]][4, c(109:120)]/fm3[4]
fs[5,c(109:120)] <-   xfs3[[1]][5, c(109:120)]/fm3[5]
fs[, 120]
fp_m <- c(1,1,1,1,1)
xxx <- HCVMSM(POC_AU, best_estimates, best_est_pop,
       disease_progress,pop_array,
       dfList, 
       param_cascade_sc = dfList_NP_2024, fib = fib, 
       modelrun="UN", proj = "POC_AU", end_Y = 20,  
       fc_sc= fs,
       fp = fp_m) 


test <- list()
cl_ext <- names(xxx)[c(13:16, 20:22)]
for(i in cl_ext){
  
  test[[i]] <- modres.flow.t(POC_AU, xxx, endYear = 20, 
                             allp = i)%>%filter(year %in% seq(7, 15, 1))%>%
    ungroup()%>%
    group_by(year, population)%>%
    summarise(best = sum(best))
  
  
}

test$newTestingPOCT_sc
test <- dplyr::bind_rows(test, .id = 'index')%>%group_by(year, population)%>%spread(index, best)


View(test)
test_fscal <- test%>%mutate(Ab = newTestingAb_sc + newTestingAb_sc_neg, 
                            RNA = (newTestingAg_sc+ newTestingAg_sc_neg + newTestingPOCT_sc + 
                                     newTestingPOCT_sc_neg),
                            treatment = newTreatment_sc)%>%select(year, population, Ab, 
                                                                  RNA, treatment)%>%
  mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), "C", "P"))%>%
  ungroup()%>%
  select(-c(population))%>%
  gather(index, value, -c(year, setting))%>%
  group_by(year, setting, index)%>%summarise(value = sum(value))%>%
  filter(year%in% c(7,8,9, 10))
test_fscal
dfList_NP$tau_ab[,, 96]
xfs[, 96]
fc_t <- xfs[[1]][, 96]/fp[1]
x<- modres.t(POC_AU, Sce_sq, endYear = 100, allp = NULL)%>%
  filter(cascade %in% c("undiag"))%>%
  group_by(timestep, population)%>%
  summarize(best = sum(best))%>%ungroup()%>%group_by(population)%>%
  filter(timestep >= 7 & timestep<8)%>%
  tail(., n = POC_AU$npops)
  
  tail(., n= POC_AU$npops)
x[1, "best"]




fc_scenario <- list()
fc_scenario[["dfList_NP_2023"]] <- fc

fc_scenario[["dfList_NP_2023"]][1,c(85:108)] <- 0.009
fc_scenario[["dfList_NP_2023"]][2,c(85:108)] <- 0.009
fc_scenario[["dfList_NP_2023"]][3,c(85:108)] <- 0.3
fc_scenario[["dfList_NP_2023"]][4,c(85:108)] <- 0.3
fc_scenario[["dfList_NP_2023"]][5,c(85:108)] <- 0.05

fc_scenario[["dfList_NP_2024"]] <- fc_scenario[["dfList_NP_2023"]]


fc_scenario[["dfList_NP_2024"]][1,c(109:120)] <- fc_scenario[["dfList_NP_2024"]][1,108]*odd_num_test[[2024]]*1.35
fc_scenario[["dfList_NP_2024"]][2,c(109:120)] <- fc_scenario[["dfList_NP_2024"]][2,108]*odd_num_test[[2024]]*1.35
fc_scenario[["dfList_NP_2024"]][3,c(109:120)] <- fc_scenario[["dfList_NP_2024"]][3,108]
fc_scenario[["dfList_NP_2024"]][4,c(109:120)] <- fc_scenario[["dfList_NP_2024"]][4,108]
fc_scenario[["dfList_NP_2024"]][5,c(109:120)] <- fc_scenario[["dfList_NP_2024"]][5,108]

fc_scenario[["dfList_NPexp_A"]] <- fc_scenario[["dfList_NP_2024"]]
fc_scenario[["dfList_NPexp_A"]][1,c(121:172)] <- fc_scenario[["dfList_NPexp_A"]][1, 120]
fc_scenario[["dfList_NPexp_A"]][2,c(121:172)] <- fc_scenario[["dfList_NPexp_A"]][2, 120]
fc_scenario[["dfList_NPexp_A"]][3,c(121:172)] <- fc_scenario[["dfList_NPexp_A"]][3, 120]
fc_scenario[["dfList_NPexp_A"]][4,c(121:172)] <- fc_scenario[["dfList_NPexp_A"]][4, 120]
fc_scenario[["dfList_NPexp_A"]][5,c(121:172)] <- fc_scenario[["dfList_NPexp_A"]][5, 120]

fc_scenario[["dfList_NPexp_B"]] <- fc_scenario[["dfList_NP_2024"]]
fc_scenario[["dfList_NPexp_B"]][1,c(121:144)] <- fc_scenario[["dfList_NPexp_B"]][1, 120]*odd_num_test[[2025]]
fc_scenario[["dfList_NPexp_B"]][2,c(121:144)] <- fc_scenario[["dfList_NPexp_B"]][2, 120]*odd_num_test[[2025]]
fc_scenario[["dfList_NPexp_B"]][3,c(121:144)] <- fc_scenario[["dfList_NPexp_B"]][3, 120]*odd_num_test[[2025]]
fc_scenario[["dfList_NPexp_B"]][4,c(121:144)] <- fc_scenario[["dfList_NPexp_B"]][4, 120]*odd_num_test[[2025]]
fc_scenario[["dfList_NPexp_B"]][5,c(121:144)] <- fc_scenario[["dfList_NPexp_B"]][5, 120]*odd_num_test[[2025]]

fc_scenario[["dfList_NPexp_B"]][1,c(145:172)] <- fc_scenario[["dfList_NPexp_B"]][1, 144]*odd_num_test[[2025]]
fc_scenario[["dfList_NPexp_B"]][2,c(145:172)] <- fc_scenario[["dfList_NPexp_B"]][2, 144]
fc_scenario[["dfList_NPexp_B"]][3,c(145:172)] <- fc_scenario[["dfList_NPexp_B"]][3, 144]
fc_scenario[["dfList_NPexp_B"]][4,c(145:172)] <- fc_scenario[["dfList_NPexp_B"]][4, 144]
fc_scenario[["dfList_NPexp_B"]][5,c(145:172)] <- fc_scenario[["dfList_NPexp_B"]][5, 144]

num_ab <- list()
ab_ndt <- matrix(0, ncol = dim(dfList$eta)[3], nrow = POC_AU$npops)

num_ab[["dfList_NP_2023"]] <- ab_ndt
num_ab[["dfList_NP_2023"]]


for(scenario in names(param_sc)){
  
  
  Sce_np[[scenario]] <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                               disease_progress,pop_array,
                               dfList, 
                               param_cascade_sc = param_sc[[scenario]] , fib = fib, 
                               modelrun="UN", proj = "POC_AU", end_Y = endY, 
                               cost = costdfList, costflow = costflow, 
                               costflow_Neg = costflow_Neg, A,
                               fp = fp)
  
}

toc <- proc.time() - tic
toc
save(Sce_sq,Sce_np,
     file = file.path(OutputFolder ,
                      paste0(project_name,"Simulations" ,".rda"))) 








Sce_np <- list()
fp <- c(1,1,1,1,1)

dim(dfList_NP_2023$eta) 
fc <- matrix(0, ncol = dim(dfList$eta)[3], nrow = POC_AU$npops)

fc[1,c(85:108)] <- 0.009
fc[2,c(85:108)] <- 0.009
fc[3,c(85:108)] <- 0.3
fc[4,c(85:108)] <- 0.3
fc[5,c(85:108)] <- 0.05
odd_num_test[[2024]]

fc[1,c(109:120)] <- fc[1,108]*odd_num_test[[2024]]*1.35
fc[2,c(109:120)] <- fc[2,108]*odd_num_test[[2024]]*1.35
fc[3,c(109:120)] <- fc[3,108]
fc[4,c(109:120)] <- fc[4,108]
fc[5,c(109:120)] <- fc[5,108]

fc[1,c(121:144)] <- fc[1,120]*odd_num_test[[2025]]
fc[2,c(121:144)] <- fc[2,120]*odd_num_test[[2025]]
fc[3,c(121:144)] <- fc[3,120]
fc[4,c(121:144)] <- fc[4,120]
fc[5,c(121:144)] <- fc[5,144]

endY <- 20
  xxx <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                        disease_progress,pop_array,
                        dfList, 
                        param_cascade_sc = dfList_NP_2023, fib = fib, 
                        modelrun="UN", proj = "POC_AU", end_Y = endY,  
                fc_sc= fc,
                   fp = fp)
  
  xxx2025 <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                disease_progress,pop_array,
                dfList, 
                param_cascade_sc = dfList_NP_2024, fib = fib, 
                modelrun="UN", proj = "POC_AU", end_Y = endY,  
                fc_sc = fc_scenario$dfList_NP_2024,
                fp = fp)
fc_scenario$dfList_NP_2024[, 120]
  
xxx$newTestingAb_sc_neg[, 81]
  

modres.flow.t(POC_AU, xxx2024, endYear = endY, 
                allp = "newTreatment_sc")%>%
    group_by(year, population)%>%
    summarize(best = sum(best))%>%filter(year>= 7 & year<10)

dfList_NP_2024$tau_ab[, , 120]
    
modres.flow.t(POC_AU, xxx2024, endYear = endY, 
              allp = "newTestingPOCT_sc_neg")%>%
  group_by(year, population)%>%
  summarize(best = sum(best))%>%filter(year>= 7 & year<10)
modres.t(POC_AU, xxx, endYear = endY, 
              allp = NULL)%>%filter(cascade == "undiag")%>%
  group_by(year, population)%>%
  summarize(best = sum(best))%>%filter(year>= 7 & year<9)

test <- list()
cl_ext <- names(xxx)[c(13:16, 20:22)]
for(i in cl_ext){
  
  test[[i]] <- modres.flow.t(POC_AU, xxx2025, endYear = endY, 
                             allp = i)%>%filter(year %in% seq(7, 15, 1))%>%
    ungroup()%>%
    group_by(year, population)%>%
    summarise(best = sum(best))
  
  
}

test$newTestingPOCT_sc
test <- dplyr::bind_rows(test, .id = 'index')%>%group_by(year, population)%>%spread(index, best)

test_fscal <- test%>%mutate(Ab = newTestingAb_sc + newTestingAb_sc_neg, 
                            RNA = (newTestingAg_sc+ newTestingAg_sc_neg + newTestingPOCT_sc + 
                                     newTestingPOCT_sc_neg),
                            treatment = newTreatment_sc)%>%select(year, population, Ab, 
                                                                  RNA, treatment)%>%
  mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), "C", "P"))%>%
  ungroup()%>%
  select(-c(population))%>%
  gather(index, value, -c(year, setting))%>%
  group_by(year, setting, index)%>%summarise(value = sum(value))%>%
  filter(year%in% c(9))
test_fscal
NP_cas_num <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/NP_cascade_num.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(year = year - 2015)%>%arrange(year, setting) 


xtest <- cbind(test_fscal%>%ungroup, dt = NP_cas_num$value)%>%
  mutate(diff = test_fscal$value - dt)%>%
  mutate(diff_sum = cumsum(diff^2))
xtest 

num_testntreatment <- ggplot(test_fscal, aes(x = year, y = value)) + 
  geom_line() + 
  geom_point(data = NP_cas_num, aes(x = year, y = value)) + 
  facet_wrap(.~setting+index, scales="free", ncol = 3) + 
  geom_smooth(data = NP_cas_num %>% group_by(setting,index),    # subset
              method = "lm", formula = y ~ x) +
  scale_x_continuous(limits = c(7, 10), breaks = seq(7,10,1), 
                     labels = seq(2022, 2025,1)) + 
  theme_bw() 
num_testntreatment <- num_testntreatment + 
  facet_custom (~setting+index,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 10000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 10000))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 500))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 2000))),
                    scale_new(5,
                              scale_y_continuous(limits = 
                                                   c(0, 10000))),
                    scale_new(6,
                              scale_y_continuous(limits = 
                                                   c(0, 1000)))
                  )) + theme_bw() 

num_testntreatment


cl_ext <- names(xxx)[c(13:16, 20:22)]
pop <- list()
for(i in cl_ext){
  pop[[i]] <- modres.flow.t(POC_AU, xxx, endYear = endY, allp = i)%>%
    group_by(year, population)%>%summarize(best = sum(best))%>%filter(year>= 7 & year<9)
}
pop <- dplyr::bind_rows(pop, .id = 'index')%>%spread(index, best)
View(pop)


modres.t(POC_AU, xxx, endYear = endY, allp = NULL)%>%filter(cascade == "undiag" & disease_prog != "a")%>%
  group_by(year, population, cascade)%>%summarize(best = sum(best))%>%filter(year>= 7 & year<9)

modres.flow.t(POC_AU, xxx, endYear = endY, 
              allp = "newTreatment")%>%
  group_by(year, population)%>%
  mutate(year  = year + POC_AU$cabY)%>%
  summarize(best = sum(best))%>%filter(year>= 2022 )%>%
  ungroup()%>%mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), "C", "P"))%>%
  group_by(year, setting)%>%summarize(best = sum(best))
