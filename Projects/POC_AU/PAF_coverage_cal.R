# this script is for calculating the effect of national program 
# 1. we first calculated the contribution of national program in treatment uptake 
#    by applying the population attribution fraction (PAF)
# 2. Then We used the PAF to estimate the coverage of national program 
# 3. We applied the coverage and PAF to estimate the parameters for HCV testing 
#    and treatment for national program scenario 
rm(list = ls())
library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)

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


# 1. PAF 
#    2022     | t_init    | w/t_init  | Total
# ----------------------------------------------------
# POC_AU      | 669 (A)   | (B)       | N1 (X)
# ----------------------------------------------------
# w/POC_AU    | 6646 (C)  | (D)       | N2(Y) 
# ----------------------------------------------------
# Total       | N3        | N4        |      
#
# Formula for PAF: (O_i - E_i)/O_i
# exmaple from the cross-table 
# PAF = [(A + C) - (X*C/Y + C)]/(A + C) = (A - X*C/Y)/(A + C)

# Since there are more treatment uptake in the prison, we set the different PAFs 
# for community and prison settings 
# national treatment uptake by C or P is extracted from the cascade report 
# we assumed all the treatment uptake in the comunities is among C_PWID and C_fPWID
# In 2022, 2563 hepatitis C treatments were initiated in prisons across 
# all Australian jurisdictions. 
# This is calculated to represent 35% (2563/7344) of all hepatitis C 
# treatment episodes in Australia in 2022
# PAF_C
#    2022     | t_init    | w/t_init  | Total
# ----------------------------------------------------
# POC_AU      | 112 (A)   | (B)       | N1 (X)
# ----------------------------------------------------
# w/POC_AU    | 4669 (C)  | (D)       | N2(Y) 
# ----------------------------------------------------
# Total       | N3        | N4        |      
# C/Y= treatment coverage in community: 0.59 
# A/X = 0.61 (0.594 in C_PWID, 0.649 in C_fPWID)
PAF_C <- (112 - (112/0.61)*0.59)/(112+4669)

OR_C <- 112/4669

# PAF_P
#    2022     | t_init    | w/t_init  | Total
# ----------------------------------------------------
# POC_AU      | 557 (A)   | (B)       | N1 (X)
# ----------------------------------------------------
# w/POC_AU    | 2560 (C)  | (D)       | N2(Y) 
# ----------------------------------------------------
# Total       | N3        | N4        |  
# C/Y= treatment coverage in community: 0.60 
# A/X = 0.6 

PAF_P <- (557 - (557/0.61)*0.6)/(557+2563)

OR_P <- 557/2563

# 2. national program coverage calculation 
# 2.1 calculation additional % of testing and treatment parameters under national program 
# X index the coverage of national program 
# RR_tauRNA = tauRNA_POCAU/tauRNA_soc
# RR_tauRNAonly = tauRNAonly_POCAU/tauRNAonly_soc
# RR_eta = eta_POCAU/eta_soc
# formula: PAF = [(RR_tauRNA - 1) + (RR_tauRNAonly - 1)]*X*(RR_eta - 1)
# X = PAF/{[(RR_tauRNA - 1) + (RR_tauRNAonly - 1)]*(RR_eta - 1)}
# ab testing in P 
abcov <- dfList$tau_ab[1,2,100]

RR_tauRNA_C <- 0.939/0.686

RR_tauRNAonly_C <- 0.709/0.33

RR_eta_C <- 0.61/0.59

RR_tauRNA_P <- 0.939/0.699

RR_tauRNAonly_P <- 0.709/0.122

RR_eta_P <- 0.61/0.6 

RRlst <- list("C" = list("tau_RNA" = RR_tauRNA_C, 
                         "tau_poct" = RR_tauRNAonly_C,
                         "eta" = RR_eta_C),
              "P" = list("tau_RNA" = RR_tauRNA_P, 
                         "tau_poct" = RR_tauRNAonly_P,
                         "eta" = RR_eta_P))

Ccal_C <- PAF_C/(((RR_tauRNA_C - 1) + (RR_tauRNAonly_C - 1))*(RR_eta_C - 1))

Ccal_P <- PAF_P/(((RR_tauRNA_P - 1) + (RR_tauRNAonly_P - 1))*(RR_eta_P - 1))

Ccal <- list("C" = Ccal_C, "P" = Ccal_P)
Ccal_scale <- lapply(Ccal, function(x) x*10)

#### RR approach #### 
# coverage = OR*[(tau_ab_soc*tau_RNA_soc + tau_poct_soc)*eta_soc]/[(tau_ab_np*tau_RNA_np + tau_poct_np)*eta_np]
Ccal_rr_C <- OR_C*((0.24*0.686 +0.33)*0.59)/((0.24*0.939+0.709)*0.62)
Ccal_rr_P <- OR_P*((0.4*0.699 +0.122)*0.60)/((0.4*0.939+0.709)*0.89)

Ccal <- list("C" = Ccal_rr_C, "P" = Ccal_rr_P)
Ccal_scale <- lapply(Ccal, function(x) x*2)





# 3. testing and treat parameters for national program 

# set up parameters for scenario 
dfList_NP <- dfList 
POC_AU$cabY
S_Yint <- 2022
S_Yend <-2024
# function to estimate the parameters of national program 
Param_cal <- function(pj, dlist, index ,S_Yint, S_Yend, ORlist, Ccal ){ 
  #proj: project name 
  # dlist: list of cascade parameters 
  # index: the parameter to calculate 
  # SYint: starting year of scenario 
  # S_Yend: end year of scenario 
  # Orlist: the oddration of the parameters in scenario (save in a list)
  # Ccal: the list of national program coverage 
  SYpoint_int <- (S_Yint - pj$cabY)/pj$timestep + 1
  SYpoint_end <- (S_Yend - pj$cabY)/POC_AU$timestep + 1
  SY_leng <- SYpoint_end - SYpoint_int
  
  intVal <- dlist[[index]][, 3, SYpoint_int] 
  
  endVal <- c() 
  endVal[1] <- intVal[1] + (ORlist[["C"]][[index]] - 1)*Ccal[["C"]]
  endVal[2] <- intVal[2] + (ORlist[["C"]][[index]] - 1)*Ccal[["C"]]
  endVal[3] <- intVal[3] + (ORlist[["P"]][[index]] - 1)*Ccal[["P"]]
  endVal[4] <- intVal[4] + (ORlist[["P"]][[index]] - 1)*Ccal[["P"]] 
  endVal[5] <- intVal[5] 
  
  endVal[1] <- ifelse(endVal[1]>1, 1, endVal[1])
  endVal[2] <- ifelse(endVal[2]>1, 1, endVal[2])
  endVal[3] <- ifelse(endVal[3]>1, 1, endVal[3])
  endVal[4] <- ifelse(endVal[4]>1, 1, endVal[4])
  
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
      intVal[5]
  }

  
  
  for ( i in 2:dim(dlist[[index]])[[2]]){
    dlist[[index]][1, i, c(SYpoint_int:pj$npts)] <- 
      c(seq(as.numeric(intVal[1]), as.numeric(endVal[1]), length = (SY_leng)), 
        rep(endVal[1], pj$npts - SYpoint_end +1))
    
    dlist[[index]][2, i, c(SYpoint_int:pj$npts)] <- 
      c(seq(as.numeric(intVal[2]), as.numeric(endVal[2]), length = (SY_leng)), 
        rep(endVal[2], pj$npts - SYpoint_end +1))
    
    dlist[[index]][3, i, c(SYpoint_int:pj$npts)] <- 
      c(seq(as.numeric(intVal[3]), as.numeric(endVal[3]), length = (SY_leng)), 
        rep(endVal[3], pj$npts - SYpoint_end +1))
    
    dlist[[index]][4, i, c(SYpoint_int:pj$npts)] <- 
      c(seq(as.numeric(intVal[4]), as.numeric(endVal[4]), length = (SY_leng)), 
        rep(endVal[4], pj$npts - SYpoint_end +1))
    
    
  }
  
  return(dlist[[index]])
}

param_var <- c("tau_RNA", "tau_poct", "eta")

for(i in param_var){ 
  dfList_NP[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP, index = i, 
                         S_Yint = 2022, S_Yend = 2024, ORlist = RRlst, 
                         Ccal = Ccal)
}

save(dfList_NP, Ccal,
     file = file.path(OutputFolder ,
                      paste0(project_name,"S_NP_test" ,".rda")))

dfList_NPscale <- dfList_NP
for(i in param_var){ 
  dfList_NPscale[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NPscale, index = i, 
                              S_Yint = 2024, S_Yend = 2027, ORlist = RRlst, 
                              Ccal = Ccal_scale)
}
save(dfList_NPscale,Ccal_scale,
     file = file.path(OutputFolder ,
                      paste0(project_name,"S_NPscale_test" ,".rda")))

tic <- proc.time()

endY <- 100
Sce_np <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                        disease_progress,pop_array,
                        dfList_NP, fib, 
                        modelrun="UN", proj = "POC_AU", end_Y = endY)

Sce_npscale <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                 disease_progress,pop_array,
                 dfList_NPscale, fib, 
                 modelrun="UN", proj = "POC_AU", end_Y = endY)


toc <- proc.time() - tic 
toc


save(Sce_np,
     file = file.path(OutputFolder ,
                      paste0(project_name,"Sce_np" ,".rda")))

save(Sce_npscale,
     file = file.path(OutputFolder ,
                      paste0(project_name,"Sce_npscale" ,".rda")))



    

