# this script is tiding up the scenario set up 
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

# run pre-national program scenario first 
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


# import the number of test 

Num_test_person_NP <- read_excel(paste0(data_path, "/01. DATA/Num_test_person_NP.xlsx"))
# assuming actively C_PWID and C_fPWID roughly equals to the pop size of CPWID
# We accounted the transition in prison setting regarding its high dynamic nature 
Num_test_person_NP <- Num_test_person_NP%>%na.omit()%>%mutate(num_pop = c(80000, 80000, 80000, 80000))

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
Ccal[[2022]]
n_ab <- read_excel(paste0(data_path, "/01. DATA/n_ab.xlsx"))

# efficacy of national program 

np_effect <- read_excel(paste0(data_path, "/01. DATA/effect_np.xlsx"))

NP_tauab_C <- unlist(as.numeric(np_effect[1,2]))

NP_tauRNA_C <- unlist(as.numeric(np_effect[2,2]))

NP_tauRNAonly_C <- unlist(as.numeric(np_effect[3,2]))

# the treatment initiation is the % of people tested RNA+ initiated DAA within 120 days 
NP_eta_C <-  1- (1- unlist(as.numeric(np_effect[4,2])))^(1/POC_AU$timestep/4)

NP_tauab_P <- unlist(as.numeric(np_effect[1,3]))

NP_tauRNA_P <- unlist(as.numeric(np_effect[2,3]))


NP_tauRNAonly_P <- unlist(as.numeric(np_effect[3,3]))

# NP_eta_P <- 0.89 for 1 month 
# turn it back to annual probability
NP_eta_P <- 1- (1-unlist(as.numeric(np_effect[4,3])))^(1/POC_AU$timestep/1)
NPlst <- list("C" = list("tau_ab" = NP_tauab_C,
                         "tau_RNA" = NP_tauRNA_C, 
                         "tau_poct" = NP_tauRNAonly_C,
                         "eta" = NP_eta_C),
              "P" = list("tau_ab" = NP_tauab_P,
                         "tau_RNA" = NP_tauRNA_P, 
                         "tau_poct" = NP_tauRNAonly_P,
                         "eta" = NP_eta_P))

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
  
  scVal_dt[1] <- ifelse(scVal_dt[1] >= 1, 1, scVal_dt[1])
  scVal_dt[2] <- ifelse(scVal_dt[2] >= 1, 1, scVal_dt[2])
  scVal_dt[3] <- ifelse(scVal_dt[3] >= 1, 1, scVal_dt[3])
  scVal_dt[4] <- ifelse(scVal_dt[4] >= 1, 1, scVal_dt[4])
  scVal_dt[5] <- ifelse(scVal_dt[5] >= 1, 1, scVal_dt[5])
  
  
  
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


# numb_ab ==0 , raplaced by immRNA
# frac_ab ==0, then frac_ab = 1 
fc <- matrix(0, ncol = dim(dfList$eta)[3], nrow = POC_AU$npops)

fs_estimate <- function(num_ab, cov_np, frac_ab,fp, year, coverage, modsim, endY){ 
  # num_ab in year <- c(num_ab_c, num_ab_c, num_ab_p, num_ab_p, num_ab_p)
  # fp is the factor multiplied the coverage of national program 
  if(is.null(endY)){ 
    
    endY <- 100}
  t_init <- (year - POC_AU$cabY)/POC_AU$timestep + 1
  t_end <- (year + 1 - POC_AU$cabY)/POC_AU$timestep
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
  
  fs[5, t_init: t_end] <- 
    rep((num_ab[2]/(cov_np[2]*frac_ab[2]) - fp[5]*unlist(undiag[2, "best"]))/unlist(s_bar[2, "best"]), 
                              t_end - t_init + 1)
  
  
  
  return(list(fs, undiag, s_bar))
  
}

dfList_NP <- lapply(dfList, function(x) x*0)

param_var <- c("tau_ab","tau_RNA", "tau_poct", "eta") 
Ccal[[2022]]
n_ab_np <- list()
n_ab_np[["2022"]] <- c(unlist(as.numeric(n_ab[n_ab$Year == 2022,"community"])), 
                       unlist(as.numeric(n_ab[n_ab$Year == 2022,"prison"])))
n_ab_np[["2023"]] <- c(unlist(as.numeric(n_ab[n_ab$Year == 2023,"community"])), 
                       unlist(as.numeric(n_ab[n_ab$Year == 2023,"prison"])))
n_ab_np[["2024"]] <- c(unlist(as.numeric(n_ab[n_ab$Year == 2024,"community"])), 
                       unlist(as.numeric(n_ab[n_ab$Year == 2024,"prison"])))

frac_ab <- list()
frac_ab[["2022"]] <- c(frac_test[[2022]]$C$reflex, frac_test[[2022]]$P$immeRNA )
frac_ab[["2023"]] <- c(frac_test[[2023]]$C$reflex, frac_test[[2023]]$P$reflex )
frac_ab[["2024"]] <- c(frac_test[[2023]]$C$reflex, frac_test[[2023]]$P$reflex )
# assuming the fraction of Ab testing in community is based on the % of people have been told the hiostory of HCv infeciton 
# the fraction of Ab testing in prison is based on the calibration the number of tests 
frac_ab[["2024"]] <- c(0.6, 0.2)

# calibrating the fm value 
fm <- list()
fm[["2022"]] <- c(1, 1, 9, 9, 1)
fm[["2023"]] <- c(1, 1, 16, 16, 1)
fm[["2024"]] <- c(1, 1, 16, 16, 1)
fm[["2025"]] <- c(0.9, 0.9, 6, 6, 1)
fm[["2026"]] <- c(0.9, 0.9, 6, 6, 1)
fm[["2027"]] <- c(0.9, 0.9, 5.6, 5.6, 1)

coverage_np <- list()
coverage_np[["2022"]] <- c(Ccal[[2022]]$C, Ccal[[2022]]$P)
coverage_np[["2023"]] <- c(Ccal[[2023]]$C, Ccal[[2023]]$P)


xfs <- list()
xfs[["2022"]] <- fs_estimate(num_ab = n_ab_np[["2022"]], 
                             cov_np = coverage_np[["2022"]], 
                             frac_ab = frac_ab[["2022"]], 
                             fp = fm[["2022"]], year = 2022,endY = 100,
                             modsim = Sce_sq)

fs <- list()
fs[["2022"]] <- xfs[["2022"]][[1]] 



fs[["2022"]][1,] <-   xfs[["2022"]][[1]][1, ]/fm[["2022"]][1]
fs[["2022"]][2,] <-   xfs[["2022"]][[1]][2, ]/fm[["2022"]][2]
fs[["2022"]][3,] <-   xfs[["2022"]][[1]][3, ]/fm[["2022"]][3]
fs[["2022"]][4,] <-   xfs[["2022"]][[1]][4, ]/fm[["2022"]][4]
fs[["2022"]][5,] <-   xfs[["2022"]][[1]][5, ]/xfs[["2022"]][[1]][5, ]

fs[["2022"]][5,91]

dfList_NP <- lapply(dfList, function(x) x*0)

param_var <- c("tau_ab","tau_RNA", "tau_poct", "eta") 
for(i in param_var){  
  dfList_NP[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP, index = i, 
                              frac_testing = frac_test[[2022]],
                              S_Yint = 2022, S_Yend = 2023, r_Yend = 2023, NPlst = NPlst, 
                              fp = c(Ccal[[2022]]$C*fm[["2022"]][1],
                                     Ccal[[2022]]$C*fm[["2022"]][2],
                                     Ccal[[2022]]$P*fm[["2022"]][3],
                                     Ccal[[2022]]$P*fm[["2022"]][4],
                                     Ccal[[2022]]$P*fm[["2022"]][5]))
  ini_dt <- (2022- POC_AU$cabY)/POC_AU$timestep + 1 
  end_dt <- ((2022 + 1 ) - POC_AU$cabY)/POC_AU$timestep
  dfList_NP[[i]][, , ini_dt:end_dt] <- dfList_NP[[i]][, ,end_dt]
}

#### #####

dfList_NP_2023 <- dfList_NP

for(i in param_var){  
  dfList_NP_2023[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP, index = i, 
                                   frac_testing = frac_test[[2023]],
                                   S_Yint = 2023, S_Yend = 2024, r_Yend = 2024, NPlst = NPlst, 
                                   fp = c(Ccal[[2023]]$C*fm[["2023"]][1], 
                                          Ccal[[2023]]$C*fm[["2023"]][2], 
                                          Ccal[[2023]]$P*fm[["2023"]][3], 
                                          Ccal[[2023]]$P*fm[["2023"]][4], 
                                          Ccal[[2023]]$P*fm[["2023"]][5]))
  ini_dt <- (2023 - POC_AU$cabY)/POC_AU$timestep + 1 
  end_dt <- ((2023 + 1 ) - POC_AU$cabY)/POC_AU$timestep
  dfList_NP_2023[[i]][, , ini_dt:end_dt] <- dfList_NP_2023[[i]][, ,end_dt]
}

for(i in param_var){
  # begining of 2024 
  b_pt <- (2024 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2023[[i]])[3]
  dfList_NP_2023[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

xfs[["2023"]] <- fs_estimate(num_ab = n_ab_np[["2023"]], 
                             cov_np = coverage_np[["2023"]], 
                             frac_ab = frac_ab[["2023"]], 
                             fp = fm[["2023"]], year = 2023, endY = 100,
                             modsim = Sce_sq)


fs[["2023"]] <- fs[["2022"]]
ini_dt <- (2023 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2023 + 1 ) - POC_AU$cabY)/POC_AU$timestep
fs[["2023"]][1, ini_dt:end_dt ] <-   xfs[["2023"]][[1]][1, ini_dt:end_dt ]/fm[["2023"]][1]
fs[["2023"]][2, ini_dt:end_dt ] <-   xfs[["2023"]][[1]][2, ini_dt:end_dt ]/fm[["2023"]][2]
fs[["2023"]][3, ini_dt:end_dt ] <-   xfs[["2023"]][[1]][3, ini_dt:end_dt ]/fm[["2023"]][3]
fs[["2023"]][4, ini_dt:end_dt ] <-   xfs[["2023"]][[1]][4, ini_dt:end_dt ]/fm[["2023"]][4]
fs[["2023"]][5, ini_dt:end_dt ] <-   xfs[["2023"]][[1]][5, ini_dt:end_dt ]/xfs[["2023"]][[1]][5, ini_dt:end_dt ] 

fm[["2022"]]
dfList_NP_2024 <- dfList_NP_2023 
Ccal[[2024]] <- lapply(Ccal[[2023]], function(x) x*2)

frac_test[[2024]] <- frac_test[[2023]]
frac_test[[2024]]$C$reflex <- frac_ab[["2024"]][1]
frac_test[[2024]]$C$immeRNA <- 1 - frac_ab[["2024"]][1]
frac_test[[2024]]$P$reflex <- frac_ab[["2024"]][2]
frac_test[[2024]]$P$immeRNA <- 1 - frac_ab[["2024"]][2]


for(i in param_var){  
  dfList_NP_2024[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2024, index = i, 
                                   frac_testing = frac_test[[2024]],
                                   S_Yint = 2024, S_Yend = 2025, r_Yend = 2025, NPlst = NPlst, 
                                   fp = c(Ccal[[2024]]$C*fm[["2024"]][1], 
                                                               Ccal[[2024]]$C*fm[["2024"]][2], 
                                                               Ccal[[2024]]$P*fm[["2024"]][3], 
                                                               Ccal[[2024]]$P*fm[["2024"]][4], 
                                                               Ccal[[2024]]$P*fm[["2024"]][5]))
  
}

for(i in param_var){
  # begining of 2024 
  b_pt <- (2025 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2024[[i]])[3]
  dfList_NP_2024[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2024"]] <- c(Ccal[[2024]]$C, Ccal[[2024]]$P)

xfs[["2024"]] <- fs_estimate(num_ab = n_ab_np[["2024"]], 
                             cov_np = coverage_np[["2024"]], 
                             frac_ab = frac_ab[["2024"]], 
                             fp = fm[["2024"]], year = 2024, endY = 100,
                             modsim = Sce_sq)


fs[["2024"]] <- fs[["2023"]]
ini_dt <- (2024 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2024 + 1 ) - POC_AU$cabY)/POC_AU$timestep
fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 


####NP expand_A ####
dfList_NPexp_A <- dfList_NP_2023
for(i in param_var){  
  dfList_NPexp_A[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2023, index = i, 
                                   frac_testing = frac_test[[2024]],
                                   S_Yint = 2024, S_Yend = 2025, r_Yend = 2028, NPlst = NPlst, 
                                   fp = c(Ccal[[2024]]$C*fm[["2024"]][1], 
                                          Ccal[[2024]]$C*fm[["2024"]][2], 
                                          Ccal[[2024]]$P*fm[["2024"]][3], 
                                          Ccal[[2024]]$P*fm[["2024"]][4], 
                                          Ccal[[2024]]$P*fm[["2024"]][5]))
  
}

for(i in param_var){
  # begining of 2024 
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2024[[i]])[3]
  dfList_NPexp_A[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

fs[["dfList_NPexp_A"]] <- fs[["2024"]]
ini_dt <- (2025 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2027 + 1 ) - POC_AU$cabY)/POC_AU$timestep
fs[["dfList_NPexp_A"]][, ini_dt:end_dt] <- fs[["dfList_NPexp_A"]][, ini_dt - 1]


####NP expand_B ####
odd_num_test <- 30000/20000
Ccal[[2025]] <- lapply(Ccal[[2024]],function(x) x*odd_num_test)

frac_test[[2025]] <- frac_test[[2024]]
frac_test[[2025]]$P$reflex <- 0.25
frac_test[[2025]]$P$immeRNA <- 0.75
frac_ab[["2025"]] <- c(unlist(as.numeric(frac_test[[2025]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2025]]$P$reflex)))
dfList_NPexp_B <- dfList_NP_2024
for(i in param_var){  
  dfList_NPexp_B[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2024, index = i, 
                                   frac_testing = frac_test[[2025]],
                                   S_Yint = 2025, S_Yend = 2026, r_Yend = 2028, NPlst = NPlst, 
                                   fp = c(Ccal[[2025]]$C*fm[["2025"]][1], 
                                          Ccal[[2025]]$C*fm[["2025"]][2], 
                                          Ccal[[2025]]$P*fm[["2025"]][3], 
                                          Ccal[[2025]]$P*fm[["2025"]][4], 
                                          Ccal[[2025]]$P*fm[["2025"]][5]))
  
}

for(i in param_var){
  # begining of 2024 
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2024[[i]])[3]
  dfList_NPexp_B[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

n_ab_np[["2025"]] <- n_ab_np[["2024"]]*odd_num_test
coverage_np[["2025"]] <- coverage_np[["2024"]]*odd_num_test

xfs[["2025"]] <- fs_estimate(num_ab = n_ab_np[["2025"]], 
                             cov_np = coverage_np[["2025"]], 
                             frac_ab = frac_ab[["2025"]], 
                             fp = fm[["2025"]], year = 2025, endY = 100,
                             modsim = Sce_sq)

fs[["dfList_NPexp_B"]] <- fs[["2024"]]
ini_dt <- (2025 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2025 + 1 ) - POC_AU$cabY)/POC_AU$timestep
fs[["dfList_NPexp_B"]][1, ini_dt:end_dt ] <-   xfs[["2025"]][[1]][1, ini_dt:end_dt ]/fm[["2025"]][1]
fs[["dfList_NPexp_B"]][2, ini_dt:end_dt ] <-   xfs[["2025"]][[1]][2, ini_dt:end_dt ]/fm[["2025"]][2]
fs[["dfList_NPexp_B"]][3, ini_dt:end_dt ] <-   xfs[["2025"]][[1]][3, ini_dt:end_dt ]/fm[["2025"]][3]
fs[["dfList_NPexp_B"]][4, ini_dt:end_dt ] <-   xfs[["2025"]][[1]][4, ini_dt:end_dt ]/fm[["2025"]][4]
fs[["dfList_NPexp_B"]][5, ini_dt:end_dt ] <-   xfs[["2025"]][[1]][5, ini_dt:end_dt ]/xfs[["2025"]][[1]][5, ini_dt:end_dt ] 

ini_dt <- (2026 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2027 + 1 ) - POC_AU$cabY)/POC_AU$timestep
fs[["dfList_NPexp_B"]][, ini_dt:end_dt] <- fs[["dfList_NPexp_B"]][, ini_dt - 1]


####NP expand_C ####
odd_num_test <- 40000/30000
Ccal[[2026]] <- lapply(Ccal[[2025]],function(x) x*odd_num_test)


dfList_NPexp_C <- dfList_NPexp_B
for(i in param_var){  
  dfList_NPexp_C[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NPexp_B, index = i, 
                                   frac_testing = frac_test[[2025]],
                                   S_Yint = 2026, S_Yend = 2027, r_Yend = 2028, NPlst = NPlst, 
                                   fp = c(Ccal[[2026]]$C*fm[["2026"]][1], 
                                          Ccal[[2026]]$C*fm[["2026"]][2], 
                                          Ccal[[2026]]$P*fm[["2026"]][3], 
                                          Ccal[[2026]]$P*fm[["2026"]][4], 
                                          Ccal[[2026]]$P*fm[["2026"]][5]))
  
}

for(i in param_var){
  # begining of 2024 
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2024[[i]])[3]
  dfList_NPexp_C[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

n_ab_np[["2026"]] <- n_ab_np[["2025"]]*odd_num_test
coverage_np[["2026"]] <- coverage_np[["2025"]]*odd_num_test

xfs[["2026"]] <- fs_estimate(num_ab = n_ab_np[["2026"]], 
                             cov_np = coverage_np[["2026"]], 
                             frac_ab = frac_ab[["2025"]], 
                             fp = fm[["2026"]], year = 2026, endY = 100,
                             modsim = Sce_sq)

fs[["dfList_NPexp_C"]] <- fs[["dfList_NPexp_B"]]
ini_dt <- (2026 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2026 + 1 ) - POC_AU$cabY)/POC_AU$timestep
fs[["dfList_NPexp_C"]][1, ini_dt:end_dt ] <-   xfs[["2026"]][[1]][1, ini_dt:end_dt ]/fm[["2026"]][1]
fs[["dfList_NPexp_C"]][2, ini_dt:end_dt ] <-   xfs[["2026"]][[1]][2, ini_dt:end_dt ]/fm[["2026"]][2]
fs[["dfList_NPexp_C"]][3, ini_dt:end_dt ] <-   xfs[["2026"]][[1]][3, ini_dt:end_dt ]/fm[["2026"]][3]
fs[["dfList_NPexp_C"]][4, ini_dt:end_dt ] <-   xfs[["2026"]][[1]][4, ini_dt:end_dt ]/fm[["2026"]][4]
fs[["dfList_NPexp_C"]][5, ini_dt:end_dt ] <-   xfs[["2026"]][[1]][5, ini_dt:end_dt ]/xfs[["2026"]][[1]][5, ini_dt:end_dt ] 

ini_dt <- (2027 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2027 + 1 ) - POC_AU$cabY)/POC_AU$timestep
fs[["dfList_NPexp_C"]][, ini_dt:end_dt] <- fs[["dfList_NPexp_C"]][, ini_dt - 1]


####NP expand_D ####
odd_num_test <- 50000/40000
Ccal[[2027]] <- lapply(Ccal[[2026]],function(x) x)

frac_test[[2027]] <- frac_test[[2025]]
frac_test[[2027]]$P$reflex <- 0.5
frac_test[[2027]]$P$immeRNA <- 0.5
dfList_NPexp_D <- dfList_NPexp_C
for(i in param_var){  
  dfList_NPexp_D[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NPexp_C, index = i, 
                                   frac_testing = frac_test[[2025]],
                                   S_Yint = 2027, S_Yend = 2028, r_Yend = 2028, NPlst = NPlst, 
                                   fp = c(Ccal[[2027]]$C*fm[["2027"]][1], 
                                          Ccal[[2027]]$C*fm[["2027"]][2], 
                                          Ccal[[2027]]$P*fm[["2027"]][3], 
                                          Ccal[[2027]]$P*fm[["2027"]][4], 
                                          Ccal[[2027]]$P*fm[["2027"]][5]))
  
}

for(i in param_var){
  # begining of 2024 
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2024[[i]])[3]
  dfList_NPexp_D[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

n_ab_np[["2027"]] <- n_ab_np[["2026"]]
coverage_np[["2027"]] <- coverage_np[["2026"]]

xfs[["2027"]] <- fs_estimate(num_ab = n_ab_np[["2027"]], 
                             cov_np = coverage_np[["2027"]], 
                             frac_ab = frac_ab[["2025"]], 
                             fp = fm[["2027"]], year = 2027, endY = 100,
                             modsim = Sce_sq)

fs[["dfList_NPexp_D"]] <- fs[["dfList_NPexp_C"]]
ini_dt <- (2027 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2027 + 1 ) - POC_AU$cabY)/POC_AU$timestep
fs[["dfList_NPexp_D"]][1, ini_dt:end_dt ] <-   xfs[["2027"]][[1]][1, ini_dt:end_dt ]/fm[["2027"]][1]
fs[["dfList_NPexp_D"]][2, ini_dt:end_dt ] <-   xfs[["2027"]][[1]][2, ini_dt:end_dt ]/fm[["2027"]][2]
fs[["dfList_NPexp_D"]][3, ini_dt:end_dt ] <-   xfs[["2027"]][[1]][3, ini_dt:end_dt ]/fm[["2027"]][3]
fs[["dfList_NPexp_D"]][4, ini_dt:end_dt ] <-   xfs[["2027"]][[1]][4, ini_dt:end_dt ]/fm[["2027"]][4]
fs[["dfList_NPexp_D"]][5, ini_dt:end_dt ] <-   xfs[["2027"]][[1]][5, ini_dt:end_dt ]/xfs[["2027"]][[1]][5, ini_dt:end_dt ] 


endY <- 20
xxx <- HCVMSM(POC_AU, best_estimates, best_est_pop,
              disease_progress,pop_array,
              dfList, 
              param_cascade_sc = dfList_NPexp_D, fib = fib, 
              modelrun="UN", proj = "POC_AU", end_Y = endY,  
              fc_sc= fs[["dfList_NPexp_D"]],
              fp = c(1,1,1,1,1))



scenario_cascade <- list("dfList_NP_2023" = dfList_NP_2023, 
                         "dfList_NP_2024" = dfList_NP_2024, 
                         "dfList_NPexp_A" = dfList_NPexp_A, 
                         "dfList_NPexp_B" = dfList_NPexp_B, 
                         "dfList_NPexp_C" = dfList_NPexp_C, 
                         "dfList_NPexp_D" = dfList_NPexp_D)

scenario_fc <- list("dfList_NP_2023" = fs[["2023"]], 
                    "dfList_NP_2024" = fs[["2024"]], 
                    "dfList_NPexp_A" = fs$dfList_NPexp_A, 
                    "dfList_NPexp_B" = fs$dfList_NPexp_B, 
                    "dfList_NPexp_C" = fs$dfList_NPexp_C, 
                    "dfList_NPexp_D" = fs$dfList_NPexp_D)



save(scenario_cascade,
     scenario_fc,
     file = file.path(OutputFolder,
                      paste0(project_name, "scenario_cascade", ".rda")))

tic <- proc.time()
endY <- 100
Sce_np <- list()

for(scenario in names(scenario_cascade)){ 
  Sce_np[[scenario]] <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                               disease_progress,pop_array,
                               dfList,  
                               param_cascade_sc = scenario_cascade[[scenario]] , fib = fib, 
                               modelrun="UN", proj = "POC_AU", end_Y = endY, 
                               cost = costdfList, costflow = costflow, 
                               costflow_Neg = costflow_Neg, 
                               fc = scenario_fc[[scenario]])
  
}

toc <- proc.time() - tic
toc
save(Sce_sq,Sce_np,
     file = file.path(OutputFolder ,
                      paste0(project_name,"Simulations" ,".rda"))) 


test <- list()
cl_ext <- names(Sce_np$dfList_NPexp_D)[c(10:22)]
for(i in cl_ext){
  
  test[[i]] <- modres.flow.t(POC_AU, Sce_np$dfList_NPexp_D, endYear = 100, 
                             allp = i)%>%
    ungroup()%>%
    group_by(year, population)%>%
    summarise(best = sum(best))
  
  
}
test_sq <- list()
for(i in cl_ext){
  
  test_sq[[i]] <- modres.flow.t(POC_AU, Sce_sq, endYear = 100, 
                             allp = i)%>%
    ungroup()%>%
    group_by(year, population)%>%
    summarise(best = sum(best))
  
  
}

test$newTestingPOCT_sc
test <- dplyr::bind_rows(test, .id = 'index')%>%group_by(year, population)%>%spread(index, best)
test_sq <- dplyr::bind_rows(test_sq, .id = 'index')%>%group_by(year, population)%>%spread(index, best)


test_fscal <- test%>%mutate(Ab = newTestingAb_sc + newTestingAb_sc_neg, 
                            RNA = (newTestingAg_sc+ newTestingAg_sc_neg + newTestingPOCT_sc + 
                                     newTestingPOCT_sc_neg ))%>%select(year, population, Ab, 
                                                                  RNA)%>%
  mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), "C", "P"))%>%
  ungroup()%>%
  select(-c(population))%>%
  gather(index, value, -c(year, setting))%>%
  group_by(year, setting, index)%>%summarise(value = sum(value))%>%
  mutate(scenario = "exp_D")
dtp <- c(6615, 11276, 19480, 30000, 40000, 50000)
View(test_fscal%>%filter(year%in% c(7,8,9,10,11,12))
              )
View(test_fscal%>%filter(year%in% c(7,8,9,10,11,12))%>%group_by(year)%>%
       summarise(tot_NP_test = sum(value))%>%
       mutate(year = year + 2015, 
              data_point = unlist(dtp))%>%
       mutate(diff = data_point- tot_NP_test))
test%>%mutate(setting = ifelse(population %in% POC_AU$popNames[3:5], "P", "C"))%>%
  group_by(setting, year)%>%summarise(best = sum(newTreatment_sc))%>%
  filter(year%in% c(7:12))

test_fscal

test_fscal_sq <- test_sq%>%mutate(Ab = newTestingAb_sc + newTestingAb_sc_neg + newTestingAb + newTestingAb_neg, 
                            RNA = (newTestingAg_sc+ newTestingAg_sc_neg + newTestingPOCT_sc + 
                                     newTestingPOCT_sc_neg + 
                                     newTestingAg+ newTestingAg_neg + newTestingPOCT + 
                                     newTestingPOCT_neg))%>%select(year, population, Ab, 
                                                                   RNA)%>%
  mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), "C", "P"))%>%
  ungroup()%>%
  select(-c(population))%>%
  gather(index, value, -c(year, setting))%>%
  group_by(year, setting, index)%>%summarise(value = sum(value))%>%
  mutate(scenario = "sq")


View(test_fscal%>%group_by(year)%>%summarise(best = sum(value)))

x <- rbind(test_fscal_sq, test_fscal)

ggplot(data = x%>%filter(index == "Ab"), aes(x = year, y = value)) + 
  geom_line(aes(colour = scenario)) + 
  facet_wrap(~setting) 


#################################################################################
Sce_np$dfList_NP_2023$newTreatment

calibrateFlow_sc <- Sce_np$dfList_NP_2023[names(Sce_np$dfList_NP_2023)%in%
                                 c("newTreatment", "newRetreat", 
                                   "newTreatment_sc")]
flow_sub <- list()


flow_sub <- lapply(names(calibrateFlow_sc), function(x){ 
  a <- indicatorResults(POC_AU, calibrateFlow_sc, x, 
                        pop=POC_AU$popNames,
                        paramR = NULL, range = NULL,
                        endY = endY)
})

names(flow_sub) <- names(calibrateFlow_sc)


flow_setting <- lapply(flow_sub, function(x){ 
  
  a <- x%>%
    mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), 
                            "commu", "prisons"))
  
  a <- a%>%group_by(year, setting)%>%summarise_at("best", sum)%>%
    mutate(population = setting)%>%select(-setting)
}) 

x <- cbind(flow_setting$newTreatment[, c(1:3)], flow_setting$newRetreat[, 2], 
      flow_setting$newTreatment_sc[,2])
View(x)
N_treatment <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                each = 2),
                     population = flow_setting$newTreatment[ ,3],
                     as.data.frame(flow_setting$newTreatment[, -c(1,3)] + 
                                     flow_setting$newRetreat[, -c(1,3)] + 
                                     flow_setting$newTreatment_sc[, -c(1,3)]))%>%
  tibble::as_tibble()%>%
  mutate(population = factor(population, 
                             levels = c("commu", "prisons"), 
                             labels = c("Community", "Prisons" )))



HCVtreatinitN_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatinitN_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = realpop,
                           up = upper,
                           low = lower,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))

endY_plot <- 2030- POC_AU$cabY

N_treatment_setting_p <- indicatorPlot(POC_AU, N_treatment, 
                                       ylabel = "N",
                                       xlimits = c(POC_AU$startYear, 
                                                   (POC_AU$startYear+endY_plot), 5),
                                       calibration_Y = POC_AU$cabY,
                                       rangeun = NULL, 
                                       groupPlot = NULL, 
                                       facetPlot = population,
                                       observationData = HCVtreatinitN_setting_fit, 
                                       simulateYear = POC_AU$simY) + 
  ggtitle("Number of treatment init") + theme_bw()

N_treatment_setting_p <- N_treatment_setting_p + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 40000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 5000)))))
N_treatment_setting_p 
N_treatment%>%filter(year %in% c(8,9))
ggplot(data = N_treatment%>%filter(year<=15 & population == "Prisons"), 
       aes(x = year, y = best)) + 
  geom_line() + geom_point(aes(x = 7, y = 2679))

unlist(c(N_treatment%>%group_by(year)%>%summarise(best= sum(best))%>%mutate(year = year + POC_AU$cabY -1)%>%
  filter(year%in% seq(2015,2023,1))%>%select(best)))

N_treatment%>%filter(year%in% c(8,9))

HCVtreatinitN_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatinitN_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = as.numeric(realpop),
                           up = as.numeric(upper),
                           low = as.numeric(lower))


N_treatment%>%group_by(year)%>%summarise(best = sum(best))%>%indicatorPlot(POC_AU,. , 
                                       ylabel = "N",
                                       xlimits = c(POC_AU$startYear, 
                                                   (POC_AU$startYear+endY_plot), 5),
                                       calibration_Y = POC_AU$cabY,
                                       rangeun = NULL, 
                                       groupPlot = NULL, 
                                       facetPlot = NULL,
                                       observationData = HCVtreatinitN_fit, 
                                       simulateYear = POC_AU$simY) + 
  ggtitle("Number of treatment init") + theme_bw()
