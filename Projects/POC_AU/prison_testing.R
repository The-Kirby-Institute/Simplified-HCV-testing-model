rm(list = ls())
library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
library(readxl)
project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Documents/Simplified-HCV-testing-model-fresh")

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

if (!dir.exists(OutputFolder)) {
  dir.create(OutputFolder, recursive = TRUE, showWarnings = FALSE)
}

load(file.path(RdaFolder, paste0(project_name, ".rda")))
load(file.path(RdaFolder, paste0(project_name, "cali.rda")))
load(file.path(RdaFolder, paste0(project_name, "cali_timev.rda")))

source(file.path(Rcode, "Functions/HCV_model.R"))

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

#### sensitivity total cost(including program cost) ####
###################### debug required ##########################################
files <- list.files(path = paste0(DataFolder, 
                                  "/cost", sep =  ""), pattern = '*.csv')

# parameter sets for cost data 
# +- 10% 
costdfList <- list()
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

###############################################################################

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

Num_test_person_NP <- read_excel(paste0(data_path, "/01. DATA/Num_test_person_NP.xlsx"))%>%
  select(year, settings, num_ab, num_RNA, num_tests, num_person)%>%
  mutate(num_RNA = as.numeric(num_RNA), 
         num_tests = as.numeric(num_tests),
         num_person = as.numeric(num_person))%>%
  filter(year%in% c(2022:2024))
# assuming actively C_PWID and C_fPWID roughly equals to the pop size of CPWID
# We accounted the transition in prison setting regarding its high dynamic nature 
Num_test_person_NP <- Num_test_person_NP%>%na.omit()%>%mutate(num_pop = c(80000, 80000, 80000, 80000, 80000,80000))

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

for(i in c(2022:2024)){ 
  
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
Ccal[[2023]]

# for year of 2022, using RNA testing in prisons for n_ab to estimate the testing coverage 
n_ab <- Num_test_person_NP%>%select(year, settings, num_ab)%>%
  spread(settings, -c(year))

n_ab[1, "prison"] <- Num_test_person_NP%>%
  filter(year == 2022 & settings == "prison")%>%select(num_RNA)%>%unlist()%>%c()

# efficacy of national program 

np_effect <- read_excel(paste0(data_path, "/01. DATA/effect_np.xlsx"))

NP_tauab_C <- unlist(as.numeric(np_effect[1,2]))

NP_tauRNA_C <- unlist(as.numeric(np_effect[2,2]))

NP_tauRNAonly_C <- unlist(as.numeric(np_effect[3,2]))

# the treatment initiation is the % of people tested RNA+ initiated DAA within 120 days 
# old 
# NP_eta_C <-  1- (1- unlist(as.numeric(np_effect[4,2])))^(1/POC_AU$timestep/4)

# test on new data informed by national program [2025/07/30]
# NP_eta_C <- 0.6
NP_eta_C <- 1- (1- unlist(as.numeric(np_effect[4,2])))

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
n_ab_np[["2022"]] <- c(unlist(as.numeric(n_ab[n_ab$year == 2022,"community"])), 
                       unlist(as.numeric(n_ab[n_ab$year == 2022,"prison"])))
n_ab_np[["2023"]] <- c(unlist(as.numeric(n_ab[n_ab$year == 2023,"community"])), 
                       unlist(as.numeric(n_ab[n_ab$year == 2023,"prison"])))
n_ab_np[["2024"]] <- c(unlist(as.numeric(n_ab[n_ab$year == 2024,"community"])), 
                       unlist(as.numeric(n_ab[n_ab$year == 2024,"prison"])))

frac_ab <- list()
frac_ab[["2022"]] <- c(frac_test[[2022]]$C$reflex, frac_test[[2022]]$P$immeRNA )
frac_ab[["2023"]] <- c(frac_test[[2023]]$C$reflex, frac_test[[2023]]$P$reflex )
frac_ab[["2024"]] <- c(frac_test[[2024]]$C$reflex, frac_test[[2024]]$P$reflex )
frac_ab[["2025"]] <- c(frac_test[[2024]]$C$reflex, frac_test[[2024]]$P$reflex )
# assuming the fraction of Ab testing in community is based on the % of people have been told the history of HCv infection 
# the fraction of Ab testing in prison is based on the calibration the number of tests 

# prisons testing scenarios 
fm <- list()
fm[["2022"]] <- c(1.1,1.1,  2.4,  2.4 , 0.01)
fm[["2023"]] <- c(1.1, 1.1, 18.5, 18.5, 1)
fm[["2024"]] <- c(1.1, 1.1, 18.5, 18.5, 1)
# fm[["2025"]] <- c(0.1, 0.1, 15, 15, 1)
fm[["2025"]] <- c(1.1, 1.1, 15, 15, 1)
fm[["2026"]] <- c(1.1, 1.1, 15, 15, 1)
#fm[["2026"]] <- c(0.1, 0.1, 15, 15, 1)
# fm[["2027"]] <- c(0.1, 0.1, 13, 13, 1)
fm[["2027"]] <- c(0.8, 0.8, 13, 13, 1)
fm[["2028"]] <- c(0.8, 0.8, 15, 15, 1)
fm[["2029"]] <- c(0.8, 0.8, 13, 13, 1)
fm[["2030"]] <- c(0.8, 0.8, 15, 15, 1)
coverage_np <- list()
coverage_np[["2022"]] <- c(Ccal[[2022]]$C, Ccal[[2022]]$P)
coverage_np[["2023"]] <- c(Ccal[[2023]]$C, Ccal[[2023]]$P)
coverage_np[["2024"]] <- c(Ccal[[2024]]$C, Ccal[[2024]]$P)

xfs <- list()
xfs[["2022"]] <- fs_estimate(num_ab = n_ab_np[["2022"]], 
                             cov_np = coverage_np[["2022"]], 
                             frac_ab = frac_ab[["2022"]], 
                             fp = fm[["2022"]], year = 2022,endY = 100,
                             modsim = Sce_sq)

xfs[["2023"]] <- fs_estimate(num_ab = n_ab_np[["2023"]], 
                             cov_np = coverage_np[["2023"]], 
                             frac_ab = frac_ab[["2023"]], 
                             fp = fm[["2023"]], year = 2023,endY = 100,
                             modsim = Sce_sq)

xfs[["2024"]] <- fs_estimate(num_ab = n_ab_np[["2024"]], 
                             cov_np = coverage_np[["2024"]], 
                             frac_ab = frac_ab[["2024"]], 
                             fp = fm[["2024"]], year = 2024,endY = 100,
                             modsim = Sce_sq)

fs <- list()
fs[["2022"]] <- xfs[["2022"]][[1]] 


dfList_NP <- lapply(dfList, function(x) x*0)

param_var <- c("tau_ab","tau_RNA", "tau_poct", "eta") 
ini_dt <- (2022 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2022 + 1 ) - POC_AU$cabY)/POC_AU$timestep
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



ini_dt <- (2023 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2023 + 1 ) - POC_AU$cabY)/POC_AU$timestep
fs[["2023"]] <- fs[["2022"]]
fs[["2023"]][1, ini_dt: end_dt] <-   xfs[["2023"]][[1]][1, ini_dt: end_dt]/fm[["2023"]][1]
fs[["2023"]][2, ini_dt: end_dt] <-   xfs[["2023"]][[1]][2, ini_dt: end_dt]/fm[["2023"]][2]
fs[["2023"]][3, ini_dt: end_dt] <-   xfs[["2023"]][[1]][3, ini_dt: end_dt]/fm[["2023"]][3]
fs[["2023"]][4, ini_dt: end_dt] <-   xfs[["2023"]][[1]][4, ini_dt: end_dt]/fm[["2023"]][4]
fs[["2023"]][5, ini_dt: end_dt] <-   xfs[["2023"]][[1]][5, ini_dt: end_dt]/xfs[["2023"]][[1]][5, ini_dt: end_dt]
# dfList_NP_2024 <- dfList_NP_2023 
# Ccal[[2024]] <- lapply(Ccal[[2023]], function(x) x*2)

# frac_test[[2024]] <- frac_test[[2023]]
# frac_test[[2024]]$C$reflex <- frac_ab[["2024"]][1]
# frac_test[[2024]]$C$immeRNA <- 1 - frac_ab[["2024"]][1]
# frac_test[[2024]]$P$reflex <- frac_ab[["2024"]][2]
# frac_test[[2024]]$P$immeRNA <- 1 - frac_ab[["2024"]][2]

dfList_NP_2024 <- dfList_NP_2023

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


# fs[["2024"]] <- fs[["2023"]]
ini_dt <- (2024 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2024 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2024"]] <- fs[["2023"]]
fs[["2024"]][1, ini_dt: end_dt] <-   xfs[["2024"]][[1]][1, ini_dt: end_dt]/fm[["2024"]][1]
fs[["2024"]][2, ini_dt: end_dt] <-   xfs[["2024"]][[1]][2, ini_dt: end_dt]/fm[["2024"]][2]
fs[["2024"]][3, ini_dt: end_dt] <-   xfs[["2024"]][[1]][3, ini_dt: end_dt]/fm[["2024"]][3]
fs[["2024"]][4, ini_dt: end_dt] <-   xfs[["2024"]][[1]][4, ini_dt: end_dt]/fm[["2024"]][4]
fs[["2024"]][5, ini_dt: end_dt] <-   xfs[["2024"]][[1]][5, ini_dt: end_dt]/xfs[["2024"]][[1]][5, ini_dt: end_dt]

# 2025

odd_num_test <- 1

Ccal[[2025]] <- lapply(Ccal[[2024]],function(x) x*odd_num_test)

frac_test[[2025]] <- frac_test[[2024]]

frac_ab[["2025"]] <- c(unlist(as.numeric(frac_test[[2025]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2025]]$P$reflex)))

dfList_NP_2025 <- dfList_NP_2024
for(i in param_var){  
  dfList_NP_2025[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2025, index = i, 
                                   frac_testing = frac_test[[2025]],
                                   S_Yint = 2025, S_Yend = 2026, r_Yend = 2026, NPlst = NPlst, 
                                   fp = c(Ccal[[2025]]$C*fm[["2025"]][1], 
                                          Ccal[[2025]]$C*fm[["2025"]][2], 
                                          Ccal[[2025]]$P*fm[["2025"]][3], 
                                          Ccal[[2025]]$P*fm[["2025"]][4], 
                                          Ccal[[2025]]$P*fm[["2025"]][5]))
  
}

for(i in param_var){
  # begining of 2025
  b_pt <- (2026 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2025[[i]])[3]
  dfList_NP_2025[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2025"]] <- c(Ccal[[2025]]$C, Ccal[[2025]]$P)
n_ab_np[["2025"]] <- n_ab_np[["2024"]]
xfs[["2025"]] <- fs_estimate(num_ab = n_ab_np[["2025"]], 
                             cov_np = coverage_np[["2025"]], 
                             frac_ab = frac_ab[["2025"]], 
                             fp = fm[["2025"]], year = 2025, endY = 100,
                             modsim = Sce_sq)


# fs[["2024"]] <- fs[["2023"]]
ini_dt <- (2025 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2025 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2025"]] <- fs[["2024"]]
fs[["2025"]][1, ini_dt: end_dt] <-   xfs[["2025"]][[1]][1, ini_dt: end_dt]/fm[["2025"]][1]
fs[["2025"]][2, ini_dt: end_dt] <-   xfs[["2025"]][[1]][2, ini_dt: end_dt]/fm[["2025"]][2]
fs[["2025"]][3, ini_dt: end_dt] <-   xfs[["2025"]][[1]][3, ini_dt: end_dt]/fm[["2025"]][3]
fs[["2025"]][4, ini_dt: end_dt] <-   xfs[["2025"]][[1]][4, ini_dt: end_dt]/fm[["2025"]][4]
fs[["2025"]][5, ini_dt: end_dt] <-   xfs[["2025"]][[1]][5, ini_dt: end_dt]/xfs[["2025"]][[1]][5, ini_dt: end_dt]

# 2026

odd_num_test <- 1.08

Ccal[[2026]] <- lapply(Ccal[[2025]],function(x) x*odd_num_test)

frac_test[[2026]] <- frac_test[[2025]]

frac_ab[["2026"]] <- c(unlist(as.numeric(frac_test[[2026]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2026]]$P$reflex)))

dfList_NPPhaseII <- dfList_NP_2025
for(i in param_var){  
  dfList_NPPhaseII[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NPPhaseII, index = i, 
                                     frac_testing = frac_test[[2026]],
                                     S_Yint = 2026, S_Yend = 2027, r_Yend = 2027, NPlst = NPlst, 
                                     fp = c(Ccal[[2026]]$C*fm[["2026"]][1], 
                                            Ccal[[2026]]$C*fm[["2026"]][2], 
                                            Ccal[[2026]]$P*fm[["2026"]][3], 
                                            Ccal[[2026]]$P*fm[["2026"]][4], 
                                            Ccal[[2026]]$P*fm[["2026"]][5]))
  
}

for(i in param_var){
  # begining of 2026
  b_pt <- (2027 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPPhaseII[[i]])[3]
  dfList_NPPhaseII[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2026"]] <- c(Ccal[[2026]]$C, Ccal[[2026]]$P)
n_ab_np[["2026"]] <- n_ab_np[["2025"]]*odd_num_test
xfs[["2026"]] <- fs_estimate(num_ab = n_ab_np[["2026"]], 
                             cov_np = coverage_np[["2026"]], 
                             frac_ab = frac_ab[["2026"]], 
                             fp = fm[["2026"]], year = 2026, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2026 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2026 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2026"]] <- fs[["2025"]]
fs[["2026"]][1, ini_dt: end_dt] <-   xfs[["2026"]][[1]][1, ini_dt: end_dt]/fm[["2026"]][1]
fs[["2026"]][2, ini_dt: end_dt] <-   xfs[["2026"]][[1]][2, ini_dt: end_dt]/fm[["2026"]][2]
fs[["2026"]][3, ini_dt: end_dt] <-   xfs[["2026"]][[1]][3, ini_dt: end_dt]/fm[["2026"]][3]
fs[["2026"]][4, ini_dt: end_dt] <-   xfs[["2026"]][[1]][4, ini_dt: end_dt]/fm[["2026"]][4]
fs[["2026"]][5, ini_dt: end_dt] <-   xfs[["2026"]][[1]][5, ini_dt: end_dt]/xfs[["2026"]][[1]][5, ini_dt: end_dt]

# 2027

odd_num_test <- 1

Ccal[[2027]] <- lapply(Ccal[[2026]],function(x) x*odd_num_test)

frac_test[[2027]] <- frac_test[[2026]]

frac_ab[["2027"]] <- c(unlist(as.numeric(frac_test[[2027]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2027]]$P$reflex)))

dfList_NP_2027 <- dfList_NPPhaseII
for(i in param_var){  
  dfList_NP_2027[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2027, index = i, 
                                   frac_testing = frac_test[[2027]],
                                   S_Yint = 2027, S_Yend = 2028, r_Yend = 2028, NPlst = NPlst, 
                                   fp = c(Ccal[[2027]]$C*fm[["2027"]][1], 
                                          Ccal[[2027]]$C*fm[["2027"]][2], 
                                          Ccal[[2027]]$P*fm[["2027"]][3], 
                                          Ccal[[2027]]$P*fm[["2027"]][4], 
                                          Ccal[[2027]]$P*fm[["2027"]][5]))
  
}

for(i in param_var){
  # begining of 2027
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPPhaseII[[i]])[3]
  dfList_NP_2027[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2027"]] <- c(Ccal[[2027]]$C, Ccal[[2027]]$P)
n_ab_np[["2027"]] <- n_ab_np[["2026"]]*odd_num_test
xfs[["2027"]] <- fs_estimate(num_ab = n_ab_np[["2027"]], 
                             cov_np = coverage_np[["2027"]], 
                             frac_ab = frac_ab[["2027"]], 
                             fp = fm[["2027"]], year = 2027, endY = 100,
                             modsim = Sce_sq)

View(xfs[["2027"]][[1]])
ini_dt <- (2027 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2027 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2027"]] <- fs[["2026"]]
fs[["2027"]][1, ini_dt: end_dt] <-   xfs[["2027"]][[1]][1, ini_dt: end_dt]/fm[["2027"]][1]
fs[["2027"]][2, ini_dt: end_dt] <-   xfs[["2027"]][[1]][2, ini_dt: end_dt]/fm[["2027"]][2]
fs[["2027"]][3, ini_dt: end_dt] <-   xfs[["2027"]][[1]][3, ini_dt: end_dt]/fm[["2027"]][3]
fs[["2027"]][4, ini_dt: end_dt] <-   xfs[["2027"]][[1]][4, ini_dt: end_dt]/fm[["2027"]][4]
fs[["2027"]][5, ini_dt: end_dt] <-   xfs[["2027"]][[1]][5, ini_dt: end_dt]/xfs[["2027"]][[1]][5, ini_dt: end_dt]




# 2028 

odd_num_test <- 1.007

Ccal[[2028]] <- lapply(Ccal[[2027]],function(x) x*odd_num_test)

frac_test[[2028]] <- frac_test[[2026]]

frac_ab[["2028"]] <- c(unlist(as.numeric(frac_test[[2028]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2028]]$P$reflex)))

dfList_NP_2028 <- dfList_NP_2027
for(i in param_var){  
  dfList_NP_2028[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2028, index = i, 
                                   frac_testing = frac_test[[2028]],
                                   S_Yint = 2028, S_Yend = 2029, r_Yend = 2029, NPlst = NPlst, 
                                   fp = c(Ccal[[2028]]$C*fm[["2028"]][1], 
                                          Ccal[[2028]]$C*fm[["2028"]][2], 
                                          Ccal[[2028]]$P*fm[["2028"]][3], 
                                          Ccal[[2028]]$P*fm[["2028"]][4], 
                                          Ccal[[2028]]$P*fm[["2028"]][5]))
  
}

for(i in param_var){
  # begining of 2026
  b_pt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2028[[i]])[3]
  dfList_NP_2028[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2028"]] <- c(Ccal[[2028]]$C, Ccal[[2028]]$P)
n_ab_np[["2028"]] <- n_ab_np[["2027"]]*odd_num_test
xfs[["2028"]] <- fs_estimate(num_ab = n_ab_np[["2028"]], 
                             cov_np = coverage_np[["2028"]], 
                             frac_ab = frac_ab[["2028"]], 
                             fp = fm[["2028"]], year = 2028, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2028 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2028"]] <- fs[["2027"]]
fs[["2028"]][1, ini_dt: end_dt] <-   xfs[["2028"]][[1]][1, ini_dt: end_dt]/fm[["2028"]][1]
fs[["2028"]][2, ini_dt: end_dt] <-   xfs[["2028"]][[1]][2, ini_dt: end_dt]/fm[["2028"]][2]
fs[["2028"]][3, ini_dt: end_dt] <-   xfs[["2028"]][[1]][3, ini_dt: end_dt]/fm[["2028"]][3]
fs[["2028"]][4, ini_dt: end_dt] <-   xfs[["2028"]][[1]][4, ini_dt: end_dt]/fm[["2028"]][4]
fs[["2028"]][5, ini_dt: end_dt] <-   xfs[["2028"]][[1]][5, ini_dt: end_dt]/xfs[["2028"]][[1]][5, ini_dt: end_dt]

# 2029

odd_num_test <- 1

Ccal[[2029]] <- lapply(Ccal[[2028]],function(x) x*odd_num_test)

frac_test[[2029]] <- frac_test[[2028]]

frac_ab[["2029"]] <- c(unlist(as.numeric(frac_test[[2029]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2029]]$P$reflex)))

dfList_NP_2029 <- dfList_NP_2028
for(i in param_var){  
  dfList_NP_2029[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2029, index = i, 
                                   frac_testing = frac_test[[2029]],
                                   S_Yint = 2029, S_Yend = 2030, r_Yend = 2030, NPlst = NPlst, 
                                   fp = c(Ccal[[2029]]$C*fm[["2029"]][1], 
                                          Ccal[[2029]]$C*fm[["2029"]][2], 
                                          Ccal[[2029]]$P*fm[["2029"]][3], 
                                          Ccal[[2029]]$P*fm[["2029"]][4], 
                                          Ccal[[2029]]$P*fm[["2029"]][5]))
  
}

for(i in param_var){
  # begining of 2029
  b_pt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2029[[i]])[3]
  dfList_NP_2029[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2029"]] <- c(Ccal[[2029]]$C, Ccal[[2029]]$P)
n_ab_np[["2029"]] <- n_ab_np[["2028"]]*odd_num_test
xfs[["2029"]] <- fs_estimate(num_ab = n_ab_np[["2029"]], 
                             cov_np = coverage_np[["2029"]], 
                             frac_ab = frac_ab[["2029"]], 
                             fp = fm[["2029"]], year = 2029, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2029 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2029"]] <- fs[["2028"]]
fs[["2029"]][1, ini_dt: end_dt] <-   xfs[["2029"]][[1]][1, ini_dt: end_dt]/fm[["2029"]][1]
fs[["2029"]][2, ini_dt: end_dt] <-   xfs[["2029"]][[1]][2, ini_dt: end_dt]/fm[["2029"]][2]
fs[["2029"]][3, ini_dt: end_dt] <-   xfs[["2029"]][[1]][3, ini_dt: end_dt]/fm[["2029"]][3]
fs[["2029"]][4, ini_dt: end_dt] <-   xfs[["2029"]][[1]][4, ini_dt: end_dt]/fm[["2029"]][4]

fs[["2029"]][5, ini_dt: end_dt] <-   xfs[["2029"]][[1]][5, ini_dt: end_dt]/xfs[["2029"]][[1]][5, ini_dt: end_dt]


# 2030

odd_num_test <- 1

Ccal[[2030]] <- lapply(Ccal[[2029]],function(x) x*odd_num_test)

frac_test[[2030]] <- frac_test[[2029]]

frac_ab[["2030"]] <- c(unlist(as.numeric(frac_test[[2030]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2030]]$P$reflex)))

dfList_NPPhaseIII_A <- dfList_NP_2029
for(i in param_var){  
  dfList_NPPhaseIII_A[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NPPhaseIII_A, index = i, 
                                        frac_testing = frac_test[[2030]],
                                        S_Yint = 2030, S_Yend = 2031, r_Yend = 2031, NPlst = NPlst, 
                                        fp = c(Ccal[[2030]]$C*fm[["2030"]][1], 
                                               Ccal[[2030]]$C*fm[["2030"]][2], 
                                               Ccal[[2030]]$P*fm[["2030"]][3], 
                                               Ccal[[2030]]$P*fm[["2030"]][4], 
                                               Ccal[[2030]]$P*fm[["2030"]][5]))
  
}

for(i in param_var){
  # begining of 2029
  b_pt <- (2031 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPPhaseIII_A[[i]])[3]
  dfList_NPPhaseIII_A[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2030"]] <- c(Ccal[[2030]]$C, Ccal[[2030]]$P)
n_ab_np[["2030"]] <- n_ab_np[["2029"]]*odd_num_test
xfs[["2030"]] <- fs_estimate(num_ab = n_ab_np[["2030"]], 
                             cov_np = coverage_np[["2030"]], 
                             frac_ab = frac_ab[["2030"]], 
                             fp = fm[["2030"]], year = 2030, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2030 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2030"]] <- fs[["2029"]]
fs[["2030"]][1, ini_dt: end_dt] <-   xfs[["2030"]][[1]][1, ini_dt: end_dt]/fm[["2030"]][1]
fs[["2030"]][2, ini_dt: end_dt] <-   xfs[["2030"]][[1]][2, ini_dt: end_dt]/fm[["2030"]][2]
fs[["2030"]][3, ini_dt: end_dt] <-   xfs[["2030"]][[1]][3, ini_dt: end_dt]/fm[["2030"]][3]
fs[["2030"]][4, ini_dt: end_dt] <-   xfs[["2030"]][[1]][4, ini_dt: end_dt]/fm[["2030"]][4]
fs[["2030"]][5, ini_dt: end_dt] <-   xfs[["2030"]][[1]][5, ini_dt: end_dt]/xfs[["2030"]][[1]][5, ini_dt: end_dt]








####NP Phase III (2027-2030: expand Phase II-B)-B #### 
# 2027

odd_num_test <- 1.112

Ccal[[2027]] <- lapply(Ccal[[2026]],function(x) x*odd_num_test)

frac_test[[2027]] <- frac_test[[2026]]

frac_ab[["2027"]] <- c(unlist(as.numeric(frac_test[[2027]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2027]]$P$reflex)))

dfList_NP_2027 <- dfList_NPPhaseII
for(i in param_var){  
  dfList_NP_2027[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2027, index = i, 
                                   frac_testing = frac_test[[2027]],
                                   S_Yint = 2027, S_Yend = 2028, r_Yend = 2028, NPlst = NPlst, 
                                   fp = c(Ccal[[2027]]$C*fm[["2027"]][1], 
                                          Ccal[[2027]]$C*fm[["2027"]][2], 
                                          Ccal[[2027]]$P*fm[["2027"]][3], 
                                          Ccal[[2027]]$P*fm[["2027"]][4], 
                                          Ccal[[2027]]$P*fm[["2027"]][5]))
  
}

for(i in param_var){
  # begining of 2027
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPPhaseII[[i]])[3]
  dfList_NP_2027[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2027"]] <- c(Ccal[[2027]]$C, Ccal[[2027]]$P)
n_ab_np[["2027"]] <- n_ab_np[["2026"]]*odd_num_test
xfs[["2027"]] <- fs_estimate(num_ab = n_ab_np[["2027"]], 
                             cov_np = coverage_np[["2027"]], 
                             frac_ab = frac_ab[["2027"]], 
                             fp = fm[["2027"]], year = 2027, endY = 100,
                             modsim = Sce_sq)

View(xfs[["2027"]][[1]])
ini_dt <- (2027 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2027 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2027"]] <- fs[["2026"]]
fs[["2027"]][1, ini_dt: end_dt] <-   xfs[["2027"]][[1]][1, ini_dt: end_dt]/fm[["2027"]][1]
fs[["2027"]][2, ini_dt: end_dt] <-   xfs[["2027"]][[1]][2, ini_dt: end_dt]/fm[["2027"]][2]
fs[["2027"]][3, ini_dt: end_dt] <-   xfs[["2027"]][[1]][3, ini_dt: end_dt]/fm[["2027"]][3]
fs[["2027"]][4, ini_dt: end_dt] <-   xfs[["2027"]][[1]][4, ini_dt: end_dt]/fm[["2027"]][4]
fs[["2027"]][5, ini_dt: end_dt] <-   xfs[["2027"]][[1]][5, ini_dt: end_dt]/xfs[["2027"]][[1]][5, ini_dt: end_dt]





# 2028 

odd_num_test <- 1.17
fm[["2028"]] <- c(0.8, 0.8, 12, 12, 1)

Ccal[[2028]] <- lapply(Ccal[[2027]],function(x) x*odd_num_test)

frac_test[[2028]] <- frac_test[[2026]]

frac_ab[["2028"]] <- c(unlist(as.numeric(frac_test[[2028]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2028]]$P$reflex)))

dfList_NP_2028 <- dfList_NP_2027
for(i in param_var){  
  dfList_NP_2028[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2028, index = i, 
                                   frac_testing = frac_test[[2028]],
                                   S_Yint = 2028, S_Yend = 2029, r_Yend = 2029, NPlst = NPlst, 
                                   fp = c(Ccal[[2028]]$C*fm[["2028"]][1], 
                                          Ccal[[2028]]$C*fm[["2028"]][2], 
                                          Ccal[[2028]]$P*fm[["2028"]][3], 
                                          Ccal[[2028]]$P*fm[["2028"]][4], 
                                          Ccal[[2028]]$P*fm[["2028"]][5]))
  
}

for(i in param_var){
  # begining of 2026
  b_pt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2028[[i]])[3]
  dfList_NP_2028[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2028"]] <- c(Ccal[[2028]]$C, Ccal[[2028]]$P)
n_ab_np[["2028"]] <- n_ab_np[["2027"]]*odd_num_test
xfs[["2028"]] <- fs_estimate(num_ab = n_ab_np[["2028"]], 
                             cov_np = coverage_np[["2028"]], 
                             frac_ab = frac_ab[["2028"]], 
                             fp = fm[["2028"]], year = 2028, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2028 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2028"]] <- fs[["2027"]]
fs[["2028"]][1, ini_dt: end_dt] <-   xfs[["2028"]][[1]][1, ini_dt: end_dt]/fm[["2028"]][1]
fs[["2028"]][2, ini_dt: end_dt] <-   xfs[["2028"]][[1]][2, ini_dt: end_dt]/fm[["2028"]][2]
fs[["2028"]][3, ini_dt: end_dt] <-   xfs[["2028"]][[1]][3, ini_dt: end_dt]/fm[["2028"]][3]
fs[["2028"]][4, ini_dt: end_dt] <-   xfs[["2028"]][[1]][4, ini_dt: end_dt]/fm[["2028"]][4]
fs[["2028"]][5, ini_dt: end_dt] <-   xfs[["2028"]][[1]][5, ini_dt: end_dt]/xfs[["2028"]][[1]][5, ini_dt: end_dt]

# 2029
fm[["2029"]] <- c(0.8, 0.8, 10, 10, 1)
fm[["2030"]] <- c(0.8, 0.8, 10, 10, 1)
odd_num_test <- 1.137

Ccal[[2029]] <- lapply(Ccal[[2028]],function(x) x*odd_num_test)

frac_test[[2029]] <- frac_test[[2028]]

frac_ab[["2029"]] <- c(unlist(as.numeric(frac_test[[2029]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2029]]$P$reflex)))

dfList_NP_2029 <- dfList_NP_2028
for(i in param_var){  
  dfList_NP_2029[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2029, index = i, 
                                   frac_testing = frac_test[[2029]],
                                   S_Yint = 2029, S_Yend = 2030, r_Yend = 2030, NPlst = NPlst, 
                                   fp = c(Ccal[[2029]]$C*fm[["2029"]][1], 
                                          Ccal[[2029]]$C*fm[["2029"]][2], 
                                          Ccal[[2029]]$P*fm[["2029"]][3], 
                                          Ccal[[2029]]$P*fm[["2029"]][4], 
                                          Ccal[[2029]]$P*fm[["2029"]][5]))
  
}


for(i in param_var){
  # begining of 2029
  b_pt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2029[[i]])[3]
  dfList_NP_2029[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2029"]] <- c(Ccal[[2029]]$C, Ccal[[2029]]$P)
n_ab_np[["2029"]] <- n_ab_np[["2028"]]*odd_num_test
xfs[["2029"]] <- fs_estimate(num_ab = n_ab_np[["2029"]], 
                             cov_np = coverage_np[["2029"]], 
                             frac_ab = frac_ab[["2029"]], 
                             fp = fm[["2029"]], year = 2029, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2029 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2029"]] <- fs[["2028"]]
fs[["2029"]][1, ini_dt: end_dt] <-   xfs[["2029"]][[1]][1, ini_dt: end_dt]/fm[["2029"]][1]
fs[["2029"]][2, ini_dt: end_dt] <-   xfs[["2029"]][[1]][2, ini_dt: end_dt]/fm[["2029"]][2]
fs[["2029"]][3, ini_dt: end_dt] <-   xfs[["2029"]][[1]][3, ini_dt: end_dt]/fm[["2029"]][3]
fs[["2029"]][4, ini_dt: end_dt] <-   xfs[["2029"]][[1]][4, ini_dt: end_dt]/fm[["2029"]][4]
fs[["2029"]][5, ini_dt: end_dt] <-   xfs[["2029"]][[1]][5, ini_dt: end_dt]/xfs[["2029"]][[1]][5, ini_dt: end_dt]


# 2030

odd_num_test <- 1.065

Ccal[[2030]] <- lapply(Ccal[[2029]],function(x) x*odd_num_test)

frac_test[[2030]] <- frac_test[[2029]]

frac_ab[["2030"]] <- c(unlist(as.numeric(frac_test[[2030]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2030]]$P$reflex)))

dfList_NPPhaseIII_B <- dfList_NP_2029
for(i in param_var){  
  dfList_NPPhaseIII_B[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NPPhaseIII_B, index = i, 
                                        frac_testing = frac_test[[2030]],
                                        S_Yint = 2030, S_Yend = 2031, r_Yend = 2031, NPlst = NPlst, 
                                        fp = c(Ccal[[2030]]$C*fm[["2030"]][1], 
                                               Ccal[[2030]]$C*fm[["2030"]][2], 
                                               Ccal[[2030]]$P*fm[["2030"]][3], 
                                               Ccal[[2030]]$P*fm[["2030"]][4], 
                                               Ccal[[2030]]$P*fm[["2030"]][5]))
  
}

for(i in param_var){
  # begining of 2029
  b_pt <- (2031 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPPhaseIII_B[[i]])[3]
  dfList_NPPhaseIII_B[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2030"]] <- c(Ccal[[2030]]$C, Ccal[[2030]]$P)
n_ab_np[["2030"]] <- n_ab_np[["2029"]]*odd_num_test
xfs[["2030"]] <- fs_estimate(num_ab = n_ab_np[["2030"]], 
                             cov_np = coverage_np[["2030"]], 
                             frac_ab = frac_ab[["2030"]], 
                             fp = fm[["2030"]], year = 2030, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2030 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["dfList_NPPhaseIII_B"]] <- fs[["2029"]]
fs[["dfList_NPPhaseIII_B"]][1, ini_dt: end_dt] <-   xfs[["2030"]][[1]][1, ini_dt: end_dt]/fm[["2030"]][1]
fs[["dfList_NPPhaseIII_B"]][2, ini_dt: end_dt] <-   xfs[["2030"]][[1]][2, ini_dt: end_dt]/fm[["2030"]][2]
fs[["dfList_NPPhaseIII_B"]][3, ini_dt: end_dt] <-   xfs[["2030"]][[1]][3, ini_dt: end_dt]/fm[["2030"]][3]
fs[["dfList_NPPhaseIII_B"]][4, ini_dt: end_dt] <-   xfs[["2030"]][[1]][4, ini_dt: end_dt]/fm[["2030"]][4]
fs[["dfList_NPPhaseIII_B"]][5, ini_dt: end_dt] <-   xfs[["2030"]][[1]][5, ini_dt: end_dt]/xfs[["2030"]][[1]][5, ini_dt: end_dt]




#### prison_testing_I: program sustained by reducing testing in community ####
# 2027

odd_num_test <- c(0.8, 1)

Ccal[[2027]] <- list("C" = Ccal[[2026]]$C*odd_num_test[1],
                     "P" = Ccal[[2026]]$P*odd_num_test[2])

frac_test[[2027]] <- frac_test[[2026]]

frac_ab[["2027"]] <- c(unlist(as.numeric(frac_test[[2027]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2027]]$P$reflex)))

dfList_NP_2027 <- dfList_NPPhaseII
for(i in param_var){  
  dfList_NP_2027[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2027, index = i, 
                                   frac_testing = frac_test[[2027]],
                                   S_Yint = 2027, S_Yend = 2028, r_Yend = 2028, NPlst = NPlst, 
                                   fp = c(Ccal[[2027]]$C*fm[["2027"]][1], 
                                          Ccal[[2027]]$C*fm[["2027"]][2], 
                                          Ccal[[2027]]$P*fm[["2027"]][3], 
                                          Ccal[[2027]]$P*fm[["2027"]][4], 
                                          Ccal[[2027]]$P*fm[["2027"]][5]))
  
}

for(i in param_var){
  # begining of 2027
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPPhaseII[[i]])[3]
  dfList_NP_2027[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2027"]] <- c(Ccal[[2027]]$C, Ccal[[2027]]$P)
n_ab_np[["2027"]] <- n_ab_np[["2026"]]*odd_num_test
xfs[["2027"]] <- fs_estimate(num_ab = n_ab_np[["2027"]], 
                             cov_np = coverage_np[["2027"]], 
                             frac_ab = frac_ab[["2027"]], 
                             fp = fm[["2027"]], year = 2027, endY = 100,
                             modsim = Sce_sq)

View(xfs[["2027"]][[1]])
ini_dt <- (2027 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2027 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2027"]] <- fs[["2026"]]
fs[["2027"]][1, ini_dt: end_dt] <-   xfs[["2027"]][[1]][1, ini_dt: end_dt]/fm[["2027"]][1]
fs[["2027"]][2, ini_dt: end_dt] <-   xfs[["2027"]][[1]][2, ini_dt: end_dt]/fm[["2027"]][2]
fs[["2027"]][3, ini_dt: end_dt] <-   xfs[["2027"]][[1]][3, ini_dt: end_dt]/fm[["2027"]][3]
fs[["2027"]][4, ini_dt: end_dt] <-   xfs[["2027"]][[1]][4, ini_dt: end_dt]/fm[["2027"]][4]
fs[["2027"]][5, ini_dt: end_dt] <-   xfs[["2027"]][[1]][5, ini_dt: end_dt]/xfs[["2027"]][[1]][5, ini_dt: end_dt]




# 2028 

odd_num_test <- c(0.9, 1.05)

Ccal[[2028]] <- list("C" = Ccal[[2027]]$C*odd_num_test[1],
                     "P" = Ccal[[2027]]$P*odd_num_test[2])

frac_test[[2028]] <- frac_test[[2026]]

frac_ab[["2028"]] <- c(unlist(as.numeric(frac_test[[2028]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2028]]$P$reflex)))

dfList_NP_2028 <- dfList_NP_2027
for(i in param_var){  
  dfList_NP_2028[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2028, index = i, 
                                   frac_testing = frac_test[[2028]],
                                   S_Yint = 2028, S_Yend = 2029, r_Yend = 2029, NPlst = NPlst, 
                                   fp = c(Ccal[[2028]]$C*fm[["2028"]][1], 
                                          Ccal[[2028]]$C*fm[["2028"]][2], 
                                          Ccal[[2028]]$P*fm[["2028"]][3], 
                                          Ccal[[2028]]$P*fm[["2028"]][4], 
                                          Ccal[[2028]]$P*fm[["2028"]][5]))
  
}

for(i in param_var){
  # begining of 2026
  b_pt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2028[[i]])[3]
  dfList_NP_2028[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2028"]] <- c(Ccal[[2028]]$C, Ccal[[2028]]$P)
n_ab_np[["2028"]] <- n_ab_np[["2027"]]*odd_num_test
xfs[["2028"]] <- fs_estimate(num_ab = n_ab_np[["2028"]], 
                             cov_np = coverage_np[["2028"]], 
                             frac_ab = frac_ab[["2028"]], 
                             fp = fm[["2028"]], year = 2028, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2028 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2028"]] <- fs[["2027"]]
fs[["2028"]][1, ini_dt: end_dt] <-   xfs[["2028"]][[1]][1, ini_dt: end_dt]/fm[["2028"]][1]
fs[["2028"]][2, ini_dt: end_dt] <-   xfs[["2028"]][[1]][2, ini_dt: end_dt]/fm[["2028"]][2]
fs[["2028"]][3, ini_dt: end_dt] <-   xfs[["2028"]][[1]][3, ini_dt: end_dt]/fm[["2028"]][3]
fs[["2028"]][4, ini_dt: end_dt] <-   xfs[["2028"]][[1]][4, ini_dt: end_dt]/fm[["2028"]][4]
fs[["2028"]][5, ini_dt: end_dt] <-   xfs[["2028"]][[1]][5, ini_dt: end_dt]/xfs[["2028"]][[1]][5, ini_dt: end_dt]

# 2029

odd_num_test <- c(0.9, 1)

Ccal[[2029]] <- list("C" = Ccal[[2028]]$C*odd_num_test[1],
                     "P" = Ccal[[2028]]$P*odd_num_test[2])

frac_test[[2029]] <- frac_test[[2028]]

frac_ab[["2029"]] <- c(unlist(as.numeric(frac_test[[2029]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2029]]$P$reflex)))

dfList_NP_2029 <- dfList_NP_2028
for(i in param_var){  
  dfList_NP_2029[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2029, index = i, 
                                   frac_testing = frac_test[[2029]],
                                   S_Yint = 2029, S_Yend = 2030, r_Yend = 2030, NPlst = NPlst, 
                                   fp = c(Ccal[[2029]]$C*fm[["2029"]][1], 
                                          Ccal[[2029]]$C*fm[["2029"]][2], 
                                          Ccal[[2029]]$P*fm[["2029"]][3], 
                                          Ccal[[2029]]$P*fm[["2029"]][4], 
                                          Ccal[[2029]]$P*fm[["2029"]][5]))
  
}

for(i in param_var){
  # begining of 2029
  b_pt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2029[[i]])[3]
  dfList_NP_2029[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2029"]] <- c(Ccal[[2029]]$C, Ccal[[2029]]$P)
n_ab_np[["2029"]] <- n_ab_np[["2028"]]*odd_num_test
xfs[["2029"]] <- fs_estimate(num_ab = n_ab_np[["2029"]], 
                             cov_np = coverage_np[["2029"]], 
                             frac_ab = frac_ab[["2029"]], 
                             fp = fm[["2029"]], year = 2029, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2029 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2029"]] <- fs[["2028"]]
fs[["2029"]][1, ini_dt: end_dt] <-   xfs[["2029"]][[1]][1, ini_dt: end_dt]/fm[["2029"]][1]
fs[["2029"]][2, ini_dt: end_dt] <-   xfs[["2029"]][[1]][2, ini_dt: end_dt]/fm[["2029"]][2]
fs[["2029"]][3, ini_dt: end_dt] <-   xfs[["2029"]][[1]][3, ini_dt: end_dt]/fm[["2029"]][3]
fs[["2029"]][4, ini_dt: end_dt] <-   xfs[["2029"]][[1]][4, ini_dt: end_dt]/fm[["2029"]][4]

fs[["2029"]][5, ini_dt: end_dt] <-   xfs[["2029"]][[1]][5, ini_dt: end_dt]/xfs[["2029"]][[1]][5, ini_dt: end_dt]


# 2030

odd_num_test <- c(0.9, 1)

Ccal[[2030]] <- list("C" = Ccal[[2029]]$C*odd_num_test[1],
                     "P" = Ccal[[2029]]$P*odd_num_test[2])

frac_test[[2030]] <- frac_test[[2029]]

frac_ab[["2030"]] <- c(unlist(as.numeric(frac_test[[2030]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2030]]$P$reflex)))

Prison_testing_I <- dfList_NP_2029
for(i in param_var){  
  Prison_testing_I[[i]] <- Param_cal(pj = POC_AU, dlist = Prison_testing_I, index = i, 
                                     frac_testing = frac_test[[2030]],
                                     S_Yint = 2030, S_Yend = 2031, r_Yend = 2031, NPlst = NPlst, 
                                     fp = c(Ccal[[2030]]$C*fm[["2030"]][1], 
                                            Ccal[[2030]]$C*fm[["2030"]][2], 
                                            Ccal[[2030]]$P*fm[["2030"]][3], 
                                            Ccal[[2030]]$P*fm[["2030"]][4], 
                                            Ccal[[2030]]$P*fm[["2030"]][5]))
  
}

for(i in param_var){
  # begining of 2029
  b_pt <- (2031 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(Prison_testing_I[[i]])[3]
  Prison_testing_I[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2030"]] <- c(Ccal[[2030]]$C, Ccal[[2030]]$P)
n_ab_np[["2030"]] <- n_ab_np[["2029"]]*odd_num_test
xfs[["2030"]] <- fs_estimate(num_ab = n_ab_np[["2030"]], 
                             cov_np = coverage_np[["2030"]], 
                             frac_ab = frac_ab[["2030"]], 
                             fp = fm[["2030"]], year = 2030, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2030 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["Prison_testing_I"]] <- fs[["2029"]]
fs[["Prison_testing_I"]][1, ini_dt: end_dt] <-   xfs[["2030"]][[1]][1, ini_dt: end_dt]/fm[["2030"]][1]
fs[["Prison_testing_I"]][2, ini_dt: end_dt] <-   xfs[["2030"]][[1]][2, ini_dt: end_dt]/fm[["2030"]][2]
fs[["Prison_testing_I"]][3, ini_dt: end_dt] <-   xfs[["2030"]][[1]][3, ini_dt: end_dt]/fm[["2030"]][3]
fs[["Prison_testing_I"]][4, ini_dt: end_dt] <-   xfs[["2030"]][[1]][4, ini_dt: end_dt]/fm[["2030"]][4]
fs[["Prison_testing_I"]][5, ini_dt: end_dt] <-   xfs[["2030"]][[1]][5, ini_dt: end_dt]/xfs[["2030"]][[1]][5, ini_dt: end_dt]



#### prison_testing_II ####
# 2027

# 2027

odd_num_test <-  c(1.03, 1.114)

Ccal[[2027]] <- list("C" = Ccal[[2026]]$C*odd_num_test[1],
                     "P" = Ccal[[2026]]$P*odd_num_test[2])

frac_test[[2027]] <- frac_test[[2026]]

frac_ab[["2027"]] <- c(unlist(as.numeric(frac_test[[2027]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2027]]$P$reflex)))

dfList_NP_2027 <- dfList_NPPhaseII
for(i in param_var){  
  dfList_NP_2027[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2027, index = i, 
                                   frac_testing = frac_test[[2027]],
                                   S_Yint = 2027, S_Yend = 2028, r_Yend = 2028, NPlst = NPlst, 
                                   fp = c(Ccal[[2027]]$C*fm[["2027"]][1], 
                                          Ccal[[2027]]$C*fm[["2027"]][2], 
                                          Ccal[[2027]]$P*fm[["2027"]][3], 
                                          Ccal[[2027]]$P*fm[["2027"]][4], 
                                          Ccal[[2027]]$P*fm[["2027"]][5]))
  
}

for(i in param_var){
  # begining of 2027
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPPhaseII[[i]])[3]
  dfList_NP_2027[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2027"]] <- c(Ccal[[2027]]$C, Ccal[[2027]]$P)
n_ab_np[["2027"]] <- n_ab_np[["2026"]]*odd_num_test
xfs[["2027"]] <- fs_estimate(num_ab = n_ab_np[["2027"]], 
                             cov_np = coverage_np[["2027"]], 
                             frac_ab = frac_ab[["2027"]], 
                             fp = fm[["2027"]], year = 2027, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2027 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2027 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2027"]] <- fs[["2026"]]
fs[["2027"]][1, ini_dt: end_dt] <-   xfs[["2027"]][[1]][1, ini_dt: end_dt]/fm[["2027"]][1]
fs[["2027"]][2, ini_dt: end_dt] <-   xfs[["2027"]][[1]][2, ini_dt: end_dt]/fm[["2027"]][2]
fs[["2027"]][3, ini_dt: end_dt] <-   xfs[["2027"]][[1]][3, ini_dt: end_dt]/fm[["2027"]][3]
fs[["2027"]][4, ini_dt: end_dt] <-   xfs[["2027"]][[1]][4, ini_dt: end_dt]/fm[["2027"]][4]
fs[["2027"]][5, ini_dt: end_dt] <-   xfs[["2027"]][[1]][5, ini_dt: end_dt]/xfs[["2027"]][[1]][5, ini_dt: end_dt]





# 2028 
odd_num_test <-  c(0.99, 1.1797)
Ccal[[2028]] <- list("C" = Ccal[[2027]]$C*odd_num_test[1],
                     "P" = Ccal[[2027]]$P*odd_num_test[2])
fm[["2028"]] <- c(0.8, 0.8, 12, 12, 1)


frac_test[[2028]] <- frac_test[[2026]]

frac_ab[["2028"]] <- c(unlist(as.numeric(frac_test[[2028]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2028]]$P$reflex)))

dfList_NP_2028 <- dfList_NP_2027
for(i in param_var){  
  dfList_NP_2028[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2028, index = i, 
                                   frac_testing = frac_test[[2028]],
                                   S_Yint = 2028, S_Yend = 2029, r_Yend = 2029, NPlst = NPlst, 
                                   fp = c(Ccal[[2028]]$C*fm[["2028"]][1], 
                                          Ccal[[2028]]$C*fm[["2028"]][2], 
                                          Ccal[[2028]]$P*fm[["2028"]][3], 
                                          Ccal[[2028]]$P*fm[["2028"]][4], 
                                          Ccal[[2028]]$P*fm[["2028"]][5]))
  
}

for(i in param_var){
  # begining of 2026
  b_pt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2028[[i]])[3]
  dfList_NP_2028[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2028"]] <- c(Ccal[[2028]]$C, Ccal[[2028]]$P)
n_ab_np[["2028"]] <- n_ab_np[["2027"]]*odd_num_test
xfs[["2028"]] <- fs_estimate(num_ab = n_ab_np[["2028"]], 
                             cov_np = coverage_np[["2028"]], 
                             frac_ab = frac_ab[["2028"]], 
                             fp = fm[["2028"]], year = 2028, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2028 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2028"]] <- fs[["2027"]]
fs[["2028"]][1, ini_dt: end_dt] <-   xfs[["2028"]][[1]][1, ini_dt: end_dt]/fm[["2028"]][1]
fs[["2028"]][2, ini_dt: end_dt] <-   xfs[["2028"]][[1]][2, ini_dt: end_dt]/fm[["2028"]][2]
fs[["2028"]][3, ini_dt: end_dt] <-   xfs[["2028"]][[1]][3, ini_dt: end_dt]/fm[["2028"]][3]
fs[["2028"]][4, ini_dt: end_dt] <-   xfs[["2028"]][[1]][4, ini_dt: end_dt]/fm[["2028"]][4]
fs[["2028"]][5, ini_dt: end_dt] <-   xfs[["2028"]][[1]][5, ini_dt: end_dt]/xfs[["2028"]][[1]][5, ini_dt: end_dt]

# 2029
fm[["2029"]] <- c(0.8, 0.8, 10, 10, 1)
fm[["2030"]] <- c(0.8, 0.8, 10, 10, 1)

odd_num_test <- c(0.95,1.14)
Ccal[[2029]] <- list("C" = Ccal[[2028]]$C*odd_num_test[1],
                     "P" = Ccal[[2028]]$P*odd_num_test[2])



frac_test[[2029]] <- frac_test[[2028]]

frac_ab[["2029"]] <- c(unlist(as.numeric(frac_test[[2029]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2029]]$P$reflex)))

dfList_NP_2029 <- dfList_NP_2028
for(i in param_var){  
  dfList_NP_2029[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2029, index = i, 
                                   frac_testing = frac_test[[2029]],
                                   S_Yint = 2029, S_Yend = 2030, r_Yend = 2030, NPlst = NPlst, 
                                   fp = c(Ccal[[2029]]$C*fm[["2029"]][1], 
                                          Ccal[[2029]]$C*fm[["2029"]][2], 
                                          Ccal[[2029]]$P*fm[["2029"]][3], 
                                          Ccal[[2029]]$P*fm[["2029"]][4], 
                                          Ccal[[2029]]$P*fm[["2029"]][5]))
  
}


for(i in param_var){
  # begining of 2029
  b_pt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2029[[i]])[3]
  dfList_NP_2029[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2029"]] <- c(Ccal[[2029]]$C, Ccal[[2029]]$P)
n_ab_np[["2029"]] <- n_ab_np[["2028"]]*odd_num_test
xfs[["2029"]] <- fs_estimate(num_ab = n_ab_np[["2029"]], 
                             cov_np = coverage_np[["2029"]], 
                             frac_ab = frac_ab[["2029"]], 
                             fp = fm[["2029"]], year = 2029, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2029 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2029"]] <- fs[["2028"]]
fs[["2029"]][1, ini_dt: end_dt] <-   xfs[["2029"]][[1]][1, ini_dt: end_dt]/fm[["2029"]][1]
fs[["2029"]][2, ini_dt: end_dt] <-   xfs[["2029"]][[1]][2, ini_dt: end_dt]/fm[["2029"]][2]
fs[["2029"]][3, ini_dt: end_dt] <-   xfs[["2029"]][[1]][3, ini_dt: end_dt]/fm[["2029"]][3]
fs[["2029"]][4, ini_dt: end_dt] <-   xfs[["2029"]][[1]][4, ini_dt: end_dt]/fm[["2029"]][4]
fs[["2029"]][5, ini_dt: end_dt] <-   xfs[["2029"]][[1]][5, ini_dt: end_dt]/xfs[["2029"]][[1]][5, ini_dt: end_dt]


# 2030

odd_num_test <- c(0.95, 1.18)
Ccal[[2030]] <- list("C" = Ccal[[2029]]$C*odd_num_test[1],
                     "P" = Ccal[[2028]]$P*odd_num_test[2])


frac_test[[2030]] <- frac_test[[2029]]

frac_ab[["2030"]] <- c(unlist(as.numeric(frac_test[[2030]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2030]]$P$reflex)))

Prison_testing_II <- dfList_NP_2029
for(i in param_var){  
  Prison_testing_II[[i]] <- Param_cal(pj = POC_AU, dlist = Prison_testing_II, index = i, 
                                      frac_testing = frac_test[[2030]],
                                      S_Yint = 2030, S_Yend = 2031, r_Yend = 2031, NPlst = NPlst, 
                                      fp = c(Ccal[[2030]]$C*fm[["2030"]][1], 
                                             Ccal[[2030]]$C*fm[["2030"]][2], 
                                             Ccal[[2030]]$P*fm[["2030"]][3], 
                                             Ccal[[2030]]$P*fm[["2030"]][4], 
                                             Ccal[[2030]]$P*fm[["2030"]][5]))
  
}

for(i in param_var){
  # begining of 2029
  b_pt <- (2031 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(Prison_testing_II[[i]])[3]
  Prison_testing_II[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2030"]] <- c(Ccal[[2030]]$C, Ccal[[2030]]$P)
n_ab_np[["2030"]] <- n_ab_np[["2029"]]*odd_num_test
xfs[["2030"]] <- fs_estimate(num_ab = n_ab_np[["2030"]], 
                             cov_np = coverage_np[["2030"]], 
                             frac_ab = frac_ab[["2030"]], 
                             fp = fm[["2030"]], year = 2030, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2030 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["Prison_testing_II"]] <- fs[["2029"]]
fs[["Prison_testing_II"]][1, ini_dt: end_dt] <-   xfs[["2030"]][[1]][1, ini_dt: end_dt]/fm[["2030"]][1]
fs[["Prison_testing_II"]][2, ini_dt: end_dt] <-   xfs[["2030"]][[1]][2, ini_dt: end_dt]/fm[["2030"]][2]
fs[["Prison_testing_II"]][3, ini_dt: end_dt] <-   xfs[["2030"]][[1]][3, ini_dt: end_dt]/fm[["2030"]][3]
fs[["Prison_testing_II"]][4, ini_dt: end_dt] <-   xfs[["2030"]][[1]][4, ini_dt: end_dt]/fm[["2030"]][4]
fs[["Prison_testing_II"]][5, ini_dt: end_dt] <-   xfs[["2030"]][[1]][5, ini_dt: end_dt]/xfs[["2030"]][[1]][5, ini_dt: end_dt]



####sustained-scale-up: community focus ####
#### prison_testing_III ####
# 2027

# 2027

odd_num_test <-  c(1.82, 1)

Ccal[[2027]] <- list("C" = Ccal[[2026]]$C*odd_num_test[1],
                     "P" = Ccal[[2026]]$P*odd_num_test[2])

frac_test[[2027]] <- frac_test[[2026]]

frac_ab[["2027"]] <- c(unlist(as.numeric(frac_test[[2027]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2027]]$P$reflex)))

dfList_NP_2027 <- dfList_NPPhaseII
for(i in param_var){  
  dfList_NP_2027[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2027, index = i, 
                                   frac_testing = frac_test[[2027]],
                                   S_Yint = 2027, S_Yend = 2028, r_Yend = 2028, NPlst = NPlst, 
                                   fp = c(Ccal[[2027]]$C*fm[["2027"]][1], 
                                          Ccal[[2027]]$C*fm[["2027"]][2], 
                                          Ccal[[2027]]$P*fm[["2027"]][3], 
                                          Ccal[[2027]]$P*fm[["2027"]][4], 
                                          Ccal[[2027]]$P*fm[["2027"]][5]))
  
}

for(i in param_var){
  # begining of 2027
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPPhaseII[[i]])[3]
  dfList_NP_2027[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2027"]] <- c(Ccal[[2027]]$C, Ccal[[2027]]$P)
n_ab_np[["2027"]] <- n_ab_np[["2026"]]*odd_num_test
xfs[["2027"]] <- fs_estimate(num_ab = n_ab_np[["2027"]], 
                             cov_np = coverage_np[["2027"]], 
                             frac_ab = frac_ab[["2027"]], 
                             fp = fm[["2027"]], year = 2027, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2027 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2027 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2027"]] <- fs[["2026"]]
fs[["2027"]][1, ini_dt: end_dt] <-   xfs[["2027"]][[1]][1, ini_dt: end_dt]/fm[["2027"]][1]
fs[["2027"]][2, ini_dt: end_dt] <-   xfs[["2027"]][[1]][2, ini_dt: end_dt]/fm[["2027"]][2]
fs[["2027"]][3, ini_dt: end_dt] <-   xfs[["2027"]][[1]][3, ini_dt: end_dt]/fm[["2027"]][3]
fs[["2027"]][4, ini_dt: end_dt] <-   xfs[["2027"]][[1]][4, ini_dt: end_dt]/fm[["2027"]][4]
fs[["2027"]][5, ini_dt: end_dt] <-   xfs[["2027"]][[1]][5, ini_dt: end_dt]/xfs[["2027"]][[1]][5, ini_dt: end_dt]





# 2028 
odd_num_test <-  c(1.73, 1)
Ccal[[2028]] <- list("C" = Ccal[[2027]]$C*odd_num_test[1],
                     "P" = Ccal[[2027]]$P*odd_num_test[2])
fm[["2028"]] <- c(0.8, 0.8, 12, 12, 1)


frac_test[[2028]] <- frac_test[[2026]]

frac_ab[["2028"]] <- c(unlist(as.numeric(frac_test[[2028]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2028]]$P$reflex)))

dfList_NP_2028 <- dfList_NP_2027
for(i in param_var){  
  dfList_NP_2028[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2028, index = i, 
                                   frac_testing = frac_test[[2028]],
                                   S_Yint = 2028, S_Yend = 2029, r_Yend = 2029, NPlst = NPlst, 
                                   fp = c(Ccal[[2028]]$C*fm[["2028"]][1], 
                                          Ccal[[2028]]$C*fm[["2028"]][2], 
                                          Ccal[[2028]]$P*fm[["2028"]][3], 
                                          Ccal[[2028]]$P*fm[["2028"]][4], 
                                          Ccal[[2028]]$P*fm[["2028"]][5]))
  
}

for(i in param_var){
  # begining of 2026
  b_pt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2028[[i]])[3]
  dfList_NP_2028[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2028"]] <- c(Ccal[[2028]]$C, Ccal[[2028]]$P)
n_ab_np[["2028"]] <- n_ab_np[["2027"]]*odd_num_test
xfs[["2028"]] <- fs_estimate(num_ab = n_ab_np[["2028"]], 
                             cov_np = coverage_np[["2028"]], 
                             frac_ab = frac_ab[["2028"]], 
                             fp = fm[["2028"]], year = 2028, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2028 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2028"]] <- fs[["2027"]]
fs[["2028"]][1, ini_dt: end_dt] <-   xfs[["2028"]][[1]][1, ini_dt: end_dt]/fm[["2028"]][1]
fs[["2028"]][2, ini_dt: end_dt] <-   xfs[["2028"]][[1]][2, ini_dt: end_dt]/fm[["2028"]][2]
fs[["2028"]][3, ini_dt: end_dt] <-   xfs[["2028"]][[1]][3, ini_dt: end_dt]/fm[["2028"]][3]
fs[["2028"]][4, ini_dt: end_dt] <-   xfs[["2028"]][[1]][4, ini_dt: end_dt]/fm[["2028"]][4]
fs[["2028"]][5, ini_dt: end_dt] <-   xfs[["2028"]][[1]][5, ini_dt: end_dt]/xfs[["2028"]][[1]][5, ini_dt: end_dt]

# 2029
fm[["2029"]] <- c(0.8, 0.8, 10, 10, 1)
fm[["2030"]] <- c(0.8, 0.8, 10, 10, 1)

odd_num_test <- c(1.04,1)
Ccal[[2029]] <- list("C" = Ccal[[2028]]$C*odd_num_test[1],
                     "P" = Ccal[[2028]]$P*odd_num_test[2])



frac_test[[2029]] <- frac_test[[2028]]

frac_ab[["2029"]] <- c(unlist(as.numeric(frac_test[[2029]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2029]]$P$reflex)))

dfList_NP_2029 <- dfList_NP_2028
for(i in param_var){  
  dfList_NP_2029[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2029, index = i, 
                                   frac_testing = frac_test[[2029]],
                                   S_Yint = 2029, S_Yend = 2030, r_Yend = 2030, NPlst = NPlst, 
                                   fp = c(Ccal[[2029]]$C*fm[["2029"]][1], 
                                          Ccal[[2029]]$C*fm[["2029"]][2], 
                                          Ccal[[2029]]$P*fm[["2029"]][3], 
                                          Ccal[[2029]]$P*fm[["2029"]][4], 
                                          Ccal[[2029]]$P*fm[["2029"]][5]))
  
}


for(i in param_var){
  # begining of 2029
  b_pt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2029[[i]])[3]
  dfList_NP_2029[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2029"]] <- c(Ccal[[2029]]$C, Ccal[[2029]]$P)
n_ab_np[["2029"]] <- n_ab_np[["2028"]]*odd_num_test
xfs[["2029"]] <- fs_estimate(num_ab = n_ab_np[["2029"]], 
                             cov_np = coverage_np[["2029"]], 
                             frac_ab = frac_ab[["2029"]], 
                             fp = fm[["2029"]], year = 2029, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2029 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2029"]] <- fs[["2028"]]
fs[["2029"]][1, ini_dt: end_dt] <-   xfs[["2029"]][[1]][1, ini_dt: end_dt]/fm[["2029"]][1]
fs[["2029"]][2, ini_dt: end_dt] <-   xfs[["2029"]][[1]][2, ini_dt: end_dt]/fm[["2029"]][2]
fs[["2029"]][3, ini_dt: end_dt] <-   xfs[["2029"]][[1]][3, ini_dt: end_dt]/fm[["2029"]][3]
fs[["2029"]][4, ini_dt: end_dt] <-   xfs[["2029"]][[1]][4, ini_dt: end_dt]/fm[["2029"]][4]
fs[["2029"]][5, ini_dt: end_dt] <-   xfs[["2029"]][[1]][5, ini_dt: end_dt]/xfs[["2029"]][[1]][5, ini_dt: end_dt]


# 2030

odd_num_test <- c(1.365, 1)
Ccal[[2030]] <- list("C" = Ccal[[2029]]$C*odd_num_test[1],
                     "P" = Ccal[[2028]]$P*odd_num_test[2])


frac_test[[2030]] <- frac_test[[2029]]

frac_ab[["2030"]] <- c(unlist(as.numeric(frac_test[[2030]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2030]]$P$reflex)))

Prison_testing_III <- dfList_NP_2029
for(i in param_var){  
  Prison_testing_III[[i]] <- Param_cal(pj = POC_AU, dlist = Prison_testing_III, index = i, 
                                       frac_testing = frac_test[[2030]],
                                       S_Yint = 2030, S_Yend = 2031, r_Yend = 2031, NPlst = NPlst, 
                                       fp = c(Ccal[[2030]]$C*fm[["2030"]][1], 
                                              Ccal[[2030]]$C*fm[["2030"]][2], 
                                              Ccal[[2030]]$P*fm[["2030"]][3], 
                                              Ccal[[2030]]$P*fm[["2030"]][4], 
                                              Ccal[[2030]]$P*fm[["2030"]][5]))
  
}

for(i in param_var){
  # begining of 2029
  b_pt <- (2031 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(Prison_testing_III[[i]])[3]
  Prison_testing_III[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2030"]] <- c(Ccal[[2030]]$C, Ccal[[2030]]$P)
n_ab_np[["2030"]] <- n_ab_np[["2029"]]*odd_num_test
xfs[["2030"]] <- fs_estimate(num_ab = n_ab_np[["2030"]], 
                             cov_np = coverage_np[["2030"]], 
                             frac_ab = frac_ab[["2030"]], 
                             fp = fm[["2030"]], year = 2030, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2030 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["Prison_testing_III"]] <- fs[["2029"]]
fs[["Prison_testing_III"]][1, ini_dt: end_dt] <-   xfs[["2030"]][[1]][1, ini_dt: end_dt]/fm[["2030"]][1]
fs[["Prison_testing_III"]][2, ini_dt: end_dt] <-   xfs[["2030"]][[1]][2, ini_dt: end_dt]/fm[["2030"]][2]
fs[["Prison_testing_III"]][3, ini_dt: end_dt] <-   xfs[["2030"]][[1]][3, ini_dt: end_dt]/fm[["2030"]][3]
fs[["Prison_testing_III"]][4, ini_dt: end_dt] <-   xfs[["2030"]][[1]][4, ini_dt: end_dt]/fm[["2030"]][4]
fs[["Prison_testing_III"]][5, ini_dt: end_dt] <-   xfs[["2030"]][[1]][5, ini_dt: end_dt]/xfs[["2030"]][[1]][5, ini_dt: end_dt]



scenario_cascade <- list(
  "Prison_testing_I" = Prison_testing_I, 
  "Prison_testing_II" = Prison_testing_II,
  "Prison_testing_III" = Prison_testing_III, 
  "dfList_NPPhaseIII_A" = dfList_NPPhaseIII_A, 
  "dfList_NPPhaseIII_B" = dfList_NPPhaseIII_B)

scenario_fc <- list( 
  "Prison_testing_I" = fs[["Prison_testing_I"]], 
  "Prison_testing_II" = fs[["Prison_testing_II"]],
  "Prison_testing_III" = fs[["Prison_testing_III"]],
  "dfList_NPPhaseIII_A" = fs[["2030"]], 
  "dfList_NPPhaseIII_B" = fs[["dfList_NPPhaseIII_B"]])



#===============================================================================
#
#                               run simulations
#
#===============================================================================
endY <- 100
Sce_np <- list()
tic <- proc.time()

for (scenario in names(scenario_cascade)) { 
  Sce_np[[scenario]] <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                               disease_progress, pop_array,
                               dfList,  
                               param_cascade_sc = scenario_cascade[[scenario]], 
                               fib = fib, 
                               modelrun = "UN", proj = "POC_AU", end_Y = endY, 
                               cost = NULL, costflow = NULL, 
                               costflow_Neg = NULL, 
                               fc = scenario_fc[[scenario]])
}

toc <- proc.time() - tic
# print(paste0("Completed: ", cost_type, " | Time: "))
print(toc)
save(Sce_sq, Sce_np,
     file = file.path(OutputFolder,
                      paste0("Simulations_",".rda")))

test <- list()
cl_ext <- names(Sce_np)[c(10:22)]
for(i in cl_ext){
  
  test[[i]] <- modres.flow.t(POC_AU, Sce_np, endYear = 100, 
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
  mutate(scenario = "prisons_testing_I")
dtp <- c(6667,6667, 11476,11476,20393, 20393, 12000, 13000, 11000, 14000, 10000, 15000,9000, 16000,8000, 17000,7000, 18000 )
# View(test_fscal%>%filter(year%in% c(7:15)))
View(test_fscal%>%filter(year%in% c(7:15))%>%group_by(year, setting)%>%
       summarise(tot_NP_test = sum(value))%>%
       mutate(year = year + 2015))
test%>%mutate(setting = ifelse(population %in% POC_AU$popNames[3:5], "P", "C"))%>%
  group_by(setting, year)%>%summarise(best = sum(newTreatment_sc))%>%
  filter(year%in% c(7:12))
