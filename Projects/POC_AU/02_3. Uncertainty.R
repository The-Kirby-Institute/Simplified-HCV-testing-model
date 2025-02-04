# uncertainty set up 
# parameter set
# Number of parameter sets to run
rm(list = ls())
gc()
project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Documents/Simplified-HCV-testing-model")

data_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                   project_name, sep = "")
# Load useful libraries
library("lhs")
library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")
library("ggplot2")
library("viridis") 
library("openxlsx")
library("gt")
library("dplyr")
Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
OutputFolder <- file.path(data_path, "02. Output")
OutputFig <- file.path(OutputFolder, "Figs")

load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "cali.rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_timev.rda")))
load(file.path(OutputFolder, paste0(project_name, "scenario_cascade.rda")))
source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotManuscript.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R")) 

source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 

urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
runSamples <- TRUE
saveAsBase <- TRUE  # if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing
number_samples <- 1000
start_lower <- 0.75
start_upper <- 1.25
end_lower <- 0.75
end_upper <- 1.25


POC_AU$numberSamples <- number_samples


# LHS sampling 
# test for the constant parameters 
## constants beta and mortality_b
param_constant <- read.csv(file.path(DataFolder, "parameters_constants.csv"), 
                           header = TRUE)%>%dplyr::select(-"X")%>%as.data.frame() 
# disease progress 

# import UL & LL datasets
# UL 
disease_progUL <- read.csv(file.path(DataFolder, "diseaseProgress_UU.csv"), 
                           header = TRUE)
disease_progressUL <- c(disease_progUL[, -1])%>%as_tibble()

disease_progressUL <- disease_progressUL%>%
  mutate(a_f0 = disease_progress$a_f0*start_upper,
         f3_hcc = disease_progress$f3_hcc*start_upper)


# LL 
disease_progLL <- read.csv(file.path(DataFolder, "diseaseProgress_LL.csv"), 
                           header = TRUE)
disease_progressLL <- c(disease_progLL[, -1])%>%as_tibble()

disease_progressLL <- disease_progressLL%>%
  mutate(a_f0 = disease_progress$a_f0*start_lower,
         f3_hcc = disease_progress$f3_hcc*start_lower)


# cured disease progress 
cured_progUL <- read.csv(file.path(DataFolder, "curedProgress_UU.csv"), 
                         header = TRUE)
cured_progressUL <- c(cured_progUL[, -1])%>%as_tibble()

cured_progressUL <- cured_progressUL%>%
  mutate(lt_plt = fib$lt_plt*start_upper,
         lt_cured_plt_cured = fib$lt_cured_plt_cured*start_upper)


# LL 
cured_progLL <- read.csv(file.path(DataFolder, "curedProgress_LL.csv"), 
                         header = TRUE)
cured_progressLL <- c(cured_progLL[, -1])%>%as_tibble()

cured_progressLL <- cured_progressLL%>%
  mutate(lt_plt = fib$lt_plt*start_lower,
         lt_cured_plt_cured = fib$lt_cured_plt_cured*start_lower)



set.seed(123456)
names(disease_progress)
lhs_samples <- randomLHS(n = POC_AU$numberSamples, 
                         k = nrow(param_constant) +  
                           dim(disease_progress)[2] + dim(fib)[2]  + 1 + POC_AU$npops)
colnames(lhs_samples) <- c(param_constant$parameter, names(disease_progress), 
                           names(fib), "poparray", 
                           "C_PWID", "C_fPWID", "P_PWID", "P_fPWID", "P_nPWID")

param_constant[is.na(param_constant$lower),"lower"] <- start_lower
param_constant[is.na(param_constant$upper),"upper"] <- start_upper

param_constant <- param_constant%>%mutate(lower = ifelse(lower == start_lower, start_lower*value, lower),
                                          upper = ifelse(upper == start_upper, start_upper*value, upper))

# paramset pop size in each compartment
init_pop <- filter(initialPops, parameter == "init_pop")$value

pop_prop <- initialPops%>%filter(parameter%in% c("pop_prop1", "pop_prop2", 
                                                 "pop_prop3", "pop_prop4", 
                                                 "pop_prop5"))%>%
  dplyr::select(value)%>%unlist()%>%as.vector()

popProp <- as.numeric(init_pop)*pop_prop 
init_prop_I <- c(best_estimates$HCVP1[1], best_estimates$HCVP2[1], 
                 best_estimates$HCVP3[1], best_estimates$HCVP4[1],
                 best_estimates$HCVP5[1])
init_prop_S <- c(1 - init_prop_I) 

estPops <- read.csv(file.path(paste0(DataFolder, "/Estimate_initial_pop.csv")), header = TRUE)

# Scale the samples to the specified ranges
x <- list()
randompara <- lhs_samples[,c(param_constant$parameter)]
parampop <- lhs_samples[, c(POC_AU$popNames)]
y <- matrix(NA, nrow = POC_AU$npops, ncol = length(names(disease_progress)))
colnames(y) <- names(disease_progress)

Param_disease_progress <- rep(list(y), POC_AU$numberSamples) 

z <- matrix(NA, nrow = POC_AU$npops, ncol = length(names(fib)))

colnames(z) <- names(fib)

Param_fib <- rep(list(z), POC_AU$numberSamples) 

disease_progparam <- list()
cured_progparam <- list()
parampopsize <- list()
for(i in 1:POC_AU$numberSamples){ 
  x[[i]] <- matrix(NA, nrow = 1, ncol = nrow(param_constant))
  
  x[[i]] <- randompara[i,] *(param_constant$upper - param_constant$lower) + param_constant$lower
  
  parampopsize[[i]] <- matrix(NA, nrow = POC_AU$npops, ncol = POC_AU$npops)
  
  parampopsize[[i]] <- parampop[i, ]*(popProp*start_upper - popProp*start_lower) +popProp*start_lower
  
  
  disease_progparam[[i]] <- do.call(rbind, replicate(POC_AU$npops, 
                                                     lhs_samples[i,c(names(disease_progress))], simplify=FALSE))
  
  cured_progparam[[i]] <- do.call(rbind, replicate(POC_AU$npops, 
                                                   lhs_samples[i,c(names(fib))], simplify=FALSE))
  
  Param_disease_progress[[i]] <- disease_progparam[[i]]*(disease_progressUL - disease_progressLL) + disease_progressLL

  Param_fib[[i]] <- cured_progparam[[i]]*(cured_progressUL - cured_progressLL) + cured_progressLL
  }

Param_Cparam <- lapply(x, function(m) as_tibble(t(m)))
Param_estimates <- lapply(Param_Cparam, function(x) x[rep(1,POC_AU$npts), ])

param_poparray <- rep(list(pop_array), POC_AU$numberSamples) 
pop_array_LL <- pop_array*start_lower
pop_array_UU <- pop_array*start_upper


Param_estPops <- list()

Param_Pops <- list()
for ( i in 1: POC_AU$numberSamples){ 
  Param_estPops[[i]] <- estPops%>%mutate(
    pop_group = rep(c(parampopsize[[i]]),dim(estPops)[1]/POC_AU$npops),
    SIprop = ifelse(estPops$SI=="S", 
                    rep(init_prop_S, POC_AU$diseaseprogress_n*POC_AU$npops),
                    rep(init_prop_I, POC_AU$ncomponent*POC_AU$npops - 
                          POC_AU$diseaseprogress_n*POC_AU$npops)),
    est_pop = value*pop_group*SIprop)
  
  Param_Pops[[i]] <- as.matrix(as.data.frame(matrix(Param_estPops[[i]]$est_pop, 
                                                    ncol = POC_AU$ncomponent,  
                                                    nrow = POC_AU$npops)))
}

Param_Pops <- lapply(Param_Pops, "colnames<-", c(POC_AU$component_name))

paramDflist <- rep(list(dfList), POC_AU$numberSamples) 
dflist_LL <- lapply(dfList, function(x) x*start_lower)
dflist_UU <- lapply(dfList, function(x) x*start_upper)
for(i in 1:POC_AU$numberSamples){ 
  param_poparray[[i]] <- lhs_samples[i,"poparray"]* (pop_array_UU - pop_array_LL) + pop_array_LL
  
  param_poparray[[i]][param_poparray[[i]]>1] <- 1 
  
}
 
save(Param_estimates,
     Param_disease_progress, 
     param_poparray, 
     parampopsize,
     Param_Pops,
     Param_fib,
     lhs_samples,
     file = file.path(OutputFolder,
                      paste0(project_name, "param", ".rda")))

rm(Param_estimates, Param_disease_progress, param_poparray, 
   parampopsize,
   Param_Pops,
   Param_fib)

gc()
# cured and lota using the duration of treatment, don't random sampling 
for(i in 1:POC_AU$numberSamples){
  paramDflist[[i]][["cured"]] <- dfList$cured

}
for(i in 1:POC_AU$numberSamples){
  paramDflist[[i]][["eta"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["eta"]] - dflist_LL[["eta"]]) + dflist_LL[["eta"]] 
  paramDflist[[i]][["eta"]][paramDflist[[i]][["eta"]]>1] <- 1                               
} 

for(i in 1:POC_AU$numberSamples){
  paramDflist[[i]][["lota"]] <- dfList$lota 

}
for(i in 1:POC_AU$numberSamples){
  paramDflist[[i]][["rho"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["rho"]] - dflist_LL[["rho"]]) + dflist_LL[["rho"]] 
  paramDflist[[i]][["rho"]][paramDflist[[i]][["rho"]]>1] <- 1  
}
for(i in 1:POC_AU$numberSamples){
  paramDflist[[i]][["tau_ab"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["tau_ab"]] - dflist_LL[["tau_ab"]]) + dflist_LL[["tau_ab"]] 
  paramDflist[[i]][["tau_ab"]][paramDflist[[i]][["tau_ab"]]>1] <- 1  
}
for(i in 1:POC_AU$numberSamples){  
  paramDflist[[i]][["tau_poct"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["tau_poct"]] - dflist_LL[["tau_poct"]]) + dflist_LL[["tau_poct"]] 
  paramDflist[[i]][["tau_poct"]][paramDflist[[i]][["tau_poct"]]>1] <- 1  
}
for(i in 1:POC_AU$numberSamples){   
  paramDflist[[i]][["tau_RNA"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["tau_RNA"]] - dflist_LL[["tau_RNA"]]) + dflist_LL[["tau_RNA"]] 
  paramDflist[[i]][["tau_RNA"]][paramDflist[[i]][["tau_RNA"]]>1] <- 1
}

save(paramDflist,
     file = file.path(OutputFolder,
                      paste0(project_name, "paramDflist", ".rda")))

save(lhs_samples, file = file.path(OutputFolder,
                                   paste0(project_name, "lhs_sampling", ".rda")))

rm(paramDflist)
gc()



files <- list.files(path = paste0(DataFolder, 
                                  "/cost/", sep =  ""), pattern = '*.csv')

# parameter sets for cost data 
# +- 10% 
costdfList <- lapply(files, function(f) {
  
  df <- read.csv(file.path(paste0(DataFolder, "/cost/", f, sep = "")), header = TRUE)
  
  df <- df[, -1]
  
  df <- df%>%as_tibble()
  
  df <- as.matrix(df, nrow = npops, ncol = length(.) + 1)
  
})
costdfList$QALYPops_LL
names(costdfList) <- c(gsub("^|.csv", "", files)) # ^: from beginning, \ end before .csv


cost_state <- costdfList$state
costflow <- list()
costflow[[1]] <- costdfList$costFlow
costflow[[2]] <- costdfList$costFlow_POCRNA

costflow_Neg <- list()
costflow_Neg[[1]] <- costdfList$costFlow_NEG
costflow_Neg[[2]] <- costdfList$`costFlow_POCRNA _NEG`



set.seed(123456) 
rand_multiply <- runif(number_samples, 0.9, 1.1)

param_cost <- lapply(rand_multiply, function(x) lapply(costdfList, function(y) y*x))

for(i in 1: length(rand_multiply)){ 
  names(param_cost[[i]]) <- names(costdfList)
  param_cost[[i]]$QALY <- lhs_samples[i,"poparray"]*(costdfList$QALYPops_UU - costdfList$QALYPops_LL) + costdfList$QALYPops_LL 
  }

param_cost_flow <- list()
param_costflow_Neg <- list()
param_QALY <- list()
for(i in 1: number_samples){ 
  param_cost_flow[[i]] <- list(param_cost[[i]]$costFlow, 
                               param_cost[[i]]$costFlow_POCRNA)
  
  param_costflow_Neg[[i]] <- list(param_cost[[i]]$costFlow_NEG,
                                  param_cost[[i]]$`costFlow_POCRNA _NEG`)
    
                                   
  param_QALY[[i]] <-  param_cost[[i]]$QALY  
  
} 


save(param_cost,
     param_cost_flow, 
     param_costflow_Neg, 
     param_QALY,
     rand_multiply ,
     file = file.path(OutputFolder,
                      paste0(project_name, "param_cost", ".rda")))
  

