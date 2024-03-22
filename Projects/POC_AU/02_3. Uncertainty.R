# uncertainty set up 
# parameter set
# Number of parameter sets to run
rm(list = ls())

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
library(gt)
library(dplyr)
Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
OutputFolder <- file.path(data_path, "02. Output")
OutputFig <- file.path(OutputFolder, "Figs")

load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "cali.rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_timev.rda")))
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
                           header = TRUE)%>%select(-"X")%>%as.data.frame() 
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
                           dim(disease_progress)[2] + dim(fib)[2]  + 1)
colnames(lhs_samples) <- c(param_constant$parameter, names(disease_progress), 
                           names(fib), "poparray")

param_constant[is.na(param_constant$lower),"lower"] <- start_lower
param_constant[is.na(param_constant$upper),"upper"] <- start_upper

param_constant <- param_constant%>%mutate(lower = ifelse(lower == start_lower, start_lower*value, lower),
                                          upper = ifelse(upper == start_upper, start_upper*value, upper))



# Scale the samples to the specified ranges
x <- list()
randompara <- lhs_samples[,c(param_constant$parameter)]

y <- matrix(NA, nrow = POC_AU$npops, ncol = length(names(disease_progress)))
colnames(y) <- names(disease_progress)

Param_disease_progress <- rep(list(y), POC_AU$numberSamples) 

z <- matrix(NA, nrow = POC_AU$npops, ncol = length(names(fib)))

colnames(z) <- names(fib)

Param_fib <- rep(list(z), POC_AU$numberSamples) 

disease_progparam <- list()
cured_progparam <- list()

for(i in 1:POC_AU$numberSamples){ 
  
  x[[i]] <- matrix(NA, nrow = 1, ncol = nrow(param_constant))
  x[[i]] <- randompara[i,] *(param_constant$upper - param_constant$lower) + param_constant$lower
  
  disease_progparam[[i]] <- do.call(rbind, replicate(POC_AU$npops, 
                                                     lhs_samples[i,c(names(disease_progress))], simplify=FALSE))
  
  cured_progparam[[i]] <- do.call(rbind, replicate(POC_AU$npops, 
                                                   lhs_samples[i,c(names(fib))], simplify=FALSE))
  
  Param_disease_progress[[i]] <- disease_progparam[[i]]*(disease_progressUL - disease_progressLL) + disease_progressLL

  Param_fib[[i]] <- cured_progparam[[i]]*(cured_progressUL - cured_progressLL) + cured_progressLL
  }

Param_Cparam <- lapply(x, function(m) as.data.frame(t(m)))
Param_estimates <- lapply(Param_Cparam, function(x) x[rep(1,POC_AU$npts), ] )

param_poparray <- rep(list(pop_array), POC_AU$numberSamples) 
pop_array_LL <- pop_array*start_lower
pop_array_UU <- pop_array*start_upper


paramDflist <- rep(list(dfList), POC_AU$numberSamples) 
dflist_LL <- lapply(dfList, function(x) x*start_lower)
dflist_UU <- lapply(dfList, function(x) x*start_upper)
for(i in 1:POC_AU$numberSamples){ 
  param_poparray[[i]] <- lhs_samples[i,"poparray"]* (pop_array_UU - pop_array_LL) + pop_array_LL
  
  param_poparray[[i]][param_poparray[[i]]>1] <- 1 
  
}

for(i in 1:1000){
  paramDflist[[i]][["cured"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["cured"]] - dflist_LL[["cured"]]) + dflist_LL[["cured"]]
  paramDflist[[i]][["eta"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["eta"]] - dflist_LL[["eta"]]) + dflist_LL[["eta"]]
  paramDflist[[i]][["lota"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["lota"]] - dflist_LL[["lota"]]) + dflist_LL[["lota"]]
  paramDflist[[i]][["rho"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["rho"]] - dflist_LL[["rho"]]) + dflist_LL[["rho"]]
  paramDflist[[i]][["tau_ab"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["tau_ab"]] - dflist_LL[["tau_ab"]]) + dflist_LL[["tau_ab"]]
  paramDflist[[i]][["tau_poct"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["tau_poct"]] - dflist_LL[["tau_poct"]]) + dflist_LL[["tau_poct"]]
  paramDflist[[i]][["tau_RNA"]] <- lhs_samples[i,"poparray"]*(dflist_UU[["tau_RNA"]] - dflist_LL[["tau_RNA"]]) + dflist_LL[["tau_RNA"]]
 
    }

save(Param_estimates,
     Param_disease_progress, 
     param_poparray, 
     paramDflist,
     Param_fib,
     lhs_samples,
     file = file.path(OutputFolder,
                      paste0(project_name, "param", ".rda")))



#save(HCV,constants ,disease_progress, fib, dfList, pop_array, 
#     constantsDf, initialPops, best_estimates, best_initial_pop, 
#     file = file.path(projectFolder,
#                      paste0(project_name, ".rda")))


#### saving parameter set as rds####
# check whether value is reasonable and debugging

library(openxlsx)
str(Param_pop_array) 
write.xlsx(Param_disease_progress, 
           file.path(basePath,"01. Data/model input/Param_disease_progress.xlsx"))

write.xlsx(Param_fib, 
           file.path(basePath,"01. Data/model input/Param_fib.xlsx"))
# not working for array 
# write.xlsx(Param_pop_array, 
#          file.path(basePath,"01. Data/model input/Param_pop_array.xlsx"))

write.xlsx(Param_estimates, 
           file.path(basePath,"01. Data/model input/Param_estimates.xlsx"))

write.xlsx(Param_Pops , 
           file.path(basePath,"01. Data/model input/Param_Pops.xlsx"))
