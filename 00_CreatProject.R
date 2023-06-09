# Creating Project

# Re-start R and set to project directory - load packages
rm(list = ls())
library(dplyr)
library(tidyverse)
library(deSolve)
library(ggplot2)
library(reshape)
library(lubridate)

# set up working directory
basePath <- getwd()

# Split equations into function and this is a script to run the functions.
project_name <- "Taiwan-MSM-HCV-model"
project_directory <- file.path(basePath)

###########################################################################
# start of the first year
startYear <-1
  
# start of the final year
endYear <-100

nyears <- endYear - startYear + 1  

# years
Time_Y <- seq(startYear, endYear, by = 1) 

# simulation timestep 
timestep <- 0.25


# population

full_pop_names <- c("HIV neg MSM", 
                    "HIV neg MSM on PrEP", 
                    "HIV pos undiag MSM", 
                    "HIV pos diag MSM")

pop_names <- c("HIV-", 
               "HIV-PrEP", 
               "HIV+", 
               "HIV+d")

npops <- length(pop_names)
###############################################################################
# full progress name 
full_progress_name <- c("suspectible", 
                        "Acute", 
                        "F0", 
                        "F1", 
                        "F2", 
                        "F3", 
                        "F4", 
                        "decompensated cirrhosis", 
                        "hepatocellular carcinoma",
                        "liver treatment", 
                        "post-liver treatment")

short_progress_name <- c("s", 
                         "a", 
                         "f0", 
                         "f1", 
                         "f2", 
                         "f3", 
                         "f4", 
                         "dc", 
                         "hcc", 
                         "lt", 
                         "plt")
# cascade
# cured: {"s0ab+", "s1ab+", "s2ab+", "s3ab+", "s4ab+", "sdcab+", "shccab+"}
cascade_name <- c("Undiagosed", 
                  "diagnosed_ab+", 
                  "diagnosed_ag/RNA", 
                  "treatment", 
                  "Treatment failed",
                  "cured") # without S

short_cascade_name <- c("undiag", 
                        "diag_ab", 
                        "diag_RNA", 
                        "treat", 
                        "treat_f", 
                        "cured")

n.cascade <- length(cascade_name)

n.progress <- length(short_progress_name)

# careful: drop "s" in progress
# overall component: 61
compon <- matrix(NA, nrow = n.cascade, ncol = n.progress-1)

compon_name<- c("s", unlist(lapply(short_progress_name[-1],
                                     FUN = function(x) { 
                                       paste0(x, sep = "_", 
                                              short_cascade_name) })))
compon_nameWOS <- compon_name[-1]


# fib name 
fibName <- c("lt_plt", 
             "f3_cured_f4_cured", 
             "f3_cured_hcc_cured", 
             "f4_cured_dc_cured", 
             "f4_cured_hcc_cured", 
             "dc_cured_hcc_cured", 
             "dc_cured_lt_cured", 
             "hcc_cured_lt_cured", 
             "lt_cured_plt_cured")



transitionName <- c("a_f0", 
                    "f0_f1", 
                    "f1_f2", 
                    "f2_f3", 
                    "f3_f4", 
                    "f3_hcc", 
                    "f4_dc", 
                    "f4_hcc", 
                    "dc_hcc", 
                    "dc_lt", 
                    "hcc_lt")
###############################################################################
# parameters

## constant parameters and adjust number of parameters 
cparam_all <- c("entry", 
                "mortality_b", 
                "MSM_pro",
                "HIV_undiag_diseaseprog",
                "HIV_diag_diseaseprog",  
                "mordc",
                "morhcc", 
                "oddsOfmordc_HIVp",
                "oddsOfmorhcc_HIVp", 
                "morlt", 
                "morplt",
                "Cure_mordc_Reduction", 
                "Cure_morhcc_Reduction")

## constant parameters for each population
cparam_pop <-  c("beta", "spc")

## number of columns of constants parameters 
numConstants <- npops * length(cparam_pop) + length(cparam_all)

# constant parameters data fram generation
constants <- tibble(parameter = character(numConstants),
                        value = 0, lower = 1, upper = 1)

# reorganized columns: cparam_pop first, order like: a1,a2,a3,a4....
# them cparam_all
for (param in 1:length(cparam_pop)) {
  popParams <- paste0(cparam_pop[param], as.character(1:npops))
  index <- (param - 1) * npops + 1
  constants$parameter[index:(index+npops-1)] <- popParams
}
constants$parameter[(index + npops):numConstants] <- cparam_all


# initial population parameters from other parameters
initialPopsNames <- c("init_pop", compon_name, "pop_prop")

numinitialPop <- npops * (length(initialPopsNames)-1) + 1 # init_pop

initial_pop <- tibble(parameter = character(numinitialPop),
                          value = 0, lower = 1, upper = 1)
dim(initial_pop)
fullInitialNames <- character()
for (param in 1:length(initialPopsNames)) {
  if (initialPopsNames[param] %in% c("init_pop")) {
    fullInitialNames <- c(fullInitialNames, initialPopsNames[param])
  } 
  else {
    initialParams <- paste0(initialPopsNames[param],
                            as.character(1:npops))
    fullInitialNames <- c(fullInitialNames, initialParams)
  }
  
}
initial_pop$parameter <- fullInitialNames

# Write constants to file
path <- paste0(project_directory, "/01. DATA/model input")

write.csv(initial_pop, file.path(path, "initial_populations.csv"))

write.csv(constants, file.path(path, "parameters_constants.csv"))


# disease progress rate file 


transition <- as.data.frame(matrix(0, nrow = npops, 
                                   ncol = length(transitionName)))%>%
  setNames(transitionName)


write.csv(transition, file.path(path, "diseaseProgress.csv"))

## transition parameters part II
## these disease progress parameter are constant overtime and populations


fib <- as.data.frame(matrix(0, nrow = npops, ncol = length(fibName)))%>%
  setNames(fibName)

write.csv(fib, file.path(path, "curedProgress.csv"))


# parameters that varies over disease stage and populations
# tau_ab, tau_ag, tau_poct, eta, lota, rho, cured
parameter_variedstage_set <- c("tau_ab", "tau_ag", "tau_poct", "eta", "lota",
                               "rho", "cured") 

param_frame <- as.data.frame(matrix(0, nrow = npops, 
                                    ncol = length(short_progress_name)-1))%>%
  setNames(short_progress_name[-1])

outlist <- list(param_frame)

outlist <- append(outlist, rep(list(param_frame), 
                               length(parameter_variedstage_set)-1))


names(outlist) <- parameter_variedstage_set

# write to seperate csv file 
for(i in names(outlist)){
  write.csv(outlist[[i]], file.path(path, paste0(i,".csv")))
}



# population transition: time-varing
# create a dataframe containing all combination : 4*4 (col), nyears (rows)
transitions <- as.data.frame(matrix(0, nrow = nyears, ncol = npops*npops ))

colnames(transitions) <- c(unlist(lapply(as.list(pop_names), 
                              function(x) paste0(x, pop_names))))


write.csv(transitions, file.path(path, "population_transitions.csv"))


# Create project specifications list

HCV <-list()
HCV$project_name <- project_name
HCV$project_directory <- project_directory

# population
HCV$npops <- npops
HCV$popNames <- pop_names

# component
HCV$component_name <- compon_name
HCV$ncomponent <- length(compon_name)

# component without s
HCV$component_nameWOS <- compon_nameWOS
HCV$ncomponentWOS <- length(compon_nameWOS)

# disease progress
HCV$progress_name <- short_progress_name[-1]
HCV$nprogress <- length(short_progress_name)-1

# disease progression rate 
HCV$diseaseprogress_Name <- transitionName
HCV$diseaseprogress_n <- length(transitionName)

# disease progression rate that are constant overtime and populations
HCV$fibName <- fibName
HCV$fib_n <-length(fibName)

# cascade
HCV$cascade_name <- short_cascade_name
HCV$ncascade <- length(short_cascade_name)
# time related
HCV$startYear <- startYear
HCV$endYear <- endYear
HCV$timestep <- timestep
HCV$nyears <- nyears
HCV$years <- startYear:endYear
  
# Create project .rda files    
save(HCV, file = file.path(project_directory,
                                    paste0(project_name, ".rda")))  



#### cost dataframe #### 

cparameter_cascade <- c("ctau_ab", "ctau_ag", "ctau_poct", "ceta", "clota",
                        "crho", "ccured") 

param_frame <- as.data.frame(matrix(0, nrow = HCV$npops, 
                                    ncol = length(cparameter_cascade)))%>%
  setNames(cparameter_cascade)

write.csv(param_frame, file.path(paste(path, 
                                    "/cost", sep = ""), "costFlow.csv"))

# write to seperate csv file 
# cascade cost
for(i in names(outlist)){
  write.csv(outlist[[i]], file.path(paste(path, 
                                          "/cost", sep = ""), paste0(i,".csv")))
}
# cost of each state 
costPops <- replace(best_initial_pop, best_initial_pop!=0,0)
write.csv(costPops, file.path(paste(path, 
                                    "/cost", sep = ""), "costPops.csv"))
# QALY of each state 
write.csv(costPops, file.path(paste(path, 
                                    "/cost", sep = ""), "QALYPops.csv"))
