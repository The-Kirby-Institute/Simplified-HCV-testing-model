#### set up projects
rm(list = ls()) 
# library
library(readxl)
library(dplyr)
library(here)
library(doMC)
library(abind)
library(tidyverse)
library(reshape)
library(lubridate)

# set up dircetory 
basePath <- here()
dir.create("01. DATA") # create subdircetory
dir.create("01. DATA/data input")
dir.create("02. Output") 
dir.create("02. Output/RDA")
project_directory <- here()
DataFolder  <- file.path(here(), "01. DATA/model input")
OutputFolder <- file.path(here(), "02. Output/RDA")
project_name <- "POC_AU"
# source above a layer, if above doulble layer: adding %>%dirname() 
# file.path(here()%>%dirname(), "TWPrisoners")

###########################################################################
# start of the first year
startYear <-1

# start of the final year
endYear <-1000

nyears <- endYear - startYear + 1  

# years
Time_Y <- seq(startYear, endYear, by = 1) 

# simulation timestep 
timestep <- 1/12

# population

full_pop_names <- c("PWID in community",  "Former PWID in community", 
                    "PWID in prisons",  "Former PWID in prisons", 
                    "nonPWID in prisons")

pop_names <- c("C_PWID", 
               "C_fPWID", 
               "P_PWID", 
               "P_fPWID",
               "P_nPWID")

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
                        "liver transplant", 
                        "post-liver transplant")

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
                  "diagnosed_RNA", 
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
                "Cure_mordc_Reduction", 
                "Cure_morhcc_Reduction")

## constant parameters for each population
cparam_pop <-  c("beta", "spc", "mortality_b", "mordc",
                 "morhcc","morlt", "morplt")

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


write.csv(initial_pop, file.path(DataFolder , "initial_populations.csv"))

write.csv(constants, file.path(DataFolder , "parameters_constants.csv"))


# disease progress rate file 


transition <- as.data.frame(matrix(0, nrow = npops, 
                                   ncol = length(transitionName)))%>%
  setNames(transitionName)


write.csv(transition, file.path(DataFolder, "diseaseProgress.csv"))

## transition parameters part II
## these disease progress parameter are constant overtime and populations


fib <- as.data.frame(matrix(0, nrow = npops, ncol = length(fibName)))%>%
  setNames(fibName)

write.csv(fib, file.path(DataFolder, "curedProgress.csv"))


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
  write.csv(outlist[[i]], file.path(DataFolder, paste0(i,".csv")))
}



# population transition: time-varing
# create a dataframe containing all combination : 4*4 (col), nyears (rows)
transitions <- as.data.frame(matrix(0, nrow = nyears, ncol = npops*npops ))

colnames(transitions) <- c(unlist(lapply(as.list(pop_names), 
                                         function(x) paste0(x, pop_names))))


write.csv(transitions, file.path(DataFolder, "population_transitions.csv"))


# Create project specifications list

POC_AU <-list()
POC_AU$project_name <- project_name


# population
POC_AU$npops <- npops
POC_AU$popNames <- pop_names

# component
POC_AU$component_name <- compon_name
POC_AU$ncomponent <- length(compon_name)

# component without s
POC_AU$component_nameWOS <- compon_nameWOS
POC_AU$ncomponentWOS <- length(compon_nameWOS)

# disease progress
POC_AU$progress_name <- short_progress_name[-1]
POC_AU$nprogress <- length(short_progress_name)-1

# disease progression rate 
POC_AU$diseaseprogress_Name <- transitionName
POC_AU$diseaseprogress_n <- length(transitionName)

# disease progression rate that are constant overtime and populations
POC_AU$fibName <- fibName
POC_AU$fib_n <-length(fibName)

# cascade
POC_AU$cascade_name <- short_cascade_name
POC_AU$ncascade <- length(short_cascade_name)
# time related
POC_AU$startYear <- startYear
POC_AU$endYear <- endYear
POC_AU$timestep <- timestep
POC_AU$nyears <- nyears
POC_AU$years <- startYear:endYear

# Create project .rda files    
save(POC_AU, file = file.path(OutputFolder,
                           paste0(project_name, ".rda")))  


