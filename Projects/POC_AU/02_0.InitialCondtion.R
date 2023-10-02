# Running 02_1.CalibrateOptimze.R and revised the input parameters then 
# Run 01.SetupModel to align the input parameters to fit the format that model need 
# Then Run this script to get the initial condition 

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
# load(file.path(OutputFolder, paste0(project_name, "cali.rda")))


urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
runSamples <- FALSE
saveAsBase <- TRUE  # if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing

source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R"))


# extend to same length 

#### finding steady stage ####

tic <- proc.time()

steady <- HCVMSM(POC_AU, best_estimates, best_initial_pop,
                 disease_progress, pop_array, dfList, fib,
                 modelrun = "steady", proj = "POC_AU")

toc <- proc.time() - tic 

toc
#### extract infected proportion as the infected population allocation #### 
# finding equilibrium

check_steady(model_result = steady, endY = POC_AU$endYear,
             timestep = POC_AU$timestep, 
             Ncomp = POC_AU$ncomponent*POC_AU$npops, 
             Tequilibrium = 800)


#extract proportion of pops for initial condition 

df_list <- lapply(steady, as.data.frame.table)

allpop <- df_list$allPops%>%mutate(time = rep(seq(1,(POC_AU$endYear - POC_AU$timestep),POC_AU$timestep), 
                                              each=POC_AU$ncomponent*POC_AU$npops),
                                   Frequency=Freq)

popPro_extract <- df_list$allPops%>%
  mutate(time = rep(seq(POC_AU$startYear, (POC_AU$endYear - POC_AU$timestep), POC_AU$timestep), 
                    each=POC_AU$ncomponent * POC_AU$npops),
         Frequency=Freq)%>%
  filter(time == 500)%>%
  mutate(cascade_status = sub("^[^_]*_", "", Var2), 
         dis_prog = sub("\\_.*", "", Var2),
         SI = ifelse(cascade_status%in%c("s", "cured"), "S","I"),
         parameter =Var2)%>%group_by(Var1 ,SI)%>%
  mutate(total = sum(Frequency),
         value = ifelse(Frequency==0, 0, Frequency/total))%>%
  ungroup()%>%group_by(Var1)%>%mutate(pop_prop = ifelse(
    Frequency==0, 0, Frequency/sum(Frequency)))%>%
  ungroup()%>%select(Var1,parameter, value, SI)


write.csv(popPro_extract, 
          file.path(DataFolder,"/Estimate_initial_pop.csv")) 

#### number of PWID/former PWID in each population ####

estPops <- read.csv(file.path(DataFolder, "Estimate_initial_pop.csv"), 
                   header = TRUE)%>%select(-"X")


init_pop <- filter(initialPops, parameter == "init_pop")$value

pop_prop <- initialPops%>%filter(parameter%in% c("pop_prop1", "pop_prop2", 
                                                 "pop_prop3", "pop_prop4", 
                                                 "pop_prop5"))%>%
  select(value)%>%unlist()%>%as.vector()

popProp <- as.numeric(init_pop)*pop_prop 


# prevalence at initial
init_prop_I <- c(constantsDf$HCVP1[1], constantsDf$HCVP2[1], 
                 constantsDf$HCVP3[1], constantsDf$HCVP4[1], constantsDf$HCVP5[1])

init_prop_S <-c(1 - init_prop_I)

estPops <- estPops%>%mutate(
  pop_group = rep(c(popProp),dim(estPops)[1]/POC_AU$npops),
  SIprop = case_when(Var1 == "C_PWID" & SI == "S" ~ init_prop_S[1],
                     Var1 == "C_PWID" & SI == "I" ~ init_prop_I[1],
                     Var1 == "C_fPWID" & SI == "S" ~ init_prop_S[2],
                     Var1 == "C_fPWID" & SI == "I" ~ init_prop_I[2],
                     Var1 == "P_PWID" & SI == "S" ~ init_prop_S[3],
                     Var1 == "P_PWID" & SI == "I" ~ init_prop_I[3],
                     Var1 == "P_fPWID" & SI == "S" ~ init_prop_S[4],
                     Var1 == "P_fPWID" & SI == "I" ~ init_prop_I[4], 
                     Var1 == "P_nPWID" & SI == "S" ~ init_prop_S[5],
                     Var1 == "P_nPWID" & SI == "I" ~ init_prop_I[5]
                     ),
  est_pop = value*pop_group*SIprop)

best_est_pop <- as.matrix(as.data.frame(matrix(estPops$est_pop, 
                                               ncol = POC_AU$ncomponent,  
                                               nrow = POC_AU$npops)))

colnames(best_est_pop) <- c(POC_AU$component_name)

save(project_name,steady, best_est_pop, 
     file = file.path(OutputFolder ,
                      paste0(project_name,"cali" ,".rda")))


# calibration 
tic <- proc.time()
endY <- 36
calibrateInit <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                       disease_progress,  pop_array, dfList, fib,
                       modelrun="UN", proj = "POC_AU", end_Y = endY)


toc <- proc.time() - tic 
toc




save(calibrateInit, 
     file = file.path(OutputFolder ,
                      paste0(project_name,"cali_init" ,".rda")))
                                
                                     
                                        