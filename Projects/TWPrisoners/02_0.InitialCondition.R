# Running 02_1.CalibrateOptimze.R and revised the input parameters then 
# Run 01.SetupModel to align the input parameters to fit the format that model need 
# Then Run this script to get the initial condition 

rm(list = ls())

library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)

project_name <- "TWPrisoners"

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
OutputFolder <- file.path(data_path, "02. Output/RDA")

load(file.path(OutputFolder, paste0(project_name, ".rda")))


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

steady <- HCVMSM(TWPrisoners, best_estimates, best_initial_pop,
                 disease_progress, pop_array, dfList, fib, end_Y = NULL, 
                 modelrun = "steady", proj = "TWPrisoners")

toc <- proc.time() - tic 
toc



check_steady(model_result = steady, endY = TWPrisoners$endYear,
             timestep = TWPrisoners$timestep, 
             Ncomp = TWPrisoners$ncomponent*TWPrisoners$npops, 
             Tequilibrium = 1500)


df_list <- lapply(steady, as.data.frame.table)

popPro_extract <- df_list$allPops%>%
  mutate(time = rep(seq(TWPrisoners$startYear, (TWPrisoners$endYear - TWPrisoners$timestep), TWPrisoners$timestep), 
                    each=TWPrisoners$ncomponent * TWPrisoners$npops))%>%
  filter(time==1500)%>%
  mutate(cascade_status = sub("^[^_]*_", "", Var2), 
         dis_prog = sub("\\_.*", "", Var2),
         SI = ifelse(cascade_status%in%c("s", "cured"), "S","I"),
         parameter =Var2)%>%group_by(Var1 ,SI)%>%
  mutate(total = sum(Freq),
         value = ifelse(Freq==0, 0, Freq/total))%>%
  ungroup()%>%group_by(Var1)%>%mutate(pop_prop = ifelse(
    Freq==0, 0, Freq/sum(Freq)))%>%
  ungroup()%>%dplyr::select(Var1,parameter, value, SI)



write.csv(popPro_extract, 
          file.path(DataFolder,"/Estimate_initial_pop.csv")) 

#### number of people in each population ####

estPops<- read.csv(file.path(DataFolder, "Estimate_initial_pop.csv"), 
                   header = TRUE)%>%dplyr::select(-"X")

init_pop <- filter(initialPops, parameter == "init_pop")$value

pop_prop <- initialPops%>%filter(parameter%in% c("pop_prop1", "pop_prop2", 
                                                 "pop_prop3", "pop_prop4", 
                                                 "pop_prop5", "pop_prop6"))%>%
  dplyr::select(value)%>%unlist()%>%as.vector()

sum(pop_prop)
popProp <- as.numeric(init_pop)*pop_prop 

sum(popProp)
# prevalence at initial
init_prop_I <- c(constantsDf$HCVP1[1], constantsDf$HCVP2[1], 
                 constantsDf$HCVP3[1], constantsDf$HCVP4[1],
                 constantsDf$HCVP5[1], constantsDf$HCVP6[1])

init_prop_S <-c(1 - init_prop_I)

estPops <- estPops%>%mutate(
  pop_group = rep(c(popProp),dim(estPops)[1]/TWPrisoners$npops),
  SIprop = case_when(Var1 == "N_inca_NCID" & SI == "S" ~ init_prop_S[1],
                     Var1 == "N_inca_NCID" & SI == "I" ~ init_prop_I[1],
                     Var1 == "N_inca_CID" & SI == "S" ~ init_prop_S[2],
                     Var1 == "N_inca_CID" & SI == "I" ~ init_prop_I[2],
                     Var1 == "E_inca_NCID" & SI == "S" ~ init_prop_S[3],
                     Var1 == "E_inca_NCID" & SI == "I" ~ init_prop_I[3],
                     Var1 == "E_inca_CID" & SI == "S" ~ init_prop_S[4],
                     Var1 == "E_inca_CID" & SI == "I" ~ init_prop_I[4],
                     Var1 == "D_inca" & SI == "S" ~ init_prop_S[5],
                     Var1 == "D_inca" & SI == "I" ~ init_prop_I[5],
                     Var1 == "P_inca" & SI == "S" ~ init_prop_S[6],
                     Var1 == "P_inca" & SI == "I" ~ init_prop_I[6]),
  
  est_pop = value*pop_group*SIprop)

best_est_pop <- as.matrix(as.data.frame(matrix(estPops$est_pop, 
                                               ncol = TWPrisoners$ncomponent,  
                                               nrow = TWPrisoners$npops)))

colnames(best_est_pop) <- c(TWPrisoners$component_name)
sum(best_est_pop)

save(TWPrisoners,steady, best_est_pop, 
     file = file.path(OutputFolder,
                      paste0(project_name,"cali" ,".rda")))




# calibration 
tic <- proc.time()
endY <- 15
calibrateInit <- HCVMSM(TWPrisoners, best_estimates, best_est_pop,
                        disease_progress,  pop_array, dfList, fib,
                        modelrun="UN", proj = "TWPrisoners", end_Y = endY)


toc <- proc.time() - tic 
toc

save(calibrateInit, 
     file = file.path(OutputFolder ,
                      paste0(project_name,"cali_init" ,".rda")))



