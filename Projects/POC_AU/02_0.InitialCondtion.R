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

load(file.path(data_path, paste0(project_name, ".rda")))


urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
runSamples <- FALSE
saveAsBase <- TRUE  # if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing

source(file.path(Rcode, "/Functions/HCV_model.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 

load(file.path(DataFolder, paste0(project_name, ".rda")))

# extend to same length 

#### finding steady stage ####

tic <- proc.time()

steady <- HCVMSM(POC_AU, best_estimates, best_initial_pop,
                 disease_progress, pop_array, dfList, fib,
                 modelrun = "steady", proj = "POC_AU")

toc <- proc.time() - tic 

toc
#### extract infected proportion as the infected population allocation #### 

df_list <- lapply(steady, as.data.frame.table)


allpop <- df_list$allPops%>%mutate(time = rep(seq(1,(POC_AU$endYear - POC_AU$timestep),POC_AU$timestep), 
                                              each=POC_AU$ncomponent*POC_AU$npops),
                                   Frequency=round(Freq, digits = 3))
allpop_cascade <- allpop%>%mutate(cascade_status = sub("^[^_]*_", "", Var2), 
                                  dis_prog = sub("\\_.*", "", Var2))%>%
  group_by(Var1, dis_prog, cascade_status, time)%>%
  summarise(statenum = sum(Frequency))

state_change <- allpop_cascade%>%
  mutate(state = ifelse(cascade_status=="cured", "cured",dis_prog))%>%
  group_by(Var1, state, time)%>%summarise(state_num = sum(statenum))%>%
  filter(time<=100)%>%
  ggplot(data = ., aes(x = time, y = state_num )) + 
  geom_line(aes(colour = state)) + 
  xlab("Year") + ylab("number") +labs(tag = "") + 
  scale_x_continuous(expand = c(0,0), limits = c(0, 100), 
                     breaks = seq(0, 100, 10))  + 
  facet_wrap(.~Var1, ncol=2, scales = "free") + 
  plotOpts 

state_change

allpop_cascade%>%
  mutate(state = ifelse(cascade_status=="cured", "cured",dis_prog))%>%
  mutate(CP = ifelse(Var1%in%c("C_PWID", "C_fPWID"), "C",
                     ifelse(Var1%in%c("P_PWID", "P_fPWID"), "P", "P_nonPWID")))%>%
  group_by(CP, state, time)%>%summarise(state_num = sum(statenum))%>%
  filter(time<=1500)%>%
  ggplot(data = ., aes(x = time, y = state_num )) + 
  geom_line(aes(colour = state)) + 
  xlab("Year") + ylab("number") +labs(tag = "") + 
  scale_x_continuous(expand = c(0,0), limits = c(0, 1500), 
                     breaks = seq(0, 1500, 100))  + 
  facet_wrap(.~CP, ncol=2, scales = "free") 




allpop_cascade%>%
  mutate(state = ifelse(cascade_status=="cured", "cured",dis_prog))%>%
  mutate(CP = ifelse(Var1%in%c("C_PWID", "C_fPWID"), "C",
                     ifelse(Var1%in%c("P_PWID", "P_fPWID"), "P", "P_nonPWID")))%>%
  group_by(state, time)%>%summarise(state_num = sum(statenum))%>%
  filter(time<=1500)%>%
  ggplot(data = ., aes(x = time, y = state_num )) + 
  geom_line(aes(colour = state)) + 
  xlab("Year") + ylab("number") +labs(tag = "") + 
  scale_x_continuous(expand = c(0,0), limits = c(0, 1500), 
                     breaks = seq(0, 1500, 100))  








popPro_extract <- df_list$allPops%>%
  mutate(time = rep(seq(POC_AU$startYear, (POC_AU$endYear - POC_AU$timestep), POC_AU$timestep), 
                    each=POC_AU$ncomponent * POC_AU$npops),
         Frequency=round(Freq, digits = 3))%>%
  filter(time==250)%>%
  mutate(cascade_status = sub("^[^_]*_", "", Var2), 
         dis_prog = sub("\\_.*", "", Var2),
         SI = ifelse(cascade_status%in%c("s", "cured"), "S","I"),
         parameter =Var2)%>%group_by(Var1 ,SI)%>%
  mutate(total = sum(Frequency),
         value = ifelse(Frequency==0, 0, round(Frequency/total, digits = 4)))%>%
  ungroup()%>%group_by(Var1)%>%mutate(pop_prop = ifelse(
    Frequency==0, 0, round(Frequency/sum(Frequency), digits = 4)))%>%
  ungroup()%>%select(Var1,parameter, value, SI)

ggplot(data = popPro_extract) + 
  geom_line()

write.csv(popPro_extract, 
          file.path(basePath,"01. Data/model input/Estimate_initial_pop.csv")) 

#### number of MSM in each population ####

estPops <- read.csv(file.path(DataFolder, "Estimate_initial_pop.csv"), 
                   header = TRUE)%>%select(-"X")

estPops <- popPro_extract 

init_pop <- filter(initialPops, parameter == "init_pop")$value

pop_prop <- initialPops%>%filter(parameter%in% c("pop_prop1", "pop_prop2", 
                                                 "pop_prop3", "pop_prop4", "pop_prop5"))%>%
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
                     Var1 == "P_PWID+" & SI == "S" ~ init_prop_S[3],
                     Var1 == "P_PWID+" & SI == "I" ~ init_prop_I[3],
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

save(HCV,steady, best_est_pop, 
     file = file.path(projectFolder,
                      paste0(project_name,"cali" ,".rda")))


