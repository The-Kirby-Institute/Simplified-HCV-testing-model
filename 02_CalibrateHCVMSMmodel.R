# expected to be a function to update the model parameters, and rerun the HCVMSM 
# model for the best estimate parameters. Then output the model plot. 
rm(list = ls())

library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
basePath <- getwd()
projectFolder <- file.path(basePath)
project_name <- "HCVModel"
DataFolder <- file.path(basePath, "01. DATA/model input" )
Rcode <- file.path(basePath, "03. Code")
urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
runSamples <- TRUE
saveAsBase <- TRUE  # if TRUE doesn't append time to results and overwrites
                    # a base file. Useful for storing main results or
                    # testing

source(file.path(Rcode, "/Functions/HCV_model.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 

load(file.path(projectFolder, paste0(project_name, ".rda")))

# project set-up



# extend to same length 

#### finding steady stage ####

tic <- proc.time()

steady <- HCVMSM(HCV, best_estimates, best_initial_pop,
                      disease_progress, pop_array, dfList, fib,
                modelrun = "steady")

toc <- proc.time() - tic 

#### extract infected proportion as the infected population allocation #### 

df_list <- lapply(steady, as.data.frame.table)

popPro_extract <- df_list$allPops%>%
  mutate(time = rep(seq(HCV$startYear, (HCV$endYear - HCV$timestep), 0.1), 
                    each=HCV$ncomponent * HCV$npops),
         Frequency=round(Freq, digits = 3))%>%
  filter(time==800)%>%
  mutate(cascade_status = sub("^[^_]*_", "", Var2), 
         dis_prog = sub("\\_.*", "", Var2),
         SI = ifelse(cascade_status%in%c("s", "cured"), "S","I"),
         parameter =Var2)%>%group_by(Var1 ,SI)%>%
  mutate(total = sum(Frequency),
         value = ifelse(Frequency==0, 0, round(Frequency/total, digits = 4)))%>%
  ungroup()%>%group_by(Var1)%>%mutate(pop_prop = ifelse(
    Frequency==0, 0, round(Frequency/sum(Frequency), digits = 4)))%>%
  ungroup()%>%select(Var1,parameter, value, SI)

basePath <- getwd()

write.csv(popPro_extract, 
          file.path(basePath,"01. Data/model input/Estimate_initial_pop.csv")) 

#### number of MSM in each population ####

estPops<- read.csv(file.path(DataFolder, "Estimate_initial_pop.csv"), 
                   header = TRUE)%>%select(-"X")

init_pop <- filter(initialPops, parameter == "init_pop")$value*constants$MSM_pro

pop_prop <- initialPops%>%filter(parameter%in% c("pop_prop1", "pop_prop2", 
                                                 "pop_prop3", "pop_prop4"))%>%
  select(value)%>%unlist()%>%as.vector()

popProp <- as.numeric(init_pop)*pop_prop 


# prevalence at initial
init_prop_I <- c(constantsDf$HCVP1[1], constantsDf$HCVP2[1], 
                 constantsDf$HCVP3[1], constantsDf$HCVP4[1])

init_prop_S <-c(1 - init_prop_I)

estPops <- estPops%>%mutate(
  pop_group = rep(c(popProp),dim(estPops)[1]/HCV$npops),
  SIprop = case_when(Var1 == "HIV-" & SI == "S" ~ init_prop_S[1],
                       Var1 == "HIV-" & SI == "I" ~ init_prop_I[1],
                       Var1 == "HIV-PrEP" & SI == "S" ~ init_prop_S[2],
                       Var1 == "HIV-PrEP" & SI == "I" ~ init_prop_I[2],
                       Var1 == "HIV+" & SI == "S" ~ init_prop_S[3],
                       Var1 == "HIV+" & SI == "I" ~ init_prop_I[3],
                       Var1 == "HIV+d" & SI == "S" ~ init_prop_S[4],
                       Var1 == "HIV+d" & SI == "I" ~ init_prop_I[4]),
  est_pop = value*pop_group*SIprop)

best_est_pop <- as.matrix(as.data.frame(matrix(estPops$est_pop, 
                                               ncol = HCV$ncomponent,  
                                               nrow = HCV$npops)))

colnames(best_est_pop) <- c(HCV$component_name)

save(HCV,steady, best_est_pop, 
     file = file.path(projectFolder,
                      paste0(project_name,"cali" ,".rda")))
