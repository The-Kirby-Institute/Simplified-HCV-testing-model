# Set up TWPrisoners model

rm(list = ls())
# Load useful libraries

library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")
library("here")

# Setup directories after setting working directory to source file 
# directory (so it is easy to move code around)
basePath <- here()


# Setup directories after setting working directory to source file 
# directory 

Rcode <- file.path(here()%>%dirname()%>%dirname(),  "03. Code/Functions")
DataFolder <- file.path(basePath, "01. DATA/model input" )
OutputFolder <- file.path(basePath, "02. Output/RDA" )
projectFolder <- file.path(basePath)

project_name <- "TWPrisoners"

# project set-up

load(file.path(OutputFolder, paste0(project_name, ".rda")))

startYear <-1
endYear <- 2000
timestep <- 1/12 
simulateY <- 2022 
calibrateY <- 2015
nyears <- length(seq(startYear,endYear,1))
numberSamples <- 1000

# timeframe 
HCV$startYear <- startYear
HCV$endYear <- endYear
HCV$timestep <- timestep
HCV$nyears <- nyears
HCV$years <- startYear:endYear

HCV$pts <- seq(HCV$startYear, HCV$endYear, by = timestep)

HCV$pts <- head(HCV$pts, -1)
HCV$npts <- length(HCV$pts)
dimNames <- list(HCV$popNames, HCV$componentName)
npops <- HCV$npops
HCV$simY <- simulateY 
HCV$cabY <- calibrateY 
HCV$numberSamples <- numberSamples


# disease progression 
disease_prog <- read.csv(file.path(DataFolder, "diseaseProgress.csv"), 
                         header = TRUE)
disease_progress <- disease_prog[, -1]

cure_prog <- read.csv(file.path(DataFolder, "curedprogress.csv"), 
                      header = TRUE)
fib <- cure_prog[, -1]




# generate () array for parameter_varied over stages
parameter_variedstage_set <- c("tau_ab", 
                               "tau_ag", 
                               "tau_poct", 
                               "eta", 
                               "lota",
                               "rho", 
                               "cured" )

testing_matrix <- array(0, c(npops, HCV$nprogress, HCV$npts + 1), 
                        dimnames = dimNames)

testinglist <- list(testing_matrix)

testinglist <- append(testinglist, rep(list(testing_matrix), 
                                       length(parameter_variedstage_set) - 1))

names(testinglist) <- parameter_variedstage_set

# read all parameter_varied over stages csv files 

files <- list.files(path = paste0(DataFolder, 
                                  "/parameter_varied_stages/", sep = ""),
                    pattern = '*.csv')


dfList <- lapply(files, function(f) {
  df <- read.csv(file.path(paste0(DataFolder, 
                                  "/parameter_varied_stages/", f, sep = "")))
  df <- df[, -1]
  df <- df%>%as_tibble()#%>%mutate_all(funs(1-exp(-.))) # probability to rate 
  df <- as.matrix(df, nrow = npops, ncol = HCV$nprogress)
  df <- replicate(HCV$npts, df)
})

names(dfList) <- c(gsub("^|.csv", "", files)) # ^: from beginning, \ end before .csv

# population trnasiiton array

## import format: long, each column represents the rate of population transition
## data wrangling a bit: drop row number, make first row as column name and 
## drop first column
## Then transpose the data: each row represent the rate of population transition
pop_transitions <- read.csv(file.path(DataFolder, "population_transitions.csv"), 
                            header = FALSE)
str(pop_transitions)
pop_transitions <- pop_transitions[ , -1] # drop the row number

colnames(pop_transitions) <- pop_transitions[1, ]
str(pop_transitions)
pop_transitions <- pop_transitions[-1, ]

# make sure the database is the numeric 
pop_transitions <- mutate_all(pop_transitions, 
                              function(x) as.numeric(as.character(x)))

# expanding pop_transitions to fit the simulation years 
pop_tran <- as.data.frame(matrix(0, nrow = HCV$nyears, 
                                 ncol = ncol(pop_transitions)))

colnames(pop_tran) <- colnames(pop_transitions)

pop_tran[c(1: nrow(pop_transitions)), ] <- pop_transitions  

pop_tran[seq(nrow(pop_transitions) + 1, nrow(pop_tran), 1), ] <- 
  pop_transitions[nrow(pop_transitions),] 



# expand to all points by liner extrapolation between years
pop_extrapolation <- as.data.frame(matrix(0, nrow = HCV$npts,
                                          ncol = ncol(pop_tran)))

colnames(pop_extrapolation) <- colnames(pop_tran)
## assumed pop_array is constant value


for (var in colnames(pop_extrapolation)) {
  
  tempValues <- pop_tran[, var] 
  
  for (year in 1:(HCV$nyears - 1)) {
    
    indices <- ((year - 1) * 1/HCV$timestep + 1): (year * 1/HCV$timestep + 1 )
    
    yearValues <- seq(tempValues[year], tempValues[year + 1], 
                      length = (1/HCV$timestep + 1))
    
    pop_extrapolation[indices, var] <- yearValues
  }
  
}


x.y <- lapply(as.list(1:dim(pop_extrapolation)[1]), 
              function(x) pop_extrapolation[x[1],])

xx <- array(unlist(x.y), c(npops, npops, HCV$npts))
pop_array <- aperm(xx, c(2, 1, 3)) # array transpose # c(2,1,3) 
# means the subscript permutation vector, 
# which must be a permutation of the integers 1:n, 
# where n is the number of dimensions of a

param_constant <- read.csv(file.path(DataFolder, "parameters_constants.csv"), 
                           header = TRUE)%>%select(-"X")%>%as.data.frame()

parameters <- param_constant$parameter 
# reshape 

param_C <- param_constant%>%gather(., "estimate", "value", -"parameter")



constants <- param_constant %>%
  select(parameter, value)
constants <- as.data.frame(t(constants[, -1]))
colnames(constants) <- parameters
rownames(constants) <- NULL

constantsDf <- constants[rep(1, HCV$npts), ]

best_estimates <- as.data.frame(matrix(0, ncol = ncol(constantsDf),
                                       nrow = HCV$npts))

colnames(best_estimates) <- colnames(constantsDf)

for (var in colnames(constantsDf)) {
  
  tempValues <- constantsDf[, var]
  
  for (year in 1:(HCV$nyears - 1)) {
    
    indices <- ((year - 1) * 1/HCV$timestep + 1): (year * 1/HCV$timestep + 1)
    
    yearValues <- seq(tempValues[year], tempValues[year + 1], 
                      length = (1 + 1/HCV$timestep))
    
    best_estimates[indices, var] <- yearValues
  }  
}
best_estimates

## initial population
initialPops<- read.csv(file.path(DataFolder, "initial_populations.csv"), 
                       header = TRUE)%>%select(-"X")

init_pop <- filter(initialPops, parameter == "init_pop")$value 
initial_pop <- initialPops %>%
  filter(parameter != "init_pop")

initialPopMat <- as.data.frame(matrix(initial_pop$value, 
                                      ncol = HCV$ncomponent + 1,  
                                      nrow = HCV$npops))

colnames(initialPopMat) <- c(HCV$component_name, "pop_prop")

rownames(initialPopMat) <- HCV$population_names

popProp <- as.numeric(init_pop) * initialPopMat$pop_prop

best_initial_pop <- apply(initialPopMat[, 1:length(HCV$component_name)], 2, 
                          function(x) x * popProp)



######
save(HCV,constants ,disease_progress, fib, dfList, pop_array, 
     constantsDf, initialPops, best_estimates, best_initial_pop, param_constant, 
     file = file.path(OutputFolder,
                      paste0(project_name, ".rda")))


load(file.path(OutputFolder, paste0(project_name, ".rda")))
source(file.path(Rcode, "/HCV_model.R"))

tic <- proc.time()

steady <- HCVMSM(HCV, best_estimates, best_initial_pop,
                 disease_progress, pop_array, dfList, fib, end_Y = NULL, 
                 modelrun = "steady")

toc <- proc.time() - tic 


df_list <- lapply(steady, as.data.frame.table)

popPro_extract <- df_list$allPops%>%
  mutate(time = rep(seq(HCV$startYear, (HCV$endYear - HCV$timestep), HCV$timestep), 
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



write.csv(popPro_extract, 
          file.path(here(),"01. Data/model input/Estimate_initial_pop.csv")) 

#### number of people in each population ####

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