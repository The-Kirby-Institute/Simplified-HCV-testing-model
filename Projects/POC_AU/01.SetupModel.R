rm(list = ls())

# Setup directories after setting working directory to source file 
# directory (so it is easy to move code around)
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

POC_AU$endYear

POC_AU$dimName <- list(POC_AU$popNames, POC_AU$progress_name)

#### disease progression ####
disease_prog <- read.csv(file.path(DataFolder, "diseaseProgress.csv"), 
                         header = TRUE)
disease_progress <- disease_prog[, -1]

rownames(disease_progress) <- POC_AU$popNames

cure_prog <- read.csv(file.path(DataFolder, "curedprogress.csv"), 
                      header = TRUE)
fib <- cure_prog[, -1]

rownames(fib) <- POC_AU$popNames



#### generate () array for parameter_varied over stages ####
parameter_variedstage_set <- c("tau_ab", 
                               "tau_RNA", 
                               "tau_poct", 
                               "eta", 
                               "lota",
                               "rho", 
                               "cured" )

testing_matrix <- array(0, c(POC_AU$npops, POC_AU$nprogress, POC_AU$npts + 1), 
                        dimnames = POC_AU$dimName)

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
  df <- df%>%as_tibble()
  df <- as.matrix(df, nrow = npops, ncol = HCV$nprogress)
  df <- replicate(POC_AU$npts, df)
  
  rownames(df) <- POC_AU$popNames
  return(df)
})

names(dfList) <- c(gsub("^|.csv", "", files)) # ^: from beginning, \ end before .csv


#### population trnasition array ####

## import format: long, each column represents the rate of population transition
## data wrangling a bit: drop row number, make first row as column name and 
## drop first column
## Then transpose the data: each row represent the rate of population transition
pop_transitions <- read.csv(file.path(DataFolder, "population_transitions.csv"), 
                            header = FALSE)

pop_transitions <- pop_transitions[ , -1] # drop the row number

colnames(pop_transitions) <- pop_transitions[1, ]

pop_transitions <- pop_transitions[-1, ]

# make sure the database is the numeric 
pop_transitions <- mutate_all(pop_transitions, 
                              function(x) as.numeric(as.character(x)))

# expanding pop_transitions to fit the simulation years 
pop_tran <- as.data.frame(matrix(0, nrow = POC_AU$nyears, 
                                 ncol = ncol(pop_transitions)))

colnames(pop_tran) <- colnames(pop_transitions)

pop_tran[c(1: nrow(pop_transitions)), ] <- pop_transitions  

pop_tran[seq(nrow(pop_transitions) + 1, nrow(pop_tran), 1), ] <- 
  pop_transitions[nrow(pop_transitions),] 



# expand to all points by liner extrapolation between years
pop_extrapolation <- as.data.frame(matrix(0, nrow = POC_AU$npts,
                                          ncol = ncol(pop_tran)))

colnames(pop_extrapolation) <- colnames(pop_tran)
## assumed pop_array is constant value


for (var in colnames(pop_extrapolation)) {
  
  tempValues <- pop_tran[, var] 
  
  for (year in 1:(POC_AU$nyears - 1)) {
    
    indices <- ((year - 1) * 1/POC_AU$timestep + 1): (year * 1/POC_AU$timestep + 1 )
    
    yearValues <- seq(tempValues[year], tempValues[year + 1], 
                      length = (1/POC_AU$timestep + 1))
    
    pop_extrapolation[indices, var] <- yearValues
  }
  
}


x.y <- lapply(as.list(1:dim(pop_extrapolation)[1]), 
              function(x) pop_extrapolation[x[1],])

xx <- array(unlist(x.y), c(POC_AU$npops, POC_AU$npops, POC_AU$npts))
pop_array <- aperm(xx, c(2, 1, 3)) # array transpose # c(2,1,3) 
# means the subscript permutation vector, 
# which must be a permutation of the integers 1:n, 
# where n is the number of dimensions of a

dimnames(pop_array) <- list(POC_AU$popNames, POC_AU$popNames) 


#### other parameters ####
param_constant <- read.csv(file.path(DataFolder, "parameters_constants.csv"), 
                           header = TRUE)%>%select(-"X")%>%as.data.frame()

parameters <- param_constant$parameter 

param_C <- param_constant%>%gather(., "estimate", "value", -"parameter")

constants <- param_constant %>%
  select(parameter, value)

constants <- as.data.frame(t(constants[, -1]))
colnames(constants) <- parameters
rownames(constants) <- NULL

# entry = male entry* MSM prop
constants$entry <- constants$entry*constants$MSM_pro

constantsDf <- constants[rep(1, POC_AU$npts), ]

best_estimates <- as.data.frame(matrix(0, ncol = ncol(constantsDf),
                                       nrow = POC_AU$npts))

colnames(best_estimates) <- colnames(constantsDf)

for (var in colnames(constantsDf)) {
  
  tempValues <- constantsDf[, var]
  
  for (year in 1:(POC_AU$nyears - 1)) {
    
    indices <- ((year - 1) * 1/POC_AU$timestep + 1): (year * 1/POC_AU$timestep + 1)
    
    yearValues <- seq(tempValues[year], tempValues[year + 1], 
                      length = (1 + 1/POC_AU$timestep))
    
    best_estimates[indices, var] <- yearValues
  }  
}


#### initial population ####
initialPops <- read.csv(file.path(DataFolder, "initial_populations.csv"), 
                       header = TRUE)%>%select(-"X")


#init_pop <- filter(initialPops, parameter == "init_pop")$value * constants$MSM_pro 

initial_pop <- initialPops %>%
  filter(parameter != "init_pop")

init_pop <- initialPops %>%
  filter(parameter == "init_pop")

initialPopMat <- as.data.frame(matrix(initial_pop$value, 
                                      ncol = POC_AU$ncomponent + 1,  
                                      nrow = POC_AU$npops))

colnames(initialPopMat) <- c(POC_AU$component_name, "pop_prop")

rownames(initialPopMat) <- POC_AU$population_names

popProp <- as.numeric(init_pop$value) * initialPopMat$pop_prop

best_initial_pop <- apply(initialPopMat[, 1:length(POC_AU$component_name)], 2, 
                          function(x) x * popProp)



######
save(POC_AU,constants ,disease_progress, fib, dfList, pop_array, 
     constantsDf, initialPops, best_estimates, best_initial_pop, param_constant, 
     file = file.path(DataFolder,
                      paste0(project_name, ".rda")))
