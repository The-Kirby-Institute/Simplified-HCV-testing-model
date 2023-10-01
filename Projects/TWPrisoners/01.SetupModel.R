# Set up TWPrisoners model

rm(list = ls())

# Setup directories after setting working directory to source file 
# directory (so it is easy to move code around)
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
TWPrisoners$startYear <- startYear
TWPrisoners$endYear <- endYear
TWPrisoners$timestep <- timestep
TWPrisoners$nyears <- nyears
TWPrisoners$years <- startYear:endYear

TWPrisoners$pts <- seq(TWPrisoners$startYear, TWPrisoners$endYear, by = timestep)

TWPrisoners$pts <- head(TWPrisoners$pts, -1)
TWPrisoners$npts <- length(TWPrisoners$pts)
dimNames <- list(TWPrisoners$popNames, TWPrisoners$componentName)
npops <- TWPrisoners$npops
TWPrisoners$simY <- simulateY 
TWPrisoners$cabY <- calibrateY 
TWPrisoners$numberSamples <- numberSamples


# disease progression 
disease_prog <- read.csv(file.path(DataFolder, "diseaseProgress.csv"), 
                         header = TRUE)
disease_progress <- disease_prog[, -1]

cure_prog <- read.csv(file.path(DataFolder, "curedprogress.csv"), 
                      header = TRUE)
fib <- cure_prog[, -1]




# generate () array for parameter_varied over stages
parameter_variedstage_set <- c("tau_ab", 
                               "tau_RNA", 
                               "tau_poct", 
                               "eta", 
                               "lota",
                               "rho", 
                               "cured" )

testing_matrix <- array(0, c(npops, TWPrisoners$nprogress, TWPrisoners$npts + 1), 
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
  df <- as.matrix(df, nrow = npops, ncol = TWPrisoners$nprogress)
  df <- replicate(TWPrisoners$npts, df)
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
pop_tran <- as.data.frame(matrix(0, nrow = TWPrisoners$nyears, 
                                 ncol = ncol(pop_transitions)))

colnames(pop_tran) <- colnames(pop_transitions)

pop_tran[c(1: nrow(pop_transitions)), ] <- pop_transitions  

pop_tran[seq(nrow(pop_transitions) + 1, nrow(pop_tran), 1), ] <- 
  pop_transitions[nrow(pop_transitions),] 



# expand to all points by liner extrapolation between years
pop_extrapolation <- as.data.frame(matrix(0, nrow = TWPrisoners$npts,
                                          ncol = ncol(pop_tran)))

colnames(pop_extrapolation) <- colnames(pop_tran)
## assumed pop_array is constant value


for (var in colnames(pop_extrapolation)) {
  
  tempValues <- pop_tran[, var] 
  
  for (year in 1:(TWPrisoners$nyears - 1)) {
    
    indices <- ((year - 1) * 1/TWPrisoners$timestep + 1): (year * 1/TWPrisoners$timestep + 1 )
    
    yearValues <- seq(tempValues[year], tempValues[year + 1], 
                      length = (1/TWPrisoners$timestep + 1))
    
    pop_extrapolation[indices, var] <- yearValues
  }
  
}


x.y <- lapply(as.list(1:dim(pop_extrapolation)[1]), 
              function(x) pop_extrapolation[x[1],])

xx <- array(unlist(x.y), c(npops, npops, TWPrisoners$npts))
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

constantsDf <- constants[rep(1, TWPrisoners$npts), ]

best_estimates <- as.data.frame(matrix(0, ncol = ncol(constantsDf),
                                       nrow = TWPrisoners$npts))

colnames(best_estimates) <- colnames(constantsDf)

for (var in colnames(constantsDf)) {
  
  tempValues <- constantsDf[, var]
  
  for (year in 1:(TWPrisoners$nyears - 1)) {
    
    indices <- ((year - 1) * 1/TWPrisoners$timestep + 1): (year * 1/TWPrisoners$timestep + 1)
    
    yearValues <- seq(tempValues[year], tempValues[year + 1], 
                      length = (1 + 1/TWPrisoners$timestep))
    
    best_estimates[indices, var] <- yearValues
  }  
}
best_estimates

## initial population
initialPops<- read.csv(file.path(DataFolder, "initial_populations.csv"), 
                       header = TRUE)%>%dplyr::select(-"X")

init_pop <- filter(initialPops, parameter == "init_pop")$value 
initial_pop <- initialPops %>%
  filter(parameter != "init_pop")



initialPopMat <- as.data.frame(matrix(initial_pop$value, 
                                      ncol = TWPrisoners$ncomponent + 1,  
                                      nrow = TWPrisoners$npops))

colnames(initialPopMat) <- c(TWPrisoners$component_name, "pop_prop")

rownames(initialPopMat) <- TWPrisoners$population_names

popProp <- as.numeric(init_pop) * initialPopMat$pop_prop

best_initial_pop <- apply(initialPopMat[, 1:length(TWPrisoners$component_name)], 2, 
                          function(x) x * popProp)


save(TWPrisoners,constants ,disease_progress, fib, dfList, pop_array, 
     constantsDf, initialPops, best_estimates, best_initial_pop, param_constant, 
     file = file.path(OutputFolder,
                      paste0(project_name, ".rda")))
