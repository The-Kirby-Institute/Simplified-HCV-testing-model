# parameter set
# Number of parameter sets to run
rm(list = ls())
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
load(file.path(projectFolder, paste0(project_name, "cali",".rda")))


number_samples <- 1000
start_lower <- 0.75
start_upper <- 1.25
end_lower <- 0.75
end_upper <- 1.25
sParam <- 250

HCV$numberSamples <- number_samples
HCV$sParam <- sParam


#### parameter set -------------------------------------------------------------

#### list apply 
# constant
randomParams <- function(param, constants) { 
  
  paramValues <- filter(constants, parameter == param)
  
  paramValues_up <-  paramValues[paramValues$estimate == "upper", "value"]
  
  paramValues_low <- paramValues[paramValues$estimate == "lower", "value"]
  
  paramSample <- runif(HCV$numberSamples, paramValues_low, 
                       paramValues_up)
  
  return(paramSample)
  
}


# parameter_varied over stages list of array
randomParams_array <- runif(number_samples, start_lower, start_upper)



# Create parameter set samples --------------------------------------------
set.seed(123456) 


## disease progress has three dataset: best estimate, Upper limit, Lower limt 
# random sampling between UL & LL 
# then generating a list of each set, which is similar to disease progress for 
# best estimate 

##### make a function ##### 

## for disease progress and cured progress 
# random sampling the number 

randomParam_progress <- function(LL, UL, best) {
  ran <- list()
  ran_w <- list()
  for (var in colnames(best)) { 
    ran[[var]] <-cbind(c(runif(HCV$numberSamples, LL[var][[1]][1], 
                               UL[var][[1]][1])),
                     c(runif(HCV$numberSamples, LL[var][[1]][2], 
                             UL[var][[1]][2])),
                     c(runif(HCV$numberSamples, LL[var][[1]][3], 
                             UL[var][[1]][3])),
                     c(runif(HCV$numberSamples, LL[var][[1]][4], 
                             UL[var][[1]][4]))) 
      # wide form: [ ,1] = 1 st sample set; [ ,2]: 2nd sample set  
      ran_w[[var]] <- t(ran[[var]]) 
      }
  m <- list()
  Rdt <- list()
  
  for ( i in 1:HCV$numberSamples) {
    m[[i]] <- data.frame(matrix(unlist(
      lapply(ran_w, function(x)  cbind(x[ ,i]))),
      nrow = HCV$npops, ncol = ncol(best)))   
  }
  
  Rdt <- lapply(m, "colnames<-", c(colnames(best)))
  
  return(Rdt)
}

# disease progress 

# import UL & LL datasets
# UL 
disease_progUL <- read.csv(file.path(DataFolder, "diseaseProgress_UL.csv"), 
                           header = TRUE)
disease_progressUL <- c(disease_progUL[, -1])

# LL 
disease_progLL <- read.csv(file.path(DataFolder, "diseaseProgress_LL.csv"), 
                           header = TRUE)
disease_progressLL <- c(disease_progLL[, -1])



Param_disease_progress <- randomParam_progress(disease_progressLL, 
                                               disease_progressUL, 
                                               disease_progress) 


# cured disease progress 
cured_progUL <- read.csv(file.path(DataFolder, "curedProgress_UL.csv"), 
                           header = TRUE)
cured_progressUL <- c(cured_progUL[, -1])

# LL 
cured_progLL <- read.csv(file.path(DataFolder, "curedProgress_LL.csv"), 
                           header = TRUE)
cured_progressLL <- c(cured_progLL[, -1])



Param_fib <- randomParam_progress(cured_progressLL, 
                                  cured_progressUL, 
                                               fib) 





################################################################################

# population array sampling parameter set 
Param_pop_array <-lapply(randomParams_array, function(x) x*pop_array)




################################### this need to revise ########################
# It's testing cascade
# Problem 1: memory is not big enough when run simulation = 1000
# Problem 2: current HCV model function seems not fit this list 
# Problem 3: need to find how to set names for list of list when the list turn inside out 
# dfList 
# increase memory limits
memory.limit(60000) 
# substr the timepoint for simulation to save the memory space 

# remember to set up endY for simulation in the 01_SetupModel
seq_time <- head(seq(1, 50, 0.1), -1) # it should be revised as seq(1, endY-1,HCV$timestep)
length(seq_time)
df_test <- lapply(dfList, function(x) x[ , , 1:length(seq_time) + 1] )



# extract name and length of dfList 

files <- list.files(path = paste0(DataFolder, 
                                  "/parameter_varied_stages/", sep = ""),
                    pattern = '*.csv')

Param_dfList <- list()
a <- list()

for (i in 1:length(files)) { 
  Param_dfList[[i]] <- lapply(randomParams_array, function(x){ 
    a[[i]] <- x * df_test[[i]] 
    ifelse(a[[i]] >= 1,0.98, a[[i]])  # probability between 0-1 if above 1 as 1 
  })
}




# turn the list inside out to align other paramsets
Param_dfList <- purrr::transpose(Param_dfList)

# extract names of cascade 



for (i in 1: length(Param_dfList)) { 
  
 names(Param_dfList[[i]]) <- c(gsub("^|.csv", "", files)) 
  
}

################################################################################

## constants beta and mortality_b
timV <- c("beta1", "beta2", "beta3", "beta4", "mortality_b")
PPAA <- list()
  for (var in timV) {
    for (i in 1: HCV$numberSamples){ 
      PPAA[[var]][[i]] <- randomParams_array[i]*best_estimates[, var]
    }
  } 




# extract the UL and  LL  of parameters from the file
param_constant <- read.csv(file.path(DataFolder, "parameters_constants.csv"), 
         header = TRUE)%>%select(-"X")%>%as.data.frame()

Param_con <- param_constant%>%
  filter(!parameter%in%timV & !parameter%in% c("HCVP1", "HCVP2", "HCVP3", 
                                               "HCVP4", "entry"))%>%
  gather("estimate", "value", -parameter)


Param_C<- list()

for (var in unique(Param_con$parameter)) {
  Param_C[[var]] <- randomParams(var, Param_con) 
}

# entry
PPAA[["entry"]] <- lapply(Param_C[["MSM_pro"]], 
                          function(x) best_estimates$entry/best_estimates$MSM_pro * x)


Param_Cparam <- lapply(Param_C, function(x) as.data.frame(t(x)))
Param_Constant <- lapply(Param_Cparam, function(x) x[rep(1,HCV$npts), ] )


################################################################################
Param_estimates <- list()
indices <- c()

yearV <- list()
for (i in 1: HCV$numberSamples) { 
  Param_estimates[[i]] <- best_estimates[1:length(seq_time) + 1, ]
  
  for (var in names(PPAA)) {
    Param_estimates[[i]][[var]] <- PPAA[[var]][[i]][1:length(seq_time) + 1]
  }
  
  for (var in c(unique(Param_con$parameter))) { 
      for ( year in 1:(50 - 1)) { 
        indices <- ((year-1) * 1/HCV$timestep + 1): (year * 1/HCV$timestep + 1)
        yearV[[var]][[i]] <- seq(Param_Constant[[var]][year, i], 
                               Param_Constant[[var]][year + 1, i], 
                               length = (1 + 1/HCV$timestep))
        Param_estimates[[i]][indices, var] <- yearV[[var]][[i]] 
    }
    }}




# initial pop 

#### number of MSM in each population ####
initialPops<- read.csv(file.path(DataFolder, "initial_populations.csv"), 
                       header = TRUE)%>%select(-"X")

estPops<- read.csv(file.path(DataFolder, "Estimate_initial_pop.csv"), 
                   header = TRUE)%>%select(-"X")

init_pop <- filter(initialPops, parameter == "init_pop")$value

param_MSMprop <- lapply(Param_C[["MSM_pro"]], function(x) x*init_pop )

pop_prop <- initialPops%>%filter(parameter%in% c("pop_prop1", "pop_prop2", 
                                                 "pop_prop3", "pop_prop4"))%>%
  select(value)%>%unlist()%>%as.vector()

popProp <- lapply(param_MSMprop, function(x) as.numeric(x)*pop_prop )


# prevalence at initial
init_prop_I <- c(constantsDf$HCVP1[1], constantsDf$HCVP2[1], 
                 constantsDf$HCVP3[1], constantsDf$HCVP4[1])
init_prop_S <-c(1-init_prop_I)
Param_estPops <- list()
Param_Pops <- list()
for ( i in 1: HCV$numberSamples){ 
  Param_estPops[[i]] <- estPops%>%mutate(
  pop_group = rep(c(popProp[[i]]),dim(estPops)[1]/HCV$npops),
  SIprop = case_when(Var1 == "HIV-" & SI == "S" ~ init_prop_S[1],
                     Var1 == "HIV-" & SI == "I" ~ init_prop_I[1],
                     Var1 == "HIV-PrEP" & SI == "S" ~ init_prop_S[2],
                     Var1 == "HIV-PrEP" & SI == "I" ~ init_prop_I[2],
                     Var1 == "HIV+" & SI == "S" ~ init_prop_S[3],
                     Var1 == "HIV+" & SI == "I" ~ init_prop_I[3],
                     Var1 == "HIV+d" & SI == "S" ~ init_prop_S[4],
                     Var1 == "HIV+d" & SI == "I" ~ init_prop_I[4]),
  est_pop = value*pop_group*SIprop)
  
  Param_Pops[[i]] <- as.matrix(as.data.frame(matrix(Param_estPops[[i]]$est_pop, 
                                               ncol = HCV$ncomponent,  
                                               nrow = HCV$npops)))
  }



Param_Pops <- lapply(Param_Pops, "colnames<-", c(HCV$component_name))


save(Param_estimates, 
     Param_Pops, 
     Param_disease_progress, 
     Param_pop_array, 
     Param_dfList,
     Param_fib,
     randomParams_array,
     file = file.path(projectFolder,
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
