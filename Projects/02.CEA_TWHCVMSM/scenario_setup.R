# This script is setting up scenarios for main analysis in project "CEA_TWHCVMSM"
# focusing on the parameters related to testing strategies and costs. 
# The code used in this script needs to source from "TWHCV-model"
rm(list = ls()) 
# library
library(readxl)
library(dplyr)
library(here)
library(doMC)
library(abind)
registerDoMC(cores = 2) 

# set up dircetory 
basePath <- here()

#file path of "TWHCV-model" project
codepath <- file.path(here() %>% dirname(), 'TWHCV-model/03. Code/Functions')

scenariopath <- file.path(here() %>% dirname(), 
                          'TWHCV-model/01. DATA/model input/Scenarios')

DataFolder <- file.path(here(), "01. DATA/model input")

# Rda file path 
# load the .rda file of base estimate 
rdapath <- file.path(here()%>%dirname(), "Taiwan-MSM-HCV-model")

load(file.path(rdapath , paste0("HCVModel", "param",".rda")))

load(file.path(rdapath , paste0("HCVModel",".rda")))

# output file path 
# dir.create("02. Output") # create subdircetory 
# dir.create("02. Output/RDA") 
outputdt <- here("02. Output/RDA")

# source 
source(file.path(codepath, "HCV_model.R"))

source(file.path(codepath, "plotFunctions.R"))

source(file.path(codepath, "plotOptions.R"))

source(file.path(codepath, "Scenarios_wrangling.R"))


# import exacel file for odds ratio and uncertainty range for scenarios 
Param_scen <- list()

shetname <- c("tau_ab", "tau_poct" , "tau_ag", "eta") 

Param_scen <- lapply(shetname, function(x) { 
  
  dt <-read_excel(paste0(scenariopath,"/Param_scenarios.xlsx"), sheet=x)
})

names(Param_scen) <- shetname

# inititaing year 
calibrateY <- 2004

scenario_Y <- 2022

uptake_Y <- 2024 

PrEPcov_Y <- 2030 

# finding time point of start of 2022
scenario_pt <- length(seq(calibrateY, scenario_Y, by = HCV$timestep))

uptake_pt <- length(seq(calibrateY, uptake_Y, by = HCV$timestep))

PrEPcov_pt <- length(seq(calibrateY, PrEPcov_Y, by = HCV$timestep))

# import exacel file for odds ratio and uncertainty range for scenarios 
Param_scen <- list()

shetname <- c("tau_ab", "tau_poct" , "tau_ag", "eta") 

Param_scen <- lapply(shetname, function(x) { 
  
  dt <-read_excel(paste0(scenariopath,"/Param_scenarios.xlsx"), sheet=x)
})

names(Param_scen) <- shetname

#### names of scenarios #### 
SceName <- Param_scen$tau_ab$scenarios

Scen_main <- list()

# parameters of stages set up for main scenarios   
for (i in 1: length(SceName)){ 
  
  Scen_main[[i]] <- dfList
}

names(Scen_main) <- SceName 


# replace the parameters in cascade 
Scen_cascade <- list()

Scen_cascade <- lapply(seq_along(SceName), 
                       function(n) lapply(c(names(Param_scen)), function(x){ 
                         a <- Scen_test(dfList, Param_scen, x , n, 
                                        scenario_pt = scenario_pt,
                                        scenario_Y = scenario_Y, 
                                        uptake_Y = uptake_Y)}))

names(Scen_cascade) <-SceName 

for( i in 1: length(Scen_cascade)){ 
  
  names(Scen_cascade[[i]]) <- c(names(Param_scen)) 
  }

# replace the parameters in scenarios set 
for ( i in 1: length(SceName)){ 
  
  for(j in c(names(Param_scen))){ 
    
    Scen_cascade[[i]][[j]] <- Scen_cascade[[i]][[j]]
    
    Scen_main[[SceName[i]]][[j]] <- Scen_cascade[[i]][[j]]
  }
}

#### uncertainty range ####
set.seed(123456) 

simY <- 100

seq_time <- head(seq(1, simY, 0.1), -1) 

sValues <- list()

sParam <- HCV$sParam

sValues <- lapply(c(names(Param_scen)), function(x){ 
  
  lapply(c(SceName), function(m) { 
    
    filter(Param_scen[[x]], scenarios == m)
  })
})

names(sValues) <- names(Param_scen)

ransParam <- lapply(c(names(Param_scen)), function(x){  
  
  SValParam(sValues, sParam, indic = x)
  })

names(ransParam) <- c(names(Param_scen))

Param_scenPara <- list()

t <- list()

for (i in 1: sParam){ 
  
  for (j in c(names(Param_scen))){ 
    
    t[[j]][[i]] <- c(ransParam[[j]]$POC_antibody[i], 
                     ransParam[[j]]$DBS[i], 
                     ransParam[[j]]$Reflex_RNA[i], 
                     ransParam[[j]]$POC_RNA[i])
  }
} 

Param_scenPara <- rep(list(Param_scen),sParam)

for (i in 1: sParam){ 
  
  for (j in c(names(Param_scen))){ 
    
    Param_scenPara[[i]][[j]]$odds <- t[[j]][[i]]
    
    Param_scenPara[[i]][[j]]$odds <- ifelse(
      
      is.nan(Param_scenPara[[i]][[j]]$odds), 0, Param_scenPara[[i]][[j]]$odds)
  }
}

# extend timepoints for parameters for running simulation for 100 years
ddn <- dim(Param_dfList[[1]]$cured)

Param_dfListExtend  <- list()

Param_dfListExtend  <- lapply(Param_dfList,function(x){
  lapply(x, function(y){ 
    tta <- array(0, c(ddn[1], ddn[2], length(seq_time)), dimnames = dimnames(y))
    tta[ , ,c(1: ddn[3])] <- y[, , c(1:ddn[3])]
    tta[, , c(ddn[3]: length(seq_time))] <- y[ , ,ddn[3]]
  
    return(tta)
    })
})



save(Param_dfListExtend, file = file.path(outputdt, "Param_dfExt.rda"))


## parameters for scenarios 
Scen_CascadeParam <- list() 

for (i in c(SceName)){ 
  
  for( x in c(names(Param_scen))){ 
    
    Scen_CascadeParam[[i]][[x]] <- lapply(seq_along(1:sParam), function(mt){ 
      
      a <- Scen_test(Param_dfListExtend[[mt]], Param_scenPara[[mt]], x , i, 
                     scenario_pt = scenario_pt, 
                     scenario_Y = scenario_Y, 
                     uptake_Y = uptake_Y)})
  }
}

# giving names for the list of list  
Scen_mainParam <- list() 

turnlst <- list()

# turn list of list inside out
 # Scen_mainParam[["testing"]][[1:1000]][["cascade"]]

for(i in c(SceName)){
  
  Scen_mainParam[[i]] <- Param_dfListExtend
  
  turnlst[[i]] <- list_transpose(Scen_CascadeParam[[i]])
  }

for(i in c(SceName)){
  
  for(n in 1: HCV$sParam){
    
    for(x in c(names(Param_scen))){
      
      Scen_mainParam[[i]][[n]][[x]] <- turnlst[[i]][[n]][[x]] 
      }
    }
  }

# extract the target coverage  
endP <- list()

endP <- lapply(seq_along(SceName), 
               
               function(n) lapply(c(names(Param_scen)),function(x) 
                 
                 Scen_main[[n]][[x]][4 , , uptake_pt]))

names(endP) <-SceName 

for( i in 1: length(endP)){ 
  
  names(endP[[i]]) <- c(names(Param_scen)) 
  }

# assign the coverage to the targeted pops 
CEA_main <- Scen_main

# subpops 1 and 3 remain as same as status quo 
# subpops 4 are the main scenario in the main analysis of the TWHCV model project
# here we uptake the coverage along hcv cascade in subpops 2
for( m in names(CEA_main)){
  
  for (x in names(Param_scen)){ 
    
    CEA_main[[m]][[x]] <- tarpop(CEA_main[[m]], x , endP[[m]])
    
    CEA_main[[m]][[x]][1, ,] <- dfList[[x]][1, ,]
    
    CEA_main[[m]][[x]][3, ,] <- dfList[[x]][3, ,]
    
    CEA_main[[m]][[x]][4, ,] <- Scen_main[[m]][[x]][4, , ]
    }
  }

#### uncertainty range ####
# uncertainty 
# extract initial % for HIV diagnosed subgroup for uncertainty 
endPinit <- list()

endPinit <- lapply(seq_along(SceName), 
                   
                   function(n) 
                     
                     lapply(seq_along(1:HCV$numberSamples), 
                            
                            function(mt) lapply(c(names(Param_scen)),function(x)  
                              
                              Scen_mainParam[[n]][[mt]][[x]][4 , , uptake_pt])))


names(endPinit) <-SceName 

for( i in 1: length(endPinit)){
  
  for(x in 1:HCV$numberSamples){
    
    names(endPinit[[i]][[x]]) <- c(names(Param_scen)) 
  }
}

endPin <- purrr::transpose(endPinit)

CEA_mainParam <- Scen_mainParam

for( m in SceName){
  
  for(u in 1: HCV$numberSamples){
    
    for (x in names(Param_scen)){ 
      
      CEA_mainParam[[m]][[u]][[x]] <- tarpop(CEA_mainParam[[m]][[u]], x , 
                 endPin[[u]][[m]])
      
      CEA_mainParam[[m]][[u]][[x]][1, ,] <- Param_dfListExtend[[u]][[x]][1, ,]
      
      CEA_mainParam[[m]][[u]][[x]][3, ,] <- Param_dfListExtend[[u]][[x]][3, ,]
      
      CEA_mainParam[[m]][[u]][[x]][4, ,] <- Scen_mainParam[[m]][[u]][[x]][4, , ]
      }
    }
  }

####save rda ####
save(CEA_main, CEA_mainParam, file = file.path(outputdt, "mainansis.rda"))

#### cost #### 
files <- list.files(path = paste0(DataFolder, 
                                  "/cost/", sep =  ""), pattern = '*.xlsx')


costdfList <- lapply(files, function(f) {
  
  df <- read_excel(file.path(paste0(DataFolder, "/cost/", f, sep = "")))
  df <- df[, -1]
  
  df <- df%>%as_tibble()
  
  df <- as.matrix(df, nrow = npops, ncol = length(.) + 1)
  })

names(costdfList) <- c(gsub("^|.xlsx", "", files)) # ^: from beginning, \ end before .csv


cost <- list()

cost[["state"]] <- costdfList$costPops

cost[["QALY"]] <- costdfList$QALYPops

cost[["flow"]] <- costdfList$costFlow 

cost[["POCab"]] <- costdfList$costFlow_POCab

cost[["Reflex"]] <- costdfList$costFlow_reflex

cost[["POCRNA"]] <- costdfList$costFlow_POCRNA



# move to csv.file with detailed items for HIV and PrEP
TestSce <- "tarpop"

if(isTRUE(TestSce=="tarpop")){ 
  
  cost[["POCab"]][1, ] <- cost[["flow"]][1, ]
  cost[["POCab"]][3, ] <- cost[["flow"]][3, ]
  
  cost[["Reflex"]][1, ] <- cost[["flow"]][1, ]
  cost[["Reflex"]][3, ] <- cost[["flow"]][3, ]
  
  cost[["POCab"]][1, ] <- cost[["flow"]][1, ]
  cost[["POCab"]][3, ] <- cost[["flow"]][3, ]
  
  cost[["POCRNA"]][1, ] <- cost[["flow"]][1, ]
  cost[["POCRNA"]][3, ] <- cost[["flow"]][3, ]
}




#### bugs #### 
cost <- lapply(cost, function(x){ 
  
  a <- replicate(HCV$npts, x)%>%as.numeric()
  })

save(cost, file = file.path(outputdt, "cost.rda"))


cost$flow

