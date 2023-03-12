# This script is preparing the parameters for scenarios 
# According to different testing, there are 4 scenarios lists 
# in each list, will incorporate the different coverage of testing 
# library
rm(list = ls(all.names = TRUE))

library(readxl)
library(dplyr)
gc()
## setup scenarios 
basePath <- getwd()
ScenarioFolder <- file.path(basePath, "01. DATA/model input/Scenarios") 
Rcode <- file.path(basePath, "03. Code")
project <- "Taiwanese MSM HCV model"
project_name <- "HCVModel"

projectFolder <- file.path(basePath)
source(file.path(Rcode, "/Functions/Scenarios_wrangling.R")) 

projectFile <- file.path(basePath,
                         paste0(project_name, ".rda"))
projectVars <- load(projectFile)
load(file.path(projectFolder, paste0(project_name, ".rda")))
load(file.path(projectFolder, paste0(project_name, "param",".rda")))


# import exacel file for odds ratio and uncertainty range for scenarios 

Param_scen <- list()

shetname <- c("tau_ab", "tau_poct" , "tau_ag", "eta") 

Param_scen <- lapply(shetname, function(x) { 
  dt <-read_excel(paste0(ScenarioFolder,"/Param_scenarios.xlsx"), sheet=x)
  })

names(Param_scen) <- shetname

# inititaing year 
calibrateY <- 2004

scenario_Y <- 2022

uptake_Y <- 2024 

# finding time point of start of 2022
scenario_pt <- length(seq(calibrateY, scenario_Y, by = HCV$timestep))

uptake_pt <- length(seq(calibrateY, uptake_Y, by = HCV$timestep)) 

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
seq_time <- head(seq(1, 50,0.1), -1) # it should be revised as seq(1, endY-1,HCV$timestep)
length(seq_time)

sValues <- list()
sParam <- 1000
sValues <- lapply(c(names(Param_scen)), function(x){ 
  lapply(c(SceName), function(m) { 
    filter(Param_scen[[x]], scenarios == m)
    })
  })

names(sValues) <- names(Param_scen)

ransParam <- lapply(c(names(Param_scen)), function(x) 
  SValParam(sValues, sParam, indic = x))

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

df_test <- lapply(dfList,function(x) x[ , , 1:length(seq_time)+1])

Scen_CascadeParam <- list() 

for (i in c(SceName)){ 
  for( x in c(names(Param_scen))){ 
    
    Scen_CascadeParam[[i]][[x]] <- lapply(seq_along(1:sParam), function(mt){ 
      a <- Scen_test(df_test, Param_scenPara[[mt]], x , i, 
                     scenario_pt = scenario_pt, 
                     scenario_Y = scenario_Y, 
                     uptake_Y = uptake_Y)})
  }
}

# giving names for the list of list  
Scen_mainParam <- list() 

for(i in c(SceName)){
  for(mt in seq_along(1:1000)){
    Scen_mainParam[[i]][[mt]] <- Param_dfList[[mt]]
  }
}

for(i in c(SceName)){
  for(mt in seq_along(1:1000)){
      for( x in c(names(Param_scen))){ 
        Scen_mainParam[[i]][[mt]][[x]][ , , c(scenario_pt:dim(Scen_mainParam[[1]][[1]][[1]])[3])] <-  
          Scen_CascadeParam[[i]][[x]][[mt]][ , , c(scenario_pt:dim(Scen_mainParam[[1]][[1]][[1]])[3])]
      }
    }
  }
rm(Scen_CascadeParam)
gc()
# extract aiming coverage 

endP <- list()
endP <- lapply(seq_along(SceName), 
               function(n) lapply(c(names(Param_scen)),function(x) 
                 Scen_main[[n]][[x]][4 , , uptake_pt]))


names(endP) <-SceName 

for( i in 1: length(endP)){ 
  
  names(endP[[i]]) <- c(names(Param_scen)) 
}


Scen_targetpop <- Scen_main

for( m in names(Scen_targetpop)){
  for (x in names(Param_scen)){ 
    
    Scen_targetpop[[m]][[x]] <- tarpop(Scen_targetpop[[m]], x , endP[[m]])
    }
  }

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

gc()

library(doMC)

registerDoMC(cores = 2) 

a <- list()

for( m in SceName){
  for(u in 1: HCV$numberSamples){
    for (x in names(Param_scen)){ 
      a <-tarpop(Scen_targetpop[[m]], x , 
                                                  endPin[[u]][[m]])
    }
  }
}

Scen_targetpopParam <-Scen_mainParam

ldt <- dim(Scen_targetpopParam[[1]][[1]][[1]])[3]

for( m in SceName){
  for(u in 1: HCV$numberSamples){
    for (x in names(Param_scen)){ 
      Scen_targetpopParam[[m]][[u]][[x]][2, , c(scenario_pt:ldt)] <-
        tarpop(Scen_targetpop[[m]], x , endPin[[u]][[m]])[2, ,c(scenario_pt:ldt)]
    }
  }
}

gc()

#### time frame: L 5y; U: 1y  ####
#### timeframe 5y #### 
Scen_timeL <- list()

# parameters of stages set up for main scenarios   
for (i in 1: length(SceName)){ 
  
  Scen_timeL[[i]] <- dfList
}

names(Scen_timeL) <- SceName 

# replace the parameters in cascade 
Scen_cascadeL <- list()

Scen_cascadeL <- lapply(seq_along(SceName), 
                       function(n) lapply(c(names(Param_scen)), function(x){ 
                         a <- Scen_test(dfList, Param_scen, x , n, 
                                        scenario_pt = scenario_pt,
                                        scenario_Y = scenario_Y, 
                                        uptake_Y = uptake_Y + 3)}))

names(Scen_cascadeL) <-SceName 

for( i in 1: length(Scen_cascadeL)){ 
  
  names(Scen_cascadeL[[i]]) <- c(names(Param_scen)) 
}

# replace the parameters in scenarios set 
for ( i in 1: length(SceName)){ 
  for(j in c(names(Param_scen))){ 
    Scen_timeL[[SceName[i]]][[j]] <- Scen_cascadeL[[i]][[j]]
  }
}

####1y ####
#### timeframe 5y #### 

Scen_timeU <- list()

# parameters of stages set up for main scenarios   
for (i in 1: length(SceName)){ 
  
  Scen_timeU[[i]] <- dfList
}

names(Scen_timeU) <- SceName 

# replace the parameters in cascade 
Scen_cascadeU <- list()

Scen_cascadeU <- lapply(seq_along(SceName), 
                        function(n) lapply(c(names(Param_scen)), function(x){ 
                          a <- Scen_test(dfList, Param_scen, x , n, 
                                         scenario_pt = scenario_pt,
                                         scenario_Y = scenario_Y, 
                                         uptake_Y = uptake_Y - 1)}))

names(Scen_cascadeL) <-SceName 

for( i in 1: length(Scen_cascadeU)){ 
  
  names(Scen_cascadeU[[i]]) <- c(names(Param_scen)) 
}

# replace the parameters in scenarios set 
for ( i in 1: length(SceName)){ 
  for(j in c(names(Param_scen))){ 
    Scen_timeU[[SceName[i]]][[j]] <- Scen_cascadeU[[i]][[j]]
  }
}

#### uncertainty #### 
tempf <- list()
Scen_CascadeParamL <- list() 

for (i in c(SceName)){ 
  for( x in c(names(Param_scen))){ 
    
    tempf[[i]][[x]] <- lapply(seq_along(1:HCV$numberSamples), function(mt){ 
      a <- Scen_test(df_test, Param_scenPara[[mt]], x , i, 
                     scenario_pt = scenario_pt, 
                     scenario_Y = scenario_Y, 
                     uptake_Y = uptake_Y + 3 )})
  }
}


# giving names for the list of list  
Scen_timeLParam <- list() 

Scen_timeLParam <- Scen_mainParam

for( m in SceName){
  for(u in 1: HCV$numberSamples){
    for (x in names(Param_scen)){
      Scen_timeLParam[[m]][[u]][[x]][, ,c(scenario_pt: ldt)] <-  
        tempf[[m]][[x]][[mt]][ , , c(scenario_pt: ldt)]
    }
  }
}


  
tempf <- list()
gc()
Scen_CascadeParamU <- list() 

for (i in c(SceName)){ 
  for( x in c(names(Param_scen))){ 
    
    tempf[[i]][[x]] <- lapply(seq_along(1:HCV$numberSamples), function(mt){ 
      a <- Scen_test(df_test, Param_scenPara[[mt]], x , i, 
                     scenario_pt = scenario_pt, 
                     scenario_Y = scenario_Y, 
                     uptake_Y = uptake_Y - 1  )})
  }
}

# giving names for the list of list  
Scen_timeUParam <- list() 
Scen_timeUParam <-  Scen_mainParam

for(i in c(SceName)){
  for(mt in seq_along(1:sParam)){
    for( x in c(names(Param_scen))){ 
      Scen_timeUParam[[i]][[mt]][[x]][ , ,c(scenario_pt: ldt)] <-  
        tempf[[i]][[x]][[mt]][ , ,c(scenario_pt: ldt)]
    }
  }
}


gc()  

#### diagnosed rate: changing screening coverage to impact diagnosed rate #### 

# diagnosis rate: varying the screening to have various diagnosis rate scenario
# POCab: (0.3, 0.8)
# DBS: (0.3, 0.8)
# Reflex: (0.3, 0.8)
# POCRNA: (0.3, 0.8)

####
Scen_diagL <-ScenCover(dt = Scen_main, endP = 0.3, scenario_pt = scenario_pt,
                  uptake_pt = uptake_pt, SceName = SceName)

Scen_diagLT <- lapply(Scen_diagL,function(x) 
  lapply(x, function(y) y[ , , 1:length(seq_time)+1]))

Scen_diagLPar <- Scen_diagL


# only tau_ab using randomParams array, other parameters same as scen_mainparam
for(i in SceName){ 
  Scen_diagLPar[[i]] <- Scen_Uncert(Scen_diagLT[[i]], "tau_ab", randomParams_array)
}


Scen_diagLParam <- Scen_mainParam

for(i in SceName){
  for(m in 1: HCV$numberSamples){
    Scen_diagLParam[[i]][[m]][["tau_ab"]][ , ,c(scenario_pt: ldt)] <-  
      Scen_diagLParam[[i]][[m]][["tau_ab"]][ , ,c(scenario_pt: ldt)]
    }
  }
 

Scen_diagU <-ScenCover(dt = Scen_main, endP = 0.8, scenario_pt = scenario_pt,
                       uptake_pt = uptake_pt, SceName = SceName)


Scen_diagUT <- lapply(Scen_diagU,function(x) 
  lapply(x, function(y) y[ , , 1:length(seq_time)+1] ))

Scen_diagUPar <- Scen_diagU


# only tau_ab using randomParams array, other parameters same as scen_mainparam
for(i in SceName){ 
  Scen_diagUPar[[i]] <- Scen_Uncert(Scen_diagUT[[i]], "tau_ab", randomParams_array)
}


Scen_diagUParam <- Scen_mainParam

for(i in SceName){
  for(m in 1: HCV$numberSamples){
    Scen_diagUParam[[i]][[m]][["tau_ab"]][ , ,c(scenario_pt: ldt)] <-  
      Scen_diagUParam[[i]][[m]][["tau_ab"]][ , ,c(scenario_pt: ldt)]
  }
}


reinfL <- matrix(0, ncol = HCV$npts + 1, nrow = HCV$npops)
dim(reinfL)

reinfL[,1:scenario_pt-1] <- 1
reinfL[,scenario_pt:dim(reinfL)[2]] <- 0


reinfU <- reinfL
reinfU[,scenario_pt:dim(reinfU)[2]] <- 1.5



####save rda ####
save(HCV, sParam, SceName,Scen_main, Scen_targetpop, Scen_timeL, Scen_timeU,
     Scen_diagL, Scen_diagU, reinfL, reinfU,
     file = file.path(projectFolder, 
                      paste0(project_name, "Scenarios", ".rda")))

save(Scen_mainParam,
     file = file.path(projectFolder, 
                      paste0(project_name, "Scenarios_mainpar", ".rda")))

save(Scen_targetpopParam,
     file = file.path(projectFolder, 
                      paste0(project_name, "Scenarios_tarpoppar", ".rda")))

save(Scen_timeLParam,
     file = file.path(projectFolder, 
                      paste0(project_name, "Scenarios_timeLpar", ".rda")))
save(Scen_timeUParam,
     file = file.path(projectFolder, 
                      paste0(project_name, "Scenarios_timeUpar", ".rda")))

save(Scen_diagLParam,
     file = file.path(projectFolder, 
                      paste0(project_name, "Scenarios_diagLpar", ".rda")))

save(Scen_diagUParam,
     file = file.path(projectFolder, 
                      paste0(project_name, "Scenarios_diagUpar", ".rda")))


