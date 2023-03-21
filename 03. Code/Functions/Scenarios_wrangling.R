# This script process the parameters related to scenarios 
# In our scenarios, there are 4 different testing strategies: 
# POC antibody testing, DBS, Reflex RNA and POC RNA testing
# Those numbers in the excel file are the odds compared to current practice


#### library ####
library("readxl")
library("dplyr")
library("tidyr")
library("purrr")
library("gridExtra")

PtoR <- function(p){ 
  
  rate <- -log(1-p)
  return(rate)
}

RtoP <- function(u){ 
  
  prob <- 1-exp(-u)
  return(prob)
}


value <- function(dt,cascade , scenario){ 
  va <- dt[[cascade]]%>%filter(scenarios == scenario)%>%
    select(odds)%>%as.numeric()
  
  return(va)
}


#### get parameters for scenarios #### 
Scen_test <- function(dfList, dt, cascadeP, scenario, scenario_pt, 
                      uptake_Y, scenario_Y){
  
  # dfList: cascade parameters dataset 
  # dt: scenario parameter set (Risk ratio) 
  # cascadeP: stage of cascade such as tau_ab, tau_ag, etc...
  # scenario: testing scenario such as "POC_antibody", etc... 
  # scenario_pt: the timepoint of scenario initiating 
  # uptake_Y: the year of end of uptake scenario 
  #         For example: a scenario scaling up from 2022 to 2024. 
  #                      2024 is the uptak_Y
  # scenario_Y: the year of scenario uptaking 
  
  # time-points 
  scenario_pt <- length(seq(HCV$cabY, scenario_Y, by = HCV$timestep))
  uptake_pt <- length(seq(HCV$cabY, uptake_Y, by = HCV$timestep))
  npt <- dim(dfList[[1]])[3]
  # extract the dimension from dfList 
  initia_value <- dfList[[cascadeP]][ , ,scenario_pt] 
  
  
  mm <- list()
  
  a <- dfList[[cascadeP]]
  
  Scenario_ab <- dt[[cascadeP]][ ,"scenarios"]%>%unlist()%>%
    as.character()%>%c()
  
  oddValue <- lapply(Scenario_ab, function(x) value(dt, cascadeP, x))
  
  names(oddValue) <- Scenario_ab
  
  # Testing scenarios impact the parameters in our model, 
  # including tau_ab, tau_ag, tau_poct (for POC RNA scenario) and eta
  
  ##### tau_ab ##### 
  # tau_ab is the RR for the parameter. 
  # Estimated as oddValue[[scenario]](RR)*dfList[[cascade]]
  if(isTRUE(cascadeP=="tau_ab")){ 
    
    if(isTRUE(oddValue[[scenario]]==0)){
      
      a0 <- dfList[["tau_poct"]][ , , 1]
      mm <- a
      mm[, , scenario_pt:npt] <- a0
      return(mm)
    }
    else if(isTRUE(oddValue[[scenario]]==1)){ 
      mm <- a
      }
    else{
      mm <-a 
      tempVal <- mm[, , scenario_pt]
      t <- tempVal
      
      for(m in 1: (uptake_pt - scenario_pt + 1)){ 
         for ( j in 2:dim(mm)[2]){
            for(i in 1:(dim(mm)[1])){
          
          t[i,j] <- oddValue[[scenario]]*tempVal[i,j]
          
          mm[i, j,(m + scenario_pt -1)] <- seq(as.numeric(a[i,j, scenario_pt]), 
                                                 as.numeric(t[i,j]), 
                              length = (uptake_pt - scenario_pt + 1))[m]
            }
        }
        
      }
      mm[ , , uptake_pt: npt] <- mm[ , , uptake_pt]
      mm <- replace(mm, mm>=0.98,0.98) # limit: 98% 
      return(mm)
      }
  }
  #####tau_ag #####
  # tau_ag is the target annual probability. Estimated as oddValue[[scenario]]
  else if(isTRUE(cascadeP== "tau_ag")){
    
    if(isTRUE(oddValue[[scenario]]==0)){
      
      a0 <- dfList[["tau_poct"]][ , , 1]
      mm <- a
      mm[, , scenario_pt:npt] <- a0
      return(mm)
      
    }
    else if(isTRUE(oddValue[[scenario]]==1)){ 
      mm <- a
    }
    else{
      mm <-a 
      tempVal <- mm[, , scenario_pt]
      t <- tempVal
      
      for(m in 1: (uptake_pt - scenario_pt + 1)){ 
        for ( j in 2:(dim(mm)[2])){
          for(i in 1:(dim(mm)[1])){
            
            t[i,j] <- oddValue[[scenario]]
            mm[i, j,(m + scenario_pt -1)] <- seq(as.numeric(tempVal[i,j]), 
                                                 as.numeric(t[i, j]), 
                                                 length = (uptake_pt - 
                                                          scenario_pt + 1))[m] 
          }
          
        }
        
      }
      
      mm[ , , uptake_pt: npt] <- mm[ , , uptake_pt]
      mm <- replace(mm, mm>=0.98,0.98)
      return(mm)
    }
  }
  ##### eta #####
  # eta is the RR for the parameter. 
  # Estimated as oddValue[[scenario]](RR)*dfList[[cascade]]
  else if(isTRUE(cascadeP== "eta")){
    if(isTRUE(oddValue[[scenario]]==0)){

      a0 <- dfList[["tau_poct"]][ , , 1]
      mm <- a
      mm[, , scenario_pt:npt] <- a0
      return(mm)
      
    }
    else if(isTRUE(oddValue[[scenario]]==1)){ 
      mm <- a
    }
    else{
        mm <-a 
        tempVal <- mm[, , scenario_pt]
        t <- tempVal
        
        for(m in 1: (uptake_pt - scenario_pt + 1)){ 
          for ( j in 2:(dim(mm)[2])){
            for(i in 1:(dim(mm)[1])){
              
              t[i,j] <- oddValue[[scenario]]*tempVal[i,j]
              mm[i, j,(m + scenario_pt -1)] <- seq(as.numeric(tempVal[i,j]), 
                                                   as.numeric(t[i, j]), 
                                                   length = (uptake_pt - 
                                                            scenario_pt + 1))[m] 
              }
            
            }
            
        }
        
        mm[ , , uptake_pt: npt] <- mm[ , , uptake_pt]
        mm <- replace(mm, mm>=0.98,0.98)
        return(mm)
        }
  }
  ##### tau_poct #####
  # tau_poct estimated as tau_ab*tau_ag*oddValue[[scenario]]
  else if(isTRUE(cascadeP=="tau_poct")){ 
    
    if(isTRUE(oddValue[[scenario]]!=0)){
      
      mm <- a
      tempVal <- mm[, , scenario_pt]
      t <-  tempVal
      mm[, , scenario_pt] <- dfList[["tau_ab"]][ , , scenario_pt]*
        dfList[["tau_ag"]][ , , scenario_pt]
      
      diagVal <- matrix()
      
      diagVal <- dfList[["tau_ab"]][ , , scenario_pt]
      
      diapop <- diagVal[ 1, 3]
      
      daipopD <- diagVal[4,3]
      
      for(m in 1: (uptake_pt - scenario_pt + 1)){ 
        for ( j in 2:dim(mm)[2]){
          for(i in 1:(dim(mm)[1])){
            
            t <- oddValue[[scenario]]*mm[, , scenario_pt]
            
            mm[i, j,(m + scenario_pt -1)] <- seq(as.numeric(diagVal[i, j]), 
                                                 as.numeric(t[i,j]), 
                                                 length = (uptake_pt - 
                                                          scenario_pt + 1))[m]
          }
          
          
        }
        
      }
      
      mm[ , , uptake_pt: npt] <- mm[ , , uptake_pt]
      mm <- replace(mm, mm>=0.98,0.98)
      return(mm) 
    
    }else{ mm <- a}
  }
    
  else{ 
    
    print(paste0("This function does not work for indicator-", cascadeP, 
                 sep = ""))
    }
}

#### uncertainty range for scenarios ####
Scen_Uncert <- function(df_test, CasName, randomParams_array){
  a <- list()
  for (i in CasName){
    for(m in 1: HCV$numberSamples){ 
      Param_dfList[[m]][[i]] <- randomParams_array[m]*df_test[[i]] 
        ifelse(Param_dfList[[m]][[i]]>=1,0.98, Param_dfList[[m]][[i]])  # probability between 0-1 if above 1 as 1
      
      }
    } 
  return(Param_dfList)
  }





#### function for get the random parameter for parameter of scenario between lower
#### bound and upper bound. 
SValParam <- function(sValues, sParam,indic){
  
  # sValues: get the odds from Param_scen dataset 
  # sParam: number of sampling 
  # indic: indicator of cascade: tau_ab, tau_ag... etc
  sValuesParam <- list()
  for (i in 1: length(sValues[[1]])){ 
    sValuesParam[[i]] <- 
      c(runif(sParam, sValues[[indic]][[i]]$lower, sValues[[indic]][[i]]$upper))
    
    
  } 
  names(sValuesParam) <- c(Param_scen$tau_ab$scenarios)
  return(sValuesParam)
}


#### sensitivity for scale up in PrEP users and HIV+d #### 
tarpop <- function(nlist, indic, endval){ 
  # nlist: the list of cascade parameters 
  # indic: c(indicator names), such as tau_ab, tau_ag...
  # endVal: the list has the parameter aim to scale up to 
  # will be like :endP$POCab[[1]]
     
      for (i in 2: dim(nlist[[indic]])[2]){ 
      nlist[[indic]][2,i,  scenario_pt:uptake_pt] <- 
        seq(as.numeric(nlist[[indic]][2, i, scenario_pt]), 
            as.numeric(unlist(endval[[indic]])[2]), 
            length = (uptake_pt - scenario_pt + 1))
      
      nlist[[indic]][2,i,  uptake_pt: dim(nlist[[indic]])[3]] <- 
        as.numeric(unlist(endval[[indic]])[i])
      
      nlist[[indic]][4,i,  scenario_pt:uptake_pt] <- 
        seq(as.numeric(nlist[[indic]][2, i, scenario_pt]), 
            as.numeric(unlist(endval[[indic]])[2]), 
            length = (uptake_pt - scenario_pt + 1))
      
      nlist[[indic]][4,i,  uptake_pt: dim(nlist[[indic]])[3]] <- 
        as.numeric(unlist(endval[[indic]])[i])
    
      
      }
    
  nlist[[indic]] <- replace(nlist[[indic]], nlist[[indic]]>=0.98,0.98)
    return(nlist[[indic]])
    
    
  }
# test 
tarpoptt <- function(nlist, indic, endval){ 
  # nlist: the list of cascade parameters 
  # indic: c(indicator names), such as tau_ab, tau_ag...
  # endVal: the list has the parameter aim to scale up to 
  # will be like :endP$POCab[[1]]
  
  for (i in 2: dim(nlist[[indic]])[2]){ 
    nlist[[indic]][,i,  scenario_pt:uptake_pt] <- 
      seq(as.numeric(nlist[[indic]][, i, scenario_pt]), 
          as.numeric(unlist(endval[[indic]])[2]), 
          length = (uptake_pt - scenario_pt + 1))
    
    nlist[[indic]][,i,  uptake_pt: dim(nlist[[indic]])[3]] <- 
      as.numeric(unlist(endval[[indic]])[i])
    
    
  }
  
  nlist[[indic]] <- replace(nlist[[indic]], nlist[[indic]]>=0.98,0.98)
  return(nlist[[indic]])
  
  
}



tarpopParm <- function(nlist, indic, endval){ 
  # nlist: the list of cascade parameters 
  # indic: c(indicator names), such as tau_ab, tau_ag...
  # endVal: the list has the parameter aim to scale up to 
  # will be like :endP$POCab[[1]]
  
  
    a <- nlist[[indic]]
    m <- as.numeric(a[2, 2, scenario_pt])
    for (i in 2: dim(a)[2]){ 
      a[2,i,  scenario_pt:uptake_pt] <- 
      seq(m, 
          as.numeric(endval), 
          length = (uptake_pt - scenario_pt + 1))
    
      a[2,i,  uptake_pt: dim(a)[3]] <- 
      as.numeric(endval)
      }
    a<- replace(a, a>=0.98,0.98) 
    a <- as_tibble()
  return(a)
  
    }
  
  
####function for modify aiming coverage for screening in sensitivity analysis####

ScenCover <-function(dt, endP, scenario_pt, uptake_pt, SceName){ 
  
  dL <- dt[[1]]$tau_ab[ , , scenario_pt -1]
  for(i in 2: dim(dL)[2]){ 
    for(j in 1: dim(dL)[1]){ 
      dL[j,i] <- endP
    }
  }
  
  Scen_diagL <- dt
 
    for(m in head(SceName, -1)){ 
      for(i in 2: dim(dL)[2]){ 
        for(j in 1: (dim(dL)[1] - 1)){ 
          Scen_diagL[[m]]$tau_ab[j,i,scenario_pt:uptake_pt] <- 
            seq(Scen_diagL[[m]]$tau_ab[j,i,scenario_pt], dL[j,i],
                length = (uptake_pt - scenario_pt + 1)) 
          Scen_diagL[[m]]$tau_ab[j,i, c(uptake_pt:dim(Scen_diagL[[m]]$tau_ab)[3])] <-
            Scen_diagL[[m]]$tau_ab[j,i, uptake_pt]
        }
      }
    }
    
    for(i in 2: dim(dL)[2]){ 
      for(j in 1: (dim(dL)[1] - 1)){ 
        Scen_diagL[["POC_RNA"]]$tau_poct[j,i,scenario_pt:uptake_pt] <- 
          seq(Scen_diagL[["POC_RNA"]]$tau_poct[j,i,scenario_pt], dL[j,i],
              length = (uptake_pt - scenario_pt + 1)) 
        
        Scen_diagL[["POC_RNA"]]$tau_poct[j,i, 
                                         c(uptake_pt:
                                             dim(Scen_diagL[["POC_RNA"]]$tau_poct)[3])] <-
          Scen_diagL[["POC_RNA"]]$tau_poct[j,i, uptake_pt]
      }
    }
  return(Scen_diagL)
  }


  

 
 