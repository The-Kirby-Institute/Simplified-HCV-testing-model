#### functions for tidying up data and plots ####
# adopted from HBV model of Dr. Richard Gray. 
# Extract annual results -------------------------------------------------
library(reshape2)
library(cowplot)
library(ggplot2)
library(here)

#### tidy results ####
StryearIndex <- function(n, timestep) {
  # This function extracts the mid year index. This used for 
  # calculating per population number annual outputs by dividing
  # annual value by mid year population size
  # 
  # Args:
  #   n: length of vector to extract midyear points
  #   timestep: modelling timestep used
  # Returns:
  #   indices associated with the mid year value
  #   
  #-----------------------------------------------------------------------
  
  seq(1, n, 1 / timestep)
} 


MidyearIndex <- function(n, timestep) {
  # This function extracts the mid year index. This used for 
  # calculating per population number annual outputs by dividing
  # annual value by mid year population size
  # 
  # Args:
  #   n: length of vector to extract midyear points
  #   timestep: modelling timestep used
  # Returns:
  #   indices associated with the mid year value
  #   
  #-----------------------------------------------------------------------
  
  seq(1, n, 1 / timestep) + round(1 / timestep / 2)
} 

# end year index
EndyearIndex <- function(n, timestep) {
  # This function extracts the end of year index. This used for 
  # calculating per population number annual outputs by dividing
  # annual value by end year population size
  # 
  # Args:
  #   n: length of vector to extract end year points
  #   timestep: modelling timestep used
  # Returns:
  #   indices associated with the ned year value
  #   
  #-----------------------------------------------------------------------
  
  seq(1, n, 1 / timestep) + round(1 /timestep)-1
} 
# Organize results -------------------------------------------------------

FactorPop <- function(df, factorLabels) {
  df$population <- factor(df$population)
  levels(df$population) <- factorLabels
  
  return(df)
}

FactorState <- function(df, factorLabels) {
  df$state <- factor(df$state)
  levels(df$state) <- factorLabels
  
  return(df)
}

FactorCascade <- function(df, factorLabels) {
  df$cascade <- factor(df$cascade)
  levels(df$cascade) <- factorLabels
  
  return(df)
}

FactorDiseaseProg <- function(df, factorLabels) {
  df$disease_prog  <- factor(df$disease_prog )
  levels(df$disease_prog ) <- factorLabels
  
  return(df)
}

MidYearDf <- function(df, timestep, npops) {
  yearFrame <-  df[seq(1, nrow(df), (1 / timestep)), ]
  yearFrame <- df[MidyearIndex(nrow(df), timestep), ]
  yearFrame$year <- rep(head(HCV$years,-1), npops)
  return(yearFrame)
}  

EndYearDf <- function(df, timestep, npops) {
  yearFrame <-  df[seq(1, nrow(df), (1 / timestep)), ]
  yearFrame <- df[EndyearIndex(nrow(df), timestep), ]
  yearFrame$year <- rep(head(HCV$years, -1), npops)
  return(yearFrame)
} 

StrYearDf <-  function(df, timestep, npops) {
  yearFrame <-  df[seq(1, nrow(df), (1 / timestep)), ]
  yearFrame <- df[StryearIndex(nrow(df), timestep), ]
  yearFrame$year <- rep(head(HCV$years, -1), npops)
  return(yearFrame)
} 


# Extract the subresults we want
popResults_MidYear <- function(pg, Best, Population = NULL, 
                               Disease_prog = NULL, Cascade = NULL, 
                               MidYear = NULL, 
                               param = NULL, endYear ,YearCut = NULL, 
                               scenario = NULL, allp = NULL) {
  # This function organizes and merges the population size results 
  # example 
  #           popResults(HCV, best, Population = c("HIV-"), Disease_prog = NULL, 
  #           Cascade = c("undiag", "diag_ab"))
  # 
  # NULL means not seperate by variable 
  # time = timestep
  
  
  
  
  # First organize the best estimates  
  df_list <- lapply(Best, as.data.frame.table)
  if(!is.null(allp)){ 
    allpop <- df_list[[allp]]
  }
  else{ 
    allpop <- df_list$allPops
    
    }
  allpop <- allpop%>%
    mutate(time = rep(seq(1.0,(endYear - pg$timestep), pg$timestep), 
                      each=pg$ncomponent*pg$npops),
           cascade = sub("^[^_]*_", "", Var2), 
           disease_prog = sub("\\_.*", "", Var2))%>%
    dplyr::select(-Var3)
  
  names(allpop) <- c("population", "state", "Frequency","timestep", "cascade",
                     "disease_prog")
  
  ## MidyearIndex
  timelong <- seq(pg$startYear, endYear, pg$timestep) 
  if (is.null(YearCut)) {
    MidY <-c(timelong[MidyearIndex(length(timelong),0.1)]) 
    allpop <- allpop%>%filter(timestep%in% MidY)%>%
      mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                    each = pg$npops*pg$ncomponent))
  
  } else if (YearCut =="Start") { 
    
    MidY <-c(timelong[StryearIndex(length(timelong),0.1)]) 
    allpop <- allpop%>%filter(timestep%in% MidY)%>%
      mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                        each = pg$npops*pg$ncomponent))
    
  } else {  
    MidY <-c(timelong[EndyearIndex(length(timelong),0.1)]) 
    allpop <- allpop%>%filter(timestep%in% c(MidY, (endYear - pg$timestep)))%>%
      mutate(year = rep(seq(pg$startYear, (endYear -1),1), 
                        each = pg$npops*pg$ncomponent))
      
    
    }
  

  
#  # now prepare to append parameter sets 
  popSizes <- allpop%>%mutate(best = Frequency)%>%
    select(-c(Frequency))
  dfparam_list <- list()
  param_pop <- list()
  if(is.null(scenario)){ 
    samp <- HCV$numberSamples
  }else{ 
    samp <- HCV$sParam
    }
  
  if (!is.null(param)){
    for (set in 1: samp){
      dfparam_list[[set]] <- as.data.frame.table(param[[set]]$allPops)
      
      param_pop[[set]] <- dfparam_list[[set]]%>%
        mutate(time = rep(seq(1.0,(endYear - pg$timestep), pg$timestep),
                          each=pg$ncomponent*pg$npops),
               cascade = sub("^[^_]*_", "", Var2), 
               disease_prog = sub("\\_.*", "", Var2))%>%
        dplyr::select(-Var3)
        
      
      names(param_pop[[set]]) <- c("population", "state", "Frequency",
                                   "timestep", "cascade","disease_prog")
      timelong <- seq(pg$startYear, endYear, pg$timestep) 
      
      if (is.null(YearCut)) {
        MidY <-c(timelong[MidyearIndex(length(timelong),0.1)]) 
        param_pop[[set]] <- param_pop[[set]]%>%filter(timestep%in% MidY)%>%
          mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                            each = pg$npops*pg$ncomponent))
        
      } else if (YearCut =="Start") { 
        
        MidY <-c(timelong[StryearIndex(length(timelong),0.1)]) 
        param_pop[[set]] <-param_pop[[set]]%>%filter(timestep%in% MidY)%>%
          mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                            each = pg$npops*pg$ncomponent))
        
      } else {  
        MidY <-c(timelong[EndyearIndex(length(timelong),0.1)]) 
        param_pop[[set]] <- param_pop[[set]]%>%
          filter(timestep%in%c(MidY, (endYear - pg$timestep)))%>%
          mutate(year = rep(seq(pg$startYear, (endYear -1),1), 
                            each = pg$npops*pg$ncomponent))
        
        
      } 
      
      
      
      popSizes[, paste0("set", toString(set))] <-  param_pop[[set]]$Frequency
    }}
  
  
  # Extract the subresults we want
  if (is.null(Disease_prog) && is.null(Cascade) && is.null(Population)){
    popSizes <- popSizes%>%
      mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                        each = pg$npops*pg$ncomponent))%>%
      gather("simulation", "popsize", best:ncol(popSizes))%>%
      group_by(year, simulation) %>%
      summarise(best = sum(popsize))%>%
      spread(simulation, best)%>%arrange(year) 
  
    }else if (is.null(Disease_prog) && is.null(Cascade) && !is.null(Population)){ 
    popSizes <- popSizes%>%
      mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                        each = pg$npops*pg$ncomponent))%>%
      gather("simulation", "popsize", 
                                  best:ncol(popSizes))%>%
      group_by(year, simulation, population)%>%
      filter(population%in%Population)%>%
      summarise(best = sum(popsize))%>%
      spread(simulation, best)%>%arrange(year, population)  
  }else if (is.null(Disease_prog) && !is.null(Cascade) && is.null(Population)) {
    popSizes <- popSizes%>%
      mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                        each = pg$npops*pg$ncomponent))%>%
      gather("simulation", "popsize", 
                                  best:ncol(popSizes))%>%
      group_by(year, simulation, cascade)%>%
      filter(cascade%in%Cascade)%>%
      summarise(best = sum(popsize))%>%
      spread(simulation, best)%>%arrange(year, cascade) 
  
    }else if (!is.null(Disease_prog) && is.null(Cascade) && is.null(Population)) {
    popSizes <- popSizes%>%
      mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                        each = pg$npops*pg$ncomponent))%>%
      gather("simulation", "popsize",best:ncol(popSizes))%>%
      group_by(year, simulation, disease_prog)%>%
      filter(disease_prog%in%Disease_prog)%>%
      summarise(best = sum(popsize))%>%
      spread(simulation, best)%>%arrange(year, disease_prog) 
  }else if (!is.null(Disease_prog) && !is.null(Cascade) && is.null(Population)) {
    popSizes <- popSizes%>%
      mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                        each = pg$npops*pg$ncomponent))%>%
      gather("simulation", "popsize", 
                                  best:ncol(popSizes))%>%
      group_by(year, simulation, state)%>%
      filter(disease_prog%in%Disease_prog & cascade %in% Cascade)%>%
      summarise(best = sum(popsize))%>%
      spread(simulation, best)%>%arrange(year, state)  
  }else if (!is.null(Disease_prog) && is.null(Cascade) && !is.null(Population)) {
    popSizes <- popSizes%>%mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                                             each = pg$npops*pg$ncomponent))%>%
      gather("simulation", "popsize", 
                              best:ncol(popSizes))%>%
      group_by(year, simulation, disease_prog, population)%>%
      filter(disease_prog%in%Disease_prog & population %in% Population)%>%
      summarise(best=sum(popsize))%>%spread(simulation, best)%>%
      arrange(year, population, disease_prog)  
  }else if (is.null(Disease_prog) && !is.null(Cascade) && !is.null(Population)) {
    popSizes <- popSizes%>%mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                                             each = pg$npops*pg$ncomponent))%>%
      gather("simulation", "popsize", 
                                  best:ncol(popSizes))%>%
      group_by(year, simulation, cascade, population)%>%
      filter(cascade %in% Cascade && population %in% Population)%>%
      summarise(best = sum(popsize))%>%
      spread(simulation, best)%>%arrange(year, population, cascade) 
  }else { popSizes <- popSizes%>%arrange(year, population, state)
  
  }
  popSizes[popSizes < 0] <- NaN
 
  return(popSizes)
  }


# new dataframe to append aggegrate columns 
popResults_range <- function(pg, data, Population = NULL, 
                             Disease_prog = NULL, Cascade = NULL, end_Y){ 
  
  # new dataframe to append aggegrate columns 
  popSizeCompo <- data
  
  if (is.null(Disease_prog) && is.null(Cascade) && is.null(Population)){
    popSizeCompo <- popSizeCompo%>%arrange(year) 
  } 
  else if (!is.null(Disease_prog) && !is.null(Cascade) && !is.null(Population)){
    popSizeCompo <- popSizeCompo%>%arrange(year, population, state)
  }
  else if (is.null(Disease_prog) && !is.null(Cascade) && !is.null(Population)) {
    popSizeCompo <- popSizeCompo%>%arrange(year, population, cascade) 
  } 
  else if (!is.null(Disease_prog) && is.null(Cascade) && !is.null(Population)) {
    popSizeCompo <- popSizeCompo%>%arrange(year, population, disease_prog)  
  }
  else if (!is.null(Disease_prog) && !is.null(Cascade) && is.null(Population)) {
    popSizeCompo <- popSizeCompo%>%arrange(year, state)  
  }
  else if (is.null(Disease_prog) && is.null(Cascade) && !is.null(Population)){ 
    popSizeCompo <- popSizeCompo%>%arrange(year, population)  
  } 
  else if (is.null(Disease_prog) && !is.null(Cascade) && is.null(Population)) {
    popSizeCompo <- popSizeCompo%>%arrange(year, cascade) 
  } 
  else if (!is.null(Disease_prog) && is.null(Cascade) && is.null(Population)) {
    popSizeCompo <- popSizeCompo%>%arrange(year, disease_prog) 
  } 
  
  
  
  
  # Extract best column location
  bestcol <- which(colnames(popSizeCompo) == "best")
  popSizeQuan <- popSizeCompo%>%gather(., "simulation", "estimate", 
                                       bestcol:ncol(popSizeCompo))%>%
    filter(simulation!= "best")
  
  
  
  if (is.null(Disease_prog) && is.null(Cascade) && is.null(Population)){
    popQ <- popSizeQuan%>%group_by(year)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year)%>% 
      select(everything())
    popSizeCompo[, c("min", "max", "Med", "Mu", 
                     "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
  } 
  else if (is.null(Disease_prog) && is.null(Cascade) && !is.null(Population)){ 
    popQ <- popSizeQuan%>%group_by(year, population)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year, population)%>% 
      select(everything()) 
    
    popSizeCompo[, c("min", "max", "Med", "Mu", 
                     "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95))  
  } 
  else if (is.null(Disease_prog) && !is.null(Cascade) && is.null(Population)) {
    popQ <- popSizeQuan%>%
      group_by(year, cascade)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year, cascade)%>%
      select(everything()) 
    
    popSizeCompo[, c("min", "max", "Med", "Mu", 
                     "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95))  
    
  } 
  else if (!is.null(Disease_prog) && is.null(Cascade) && is.null(Population)) {
    popQ <- popSizeQuan%>%
      group_by(year, disease_prog)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year, disease_prog)%>%
      select(everything()) 
    
    popSizeCompo[, c("min", "max", "Med", "Mu", 
                     "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
  } 
  else if (!is.null(Disease_prog) && !is.null(Cascade) && is.null(Population)) {
    popQ <- popSizeQuan%>%
      group_by(year, state)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year, state)%>%
      select(everything()) 
    
    popSizeCompo[, c("min", "max", "Med", "Mu", 
                     "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95))
  } 
  else if (!is.null(Disease_prog) && is.null(Cascade) && !is.null(Population)) {
    popQ <- popSizeQuan%>%
      group_by(year, disease_prog, population)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year, population, disease_prog)%>%
      select(everything()) 
    
    popSizeCompo[, c("min", "max", "Med", "Mu", 
                     "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
  }
  else if (is.null(Disease_prog) && !is.null(Cascade) && !is.null(Population)) {
    popQ <- popSizeQuan%>%
      group_by(year, cascade, population)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year, population, cascade)%>%
      select(everything()) 
    
    popSizeCompo[, c("min", "max", "Med", "Mu", 
                     "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
  }
  else if (!is.null(Disease_prog) && !is.null(Cascade) && !is.null(Population)) {
    popQ <- popSizeQuan%>%
      group_by(year, population, state)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year, population, state)%>%
      select(everything()) 
    
    popSizeCompo[, c("min", "max", "Med", "Mu", 
                     "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
    
  }
  
  
  return(popSizeCompo)
  
}
  
  
  


### incidence #### 
  indicatorResult_uno<- function(pg, Best, indicator,
                                 populations = NULL, endYear ) {
    
   
    # remove the all pops and newpop in  the best result of calibration
    # becuz different structure of dataset in the list
    df_list <- within(Best, rm(allPops,newpop_tran,death_hcv, HCVdeathState, 
                               newDeathState, costPops,QALYPops))
    
    # convert to data.frame table format 
    df_list <- lapply(df_list, as.data.frame.table)
    
    result_list<- df_list%>%lapply(., function(x) {
      mutate(x, 
             timestep = c(rep(seq(HCV$startYear,
                                  endYear - HCV$timestep,HCV$timestep), 
                              each = HCV$npops)))%>%dplyr::select(-Var2)}) 
    
    names(result_list[[indicator]]) <- c("population", "total_pop", "timestep") 
    
    
    # update the number of first timepoint as the average number in the year
    
    impute<- result_list[[indicator]]%>%filter(timestep >1 &timestep <2)%>%
      group_by(population)%>%
      mutate(total_pop = sum(total_pop)/length(total_pop))%>%
      select(total_pop)%>%ungroup()
    impute <- impute[c(1:HCV$npops),]
    
    result_list[[indicator]][c(1:HCV$npops), "total_pop"] <- impute$total_pop 
    
    
    # ready for appending parameter set 
    indicatorEstimate <- result_list[[indicator]]%>%
      relocate(.,population, timestep, total_pop)
    
    
    if (length(populations) == 1 && populations == "all") {
      indicatorEstimate <- indicatorEstimate%>%
        arrange(population, timestep)%>%
        dplyr::mutate(year = rep(rep(seq(1, endYear-HCV$timestep), 
                                     each=1/HCV$timestep),HCV$npops))%>%
        group_by(year)%>%summarise(best = sum(total_pop))%>%ungroup()
    } else if (is.null(populations)){
      indicatorEstimate<- indicatorEstimate%>%
        arrange(population, timestep)%>%
        dplyr::mutate(year = rep(rep(seq(1, endYear-HCV$timestep), 
                                     each=1/HCV$timestep),HCV$npops))%>%
        group_by(year, population)%>%
        summarise(best = sum(total_pop))%>%ungroup()
    } else if (length(populations) > 1) {
      indicatorEstimate <- indicatorEstimate%>%
        dplyr::filter(population %in% populations)%>%
        arrange(population, timestep)%>%
        dplyr::mutate(year = rep(rep(seq(1, endYear-HCV$timestep), 
                                     each=1/HCV$timestep),length(populations)))%>%
        group_by(year, population)%>%
        summarise(best = sum(total_pop))%>%ungroup()
      
      
    } else{ 
      indicatorEstimate <- indicatorEstimate%>%
        dplyr::filter(population %in% populations)%>%
        arrange(population, timestep)%>%
        dplyr::mutate(year = rep(rep(seq(1, endYear-HCV$timestep), 
                                     each=1/HCV$timestep),1))%>%
        group_by(year)%>%
        summarise(best = sum(total_pop))%>%ungroup()
      
      }
    
    
    
    return(indicatorEstimate)
    
    
  }
  
  
  # function to manipulte data wt/wo parameter sets
  
indicatorResults <- function(pg, Best, indicator, paramR = NULL,
                             pop = NULL, range = NULL, endY, scenario = NULL){
  
  # Best estimate 
  indicatorEstimate <- indicatorResult_uno(pg, Best, indicator,
                                           populations = pop, endYear = endY)
  
  if (!is.null(paramR)){ 
    # apply function to the list of parameter set
    paramlst <- list()
    
    paramlst <- lapply(paramR, function(x) indicatorResult_uno(pg, x, indicator,
                                                               populations = pop, 
                                                               endYear = endY))
    
    #if (length(pop) == 1 && pop == "all"){
    #  paramlst <- lapply(paramlst, function(x) arrange(x, year)) }
    #else{ 
    #  
    #  paramlst <- lapply(paramlst, function(x) arrange(x, year, populations))
    #}
    # make value of each parameter as columns in a dataframe 
    
    paramDf <- as.data.frame(lapply(lapply(paramlst, `[`, "best"), unlist))
    
    indicatorEstimate <- cbind(indicatorEstimate, paramDf) 
    if(is.null(scenario)){ 
      
      samp <- HCV$numberSamples
    }
    else{ samp <- HCV$sParam}
    
    # rename columns 
    if (length(pop) == 1){  
      colnames(indicatorEstimate) <-  c("year","best",
                                        c(paste0("set",
                                                 seq(1,samp, 1 ), 
                                                 sep="")))
      
      
    } else{ 
      colnames(indicatorEstimate) <-  c("year","population","best",
                                        c(paste0("set",
                                                 seq(1,samp, 1 ), 
                                                 sep="")))}
  }
  bestcol <- which(colnames(indicatorEstimate) == "best")
  
  if (is.null(pop)) {
    bestValues <- indicatorEstimate%>%
      group_by(year, population)%>%select(best)
  } else { 
    
    bestValues <- indicatorEstimate%>%group_by(year)%>%select(best)
    
  }
  
  
  if (!is.null(paramR)) {
    
    if (!is.null(range)) {
      indicatorQ <- gather(indicatorEstimate, "sim", "estimate", 
                           bestcol:ncol(indicatorEstimate))%>%
        filter(sim != "best") # aggegrate parameter set only 
      
      if (is.null(pop)|| length(pop)>1) {
        indicatorQuan <- indicatorQ %>%
          group_by(year, population) %>%
          summarise(min = min(estimate, na.rm = TRUE),
                    max = max(estimate, na.rm = TRUE),
                    Med = median(estimate, na.rm = TRUE),
                    Mu = mean(estimate, na.rm = TRUE),
                    q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                    q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                    q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                    q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
          ungroup() %>%
          arrange(year, population) %>%
          mutate(best = bestValues) %>%
          select(year, population, 
                 best, everything())
        
        indicatorEstimate <- indicatorEstimate%>%arrange(year, population)
        
        indicatorEstimate[, c("min", "max", "Med", "Mu", 
                              "q5", "q25", "q75", "q95")] <- 
          indicatorQuan%>%select(., -c("year", "population", "best"))
        
      } else {
        indicatorQuan <- indicatorQ %>%
          group_by(year) %>%
          summarise(min = min(estimate, na.rm = TRUE),
                    max = max(estimate, na.rm = TRUE),
                    Med = median(estimate, na.rm = TRUE),
                    Mu = mean(estimate, na.rm = TRUE),
                    q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                    q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                    q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                    q95 = quantile(estimate, prob = 0.975, na.rm = TRUE)) %>%
          ungroup() %>%arrange(year)%>%
          mutate(best = bestValues) %>%
          select(year, best, everything())
        
        indicatorEstimate <- indicatorEstimate%>%arrange(year)
        indicatorEstimate[, c("min", "max", "Med", "Mu", 
                              "q5", "q25", "q75", "q95")] <- 
          indicatorQuan%>%select(., -c("year", "best"))
      }
    }
    
    
    
  }
  
  return(indicatorEstimate) 
}
  

# function of modifying y-limits 

expandy <- function(data, YE){ 
  # data: dataset  for ggplot 
  # YE: the end of year in the ggplot 
  
  t <- data%>%filter(year <= YE)
  
  max.y = max(t$q95, na.rm=TRUE)
  min.log = floor(log10(max.y))
  
  min.y = min(t$q5, na.rm=TRUE)
  min.log_min = floor(log10(min.y))
  if(max.y!=0 && min.log_min <4){ 
    scale_y_continuous(limits =c(0,ceiling(max.y/10^min.log)*10^min.log))
  } else if (max.y!=0 && min.log_min >=4){ 
    scale_y_continuous(limits =c(1*10^min.log_min,ceiling(max.y/10^min.log)*10^min.log))
  } else {scale_y_continuous(limits =c(-1,1)) }

  }



## not working ## 


################################################################################
#### modify facet wrap axis ####
## reference: https://fishandwhistle.net/post/2018/
##            modifying-facet-scales-in-ggplot2/
## Example: "**" == modified area
##         p + facet_custom(~"*facet_name*", 
##                          scales = "free", ncol = "4", 
##                          scale_overrides = list(
##            scale_override("*1*", scale_y_continuous(
##                                  limits = c(limit_facet(DF)[1], 
##                                          limit_facet(DF)[1+npops])))))
################################################################################

# get dataset to find max and min of y value 
limit_y <- function(data, plotyear){ 
  dt <- data%>%filter(year<plotyear + 1)
  return(dt)
  }



Facet <- ggproto(               
  init_scales = function(layout, x_scale = NULL, y_scale = NULL, params) {
    scales <- list()
    if (!is.null(x_scale)) {
      scales$x <- plyr::rlply(max(layout$SCALE_X), x_scale$clone())
    }
    if (!is.null(y_scale)) {
      scales$y <- plyr::rlply(max(layout$SCALE_Y), y_scale$clone())
    }
    scales
  }
)
scale_new <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

facet_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}

#### finding limits for facet_wrap #### 

limit_facet <- function(DF, YE){ 
  df <- list()
  break_list <- list()
  breaks_min <- list()
  breaks_max <- list()
  pop <- c(unique(DF$population))
  
  for (i in 1: length(unique(DF$population))){ 
    
    df[[i]] <- DF%>%filter(population == pop[i] & year <= YE)
  }
  
  max.y = lapply(df, function(x){max(x$q95, na.rm=TRUE)})
  
  min.log = lapply(max.y, function(x){floor(log10(x))})
  
  min.y = lapply(df, function(x){min(x$q5, na.rm=TRUE)})
  
  min.log_min = lapply(min.y, function(x){floor(log10(x))})
  
  for (i in 1: length(unique(DF$population))) {
    if(max.y[[i]]==0) { 
      breaks_min[[i]] <- -1
      breaks_max[[i]] <- 1
    } else if (min.y[[i]]==0 & max.y[[i]]<=25){ 
      breaks_min[[i]] <- 0
      breaks_max[[i]] <- 30
    } else if (min.y[[i]]==0 & max.y[[i]]<=50){ 
      breaks_min[[i]] <- 0
      breaks_max[[i]] <- 60
    } else if (min.y[[i]]==0 & max.y[[i]]<=120){ 
      breaks_min[[i]] <- 0
      breaks_max[[i]] <- 150
    } else if (min.y[[i]]==0 & max.y[[i]]<=250){ 
      breaks_min[[i]] <- 0
      breaks_max[[i]] <- 300
    } else if (min.y[[i]]==0 & max.y[[i]]<=450){ 
      breaks_min[[i]] <- 0
      breaks_max[[i]] <- 500
    } else if (min.log_min[[i]]!= -Inf &&min.log_min[[i]]>=4) { 
      breaks_min[[i]] <- 10^min.log_min[[i]]
      breaks_max[[i]] <- ceiling(max.y[[i]]/10^min.log[[i]])*10^min.log[[i]]
    } else{ 
      breaks_min[[i]] <- 0
      breaks_max[[i]] <- ceiling(max.y[[i]]/10^min.log[[i]])*10^min.log[[i]]
      
      
    }
    
  } 
  return (c(unlist(breaks_min), unlist(breaks_max)))
  
  
} 

#### PlotOpts
plotOpts <- theme_bw() + theme(plot.title = element_text(face = "bold",
                                                         size = rel(1.2), 
                                                         hjust = 0.5),
                               text = element_text(face = "bold",size=14,
                                                   colour="black"),
                               axis.text.x = element_text(face = "bold", 
                                                          angle = 90,
                                                          size=14,
                                                          colour="black", 
                                                          hjust = 0, vjust = 0.5),
                               axis.text.y = element_text(face = "bold",
                                                          size=14,
                                                          colour="black"),
                               axis.line.x = element_line(colour="black"),
                               axis.line.y = element_line(colour="black"),
                               axis.ticks = element_line(colour="black"),
                               legend.background = element_rect(),
                               legend.key = element_blank(),
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(), 
                               panel.background = element_blank(), 
                               panel.border = element_blank(),
                               axis.line = element_line(colour = "black"),
                               strip.background = element_blank(),
                               plot.margin =unit(c(2.5,2.5,2.5,2.5),"mm"),
                               strip.text = element_text(face="bold")) +
  theme(plot.title = element_text(hjust = 0.5))





theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size) +
      theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(angle=90, vjust = 0.5,hjust = -1,
                                       face = "bold",size = rel(1)), 
            axis.text.y = element_text(
              face = "bold",size = rel(1)),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = rel(1)),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="bold", size= rel(1)),
            plot.margin =unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
      ))
  
}


theme_Publication_facet <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size) + 
      theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border     = element_rect(fill = NA, colour = "black", 
                                            size = NA),
            panel.grid       = element_line(colour = "black"),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(angle=90, vjust = 0.5,hjust = -1,
                                       face = "bold",size = rel(1)), 
            axis.text.y = element_text(
              face = "bold",size = rel(1)),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = rel(1)),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="bold", size= rel(1)),
            plot.margin = unit(c(10,5,5,5),"mm"),
            strip.background = element_blank(),
            strip.text = element_text(face="bold")
            
      ))
  
}
# 
#strip.background = element_blank(),
#strip.text = element_text(face="bold")
#####
 

#### plot generation #### 
# becareful about name of x-axis and y-axis 

indicatorPlot <- function(data,  ylabel = NULL, 
                           rangeun = FALSE,
                           xlimits = NULL, facetPlot = NULL, 
                           calibration_Y,
                           groupPlot = NULL,
                           observationData = NULL, simulateYear){
  
  
  
  # Set defaults
  
  if (is.null(ylabel)) {
    ylabel = indicator
  }
  
  if (is.null(xlimits)) {
    xlimits = c(min(data$year), max(data$year), 20)
  }
  
  if (is.null(facetPlot)) {
    plotFacets = FALSE 
  } else {
    plotFacets = TRUE
  }
  
  if (is.null(groupPlot)) {
    plotGroups = FALSE
  } else {
    plotGroups = TRUE
  }
  
  
  # Create the plot 
  # long form and remove the q25 & q75
  if (is.null(rangeun)){ 
    data<- data%>%filter(year<= xlimits[2])
    
    data[is.na(data)] <- 0 
  } else {
    data_long <- data%>%pivot_longer(
      cols = c(paste0("set", seq(1,HCV$numberSamples,1), sep ="")), 
                                     names_to = "sim", values_to = "val")
    
    }
  # turn na as 0 
  
  
  if (!is.null(observationData)){
    if (!is.null(rangeun)) {
      
      
      if (plotGroups) {
        plot <- ggplot(data = data, aes_string(x = "year", 
                                               group = groupPlot)) +
          
          geom_line(aes(y = best, colour = groupPlot), size = 1) +
          #geom_line(aes(y =Med),  colour = "red", size = 1) + 
          #geom_line(aes_string(y = "q5"),linetype ="dashed") +
          #geom_line(aes_string(y = "q95"),linetype ="dashed") +
          #geom_line(aes_string(y = "val", group = "sim"),  
          #          colour = "gray48", alpha = 0.4) +
          
          geom_ribbon(aes(ymin = q5, ymax = q95, 
                                 fill = groupPlot), alpha = 0.4) +
          xlab("Year") + ylab(ylabel) + 
          plotOpts + 
          geom_point(data=observationData, 
                     aes(y=realPop, x = time), 
                     colour = "black") + 
          geom_segment(data = observationData, 
                       aes ( y = low, yend = up, x = time, xend = time)) 
      } else {
        plot <- ggplot(data = data, aes(x = year)) +
          geom_ribbon(aes(ymin = q5, ymax = q95), 
                                 fill = "skyblue", alpha = 0.4) +
          geom_line(aes(y =best), colour = "blue", size =1) +
          #geom_line(aes(y =Med),  colour = "red", size = 1) + 
          
          #geom_line(aes_string(y = "q5"),linetype ="dashed") +
          #geom_line(aes_string(y = "q95"),linetype ="dashed") +
          #geom_line(aes_string(y = "val", group = "sim"),  
          #          colour = "gray48", alpha = 0.4) +
          
          xlab("Year") + ylab(ylabel) + 
          plotOpts + 
          geom_point(data=observationData, aes(y=realPop, x = time), 
                     colour = "black") +
          geom_segment(data = observationData, 
                       aes ( y = low, yend = up, x = time, xend = time)) }
    
      } else { 
        if (plotGroups) {
        plot <- ggplot(data = data, aes_string(x = "year", 
                                               group = groupPlot)) +
          geom_line(aes(y = best, colour = groupPlot), size =1) +
          xlab("Year") + ylab(ylabel) + 
          plotOpts + 
          
          geom_point(data=observationData, aes(y=realPop, x = time), 
                     colour = "black") +
          geom_segment(data = observationData, 
                       aes ( y = low, yend = up, x = time, xend = time)) 
      } else {
        plot <- ggplot(data = data, aes_string(x = "year")) +
          geom_line(aes(y = best), colour = "blue" , size =1) +
          xlab("Year") + ylab(ylabel) + 
          plotOpts + 
          geom_point(data=observationData, aes(y=realPop, x = time), 
                     colour = "black") +
          geom_segment(data = observationData, 
                       aes ( y = low, yend = up, x = time, xend = time)) }
      
    }
  } else { 
    if (is.null(rangeun)) {
      if (plotGroups) {
        plot <- ggplot(data = data, aes_string(x = "year", 
                                               group = groupPlot)) +
          geom_line(aes(y = best, colour = groupPlot), size =1 ) +
          xlab("Year") + ylab(ylabel)  +
          plotOpts 
      } else {
        plot <- ggplot(data = data, aes_string(x = "year")) +
          geom_line(aes(y = best), colour = "blue") +
          xlab("Year") + ylab(ylabel) +
          
          plotOpts 
        }
    }
    else { 
      if (plotGroups) {
        plot <- ggplot(data = data, aes_string(x = "year", 
                                               group = groupPlot)) +
          geom_ribbon(aes(ymin = q5, ymax = q95, 
                                 fill = groupPlot), alpha = 0.4) +
          geom_line(aes(y = best, colour = groupPlot)) +
          geom_line(aes(y =Med),  colour = "red", size = 1) + 
          
          #geom_line(aes_string(y = "q5"),linetype ="dashed") +
          #geom_line(aes_string(y = "q95"),linetype ="dashed") +
          #geom_line(aes_string(y = "val", group = "sim"),  
          #          colour = "gray48", alpha = 0.4) + 
         
          xlab("Year") + ylab(ylabel) + 
          plotOpts 
      } else {
        plot <- ggplot(data = data, aes(x = year)) +
          geom_line(aes(y = best), colour = "blue", size =1) +
          geom_ribbon(aes(ymin = q5, ymax = q95), fill = "skyblue", 
                      alpha = 0.4) +
          #geom_line(aes(y =Med),  colour = "red", size = 1) + 
          
          #geom_line(aes_string(y = "q5"),linetype ="dashed") +
          #geom_line(aes_string(y = "q95"),linetype ="dashed") +
          #geom_line(aes_string(y = "val", group = "sim"),  
          #          colour = "gray48", alpha = 0.4) +
          xlab("Year") + ylab(ylabel) + 
          plotOpts  
        }
    }
  }
  
  
  
  
  plot <- plot + coord_cartesian(xlim = xlimits[1:2]) +
    scale_x_continuous(limits =c(xlimits[1], xlimits[2]) ,
                       breaks = seq(xlimits[1], xlimits[2], 
                                    by = xlimits[3]),
                       labels = seq(calibration_Y, 
                                    (calibration_Y+xlimits[2]-xlimits[1]),
                                    xlimits[3])) + 
    geom_vline(xintercept=simulateYear, 
               linetype = "dashed", size = 1)
  
  if (plotFacets) { 
    
    
    # converrt variabel to string for facet_wrap function
    facetPlot <- deparse(substitute(facetPlot))
    
    plot <- plot + facet_wrap(as.formula(paste0("~", facetPlot, collapse = "")),
                              ncol = 2, scales = "free") +
      theme(panel.spacing = unit(2, "lines"))
  }
  return(plot)
} 


  
#### Involving Scenarios ####
# not finished yet #
# The functions for generation results and plots for scenarios are similar to 
# the functions for calibration. I chose to make it as a seperate function 
# becuz easy to modify and increase readability. 
# However, this needs to improve for flexibility. 


# functions combined results of scenarios and best result 

indicatorResultsScen <- function(pg, Best, indicator, ScenaR = NULL,
                                 pop = NULL, range = NULL, endY){
  
  # Best estimate 
  indicatorEstimate <- indicatorResult_uno(pg, Best, indicator,
                                           populations = pop, endYear = endY)
  
  if (!is.null(ScenaR)){ 
    # apply function to the list of parameter set
    paramlst <- list()
    
    paramlst <- lapply(ScenaR, function(x) indicatorResult_uno(pg, x, indicator,
                                                               populations = pop, 
                                                               endYear = endY))
    
    #if (length(pop) == 1 && pop == "all"){
    #  paramlst <- lapply(paramlst, function(x) arrange(x, year)) }
    #else{ 
    #  
    #  paramlst <- lapply(paramlst, function(x) arrange(x, year, populations))
    #}
    # make value of each parameter as columns in a dataframe 
    
    paramDf <- as.data.frame(lapply(lapply(paramlst, `[`, "best"), unlist))
    
    indicatorEstimate <- cbind(indicatorEstimate, paramDf)
    
    # rename columns 
    #if (length(pop) == 1){  
    #  colnames(indicatorEstimate) <-  c("year","best",
    #                                    c("POCab", "DBS", "Reflex RNA", "POCRNA"))
      
      
    #} else{ 
    #  colnames(indicatorEstimate) <-  c("year","population","best",
    #                                    c("POCab", "DBS", "Reflex RNA", "POCRNA"))
    #}
  }
  
  
  
  return(indicatorEstimate) 
}

# extract data from scenariosResults - indicators 


Scenario_R <- function( SRdlst, indicator){ 
  # SRdlst: scenariosResult
  # dLst: extraction of coverage 
  # indicators: the result in the SRdlst
  SceName <- names(SRdlst)
  Scenarios_R <- list()
  for (i in 1:length(SceName)){
    Scenarios_R[[i]] <- lapply(SRdlst[[i]], function(x) 
      indicatorResults(HCV, x, indicator, 
                     pop=HCV$popNames, paramR = NULL ,range = "y",
                     endY = 50))
  }
  
  return(Scenarios_R)
} 



popResultsSce_MidYear <- function(pg, Best, Population = NULL, 
                               Disease_prog = NULL, Cascade = NULL, 
                               MidYear = NULL, 
                               endYear ,
                               YearCut = NULL) {
  # This function organizes and merges the population size results 
  # example 
  #           popResults(HCV, best, Population = c("HIV-"), Disease_prog = NULL, 
  #           Cascade = c("undiag", "diag_ab"))
  # 
  # NULL means not seperate by variable 
  # time = timestep
  
  
  
  
  # First organize the best estimates  
  df_list <- lapply(Best, as.data.frame.table)
  allpop <- df_list$allPops%>%
    mutate(time = rep(seq(1.0,(endYear - pg$timestep), pg$timestep), 
                      each=pg$ncomponent*pg$npops),
           cascade = sub("^[^_]*_", "", Var2), 
           disease_prog = sub("\\_.*", "", Var2))%>%
    dplyr::select(-Var3)
  
  names(allpop) <- c("population", "state", "Frequency","timestep", "cascade",
                     "disease_prog")
  
  ## MidyearIndex
  timelong <- seq(pg$startYear, endYear, pg$timestep) 
  if (is.null(YearCut)) {
    MidY <-c(timelong[MidyearIndex(length(timelong),0.1)]) 
    allpop <- allpop%>%filter(timestep%in% MidY)%>%
      mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                        each = pg$npops*pg$ncomponent))
    
  } else if (YearCut =="Start") { 
    
    MidY <-c(timelong[StryearIndex(length(timelong),0.1)]) 
    allpop <- allpop%>%filter(timestep%in% MidY)%>%
      mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                        each = pg$npops*pg$ncomponent))
    
  } else {  
    MidY <-c(timelong[EndyearIndex(length(timelong),0.1)]) 
    allpop <- allpop%>%filter(timestep%in% MidY)%>%
      mutate(year = rep(seq(pg$startYear, (endYear-1),1), 
                        each = pg$npops*pg$ncomponent))
    
    
  }
  
  
  
  #  # now prepare to append parameter sets 
  popSizes <- allpop%>%mutate(best = Frequency)%>%
    select(-c(Frequency))
  
  return(popSizes)
}
  
 

indicatorPlotSce <- function(data,  ylabel = NULL, 
                          rangeun = FALSE,
                          xlimits = NULL, facetPlot = NULL, 
                          calibration_Y,
                          groupPlot = NULL,
                          observationData = NULL, simulateYear, 
                          Scen = NULL){
  
  
  
  # Set defaults
  
  if (is.null(ylabel)) {
    ylabel = indicator
  }
  
  if (is.null(xlimits)) {
    xlimits = c(min(data$year), max(data$year), 20)
  }
  
  if (is.null(facetPlot)) {
    plotFacets = FALSE 
  } else {
    plotFacets = TRUE
  }
  
  if (is.null(groupPlot)) {
    plotGroups = FALSE
  } else {
    plotGroups = TRUE
  }
  
  
  # Create the plot 
  # long form and remove the q25 & q75
  if (is.null(rangeun)){ 
    data<- data%>%filter(year<= xlimits[2])
    
    data[is.na(data)] <- 0 
  } else {
    data_long <- data%>%pivot_longer(
      cols = c(paste0("set", seq(1,HCV$numberSamples,1), sep ="")), 
      names_to = "sim", values_to = "val")
    
    
  }
  # turn na as 0 
  
  
  if (!is.null(observationData)){
    if (!is.null(rangeun)) {
      
      
      if (plotGroups) {
        plot <- ggplot(data = data_long, aes_string(x = "year", 
                                                    group = groupPlot)) +
          geom_ribbon(aes_string(ymin = "min", ymax = "max", 
                                 fill = groupPlot), alpha = 0.4) +
          geom_line(aes_string(y = "best", colour = groupPlot), size = 1) +
          geom_line(aes_string(y ="Mu"),  colour = "red", size = 1) + 
          #geom_line(aes_string(y = "q5"),linetype ="dashed") +
          #geom_line(aes_string(y = "q95"),linetype ="dashed") +
          #geom_line(aes_string(y = "val", group = "sim"),  
          #          colour = "gray48", alpha = 0.4) +
          
          geom_ribbon(aes_string(ymin = "q5", ymax = "q95", 
                                 fill = groupPlot), alpha = 0.4) +
          xlab("Year") + ylab(ylabel) + 
          plotOpts + 
          geom_point(data=observationData, 
                     aes(y=realPop, x = time), 
                     colour = "black") + 
          geom_segment(data = observationData, 
                       aes ( y = low, yend = up, x = time, xend = time)) 
      } else {
        plot <- ggplot(data = data_long, aes(x = year)) +
          geom_ribbon(aes(ymin = q5, ymax = q95), 
                      fill = "skyblue", alpha = 0.4) +
          geom_line(aes(y =best), colour = "blue", size =1) +
          geom_line(aes_string(y ="Mu"),  colour = "red", size = 1) + 
          
          #geom_line(aes_string(y = "q5"),linetype ="dashed") +
          #geom_line(aes_string(y = "q95"),linetype ="dashed") +
          #geom_line(aes_string(y = "val", group = "sim"),  
          #          colour = "gray48", alpha = 0.4) +
          
          xlab("Year") + ylab(ylabel) + 
          plotOpts + 
          geom_point(data=observationData, aes(y=realPop, x = time), 
                     colour = "black") +
          geom_segment(data = observationData, 
                       aes ( y = low, yend = up, x = time, xend = time)) }
      
    } else { 
      if (plotGroups) {
        plot <- ggplot(data = data, aes_string(x = "year", 
                                               group = groupPlot)) +
          geom_line(aes_string(y = "best", colour = groupPlot), size =1) +
          xlab("Year") + ylab(ylabel) + 
          plotOpts + 
          
          geom_point(data=observationData, aes(y=realPop, x = time), 
                     colour = "black") +
          geom_segment(data = observationData, 
                       aes ( y = low, yend = up, x = time, xend = time)) 
      } else {
        plot <- ggplot(data = data, aes_string(x = "year")) +
          geom_line(aes_string(y = "best"), colour = "blue" , size =1) +
          xlab("Year") + ylab(ylabel) + 
          plotOpts + 
          geom_point(data=observationData, aes(y=realPop, x = time), 
                     colour = "black") +
          geom_segment(data = observationData, 
                       aes ( y = low, yend = up, x = time, xend = time)) }
      
    }
  } else { 
    if (is.null(rangeun)) {
      if (plotGroups) {
        plot <- ggplot(data = data, aes_string(x = "year", 
                                               group = groupPlot)) +
          geom_line(aes_string(y = "best", colour = groupPlot), size =1 ) +
          xlab("Year") + ylab(ylabel)  +
          plotOpts 
      } else {
        plot <- ggplot(data = data, aes_string(x = "year")) +
          geom_line(aes_string(y = "best"), colour = "blue") +
          xlab("Year") + ylab(ylabel) +
          
          plotOpts 
      }
    }
    else { 
      if (plotGroups) {
        plot <- ggplot(data = data_long, aes_string(x = "year", 
                                                    group = groupPlot)) +
          geom_ribbon(aes_string(ymin = "q5", ymax = "q95", 
                                 fill = groupPlot), alpha = 0.4) +
          geom_line(aes_string(y = "best", colour = groupPlot)) +
          geom_line(aes_string(y ="Mu"),  colour = "red", size = 1) + 
          
          #geom_line(aes_string(y = "q5"),linetype ="dashed") +
          #geom_line(aes_string(y = "q95"),linetype ="dashed") +
          #geom_line(aes_string(y = "val", group = "sim"),  
          #          colour = "gray48", alpha = 0.4) + 
          
          xlab("Year") + ylab(ylabel) + 
          plotOpts 
      } else {
        plot <- ggplot(data = data_long, aes(x = year)) +
          geom_line(aes(y = best), colour = "blue", size =1) +
          geom_ribbon(aes(ymin = q5, ymax = q95), fill = "skyblue", 
                      alpha = 0.4) +
          geom_line(aes_string(y ="Mu"),  colour = "red", size = 1) + 
          
          #geom_line(aes_string(y = "q5"),linetype ="dashed") +
          #geom_line(aes_string(y = "q95"),linetype ="dashed") +
          #geom_line(aes_string(y = "val", group = "sim"),  
          #          colour = "gray48", alpha = 0.4) +
          xlab("Year") + ylab(ylabel) + 
          plotOpts  
      }
    }
  }
  
  
  
  
  plot <- plot + coord_cartesian(xlim = xlimits[1:2]) +
    scale_x_continuous(limits =c(xlimits[1], xlimits[2]) ,
                       breaks = seq(xlimits[1], xlimits[2], 
                                    by = xlimits[3]),
                       labels = seq(calibration_Y, 
                                    (calibration_Y+xlimits[2]-xlimits[1]),
                                    xlimits[3])) + 
    geom_vline(xintercept=simulateYear-calibration_Y +1, 
               linetype = "dashed", size = 1)
  
  if (plotFacets) { 
    
    
    # converrt variabel to string for facet_wrap function
    facetPlot <- deparse(substitute(facetPlot))
    
    plot <- plot + facet_wrap(as.formula(paste0("~", facetPlot, collapse = "")),
                              ncol = 2, scales = "free") +
      theme(panel.spacing = unit(2, "lines"))
  }
  return(plot)
} 




#### finding limits for facet_wrap: manuscript plots#### 
expandy_man <- function(data, YE){ 
  # data: dataset  for ggplot 
  # YE: the end of year in the ggplot 
  
  t <- data%>%filter(year <= YE)
  
  max.y = max(t, na.rm=TRUE)
  min.log = floor(log10(max.y))
  
  min.y = min(t, na.rm=TRUE)
  min.log_min = floor(log10(min.y))
  if(max.y!=0 && min.log_min <4){ 
    scale_y_continuous(limits =c(0,ceiling(max.y/10^min.log)*10^min.log))
  } else if (max.y!=0 && min.log_min >=4){ 
    scale_y_continuous(limits =c(1*10^min.log_min,ceiling(max.y/10^min.log)*10^min.log))
  } else {scale_y_continuous(limits =c(-1,1)) }
  
}






limit_facet_man <- function(DF, YE){ 
  df <- list()
  break_list <- list()
  breaks_min <- list()
  breaks_max <- list()
  pop <- c(unique(DF$population))
  
  for (i in 1: length(unique(DF$population))){ 
    
    df[[i]] <- DF%>%filter(population == pop[i] & year< YE)
  }
  
  max.y = lapply(df, function(x){max(x$best, na.rm=TRUE)})
  min.log = lapply(max.y, function(n){floor(log10(n))})
  

  for (i in 1: length(unique(DF$population))) {
      
      breaks_min[[i]] <- 0
      breaks_max[[i]] <- ceiling(max.y[[i]]/10^min.log[[i]])*10^min.log[[i]]
    
    
  } 
  return (c(unlist(breaks_min), unlist(breaks_max)))
  
  
} 



