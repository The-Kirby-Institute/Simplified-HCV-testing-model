#### function for plots generation in manuscripts #### 
# colorblind friendly palatte 
# remotes::install_github("clauswilke/colorblindr")
library(colorblindr)
library(gridExtra)
library(grid)
library(ggplot2)
library(xtable)
library(ggpubr)
library(here())
basePath <- here()


#### tidy up the datasets for plot #### 
# prevalence 
scenario_Prevalence <- function(srlst, sparamlst, pop = NULL, statusQ = NULL){ 
  if(is.null(statusQ)){
    nn <- names(srlst)
  }else{ 
    nn <- "Status Quo"
    
    }
  
  if(is.null(pop)){ 
    StempNOTInfected_all <- list()
    
    StempNOTInfected_all <- 
      lapply(seq_along(srlst), function(x){
        
        a <- popResults_MidYear(HCV, srlst[[x]],
                                Population = NULL,Disease_prog = NULL,
                                Cascade = c("s", "cured"), 
                                param = sparamlst[[x]],
                                endYear = 30, scenario = "Main")%>%
          as.data.frame()%>%group_by(year)%>%
          summarise_at(.vars = c(colnames(.)[-c(1,2)]),sum)
      })
    
    
    
  
    # number of all population 
    StempTotal_all <- list()
    StempTotal_all <- 
      lapply(seq_along(srlst), function(x){
        
        a <- popResults_MidYear(HCV, srlst[[x]],
                                Population = NULL,
                                Disease_prog = NULL, 
                                Cascade = NULL, 
                                param = sparamlst[[x]], 
                                endYear = 30, scenario = "Main")
      })
    # finding columns with na value 
    
      StempPrev_all <- list()
      StempPrev_all <-  lapply(seq_along(srlst), function(x){
        a <-cbind(year = seq(HCV$startYear , 30-1 ,1),
                  as.data.frame((StempTotal_all[[x]][, -1] - 
                                 StempNOTInfected_all[[x]][ ,-1])/ 
                                  StempTotal_all[[x]][ ,-1]*100))%>%
        tibble::as_tibble() })
    
    
    
    # summarise q5/q95 and median for prevalence 
    
      StempPrev_all_range <- list() 
      StempPrev_all_range <- lapply(seq_along(StempPrev_all), function(x){
        a <-popResults_range(HCV, StempPrev_all[[x]], Population = NULL,
                             Disease_prog = NULL , 
                             Cascade = NULL, end_Y = 30)%>%as.data.frame()%>%
          mutate(testing = nn[x])%>%
          select(year, testing, best, Mu, Med, max, min, q5, q25, q75, q95 )
    }) 
    
    
  }
  else{ 
    StempNOTInfected_all <- list()
    StempNOTInfected_all <- 
      lapply(seq_along(srlst), function(x){
        
        a <- popResults_MidYear(HCV, srlst[[x]],
                                Population = HCV$popNames, Disease_prog = NULL,
                                Cascade = c("s", "cured"), 
                                param = sparamlst[[x]],
                                endYear = 30, scenario = "Main")%>%
          as.data.frame()%>%ungroup()%>%group_by(year, population)%>%
          summarise_at(.vars = c(colnames(.)[-c(1,2,3)]),sum)
      })
    
    StempTotal_all <- list()
    StempTotal_all <- 
      lapply(seq_along(srlst), function(x){
        
        a <- popResults_MidYear(HCV, srlst[[x]],
                                Population = HCV$popNames,
                                Disease_prog = NULL, 
                                Cascade = NULL, 
                                param = sparamlst[[x]], 
                                endYear = 30, scenario = "Main")
      })
      
    StempPrev_all <- list()
      
    StempPrev_all <-  lapply(seq_along(srlst), function(x){
        
        a <-cbind(year = rep(seq(HCV$startYear , 30-1 ,1), each = HCV$npops),
                  population = HCV$popNames,
                  as.data.frame((StempTotal_all[[x]][, -c(1,2)] - 
                                   StempNOTInfected_all[[x]][ ,-c(1,2)])/ 
                                  StempTotal_all[[x]][ ,-c(1,2)]*100))%>%
          tibble::as_tibble() })
      
      # summarise q5/q95 and median for prevalence 
      StempPrev_all_range <- list() 
      StempPrev_all_range <- lapply(seq_along(StempPrev_all), function(x){
        
        a <-popResults_range(HCV, StempPrev_all[[x]], Population = HCV$popNames,
                             Disease_prog = NULL , 
                             Cascade = NULL, end_Y = 30)%>%
          as.data.frame()%>%
          mutate(testing = nn[x])%>%
          select(year, population, testing, best, Mu, Med, max, min, q5, q25, 
                 q75, q95 )
    })
      
    
  
   
    }
    
    StempPrev_all_range <- do.call("rbind", StempPrev_all_range)
  
    return(StempPrev_all_range) 
  }
# number living with chronic HCV
scenario_CHCVnum <- function(srlst, sparamlst, pop = NULL, statusQ = NULL){ 
  if(is.null(statusQ)){
    nn <- names(srlst)
  }else{ 
    nn <- "Status Quo"
    
  }
  
  StempTotal_all <- list()
  StempTotal_all <- 
    lapply(seq_along(srlst), function(x){
      
      a <- popResults_MidYear(HCV, srlst[[x]],
                              Population = NULL,
                              Disease_prog = NULL, 
                              Cascade = NULL, 
                              param = sparamlst[[x]], 
                              endYear = 30, scenario = "Main")
    })
  
  tempTotal_all_nadel <- list()
  tempTotal_all_nadel <- lapply(seq_along(StempTotal_all), function(x){
    names(which(colSums(is.na(StempTotal_all[[x]]))>0))%>%c()})
  
  
  if(is.null(pop)){ 
    Sce_statepop <- list()
    Sce_statepop <- 
      lapply(seq_along(srlst), function(x){
        
        a <- popResults_MidYear(HCV, srlst[[x]],
                                Population = NULL,
                                Disease_prog = c("s",HCV$progress_name), 
                                Cascade = c("s",HCV$cascade_name),
                                param = sparamlst[[x]],
                                endYear = 30, scenario = "Main")%>%
          as.data.frame()%>%mutate(cascade = sub("^[^_]*_", "", state),
                        disease_prog = sub("\\_.*", "", state))%>%
          filter(!cascade%in%c("s", "cured") & !disease_prog%in%c("a"))%>%
          group_by(year)%>%summarise(across(
            c(best,paste0("set", seq(1, HCV$sParam,1))), ~ sum(.x, na.rm = TRUE)))%>%
          as.data.frame()%>%
          dplyr::select(!c(tempTotal_all_nadel[[x]]))
      })
    
    
    
    Sce_statepopQ <- list()
    Sce_statepopQ <- lapply(seq_along(Sce_statepop), function(x){ 
      
      
      a <- popResults_range(HCV, Sce_statepop[[x]], Population = NULL,
                            Disease_prog = NULL , 
                            Cascade = NULL, end_Y = 30)%>%
        as.data.frame()%>%
        mutate(testing = nn[x])%>%
        select(year, testing, best, Mu, Med, max, min, q5, q25, q75, q95 )
    })
    
    Sce_chcQ <- do.call("rbind", Sce_statepopQ)
    
    
  }
  else{ 
    Sce_statepop <- list()
    Sce_statepop <- 
      lapply(seq_along(srlst), function(x){
        
        a <- popResults_MidYear(HCV, srlst[[x]],
                                Population = HCV$popNames,
                                Disease_prog = c("s",HCV$progress_name), 
                                Cascade = c("s",HCV$cascade_name),
                                param = sparamlst[[x]],
                                endYear = 30, scenario = "Main")%>%
          as.data.frame()
        
        
        a <- a%>%mutate(cascade = sub("^[^_]*_", "", state),
                        disease_prog = sub("\\_.*", "", state))%>%
          filter(!cascade%in%c("s", "cured") & !disease_prog%in%c("a"))%>%
          group_by(year, population)%>%summarise(across(
            c(best,paste0("set", seq(1, HCV$sParam,1))), ~ sum(.x, na.rm = TRUE)))%>%
          as.data.frame()%>%dplyr::select(!c(tempTotal_all_nadel[[x]]))
      })
    
    
    Sce_statepopQ <- list()
    Sce_statepopQ <- lapply(seq_along(Sce_statepop), function(x){ 
      
      a <- popResults_range(HCV, Sce_statepop[[x]], Population = HCV$popNames,
                            Disease_prog = NULL , 
                            Cascade = NULL, end_Y = 30)%>%
        as.data.frame()%>%
        mutate(testing = nn[x])%>%
        select(year,population, testing, best, Mu, Med, max, min, q5, q25, q75, q95 )
    })
    
    Sce_chcQ <- do.call("rbind", Sce_statepopQ)
  }
  
  
  return(Sce_chcQ) 
}
# incidence
scenario_Incidence  <- function(srlst, sparamlst, pop = NULL, 
                                indicator = NULL, statusQ = NULL){ 
  if(is.null(statusQ)){
    nn <- names(srlst)
  }else{ 
    nn <- "Status Quo"
    
  }
  
  totalpop <- list()
  totalpop <- lapply(seq_along(srlst), function(x) { 
    
    a <- popResults_MidYear(HCV, srlst[[x]], Population = NULL,
                            Disease_prog = NULL , 
                            Cascade = NULL, param = sparamlst[[x]], 
                            endYear = 30, 
                            scenario = "main")%>%
      as.data.frame()
  }) 
  
  
 
  
  
  if(is.null(indicator)){
    if(is.null(pop)){ 
      
      popS <- list()
      popS <- lapply(seq_along(srlst), function(x) { 
        
        a <- popResults_MidYear(HCV, srlst[[x]], Population = NULL,
                                Disease_prog = NULL , 
                                Cascade = NULL, param = sparamlst[[x]], 
                                endYear = 30, 
                                scenario = "main")%>%
          as.data.frame()
      }) 
      
      
      
      HCVInfect_all <- list()
      HCVInfect_all <- lapply(seq_along(srlst), function(x){ 
        
        a <- indicatorResults(HCV, srlst[[x]], "newInfections", 
                              pop="all",
                              paramR = sparamlst[[x]], 
                              range = "y",
                              endY = 30, scenario = "main")
      })
      
      # reorder columns and only keep best and parameterset 
      colset <- HCVInfect_all[[1]]%>%select(., contains("Set"))%>%
        select(., str_sort(paste0("set", seq(1, HCV$sParam,1)),numeric = FALSE))
      
      HCVInfect <- lapply(HCVInfect_all, function(x){ 
        
        
        a <- x%>%select(year, best, colnames(colset))
      })
      
      
      HCVInc <- list()
      HCVInc <- lapply(seq_along(HCVInfect), function(x){ 
        a <- cbind(year = seq(HCV$startYear , 30-1 ,1),
                   as.data.frame(HCVInfect[[x]][, -1] / 
                                   popS[[x]][ ,-1]*1000))%>%
          tibble::as_tibble()
      })
      
      # get aggregate columns 
      
      HCVInc_range <- list()
      HCVInc_range <- lapply(seq_along(HCVInc), function(x){ 
        
        a <- popResults_range(HCV, HCVInc[[x]], Population = NULL,
                              Disease_prog = NULL , 
                              Cascade = NULL, end_Y = 30)%>%
          as.data.frame()%>%
          mutate(testing = nn[x])%>%
          select(year, testing, best, Mu, Med, max, min, q5, q25, q75, q95 )
      })
      
      StempInc_all_range <- do.call("rbind", HCVInc_range) 
      
      
      
    }
    else{ 
      popS <- list()
      
      
      popS <- lapply(seq_along(srlst), function(x) { 
        
        a <- popResults_MidYear(HCV, srlst[[x]], Population = HCV$popNames,
                                Disease_prog = NULL , 
                                Cascade = NULL, param = sparamlst[[x]], 
                                endYear = 30, 
                                scenario = "main")%>%as.data.frame()
        
        
      })
      
      
      HCVInfect_all <- list()
      HCVInfect_all <- lapply(seq_along(srlst), function(x){ 
        
        a <- indicatorResults(HCV, srlst[[x]], "newInfections", 
                              pop=HCV$popNames,
                              paramR = sparamlst[[x]], 
                              range = "y",
                              endY = 30, scenario = "main")
      })
      
      # reorder columns and only keep best and parameterset 
      colset <- HCVInfect_all[[1]]%>%select(., contains("Set"))%>%
        select(., str_sort(paste0("set", seq(1, HCV$sParam,1)),numeric = FALSE))
      
      HCVInfect <- lapply(HCVInfect_all, function(x){ 
        
        
        a <- x%>%select(year, population, best, colnames(colset))
      })
      
      
      HCVInc <- list()
      HCVInc <- lapply(seq_along(HCVInfect), function(x){ 
        a <- cbind(year = rep(seq(HCV$startYear , 30-1 ,1), each = HCV$npops),
                   population = HCV$popNames,
                   as.data.frame(HCVInfect[[x]][, -c(1,2)] / 
                                   popS[[x]][ ,-c(1,2)]*1000))%>%
          tibble::as_tibble()
      })
      
      # get aggregate columns 
      
      HCVInc_range <- list()
      HCVInc_range <- lapply(seq_along(HCVInc), function(x){ 
        
        a <- popResults_range(HCV, HCVInc[[x]], Population = HCV$popNames,
                              Disease_prog = NULL , 
                              Cascade = NULL, end_Y = 30)%>%
          as.data.frame()%>%
          mutate(testing = nn[x])%>%
          select(year, population, testing, best, Mu, Med, max, min, q5, 
                 q25, q75, q95 )
      })
      
      StempInc_all_range <- do.call("rbind", HCVInc_range) 
    }
    return(StempInc_all_range)
  }
  else if(indicator =="cum"){
    if(is.null(pop)){ 
      
      HCVInfect_all <- list()
      HCVInfect_all <- lapply(seq_along(srlst), function(x){ 
        
        a <- indicatorResults(HCV, srlst[[x]], "newInfections", 
                              pop="all",
                              paramR = sparamlst[[x]], 
                              range = "y",
                              endY = 30, scenario = "main")%>%
          mutate(testing = nn[x])
        
        
        a <-a%>%mutate(across(c(best, paste0("set", seq(1, HCV$numberSamples,1))),
                            cumsum, .names = "set{col}"))%>%
          select(year, testing, c(best, 
                                  paste0("set", seq(1, HCV$numberSamples,1))))
        
        a <-popResults_range(HCV, a, Population = NULL,
                         Disease_prog = NULL , 
                         Cascade = NULL, end_Y = 30)%>%
          as.data.frame()%>%
          select(year, testing, best, Mu, Med, max, min, q5, 
                 q25, q75, q95 )
      })
      
      StempInc_all_range <- do.call("rbind", HCVInfect_all)
      
    }
    else{ 
      HCVInfect_all <- list()
      HCVInfect_all <- lapply(seq_along(srlst), function(x){ 
        
        a <- indicatorResults(HCV, srlst[[x]], "newInfections", 
                              pop=HCV$popNames,
                              paramR = sparamlst[[x]], 
                              range = "y",
                              endY = 30, scenario = "main")%>%
          mutate(testing = nn[x])
        a <-a%>%mutate(across(c(best, paste0("set", seq(1, HCV$numberSamples,1))),
                              cumsum, .names = "set{col}"))%>%
          select(year, population, testing, 
                 c(best, paste0("set", seq(1, HCV$numberSamples,1))))
        
        a <-popResults_range(HCV, a, Population = HCV$popNames,
                             Disease_prog = NULL , 
                             Cascade = NULL, end_Y = 30)%>%
          as.data.frame()%>%
          select(year, population, testing, best, Mu, Med, max, min, q5, 
                 q25, q75, q95 )
        
        
        
        
        
        })
      
      StempInc_all_range <- do.call("rbind", HCVInfect_all) 
    }
    return(StempInc_all_range)
    
  }
  else{ 
    if(is.null(pop)){ 
      
      HCVInfect_all <- list()
      HCVInfect_all <- lapply(seq_along(srlst), function(x){ 
        
        a <- indicatorResults(HCV, srlst[[x]], indicator, 
                              pop="all",
                              paramR = sparamlst[[x]], 
                              range = "y",
                              endY = 30, scenario = "main")%>%
          mutate(testing = nn[x])%>%
          select(year, testing, best, Mu, Med, max, min, q5, 
                                        q25, q75, q95 )
      })
      
      StempInc_all_range <- do.call("rbind", HCVInfect_all)
      
    }
    else{ 
      HCVInfect_all <- list()
      HCVInfect_all <- lapply(seq_along(srlst), function(x){ 
        
        a <- indicatorResults(HCV, srlst[[x]], indicator, 
                              pop=HCV$popNames,
                              paramR = sparamlst[[x]], 
                              range = "y",
                              endY = 30, scenario = "main")%>%
          mutate(testing = nn[x])%>%
          select(year, population, testing, best, Mu, Med, max, min, q5, 
                 q25, q75, q95 )
      })
      
      StempInc_all_range <- do.call("rbind", HCVInfect_all) 
    }
    return(StempInc_all_range)
  }
  
  return(StempInc_all_range)
}

scenario_cumInci  <- function(srlst, sparamlst, pop = NULL, statusQ = NULL){ 
  if(is.null(statusQ)){
    nn <- names(srlst)
  }else{ 
    nn <- "Status Quo"
    
  }
    if(is.null(pop)){ 
      
      HCVInfect_all <- list()
      HCVInfect_all <- lapply(seq_along(srlst), function(x){ 
        
        a <- indicatorResults(HCV, srlst[[x]], "newInfections", 
                              pop="all",
                              paramR = sparamlst[[x]], 
                              range = "y",
                              endY = 30, scenario = "main")%>%
          mutate(testing = nn[x])
        
        
        a <-a%>%mutate(across(c(best, paste0("set", seq(1, HCV$numberSamples,1))),
                              cumsum, .names = "set{col}"))%>%
          select(year, testing, c(best, paste0("set", seq(1, HCV$numberSamples,1))))
        
        a <-popResults_range(HCV, a, Population = NULL,
                             Disease_prog = NULL , 
                             Cascade = NULL, end_Y = 30)%>%
          as.data.frame()
      })
      
      StempInc_all_range <- do.call("rbind", HCVInfect_all)
      
    }
    else{ 
      HCVInfect_all <- list()
      HCVInfect_all <- lapply(seq_along(srlst), function(x){ 
        
        a <- indicatorResults(HCV, srlst[[x]], "newInfections", 
                              pop=HCV$popNames,
                              paramR = sparamlst[[x]], 
                              range = "y",
                              endY = 30, scenario = "main")%>%
          mutate(testing = nn[x])
        a <-a%>%mutate(across(c(best, paste0("set", seq(1, HCV$numberSamples,1))),
                              cumsum, .names = "set{col}"))%>%
          select(year, population, testing, 
                 c(best, paste0("set", seq(1, HCV$numberSamples,1))))
        
        a <-popResults_range(HCV, a, Population = HCV$popNames,
                             Disease_prog = NULL , 
                             Cascade = NULL, end_Y = 30)%>%
          as.data.frame()
        
        
        
        
        
      })
      
      StempInc_all_range <- do.call("rbind", HCVInfect_all) 
    }
    return(StempInc_all_range)
    
  }
 

#### outcome combine ####

Restbl <- function(srlst, sparamlst, pop = NULL, statusQ = NULL){
  rtbl <- list()
  if(is.null(statusQ)){
    if(is.null(pop)){
      scenPrev <- scenario_Prevalence(srlst, sparamlst, pop = NULL)%>%
        mutate(indicator = "Prevalence")
      
      scenInci <- scenario_Incidence(srlst, sparamlst, pop = NULL)%>%
        mutate(indicator = "Incidence")
      
      scenCHCVnum <- scenario_CHCVnum(srlst, sparamlst, pop = NULL)%>%
        mutate(indicator = "CHCVNum")
      
      scenIncinum <- scenario_Incidence(srlst, sparamlst, pop = NULL,
                                        indicator = "newInfections")%>%
        mutate(indicator = "InciNum")
      
      scenHCVdeath <- scenario_Incidence(srlst,sparamlst, pop = NULL, 
                                         indicator = "newHCVdeaths")%>%
        mutate(indicator = "HCVdeath")
    }
    else{ 
      scenPrev <- scenario_Prevalence(srlst, sparamlst, pop = HCV$popNames)%>%
        mutate(indicator = "Prevalence")
      
      scenInci <- scenario_Incidence(srlst, sparamlst, pop = HCV$popNames)%>%
        mutate(indicator = "Incidence")
      
      scenCHCVnum <- scenario_CHCVnum(srlst, sparamlst, pop = HCV$popNames)%>%
        mutate(indicator = "CHCVNum")
      
      scenIncinum <- scenario_Incidence(srlst, sparamlst, pop = HCV$popNames,
                                        indicator = "newInfections")%>%
        mutate(indicator = "InciNum")
      
      scenHCVdeath <- scenario_Incidence(srlst,sparamlst, pop = HCV$popNames, 
                                         indicator = "newHCVdeaths")%>%
        mutate(indicator = "HCVdeath")
  }
  
    
    
    
  }
  else{ 
    if(is.null(pop)){
      scenPrev <- scenario_Prevalence(srlst, sparamlst, pop = NULL, 
                                      statusQ = "y")%>%
        mutate(indicator = "Prevalence")
      
      scenInci <- scenario_Incidence(srlst, sparamlst, pop = NULL, 
                                     statusQ = "y")%>%
        mutate(indicator = "Incidence")
      
      scenCHCVnum <- scenario_CHCVnum(srlst, sparamlst, pop = NULL, 
                                      statusQ = "y")%>%
        mutate(indicator = "CHCVNum")
      
      scenIncinum <- scenario_Incidence(srlst, sparamlst, pop = NULL,
                                        indicator = "newInfections", 
                                        statusQ = "y")%>%
        mutate(indicator = "InciNum")
      
      scenHCVdeath <- scenario_Incidence(srlst,sparamlst, pop = NULL, 
                                         indicator = "newHCVdeaths", 
                                         statusQ = "y")%>%
        mutate(indicator = "HCVdeath")
    }
    else{ 
      scenPrev <- scenario_Prevalence(srlst, sparamlst, pop = HCV$popNames, 
                                      statusQ = "y")%>%
        mutate(indicator = "Prevalence")
      
      scenInci <- scenario_Incidence(srlst, sparamlst, pop = HCV$popNames, 
                                     statusQ = "y")%>%
        mutate(indicator = "Incidence")
      
      scenCHCVnum <- scenario_CHCVnum(srlst, sparamlst, pop = HCV$popNames, 
                                      statusQ = "y")%>%
        mutate(indicator = "CHCVNum")
      
      scenIncinum <- scenario_Incidence(srlst, sparamlst, pop = HCV$popNames,
                                        indicator = "newInfections", 
                                        statusQ = "y")%>%
        mutate(indicator = "InciNum")
      
      scenHCVdeath <- scenario_Incidence(srlst,sparamlst, pop = HCV$popNames, 
                                         indicator = "newHCVdeaths", 
                                         statusQ = "y")%>%
        mutate(indicator = "HCVdeath")
    }
    }
    
    
    rtbl <- rbind(scenPrev, scenInci, scenCHCVnum, scenIncinum, scenHCVdeath)
    
    return(rtbl)
    }
  

#### scenario plots 
ScePlot <- function(dtt, ylabel = NULL,
                    ribbonarea = FALSE,
                    xlimits = NULL, facetPlot = NULL,
                    labelN = NULL){
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
  
  if(is.null(labelN)){ 
    labelN = list(c("Testing scenarios"), c("Status Quo", 
                                            "Point-of-care antibody testing", 
                                            "Dried blood spot testing", 
                                            "Clinic-based reflex RNA testing", 
                                            "Single visit point-of-care RNA testing"))
  }
  
  dtt <- dtt%>%mutate(testing = factor(testing, 
                                       levels = c(unique(dtt$testing))),
                      year = year + HCV$cabY -1)
  
  if(!is.null(ribbonarea)){
    plot <- ggplot(data = dtt, aes_string(x = "year", 
                                          group = "testing")) +
      
      geom_line(aes_string(y = "Med", colour = "testing"), size = 2) +
      
      geom_ribbon(aes_string(ymin = "q25", ymax = "q75", 
                      fill = "testing"), alpha = 0.4) +
      xlab("Year") + ylab(ylabel) + 
      plotOpts + 
      scale_color_OkabeIto(name  =labelN[[1]],
                           breaks=c(unique(dtt$testing)),
                           labels=labelN[[2]]) +
      scale_fill_OkabeIto(name  =labelN[[1]],
                          breaks=c(unique(dtt$testing)),
                          labels=labelN[[2]])
    
    
  }
  else{ 
    plot <- ggplot(data = dtt, aes_string(x = "year", 
                                          group = "testing")) +
      
      geom_line(aes_string(y = "Med", colour = "testing"), size = 2) +
      xlab("Year") + ylab(ylabel) + 
      plotOpts + 
      scale_color_OkabeIto(name  =labelN[[1]],
                           breaks=c(unique(dtt$testing)),
                           labels=c(labelN[[2]]))
  }
  
  
  
  plot <- plot + coord_cartesian(xlim = xlimits[1:2]) +
    scale_x_continuous(expand = c(0,0), limits =c(xlimits[1], xlimits[2]) ,
                       breaks = seq(xlimits[1], xlimits[2], 
                                    by = xlimits[3]),
                       labels = seq(xlimits[1], xlimits[2], 
                                    by = xlimits[3])) + 
    scale_y_continuous(expand = c(0,0))
  
  if (plotFacets) { 
    
    
    # converrt variabel to string for facet_wrap function
    facetPlot <- deparse(substitute(facetPlot))
    
    plot <- plot + facet_wrap(as.formula(paste0("~", facetPlot, collapse = "")),
                              ncol = 2, scales = "free") +
      theme(panel.spacing = unit(2, "lines"))
  }
  return(plot)
} 

#### cascade calculation #### 

# coverage uncertainty 
CosCovcal <- function(dt, pop){
  bestcol <- which(colnames(dt) == "cum_best")
  popSizeQuan <- dt%>%gather(., "simulation", "estimate", 
                             bestcol:ncol(dt))%>%
    filter(simulation!= "best")
  
  
  
  if (is.null(pop)){
    popQ <- popSizeQuan%>%group_by(year) %>%
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
    dt[, c("min", "max", "Med", "Mu", 
           "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
  }
  else{
    
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
    
    dt[, c("min", "max", "Med", "Mu", 
           "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95))  
  }
  
  return(dt)
}

# functions to generate the HCV cascade of care 
HCVCas <- function(HCV, dt, dtpar, pop, scenario, initY){
  # HCV: the rda file save the detail of project 
  # dt: the scenario results file
  # dtpar: the result file of uncertainty of scenario
  # pop: overall==null/subgroups
  # scenario: scenarios results or not, this is related to the uncertainty range 
  # bestResults was 1000 simulation; other scenarios were only 250 simulations.
  # initY: initiation year for coverage 
  if(is.null(scenario)){ 
    samp <- HCV$numberSamples
  }
  else{ 
    samp <- HCV$sParam
  }
  
  if(is.null(pop)){ 
    statepop <- popResults_MidYear(HCV, dt, Population = NULL,
                                   Disease_prog = c("s",HCV$progress_name), 
                                   Cascade = c("s",HCV$cascade_name), 
                                   param = dtpar, endYear = 30, 
                                   YearCut = "End", scenario = scenario)%>%
      as.data.frame()
    # people living with chronic HCV (remove "s", "cured", "acute") 
    statepopa <- statepop%>%mutate(cascade = sub("^[^_]*_", "", state),
                                   disease_prog = sub("\\_.*", "", state))%>%
      filter(!cascade%in%c("s", "cured") & !disease_prog%in%c("a"))%>%
      group_by(year)%>%
      summarise(across(c(best,paste0("set",seq(1, samp,1))), 
                       ~ sum(.x, na.rm = TRUE)))%>%
      mutate(year = year + HCV$cabY - 1)%>%filter(year>=(initY - 1))
    rm(statepop)
    # cumulative diagnosed
    # diagnosed has two pathway: diagnosed by RNA or POCT 
    HCVscreen_all <- indicatorResults(HCV, dt, "newTestingAb", pop="all",
                                    paramR = dtpar, range = 'y',
                                    endY = 30, scenario = scenario)%>%
      select(year, best, c(paste0("set", seq(1, samp,1))))%>%filter(best!=0)%>%
      slice(n())
    
    HCVdiag_all <- indicatorResults(HCV, dt, "newTestingAg", pop="all",
                                    paramR = dtpar, range = 'y',
                                    endY = 30, scenario = scenario)%>%
      select(year, best, c(paste0("set", seq(1, samp,1))))
    
    HCVdiagPOCT_all <- indicatorResults(HCV, dt, "newTestingPOCT", pop="all",
                                        paramR = dtpar,  range = 'y',
                                        endY = 30, scenario = scenario)%>%
      select(year, best, c(paste0("set", seq(1, samp,1))))
    
    HCVdiag <- cbind(year = HCVdiag_all$year + (HCV$cabY - 1),
                     as.data.frame(HCVdiag_all[ , -1] + 
                                     HCVdiagPOCT_all[ ,-1]))%>%
      tibble::as_tibble()%>%as.data.frame()%>%
      filter(year>=(initY - 1))%>%
      mutate(across(c(best, paste0("set", seq(1, samp,1))),
                    cumsum, .names = "cum_{col}"))%>%
      select(year, cum_best, c(paste0("cum_set", seq(1, samp,1))))
    
    rm(HCVdiag_all, HCVdiagPOCT_all)
    # calculate cumulative number of diagnosed since the initiation Year 
    Fir <- HCVdiag[1,2:ncol(HCVdiag)]
    ltp <- HCVscreen_all[1,2:ncol(HCVscreen_all)]
    
    for (i in 1: nrow(HCVdiag)) {
      
      HCVdiag[i, 2:ncol(HCVdiag)] <- HCVdiag[i, 2:ncol(HCVdiag)] - Fir + 
        (ltp*HCV$timestep*5)
    } 
    
    # cumulative treatment initiation 
    HCVtreat_all <- indicatorResults(HCV, dt, "newTreatment", pop="all",
                                     paramR = dtpar, range = 'y',
                                     endY = 30, scenario = scenario)%>%
      mutate(year = year + HCV$cabY - 1)%>%filter(year>=(initY - 1))%>%
      mutate(across(c(best, paste0("set", seq(1, samp,1))),
                    cumsum, .names = "cum_{col}"))%>%
      select(year, c(cum_best, paste0("cum_set", seq(1, samp,1))))
    
    # calculate cumulative number of treatment initiation since the initiation Year
    Fir <- HCVtreat_all[1,2:ncol(HCVtreat_all)]
    
    for (i in 1: nrow(HCVtreat_all)) {
      HCVtreat_all[i, 2:ncol(HCVtreat_all)] <- 
        HCVtreat_all[i, 2:ncol(HCVtreat_all)] - Fir
    } 
    
    # cumulative cured
    HCVCure_all <- indicatorResults(HCV, dt, "newCured", pop="all",
                                    paramR = dtpar, range = 'y',
                                    endY = 30, scenario = scenario)%>%
      mutate(year = year + HCV$cabY - 1)%>%filter(year>=(initY - 1))%>%
      mutate(across(c(best, paste0("set", seq(1, samp,1))),
                    cumsum, .names = "cum_{col}"))%>%
      select(year, c(cum_best, paste0("cum_set", seq(1, samp,1))))
    
    # calculate cumulative number of treatment initiation since the initiation Year
    Fir <- HCVCure_all[1,2:ncol(HCVCure_all)]
    
    for (i in 1: nrow(HCVCure_all)) {
      HCVCure_all[i, 2:ncol(HCVCure_all)] <- 
        HCVCure_all[i, 2:ncol(HCVCure_all)] - Fir
    } 
    
    
    
    # cascade coverage% calculation 
    # take diagnosed coverage as example: 
    #  cumulative diagnosed_y/ (people living with chronic HCV_end of (y-1) + cumulative cured_y)
    
    # people living with HCV (y- 1) 
    statepopa <- slice(statepopa, 1:(n() - 1))
    
    HCVicover <- cbind(year = HCVCure_all[-1, 1],
                       as.data.frame((HCVCure_all[-1, -1] + 
                                        statepopa[, -1]) / 
                                       (HCVCure_all[-1, -1] + 
                                          statepopa[, -1])*100))%>%
      tibble::as_tibble()%>%as.data.frame()%>%
      CosCovcal(., pop = NULL)
    # diagnosed coverage
    HCVdcover <- cbind(year = HCVCure_all[-1, 1],
                       as.data.frame(HCVdiag[-1, -1] / 
                                       (HCVCure_all[-1, -1] + 
                                          statepopa[, -1])*100))%>%
      tibble::as_tibble()%>%as.data.frame()%>%
      CosCovcal(., pop = NULL)
    
    
    HCVtcover <- cbind(year = HCVCure_all[-1, 1],
                       as.data.frame(HCVtreat_all[-1, -1] / 
                                       (HCVCure_all[-1, -1] + 
                                          statepopa[, -1])*100))%>%
      tibble::as_tibble()%>%as.data.frame()%>%
      CosCovcal(., pop = NULL)
    
    HCVccover <- cbind(year = HCVCure_all[-1, 1],
                       as.data.frame(HCVCure_all[-1, -1] / 
                                       (HCVCure_all[-1, -1] + 
                                          statepopa[, -1])*100))%>%
      tibble::as_tibble()%>%as.data.frame()%>%
      CosCovcal(., pop = NULL)
    
    
    Coverage <- list()
    Coverage[["infections"]] <-HCVicover
    Coverage[["diagnosed"]] <-HCVdcover
    Coverage[["treatmentinit"]] <- HCVtcover
    Coverage[["Cured"]] <- HCVccover
    
  }
  else{ statepop <- popResults_MidYear(HCV, dt, Population = pop,
                                       Disease_prog = c("s",HCV$progress_name), 
                                       Cascade = c("s",HCV$cascade_name), 
                                       param = dtpar, endYear = 30, 
                                       YearCut = "End", scenario = scenario)%>%
    as.data.frame()
  # people living with chronic HCV (remove "s", "cured", "acute") 
  statepopa <- statepop%>%mutate(cascade = sub("^[^_]*_", "", state),
                                 disease_prog = sub("\\_.*", "", state))%>%
    filter(!cascade%in%c("s", "cured") & !disease_prog%in%c("a"))%>%
    group_by(year, population)%>%
    summarise(across(c(best,paste0("set",seq(1, samp,1))), 
                     ~ sum(.x, na.rm = TRUE)))%>%
    mutate(year = year + HCV$cabY - 1, 
           population = pop)%>%filter(year>=(initY - 1))%>%
    ungroup()%>%arrange(year, population)
  rm(statepop)
  # cumulative diagnosed
  # diagnosed has two pathway: diagnosed by RNA or POCT 
  HCVscreen_all <- indicatorResults(HCV, dt, "newTestingAb", pop=pop,
                                    paramR = dtpar, range = 'y',
                                    endY = 30, scenario = scenario)%>%
    select(year,population, best, c(paste0("set", seq(1, samp,1))))%>%
    filter(best!=0)%>%arrange(year, population)%>%
    slice(n())
  HCVdiag_all <- indicatorResults(HCV, dt, "newTestingAg", pop=pop,
                                  paramR = dtpar, range = 'y',
                                  endY = 30, scenario = scenario)%>%
    select(year,population, best, c(paste0("set", seq(1, samp,1))))%>%
    arrange(year, population)
  
  HCVdiagPOCT_all <- indicatorResults(HCV, dt, "newTestingPOCT", pop=pop,
                                      paramR = dtpar,  range = 'y',
                                      endY = 30, scenario = scenario)%>%
    select(year,population, best, c(paste0("set", seq(1, samp,1))))%>%
    arrange(year, population)
  
  HCVdiag <- cbind(year = HCVdiag_all$year + (HCV$cabY - 1), 
                   population = HCVdiag_all$population,
                   as.data.frame(HCVdiag_all[ , -c(1,2)] + 
                                   HCVdiagPOCT_all[ ,-c(1,2)]))%>%
    tibble::as_tibble()%>%as.data.frame()%>%
    filter(year>=(initY - 1))%>%
    group_by(population)%>%
    mutate(across(c(best, paste0("set", seq(1, samp,1))),
                  cumsum, .names = "cum_{col}"))%>%
    select(year,population, cum_best, c(paste0("cum_set", seq(1, samp,1))))%>%
    arrange(year, population)
  
  rm(HCVdiag_all, HCVdiagPOCT_all)
  
  # calculate cumulative number of diagnosed since the initiation Year 
  Fir <- HCVdiag[1:HCV$npops,3:ncol(HCVdiag)]
  ltp <- HCVscreen_all[1:HCV$npops, 3:ncol(HCVscreen_all)]
  for (i in 1: nrow(HCVdiag)) {
    if(i%%HCV$npops==1){
      HCVdiag[i, 3:ncol(HCVdiag)] <- HCVdiag[i, 3:ncol(HCVdiag)] - Fir[1,] + 
        (ltp[1,]*HCV$timestep)
    }
    else if(i%%HCV$npops==2){ 
      HCVdiag[i, 3:ncol(HCVdiag)] <- HCVdiag[i, 3:ncol(HCVdiag)] - Fir[2,] + 
        (ltp[2,]*HCV$timestep)
    }
    else if(i%%HCV$npops==3){ 
      HCVdiag[i, 3:ncol(HCVdiag)] <- HCVdiag[i, 3:ncol(HCVdiag)] - Fir[3,] + 
        (ltp[3,]*HCV$timestep)
    }
    else if(i%%HCV$npops==0){ 
      HCVdiag[i, 3:ncol(HCVdiag)] <- HCVdiag[i, 3:ncol(HCVdiag)] - Fir[4,] + 
        (ltp[4,]*HCV$timestep)
      }
    
  } 
  
  # cumulative treatment initiation 
  HCVtreat_all <- indicatorResults(HCV, dt, "newTreatment", pop=pop,
                                   paramR = dtpar, range = 'y',
                                   endY = 30, scenario = scenario)%>%
    mutate(year = year + HCV$cabY - 1)%>%filter(year>=(initY - 1))%>%
    group_by(population)%>%
    mutate(across(c(best, paste0("set", seq(1, samp,1))),
                  cumsum, .names = "cum_{col}"))%>%
    select(year, population, c(cum_best, paste0("cum_set", seq(1, samp,1))))%>%
    arrange(year, population)
  
  # calculate cumulative number of treatment initiation since the initiation Year
  Fir <- HCVtreat_all[1:HCV$npops,3:ncol(HCVtreat_all)]
  
  for (i in 1: nrow(HCVtreat_all)) {
    if(i%%HCV$npops==1){
      HCVtreat_all[i, 3:ncol(HCVtreat_all)] <- 
        HCVtreat_all[i, 3:ncol(HCVtreat_all)] - Fir[1,]
      }
    else if(i%%HCV$npops==2){
      HCVtreat_all[i, 3:ncol(HCVtreat_all)] <- 
        HCVtreat_all[i, 3:ncol(HCVtreat_all)] - Fir[2,]
    }
    else if(i%%HCV$npops==3){ 
      HCVtreat_all[i, 3:ncol(HCVtreat_all)] <- 
        HCVtreat_all[i, 3:ncol(HCVtreat_all)] - Fir[3,]
    }
    else if(i%%HCV$npops==0){ 
      HCVtreat_all[i, 3:ncol(HCVtreat_all)] <- 
        HCVtreat_all[i, 3:ncol(HCVtreat_all)] - Fir[4,]
      }
   
  } 
  
  # cumulative cured
  HCVCure_all <- indicatorResults(HCV, dt, "newCured", pop=pop,
                                  paramR = dtpar, range = 'y',
                                  endY = 30, scenario = scenario)%>%
    mutate(year = year + HCV$cabY - 1)%>%filter(year>=(initY - 1))%>%
    group_by(population)%>%
    mutate(across(c(best, paste0("set", seq(1, samp,1))),
                  cumsum, .names = "cum_{col}"))%>%
    select(year,population, c(cum_best, paste0("cum_set", seq(1, samp,1))))%>%
    arrange(year, population)
  
  # calculate cumulative number of treatment initiation since the initiation Year
  Fir <- HCVCure_all[1:HCV$npops,3:ncol(HCVCure_all)]
  
  for (i in 1: nrow(HCVCure_all)) {
    if(i%%HCV$npops==1){
      HCVCure_all[i, 3:ncol(HCVCure_all)] <- 
        HCVCure_all[i, 3:ncol(HCVCure_all)] - Fir[1, ]
    }
    else if(i%%HCV$npops==2){ 
      HCVCure_all[i, 3:ncol(HCVCure_all)] <- 
        HCVCure_all[i, 3:ncol(HCVCure_all)] - Fir[2, ]
    }
    else if(i%%HCV$npops==3){
      HCVCure_all[i, 3:ncol(HCVCure_all)] <- 
        HCVCure_all[i, 3:ncol(HCVCure_all)] - Fir[3, ]
      
    }
    else if(i%%HCV$npops==0){
      HCVCure_all[i, 3:ncol(HCVCure_all)] <- 
        HCVCure_all[i, 3:ncol(HCVCure_all)] - Fir[4, ]
      
    }
  }
  
  
  # cascade coverage% calculation 
  # take diagnosed coverage as example: 
  #  cumulative diagnosed_y/ (people living with chronic HCV_end of (y-1) + cumulative cured_y)
  
  # people living with HCV (y- 1) 
  statepopa <- slice(statepopa, 1:(n() - HCV$npops))
  
  
  # infection 
  HCVicover <- cbind(year = HCVCure_all[-c(1: HCV$npops), 1],
                     population = HCVCure_all[-c(1: HCV$npops), 2 ], 
                     as.data.frame((HCVCure_all[-c(1: HCV$npops), -c(1,2)] + 
                                      statepopa[, -c(1,2)]) / 
                                     (HCVCure_all[-c(1: HCV$npops), -c(1,2)] + 
                                        statepopa[, -c(1,2)])*100))%>%
    tibble::as_tibble()%>%as.data.frame()%>%
    CosCovcal(., pop = pop)
  
  # diagnosed coverage
  HCVdcover <- cbind(year = HCVCure_all[-c(1: HCV$npops), 1],
                     population = HCVCure_all[-c(1: HCV$npops), 2 ], 
                     as.data.frame(HCVdiag[-c(1: HCV$npops), -c(1,2)] / 
                                     (HCVCure_all[-c(1: HCV$npops), -c(1,2)] + 
                                        statepopa[, -c(1,2)])*100))%>%
    tibble::as_tibble()%>%as.data.frame()%>%
    CosCovcal(., pop = pop)
  
  
  HCVtcover <- cbind(year = HCVCure_all[-c(1: HCV$npops), 1],
                     population = HCVCure_all[-c(1: HCV$npops), 2 ],
                     as.data.frame(HCVtreat_all[-c(1: HCV$npops), -c(1,2)] / 
                                     (HCVCure_all[-c(1: HCV$npops),-c(1,2)] + 
                                        statepopa[, -c(1,2)])*100))%>%
    tibble::as_tibble()%>%as.data.frame()%>%
    CosCovcal(., pop = pop)
  
  HCVccover <- cbind(year = HCVCure_all[-c(1: HCV$npops), 1],
                     population = HCVCure_all[-c(1: HCV$npops), 2 ],
                     as.data.frame(HCVCure_all[-c(1: HCV$npops), -c(1,2)] / 
                                     (HCVCure_all[-c(1: HCV$npops), -c(1,2)] + 
                                        statepopa[, -c(1,2)])*100))%>%
    tibble::as_tibble()%>%as.data.frame()%>%
    CosCovcal(., pop = pop)
  
  
  Coverage <- list()
  Coverage[["infection"]] <-HCVicover
  Coverage[["diagnosed"]] <-HCVdcover
  Coverage[["treatmentinit"]] <- HCVtcover
  Coverage[["Cured"]] <- HCVccover
  }
  
  return(Coverage)
  
}


#### combine the coverage from status quo and scenarios #### 

coverbind <- function(bestcover, scecover, pop){ 
  
  if(is.null(pop)){
    cascade <- c(names(bestcover))
    
    xa <- lapply(seq_along(bestcover), function(x){ 
      a <- bestcover[[x]]%>%
        mutate(scenario = "Status Quo",cascade = names(bestcover)[x])
    })
    
    ttt <- do.call("rbind", xa)%>%
      select(year, scenario, cascade, 
             c(cum_best,min, max, Med,Mu,q5, q25, q75, q95))
    
    xb <- lapply(seq_along(scecover), function(x){
      m <- lapply(seq_along(scecover[[x]]), function(y){
        a <- scecover[[x]][[y]]%>%as.data.frame()%>%
          mutate(scenario = names(scecover)[x],cascade = cascade[y])
      })
      
      aa <- do.call("rbind", m)
    })
    
    uuu <- do.call("rbind", xb)%>%as.data.frame()%>%
      select(year, scenario, cascade, 
             c(cum_best,min, max, Med,Mu, q5, q25, q75, q95))
    
    vvv <- rbind(ttt, uuu)
    
    
    forplot <- vvv%>%gather(quant, value, -c(year, scenario, cascade, cum_best))%>%
      mutate(cascade = factor(cascade, levels = c(unique(vvv$cascade)), 
                              labels = c("Infected","Diagnosed", 
                                         "Treatment initiation", 
                                         "Cured")),
             scenario = factor(scenario, 
                               levels = c(unique(vvv$scenario))))
  }
  else{ 
    cascade <- c(names(bestcover))
    
    xa <- lapply(seq_along(bestcover), function(x){ 
      a <- bestcover[[x]]%>%
        mutate(scenario = "Status Quo",cascade = names(bestcover)[x])
    })
    
    ttt <- do.call("rbind", xa)%>%
      select(year, population, scenario, cascade, 
             c(cum_best,min, max, Med,Mu,q5, q25, q75, q95))
    
    xb <- lapply(seq_along(scecover), function(x){
      m <- lapply(seq_along(scecover[[x]]), function(y){
        a <- scecover[[x]][[y]]%>%as.data.frame()%>%
          mutate(scenario = names(scecover)[x],cascade = cascade[y])
      })
      
      aa <- do.call("rbind", m)
    })
    
    uuu <- do.call("rbind", xb)%>%as.data.frame()%>%
      select(year,population, scenario, cascade, 
             c(cum_best,min, max, Med,Mu, q5, q25, q75, q95))
    
    vvv <- rbind(ttt, uuu)
    
    
    forplot <- vvv%>%gather(quant, value, 
                            -c(year, population, scenario, cascade, cum_best))%>%
      mutate(cascade = factor(cascade, levels = c(unique(vvv$cascade)), 
                              labels = c("Infected", "Diagnosed", 
                                         "Treatment initiated", "Cured")),
             scenario = factor(scenario, 
                               levels = c(unique(vvv$scenario))))
  }
  
  return(forplot)
  
  
}


#### Cascade plot ####






### forest function revised 
## fork from github with a twist on label of x_axis and name for estimate 
# citiation
# citation("forester")
#> 
#> To cite forester in publications use:
#> 
#>   Boyes, Randy (2021). Forester: An R package for creating
#>   publication-ready forest plots. R package version 0.3.0. Available
#>   at: https://github.com/rdboyes/forester
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {Forester: An R package for creating publication-ready forest plots.},
#>     author = {Randy Boyes},
#>     year = {2021},
#>     note = {R package version 0.3.0},
#>     url = {https://github.com/rdboyes/forester},
#>   }



foresterR <- function(left_side_data,
                      estimate,
                      ci_low,
                      ci_high,
                      ci_sep = " to ",
                      right_side_data = NULL,
                      estimate_precision = 1,
                      ggplot_width = 30,
                      null_line_at = 0,
                      file_path = file.path(tempdir(), paste0("forester_plot.", render_as)),
                      dpi = 600,
                      display = TRUE,
                      font_family = "mono",
                      estimate_col_name,
                      stripe_colour = "#eff3f2",
                      background_colour = "white",
                      x_scale_linear = TRUE,
                      xlim = NULL,
                      xbreaks = NULL,
                      nudge_y = 0,
                      nudge_x = 1,
                      nudge_height = 0,
                      nudge_width = 0,
                      justify = 0,
                      arrows = FALSE,
                      arrow_labels = c("Lower", "Higher"),
                      add_plot = NULL,
                      add_plot_width = 1,
                      add_plot_gap = FALSE,
                      point_sizes = 3,
                      point_shapes = 16,
                      center_ggplot = NULL,
                      lower_header_row = FALSE,
                      render_as = "png"){
  
  if(!length(justify) == 1){
    justify <- c(justify[1:(length(justify) - 1)], 0, justify[length(justify)])
    justify <- matrix(justify, ncol=length(justify), nrow=nrow(left_side_data) + 3, byrow=TRUE)
    justify <- as.vector(justify)
  }
  
  if(lower_header_row == FALSE){
    theme <- gridExtra::ttheme_minimal(core=list(
      fg_params = list(hjust = justify, x = (0.05 + (0.45/0.5) * justify), fontfamily = font_family),
      bg_params = list(fill=c(rep(c(stripe_colour, background_colour), length.out=nrow(left_side_data)), background_colour, background_colour, background_colour))
    ),
    colhead = list(fg_params = list(hjust = 0, x = 0.05,
                                    fontfamily = font_family),
                   bg_params = list(fill = background_colour))
    )
  }else{
    theme <- gridExtra::ttheme_minimal(core=list(
      fg_params = list(hjust = justify, x = (0.05 + (0.45/0.5) * justify), fontfamily = font_family),
      bg_params = list(fill=c(rep(c(background_colour, stripe_colour), length.out=nrow(left_side_data)), background_colour, background_colour, background_colour))
    ),
    colhead = list(fg_params = list(hjust = 0, x = 0.05,
                                    fontfamily = font_family,
                                    fill = background_colour),
                   bg_params = list(fill = background_colour))
    )
  }
  
  gdata <- data.frame(estimate = estimate,
                      ci_low = ci_low,
                      ci_high = ci_high)
  
  if(lower_header_row){
    gdata <- tibble::add_row(gdata, .before = 1)
  }
  
  if(is.null(right_side_data)){
    tdata <- gdata
    
    tdata <- dplyr::mutate_all(tdata, ~sprintf(.,
                                               fmt = paste0('%#.', estimate_precision,'f')
    ))
    
    tdata[tdata == "NA"] <- " "
    # pretty formatting for confidence intervals
    right_side_data <- data.frame(Estimate = ifelse(tdata$estimate == " ",
                                                    " ", paste0(tdata$estimate, " (", tdata$ci_low,
                                                                ci_sep, tdata$ci_high, ")")))
    
    colnames(right_side_data) <- estimate_col_name
    
  }
  
  # finds width in number of characters for monospaced font
  
  find_width_mono <- function(data){
    num_of_rows <- nrow(data)
    num_of_cols <- ncol(data)
    
    print_data <- dplyr::mutate_all(data, as.character)
    
    num_char_across <- 0
    width <- 0
    
    for(i in 1:num_of_cols){
      for(j in 1:num_of_rows){
        num_char_across[j] <- nchar(print_data[j, i])
      }
      width[i] <- max(max(num_char_across, na.rm = TRUE),
                      nchar(colnames(print_data)[i]), na.rm = TRUE)
    }
    return(sum(width, na.rm = TRUE))
  }
  
  # finds width using shape_string from the systemfonts package
  # if not using monospaced font
  
  find_width <- function(data){
    num_of_rows <- nrow(data)
    num_of_cols <- ncol(data)
    
    print_data <- dplyr::mutate_all(data, as.character)
    
    width <- 0
    
    names <- colnames(print_data)
    
    for (i in 1:num_of_cols){
      temp <- systemfonts::shape_string(print_data[[names[i]]], family = font_family)
      temp_col <- systemfonts::shape_string(names[i], family = font_family)
      width[i] <- max(max(temp$metrics$width, na.rm = TRUE),
                      temp_col$metrics$width, na.rm = TRUE)
    }
    return(sum(width, na.rm = TRUE)/7.2)
  }
  
  # calculate widths for each side with the appropriate function
  
  if(font_family == "mono"){
    left_width <- find_width_mono(left_side_data) 
    right_width <- find_width_mono(right_side_data) + 2
  }else{
    left_width <- find_width(left_side_data) 
    right_width <- find_width(right_side_data) + 2
  }
  
  # insert a blank column so we can put the ggplot object on top
  # and correctly order columns
  
  total_width <- left_width + right_width + 2 + ggplot_width
  
  tdata_print <- left_side_data 
  
  if(lower_header_row){rbind.data.frame(colnames(tdata_print), tdata_print)}
  
  tdata_print$` ` <- paste(rep(" ", times = round(ggplot_width, 0)),
                           collapse = '')
  tdata_print <- cbind(tdata_print, right_side_data)
  
  tdata_print <- tibble::add_row(tdata_print)
  tdata_print <- tibble::add_row(tdata_print)
  tdata_print <- tibble::add_row(tdata_print)
  
  tdata_print <- dplyr::mutate_all(tdata_print, as.character)
  tdata_print[is.na(tdata_print)] <- " "
  
  ## formatting functions
  
  mono_column <- function(table, col){
    col_indexes <- function(table, col, name="core-fg"){
      l <- table$layout
      which(l$l==col & l$name==name)
    }
    
    ind <- col_indexes(table, col, "core-fg")
    
    for(i in ind){
      table$grobs[i][[1]][["gp"]] <- grid::gpar(fontfamily = "mono")
    }
    return(table)
  }
  
  white_column <- function(table, col){
    col_indexes <- function(table, col, name="core-bg"){
      l <- table$layout
      which(l$l==col & l$name==name)
    }
    
    ind <- col_indexes(table, col, "core-bg")
    ind_fg <- col_indexes(table, col, "core-fg")
    
    for(i in ind){
      table$grobs[i][[1]][["gp"]] <- grid::gpar(fill = background_colour, col = background_colour)
    }
    
    for(i in ind_fg){
      table$grobs[i][[1]][["gp"]] <- grid::gpar(fontfamily = "mono")
    }
    return(table)
  }
  
  ######## calculations for the top and bottom of the plot
  
  gdata$row_num <- (nrow(gdata) - 1):0
  
  
  h_adj <- dplyr::case_when(
    font_family == "mono" ~ 0.2,
    font_family == "serif" ~ .43,
    font_family == "sans" ~ .37,
    TRUE ~ 0
  )
  
  h_adj <- nudge_y + h_adj
  
  
  slope_adj <- dplyr::case_when(
    font_family == "mono" ~ -0.175,
    font_family == "serif" ~ -.19,
    font_family == "sans" ~ -.16,
    TRUE ~ 0
  )
  
  font_adj <- 0.3 + h_adj + log(nrow(gdata)) * slope_adj
  
  y_low <- -.5 + font_adj + -.1381 * log(nrow(gdata))
  y_high <- 1.017 * nrow(gdata) - 0.6
  
  #### add shapes and sizes to gdata ########
  
  gdata$shape <- point_shapes
  gdata$sizes <- point_sizes
  
  #### if a ci will be out of bounds, add arrow on the oob side  ###############
  
  g_oob <- tibble::tibble()
  
  if(!is.null(xlim)){
    oob_arrows <- gdata
    
    oob_arrows$x_low <- xlim[1]
    oob_arrows$x_high <- xlim[2]
    
    ra <- sum(oob_arrows$ci_high > oob_arrows$x_high, na.rm = T) > 0
    la <- sum(oob_arrows$ci_low < oob_arrows$x_low, na.rm = T) > 0
    
    if(ra){
      right_arrows <- dplyr::select(dplyr::filter(oob_arrows, ci_high > .data$x_high), start = estimate, end = .data$x_high, y = .data$row_num)
    }
    if(la){
      left_arrows <- dplyr::select(dplyr::filter(oob_arrows, ci_low < .data$x_low), start = estimate, end = .data$x_low, y = .data$row_num)
    }
    
    if(ra && !la){
      g_oob <- right_arrows
    }else if(!ra && la){
      g_oob <- left_arrows
    }else if(ra && la){
      g_oob <- rbind.data.frame(right_arrows, left_arrows)
    }
  }
  
  ########## the main figure - this will be overlaid on the table ##############
  
  center <- ggplot2::ggplot() +
    ggplot2::geom_point(data = gdata, ggplot2::aes(y = row_num, x = estimate, size = sizes, shape = shape), na.rm = TRUE) +
    ggplot2::geom_errorbarh(data = gdata, ggplot2::aes(y = row_num,
                                                       xmin = ci_low,
                                                       xmax = ci_high),
                            height = .25,
                            na.rm = TRUE) +
    ggplot2::theme_classic() + # base theme
    ggplot2::theme(axis.title.y = ggplot2::element_blank(), # remove axis, make bg transparent
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank(),
                   axis.ticks.length.x = grid::unit(.07, "in"),
                   text = ggplot2::element_text(family = font_family, size = 12),
                   panel.background = ggplot2::element_rect(fill = "transparent"),
                   plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(fill = "transparent"),
                   legend.box.background = ggplot2::element_rect(fill = "transparent")) +
    ggplot2::geom_vline(xintercept = null_line_at, linetype = "dashed") +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_shape_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::xlab("")
  
  ### add oob arrows if required ###
  
  if(nrow(g_oob) > 0){
    center <- center +
      ggplot2::geom_segment(data = g_oob,
                            ggplot2::aes(x = start,
                                         xend = end,
                                         y = y,
                                         yend = y),
                            arrow = ggplot2::arrow(angle = 15,
                                                   type = "closed",
                                                   length = grid::unit(0.1, "in")))
  }
  
  ####### fix plot zoom ######
  
  if(is.null(xlim)){
    center <- center + ggplot2::coord_cartesian(ylim = c(y_low, y_high))
  }else{
    center <- center + ggplot2::coord_cartesian(ylim = c(y_low, y_high), xlim = xlim)
  }
  
  ######## handle breaks, log vs linear scales ########
  
  if(x_scale_linear){
    if(is.null(xbreaks)){
      center <- center + ggplot2::scale_x_continuous(labels = scales::number_format(accuracy = 0.1),
                                                     expand = c(0,0))
    }else{
      center <- center + ggplot2::scale_x_continuous(labels = c(paste0(xbreaks, "%", sep = "")),
                                                     breaks = xbreaks,
                                                     expand = c(0,0))
    }
  }else{
    if(is.null(xbreaks)){
      center <- center + ggplot2::scale_x_log10(labels = scales::number_format(accuracy = 0.1),
                                                expand = c(0,0))
    }else{
      center <- center + ggplot2::scale_x_log10(labels = c(paste0(xbreaks, "%", sep = "")),
                                                breaks = xbreaks,
                                                expand = c(0,0))
    }
  }
  
  #### allow overwrite of central plot #######################
  
  if(!is.null(center_ggplot)){center <- center_ggplot}
  
  ######################## Arrows ##############################
  if(arrows == TRUE){
    
    # this df has the text labels
    xlab_df <- data.frame(text = arrow_labels,
                          x = xlim,
                          y = c(0, 0),
                          hjust = c(0, 1))
    
    a_small_amount <- abs(xlim[1] - xlim[2])/35
    
    # this df has the arrows
    if(x_scale_linear == TRUE){
      arrow_df <- data.frame(id = c(1,2),
                             xstart = c(null_line_at - a_small_amount, null_line_at + a_small_amount),
                             xend = c(xlim[1] + a_small_amount, xlim[2] - a_small_amount),
                             y = c(1, 1))
    }else{
      arrow_df <- data.frame(id = c(1,2),
                             xstart = c(null_line_at - a_small_amount, null_line_at + a_small_amount),
                             xend = c(xlim[1], xlim[2]),
                             y = c(1, 1))
    }
    
    # create the arrow/label ggplot object
    arrows_plot <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = arrow_df, ggplot2::aes(x = .data$xstart, xend = .data$xend, y = .data$y, yend = .data$y),
                            arrow = ggplot2::arrow(angle = 15, type = "closed", length = grid::unit(0.1, "in"))) +
      ggplot2::geom_text(data = xlab_df, ggplot2::aes(x = .data$x, y = .data$y, label = .data$text, hjust = .data$hjust),
                         family = font_family, size = 3) +
      ggplot2::scale_y_continuous(expand = c(0,0), limits = c(-0.5, 1.75)) +
      ggplot2::scale_x_continuous(expand = c(0,0), limits = xlim) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "transparent"),
                     plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(fill = "transparent"),
                     legend.box.background = ggplot2::element_rect(fill = "transparent"),
                     panel.border = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.line.y = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.line.x = ggplot2::element_blank())
    
    if(x_scale_linear == FALSE){
      arrows_plot <- arrows_plot + ggplot2::scale_x_log10(expand = c(0,0), limits = xlim)
    }
    
  }
  
  
  ######### using patchwork, overlay the ggplot on the table ###################
  
  png_width <- total_width/10 + nudge_x
  png_height <- (nrow(gdata) + 3)/3.8
  
  if(is.null(add_plot)){
    
    table_final <- mono_column(gridExtra::tableGrob(tdata_print, theme = theme, rows = NULL), ncol(left_side_data) + 1)
    
    table_final$widths[ncol(left_side_data) + 1] <- grid::unit(ggplot_width/10, "in")
    
    table_final$heights <- grid::unit(rep(0.255, times = length(table_final$heights)), "in")
    
    final <- patchwork::wrap_elements(table_final) +
      patchwork::inset_element(center,
                               align_to = "full",
                               left = (left_width/total_width),
                               right = ((ggplot_width + left_width)/total_width),
                               top = 1,
                               bottom = 0.35/nrow(gdata))
    
    if(arrows == TRUE){
      final <- final + patchwork::inset_element(arrows_plot,
                                                align_to = "full",
                                                left = (left_width/total_width),
                                                right = ((ggplot_width + left_width)/total_width),
                                                top = 1.5/nrow(gdata),
                                                bottom = 0)
    }
    
  }else{
    tdata_print$`  ` <- paste(rep(" ", times = round(ggplot_width, 0)),
                              collapse = '')
    
    table_final <- mono_column(gridExtra::tableGrob(tdata_print, theme = theme, rows = NULL), ncol(left_side_data) + 1)
    
    table_final <- white_column(table_final, ncol(table_final))
    
    table_final$widths[ncol(left_side_data) + 1] <- grid::unit(ggplot_width/10, "in")
    table_final$widths[ncol(table_final)] <- grid::unit(ggplot_width/10, "in")
    
    table_final$heights <- grid::unit(rep(0.255, times = length(table_final$heights)), "in")
    
    new_full_width <- total_width + ggplot_width
    
    png_width <- new_full_width/10 + nudge_x
    
    if (add_plot_gap){
      add_plot <- add_plot + ggplot2::scale_y_continuous(limits = c(y_low, y_high), expand = c(0,0)) +
        ggplot2::theme_classic() + # base theme
        ggplot2::theme(axis.title.y = ggplot2::element_text(colour = "transparent"), # make axis transparent rather than removing it
                       axis.text.y = ggplot2::element_text(colour = "transparent"), # this makes alignment much easier
                       axis.ticks.y = ggplot2::element_line(colour = "transparent"),
                       axis.line.y = ggplot2::element_line(colour = "transparent"),
                       axis.title.x = ggplot2::element_text(colour = "transparent"),
                       axis.text.x = ggplot2::element_text(colour = "transparent"),
                       axis.ticks.x = ggplot2::element_line(colour = "transparent"),
                       axis.line.x = ggplot2::element_line(colour = "transparent"),
                       panel.background = ggplot2::element_rect(fill = "transparent"),
                       plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.background = ggplot2::element_rect(fill = "transparent"),
                       legend.box.background = ggplot2::element_rect(fill = "transparent"),
                       legend.position = "none")
    }else{
      add_plot <- add_plot + ggplot2::scale_y_continuous(limits = c(y_low, y_high), expand = c(0,0)) +
        ggplot2::theme_classic() + # base theme
        ggplot2::theme(axis.title.x = ggplot2::element_text(colour = "transparent"), # make x axis (only) transparent
                       axis.text.x = ggplot2::element_text(colour = "transparent"),
                       axis.ticks.x = ggplot2::element_line(colour = "transparent"),
                       axis.line.x = ggplot2::element_line(colour = "transparent"),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "transparent"),
                       plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.background = ggplot2::element_rect(fill = "transparent"),
                       legend.box.background = ggplot2::element_rect(fill = "transparent"),
                       legend.position = "none")
    }
    
    final <- patchwork::wrap_elements(table_final) +
      patchwork::inset_element(center,
                               align_to = "full",
                               left = (left_width/new_full_width),
                               right = ((ggplot_width + left_width)/new_full_width + 1),
                               top = 1,
                               bottom = 0.35/nrow(gdata)) +
      patchwork::inset_element(add_plot,
                               align_to = "full",
                               left = total_width/new_full_width,
                               right = 1,
                               top = 1,
                               bottom = 0.35/nrow(gdata))
    
    if(arrows == TRUE){
      final <- final + patchwork::inset_element(arrows_plot,
                                                align_to = "full",
                                                left = (left_width/new_full_width),
                                                right = ((ggplot_width + left_width)/new_full_width + 1),
                                                top = 1.5/nrow(gdata),
                                                bottom = 0)
    }
  }
  
  ######### save the plot as a png, then display it with magick ################
  if(!(render_as == "rmarkdown")){
    ggplot2::ggsave(
      dpi = dpi,
      height = png_height + nudge_height,
      width = png_width + nudge_width,
      units = "in",
      filename = file_path,
      device = render_as
    )
    
    if(display == TRUE){
      system(paste0('open "', file_path, '"'))
    }
  }else{
    final
  }
}












#### sensitivity: Tornado plot #### 
# Step1: create a dataframe to save the results 
# Testing   Scenarios         Parameters    Incidence   Prevalence
#---------|----------------|-------------|-----------|------------|
# Base_POC    X                   XX%         XX%         XX%
# POC ab    diagnosed rate      30%           YY%          YY% 
# POC ab    diagnosed rate      50%           YY%          YY% 
# ......
#### import base scenario dataframe 
### import dataset 






Sensi_dlst <- function(Res, dtbase, scen){ 
  
  
  
  ScentempTotal_all <- list()
  ScentempTotal_all <- lapply(Res, function(x){ 
    popResults_MidYear(HCV, x,Population = NULL,
                       Disease_prog = NULL, 
                       Cascade = NULL, param = NULL, 
                       endYear = 50)})
  
  
  ScentempNOTInfected_all <- list()
  
  ScentempNOTInfected_all <- lapply(Res, function(x){ 
    popResults_MidYear(HCV, x,Population = NULL,
                       Disease_prog = NULL, 
                       Cascade = c("s", "cured"), 
                       param = NULL ,endYear = 50)%>%
      as.data.frame()%>%group_by(year)%>%dplyr::summarize(best = sum(best))})
  
  Sce_all <- list()
  
  for (i in 1: length(ScentempNOTInfected_all)){ 
    Sce_all[[i]] <- cbind(year = seq(HCV$startYear , 50-1 ,1),
                          as.data.frame((ScentempTotal_all[[i]][, -1] - 
                                           ScentempNOTInfected_all[[i]][ ,-1])/ 
                                          ScentempTotal_all[[i]][ ,-1]*100))%>%
      tibble::as_tibble() 
  }
  
  aa <- list()
  
  aa <- list()
  aa <- lapply(Sce_all, function(x) cbind(x$best))%>%as.data.frame()
  aa <- cbind(Sce_all[[1]]$year,rep("Overall", length(Sce_all[[1]]$year)),aa)%>%
    as.data.frame()
  
  names(aa) <- c("year","population",c(names(Res)))
  
  aalng <- aa%>%pivot_longer(!c(year, population), 
                           names_to = "testing",
                           values_to = "prevalence") 
  
  
  ScentempTotal_group <- list()
  ScentempTotal_group <- lapply(Res, function(x){ 
    popResults_MidYear(HCV, x,Population = HCV$popNames,
                       Disease_prog = NULL, 
                       Cascade = NULL, param = NULL, 
                       endYear = 50)})
  
  ScentempNOTInfected_group <- list()
  
  ScentempNOTInfected_group <- lapply(Res, function(x){ 
    popResults_MidYear(HCV, x,Population = HCV$popNames,
                       Disease_prog = NULL, 
                       Cascade = c("s", "cured"), 
                       param = NULL ,endYear = 50)%>%
      as.data.frame()%>%group_by(year, population)%>%
      dplyr::summarize(best = sum(best))})
  
  Sce_group <- list()
  
  for (i in 1: length(ScentempNOTInfected_group)){ 
    Sce_group[[i]] <- cbind(year = rep(seq(HCV$startYear , 50-1 ,1), each = HCV$npops),
                            population = HCV$popNames,
                            as.data.frame((ScentempTotal_group[[i]][, -c(1,2)] - 
                                             ScentempNOTInfected_group[[i]][ ,-c(1,2)])/ 
                                            ScentempTotal_group[[i]][ ,-c(1,2)]*100))%>%
      tibble::as_tibble() 
  }
  
  
  bb <- list()
  bb <- lapply(Sce_group, function(x) cbind(x$best))%>%as.data.frame()
  bb <- cbind(Sce_group[[1]]$year, Sce_group[[1]]$population,bb)%>%as.data.frame()
  names(bb) <- c("year","population",c(names(Res)))
  
  bblng <- bb%>%pivot_longer(!c(year, population), 
                             names_to = "testing",
                             values_to = "prevalence") 
  bblng<- bblng%>%mutate(year = as.numeric(year),
                         prevalence = as.numeric(prevalence))
  bblng <- rbind(aalng,bblng)%>%as.data.frame()%>%arrange(year, population)%>%
    mutate(scenario =scen)%>%
    select(year, population, testing, scenario, prevalence)

  #### incidence #### 
  Sce_newInf <- lapply(Res, 
                       function(x) indicatorResult_uno(HCV, x, "newInfections", 
                                                       endYear = 50))%>%
    lapply(., function(x){
      
      a <- x%>%group_by(year)%>%dplyr::summarise(best = sum(best))
      
      
      
    }) 
  names(Sce_newInf) <- names(Res)
  
  SceInc_all <- list()
  
  for (i in 1: length(Sce_newInf)){ 
    
    SceInc_all[[i]] <- cbind(year = seq(HCV$startYear , 50-1 ,1),
                             as.data.frame(Sce_newInf[[i]][, -1]/ 
                                             ScentempTotal_all[[i]][ ,-1]*1000))%>%
      tibble::as_tibble() 
    
  }
  
  
  Sce_newInf_g <- lapply(Res, 
                         function(x) indicatorResult_uno(HCV, x, "newInfections", 
                                                         populations = HCV$popNames,
                                                         endYear = 50))%>%
    lapply(., function(x){
      
      a <- x%>%group_by(year, population)%>%dplyr::summarise(best = sum(best))
      
      
      
    }) 
  
  names(Sce_newInf_g) <- names(Res) 
  
  ScetempInc_g <- list()
  
  for (i in 1: length(Sce_newInf_g)){ 
    
    ScetempInc_g[[i]] <- cbind(year = rep(seq(HCV$startYear , 50-1 ,1), each = HCV$npops),
                               population = HCV$popNames,
                               as.data.frame(Sce_newInf_g[[i]][, -c(1,2)]/ 
                                               ScentempTotal_group[[i]][ ,-c(1,2)]*1000))%>%
      tibble::as_tibble() 
    
  }
  
  
  
  
  
  
  aaInc <- list()
  aaInc <- lapply(SceInc_all, function(x) cbind(x$best))%>%as.data.frame()
  aaInc <- cbind(SceInc_all[[1]]$year,
                 rep("Overall", length(SceInc_all[[1]]$year)),aaInc)%>%
    as.data.frame()
  
  names(aaInc) <- c("year","population",c(names(Res)))
  aaInclng <- aaInc%>%pivot_longer(!c(year, population), 
                             names_to = "testing",
                             values_to = "incidence") 
  
  bbInc <- list()
  bbInc <- lapply(ScetempInc_g, function(x) cbind(x$best))%>%as.data.frame()
  bbInc <- cbind(ScetempInc_g[[1]]$year, ScetempInc_g[[1]]$population,bbInc)%>%as.data.frame()
  names(bbInc) <- c("year","population",c(names(Res)))
  bbInclng <- bbInc%>%pivot_longer(!c(year, population), 
                             names_to = "testing",
                             values_to = "incidence") 
  bbInclng<- bbInclng%>%mutate(year = as.numeric(year),
                         incidence = as.numeric(incidence))
  bbInclng <- rbind(aaInclng,bbInclng)%>%as.data.frame()%>%
    arrange(year, population)
  
  
  # number of HCV new infection
  aaIncNum <- list()
  aaIncNum <- lapply(Sce_newInf, function(x) cbind(x$best))%>%as.data.frame()
  aaIncNum <- cbind(Sce_newInf[[1]]$year,
                 rep("Overall", length(Sce_newInf[[1]]$year)),aaIncNum)%>%
    as.data.frame()
  
  names(aaIncNum) <- c("year","population",c(names(Res)))
  aaIncNumlng <- aaIncNum%>%pivot_longer(!c(year, population), 
                                   names_to = "testing",
                                   values_to = "numberInf") 
  
  bbIncNum <- list()
  bbIncNum <- lapply(Sce_newInf_g, function(x) cbind(x$best))%>%as.data.frame()
  bbIncNum <- cbind(Sce_newInf_g[[1]]$year, Sce_newInf_g[[1]]$population,bbIncNum)%>%as.data.frame()
  names(bbIncNum) <- c("year","population",c(names(Res)))
  bbIncNumlng <- bbIncNum%>%pivot_longer(!c(year, population), 
                                   names_to = "testing",
                                   values_to = "numberInf") 
  bbIncNumlng <- bbIncNumlng%>%mutate(year = as.numeric(year),
                               numberInf = as.numeric(numberInf))
  bbIncNumlng <- rbind(aaIncNumlng,bbIncNumlng)%>%as.data.frame()%>%
    arrange(year, population)
  
  
  
  
  
  
  sensi_dt <- bblng%>%as.data.frame()%>%
    mutate(incidence = bbInclng$incidence,
           numberInfect = bbIncNumlng$numberInf)
  sensi_dt <- rbind(dtBase, sensi_dt)%>%as.data.frame()
  return(sensi_dt)
  }


#### Scenario 95% CR #### 
# indicator: c("incidence", "newInfections", "newHCVdeaths", "HCVMortality", 
#              "newCured", "RSVR")
#### ####
sResCR <- function(sRes,sResPar, numberSamples, indicator, pop){
  
  # popsize all#
  popSPar <- popResults_MidYear(HCV, sRes[[1]], Population = NULL,
                                Disease_prog = NULL, 
                                Cascade = NULL, 
                                param = sResPar, 
                                endYear = 50)%>%as.data.frame()
  # new infections all #
  HCVInfect_SPar <- indicatorResults(HCV, sRes[[1]], 
                                     "newInfections", pop="all",
                                     paramR = sResPar, range = "y",
                                     endY = 50)
  
  # reset colnames #
  colsetPar <- HCVInfect_SPar%>%select(., contains("Set"))%>%
    select(., str_sort(paste0("set", seq(1, numberSamples,1)),
                       numeric = FALSE))
  
  HCVInfectPar <- HCVInfect_SPar%>%select(year, best, colnames(colsetPar))
  
  # popsize pop # 
  popSPar_pop <- popResults_MidYear(HCV, sRes[[1]], Population = HCV$popNames,
                                Disease_prog = NULL, 
                                Cascade = NULL, 
                                param = sResPar, 
                                endYear = 50)%>%as.data.frame()
  # new infections pop #
  HCVInfect_SPar_pop <- indicatorResults(HCV, sRes[[1]], "newInfections", 
                                     pop=HCV$popNames, paramR = sResPar, 
                                     range = "y",endY = 50)%>%
    select(year, population, best, colnames(colsetPar))
  
  
  # newHCVdeath #
  
  HCVD_SPar <- indicatorResults(HCV, sRes[[1]], 
                                     "newHCVdeaths", pop="all",
                                     paramR = sResPar, range = "y",
                                     endY = 50)
  
  # reset colnames #
  colsetPar <- HCVD_SPar%>%select(., contains("Set"))%>%
    select(., str_sort(paste0("set", seq(1, numberSamples,1)),
                       numeric = FALSE))
  
  HCVDPar <- HCVD_SPar%>%select(year, best, colnames(colsetPar))
  
  
  # new HCVdeath pop #
  HCVD_SPar_pop <- indicatorResults(HCV, sRes[[1]], "newHCVdeaths", 
                                         pop=HCV$popNames, paramR = sResPar, 
                                         range = "y",endY = 50)%>%
    select(year, population, best, colnames(colsetPar))
  
  
  # S and cured population all #
  
  StempNOTInfected_all <- popResults_MidYear(HCV, sRes[[1]],Population = NULL,
                                            Disease_prog = NULL, 
                                            Cascade = c("s", "cured"), 
                                            param = sResPar ,endYear = 50)%>%
    as.data.frame()%>%group_by(year)%>%
    summarise_at(.vars = c(colnames(popSPar)[-1]),sum)%>%
    select(year, best, colnames(colsetPar))
  
  
  # S and cured population all pop #
  StempNOTInfected <- popResults_MidYear(HCV, sRes[[1]],
                                        Population = HCV$popNames,
                                        Disease_prog = NULL, 
                                        Cascade = c("s", "cured"), 
                                        param = sResPar,
                                        endYear = 50)%>%
    group_by(year, population)%>%
    summarise_at(.vars = c(colnames(popSPar)[-1]), sum)%>%
    select(year, population, best, colnames(colsetPar))
  
  
  # cure
  SCureInit <- indicatorResults(HCV, sRes[[1]], "newCured", 
                                pop="all",
                                paramR = sResPar, range = "y",
                                endY = 50) 
  
  SCureInitPop <-indicatorResults(HCV, sRes[[1]], "newCured", 
                                  pop= HCV$popNames,
                                  paramR = sResPar, range = "y",
                                  endY = 50) 
  
  # incidence 
  if(isTRUE(indicator =="incidence")){
    if(is.null(pop)==TRUE){ 
      
      popSPar <- popSPar%>%select(year, best, colnames(colsetPar))
      HCVIncPar <- cbind(year = seq(HCV$startYear , 50-1 ,1),
                         as.data.frame(HCVInfectPar[, -1]/popSPar[ ,-1]*1000))%>%
        tibble::as_tibble() 
      # get aggregate columns 
      
      HCVInc_rangePar <- popResults_range(HCV, HCVIncPar, Population = NULL,
                                          Disease_prog = NULL , 
                                          Cascade = NULL, end_Y = 50)
    
    }else{ 
      
      popSPar_pop <- popSPar_pop%>%select(year, population, best, 
                                          colnames(colsetPar))
      
      HCVIncPar <- cbind(year = rep(seq(HCV$startYear , 50-1 ,1), 
                                    each = HCV$npops),
                         population = popSPar_pop$population,
                         as.data.frame(HCVInfect_SPar_pop[, -c(1,2)] / 
                                         popSPar_pop[ ,-c(1,2)]*1000))%>%
        tibble::as_tibble() 
      
      
      # get aggregate columns 
      
      HCVInc_rangePar <- popResults_range(HCV, HCVIncPar, 
                                          Population = HCV$popNames,
                                          Disease_prog = NULL , 
                                          Cascade = NULL, end_Y = 50)
    }
    return(HCVInc_rangePar)
    }
  else if(isTRUE(indicator =="newInfections")){
    if(is.null(pop)==TRUE){ 
      NewInf <- indicatorResults(HCV, sRes[[1]], "newInfections", 
                                 pop="all", paramR = sResPar, 
                                 range = "y",endY = 50)
    }else{
      NewInf <- indicatorResults(HCV, sRes[[1]], "newInfections", 
                                 pop=HCV$popNames, paramR = sResPar, 
                                 range = "y",endY = 50) 
      }
    return(NewInf)
  }
  else if(isTRUE(indicator =="newHCVdeaths")){ 
    if(is.null(pop)==TRUE){ NewHCVD <- HCVD_SPar}
    else{ NewHCVD <- indicatorResults(HCV, sRes[[1]], "newHCVdeaths", 
                                      pop=HCV$popNames, paramR = sResPar, 
                                      range = "y",endY = 50)
    }
    return(NewHCVD)
  }
  else if(isTRUE(indicator =="HCVMortality")){
    if(is.null(pop)==TRUE){
      
      popSPar <- popSPar%>%select(year, best, colnames(colsetPar))
      HCVMor <- cbind(year = seq(HCV$startYear , 50-1 ,1),
                      as.data.frame(100*HCVDPar[, -1]/popSPar[, -1] ))%>%
        tibble::as_tibble()  
      HCVMor_rangePar <- popResults_range(HCV, HCVMor, Population = NULL,
                                          Disease_prog = NULL , 
                                          Cascade = NULL, end_Y = 50)
      
    } else{ 
      popSPar_pop <- popSPar_pop%>%select(year, population, best, 
                                      colnames(colsetPar))
      HCVMor <- cbind(year = rep(seq(HCV$startYear , 50-1 ,1), each = HCV$npops),
                      population = HCV$popNames,
                      as.data.frame(100*HCVD_SPar_pop[, -c(1,2)]/
                                      popSPar_pop[, -c(1,2)]))%>%
        tibble::as_tibble() 
      
      
      HCVMor_rangePar <- popResults_range(HCV, HCVMor, Population = HCV$popNames,
                                          Disease_prog = NULL , 
                                          Cascade = NULL, end_Y = 50)
      }
    return(HCVMor_rangePar )
    
  }
  else if(isTRUE(indicator =="newCured")){ 
    if(is.null(pop)==TRUE){
      Scure <- SCureInit
      }
    else{ Scure <- SCureInitPop }
   
    return(Scure) 
  }
  else if(isTRUE(indicator =="RSVR")){ 
    if(is.null(pop)==TRUE){
      scure <- SCureInit%>%select(year, best, colnames(colsetPar))
      popSPar <- popSPar%>%select(year, best, colnames(colsetPar))
      
      
      Svr <-cbind(year = seq(HCV$startYear , 50-1 ,1),
                    as.data.frame(100*scure[, -1]/
                                    (popSPar[, -1] - StempNOTInfected_all[ ,-1])))%>%
        tibble::as_tibble()
      
      Svr_rangePar <- popResults_range(HCV, Svr, Population = NULL,
                                       Disease_prog = NULL , 
                                       Cascade = NULL, end_Y = 50)
      
      }
    else{ 
      scurepop <- SCureInitPop%>%select(year, population, best, 
                                        colnames(colsetPar))
      popSPar_pop <- popSPar_pop%>%select(year, population, best, 
                                          colnames(colsetPar))
      Svr <- cbind(year = rep(seq(HCV$startYear , 50-1 ,1), each = HCV$npops),
                      population = HCV$popNames,
                      as.data.frame(100*scurepop[-c(1,2)]/
                                      (popSPar_pop[, -c(1,2)] - 
                                         StempNOTInfected[ ,-c(1,2)])))%>%
      tibble::as_tibble() 
      
      Svr_rangePar <- popResults_range(HCV, Svr, Population = HCV$popNames,
                                     Disease_prog = NULL , 
                                     Cascade = NULL, end_Y = 50)
    
    }
    
    return(Svr_rangePar) 
  }
}


#### combine col of scenarioResults step1: ### 
SceTable <- function(dt){ 
  
  yt <- lapply(dt, function(x){ 
    
    tt <- x%>%select(year, best,Mu, Med, q5, q95)%>%mutate(year = year + 2003)
    
    tt <- tt%>%filter(year%in%c(2015,seq(2022,2030)))
    
    
    })
  bas <- list()
  bas <- lapply(seq_along(yt),function(x) { 
    
    a <- yt[[x]]%>%filter(year == 2015)%>%select(best)%>%as.numeric()
    })
  
  yt <- lapply(seq_along(yt), function(x)
    
    yt[[x]] <- yt[[x]]%>%
      mutate(eli_point = (best - bas[[x]])/bas[[x]]*100,
             eli_low = (q5 - bas[[x]])/bas[[x]]*100,
             eli_up = (q95 - bas[[x]])/bas[[x]]*100
             )%>%
      mutate(across(everything(),~format(round(., digits = 2), nsmall = 2)))%>%
      mutate(across(everything(),as.numeric))%>%
      mutate(eli90 = ifelse(eli_point<= -90, year,0),
             eli65 = ifelse(eli_point<= -65, year,0),
             eli80 = ifelse(eli_point <=-80, year,0))
  )
  names(yt) <- names(dt)
  return(yt)
}

####tablegeneration by scenario ####


Stable <- function(dt, paramDt, numbersample){
  HCV$numberSamples <- numbersample
  HCVInc_rangePar <- list()
  
  HCVInc_rangePar <- lapply(c(names(dt)), function(x) 
    sResCR(dt[[x]],paramDt[[x]], numberSamples = numbersample, 
           indicator = "incidence", pop = NULL))
  
  names(HCVInc_rangePar) <- names(scenarioResults) 
  
  #### newinfections #### 
  
  HCVInf_rangePar <- list()
  
  HCVInf_rangePar <- lapply(c(names(dt)), function(x) 
    sResCR(dt[[x]],paramDt[[x]], numberSamples = numbersample, 
           indicator = "newInfections", pop = NULL))
  
  names(HCVInf_rangePar) <- names(scenarioResults)
  
  ####HCV-related deaths ####
  HCVD_rangePar <- list()
  
  HCVD_rangePar <- lapply(c(names(dt)), function(x) 
    sResCR(dt[[x]],paramDt[[x]], numberSamples = numbersample, 
           indicator = "newHCVdeaths", pop = NULL))
  
  names(HCVD_rangePar) <- names(scenarioResults)
  
  #### HCVmor ####
  
  HCVMor_rangePar <- list()
  
  HCVMor_rangePar <- lapply(c(names(dt)), function(x) 
    sResCR(dt[[x]],paramDt[[x]], numberSamples = numbersample, 
           indicator = "HCVMortality", pop = NULL))
  
  names(HCVMor_rangePar) <- names(scenarioResults)
  
  ####numCure ####
  HCVCur_rangePar <- list()
  
  HCVCur_rangePar <- lapply(c(names(dt)), function(x) 
    sResCR(dt[[x]], paramDt[[x]], numberSamples = numbersample, 
           indicator = "newCured", pop = NULL))
  
  names(HCVCur_rangePar) <- names(scenarioResults)
  
  ####SVR ####
  HCVSVR_rangePar <- list()
  
  HCVSVR_rangePar <- lapply(c(names(dt)), function(x) 
    sResCR(dt[[x]], paramDt[[x]], numberSamples = numbersample, 
           indicator = "RSVR", pop = NULL))
  
  names(HCVSVR_rangePar) <- names(scenarioResults)
  
  #### combine #### 
  IncS <- SceTable(HCVInc_rangePar)
  
  numInfS <- SceTable(HCVInf_rangePar)
  numHCVDS <- SceTable(HCVD_rangePar)
  numHCVMorS <- SceTable(HCVMor_rangePar)
  numHCVCurS <- SceTable(HCVCur_rangePar)
  HCVSVRS <- SceTable(HCVSVR_rangePar)
  
  numHCVDS
  
  
  Target <- list()
  
  WHOTarget1 <- list()
  WHOTarget2 <- list()
  WHOTarget3 <- list()
  
  for (i in 1: length(numInfS)){ 
    
    WHOTarget1[[i]] <- IncS[[i]]%>%mutate(
      eli80 = ifelse(IncS[[i]]$eli80>0, eli80,0))%>%
      filter(eli80 != 0)%>%slice(1)%>%select(year)%>%as.numeric()
    
    
    WHOTarget2[[i]] <- numHCVDS[[i]]%>%mutate(
      eli65 = ifelse(numHCVDS[[i]]$eli65>0, eli65,0))%>%
      filter(eli65 != 0)%>%slice(1)%>%select(year)%>%as.numeric()
    
    
    WHOTarget3[[i]] <- HCVSVRS[[i]]%>%mutate(
      eli80 = ifelse(HCVSVRS[[i]]$best>=80, 1,0))%>%
      filter(eli80 != 0)%>%slice(1)%>%select(year)%>%as.numeric()
    
    
    Target[[i]] <- cbind(year = round(as.numeric(numInfS[[i]]$year),digits = 0),
                         scenario = names(numInfS)[i],
                         NumInf = paste0(numInfS[[i]]$best, 
                                         " (", numInfS[[i]]$q5, "-", 
                                         numInfS[[i]]$q95, ")"),
                         HCVInc = paste0(IncS[[i]]$best, 
                                         " (",IncS[[i]]$q5, "-", 
                                         IncS[[i]]$q95, ")"),
                         HCVIncReduction =  paste0(IncS[[i]]$eli_point, 
                                                   " (",IncS[[i]]$eli_low, "-", 
                                                   IncS[[i]]$eli_up, ")"),
                         WHOtar1 = WHOTarget1[[i]], 
                         
                         NumHCVD = paste0(numHCVDS[[i]]$best, 
                                          " (", numHCVDS[[i]]$q5, "-", 
                                          numHCVDS[[i]]$q95, ")"),
                         HCVDPrevent = paste0(numHCVDS[[i]]$eli_point, 
                                              " (",numHCVDS[[i]]$eli_low, "-", 
                                              numHCVDS[[i]]$eli_up, ")"),
                         WHOtar2 = WHOTarget2[[i]], 
                         
                         Cure = paste0(numHCVCurS[[i]]$best, 
                                       " (", numHCVCurS[[i]]$q5, "-", 
                                       numHCVCurS[[i]]$q95, ")"),
                         SVR = paste0(HCVSVRS[[i]]$best, 
                                      " (", HCVSVRS[[i]]$q5, "-", 
                                      HCVSVRS[[i]]$q95, ")"),
                         WHOtar3 = WHOTarget3[[i]]
                         
                         
    )%>%
      as.data.frame()%>%filter(year =="2030")
  
    }
  
  
  tt <- do.call(rbind, Target)
  return(tt)
} 
  
#### factor population names #### 

factorpop <- function(dt){ 
  dt$population <- factor(dt$population, 
         levels = c("HIV-", "HIV-PrEP", "HIV+", "HIV+d"), 
         labels = c("HIV-negative not on PrEP", 
                    "HIV-negative on PrEP", 
                    "HIV-positive undiagnosed",
                    "HIV-positive diagnosed (on treatment)"))
  return(dt)
  
  
  }



