#### function to create aggegrate data columns: q5, q95, etc... ####
### drop it ###
popRange <- function(dt, target_year = NULL, pop = NULL, test= NULL, Casca = NULL){ 
  testy <- head(unique(dt$year), 1)
  
  if(isTRUE(testy < HCV$simY)){ 
    tt <- dt%>%mutate(year = year + HCV$cabY - 1)
    
  } else{ 
    
    tt <- dt
    
  }
  
  if(!is.null(target_year)){ 
    
    tt <- tt%>%
      filter(year == target_year)
    }
  
  if (is.null(pop) && is.null(test) && is.null(Casca)){
    tt <- tt%>%arrange(year, population, testing) 
  }else if(is.null(pop) && is.null(test) && !is.null(Casca)){
    tt <- tt%>%arrange(year, population, testing, Casca) 
  }else if (!is.null(pop) && is.null(test)){
    tt <- tt%>%arrange(year, testing)
  }else if (is.null(pop) && !is.null(test)) {
    tt <- tt%>%arrange(year, population) 
  } else if (!is.null(pop) && !is.null(test)) {
    tt <- tt%>%arrange(year)  
  }
  
  bestcol <- which(colnames(tt) == "best")
  popSizeQuan <- tt%>%gather(., "simulation", "estimate", 
                             bestcol:ncol(tt))%>%
    filter(simulation!= "best")
  
  if (is.null(pop) && is.null(test) && is.null(Casca)){
    popQ <- popSizeQuan%>%group_by(year, population, testing)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year, population, testing)
    
    tt[, c("min", "max", "Med", "Mu", "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
  } 
  else if (is.null(pop) && is.null(test) && !is.null(Casca)){
    popQ <- popSizeQuan%>%group_by(year, population, testing, Casca)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year, population, testing, Casca)%>% 
      select(everything())
    tt[, c("min", "max", "Med", "Mu", 
           "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
    }else if (!is.null(pop) && is.null(test) && !is.null(Casca)){
    popQ <- popSizeQuan%>%group_by(year, testing, Casca)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year, testing, Casca)%>% 
      select(everything())
    tt[, c("min", "max", "Med", "Mu", 
           "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
    }else if (!is.null(pop) && is.null(test)){
    popQ <- popSizeQuan%>%group_by(year, testing)%>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE))%>%
      ungroup()%>%arrange(year, testing)%>% 
      select(everything())
    tt[, c("min", "max", "Med", "Mu", 
           "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
  } 
  else if (is.null(pop) && !is.null(test)){
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
    tt[, c("min", "max", "Med", "Mu", 
           "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
    
  } 
  else if (!is.null(pop) && !is.null(test)){
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
    tt[, c("min", "max", "Med", "Mu", 
           "q5", "q25", "q75", "q95")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95)) 
  } 
  return(tt)
}


# number prevented calculation 
Numavert <- function(dt, pop = NULL, dt_bench){ 
  
  if(is.null(pop)){
    dt_avr <-  cbind(year = dt$year, 
                     population = dt$population,
                     testing = dt$testing, 
                     as.data.frame(
                       (dt[, -c(1:3)] - 
                          dt_bench[ ,-c(1:3)])))%>%as_tibble()
    
    
  }
  else{ 
    
    dt_all <- dt%>%group_by(year, testing)%>%
      summarize(across(c(best, colnames(colset)), ~sum(.x)))
    
    dt_bench_all <- dt_bench%>%group_by(year, testing)%>%
      summarize(across(c(best, colnames(colset)), ~sum(.x)))%>%
      arrange(year)
    
    dt_avr <-  cbind(year = dt_all$year, 
                     testing = dt_all$testing, 
                     as.data.frame(
                       (dt_all[, -c(1:2)] - 
                          dt_bench_all[ ,-c(1:2)])))%>%as_tibble()
    
    }
  return(dt_avr)
}



# function for combine parameter for sensitivity analysis 
## this is for constant parameters (matrix format)

combinParam <- function(dt, HCV, y_pt = NULL){ 
  
  x <- list()
  
    for(i in colnames(dt[[1]])){ 
      
      x[[i]] <- paste0(i, seq(1,HCV$npops,1))
      }

  
  
  tttt <- list()
  
  for(set in 1:HCV$numberSamples){
  
      tttt[[set]] <- dt[[set]]%>%as.matrix()%>%as.vector()
      
    
    
    
    
  }
  
  t_bind <- do.call("rbind", tttt)
  colnames(t_bind) <- do.call("cbind",x)%>%as.vector()
  
  return(t_bind)
  
}


#### result summary  ####
# epi 
Res_summary_epi <- function(pj_obj, dt, dt_param, endY, base_case){ 
  
  # @pjname: project objects
  # @dt: the outcome dataframe for results outputs 
  # @dt_param: the uncertainty range outcome dataframe for results outpots 
  # @endY: how many years results going to output
  # @base_case: if !is.null(base_case) run the output for base scenarios
  # such as HIVPrEPRate, HIVDRate 
  
  # this function will generate the epi results dataframe :
  # datafrmaes 
  # 1) Epidemiological outcomes (by subpops and overall pop): 
  #    A. Incidence (primary and reinfection)
  #    B. number of newinfections (primary and reinfection)
  #    C. number of primary newinfections 
  #    D. number of reinfections 
  # if base_case, results will generate additional populations results 
  #    A. % of subpops 
  #    B. number of total population yearly 
  #    C. 
  
  
  if(!is.null(base_case)){ 
    # populations 
    
    # number of total population 
    
    Allpop <- popResults_MidYear(pj_obj, dt ,Population = NULL,
                                 Disease_prog = NULL , 
                                 Cascade = NULL, param = dt_param, 
                                 endYear = endY)%>%as.data.frame()
    
    popS_range <- popResults_range(pj_obj, Allpop, Population = NULL,
                                   Disease_prog = NULL , 
                                   Cascade = NULL, end_Y = endY)
    
    # number of HIV infected 
    HIVInfected <- popResults_MidYear(pj_obj, dt,
                                      Population = c("HIV+", "HIV+d"),
                                      Disease_prog = NULL, 
                                      Cascade = NULL, param = dt_param, 
                                      endYear = endY)
    
    HIVInf <- HIVInfected%>%ungroup()%>% dplyr::group_by(year)%>%
      summarise_at(.vars = c(colnames(Allpop)[-1]),sum)
    
    
    HIVInfectedD <-  HIVInfected%>%filter(population == "HIV+d")
    
    HIVInfectedUND <-  HIVInfected%>%filter(population == "HIV+")
    
    HIVDiagRate <- cbind(year = seq(pj_obj$startYear , endY-1 ,1),
                         as.data.frame(HIVInfectedD[, -c(1,2)]/ HIVInf[ ,-1])*100)%>%
      tibble::as_tibble() 
    
    
    HIVDiagRate_range <- popResults_range(pj_obj, HIVDiagRate, Population = NULL,
                                          Disease_prog = NULL , 
                                          Cascade = NULL, end_Y = endY) 
    subpop_percent <- list()
    
    subpop_percent_range <- list()
    
    subpop_percent[["HIV+d"]] <- cbind(year = seq(pj_obj$startYear , endY-1 ,1),
                                       as.data.frame(HIVInfectedD[, -c(1,2)]/ Allpop[ ,-1])*100)%>%
      tibble::as_tibble()
    
    subpop_percent_range[["HIV+d"]] <- popResults_range(pj_obj, 
                                                        subpop_percent[["HIV+d"]], 
                                                        Population = NULL,
                                                        Disease_prog = NULL , 
                                                        Cascade = NULL, end_Y = endY) 
    
    subpop_percent[["HIV+"]] <- cbind(year = seq(pj_obj$startYear , endY-1 ,1),
                                      as.data.frame(HIVInfectedUND[, -c(1,2)]/ Allpop[ ,-1])*100)%>%
      tibble::as_tibble()
    
    subpop_percent_range[["HIV+"]] <- popResults_range(pj_obj, 
                                                       subpop_percent[["HIV+"]], 
                                                       Population = NULL,
                                                       Disease_prog = NULL , 
                                                       Cascade = NULL, 
                                                       end_Y = endY)
    
    # % of HIV negative 
    HIVN <- popResults_MidYear(pj_obj, dt,
                               Population = c("HIV-", "HIV-PrEP"),
                               Disease_prog = NULL, 
                               Cascade = NULL, param = dt_param, 
                               endYear = endY) 
    
    HIVNeg <- HIVN%>%ungroup()%>% dplyr::group_by(year)%>%
      summarise_at(.vars = c(colnames(Allpop)[-1]),sum)
    
    HIVPrEP <- HIVN%>%filter(population =="HIV-PrEP")
    
    HIVnotPrEP <- HIVN%>%filter(population =="HIV-")
    
    HIVPrEPRate <- cbind(year = seq(pj_obj$startYear , endY - 1 ,1),
                         as.data.frame(HIVPrEP[, -c(1,2)]/HIVNeg[ ,-1])*100)%>%
      tibble::as_tibble() 
    
    
    HIVPrEPRate_range <- popResults_range(pj_obj, HIVPrEPRate, Population = NULL,
                                          Disease_prog = NULL , 
                                          Cascade = NULL, end_Y = endY) 
    
    subpop_percent[["HIV-"]] <- cbind(year = seq(pj_obj$startYear , endY-1 ,1),
                                      as.data.frame(HIVnotPrEP[, -c(1,2)]/ Allpop[ ,-1])*100)%>%
      tibble::as_tibble()
    
    subpop_percent_range[["HIV-"]] <- popResults_range(pj_obj, 
                                                       subpop_percent[["HIV-"]], 
                                                       Population = NULL,
                                                       Disease_prog = NULL , 
                                                       Cascade = NULL, end_Y = endY) 
    
    subpop_percent[["HIV-PrEP"]] <- cbind(year = seq(pj_obj$startYear , endY - 1 ,1),
                                          as.data.frame(HIVPrEP[, -c(1,2)]/ Allpop[ ,-1])*100)%>%
      tibble::as_tibble()
    
    subpop_percent_range[["HIV-PrEP"]] <- popResults_range(pj_obj, 
                                                           subpop_percent[["HIV-PrEP"]], 
                                                           Population = NULL,
                                                           Disease_prog = NULL, 
                                                           Cascade = NULL, 
                                                           end_Y = endY)
    
    # Epi outcomes 
    # incidence 
    
    Inci <- scenario_Incidence_allset(list(dt), list(dt_param), 
                                      pop = NULL, statusQ = "y", 
                                      endY = endY)
    
    thres <- Inci%>%mutate(year = pj_obj$cabY - 1 + year)%>%
      filter(year == 2015)%>%
      select(best)%>%mutate(thres = best*0.2)
    
    Inci_pop <- scenario_Incidence_allset(list(dt), list(dt_param),
                                          pop = "pop", statusQ = "y", 
                                          endY = endY)
    
    threspop <- Inci_pop%>%mutate(year = pj_obj$cabY - 1  + year)%>%
      filter(year == 2015)%>%
      select(best)%>%mutate(thres = best*0.2)
    
    threspop[2,] <- threspop[1,]
    
    infnum_pop <- scenario_Incidence_allset(list(dt), list(dt_param),
                                            pop = "pop", statusQ = "y", 
                                            indicator = "newInfections", 
                                            endY = endY)
    
    reinfnum_pop <- scenario_Incidence_allset(list(dt), list(dt_param),
                                              pop = "pop", statusQ = "y", 
                                              indicator = "newreinfection", 
                                              endY = endY)
    
    
  }
  else{ 
    Inci_pop <- scenario_Incidence_allset(dt, 
                                          dt_param, 
                                          pop = "pop", endY = endY) 
    
    Inci <- scenario_Incidence_allset(dt, 
                                      dt_param,
                                      pop = NULL, endY = endY) 
    
    infnum_pop <- scenario_Incidence_allset(dt, 
                                            dt_param,
                                            pop = "pop", 
                                            indicator = "newInfections", 
                                            endY = endY)
    
    reinfnum_pop <- scenario_Incidence_allset(dt, dt_param,
                                              pop = "pop",
                                              indicator = "newreinfection", 
                                              endY = endY)
  }
  if(!is.null(base_case)){ 
    result = list(
      Inci_pop = Inci_pop,
      Inci = Inci,
      infnum_pop = infnum_pop,
      reinfnum_pop = reinfnum_pop,
      HIVDiagRate_range = HIVDiagRate_range,
      HIVPrEPRate = HIVPrEPRate,
      subpop_percent_range = subpop_percent_range)
  }else{ 
    result = list(
      Inci_pop = Inci_pop,
      Inci = Inci,
      infnum_pop = infnum_pop,
      reinfnum_pop = reinfnum_pop)
  }
  
}

# cost and qaly
Res_summary_costqaly <- function(pj_obj, dt, dt_param, endY){ 
  
    indicators <- c("costTestingAb", "costTestingAg", "costnewTestingPOCT", 
                    "costTreatment", "costCured", "costRetreat") 
    
    cost_flow <- list()
    
    cost_box <- list()
    
    cost_box_Range <-list() 
    # cost of flows
    for(i in indicators){
      cost_flow[[i]] <- indicatorResults(pj_obj, 
                                              dt,
                                              i, 
                                              dt_param,
                                              pop = NULL, 
                                              range = "y", 
                                              endY, scenario = NULL)
    }
    
    # cost and qaly in the components 
    indicators_allp <- c("QALYPops", "costPops")
    
    for(i in indicators_allp){
      cost_box[[i]] <- 
        popResults_MidYear(pj_obj, dt, 
                           Population = pj_obj$popNames,
                           Disease_prog = pj_obj$progress_name, 
                           Cascade = c("s",pj_obj$cascade_name), 
                           param = dt_param,
                           endYear = endY, 
                           allp = i)%>%tibble::as.tibble() 
      
      
     
    }
    
    # ouptut as list
    result = list(
      cost_flow = cost_flow,
      cost_box =  cost_box)
    
    
    
    }
  
  
  
