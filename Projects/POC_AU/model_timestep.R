# each timesteps 
# compartments 
modres.t <- function(pg, Best, endYear, allp = NULL) {
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
    dplyr::select(-Var3)%>%ungroup()
  
  allpop <- allpop%>%
    filter(time!= 1)%>%
    mutate(time = c(rep(seq(pg$startYear, endYear- 2*pg$timestep,
                            pg$timestep),each = pg$npops*pg$ncomponent)))
  
  names(allpop) <- c("population", "state", "best","timestep", "cascade",
                     "disease_prog")
  
  ## MidyearIndex
  timelong <- seq(pg$startYear, endYear, pg$timestep) 
  
  allpop <- allpop%>%filter(timestep%in% timelong)%>%
    mutate(year = timestep%/%1)
}

# for flows 
modres.flow.t <- function(pg, Best, endYear, allp = NULL) {
  flowpop <- as.data.frame.table(Best[[allp]])
  
  flowpop  <- flowpop %>%
    mutate(dt = rep(seq(1.0,(endYear - pg$timestep), pg$timestep), 
                      each=pg$npops))%>%ungroup()%>%
    mutate(year = dt%/%1 - 1, 
           population = Var1, 
           best = Freq)%>%
    select(year, dt, population, best)
  return(flowpop)

  }
  
modres.flow.t (POC_AU, Sce[[1]], endYear = 100, allp = "newTreatment")



subpop_N <- modres.t(POC_AU, calibrateInit,endYear = endY)%>% as_tibble()%>%
  group_by(timestep, population)%>%
  summarise(best = sum(best))%>%arrange(timestep, population)
# all subpop in one list 

tempNOTInfectedRNA_subpop <- modres.t(POC_AU, calibrateInit,Population = POC_AU$popNames,
                                      Disease_prog = NULL, 
                                      Cascade = c("s", "cured"), 
                                      param = NULL ,endYear = endY)%>%
  as_tibble()%>%filter(cascade%in%c("s", "cured"))%>%group_by(timestep, population)%>%
  summarise(best = sum(best))%>%arrange(timestep, population)



tempPrevRNA_subpop <- cbind(timestep = subpop_N$timestep,
                            population = subpop_N$population,
                            
                            as.data.frame(100*(subpop_N[, -c(1,2)] - 
                                                 tempNOTInfectedRNA_subpop[ ,-c(1,2)])/ 
                                            subpop_N[ ,-c(1,2)]))%>%
  tibble::as_tibble()%>%mutate(year = timestep%/%1)%>%
  mutate(year = year + 2015 - 1)

tempPrevRNA_subpop <- tempPrevRNA_subpop%>%group_by(year, population)%>%
  summarize(mean(best)) 



subpop_N <- lapply(POC_AU$popNames, function(x){ 
  
  a <- popResults_MidYear(POC_AU, calibrateInit,
                          Population = x,
                          Disease_prog = NULL, 
                          Cascade = NULL, param = NULL, 
                          endYear = endY)%>%ungroup()
})

names(subpop_N) <- POC_AU$popNames

# all subpop in one list 
pop_N <- dplyr::bind_rows(subpop_N, .id = 'population')
pop_N <- pop_N%>%arrange(year, population)

tempNOTInfectedRNA_subpop <- popResults_MidYear(POC_AU, calibrateInit,Population = POC_AU$popNames,
                                                Disease_prog = NULL, 
                                                Cascade = c("s", "cured"), 
                                                param = NULL ,endYear = endY)%>%
  as_tibble()%>%group_by(year, population)%>%
  summarise(best = sum(best))%>%arrange(year, population)



tempPrevRNA_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                            population = POC_AU$popNames,
                            
                            as.data.frame(100*(pop_N[, -c(1,2)] - 
                                                 tempNOTInfectedRNA_subpop[ ,-c(1,2)])/ 
                                            pop_N[ ,-c(1,2)]))%>%
  tibble::as_tibble()
pop_labname <- c("PWID in community",  "Former PWID in community", 
                 "PWID in prisons",  "Former PWID in prisons", 
                 "nonPWID in prisons")
endY_plot <- 2030- POC_AU$cabY
tempPrevRNA_subpop <- tempPrevRNA_subpop%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))
HCVPrevRNA <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVPrevRNA_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCV.RNA.prevalence*100,
                           up = upper*100,
                           low = lower*100,
                           population = factor(population, 
                                               levels = POC_AU$popNames, 
                                               labels = pop_labname ))
popPrevRNAPlot <- indicatorPlot(POC_AU, tempPrevRNA_subpop, 
                                ylabel = "HCV prevalence (%)",
                                xlimits = c(POC_AU$startYear, 
                                            (POC_AU$startYear+endY_plot), 5),
                                calibration_Y = POC_AU$cabY,
                                rangeun = NULL, 
                                groupPlot = NULL, 
                                facetPlot = population,
                                observationData = HCVPrevRNA , 
                                simulateYear = POC_AU$simY) +
  ggtitle("HCV RNA prevalence by population") 

popPrevRNAPlot <- popPrevRNAPlot + 
  facet_custom (~population,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 60))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 50))),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 100))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 50))),
                    scale_new(5,
                              scale_y_continuous(limits = 
                                                   c(0, 5)))
                  )) + theme_bw()

popPrevRNAPlot 

varying_Yint <- 2016
varying_Yfir <- 2019
varying_Ymid <- 2020
varying_Yend <- 2021
calibration_Y <- 2015
varyingYpoint_int <- (varying_Yint - calibration_Y)/POC_AU$timestep 
varyingYpoint_fir <- (varying_Yfir - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_mid <- (varying_Ymid - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_end <- (varying_Yend - calibration_Y)/POC_AU$timestep + 1

for ( i in 2:dim(dfList$eta)[[2]]){  
  dfList$eta[3, i, c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[3]), 0.2, length = 1/POC_AU$timestep),
      seq(0.2, 0.999999, length = (varyingYpoint_fir - varyingYpoint_int)),
      seq(0.999999, as.numeric(TreatInit[3]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(TreatInit[3]), as.numeric(TreatInit[3]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(TreatInit[3]), POC_AU$npts - varyingYpoint_end))
  
  dfList$eta[4, i, c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[4]), 0.2, length = 1/POC_AU$timestep),
      seq(0.2, 0.999999, length = (varyingYpoint_fir - varyingYpoint_int)),
      seq(0.999999, as.numeric(TreatInit[4]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(TreatInit[4]), as.numeric(TreatInit[4]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(TreatInit[4]), POC_AU$npts - varyingYpoint_end))
  
  
  
  dfList$eta[5, i, c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[5]), 0.2, length = 1/POC_AU$timestep),
      seq(0.2, 0.999999, length = (varyingYpoint_fir - varyingYpoint_int )),
      seq(0.999999, as.numeric(TreatInit[5]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(TreatInit[5]), as.numeric(TreatInit[5]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(TreatInit[5]), POC_AU$npts - varyingYpoint_end))
  
  
  dfList$rho[3, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_rho[3]), 0.999999, length = (varyingYpoint_fir - varyingYpoint_int + 1)),
      seq(0.999999, as.numeric(Retreat[3]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(Retreat[3]), as.numeric(Retreat[3]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(Retreat[3]), POC_AU$npts - varyingYpoint_end))
  
  dfList$rho[4, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_rho[4]), 0.999999, length = (varyingYpoint_fir - varyingYpoint_int + 1)),
      seq(0.999999, as.numeric(Retreat[4]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(Retreat[4]), as.numeric(Retreat[4]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(Retreat[4]), POC_AU$npts - varyingYpoint_end))
  
  
  
  dfList$rho[5, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_rho[5]), 0.999999, length = (varyingYpoint_fir - varyingYpoint_int +1)),
      seq(0.999999, as.numeric(Retreat[5]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(Retreat[5]), as.numeric(Retreat[5]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(Retreat[5]), POC_AU$npts - varyingYpoint_end))
  
}




tic <- proc.time()
endY <- 100
calibrateInit <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                        disease_progress,pop_array,
                        dfList, fib, 
                        modelrun="UN", proj = "POC_AU", end_Y = endY)
calibrateFlow <- calibrateInit[!names(calibrateInit)%in%
                                 c("allPops", "newpop_tran", "newpop_tranState", "HCVdeathState",
                                   "newDeathState", "death_hcv")]
endY <- 100
flow_sub <- list()
flow_all <- list()

flow_sub <- lapply(names(calibrateFlow), function(x){ 
  a <- indicatorResults(POC_AU, calibrateFlow, x, 
                        pop=POC_AU$popNames,
                        paramR = NULL, range = NULL,
                        endY = endY)
})

names(flow_sub) <- names(calibrateFlow)

flow_setting <- lapply(flow_sub, function(x){ 
  
  a <- x%>%
    mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), 
                            "commu", "prisons"))
  
  a <- a%>%group_by(year)%>%summarise_at("best", sum)
}) 

N_treatment <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                each = 1),
                     as.data.frame(flow_setting$newTreatment[, -c(1)] + 
                                     flow_setting$newRetreat[, -c(1)] ))%>%
  tibble::as_tibble()



HCVtreatinitN_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatinitN_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = realpop,
                           up = upper,
                           low = lower)

HCVtreatinitN_setting_fit <- HCVtreatinitN_setting_fit%>%group_by(time)%>%
  summarize(realPop = sum(realpop), 
            up = sum(up), 
            low = sum(low))


N_treatment_setting_p <- indicatorPlot(POC_AU, N_treatment, 
                                       ylabel = "N",
                                       xlimits = c(POC_AU$startYear, 
                                                   (POC_AU$startYear+30), 5),
                                       calibration_Y = POC_AU$cabY,
                                       rangeun = NULL, 
                                       groupPlot = NULL, 
                                       facetPlot = NULL,
                                       observationData = HCVtreatinitN_setting_fit, 
                                       simulateYear = POC_AU$simY) + 
  ggtitle("Number of treatment init") + theme_bw()

N_treatment_setting_p <- N_treatment_setting_p + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 35000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 4000)))))
N_treatment_setting_p
