library(dplyr)
library(stringr)
# simulate the HCV model equations

HCVMSM <- function(HCV, parama, initialPop, disease_progress,
                   pop_array, param_cascade, fib,  end_Y = NULL, modelrun=NULL,
                   scenario = NULL, cost = NULL, costflow = NULL){
  
  # Args:
  #       HCV: a list containing project specifications 
  #       tp: the simulation t points in year
  #       parama: paramaeters did not change by time
  #       initialPop: initial population size
  #       
  #        
  #     pop_array:  
  #     param_cascade: a list of array includes those paramaeters that varies across 
  #               population and disease progress, (tau_ab, tau_ag, tau_poct, 
  #               eta, lota, rho, cured) 
  #     scenario: testing the reinfection reduction related to behavior intervnetions
  #                            
  # cost data (apply discount rate already)
  # costdt: cost in each state in each timestep 
  # costcas: cost in the cascade (e.g. cost of testing ab)
  
  # simulation time period 
  dt <- HCV$timestep
  
  
  if (is.null(end_Y)){ 
    
    years <- HCV$nyears
    
    npts <-HCV$npts
    
  } else if (!is.data.frame(parama)) {
    
    years <- HCV$startYear:end_Y
    pts <- seq(HCV$startYear, end_Y, by = dt)
    
    pts <- head(pts, -1)
    
    npts <- length(pts)
    
    parama <- lapply(parama, function(x) x[1:(npts + 1)])
    
    pop_array <- pop_array[ , , 1:npts] 
    
    # cost <- lapply(cost, function(x) x[, , 1:(npts+1)])
    
    # costflow <- costflow[ , , 1:npts]
    
    
  }else{
    
    years <- HCV$startYear:end_Y
    
    pts <- seq(HCV$startYear, end_Y, by = dt)
    
    pts <- head(pts, -1)
    
    npts <- length(pts)
    
    parama <- parama[1:(npts+1),]
    
    pop_array <- pop_array[ , , 1:npts] 
    
    # cost <- lapply(cost, function(x) x[, , 1:(npts+1)])
    
    # costflow <- costflow[ , , 1:npts]
    
    
  }
  
  
  
  npops <- HCV$npops
  
  ncomponent <- HCV$ncomponent
  
  ncomponentWOS <- HCV$ncomponentWOS
  
  nprogress <- HCV$nprogress
  
  ncascade <- HCV$ncascade
  
  componentName<- HCV$component_name
  
  componentName_late <- HCV$component_name[-c(1:37)] # from Dc stage
  
  ncomponentName_late <- length(componentName_late)
  
  componentNameWOS<- HCV$component_nameWOS
  
  popNames<- HCV$popNames
  
  dimNames <- list(popNames, componentName)
  
  progressName <- HCV$progress_name
  
  cascadeName <- HCV$cascade_name 
  
  diseaseprogress_n <- HCV$diseaseprogress_n
  
  diseaseprogress_Name <- HCV$diseaseprogress_Name 
  
  fib_n <- HCV$fib_n
  
  fib_name <-HCV$fibName
  
  allPops <- array(0, c(npops, ncomponent, npts), dimnames = dimNames) 
  
  newDeathState <- array(0, c(npops, ncomponent, npts), dimnames = dimNames)
  
  HCVdeathState <-  array(0,  c(npops, 20, npts), 
                          dimnames = list(popNames, c("dc_undiag", 
                                                      "dc_diag_ab",
                                                      "dc_diag_RNA", 
                                                      "dc_treat",
                                                      "dc_treat_f",
                                                      "hcc_undiag",
                                                      "hcc_diag_ab",
                                                      "hcc_diag_RNA",
                                                      "hcc_treat",
                                                      "hcc_treat_f",
                                                      "lt_undiag", 
                                                      "lt_diag_ab",
                                                      "lt_diag_RNA",
                                                      "lt_treat",
                                                      "lt_treat_f",
                                                      "plt_undiag", 
                                                      "plt_diag_ab",
                                                      "plt_diag_RNA",
                                                      "plt_treat",
                                                      "plt_treat_f")))
  
  allPops [ , , 1] <- initialPop 
  
  
  #### reinfection factor ####
  #### Reinfection factor ####
  
  reinfP <- matrix(0, ncol = npts + 1, nrow = npops)
  if(is.null(scenario)){ 
    reinfP[1,] <- 1
    reinfP[2,] <- 1
    reinfP[3,] <- 1
    reinfP[4,] <- 1
    
  } else{ 
    reinfP <- scenario
    
  }
  
  #### entry ####
  entry1 <- matrix(0, ncol = 1, nrow = 4)
  death <- matrix(0, ncol = ncomponent, nrow = npops, dimnames = dimNames)
  death_hcv <- matrix(0, ncol = ncomponentName_late, nrow = npops, 
                      dimnames = list(HCV$popNames, componentName_late))
  
  ####create result matrix####
  ResultMatrix <- matrix(0, ncol = npts, nrow = npops)
  rownames(ResultMatrix) <- popNames
  
  newS <- ResultMatrix
  
  newEntry <- ResultMatrix
  
  newDeath <-ResultMatrix
  
  newInfections <- ResultMatrix
  
  newHCVdeaths <- ResultMatrix
  
  newHCVdeathsState <- ResultMatrix
  
  newTreatment <- ResultMatrix
  
  newRetreat <- ResultMatrix
  
  newTestingAb <- ResultMatrix
  
  newTestingAg <- ResultMatrix
  
  newTestingPOCT <- ResultMatrix
  
  newCured <- ResultMatrix
  
  newreinfection <- ResultMatrix
  
  inflow <-ResultMatrix
  
  outflow <-ResultMatrix 
  
  if (!is.null(cost)){ 
    
    costTestingAb <- ResultMatrix
    costTestingAg <- ResultMatrix
    costnewTestingPOCT <- ResultMatrix
    costTreatment <- ResultMatrix
    costCured <- ResultMatrix
    costRetreat <- ResultMatrix 

    costPops <- array(0, c(npops, ncomponent, npts), dimnames = dimNames) 
    
    QALYPops <- array(0, c(npops, ncomponent, npts), dimnames = dimNames) 
    
    
    }
  
  #### pop_array ####
  popArray <- array(0, c(npops, npops, npts), dimnames = list(popNames, popNames))
  
  newarray <- popArray
  
  ####demographic paramaeters####
  # entry 
  
  entry <- matrix(0, ncol = npts + 1, nrow = npops)
  entry[1,] <- parama$entry
  entry[2,] <- 0
  entry[3,] <- 0
  entry[4,] <- 0
  entry_dt <- entry*dt  
  
  ####background mortality =rate to leave the model####
  ## 
  morb <- matrix(0, ncol = npts + 1, nrow = npops)
  morb[1, ] <- parama$mortality_b
  morb[2, ] <- parama$mortality_b
  morb[3, ] <- parama$HIV_mor
  morb[4, ] <- parama$HIVd_mor
  morb_dt <-1-(1-morb)^dt
  
  ####HCV_related mortality####
  ## DC, HCC, LT, PLT, DCcured, HCCcured, LTcured, PLTcured
  
  #####DC#####
  mordc <- matrix(0, ncol = npts + 1, nrow = npops)
  mordc[1, ] <- parama$mordc
  mordc[2, ] <- parama$mordc
  mordc[3, ] <- parama$mordc*parama$oddsOfmordc_HIVp
  mordc[4, ] <- parama$mordc*parama$oddsOfmordc_HIVp
  mordc_dt <- 1-(1-mordc)^dt
  
  #####HCC#####
  morhcc <- matrix(0, ncol = npts + 1, nrow = npops)
  morhcc[1, ] <- parama$morhcc
  morhcc[2, ] <- parama$morhcc
  morhcc[3, ] <- parama$morhcc*parama$oddsOfmorhcc_HIVp
  morhcc[4, ] <- parama$morhcc*parama$oddsOfmorhcc_HIVp
  morhcc_dt <- 1-(1-morhcc)^dt
  
  #####LT####
  morlt <- matrix(0, ncol = npts + 1, nrow = npops)
  morlt[1, ] <- parama$morlt
  morlt[2, ] <- parama$morlt
  morlt[3, ] <- parama$morlt
  morlt[4, ] <- parama$morlt
  morlt_dt <- 1-(1-morlt)^dt
  
  #####PLT#####
  morplt <- matrix(0, ncol = npts + 1, nrow = npops)
  morplt[1, ] <- parama$morplt
  morplt[2, ] <- parama$morplt
  morplt[3, ] <- parama$morplt
  morplt[4, ] <- parama$morplt
  morplt_dt <- 1-(1-morplt)^dt
  
  #####DC cured####
  mordcCure <- matrix(0, ncol = npts + 1, nrow = npops)
  mordcCure[1, ] <- mordc[1,]*parama$Cure_mordc_Reduction
  mordcCure[2, ] <- mordc[2,]*parama$Cure_mordc_Reduction
  mordcCure[3, ] <- mordc[3,]*parama$Cure_mordc_Reduction
  mordcCure[4, ] <- mordc[4,]*parama$Cure_mordc_Reduction
  mordcCure_dt <- 1-(1-mordcCure)^dt
  
  #####HCC cured####
  morhccCure <- matrix(0, ncol = npts + 1, nrow = npops)
  morhccCure[1, ]<- morhcc[1,]*(1-parama$Cure_morhcc_Reduction)
  morhccCure[2, ]<- morhcc[2,]*(1-parama$Cure_morhcc_Reduction)
  morhccCure[3, ]<- morhcc[3,]*(1-parama$Cure_morhcc_Reduction)
  morhccCure[4, ]<- morhcc[4,]*(1-parama$Cure_morhcc_Reduction)
  morhccCure_dt <- 1-(1-morhccCure)^dt
  
  ## temporarily ignore HIV excess HIV  
  
  
  #### disease progress ####
  
  transition <-matrix(0,ncol = diseaseprogress_n, nrow = npops)
  colnames(transition) <- diseaseprogress_Name
  
  trans <- as.matrix(disease_progress, nrow = npops, ncol= diseaseprogress_n)
  
  transition[1, ]<- trans[1, ]
  transition[2, ]<- trans[2, ]
  transition[3, ]<- trans[3, ]
  transition[4, ]<- trans[4, ] 
  
  transition[3, "f3_hcc"] <- trans[3 ,"f3_hcc"]*parama$HIV_undiag_diseaseprog[1]
  transition[4, "f3_hcc"] <- trans[4,"f3_hcc"]*parama$HIV_diag_diseaseprog[1]
  
  transition[3, "f4_hcc"] <- trans[3 ,"f4_hcc"]*parama$HIV_undiag_diseaseprog[1]
  transition[4, "f4_hcc"] <- trans[4,"f4_hcc"]*parama$HIV_diag_diseaseprog[1]
  
  transition[3, "dc_hcc"] <- trans[3 ,"dc_hcc"]*parama$HIV_undiag_diseaseprog[1]
  transition[4, "dc_hcc"] <- trans[4,"dc_hcc"]*parama$HIV_diag_diseaseprog[1]
  
  transition[3, "dc_lt"] <- trans[3 ,"dc_lt"]*parama$HIV_undiag_diseaseprog[1]
  transition[4, "dc_lt"] <- trans[4,"dc_lt"]*parama$HIV_diag_diseaseprog[1]
  
  transition[3, "hcc_lt"] <- trans[3 ,"hcc_lt"]*parama$HIV_undiag_diseaseprog[1]
  transition[4, "hcc_lt"] <- trans[4,"hcc_lt"]*parama$HIV_diag_diseaseprog[1]
  
  transition_dt <- 1-(1-transition)^dt 
  
  ####Lt>PLT & cured stage disease progression####
  fibprog <-matrix(0,ncol = fib_n, nrow = npops)
  colnames(fibprog) <- fib_name
  fibprog <- as.matrix(fib, nrow = npops, ncol= fib_n)
  fibprog_dt <- 1-(1-fibprog)^dt 
  
  
  # care cascade; please see Fig "B. HCV care cascade progress" in 
  # 02. document/ 01. Model structure 
  # to better understand the paramaeter as following:
  
  
  ####spontaneously clearance####
  
  spc1 <-  matrix(0, ncol = npts + 1, nrow = npops)
  spc1[1, ] <- parama$spc1
  spc1[2, ] <- parama$spc2
  spc1[3, ] <- parama$spc3
  spc1[4, ] <- parama$spc4
  spc1_dt <- 1-(1-spc1)^dt 
  
  ####testing rate#### 
  #####antibody####
  dimN <- list(popNames, progressName) 
  tau_ab <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  tau_ab <- param_cascade$tau_ab
  tau_ab_dt <- 1-(1-tau_ab)^dt
  
  #####ag/RNA (second step testing)####
  tau_ag <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  tau_ag <- param_cascade$tau_ag
  tau_ag_dt <- 1-(1-tau_ag)^dt
  
  #####POCT: undiag>>diag_RNA####
  tau_poct <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  tau_poct <- param_cascade$tau_poct
  tau_poct_dt <- 1-(1-tau_poct)^dt
  
  ####treatment####
  eta <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  eta<- param_cascade$eta
  eta_dt <- 1-(1-eta)^dt
  
  ####treatment failed#### 
  lota <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  lota <- param_cascade$lota
  lota_dt <- 1-(1-lota)^dt
  
  ####re-treated####
  rho <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  rho <- param_cascade$rho
  rho_dt <- 1-(1-rho)^dt
  
  ####cured#### 
  cure <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  cure <- param_cascade$cured
  cure_dt <- 1-(1-cure)^dt
  
  ####population transition#### 
  pop_array <- 1-(1-pop_array)^dt 
  
  
  
  #------------------------------------------------------------------------------#
  ####                     Equations                                        ####   
  #------------------------------------------------------------------------------#
  for (t in 2: npts){ 
    
    oldPop <- allPops[, , t-1]
    newPop <- allPops[, , t]
    
    #####FOI#####
    
    foi <- matrix(0, ncol = 1, nrow = npops)
    
    foi[1,] <- parama$beta1[t]
    
    foi[2, ] <- parama$beta2[t]
    
    
    foi[3, ] <- parama$beta3[t]
    
    
    foi[4, ] <- parama$beta4[t]
    
    
    foi_dt <- foi*dt 
    
    
   
   
    
    ##### mortality ####
    death[,] <- morb_dt[,t]*oldPop[,] 
    
    death_hcv[,"dc_undiag"] <- mordc_dt[, t]*oldPop[, "dc_undiag"]
    death_hcv[,"dc_diag_ab"] <- mordc_dt[, t]*oldPop[, "dc_diag_ab"]
    death_hcv[,"dc_diag_RNA"] <- mordc_dt[, t]*oldPop[, "dc_diag_RNA"]
    death_hcv[,"dc_treat"] <- mordc_dt[,t]*oldPop[, "dc_treat"]
    death_hcv[,"dc_treat_f"] <- mordc_dt[,t]*oldPop[, "dc_treat_f"]
    death_hcv[,"dc_cured"] <- mordcCure_dt[,t]*oldPop[, "dc_cured"]                      
    death_hcv[,"hcc_undiag"] <- morhcc_dt[, t]*oldPop[, "hcc_undiag"]                       
    death_hcv[,"hcc_diag_ab"] <- morhcc_dt[, t]*oldPop[, "hcc_diag_ab"]                       
    death_hcv[,"hcc_diag_RNA"] <- morhcc_dt[,t]*oldPop[, "hcc_diag_RNA"]
    death_hcv[,"hcc_treat"] <- morhcc_dt[,t]*oldPop[, "hcc_treat"]
    death_hcv[,"hcc_treat_f"] <- morhcc_dt[,t]*oldPop[, "hcc_treat_f"]
    death_hcv[,"hcc_cured"] <- morhccCure_dt[,t]*oldPop[, "hcc_cured"]
    death_hcv[,"lt_undiag"] <- morlt_dt[, t]*oldPop[, "lt_undiag"]
    death_hcv[,"lt_diag_ab"] <- morlt_dt[, t]*oldPop[, "lt_diag_ab"]
    death_hcv[,"lt_diag_RNA"] <- morlt_dt[,t]*oldPop[, "lt_diag_RNA"]
    death_hcv[,"lt_treat"] <- morlt_dt[,t]*oldPop[, "lt_treat"]
    death_hcv[,"lt_treat_f"] <- morlt_dt[,t]*oldPop[, "lt_treat_f"]
    death_hcv[,"lt_cured"] <- morlt_dt[,t]*oldPop[, "lt_cured"]
    death_hcv[,"plt_undiag"] <- morplt_dt[, t]*oldPop[, "plt_undiag"]
    death_hcv[,"plt_diag_ab"] <- morplt_dt[, t]*oldPop[, "plt_diag_ab"]
    death_hcv[,"plt_diag_RNA"] <- morplt_dt[,t]*oldPop[, "plt_diag_RNA"] 
    death_hcv[,"plt_treat"] <- morplt_dt[,t]*oldPop[, "plt_treat"]
    death_hcv[,"plt_treat_f"] <- morplt_dt[,t]*oldPop[, "plt_treat_f"]
    death_hcv[,"plt_cured"] <- morplt_dt[,t]*oldPop[, "plt_cured"]
    
    
    
    #### proportion of cured ####
    
    # total num
    N <- matrix(0, ncol = 1, nrow = npops)
    N[1, ] <- sum(oldPop[1, ], na.rm = TRUE)
    N[2, ] <- sum(oldPop[2, ], na.rm = TRUE)
    N[3, ] <- sum(oldPop[3, ], na.rm = TRUE)
    N[4, ] <- sum(oldPop[4, ], na.rm = TRUE)
    
   
    # S 
    S <- matrix(0, ncol = 1, nrow = npops)
    S[1, ] <- (oldPop[1 , "s"] + 
                 oldPop[1 , "a_cured"] +
                 oldPop[1 , "f0_cured"] +
                 oldPop[1 , "f1_cured"] +
                 oldPop[1 , "f2_cured"] +
                 oldPop[1 , "f3_cured"] +
                 oldPop[1 , "f4_cured"] +
                 oldPop[1 , "dc_cured"] +
                 oldPop[1 , "hcc_cured"] +
                 oldPop[1 , "lt_cured"] +
                 oldPop[1 , "plt_cured"])
    S[2, ] <- (oldPop[2 , "s"] + 
                    oldPop[2 , "a_cured"] +
                    oldPop[2 , "f0_cured"] +
                    oldPop[2 , "f1_cured"] +
                    oldPop[2 , "f2_cured"] +
                    oldPop[2 , "f3_cured"] +
                    oldPop[2 , "f4_cured"] +
                    oldPop[2 , "dc_cured"] +
                    oldPop[2 , "hcc_cured"] +
                    oldPop[2 , "lt_cured"] +
                    oldPop[2 , "plt_cured"])
    S[3, ] <- (oldPop[3 , "s"] + 
                    oldPop[3 , "a_cured"] +
                    oldPop[3 , "f0_cured"] +
                    oldPop[3 , "f1_cured"] +
                    oldPop[3 , "f2_cured"] +
                    oldPop[3 , "f3_cured"] +
                    oldPop[3 , "f4_cured"] +
                    oldPop[3 , "dc_cured"] +
                    oldPop[3 , "hcc_cured"] +
                    oldPop[3 , "lt_cured"] +
                    oldPop[3 , "plt_cured"])
    S[4, ] <- (oldPop[4 , "s"] + 
                    oldPop[4 , "a_cured"] +
                    oldPop[4 , "f0_cured"] +
                    oldPop[4 , "f1_cured"] +
                    oldPop[4 , "f2_cured"] +
                    oldPop[4 , "f3_cured"] +
                    oldPop[4 , "f4_cured"] +
                    oldPop[4 , "dc_cured"] +
                    oldPop[4 , "hcc_cured"] +
                    oldPop[4 , "lt_cured"] +
                    oldPop[4 , "plt_cured"])
    
    Sall <- matrix(0, ncol = 1, nrow = 1)
    Sall <- sum(S[1,], S[2,], S[3,], S[4,])
    
    
    Ps<- matrix(0, ncol = 1, nrow = npops)
    Pca <- matrix(0, ncol = 1, nrow = npops)
    Pcf0 <- matrix(0, ncol = 1, nrow = npops)
    Pcf1 <- matrix(0, ncol = 1, nrow = npops)
    Pcf2 <- matrix(0, ncol = 1, nrow = npops)
    Pcf3 <- matrix(0, ncol = 1, nrow = npops)
    Pcf4 <- matrix(0, ncol = 1, nrow = npops)
    Pcdc <- matrix(0, ncol = 1, nrow = npops)
    Pchcc <- matrix(0, ncol = 1, nrow = npops)
    Pclt <- matrix(0, ncol = 1, nrow = npops)
    Pcplt <- matrix(0, ncol = 1, nrow = npops)
    
    # S should be all S population 
    Ps[, ] <- oldPop[, "s"]/S[,]
    Pca[,] <- oldPop[, "a_cured"]/S[,]
    Pcf0[,] <- oldPop[, "f0_cured"]/S[,]
    Pcf1[,] <- oldPop[, "f1_cured"]/S[,]
    Pcf2[,] <- oldPop[, "f2_cured"]/S[,]
    Pcf3[,] <- oldPop[, "f3_cured"]/S[,]
    Pcf4[,] <- oldPop[, "f4_cured"]/S[,]
    Pcdc[,] <- oldPop[, "dc_cured"]/S[,]
    Pchcc[,] <- oldPop[, "hcc_cured"]/S[,]
    Pclt[,] <- oldPop[, "lt_cured"]/S[,]
    Pcplt[,] <- oldPop[, "plt_cured"]/S[,]
    
    Ps[is.nan(Ps)] <- 0
    Pca[is.nan(Pca)] <- 0
    Pcf0[is.nan(Pcf0)] <- 0
    Pcf1[is.nan(Pcf1)] <- 0
    Pcf2[is.nan(Pcf2)] <- 0
    Pcf3[is.nan(Pcf3)] <- 0
    Pcf4[is.nan(Pcf4)] <- 0
    Pcdc[is.nan(Pcdc)] <- 0
    Pchcc[is.nan(Pchcc)] <- 0
    Pclt[is.nan(Pclt)] <- 0
    Pcplt[is.nan(Pcplt)] <- 0
    
    
    # I
    I <- matrix(0, ncol = 1, nrow = npops)
    
    I[1, ] <- N[1, ] - S[1,]
    I[2, ] <- N[2, ] - S[2,]
    I[3, ] <- N[3, ] - S[3,]
    I[4, ] <- N[4, ] - S[4,]
    
    
    Iall <- matrix(0, ncol = 1, nrow = 1)
    Iall <- sum(I[1,], I[2,], I[3,], I[4,])
    
    II <- Iall/Sall
    II[is.nan(II)] <- 0
    
    
    #####Entry of new people####
    
    if(modelrun == "steady"){ 
      entry1[1, ] <- sum(death)+sum(death_hcv)
      
    } else { entry1[1,] <- entry_dt[1,t]
    
    
    }
    
    #if(!is.null(scenario)){ reinfP[ , ] <- reinfP[ ,t]}
    #else{ reinfP[ ,t] <- c(1,1,1,1)}
    
    
    #####S#### 
    newPop[,"s"] <- entry1[,] + oldPop[,"s"] - death[,"s"] + 
      (colSums(pop_array[, , t]*oldPop[, "s"]) -
         rowSums(pop_array[, , t]*oldPop[, "s"])) -
      #foi_dt[, ]*I[,]*Ps[,]
      foi_dt[, ]*II*oldPop[, "s"]
    
    ##### a ####
    # a_infection 
    
    newPop[, "a_undiag"] <- oldPop[ ,"a_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "a_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_undiag"])) - 
      death[,"a_undiag"] -
      transition_dt[, "a_f0"]*oldPop[, "a_undiag"] - 
      spc1_dt[, t]*oldPop[ ,"a_undiag"] -
      tau_ab_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[, "a_undiag"] - 
      tau_poct_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[, "a_undiag"] + 
      #foi_dt[, ]*I[,]*Ps[,] +
      #reinfP[, t]*foi_dt[, ]*I[,]*Pca[,] + 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf0[,]
      foi_dt[, ]*II*oldPop[, "s"]+ 
      reinfP[, t]*foi_dt[, ]*II*oldPop[, "a_cured"] +
      reinfP[, t]*foi_dt[, ]*II*oldPop[, "f0_cured"]
    
    
    
    ## a_testing, ab+
    newPop[, "a_diag_ab"] <- oldPop[,"a_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "a_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_diag_ab"])) -
      death[,"a_diag_ab"] - 
      transition_dt[, "a_f0"]*oldPop[,"a_diag_ab"] + 
      tau_ab_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[,"a_undiag"] - 
      tau_ag_dt[, "a", t]*oldPop[,"a_diag_ab"]
    
    ## a_testing, ag+
    newPop[, "a_diag_RNA"] <- oldPop[ ,"a_diag_RNA"] +
      (colSums(pop_array[, , t]*oldPop[, "a_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_diag_RNA"])) -
      death[,"a_diag_RNA"]  -
      transition_dt[, "a_f0"]*oldPop[, "a_diag_RNA"] + 
      tau_ag_dt[, "a", t]*oldPop[, "a_diag_ab"] +
      tau_poct_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[, "a_undiag"] -
      eta_dt[, "a", t]*oldPop[,"a_diag_RNA"]
    
    ## a_testing, treat
    newPop[, "a_treat"] <- oldPop[,"a_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "a_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_treat"])) - 
      death[,"a_treat"] -
      transition_dt[, "a_f0"]*oldPop[,"a_treat"] +
      eta_dt[, "a", t]*oldPop[,"a_diag_RNA"] - 
      lota_dt[, "a", t]*(1-cure_dt[, "a", t])*oldPop[,"a_treat"] -
      cure_dt[, "a", t]*oldPop[,"a_treat"] +
      rho_dt[, "a", t]*oldPop[,"a_treat_f"]
    
    # a_treat_failed  
    newPop[, "a_treat_f"] <- oldPop[,"a_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "a_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_treat_f"])) - 
      death[,"a_treat_f"] - 
      transition_dt[, "a_f0"]*oldPop[,"a_treat_f"] -
      rho_dt[, "a", t]*oldPop[,"a_treat_f"] +
      lota_dt[, "a", t]*(1-cure_dt[, "a", t])*oldPop[,"a_treat"]
    
    # a_cured 
    newPop[, "a_cured"] <- oldPop[,"a_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "a_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_cured"])) - 
      death[,"a_cured"] + 
      cure_dt[, "a", t]*oldPop[, "a_treat"] + 
      spc1_dt[, t]*oldPop[ ,"a_undiag"] -
      #reinfP[, t]*foi_dt[, ]*I[,]*Pca[,]
      reinfP[, t]*foi_dt[, ]*II*oldPop[, "a_cured"]
    
    
    ##### f0 #####
    ## f0 infection
    newPop[, "f0_undiag"] <- oldPop[ ,"f0_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_undiag"])) - 
      death[,"f0_undiag"] +
      transition_dt[, "a_f0"]*oldPop[, "a_undiag"] -
      transition_dt[, "f0_f1"]*oldPop[, "f0_undiag"] -
      tau_ab_dt[, "f0", t]*oldPop[, "f0_undiag"] -  
      tau_poct_dt[, "f0", t]*oldPop[, "f0_undiag"] 
    
    
    ## f0 testing, ab+
    newPop[, "f0_diag_ab"] <- oldPop[,"f0_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_diag_ab"])) - 
      death[,"f0_diag_ab"] + 
      transition_dt[, "a_f0"]*oldPop[, "a_diag_ab"] -
      transition_dt[, "f0_f1"]*oldPop[, "f0_diag_ab"] +
      tau_ab_dt[, "f0", t]*oldPop[, "f0_undiag"] -
      tau_ag_dt[, "f0", t]*oldPop[, "f0_diag_ab"]    
    
    
    ## f0 testing, ag+/RNA
    newPop[, "f0_diag_RNA"] <- oldPop[ ,"f0_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_diag_RNA"])) -
      death[,"f0_diag_RNA"] +
      transition_dt[, "a_f0"]*oldPop[, "a_diag_RNA"] -
      transition_dt[, "f0_f1"]*oldPop[, "f0_diag_RNA"] + 
      tau_ag_dt[, "f0", t]*oldPop[, "f0_diag_ab"] +
      tau_poct_dt[, "f0", t]*oldPop[, "f0_undiag"] -
      eta_dt[, "f0", t]*oldPop[,"f0_diag_RNA"]    
    
    
    ## f0 treat
    newPop[, "f0_treat"] <- oldPop[,"f0_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_treat"])) - 
      death[,"f0_treat"] +
      transition_dt[, "a_f0"]*oldPop[, "a_treat"] - 
      transition_dt[, "f0_f1"]*oldPop[, "f0_treat"] +  
      eta_dt[, "f0", t]*oldPop[,"f0_diag_RNA"] - 
      lota_dt[, "f0", t]*(1-cure_dt[, "f0", t])*oldPop[,"f0_treat"] -  
      cure_dt[, "f0", t]*oldPop[,"f0_treat"] +  
      rho_dt[, "f0", t]*oldPop[,"f0_treat_f"]
    
    ## f0 treat failed
    newPop[, "f0_treat_f"] <- oldPop[,"f0_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_treat_f"])) -
      death[,"f0_treat_f"] + 
      transition_dt[, "a_f0"]*oldPop[, "a_treat_f"] - 
      transition_dt[, "f0_f1"]*oldPop[, "f0_treat_f"] + 
      lota_dt[, "f0", t]*(1-cure_dt[, "f0", t])*oldPop[,"f0_treat"] - 
      rho_dt[, "f0", t]*oldPop[,"f0_treat_f"]
    
    ## f0 cured
    newPop[, "f0_cured"] <- oldPop[,"f0_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_cured"]))  - 
      death[,"f0_cured"] + 
      cure_dt[, "f0", t]*oldPop[, "f0_treat"] -
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf0[,]
      reinfP[, t]*foi_dt[, ]*II*oldPop[, "f0_cured"]
    
    ##### f1 #####
    ## f1 infection
    newPop[, "f1_undiag"] <- oldPop[,"f1_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_undiag"])) - 
      death[,"f1_undiag"]  +                                           
      transition_dt[, "f0_f1"]*oldPop[, "f0_undiag"] - 
      transition_dt[, "f1_f2"]*oldPop[, "f1_undiag"] - 
      tau_ab_dt[, "f1", t]*oldPop[, "f1_undiag"] - 
      tau_poct_dt[, "f2", t]*oldPop[, "f1_undiag"] +
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf1[,]
      reinfP[, t]*foi_dt[, ]*II*oldPop[, "f1_cured"]
    
    
    ## f1 testing, ab+
    newPop[, "f1_diag_ab"] <- oldPop[ ,"f1_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_diag_ab"])) - 
      death[,"f1_diag_ab"] + 
      transition_dt[, "f0_f1"]*oldPop[, "f0_diag_ab"] - 
      transition_dt[, "f1_f2"]*oldPop[, "f1_diag_ab"] + 
      tau_ab_dt[, "f1", t]*oldPop[, "f1_undiag"] - 
      tau_ag_dt[, "f1", t]*oldPop[, "f1_diag_ab"]    
    
    
    ## f1 testing, ag+/RNA
    newPop[, "f1_diag_RNA"] <-  oldPop[ ,"f1_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_diag_RNA"])) - 
      death[,"f1_diag_RNA"] + 
      transition_dt[, "f0_f1"]*oldPop[, "f0_diag_RNA"] - 
      transition_dt[, "f1_f2"]*oldPop[, "f1_diag_RNA"] + 
      tau_ag_dt[, "f1", t]*oldPop[, "f1_diag_ab"] + 
      tau_poct_dt[, "f1", t]*oldPop[, "f1_undiag"] - 
      eta_dt[, "f1", t]*oldPop[,"f1_diag_RNA"]    
    
    
    ## f1 treat
    newPop[, "f1_treat"] <- oldPop[ ,"f1_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_treat"])) - 
      death[,"f1_treat"] + 
      transition_dt[, "f0_f1"]*oldPop[, "f0_treat"] - 
      transition_dt[, "f1_f2"]*oldPop[, "f1_treat"] +  
      eta_dt[, "f1", t]*oldPop[,"f1_diag_RNA"]  - 
      lota_dt[, "f1", t]*(1-cure_dt[, "f1", t])*oldPop[,"f1_treat"] -  
      cure_dt[, "f1", t]*oldPop[,"f1_treat"] +  
      rho_dt[, "f1", t]*oldPop[,"f1_treat_f"]
    
    ## f1 treat failed
    newPop[, "f1_treat_f"] <- oldPop[,"f1_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_treat_f"])) - 
      death[,"f1_treat_f"] + 
      transition_dt[, "f0_f1"]*oldPop[, "f0_treat_f"] - 
      transition_dt[, "f1_f2"]*oldPop[, "f1_treat_f"] + 
      lota_dt[, "f1", t]*(1-cure_dt[, "f1", t])*oldPop[,"f1_treat"] - 
      rho_dt[, "f1", t]*oldPop[,"f1_treat_f"]
    
    ## f1 cured
    newPop[, "f1_cured"] <- oldPop[,"f1_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_cured"])) - 
      death[,"f1_cured"] + 
      cure_dt[, "f1", t]*oldPop[, "f1_treat"] -
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf1[,] 
      reinfP[, t]*foi_dt[, ]*II*oldPop[, "f1_cured"]
    
    ##### f2 #####
    ## f2 infection
    newPop[, "f2_undiag"] <- oldPop[,"f2_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_undiag"])) - 
      death[,"f2_undiag"]  +                                                                                             
      transition_dt[, "f1_f2"]*oldPop[, "f1_undiag"] - 
      transition_dt[, "f2_f3"]*oldPop[, "f2_undiag"] - 
      tau_ab_dt[, "f2", t]*oldPop[, "f2_undiag"] -  
      tau_poct_dt[, "f2", t]*oldPop[, "f2_undiag"] +
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf2[,]
      reinfP[, t]*foi_dt[, ]*II*oldPop[, "f2_cured"]
    
    ## f2 testing, ab+
    newPop[, "f2_diag_ab"] <- oldPop[,"f2_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_diag_ab"])) - 
      death[,"f2_diag_ab"] + 
      transition_dt[, "f1_f2"]*oldPop[, "f1_diag_ab"] - 
      transition_dt[, "f2_f3"]*oldPop[, "f2_diag_ab"] + 
      tau_ab_dt[, "f2", t]*oldPop[, "f2_undiag"] - 
      tau_ag_dt[, "f2", t]*oldPop[, "f2_diag_ab"]    
    
    
    ## f2 testing, ag+/RNA
    newPop[, "f2_diag_RNA"] <- oldPop[,"f2_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_diag_RNA"])) - 
      death[,"f2_diag_RNA"] + 
      transition_dt[, "f1_f2"]*oldPop[, "f1_diag_RNA"] - 
      transition_dt[, "f2_f3"]*oldPop[, "f2_diag_RNA"] + 
      tau_ag_dt[, "f2", t]*oldPop[, "f2_diag_ab"] + 
      tau_poct_dt[, "f2", t]*oldPop[, "f2_undiag"] - 
      eta_dt[, "f2", t]*oldPop[,"f2_diag_RNA"]    
    
    
    ## f2 treat
    newPop[, "f2_treat"] <- oldPop[,"f2_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_treat"])) - 
      death[,"f2_treat"] + 
      transition_dt[, "f1_f2"]*oldPop[, "f1_treat"] - 
      transition_dt[, "f2_f3"]*oldPop[, "f2_treat"] + 
      eta_dt[, "f2", t]*oldPop[,"f2_diag_RNA"] - 
      lota_dt[, "f2", t]*(1-cure_dt[, "f2", t])*oldPop[,"f2_treat"] - 
      cure_dt[, "f2", t]*oldPop[,"f2_treat"] + 
      rho_dt[, "f2", t]*oldPop[,"f2_treat_f"]
    
    ## f2 treat failed
    newPop[, "f2_treat_f"] <- oldPop[,"f2_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_treat_f"]))  - 
      death[,"f2_treat_f"] + 
      transition_dt[, "f1_f2"]*oldPop[, "f1_treat_f"] - 
      transition_dt[, "f2_f3"]*oldPop[, "f2_treat_f"] + 
      lota_dt[, "f2", t]*(1-cure_dt[, "f2", t])*oldPop[,"f2_treat"] - 
      rho_dt[, "f2", t]*oldPop[,"f2_treat_f"]
    
    ## f2 cured
    newPop[, "f2_cured"] <- oldPop[,"f2_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_cured"])) - 
      death[,"f2_cured"] + 
      cure_dt[, "f2", t]*oldPop[, "f2_treat"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf2[,] 
      reinfP[, t]*foi_dt[, ]*II*oldPop[, "f2_cured"]
    
    ##### f3 #####
    ## f3 infection
    newPop[, "f3_undiag"] <- oldPop[,"f3_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_undiag"])) - 
      death[,"f3_undiag"] + 
      
      transition_dt[, "f2_f3"]*oldPop[, "f2_undiag"] - 
      transition_dt[, "f3_f4"]*oldPop[, "f3_undiag"] - 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_undiag"] - 
      tau_ab_dt[, "f3", t]*oldPop[, "f3_undiag"] - 
      tau_poct_dt[, "f3", t]*oldPop[, "f3_undiag"] + 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf3[,]
      reinfP[, t]*foi_dt[, ]*II*oldPop[, "f3_cured"]
    
    
    ## f3 testing, ab+
    newPop[, "f3_diag_ab"] <- oldPop[,"f3_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_diag_ab"]))  - 
      death[,"f3_diag_ab"]  + 
      transition_dt[, "f2_f3"]*oldPop[, "f2_diag_ab"] - 
      transition_dt[, "f3_f4"]*oldPop[, "f3_diag_ab"] - 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_diag_ab"] + 
      tau_ab_dt[, "f3", t]*oldPop[, "f3_undiag"] - 
      tau_ag_dt[, "f3", t]*oldPop[, "f3_diag_ab"]    
    
    
    ## f3 testing, ag+/RNA
    newPop[, "f3_diag_RNA"] <- oldPop[,"f3_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_diag_RNA"])) - 
      death[,"f3_diag_RNA"] +
      transition_dt[, "f2_f3"]*oldPop[, "f2_diag_RNA"] - 
      transition_dt[, "f3_f4"]*oldPop[, "f3_diag_RNA"] - 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_diag_RNA"] + 
      tau_ag_dt[, "f3", t]*oldPop[, "f3_diag_ab"] + 
      tau_poct_dt[, "f3", t]*oldPop[, "f3_undiag"] - 
      eta_dt[, "f3", t]*oldPop[,"f3_diag_RNA"]    
    
    
    ## f3 treat
    newPop[, "f3_treat"] <- oldPop[,"f3_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_treat"])) - 
      death[,"f3_treat"]  + 
      transition_dt[, "f2_f3"]*oldPop[, "f2_treat"] - 
      transition_dt[, "f3_f4"]*oldPop[, "f3_treat"] - 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_treat"] + 
      eta_dt[, "f3", t]*oldPop[,"f3_diag_RNA"] - 
      lota_dt[, "f3", t]*(1-cure_dt[, "f3", t])*oldPop[,"f3_treat"] - 
      cure_dt[, "f3", t]*oldPop[,"f3_treat"] + 
      rho_dt[, "f3", t]*oldPop[,"f3_treat_f"]
    
    ## f3 treat failed
    newPop[, "f3_treat_f"] <- oldPop[,"f3_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_treat_f"])) - 
      death[,"f3_treat_f"] + 
      transition_dt[, "f2_f3"]*oldPop[, "f2_treat_f"] - 
      transition_dt[, "f3_f4"]*oldPop[, "f3_treat_f"] - 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_treat_f"] + 
      lota_dt[, "f3", t]*(1-cure_dt[, "f3", t])*oldPop[,"f3_treat"] - 
      rho_dt[, "f3", t]*oldPop[,"f3_treat_f"]
    
    ## f3 cured
    newPop[, "f3_cured"] <- oldPop[,"f3_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_cured"])) - 
      death[,"f3_cured"] + 
      cure_dt[, "f3", t]*oldPop[, "f3_treat"]  -
      fibprog_dt[, "f3_cured_f4_cured"]*oldPop[, "f3_cured"] - 
      fibprog_dt[, "f3_cured_hcc_cured"]*oldPop[, "f3_cured"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf3[,]
      reinfP[, t]*foi_dt[, ]*II*oldPop[, "f3_cured"]
    
    ##### f4 ##### 
    ## f4 infection
    newPop[, "f4_undiag"] <- oldPop[,"f4_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_undiag"])) - 
      death[,"f4_undiag"] + 
      
      transition_dt[, "f3_f4"]*oldPop[, "f3_undiag"] - 
      transition_dt[, "f4_dc"]*oldPop[, "f4_undiag"] - 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_undiag"] - 
      tau_ab_dt[, "f4", t]*oldPop[, "f4_undiag"] - 
      tau_poct_dt[, "f4", t]*oldPop[, "f4_undiag"] +
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf4[,]
      0*reinfP[, t]*foi_dt[, ]*II*oldPop[, "f4_cured"]
    
    ## f4 testing, ab+
    newPop[, "f4_diag_ab"] <- oldPop[,"f4_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_diag_ab"])) - 
      death[,"f4_diag_ab"] + 
      transition_dt[, "f3_f4"]*oldPop[, "f3_diag_ab"] - 
      transition_dt[, "f4_dc"]*oldPop[, "f4_diag_ab"] - 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_diag_ab"] + 
      tau_ab_dt[, "f4", t]*oldPop[, "f4_undiag"] - 
      tau_ag_dt[, "f4", t]*oldPop[, "f4_diag_ab"]    
    
    
    ## f4 testing, ag+/RNA
    newPop[, "f4_diag_RNA"] <- oldPop[,"f4_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_diag_RNA"])) - 
      death[,"f4_diag_RNA"] + 
      transition_dt[, "f3_f4"]*oldPop[, "f3_diag_RNA"] - 
      transition_dt[, "f4_dc"]*oldPop[, "f4_diag_RNA"] - 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_diag_RNA"] + 
      tau_ag_dt[, "f4", t]*oldPop[, "f4_diag_ab"] + 
      tau_poct_dt[, "f4", t]*oldPop[, "f4_undiag"] - 
      eta_dt[, "f4", t]*oldPop[,"f4_diag_RNA"]    
    
    
    ## f4 treat
    newPop[, "f4_treat"] <- oldPop[,"f4_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_treat"])) - 
      death[,"f4_treat"]  + 
      transition_dt[, "f3_f4"]*oldPop[, "f3_treat"] - 
      transition_dt[, "f4_dc"]*oldPop[, "f4_treat"] - 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_treat"] + 
      eta_dt[, "f4", t]*oldPop[,"f4_diag_RNA"] - 
      lota_dt[, "f4", t]*(1-cure_dt[, "f4", t])*oldPop[,"f4_treat"] - 
      cure_dt[, "f4", t]*oldPop[,"f4_treat"] + 
      rho_dt[, "f4", t]*oldPop[,"f4_treat_f"]
    
    ## f4 treat failed
    newPop[, "f4_treat_f"] <- oldPop[,"f4_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_treat_f"])) - 
      death[,"f4_treat_f"]  + 
      transition_dt[, "f3_f4"]*oldPop[, "f3_treat_f"] - 
      transition_dt[, "f4_dc"]*oldPop[, "f4_treat_f"] - 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_treat_f"] + 
      lota_dt[, "f4", t]*(1-cure_dt[, "f4", t])*oldPop[,"f4_treat"] - 
      rho_dt[, "f4", t]*oldPop[,"f4_treat_f"]
    
    ## f4 cured
    newPop[, "f4_cured"] <- oldPop[,"f4_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_cured"])) - 
      death[,"f4_cured"] + 
      cure_dt[, "f4", t]*oldPop[, "f4_treat"]  + 
      fibprog_dt[, "f3_cured_f4_cured"]*oldPop[, "f3_cured"] - 
      fibprog_dt[, "f4_cured_dc_cured"]*oldPop[, "f4_cured"] - 
      fibprog_dt[, "f4_cured_hcc_cured"]*oldPop[, "f4_cured"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf4[,]
      0*reinfP[, t]*foi_dt[, ]*II*oldPop[, "f4_cured"]
    
    
    
    
    ##### decompensated cirrhosis  (dc) #####
    ## DAA immediately (w/wo liver treatment)
    ## dc infection
    ## dc testing, ab+
    newPop[, "dc_undiag"] <- oldPop[,"dc_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_undiag"])) - 
      death[,"dc_undiag"]  - 
      death_hcv[,"dc_undiag"]  + 
      transition_dt[, "f4_dc"]*oldPop[, "f4_undiag"] - 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_undiag"] - 
      transition_dt[, "dc_lt"]*oldPop[, "dc_undiag"] - 
      tau_ab_dt[, "dc", t]*oldPop[, "dc_undiag"] - 
      tau_poct_dt[, "dc", t]*oldPop[, "dc_undiag"] +
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcdc[,]
      0*reinfP[, t]*foi_dt[, ]*II*oldPop[, "dc_cured"]
    
    newPop[, "dc_diag_ab"] <- oldPop[,"dc_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_diag_ab"])) - 
      death[,"dc_diag_ab"] - 
      death_hcv[,"dc_diag_ab"] + 
      transition_dt[, "f4_dc"]*oldPop[, "f4_diag_ab"] - 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_diag_ab"] - 
      transition_dt[, "dc_lt"]*oldPop[, "dc_diag_ab"] + 
      tau_ab_dt[, "dc", t]*oldPop[, "dc_undiag"] - 
      tau_ag_dt[, "dc", t]*oldPop[, "dc_diag_ab"]    
    
    
    ## dc testing, ag+/RNA
    newPop[, "dc_diag_RNA"] <- oldPop[,"dc_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_diag_RNA"])) - 
      death[,"dc_diag_RNA"] - 
      death_hcv[,"dc_diag_RNA"] + 
      transition_dt[, "f4_dc"]*oldPop[, "f4_diag_RNA"] - 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_diag_RNA"] - 
      transition_dt[, "dc_lt"]*oldPop[, "dc_diag_RNA"] + 
      tau_ag_dt[, "dc", t]*oldPop[, "dc_diag_ab"] + 
      tau_poct_dt[, "dc", t]*oldPop[, "dc_undiag"] - 
      eta_dt[, "dc", t]*oldPop[,"dc_diag_RNA"]    
    
    
    ## dc treat
    newPop[, "dc_treat"] <- oldPop[,"dc_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_treat"])) - 
      death[,"dc_treat"]  - 
      death_hcv[,"dc_treat"] + 
      transition_dt[, "f4_dc"]*oldPop[, "f4_treat"] - 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_treat"] - 
      transition_dt[, "dc_lt"]*oldPop[, "dc_treat"] + 
      eta_dt[, "dc", t]*oldPop[,"dc_diag_RNA"] - 
      lota_dt[, "dc", t]*(1-cure_dt[, "dc", t])*oldPop[,"dc_treat"] - 
      cure_dt[, "dc", t]*oldPop[,"dc_treat"] + 
      rho_dt[, "dc", t]*oldPop[,"dc_treat_f"]
    
    ## dc treat failed
    newPop[, "dc_treat_f"] <- oldPop[,"dc_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_treat_f"])) - 
      death[,"dc_treat_f"]  - 
      death_hcv[,"dc_treat_f"] + 
      transition_dt[, "f4_dc"]*oldPop[, "f4_treat_f"] - 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_treat_f"] - 
      transition_dt[, "dc_lt"]*oldPop[, "dc_treat_f"] + 
      lota_dt[, "dc", t]*(1-cure_dt[, "dc", t])*oldPop[,"dc_treat"] - 
      rho_dt[, "dc", t]*oldPop[,"dc_treat_f"]
    
    ## dc cured
    newPop[, "dc_cured"] <- oldPop[,"dc_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_cured"])) - 
      death[,"dc_cured"] - 
      death_hcv[,"dc_cured"] + 
      cure_dt[, "dc", t]*oldPop[, "dc_treat"]  - 
      fibprog_dt[, "dc_cured_lt_cured"]*oldPop[, "dc_cured"] - 
      fibprog_dt[, "dc_cured_hcc_cured"]*oldPop[, "dc_cured"] + 
      fibprog_dt[, "f4_cured_dc_cured"]*oldPop[, "f4_cured"] -
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcdc[,]
      0*reinfP[, t]*foi_dt[, ]*II*oldPop[, "dc_cured"]
    
    
    
    
    ##### hep carcinoma (HCC) #####
    ## guideline recommend treat HCC first then giving DAA 
    
    ## HCC infection
    newPop[, "hcc_undiag"] <- oldPop[,"hcc_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_undiag"])) - 
      death[,"hcc_undiag"] - 
      death_hcv[,"hcc_undiag"]  + 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_undiag"] + 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_undiag"] + 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_undiag"] - 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_undiag"] - 
      tau_ab_dt[, "hcc", t]*oldPop[, "hcc_undiag"] -  
      tau_poct_dt[, "hcc", t]*oldPop[, "hcc_undiag"] + 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pchcc[,] 
      0*reinfP[, t]*foi_dt[, ]*II*oldPop[, "hcc_cured"]
    
    
    ## hcc testing, ab+
    newPop[, "hcc_diag_ab"] <- oldPop[,"hcc_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_diag_ab"])) - 
      death[,"hcc_diag_ab"] - 
      death_hcv[,"hcc_diag_ab"] + 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_diag_ab"] + 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_diag_ab"] + 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_diag_ab"] - 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_diag_ab"] + 
      tau_ab_dt[, "hcc", t]*oldPop[, "hcc_undiag"] - 
      tau_ag_dt[, "hcc", t]*oldPop[, "hcc_diag_ab"]    
    
    
    ## hcc testing, ag+/RNA
    newPop[, "hcc_diag_RNA"] <- oldPop[,"hcc_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_diag_RNA"])) - 
      death[,"hcc_diag_RNA"] - 
      death_hcv[,"hcc_diag_RNA"] + 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_diag_RNA"] +
      transition_dt[, "dc_hcc"]*oldPop[, "dc_diag_RNA"] + 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_diag_RNA"] - 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_diag_RNA"] + 
      tau_ag_dt[, "hcc", t]*oldPop[, "hcc_diag_ab"] +
      tau_poct_dt[, "hcc", t]*oldPop[, "hcc_undiag"] - 
      eta_dt[, "hcc", t]*oldPop[,"hcc_diag_RNA"]
    
    
    ## hcc treat
    newPop[, "hcc_treat"] <- oldPop[,"hcc_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_treat"])) - 
      death[,"hcc_treat"] - 
      death_hcv[,"hcc_treat"] + 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_treat"] + 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_treat"] + 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_treat"] - 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_treat"] + 
      eta_dt[, "hcc", t]*oldPop[,"hcc_diag_RNA"] - 
      lota_dt[, "hcc", t]*(1-cure_dt[, "hcc", t])*oldPop[,"hcc_treat"] -  
      cure_dt[, "hcc", t]*oldPop[,"hcc_treat"] +
      rho_dt[, "hcc", t]*oldPop[,"hcc_treat_f"]
    
    ## hcc treat failed
    newPop[, "hcc_treat_f"] <- oldPop[,"hcc_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_treat_f"])) - 
      death[,"hcc_treat_f"] - 
      death_hcv[,"hcc_treat_f"]  + 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_treat_f"] + 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_treat_f"] + 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_treat_f"] - 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_treat_f"] + 
      lota_dt[, "hcc", t]*(1-cure_dt[, "hcc", t])*oldPop[,"hcc_treat"] - 
      rho_dt[, "hcc", t]*oldPop[,"hcc_treat_f"]
    
    ## hcc cured
    newPop[, "hcc_cured"] <- oldPop[,"hcc_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_cured"])) - 
      death[,"hcc_cured"] - 
      death_hcv[,"hcc_cured"] + 
      cure_dt[, "hcc", t]*oldPop[, "hcc_treat"]  + 
      fibprog_dt[, "f4_cured_hcc_cured"]*oldPop[, "f4_cured"] + 
      fibprog_dt[, "dc_cured_hcc_cured"]*oldPop[, "dc_cured"] +
      fibprog_dt[, "f3_cured_hcc_cured"]*oldPop[, "f3_cured"] - 
      fibprog_dt[, "hcc_cured_lt_cured"]*oldPop[, "hcc_cured"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pchcc[,] 
      0*reinfP[, t]*foi_dt[, ]*II*oldPop[, "hcc_cured"]
    
    
    
    #####  liver failure needs to liver transplant ##### 
    ### LT infection
    newPop[, "lt_undiag"] <- oldPop[,"lt_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_undiag"]))  - 
      death[,"lt_undiag"] - 
      death_hcv[,"lt_undiag"] + 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_undiag"] + 
      transition_dt[, "dc_lt"]*oldPop[, "dc_undiag"]  -
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_undiag"] - 
      tau_ab_dt[, "lt", t]*oldPop[, "lt_undiag"] - 
      tau_poct_dt[, "lt", t]*oldPop[, "lt_undiag"] + 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pclt[,] 
      0*reinfP[, t]*foi_dt[, ]*II*oldPop[, "lt_cured"]
    
    ### LT testing, ab+
    newPop[, "lt_diag_ab"] <- oldPop[,"lt_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_diag_ab"]))- 
      death[,"lt_diag_ab"] - 
      death_hcv[,"lt_diag_ab"] + 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_diag_ab"] + 
      transition_dt[, "dc_lt"]*oldPop[, "dc_diag_ab"]  - 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_diag_ab"] + 
      tau_ab_dt[, "lt", t]*oldPop[, "lt_undiag"] - 
      tau_ag_dt[, "lt", t]*oldPop[, "lt_diag_ab"]
    
    
    ### LT testing, ag+/RNA
    newPop[, "lt_diag_RNA"] <- oldPop[,"lt_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_diag_RNA"])) - 
      death[,"lt_diag_RNA"] - 
      death_hcv[,"lt_diag_RNA"] + 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_diag_RNA"] + 
      transition_dt[, "dc_lt"]*oldPop[, "dc_diag_RNA"]  - 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_diag_RNA"] +
      tau_ag_dt[, "lt", t]*oldPop[, "lt_diag_ab"] + 
      tau_poct_dt[, "lt", t]*oldPop[, "lt_undiag"] - 
      eta_dt[, "lt", t]*oldPop[,"lt_diag_RNA"]
    
    ### LT treat
    newPop[, "lt_treat"] <- oldPop[,"lt_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_treat"]))  -  
      death[,"lt_treat"] - 
      death_hcv[,"lt_treat"] + 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_treat"] + 
      transition_dt[, "dc_lt"]*oldPop[, "dc_treat"] + 
      eta_dt[, "lt", t]*oldPop[,"lt_diag_RNA"] - 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_treat"] - 
      lota_dt[, "lt", t]*(1-cure_dt[, "lt", t])*oldPop[,"lt_treat"] - 
      cure_dt[, "lt", t]*oldPop[,"lt_treat"] + 
      rho_dt[, "lt", t]*oldPop[,"lt_treat_f"]
    
    ### LT treat failed
    newPop[, "lt_treat_f"] <- oldPop[,"lt_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_treat_f"])) - 
      death[,"lt_treat_f"] - 
      death_hcv[,"lt_treat_f"]  + 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_treat_f"] + 
      transition_dt[, "dc_lt"]*oldPop[, "dc_treat_f"]  + 
      lota_dt[, "lt", t]*(1-cure_dt[, "lt", t])*oldPop[,"lt_treat"] - 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_treat_f"] - 
      rho_dt[, "lt", t]*oldPop[,"lt_treat_f"]
    
    ### LT cured
    newPop[, "lt_cured"] <- oldPop[,"lt_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_cured"]))  - 
      death[,"lt_cured"] - 
      death_hcv[,"lt_cured"] + 
      cure_dt[, "lt", t]*oldPop[, "lt_treat"]  + 
      fibprog_dt[, "hcc_cured_lt_cured"]*oldPop[, "hcc_cured"] + 
      fibprog_dt[, "dc_cured_lt_cured"]*oldPop[, "dc_cured"] - 
      fibprog_dt[, "lt_cured_plt_cured"]*oldPop[, "lt_cured"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pclt[,]
      0*reinfP[, t]*foi_dt[, ]*II*oldPop[, "lt_cured"]
    
    ##### PLT #####
    ### PLT infection
    newPop[, "plt_undiag"] <- oldPop[,"plt_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_undiag"])) - 
      death[,"plt_undiag"]  - 
      death_hcv[,"plt_undiag"]  - 
      tau_ab_dt[, "plt", t]*oldPop[, "plt_undiag"] -
      tau_poct_dt[, "plt", t]*oldPop[, "plt_undiag"] + 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_undiag"] +
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcplt[,]
      0*reinfP[, t]*foi_dt[, ]*II*oldPop[, "plt_cured"]
    
    ### PLT testing, ab+
    newPop[, "plt_diag_ab"] <- oldPop[,"plt_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_diag_ab"])) - 
      death[,"plt_diag_ab"] - 
      death_hcv[,"plt_diag_ab"] + 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_diag_ab"] + 
      tau_ab_dt[, "plt", t]*oldPop[, "plt_undiag"] - 
      tau_ag_dt[, "plt", t]*oldPop[, "plt_diag_ab"]
    
    ### PLT testing, ag+/RNA
    newPop[, "plt_diag_RNA"] <- oldPop[,"plt_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_diag_RNA"])) - 
      death[,"plt_diag_RNA"] - 
      death_hcv[,"plt_diag_RNA"]  + 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_diag_RNA"] + 
      tau_ag_dt[, "plt", t]*oldPop[, "plt_diag_ab"] + 
      tau_poct_dt[, "plt", t]*oldPop[, "plt_undiag"] - 
      eta_dt[, "plt", t]*oldPop[,"plt_diag_RNA"]
    
    ### PLT treat
    newPop[, "plt_treat"] <- oldPop[,"plt_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_treat"]))- 
      death[,"plt_treat"] - 
      death_hcv[,"plt_treat"] + 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_treat"] + 
      eta_dt[, "plt", t]*oldPop[,"plt_diag_RNA"] - 
      lota_dt[, "plt", t]*(1-cure_dt[, "plt", t])*oldPop[,"plt_treat"] - 
      cure_dt[, "plt", t]*oldPop[,"plt_treat"] +
      rho_dt[, "plt", t]*oldPop[,"plt_treat_f"]
    
    ### PLT treat failed
    newPop[, "plt_treat_f"] <- oldPop[,"plt_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_treat_f"])) - 
      death[,"plt_treat_f"] - 
      death_hcv[,"plt_treat_f"] + 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_treat_f"] + 
      lota_dt[, "plt", t]*(1-cure_dt[, "plt", t])*oldPop[,"plt_treat"] - 
      rho_dt[, "plt", t]*oldPop[,"plt_treat_f"]
    
    ### PLT cured
    newPop[, "plt_cured"] <- oldPop[,"plt_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_cured"])) - 
      death[,"plt_cured"] - 
      death_hcv[,"plt_cured"] + 
      cure_dt[, "plt", t]*oldPop[, "plt_treat"] + 
      fibprog_dt[, "lt_cured_plt_cured"]*oldPop[, "lt_cured"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcplt[,] 
      0*reinfP[, t]*foi_dt[, ]*II*oldPop[, "plt_cured"]
    
    
    ifelse (newPop[,] < 0, 0, newPop) 
    
    
    
    
    #### uncertainty 
   
    
    #### result  aggregate ####  
    ##### population #####
    allPops[, , t] <- newPop

    ##### test S #####
    
    newS[, t] <- newPop[ , "s"]
    
    
    ##### new infection #####
    newInfections[, t] <-  foi_dt[, ]*II*oldPop[, "s"]+ 
      reinfP[, t]*foi_dt[, ]*II*(oldPop[, "a_cured"] + 
                                        oldPop[, "f0_cured"] + 
                                        oldPop[, "f1_cured"] + 
                                        oldPop[, "f2_cured"] + 
                                        oldPop[, "f3_cured"] + 
                                        0*oldPop[, "f4_cured"] + 
                                        0*oldPop[, "dc_cured"] + 
                                        0*oldPop[, "hcc_cured"] + 
                                        0*oldPop[, "lt_cured"] + 
                                        0*oldPop[, "plt_cured"])
    
    
    newEntry[1,t] <- sum(entry1)
    
    ##### background death ##### 
    newDeath[1, t] <- sum(death[1,])
    newDeath[2, t] <- sum(death[2, ])
    newDeath[3, t] <- sum(death[3,])
    newDeath[4, t] <- sum(death[4,])
    
    newDeathState[1, ,t] <- morb_dt[1,t]*oldPop[1,] 
    newDeathState[2, ,t] <- morb_dt[2,t]*oldPop[2,] 
    newDeathState[3, ,t] <- morb_dt[3,t]*oldPop[3,]
    newDeathState[4, ,t] <- morb_dt[4,t]*oldPop[4,]
    ##### new HCV deaths in this timestep ##### 
    ## question: should post-liver transplant count in the new HCV death??
    # newHCVdeaths[1, t] <- sum(death_hcv[1, ])
    # newHCVdeaths[2, t] <- sum(death_hcv[2, ])
    # newHCVdeaths[3, t] <- sum(death_hcv[3, ])
    # newHCVdeaths[4, t] <- sum(death_hcv[4, ])
    
    HCVdeathState[,"dc_undiag" ,t] <- mordc_dt[, t]*oldPop[, "dc_undiag"]
    HCVdeathState[,"dc_diag_ab" ,t] <- mordc_dt[, t]*oldPop[, "dc_diag_ab"]
    HCVdeathState[,"dc_diag_RNA" ,t] <- mordc_dt[, t]*oldPop[, "dc_diag_RNA"]
    HCVdeathState[,"dc_treat" ,t] <- mordc_dt[, t]*oldPop[, "dc_treat"]
    HCVdeathState[,"dc_treat_f" ,t] <- mordc_dt[, t]*oldPop[, "dc_treat_f"]
    
    HCVdeathState[,"hcc_undiag" ,t] <- mordc_dt[, t]*oldPop[, "hcc_undiag"]
    HCVdeathState[,"hcc_diag_ab" ,t] <- mordc_dt[, t]*oldPop[, "hcc_diag_ab"]
    HCVdeathState[,"hcc_diag_RNA" ,t] <- mordc_dt[, t]*oldPop[, "hcc_diag_RNA"]
    HCVdeathState[,"hcc_treat" ,t] <- mordc_dt[, t]*oldPop[, "hcc_treat"]
    HCVdeathState[,"hcc_treat_f" ,t] <- mordc_dt[, t]*oldPop[, "hcc_treat_f"]
    
    HCVdeathState[,"lt_undiag" ,t] <- mordc_dt[, t]*oldPop[, "lt_undiag"]
    HCVdeathState[,"lt_diag_ab" ,t] <- mordc_dt[, t]*oldPop[, "lt_diag_ab"]
    HCVdeathState[,"lt_diag_RNA" ,t] <- mordc_dt[, t]*oldPop[, "lt_diag_RNA"]
    HCVdeathState[,"lt_treat" ,t] <- mordc_dt[, t]*oldPop[, "lt_treat"]
    HCVdeathState[,"lt_treat_f" ,t] <- mordc_dt[, t]*oldPop[, "lt_treat_f"]
    
    HCVdeathState[,"plt_undiag" ,t] <- mordc_dt[, t]*oldPop[, "plt_undiag"]
    HCVdeathState[,"plt_diag_ab" ,t] <- mordc_dt[, t]*oldPop[, "plt_diag_ab"]
    HCVdeathState[,"plt_diag_RNA" ,t] <- mordc_dt[, t]*oldPop[, "plt_diag_RNA"]
    HCVdeathState[,"plt_treat" ,t] <- mordc_dt[, t]*oldPop[, "plt_treat"]
    HCVdeathState[,"plt_treat_f" ,t] <- mordc_dt[, t]*oldPop[, "plt_treat_f"]
    
    newHCVdeaths[1, t] <- sum(mordc_dt[1, t]*oldPop[1, "dc_undiag"],
                              mordc_dt[1, t]*oldPop[1, "dc_diag_ab"],
                              mordc_dt[1, t]*oldPop[1, "dc_diag_RNA"],
                              mordc_dt[1,t]*oldPop[1, "dc_treat"],
                              mordc_dt[1,t]*oldPop[1, "dc_treat_f"],
                              mordcCure_dt[1,t]*oldPop[1, "dc_cured"],                      
                              morhcc_dt[1, t]*oldPop[1, "hcc_undiag"],                        
                              morhcc_dt[1, t]*oldPop[1, "hcc_diag_ab"],                       
                              morhcc_dt[1,t]*oldPop[1, "hcc_diag_RNA"],
                              morhcc_dt[1,t]*oldPop[1, "hcc_treat"],
                              morhcc_dt[1,t]*oldPop[1, "hcc_treat_f"],
                              morhccCure_dt[1,t]*oldPop[1, "hcc_cured"],
                              morlt_dt[1, t]*oldPop[1, "lt_undiag"],
                              morlt_dt[1, t]*oldPop[1, "lt_diag_ab"],
                              morlt_dt[1,t]*oldPop[1, "lt_diag_RNA"],
                              morlt_dt[1,t]*oldPop[1, "lt_treat"],
                              morlt_dt[1,t]*oldPop[1, "lt_treat_f"],
                              morlt_dt[1,t]*oldPop[1, "lt_cured"],
                              morplt_dt[1, t]*oldPop[1, "plt_undiag"],
                              morplt_dt[1, t]*oldPop[1, "plt_diag_ab"],
                              morplt_dt[1,t]*oldPop[1, "plt_diag_RNA"],
                              morplt_dt[1,t]*oldPop[1, "plt_treat"],
                              morplt_dt[1,t]*oldPop[1, "plt_treat_f"],
                              morplt_dt[1,t]*oldPop[1, "plt_cured"])
    
    newHCVdeaths[2, t] <- sum(mordc_dt[2, t]*oldPop[2, "dc_undiag"],
                              mordc_dt[2, t]*oldPop[2, "dc_diag_ab"],
                              mordc_dt[2, t]*oldPop[2, "dc_diag_RNA"],
                              mordc_dt[2,t]*oldPop[2, "dc_treat"],
                              mordc_dt[2,t]*oldPop[2, "dc_treat_f"],
                              mordcCure_dt[2,t]*oldPop[2, "dc_cured"],                      
                              morhcc_dt[2, t]*oldPop[2, "hcc_undiag"],                        
                              morhcc_dt[2, t]*oldPop[2, "hcc_diag_ab"],                       
                              morhcc_dt[2,t]*oldPop[2, "hcc_diag_RNA"],
                              morhcc_dt[2,t]*oldPop[2, "hcc_treat"],
                              morhcc_dt[2,t]*oldPop[2, "hcc_treat_f"],
                              morhccCure_dt[2,t]*oldPop[2, "hcc_cured"],
                              morlt_dt[2, t]*oldPop[2, "lt_undiag"],
                              morlt_dt[2, t]*oldPop[2, "lt_diag_ab"],
                              morlt_dt[2,t]*oldPop[2, "lt_diag_RNA"],
                              morlt_dt[2,t]*oldPop[2, "lt_treat"],
                              morlt_dt[2,t]*oldPop[2, "lt_treat_f"],
                              morlt_dt[2,t]*oldPop[2, "lt_cured"],
                              morplt_dt[2, t]*oldPop[2, "plt_undiag"],
                              morplt_dt[2, t]*oldPop[2, "plt_diag_ab"],
                              morplt_dt[2,t]*oldPop[2, "plt_diag_RNA"],
                              morplt_dt[2,t]*oldPop[2, "plt_treat"],
                              morplt_dt[2,t]*oldPop[2, "plt_treat_f"],
                              morplt_dt[2,t]*oldPop[2, "plt_cured"])
    
    newHCVdeaths[3, t] <- sum(mordc_dt[3, t]*oldPop[3, "dc_undiag"],
                              mordc_dt[3, t]*oldPop[3, "dc_diag_ab"],
                              mordc_dt[3, t]*oldPop[3, "dc_diag_RNA"],
                              mordc_dt[3,t]*oldPop[3, "dc_treat"],
                              mordc_dt[3,t]*oldPop[3, "dc_treat_f"],
                              mordcCure_dt[3,t]*oldPop[3, "dc_cured"],                      
                              morhcc_dt[3, t]*oldPop[3, "hcc_undiag"],                        
                              morhcc_dt[3, t]*oldPop[3, "hcc_diag_ab"],                       
                              morhcc_dt[3,t]*oldPop[3, "hcc_diag_RNA"],
                              morhcc_dt[3,t]*oldPop[3, "hcc_treat"],
                              morhcc_dt[3,t]*oldPop[3, "hcc_treat_f"],
                              morhccCure_dt[3,t]*oldPop[3, "hcc_cured"],
                              morlt_dt[3, t]*oldPop[3, "lt_undiag"],
                              morlt_dt[3, t]*oldPop[3, "lt_diag_ab"],
                              morlt_dt[3,t]*oldPop[3, "lt_diag_RNA"],
                              morlt_dt[3,t]*oldPop[3, "lt_treat"],
                              morlt_dt[3,t]*oldPop[3, "lt_treat_f"],
                              morlt_dt[3,t]*oldPop[3, "lt_cured"],
                              morplt_dt[3, t]*oldPop[3, "plt_undiag"],
                              morplt_dt[3, t]*oldPop[3, "plt_diag_ab"],
                              morplt_dt[3,t]*oldPop[3, "plt_diag_RNA"],
                              morplt_dt[3,t]*oldPop[3, "plt_treat"],
                              morplt_dt[3,t]*oldPop[3, "plt_treat_f"],
                              morplt_dt[3,t]*oldPop[3, "plt_cured"])
    
    newHCVdeaths[4, t] <- sum(mordc_dt[4, t]*oldPop[4, "dc_undiag"],
                              mordc_dt[4, t]*oldPop[4, "dc_diag_ab"],
                              mordc_dt[4, t]*oldPop[4, "dc_diag_RNA"],
                              mordc_dt[4,t]*oldPop[4, "dc_treat"],
                              mordc_dt[4,t]*oldPop[4, "dc_treat_f"],
                              mordcCure_dt[4,t]*oldPop[4, "dc_cured"],                      
                              morhcc_dt[4, t]*oldPop[4, "hcc_undiag"],                        
                              morhcc_dt[4, t]*oldPop[4, "hcc_diag_ab"],                       
                              morhcc_dt[4,t]*oldPop[4, "hcc_diag_RNA"],
                              morhcc_dt[4,t]*oldPop[4, "hcc_treat"],
                              morhcc_dt[4,t]*oldPop[4, "hcc_treat_f"],
                              morhccCure_dt[4,t]*oldPop[4, "hcc_cured"],
                              morlt_dt[4, t]*oldPop[4, "lt_undiag"],
                              morlt_dt[4, t]*oldPop[4, "lt_diag_ab"],
                              morlt_dt[4,t]*oldPop[4, "lt_diag_RNA"],
                              morlt_dt[4,t]*oldPop[4, "lt_treat"],
                              morlt_dt[4,t]*oldPop[4, "lt_treat_f"],
                              morlt_dt[4,t]*oldPop[4, "lt_cured"],
                              morplt_dt[4, t]*oldPop[4, "plt_undiag"],
                              morplt_dt[4, t]*oldPop[4, "plt_diag_ab"],
                              morplt_dt[4,t]*oldPop[4, "plt_diag_RNA"],
                              morplt_dt[4,t]*oldPop[4, "plt_treat"],
                              morplt_dt[4,t]*oldPop[4, "plt_treat_f"],
                              morplt_dt[4,t]*oldPop[4, "plt_cured"])
    
    ##### new treat #####
    
    newTreatment[1, t] <- sum(eta_dt[1, "a", t]*oldPop[1,"a_diag_RNA"],
                              eta_dt[1, "f0", t]*oldPop[1,"f0_diag_RNA"],
                              eta_dt[1, "f1", t]*oldPop[1,"f1_diag_RNA"],
                              eta_dt[1, "f2", t]*oldPop[1,"f2_diag_RNA"],
                              eta_dt[1, "f3", t]*oldPop[1,"f3_diag_RNA"],
                              eta_dt[1, "f4", t]*oldPop[1,"f4_diag_RNA"],
                              eta_dt[1, "dc", t]*oldPop[1,"dc_diag_RNA"],
                              eta_dt[1, "hcc", t]*oldPop[1,"hcc_diag_RNA"],
                              eta_dt[1, "lt", t]*oldPop[1,"lt_diag_RNA"],
                              eta_dt[1, "plt", t]*oldPop[1,"plt_diag_RNA"]) 
    
    
    newTreatment[2, t] <- sum(eta_dt[2, "a", t]*oldPop[2,"a_diag_RNA"],
                              eta_dt[2, "f0", t]*oldPop[2,"f0_diag_RNA"],
                              eta_dt[2, "f1", t]*oldPop[2,"f1_diag_RNA"],
                              eta_dt[2, "f2", t]*oldPop[2,"f2_diag_RNA"],
                              eta_dt[2, "f3", t]*oldPop[2,"f3_diag_RNA"],
                              eta_dt[2, "f4", t]*oldPop[2,"f4_diag_RNA"],
                              eta_dt[2, "dc", t]*oldPop[2,"dc_diag_RNA"],
                              eta_dt[2, "hcc", t]*oldPop[2,"hcc_diag_RNA"],
                              eta_dt[2, "lt", t]*oldPop[2,"lt_diag_RNA"],
                              eta_dt[2, "plt", t]*oldPop[2,"plt_diag_RNA"])
    
    newTreatment[3, t] <- sum(eta_dt[3, "a", t]*oldPop[3,"a_diag_RNA"],
                              eta_dt[3, "f0", t]*oldPop[3,"f0_diag_RNA"],
                              eta_dt[3, "f1", t]*oldPop[3,"f1_diag_RNA"],
                              eta_dt[3, "f2", t]*oldPop[3,"f2_diag_RNA"],
                              eta_dt[3, "f3", t]*oldPop[3,"f3_diag_RNA"],
                              eta_dt[3, "f4", t]*oldPop[3,"f4_diag_RNA"],
                              eta_dt[3, "dc", t]*oldPop[3,"dc_diag_RNA"],
                              eta_dt[3, "hcc", t]*oldPop[3,"hcc_diag_RNA"],
                              eta_dt[3, "lt", t]*oldPop[3,"lt_diag_RNA"],
                              eta_dt[3, "plt", t]*oldPop[3,"plt_diag_RNA"])
    
    
    newTreatment[4, t] <- sum(eta_dt[4, "a", t]*oldPop[4,"a_diag_RNA"],
                              eta_dt[4, "f0", t]*oldPop[4,"f0_diag_RNA"],
                              eta_dt[4, "f1", t]*oldPop[4,"f1_diag_RNA"],
                              eta_dt[4, "f2", t]*oldPop[4,"f2_diag_RNA"],
                              eta_dt[4, "f3", t]*oldPop[4,"f3_diag_RNA"],
                              eta_dt[4, "f4", t]*oldPop[4,"f4_diag_RNA"],
                              eta_dt[4, "dc", t]*oldPop[4,"dc_diag_RNA"],
                              eta_dt[4, "hcc", t]*oldPop[4,"hcc_diag_RNA"],
                              eta_dt[4, "lt", t]*oldPop[4,"lt_diag_RNA"],
                              eta_dt[4, "plt", t]*oldPop[4,"plt_diag_RNA"])
    
    ##### new retreat #####  
    
    newRetreat[1 ,t] <- sum(rho_dt[1, "a", t]*oldPop[1,"a_treat_f"],
                            rho_dt[1, "f0", t]*oldPop[1,"f0_treat_f"],
                            rho_dt[1, "f1", t]*oldPop[1,"f1_treat_f"],
                            rho_dt[1, "f2", t]*oldPop[1,"f2_treat_f"],
                            rho_dt[1, "f3", t]*oldPop[1,"f3_treat_f"],
                            rho_dt[1, "f4", t]*oldPop[1,"f4_treat_f"],
                            rho_dt[1, "dc", t]*oldPop[1,"dc_treat_f"],
                            rho_dt[1, "hcc", t]*oldPop[1,"hcc_treat_f"],
                            rho_dt[1, "lt", t]*oldPop[1,"lt_treat_f"],
                            rho_dt[1, "plt", t]*oldPop[1,"plt_treat_f"])
    
    newRetreat[2, t] <- sum(rho_dt[2, "a", t]*oldPop[2,"a_treat_f"],
                            rho_dt[2, "f0", t]*oldPop[2,"f0_treat_f"],
                            rho_dt[2, "f1", t]*oldPop[2,"f1_treat_f"],
                            rho_dt[2, "f2", t]*oldPop[2,"f2_treat_f"],
                            rho_dt[2, "f3", t]*oldPop[2,"f3_treat_f"],
                            rho_dt[2, "f4", t]*oldPop[2,"f4_treat_f"],
                            rho_dt[2, "dc", t]*oldPop[2,"dc_treat_f"],
                            rho_dt[2, "hcc", t]*oldPop[2,"hcc_treat_f"],
                            rho_dt[2, "lt", t]*oldPop[2,"lt_treat_f"],
                            rho_dt[2, "plt", t]*oldPop[2,"plt_treat_f"])
    
    newRetreat[3, t] <- sum(rho_dt[3, "a", t]*oldPop[3,"a_treat_f"],
                            rho_dt[3, "f0", t]*oldPop[3,"f0_treat_f"],
                            rho_dt[3, "f1", t]*oldPop[3,"f1_treat_f"],
                            rho_dt[3, "f2", t]*oldPop[3,"f2_treat_f"],
                            rho_dt[3, "f3", t]*oldPop[3,"f3_treat_f"],
                            rho_dt[3, "f4", t]*oldPop[3,"f4_treat_f"],
                            rho_dt[3, "dc", t]*oldPop[3,"dc_treat_f"],
                            rho_dt[3, "hcc", t]*oldPop[3,"hcc_treat_f"],
                            rho_dt[3, "lt", t]*oldPop[3,"lt_treat_f"],
                            rho_dt[3, "plt", t]*oldPop[3,"plt_treat_f"])
    
    newRetreat[4, t] <- sum(rho_dt[4, "a", t]*oldPop[4,"a_treat_f"],
                            rho_dt[4, "f0", t]*oldPop[4,"f0_treat_f"],
                            rho_dt[4, "f1", t]*oldPop[4,"f1_treat_f"],
                            rho_dt[4, "f2", t]*oldPop[4,"f2_treat_f"],
                            rho_dt[4, "f3", t]*oldPop[4,"f3_treat_f"],
                            rho_dt[4, "f4", t]*oldPop[4,"f4_treat_f"],
                            rho_dt[4, "dc", t]*oldPop[4,"dc_treat_f"],
                            rho_dt[4, "hcc", t]*oldPop[4,"hcc_treat_f"],
                            rho_dt[4, "lt", t]*oldPop[4,"lt_treat_f"],
                            rho_dt[4, "plt", t]*oldPop[4,"plt_treat_f"])
    
    ##### antibody test ##### 
    
    newTestingAb[1 , t] <- sum( tau_ab_dt[1,"a", t]*(1-spc1_dt[1, t])*oldPop[1,"a_undiag"],
                                tau_ab_dt[1, "f0", t]*oldPop[1,"f0_undiag"],
                                tau_ab_dt[1, "f1", t]*oldPop[1,"f1_undiag"],
                                tau_ab_dt[1, "f2", t]*oldPop[1,"f2_undiag"],
                                tau_ab_dt[1, "f3", t]*oldPop[1,"f3_undiag"],
                                tau_ab_dt[1, "f4", t]*oldPop[1,"f4_undiag"],
                                tau_ab_dt[1, "dc", t]*oldPop[1,"dc_undiag"],
                                tau_ab_dt[1, "hcc", t]*oldPop[1,"hcc_undiag"],
                                tau_ab_dt[1, "lt", t]*oldPop[1,"lt_undiag"],
                                tau_ab_dt[1, "plt", t]*oldPop[1,"plt_undiag"])
    
    newTestingAb[2, t] <- sum( tau_ab_dt[2,"a", t]*(1-spc1_dt[2, t])*oldPop[2,"a_undiag"],
                               tau_ab_dt[2, "f0", t]*oldPop[2,"f0_undiag"],
                               tau_ab_dt[2, "f1", t]*oldPop[2,"f1_undiag"],
                               tau_ab_dt[2, "f2", t]*oldPop[2,"f2_undiag"],
                               tau_ab_dt[2, "f3", t]*oldPop[2,"f3_undiag"],
                               tau_ab_dt[2, "f4", t]*oldPop[2,"f4_undiag"],
                               tau_ab_dt[2, "dc", t]*oldPop[2,"dc_undiag"],
                               tau_ab_dt[2, "hcc", t]*oldPop[2,"hcc_undiag"],
                               tau_ab_dt[2, "lt", t]*oldPop[2,"lt_undiag"],
                               tau_ab_dt[2, "plt", t]*oldPop[2,"plt_undiag"])
    
    newTestingAb[3, t] <- sum( tau_ab_dt[3,"a", t]*(1-spc1_dt[3, t])*oldPop[3,"a_undiag"],
                               tau_ab_dt[3, "f0", t]*oldPop[3,"f0_undiag"],
                               tau_ab_dt[3, "f1", t]*oldPop[3,"f1_undiag"],
                               tau_ab_dt[3, "f2", t]*oldPop[3,"f2_undiag"],
                               tau_ab_dt[3, "f3", t]*oldPop[3,"f3_undiag"],
                               tau_ab_dt[3, "f4", t]*oldPop[3,"f4_undiag"],
                               tau_ab_dt[3, "dc", t]*oldPop[3,"dc_undiag"],
                               tau_ab_dt[3, "hcc", t]*oldPop[3,"hcc_undiag"],
                               tau_ab_dt[3, "lt", t]*oldPop[3,"lt_undiag"],
                               tau_ab_dt[3, "plt", t]*oldPop[3,"plt_undiag"])
    
    newTestingAb[4, t] <- sum( tau_ab_dt[4,"a", t]*(1-spc1_dt[4, t])*oldPop[4,"a_undiag"],
                               tau_ab_dt[4, "f0", t]*oldPop[4,"f0_undiag"],
                               tau_ab_dt[4, "f1", t]*oldPop[4,"f1_undiag"],
                               tau_ab_dt[4, "f2", t]*oldPop[4,"f2_undiag"],
                               tau_ab_dt[4, "f3", t]*oldPop[4,"f3_undiag"],
                               tau_ab_dt[4, "f4", t]*oldPop[4,"f4_undiag"],
                               tau_ab_dt[4, "dc", t]*oldPop[4,"dc_undiag"],
                               tau_ab_dt[4, "hcc", t]*oldPop[4,"hcc_undiag"],
                               tau_ab_dt[4, "lt", t]*oldPop[4,"lt_undiag"],
                               tau_ab_dt[4, "plt", t]*oldPop[4,"plt_undiag"])
    
    ##### antigen test ##### 
    
    newTestingAg[1 ,t] <- sum(tau_ag_dt[1, "a", t]*oldPop[1,"a_diag_ab"],
                              tau_ag_dt[1, "f0", t]*oldPop[1,"f0_diag_ab"],
                              tau_ag_dt[1, "f1", t]*oldPop[1,"f1_diag_ab"],
                              tau_ag_dt[1, "f2", t]*oldPop[1,"f2_diag_ab"],
                              tau_ag_dt[1, "f3", t]*oldPop[1,"f3_diag_ab"],
                              tau_ag_dt[1, "f4", t]*oldPop[1,"f4_diag_ab"],
                              tau_ag_dt[1, "dc", t]*oldPop[1,"dc_diag_ab"],
                              tau_ag_dt[1, "hcc", t]*oldPop[1,"hcc_diag_ab"],
                              tau_ag_dt[1, "lt", t]*oldPop[1,"lt_diag_ab"],
                              tau_ag_dt[1, "plt", t]*oldPop[1,"plt_diag_ab"])
    
    newTestingAg[2, t] <- sum(tau_ag_dt[2, "a", t]*oldPop[2,"a_diag_ab"],
                              tau_ag_dt[2, "f0", t]*oldPop[2,"f0_diag_ab"],
                              tau_ag_dt[2, "f1", t]*oldPop[2,"f1_diag_ab"],
                              tau_ag_dt[2, "f2", t]*oldPop[2,"f2_diag_ab"],
                              tau_ag_dt[2, "f3", t]*oldPop[2,"f3_diag_ab"],
                              tau_ag_dt[2, "f4", t]*oldPop[2,"f4_diag_ab"],
                              tau_ag_dt[2, "dc", t]*oldPop[2,"dc_diag_ab"],
                              tau_ag_dt[2, "hcc", t]*oldPop[2,"hcc_diag_ab"],
                              tau_ag_dt[2, "lt", t]*oldPop[2,"lt_diag_ab"],
                              tau_ag_dt[2, "plt", t]*oldPop[2,"plt_diag_ab"])
    
    newTestingAg[3, t] <- sum(tau_ag_dt[3, "a", t]*oldPop[3,"a_diag_ab"],
                              tau_ag_dt[3, "f0", t]*oldPop[3,"f0_diag_ab"],
                              tau_ag_dt[3, "f1", t]*oldPop[3,"f1_diag_ab"],
                              tau_ag_dt[3, "f2", t]*oldPop[3,"f2_diag_ab"],
                              tau_ag_dt[3, "f3", t]*oldPop[3,"f3_diag_ab"],
                              tau_ag_dt[3, "f4", t]*oldPop[3,"f4_diag_ab"],
                              tau_ag_dt[3, "dc", t]*oldPop[3,"dc_diag_ab"],
                              tau_ag_dt[3, "hcc", t]*oldPop[3,"hcc_diag_ab"],
                              tau_ag_dt[3, "lt", t]*oldPop[3,"lt_diag_ab"],
                              tau_ag_dt[3, "plt", t]*oldPop[3,"plt_diag_ab"])
    
    newTestingAg[4, t] <- sum(tau_ag_dt[4, "a", t]*oldPop[4,"a_diag_ab"],
                              tau_ag_dt[4, "f0", t]*oldPop[4,"f0_diag_ab"],
                              tau_ag_dt[4, "f1", t]*oldPop[4,"f1_diag_ab"],
                              tau_ag_dt[4, "f2", t]*oldPop[4,"f2_diag_ab"],
                              tau_ag_dt[4, "f3", t]*oldPop[4,"f3_diag_ab"],
                              tau_ag_dt[4, "f4", t]*oldPop[4,"f4_diag_ab"],
                              tau_ag_dt[4, "dc", t]*oldPop[4,"dc_diag_ab"],
                              tau_ag_dt[4, "hcc", t]*oldPop[4,"hcc_diag_ab"],
                              tau_ag_dt[4, "lt", t]*oldPop[4,"lt_diag_ab"],
                              tau_ag_dt[4, "plt", t]*oldPop[4,"plt_diag_ab"])
    
    ##### POCT test #####
    
    newTestingPOCT[1, t] <- sum( tau_poct_dt[1, "a", t]*(1-spc1_dt[1, t])*oldPop[1,"a_undiag"],
                                 tau_poct_dt[1, "f0", t]*oldPop[1,"f0_undiag"],
                                 tau_poct_dt[1, "f1", t]*oldPop[1,"f1_undiag"],
                                 tau_poct_dt[1, "f2", t]*oldPop[1,"f2_undiag"],
                                 tau_poct_dt[1, "f3", t]*oldPop[1,"f3_undiag"],
                                 tau_poct_dt[1, "f4", t]*oldPop[1,"f4_undiag"],
                                 tau_poct_dt[1, "dc", t]*oldPop[1,"dc_undiag"],
                                 tau_poct_dt[1, "hcc", t]*oldPop[1,"hcc_undiag"],
                                 tau_poct_dt[1, "lt", t]*oldPop[1,"lt_undiag"],
                                 tau_poct_dt[1, "plt", t]*oldPop[1,"plt_undiag"])
    
    newTestingPOCT[2,t] <- sum( tau_poct_dt[2, "a", t]*(1-spc1_dt[2, t])*oldPop[2,"a_undiag"],
                                tau_poct_dt[2, "f0", t]*oldPop[2,"f0_undiag"],
                                tau_poct_dt[2, "f1", t]*oldPop[2,"f1_undiag"],
                                tau_poct_dt[2, "f2", t]*oldPop[2,"f2_undiag"],
                                tau_poct_dt[2, "f3", t]*oldPop[2,"f3_undiag"],
                                tau_poct_dt[2, "f4", t]*oldPop[2,"f4_undiag"],
                                tau_poct_dt[2, "dc", t]*oldPop[2,"dc_undiag"],
                                tau_poct_dt[2, "hcc", t]*oldPop[2,"hcc_undiag"],
                                tau_poct_dt[2, "lt", t]*oldPop[2,"lt_undiag"],
                                tau_poct_dt[2, "plt", t]*oldPop[2,"plt_undiag"])
    
    newTestingPOCT[3,t] <- sum( tau_poct_dt[3, "a", t]*(1-spc1_dt[3, t])*oldPop[3,"a_undiag"],
                                tau_poct_dt[3, "f0", t]*oldPop[3,"f0_undiag"],
                                tau_poct_dt[3, "f1", t]*oldPop[3,"f1_undiag"],
                                tau_poct_dt[3, "f2", t]*oldPop[3,"f2_undiag"],
                                tau_poct_dt[3, "f3", t]*oldPop[3,"f3_undiag"],
                                tau_poct_dt[3, "f4", t]*oldPop[3,"f4_undiag"],
                                tau_poct_dt[3, "dc", t]*oldPop[3,"dc_undiag"],
                                tau_poct_dt[3, "hcc", t]*oldPop[3,"hcc_undiag"],
                                tau_poct_dt[3, "lt", t]*oldPop[3,"lt_undiag"],
                                tau_poct_dt[3, "plt", t]*oldPop[3,"plt_undiag"])
    
    newTestingPOCT[4,t] <- sum( tau_poct_dt[4, "a", t]*(1-spc1_dt[4, t])*oldPop[4,"a_undiag"],
                                tau_poct_dt[4, "f0", t]*oldPop[4,"f0_undiag"],
                                tau_poct_dt[4, "f1", t]*oldPop[4,"f1_undiag"],
                                tau_poct_dt[4, "f2", t]*oldPop[4,"f2_undiag"],
                                tau_poct_dt[4, "f3", t]*oldPop[4,"f3_undiag"],
                                tau_poct_dt[4, "f4", t]*oldPop[4,"f4_undiag"],
                                tau_poct_dt[4, "dc", t]*oldPop[4,"dc_undiag"],
                                tau_poct_dt[4, "hcc", t]*oldPop[4,"hcc_undiag"],
                                tau_poct_dt[4, "lt", t]*oldPop[4,"lt_undiag"],
                                tau_poct_dt[4, "plt", t]*oldPop[4,"plt_undiag"])
    ##### Cured ##### 
    newCured[1, t] <- sum( cure_dt[1, "a", t]*oldPop[1,"a_treat"],
                           cure_dt[1, "f0", t]*oldPop[1,"f0_treat"],
                           cure_dt[1, "f1", t]*oldPop[1,"f1_treat"],
                           cure_dt[1, "f2", t]*oldPop[1,"f2_treat"],
                           cure_dt[1, "f3", t]*oldPop[1,"f3_treat"],
                           cure_dt[1, "f4", t]*oldPop[1,"f4_treat"],
                           cure_dt[1, "dc", t]*oldPop[1,"dc_treat"],
                           cure_dt[1, "hcc", t]*oldPop[1,"hcc_treat"],
                           cure_dt[1, "lt", t]*oldPop[1,"lt_treat"],
                           cure_dt[1, "plt", t]*oldPop[1,"plt_treat"])
    
    newCured[2, t] <- sum( cure_dt[2, "a", t]*oldPop[2,"a_treat"],
                           cure_dt[2, "f0", t]*oldPop[2,"f0_treat"],
                           cure_dt[2, "f1", t]*oldPop[2,"f1_treat"],
                           cure_dt[2, "f2", t]*oldPop[2,"f2_treat"],
                           cure_dt[2, "f3", t]*oldPop[2,"f3_treat"],
                           cure_dt[2, "f4", t]*oldPop[2,"f4_treat"],
                           cure_dt[2, "dc", t]*oldPop[2,"dc_treat"],
                           cure_dt[2, "hcc", t]*oldPop[2,"hcc_treat"],
                           cure_dt[2, "lt", t]*oldPop[2,"lt_treat"],
                           cure_dt[2, "plt", t]*oldPop[2,"plt_treat"])
    
    newCured[3, t] <- sum( cure_dt[3, "a", t]*oldPop[3,"a_treat"],
                           cure_dt[3, "f0", t]*oldPop[3,"f0_treat"],
                           cure_dt[3, "f1", t]*oldPop[3,"f1_treat"],
                           cure_dt[3, "f2", t]*oldPop[3,"f2_treat"],
                           cure_dt[3, "f3", t]*oldPop[3,"f3_treat"],
                           cure_dt[3, "f4", t]*oldPop[3,"f4_treat"],
                           cure_dt[3, "dc", t]*oldPop[3,"dc_treat"],
                           cure_dt[3, "hcc", t]*oldPop[3,"hcc_treat"],
                           cure_dt[3, "lt", t]*oldPop[3,"lt_treat"],
                           cure_dt[3, "plt", t]*oldPop[3,"plt_treat"])
    
    newCured[4, t] <- sum( cure_dt[4, "a", t]*oldPop[4,"a_treat"],
                           cure_dt[4, "f0", t]*oldPop[4,"f0_treat"],
                           cure_dt[4, "f1", t]*oldPop[4,"f1_treat"],
                           cure_dt[4, "f2", t]*oldPop[4,"f2_treat"],
                           cure_dt[4, "f3", t]*oldPop[4,"f3_treat"],
                           cure_dt[4, "f4", t]*oldPop[4,"f4_treat"],
                           cure_dt[4, "dc", t]*oldPop[4,"dc_treat"],
                           cure_dt[4, "hcc", t]*oldPop[4,"hcc_treat"],
                           cure_dt[4, "lt", t]*oldPop[4,"lt_treat"],
                           cure_dt[4, "plt", t]*oldPop[4,"plt_treat"])
    
    #####reinfection##### 
    
    
    newreinfection[, t] <- reinfP[, t]*foi_dt[, ]*II*(oldPop[, "a_cured"] + 
                                                             oldPop[, "f0_cured"] + 
                                                             oldPop[, "f1_cured"] + 
                                                             oldPop[, "f2_cured"] + 
                                                             oldPop[, "f3_cured"] + 
                                                             0*oldPop[, "f4_cured"] + 
                                                             0*oldPop[, "dc_cured"] + 
                                                             0*oldPop[, "hcc_cured"] + 
                                                             0*oldPop[, "lt_cured"] + 
                                                             0*oldPop[, "plt_cured"])
      
      
      
      
      #reinfP[, t]*foi_dt[,]*I[, ]*(Pca[,]  + Pcf0[,] + 
       #                                                   Pcf1[,] + Pcf2[,] + 
        #                                                  Pcf3[,] + Pcf4[,] + 
         #                                                 Pcdc[,] + Pchcc[,] + 
          #                                                Pclt[,] + Pcplt[,])
    
    # newreinfection[1, t] <- sum(reinfP[1 , t]*foi_dt[1, ]*oldPop[1 ,"a_cured"],
    #                            reinfP[ 1, t]*foi_dt[1, ]*oldPop[ 1,"f0_cured"],
    #                            reinfP[ 1, t]*foi_dt[1, ]*oldPop[ 1,"f1_cured"],
    #                            reinfP[ 1, t]*foi_dt[1, ]*oldPop[1 ,"f2_cured"],
    #                            reinfP[ 1, t]*foi_dt[1, ]*oldPop[ 1,"f3_cured"],
    #                            reinfP[ 1, t]*foi_dt[1, ]*oldPop[1 ,"f4_cured"],
    #                            reinfP[ 1, t]*foi_dt[1, ]*oldPop[1 ,"dc_cured"],
    #                            reinfP[ 1, t]*foi_dt[1, ]*oldPop[1 ,"hcc_cured"],
    #                            reinfP[ 1, t]*foi_dt[1, ]*oldPop[1 ,"lt_cured"],
    #                            reinfP[ 1, t]*foi_dt[1, ]*oldPop[1 ,"plt_cured"])
    
    # newreinfection[2, t] <- sum(reinfP[ 2, t]*foi_dt[2, ]*oldPop[2 ,"a_cured"],
    #                             reinfP[ 2, t]*foi_dt[2, ]*oldPop[2,"f0_cured"],
    #                            reinfP[ 2, t]*foi_dt[2, ]*oldPop[2,"f1_cured"],
    #                             reinfP[ 2, t]*foi_dt[2, ]*oldPop[2,"f2_cured"],
    #                            reinfP[ 2, t]*foi_dt[2, ]*oldPop[2,"f3_cured"],
    #                             reinfP[ 2, t]*foi_dt[2, ]*oldPop[2 ,"f4_cured"],
    #                            reinfP[ 2, t]*foi_dt[2, ]*oldPop[2 ,"dc_cured"],
    #                              reinfP[ 2, t]*foi_dt[2, ]*oldPop[2 ,"hcc_cured"],
    #                             reinfP[ 2, t]*foi_dt[2, ]*oldPop[2 ,"lt_cured"],
    #                              reinfP[ 2, t]*foi_dt[2, ]*oldPop[2 ,"plt_cured"])
    
    # newreinfection[3, t] <- sum(reinfP[ 3, t]*foi_dt[3, ]*oldPop[3 ,"a_cured"],
    #                             reinfP[ 3, t]*foi_dt[3, ]*oldPop[3,"f0_cured"],
    #                             reinfP[ 3, t]*foi_dt[3, ]*oldPop[3,"f1_cured"],
    #                             reinfP[ 3, t]*foi_dt[3, ]*oldPop[3,"f2_cured"],
    #                             reinfP[ 3, t]*foi_dt[3, ]*oldPop[3,"f3_cured"],
    #                             reinfP[ 3, t]*foi_dt[3, ]*oldPop[3 ,"f4_cured"],
    #                             reinfP[ 3, t]*foi_dt[3, ]*oldPop[3 ,"dc_cured"],
    #                             reinfP[ 3, t]*foi_dt[3, ]*oldPop[3 ,"hcc_cured"],
    #                             reinfP[ 3, t]*foi_dt[3, ]*oldPop[3 ,"lt_cured"],
    #                             reinfP[ 3, t]*foi_dt[3, ]*oldPop[3 ,"plt_cured"])
    
    # newreinfection[4, t] <- sum(reinfP[ 4, t]*foi_dt[4, ]*oldPop[4 ,"a_cured"],
    #                            reinfP[ 4, t]*foi_dt[4, ]*oldPop[4,"f0_cured"],
    #                            reinfP[ 4, t]*foi_dt[4, ]*oldPop[4,"f1_cured"],
    #                            reinfP[ 4, t]*foi_dt[4, ]*oldPop[4,"f2_cured"],
    #                            reinfP[ 4, t]*foi_dt[4, ]*oldPop[4,"f3_cured"],
    #                            reinfP[ 4, t]*foi_dt[4, ]*oldPop[4 ,"f4_cured"],
    #                            reinfP[ 4, t]*foi_dt[4, ]*oldPop[4 ,"dc_cured"],
    #                            reinfP[ 4, t]*foi_dt[4, ]*oldPop[4 ,"hcc_cured"],
    #                            reinfP[ 4, t]*foi_dt[4, ]*oldPop[4 ,"lt_cured"],
    #                            reinfP[ 4, t]*foi_dt[4, ]*oldPop[4 ,"plt_cured"]) 
    
    # pop-array 
    newarray[ , , t]<- 
      (pop_array[, , t]*oldPop[, "s"] + 
         pop_array[, , t]*oldPop[, "a_undiag"] + 
         pop_array[, , t]*oldPop[, "f0_undiag"] + 
         pop_array[, , t]*oldPop[, "f1_undiag"] + 
         pop_array[, , t]*oldPop[, "f2_undiag"] + 
         pop_array[, , t]*oldPop[, "f3_undiag"] + 
         pop_array[, , t]*oldPop[, "f4_undiag"] + 
         pop_array[, , t]*oldPop[, "dc_undiag"] +  
         pop_array[, , t]*oldPop[, "hcc_undiag"] + 
         pop_array[, , t]*oldPop[, "lt_undiag"] +  
         pop_array[, , t]*oldPop[, "plt_undiag"] + 
         pop_array[, , t]*oldPop[, "a_diag_ab"] + 
         pop_array[, , t]*oldPop[, "f0_diag_ab"] + 
         pop_array[, , t]*oldPop[, "f1_diag_ab"] + 
         pop_array[, , t]*oldPop[, "f2_diag_ab"] + 
         pop_array[, , t]*oldPop[, "f3_diag_ab"] + 
         pop_array[, , t]*oldPop[, "f4_diag_ab"] + 
         pop_array[, , t]*oldPop[, "dc_diag_ab"] + 
         pop_array[, , t]*oldPop[, "hcc_diag_ab"]+ 
         pop_array[, , t]*oldPop[, "lt_diag_ab"] + 
         pop_array[, , t]*oldPop[, "plt_diag_ab"] + 
         pop_array[, , t]*oldPop[, "a_diag_RNA"] + 
         pop_array[, , t]*oldPop[, "f0_diag_RNA"] + 
         pop_array[, , t]*oldPop[, "f1_diag_RNA"] + 
         pop_array[, , t]*oldPop[, "f2_diag_RNA"] + 
         pop_array[, , t]*oldPop[, "f3_diag_RNA"] + 
         pop_array[, , t]*oldPop[, "f4_diag_RNA"] + 
         pop_array[, , t]*oldPop[, "dc_diag_RNA"] + 
         pop_array[, , t]*oldPop[, "hcc_diag_RNA"] + 
         pop_array[, , t]*oldPop[, "lt_diag_RNA"] + 
         pop_array[, , t]*oldPop[, "plt_diag_RNA"] +
         pop_array[, , t]*oldPop[, "a_treat"] + 
         pop_array[, , t]*oldPop[, "f0_treat"] + 
         pop_array[, , t]*oldPop[, "f1_treat"] + 
         pop_array[, , t]*oldPop[, "f2_treat"] + 
         pop_array[, , t]*oldPop[, "f3_treat"] + 
         pop_array[, , t]*oldPop[, "f4_treat"] + 
         pop_array[, , t]*oldPop[, "dc_treat"] + 
         pop_array[, , t]*oldPop[, "hcc_treat"] + 
         pop_array[, , t]*oldPop[, "lt_treat"] + 
         pop_array[, , t]*oldPop[, "plt_treat"] + 
         pop_array[, , t]*oldPop[, "a_treat_f"] + 
         pop_array[, , t]*oldPop[, "f0_treat_f"] + 
         pop_array[, , t]*oldPop[, "f1_treat_f"] + 
         pop_array[, , t]*oldPop[, "f2_treat_f"] + 
         pop_array[, , t]*oldPop[, "f3_treat_f"] + 
         pop_array[, , t]*oldPop[, "f4_treat_f"] + 
         pop_array[, , t]*oldPop[, "dc_treat_f"] + 
         pop_array[, , t]*oldPop[, "hcc_treat_f"] + 
         pop_array[, , t]*oldPop[, "lt_treat_f"] + 
         pop_array[, , t]*oldPop[, "plt_treat_f"] + 
         pop_array[, , t]*oldPop[, "a_cured"] + 
         pop_array[, , t]*oldPop[, "f0_cured"] + 
         pop_array[, , t]*oldPop[, "f1_cured"] + 
         pop_array[, , t]*oldPop[, "f2_cured"] + 
         pop_array[, , t]*oldPop[, "f3_cured"] + 
         pop_array[, , t]*oldPop[, "f4_cured"] + 
         pop_array[, , t]*oldPop[, "dc_cured"] + 
         pop_array[, , t]*oldPop[, "hcc_cured"] + 
         pop_array[, , t]*oldPop[, "lt_cured"] + 
         pop_array[, , t]*oldPop[, "plt_cured"])
    
    inflow[, t]<- 
      colSums(pop_array[, , t]*oldPop[, "s"])  + 
      colSums(pop_array[, , t]*oldPop[, "a_undiag"])  + 
      colSums(pop_array[, , t]*oldPop[, "f0_undiag"]) +
      colSums(pop_array[, , t]*oldPop[, "f1_undiag"]) + 
      colSums(pop_array[, , t]*oldPop[, "f2_undiag"]) + 
      colSums(pop_array[, , t]*oldPop[, "f3_undiag"]) + 
      colSums(pop_array[, , t]*oldPop[, "f4_undiag"]) + 
      colSums(pop_array[, , t]*oldPop[, "dc_undiag"]) + 
      colSums(pop_array[, , t]*oldPop[, "hcc_undiag"]) + 
      colSums(pop_array[, , t]*oldPop[, "lt_undiag"]) + 
      colSums(pop_array[, , t]*oldPop[, "plt_undiag"]) + 
      colSums(pop_array[, , t]*oldPop[, "a_diag_ab"]) + 
      colSums(pop_array[, , t]*oldPop[, "f0_diag_ab"]) +
      colSums(pop_array[, , t]*oldPop[, "f1_diag_ab"]) + 
      colSums(pop_array[, , t]*oldPop[, "f2_diag_ab"]) + 
      colSums(pop_array[, , t]*oldPop[, "f3_diag_ab"]) + 
      colSums(pop_array[, , t]*oldPop[, "f4_diag_ab"]) + 
      colSums(pop_array[, , t]*oldPop[, "dc_diag_ab"]) + 
      colSums(pop_array[, , t]*oldPop[, "hcc_diag_ab"]) + 
      colSums(pop_array[, , t]*oldPop[, "lt_diag_ab"]) + 
      colSums(pop_array[, , t]*oldPop[, "plt_diag_ab"]) + 
      colSums(pop_array[, , t]*oldPop[, "a_diag_RNA"]) + 
      colSums(pop_array[, , t]*oldPop[, "f0_diag_RNA"]) +
      colSums(pop_array[, , t]*oldPop[, "f1_diag_RNA"]) + 
      colSums(pop_array[, , t]*oldPop[, "f2_diag_RNA"]) + 
      colSums(pop_array[, , t]*oldPop[, "f3_diag_RNA"]) + 
      colSums(pop_array[, , t]*oldPop[, "f4_diag_RNA"]) + 
      colSums(pop_array[, , t]*oldPop[, "dc_diag_RNA"]) + 
      colSums(pop_array[, , t]*oldPop[, "hcc_diag_RNA"]) + 
      colSums(pop_array[, , t]*oldPop[, "lt_diag_RNA"]) + 
      colSums(pop_array[, , t]*oldPop[, "plt_diag_RNA"]) + 
      colSums(pop_array[, , t]*oldPop[, "a_treat"]) + 
      colSums(pop_array[, , t]*oldPop[, "f0_treat"]) +
      colSums(pop_array[, , t]*oldPop[, "f1_treat"]) + 
      colSums(pop_array[, , t]*oldPop[, "f2_treat"]) + 
      colSums(pop_array[, , t]*oldPop[, "f3_treat"]) + 
      colSums(pop_array[, , t]*oldPop[, "f4_treat"]) + 
      colSums(pop_array[, , t]*oldPop[, "dc_treat"]) + 
      colSums(pop_array[, , t]*oldPop[, "hcc_treat"]) + 
      colSums(pop_array[, , t]*oldPop[, "lt_treat"]) + 
      colSums(pop_array[, , t]*oldPop[, "plt_treat"]) + 
      colSums(pop_array[, , t]*oldPop[, "a_treat_f"]) + 
      colSums(pop_array[, , t]*oldPop[, "f0_treat_f"]) +
      colSums(pop_array[, , t]*oldPop[, "f1_treat_f"]) + 
      colSums(pop_array[, , t]*oldPop[, "f2_treat_f"]) + 
      colSums(pop_array[, , t]*oldPop[, "f3_treat_f"]) + 
      colSums(pop_array[, , t]*oldPop[, "f4_treat_f"]) + 
      colSums(pop_array[, , t]*oldPop[, "dc_treat_f"]) + 
      colSums(pop_array[, , t]*oldPop[, "hcc_treat_f"]) + 
      colSums(pop_array[, , t]*oldPop[, "lt_treat_f"]) + 
      colSums(pop_array[, , t]*oldPop[, "plt_treat_f"]) + 
      colSums(pop_array[, , t]*oldPop[, "a_cured"]) + 
      colSums(pop_array[, , t]*oldPop[, "f0_cured"]) +
      colSums(pop_array[, , t]*oldPop[, "f1_cured"]) + 
      colSums(pop_array[, , t]*oldPop[, "f2_cured"]) + 
      colSums(pop_array[, , t]*oldPop[, "f3_cured"]) + 
      colSums(pop_array[, , t]*oldPop[, "f4_cured"]) + 
      colSums(pop_array[, , t]*oldPop[, "dc_cured"]) + 
      colSums(pop_array[, , t]*oldPop[, "hcc_cured"]) + 
      colSums(pop_array[, , t]*oldPop[, "lt_cured"]) + 
      colSums(pop_array[, , t]*oldPop[, "plt_cured"]) 
    
    
    outflow[, t]<- 
      rowSums(pop_array[, , t]*oldPop[, "s"])  + 
      rowSums(pop_array[, , t]*oldPop[, "a_undiag"])  + 
      rowSums(pop_array[, , t]*oldPop[, "f0_undiag"]) +
      rowSums(pop_array[, , t]*oldPop[, "f1_undiag"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f2_undiag"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f3_undiag"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f4_undiag"]) + 
      rowSums(pop_array[, , t]*oldPop[, "dc_undiag"]) + 
      rowSums(pop_array[, , t]*oldPop[, "hcc_undiag"]) + 
      rowSums(pop_array[, , t]*oldPop[, "lt_undiag"]) + 
      rowSums(pop_array[, , t]*oldPop[, "plt_undiag"]) + 
      rowSums(pop_array[, , t]*oldPop[, "a_diag_ab"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f0_diag_ab"]) +
      rowSums(pop_array[, , t]*oldPop[, "f1_diag_ab"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f2_diag_ab"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f3_diag_ab"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f4_diag_ab"]) + 
      rowSums(pop_array[, , t]*oldPop[, "dc_diag_ab"]) + 
      rowSums(pop_array[, , t]*oldPop[, "hcc_diag_ab"]) + 
      rowSums(pop_array[, , t]*oldPop[, "lt_diag_ab"]) + 
      rowSums(pop_array[, , t]*oldPop[, "plt_diag_ab"]) + 
      rowSums(pop_array[, , t]*oldPop[, "a_diag_RNA"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f0_diag_RNA"]) +
      rowSums(pop_array[, , t]*oldPop[, "f1_diag_RNA"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f2_diag_RNA"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f3_diag_RNA"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f4_diag_RNA"]) + 
      rowSums(pop_array[, , t]*oldPop[, "dc_diag_RNA"]) + 
      rowSums(pop_array[, , t]*oldPop[, "hcc_diag_RNA"]) + 
      rowSums(pop_array[, , t]*oldPop[, "lt_diag_RNA"]) + 
      rowSums(pop_array[, , t]*oldPop[, "plt_diag_RNA"]) + 
      rowSums(pop_array[, , t]*oldPop[, "a_treat"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f0_treat"]) +
      rowSums(pop_array[, , t]*oldPop[, "f1_treat"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f2_treat"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f3_treat"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f4_treat"]) + 
      rowSums(pop_array[, , t]*oldPop[, "dc_treat"]) + 
      rowSums(pop_array[, , t]*oldPop[, "hcc_treat"]) + 
      rowSums(pop_array[, , t]*oldPop[, "lt_treat"]) + 
      rowSums(pop_array[, , t]*oldPop[, "plt_treat"]) + 
      rowSums(pop_array[, , t]*oldPop[, "a_treat_f"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f0_treat_f"]) +
      rowSums(pop_array[, , t]*oldPop[, "f1_treat_f"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f2_treat_f"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f3_treat_f"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f4_treat_f"]) + 
      rowSums(pop_array[, , t]*oldPop[, "dc_treat_f"]) + 
      rowSums(pop_array[, , t]*oldPop[, "hcc_treat_f"]) + 
      rowSums(pop_array[, , t]*oldPop[, "lt_treat_f"]) + 
      rowSums(pop_array[, , t]*oldPop[, "plt_treat_f"]) + 
      rowSums(pop_array[, , t]*oldPop[, "a_cured"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f0_cured"]) +
      rowSums(pop_array[, , t]*oldPop[, "f1_cured"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f2_cured"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f3_cured"]) + 
      rowSums(pop_array[, , t]*oldPop[, "f4_cured"]) + 
      rowSums(pop_array[, , t]*oldPop[, "dc_cured"]) + 
      rowSums(pop_array[, , t]*oldPop[, "hcc_cured"]) + 
      rowSums(pop_array[, , t]*oldPop[, "lt_cured"]) + 
      rowSums(pop_array[, , t]*oldPop[, "plt_cured"]) 
    
    
    
    #### cost ####
    if (!is.null(cost)){ 
    costPops[, , t] <- allPops[, , t]*cost$state[,,t]
   
    
    costTestingAb[, t] <- costflow[,"ctau_ab", t]*
      (newTestingAb[, t] + 
         ((oldPop[, "s"] + oldPop[, "a_cured"])*tau_ab_dt[, "f0", t]))
    
    if(tau_poct_dt[2, "f0", t]==0){
      costTestingAg[, t] <- costflow[,"ctau_ag", t]*
      (newTestingAg[, t] + 
         ((S[,] - oldPop[,"s"] - oldPop[, "a_cured"])*tau_ag_dt[, "f0", t]))
    }else if(tau_poct_dt[2, "f0", t]!=0){ 
      costTestingAg[, t] <- costflow[,"ctau_ag", t]*newTestingAg[, t]
      
      }
    
    costnewTestingPOCT[, t] <- costflow[,"ctau_poct",t ]*
      (newTestingPOCT[, t] + S[,]*tau_poct_dt[, "f0", t])
    
    
    
    costTreatment[, t] <- newTreatment[, t]*costflow[,"ceta", t]
    
    
    costCured[, t] <- newCured[, t]*costflow[,"ccured", t]
    
    
    costRetreat[, t] <- newRetreat[, t]*costflow[,"crho", t]
    
    
    QALYPops[, , t] <- allPops[, , t]*cost$QALY[,,t]
    
    }
    
  }
  
  #------------------------- results output ---------------------------------#
  
  if (!is.null(cost)){
    results <- list(allPops = allPops,
                    newS = newS,
                    newEntry = newEntry,
                    newDeath = newDeath,
                    newInfections = newInfections, 
                    newHCVdeaths = newHCVdeaths, 
                    newTreatment = newTreatment, 
                    newRetreat = newRetreat, 
                    newTestingAb = newTestingAb, 
                    newTestingAg = newTestingAg,
                    newTestingPOCT = newTestingPOCT,
                    newCured = newCured,
                    newreinfection = newreinfection,
                    newpop_tran = newarray,
                    inflow = inflow,
                    outflow = outflow,
                    death_hcv = death_hcv,
                    HCVdeathState = HCVdeathState,
                    newDeathState = newDeathState,
                    costPops = costPops,
                    QALYPops = QALYPops,
                    costTestingAb = costTestingAb,
                    costTestingAg = costTestingAg,
                    costnewTestingPOCT = costnewTestingPOCT,
                    costTreatment = costTreatment,
                    costCured = costCured,
                    costRetreat = costRetreat)
  }else{ 
    results <- list(allPops = allPops,
                    newS = newS,
                    newEntry = newEntry,
                    newDeath = newDeath,
                    newInfections = newInfections, 
                    newHCVdeaths = newHCVdeaths, 
                    newTreatment = newTreatment, 
                    newRetreat = newRetreat, 
                    newTestingAb = newTestingAb, 
                    newTestingAg = newTestingAg,
                    newTestingPOCT = newTestingPOCT,
                    newCured = newCured,
                    newreinfection = newreinfection,
                    newpop_tran = newarray,
                    inflow = inflow,
                    outflow = outflow,
                    death_hcv = death_hcv,
                    HCVdeathState = HCVdeathState,
                    newDeathState = newDeathState)
    
    
    }
  
  
  
}