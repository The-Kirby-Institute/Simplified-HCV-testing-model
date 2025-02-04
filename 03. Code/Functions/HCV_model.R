library(dplyr)
library(stringr)
# simulate the HCV model equations

HCVMSM <- function(HCV, parama, initialPop, disease_progress,
                   pop_array, param_cascade, param_cascade_sc, fib,  end_Y = NULL, modelrun=NULL,
                   scenario = NULL, cost = NULL, costflow = NULL, 
                   costflow_Neg = NULL, proj=NULL, fc_sc = NULL, fp = NULL){
  
  # Args:
  #       HCV: a list containing project specifications 
  #       tp: the simulation t points in year
  #       parama: paramaeters did not change by time
  #       initialPop: initial population size
  #       
  #        
  #     pop_array:  
  #     param_cascade: a list of array includes those paramaeters that varies across 
  #               population and disease progress, (tau_ab, tau_RNA, tau_poct, 
  #               eta, lota, rho, cured)  
  #              * if intervention applied, the parameters equals to base case + scenarios (total effect)
  #    param_cascade_sc:  intervention parameters 
  #     
  #    scenario: testing the reinfection reduction related to behavior intervnetions
  #                            
  # cost: cost attached to compartments 
  # costflow: cost attached to flows 
  # costflow_Neg: cost for those not living with HCV 
  # Proj: project name 
  # simulation time period 
  # fc: a fraction for the people living without HCV received test 
  dt <- HCV$timestep
  
  
  if (is.null(end_Y)){ 
    
    years <- HCV$nyears
    
    npts <-HCV$npts - 1 
    
  } else if (!is.data.frame(parama)) {
    
    years <- HCV$startYear:end_Y
    pts <- seq(HCV$startYear, end_Y, by = dt)
    
    #pts <- head(pts, -1)
    pts <- pts[-1]
    npts <- length(pts)
    
    parama <- lapply(parama, function(x) x[1:(npts + 1)])
    
    pop_array <- pop_array[ , , 1:npts] 
    
    # cost <- lapply(cost, function(x) x[, , 1:(npts)])
    
    # costflow <- lapply(costflow, function(x) x[ , , 1:npts])
    
    
  }else{
    
    years <- HCV$startYear:end_Y
    
    pts <- seq(HCV$startYear, end_Y, by = dt)
    
    #pts <- head(pts, -1)
    pts <- pts[-1]
    npts <- length(pts)
    
    parama <- parama[1:(npts+1),]
    
    pop_array <- pop_array[ , , 1:npts] 
    
    # cost <- lapply(cost, function(x) x[, , 1:(npts)])
    
    # costflow <- lapply(costflow, function(x) x[ , , 1:npts])
    
    
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
    reinfP[ , ] <- 1
  
    
  } else{ 
    reinfP <- scenario
    
  }
  
  
  # testing distribution weigthed 
  fc <- matrix(0, ncol = npts + 1, nrow = npops)
  
  if (is.null(fc_sc)){
    fc[i, ] <- rep(1, npts + 1)
    
  } 
  else{ 
    
    fc <- fc_sc
    
    }
  
  
  
  if (is.null(fp)){
    fp <- 0
    
  }
  #### entry ####
  entry1 <- matrix(0, ncol = ncomponent, nrow = npops, dimnames = dimNames)
  death <- matrix(0, ncol = ncomponent, nrow = npops, dimnames = dimNames)
  death_hcv <- matrix(0, ncol = ncomponent, nrow = npops, 
                      dimnames = dimNames)
  
  leave <- matrix(0, ncol = ncomponent, nrow = npops, dimnames = dimNames)
  ####create result matrix####
  ResultMatrix <- matrix(0, ncol = npts, nrow = npops)
  rownames(ResultMatrix) <- popNames
  
  newS <- ResultMatrix
  
  newEntry <- ResultMatrix
  
  newDeath <-ResultMatrix
  
  newLeave <-ResultMatrix
  
  newInfections <- ResultMatrix
  
  newHCVdeaths <- ResultMatrix
  
  newHCVdeathsState <- ResultMatrix
  
  newTreatment <- ResultMatrix 
  
  newRetreat <- ResultMatrix
  
  newTestingAb <- ResultMatrix
  
  newTestingAb_neg <- ResultMatrix
  
  newTestingAg <- ResultMatrix
  
  newTestingAg_neg <- ResultMatrix
  
  newTestingPOCT <- ResultMatrix
  
  newTestingPOCT_neg <- ResultMatrix
  
  newCured <- ResultMatrix
  
  newtreatfailed <- ResultMatrix
  
  newreinfection <- ResultMatrix
  
  newreinfection_chronic <- ResultMatrix
  
  inflow <- ResultMatrix
  
  inflow_hcv <- ResultMatrix
  
  outflow <-ResultMatrix 
  
  outflow_hcv <- ResultMatrix
  
  # intervention 
  newTestingAb_sc <- ResultMatrix
  
  newTestingAg_sc <- ResultMatrix
  
  newTestingPOCT_sc <- ResultMatrix 
  
  newTreatment_sc <- ResultMatrix
  
  newTestingAb_sc_neg <- ResultMatrix
  
  newTestingAg_sc_neg <- ResultMatrix
  
  newTestingPOCT_sc_neg <- ResultMatrix 
  

  
  if (!is.null(cost)){ 
    costTestingAb <- ResultMatrix
    costTestingAg <- ResultMatrix
    costTestingPOCT <- ResultMatrix
    costTreatment <- ResultMatrix
    costCured <- ResultMatrix
    costRetreat <- ResultMatrix 
    
    costTestingAb_sc <- ResultMatrix
    costTestingAg_sc <- ResultMatrix
    costTestingPOCT_sc <- ResultMatrix
    costTreatment_sc <- ResultMatrix
    
    costPops <- array(0, c(npops, ncomponent, npts), dimnames = dimNames) 
    
    QALYPops <- array(0, c(npops, ncomponent, npts), dimnames = dimNames) 
    
    
    
    
    }
  
  #### pop_array ####
  popArray <- array(0, c(npops, npops, npts), dimnames = list(popNames, popNames))
  
  # pop transition in each compartment 
  newarray_state <- array(0, c(ncomponent, npops, npops, npts), 
                          dimnames = list(componentName, popNames, popNames))
  
  
  newarray <- popArray
  
  ####demographic paramaeters####
  # entry 
  
  
  ####mortality =rate to leave the model####
  ## 
  morb <- matrix(0, ncol = npts + 1, nrow = npops)
  
  mordc <- matrix(0, ncol = npts + 1, nrow = npops)
  
  morhcc <- matrix(0, ncol = npts + 1, nrow = npops)
  
  morlt <- matrix(0, ncol = npts + 1, nrow = npops)
  
  morplt <- matrix(0, ncol = npts + 1, nrow = npops)
  
  l <- matrix(0, ncol = npts + 1, nrow = npops)

  
  for(i in 1:npops){
    morb[i, ] <- parama[, paste0("morb", i)]
    
    mordc[i, ] <- parama[ , paste0("mordc", i)]
    
    morhcc[i, ] <- parama[ , paste0("morhcc", i)]
    
    morlt[i, ] <- parama[ , paste0("morlt", i)]
    
    morplt[i, ] <- parama[ , paste0("morplt", i)]
    
    l[i, ] <- parama[ , paste0("leave", i)]
  }
  
  morb_dt <-1-(1-morb)^dt
  
  mordc_dt <- 1-(1-mordc)^dt
  
  morhcc_dt <- 1-(1-morhcc)^dt
  
  morlt_dt <- 1-(1-morlt)^dt
  
  morplt_dt <- 1-(1-morplt)^dt
  
  leave_dt <- 1-(1-l)^dt

  #####cured####
  mordcCure <- matrix(0, ncol = npts + 1, nrow = npops)
  
  morhccCure <- matrix(0, ncol = npts + 1, nrow = npops)
  
  for(i in 1:npops){
    mordcCure[i, ] <- mordc[i,]*parama$Cure_mordc_Reduction
    morhccCure[i, ]<- morhcc[i,]*(1-parama$Cure_morhcc_Reduction)
  }
   
  mordcCure_dt <- 1-(1-mordcCure)^dt
  
  morhccCure_dt <- 1-(1-morhccCure)^dt
  
  ## temporarily ignore HIV excess HIV  
  
  
  #### disease progress ####
  
  transition <-matrix(0,ncol = diseaseprogress_n, nrow = npops)
  colnames(transition) <- diseaseprogress_Name
  
  trans <- as.matrix(disease_progress, nrow = npops, ncol= diseaseprogress_n)
  
  transition[, ]<- trans[, ]
  
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
  
  for(i in 1:npops){
    spc1[i, ] <- parama[ , paste0("spc", i)]
  }
  
  spc1_dt <- 1-(1-spc1)^dt 
  
  ####testing rate#### 
  #####antibody####
  dimN <- list(popNames, progressName) 
  tau_ab <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  tau_ab <- param_cascade$tau_ab
  tau_ab_dt <- 1-(1-tau_ab)^dt
  
  #####ag/RNA (second step testing)####
  tau_RNA <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  tau_RNA <- param_cascade$tau_RNA
  tau_RNA_dt <- 1-(1-tau_RNA)^dt
  
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

  lota_dt[, c("f0", "f1", "f2", "f3"), ] <- lota_dt[, c("f0", "f1", "f2", "f3"), ]*(1-parama$SVR[1])
  lota_dt[, c("f4", "dc", "hcc", "lt", "plt"), ] <- lota_dt[, c("f4", "dc", "hcc", "lt", "plt"), ]*(1-parama$SVRf4[1])

  
  ####re-treated####
  rho <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  rho <- param_cascade$rho
  rho_dt <- 1-(1-rho)^dt
  
  ####cured#### 
  cure <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  cure <- param_cascade$cured
  cure_dt <- 1-(1-cure)^dt
  cure_dt[, c("f0", "f1", "f2", "f3"), ] <- cure_dt[, c("f0", "f1", "f2", "f3"), ]*parama$SVR[1]
  cure_dt[, c("f4", "dc", "hcc", "lt", "plt"), ] <- cure_dt[, c("f4", "dc", "hcc", "lt", "plt"), ]*parama$SVRf4[1]
  
  # intervention 
  ####testing rate#### 
  #####antibody####
  dimN <- list(popNames, progressName) 
  tau_ab_sc <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  tau_ab_sc <- param_cascade_sc$tau_ab
  tau_ab_sc_dt <- 1-(1-tau_ab_sc)^dt
  
  #####ag/RNA (second step testing)####
  tau_RNA_sc <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  tau_RNA_sc <- param_cascade_sc$tau_RNA
  tau_RNA_sc_dt <- 1-(1-tau_RNA_sc)^dt
  
  #####POCT: undiag>>diag_RNA####
  tau_poct_sc <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  tau_poct_sc <- param_cascade_sc$tau_poct
  tau_poct_sc_dt <- 1-(1-tau_poct_sc)^dt
  
  ####treatment####
  eta_sc <- array(0, c(npops, nprogress, npts +1), dimnames = dimN)
  eta_sc <- param_cascade_sc$eta
  eta_sc_dt <- 1-(1-eta_sc)^dt
  
  
  ####population transition#### 
  pop_array <- pop_array*dt 
  
  #### cost dt ####
  if(!is.null(cost) &!is.null(costflow_Neg)){
    cost <- lapply(cost, function(x){ 
      
      a <- x*dt
    })
    
    costflow <- costflow
    costflow_Neg <- costflow_Neg
  }
  else if(!is.null(cost) & is.null(costflow_Neg)){
    cost <- lapply(cost, function(x){ 
      
      a <- x*dt
    })
    
    costflow <- costflow
    costflow_Neg <- lapply(costflow, function(x){ x*0})
  }
  
  entry <- matrix(0, ncol = npts + 1, nrow = npops)
  
  if(isTRUE(!proj%in% c("POC_AU", "TWPrisoners")) & isTRUE(modelrun != "steady")){
    entry[1,] <- parama$entry
    
    
    for(i in 2: npops){ 
      entry[i,] <- 0
    }
    entry_dt <- entry*dt  
  }
  
  
  #----------------------------------------------------------------------------#
  ####                     Equations                                        ####   
  #----------------------------------------------------------------------------#
  for (t in 2: npts){ 
    
    oldPop <- allPops[, , t-1]
    newPop <- allPops[, , t]
    
    #####FOI#####
    
    foi <- matrix(0, ncol = 1, nrow = npops)
    for(i in 1:npops){
      foi[i, ] <- parama[ , paste0("beta", i)][t]
    }
    
    foi_dt <- foi*dt 
    
    
    ####leave model  ####
    leave[,] <- leave_dt[,t]*oldPop[,] 
    ##### mortality ####
    death[,] <- morb_dt[,t]*oldPop[,] 
    
    death_hcv[,"s"] <- 0*oldPop[, "s"]
    death_hcv[,"a_undiag"] <- 0*oldPop[, "a_undiag"]
    death_hcv[,"a_diag_ab"] <- 0*oldPop[, "a_diag_ab"]
    death_hcv[,"a_diag_RNA"] <- 0*oldPop[, "a_diag_RNA"]
    death_hcv[,"a_treat"] <- 0*oldPop[, "a_treat"]
    death_hcv[,"a_treat_f"] <- 0*oldPop[, "a_treat_f"]
    death_hcv[,"a_cured"] <- 0*oldPop[, "a_cured"]
    
    death_hcv[,"f0_undiag"] <- 0*oldPop[, "f0_undiag"]
    death_hcv[,"f0_diag_ab"] <- 0*oldPop[, "f0_diag_ab"]
    death_hcv[,"f0_diag_RNA"] <- 0*oldPop[, "f0_diag_RNA"]
    death_hcv[,"f0_treat"] <- 0*oldPop[, "f0_treat"]
    death_hcv[,"f0_treat_f"] <- 0*oldPop[, "f0_treat_f"]
    death_hcv[,"f0_cured"] <- 0*oldPop[, "f0_cured"]
    
    death_hcv[,"f1_undiag"] <- 0*oldPop[, "f1_undiag"]
    death_hcv[,"f1_diag_ab"] <- 0*oldPop[, "f1_diag_ab"]
    death_hcv[,"f1_diag_RNA"] <- 0*oldPop[, "f1_diag_RNA"]
    death_hcv[,"f1_treat"] <- 0*oldPop[, "f1_treat"]
    death_hcv[,"f1_treat_f"] <- 0*oldPop[, "f1_treat_f"]
    death_hcv[,"f1_cured"] <- 0*oldPop[, "f1_cured"]
    
    death_hcv[,"f2_undiag"] <- 0*oldPop[, "f2_undiag"]
    death_hcv[,"f2_diag_ab"] <- 0*oldPop[, "f2_diag_ab"]
    death_hcv[,"f2_diag_RNA"] <- 0*oldPop[, "f2_diag_RNA"]
    death_hcv[,"f2_treat"] <- 0*oldPop[, "f2_treat"]
    death_hcv[,"f2_treat_f"] <- 0*oldPop[, "f2_treat_f"]
    death_hcv[,"f2_cured"] <- 0*oldPop[, "f2_cured"]
    
    death_hcv[,"f3_undiag"] <- 0*oldPop[, "f3_undiag"]
    death_hcv[,"f3_diag_ab"] <- 0*oldPop[, "f3_diag_ab"]
    death_hcv[,"f3_diag_RNA"] <- 0*oldPop[, "f3_diag_RNA"]
    death_hcv[,"f3_treat"] <- 0*oldPop[, "f3_treat"]
    death_hcv[,"f3_treat_f"] <- 0*oldPop[, "f3_treat_f"]
    death_hcv[,"f3_cured"] <- 0*oldPop[, "f3_cured"]
    
    death_hcv[,"f4_undiag"] <- 0*oldPop[, "f4_undiag"]
    death_hcv[,"f4_diag_ab"] <- 0*oldPop[, "f4_diag_ab"]
    death_hcv[,"f4_diag_RNA"] <- 0*oldPop[, "f4_diag_RNA"]
    death_hcv[,"f4_treat"] <- 0*oldPop[, "f4_treat"]
    death_hcv[,"f4_treat_f"] <- 0*oldPop[, "f4_treat_f"]
    death_hcv[,"f4_cured"] <- 0*oldPop[, "f4_cured"]
    
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
    
    S <- matrix(0, ncol = 1, nrow = npops) 
    
    Sc <- matrix(0, ncol = 1, nrow = npops)
    
    # I
    I <- matrix(0, ncol = 1, nrow = npops)
    
    for(i in 1: npops){ 
      N[i, ] <- sum(oldPop[i, ], na.rm = TRUE)
      S[i, ] <- (oldPop[i , "s"] + 
                   oldPop[i , "a_cured"] +
                   oldPop[i , "f0_cured"] +
                   oldPop[i , "f1_cured"] +
                   oldPop[i , "f2_cured"] +
                   oldPop[i , "f3_cured"] +
                   oldPop[i , "f4_cured"] +
                   oldPop[i , "dc_cured"] +
                   oldPop[i , "hcc_cured"] +
                   oldPop[i , "lt_cured"] +
                   oldPop[i , "plt_cured"])
      
      Sc[i, ] <- (oldPop[i , "a_cured"] +
                   oldPop[i , "f0_cured"] +
                   oldPop[i , "f1_cured"] +
                   oldPop[i , "f2_cured"] +
                   oldPop[i , "f3_cured"] +
                   oldPop[i , "f4_cured"] +
                   oldPop[i , "dc_cured"] +
                   oldPop[i , "hcc_cured"] +
                   oldPop[i , "lt_cured"] +
                   oldPop[i , "plt_cured"])
      }
    
    for(i in 1: npops){ 
      
      I[i, ] <- N[i, ] - S[i,]
      
      }
    
    Sall <- matrix(0, ncol = 1, nrow = 1)
    
    Iall <- matrix(0, ncol = 1, nrow = 1)
    
    Sall <- colSums(S)
    
    Iall <- colSums(I)
    
    II <- matrix(0, ncol = 1, nrow = npops)
    
    if(proj %in% c("POC_AU", "TWPrisoners")){ 
      II[1,] <- (I[1,] + I[2,])/(N[1, ] + N[2, ])
      II[2,] <- (I[1,] + I[2,])/(N[1, ] + N[2, ])
      II[3,] <- (I[3,] + I[4,] + 0*I[5, ])/(N[3, ] + N[4, ] + 0*N[5, ])
      II[4,] <- (I[3,] + I[4,] + 0*I[5, ])/(N[3, ] + N[4, ] + 0*N[5, ])
      II[5,] <- (0*I[3,] + 0*I[4,] + I[5, ])/(0*N[3, ] + 0*N[4, ] + N[5, ])
      
    }else{ 
      for(i in 1: npops){   
        II[i, ] <-  Iall/(Sall + Iall)
      }
    }
   
     
     

  
    II[is.nan(II)] <- 0
    
    
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
    
    #####Entry of new people####
    
    if(isTRUE(proj == "POC_AU")){ # assuming no HCV death in nonPWID in prison 
      
      
      entry1[1,"s" ] <- sum(death[1, ])  +  sum(death[3, ])  + 
        sum(death_hcv[1, ]) +  sum(death_hcv[3, ]) + 
        sum(leave[1, ])  + sum(leave[3, ]) + 
        sum(death[2, ]) + sum(death[4, ]) + 
        sum(death_hcv[2, ]) + sum(leave[2, ]) + 
        sum(death_hcv[4, ])  + sum(leave[4, ])
    
     
      entry1[5, ] <- sum(leave[5, ])*oldPop[5, ]/sum(oldPop[5, ])

      

      
      }
    else if(isTRUE(proj == "TWPrisoners")){
      entry1[1, "s"] <- sum(death[1, ])  + sum(death[2, ]) + 
        sum(death[3, ])  + sum(death[4, ]) + 
        sum(death[5, ])  + sum(death[6, ]) +
        sum(death_hcv[1, ]) + sum(death_hcv[2, ]) + 
        sum(death_hcv[3, ]) + sum(death_hcv[4, ]) +
        sum(death_hcv[5, ]) + sum(death_hcv[6, ]) + 
        sum(leave[1, ]) + sum(leave[2, ]) + 
        sum(leave[3, ]) + sum(leave[4, ]) +
        sum(leave[5, ]) + sum(leave[6, ]) 
      
    }
    else { entry1[1,"s"] <- entry_dt[1,t]
      
      
      }
 
    #if(!is.null(scenario)){ reinfP[ , ] <- reinfP[ ,t]}
    #else{ reinfP[ ,t] <- c(1,1,1,1)}
    
    
    #####S#### 
    newPop[,"s"] <- entry1[,"s"] + oldPop[,"s"] - 
      death[,"s"] -leave[,"s"] - death_hcv[, "s"] + 
      (colSums(pop_array[, , t]*oldPop[, "s"]) -
         rowSums(pop_array[, , t]*oldPop[, "s"])) - 
      #foi_dt[, ]*I[,]*Ps[,]
      foi_dt[, ]*II[,]*oldPop[, "s"]
    
    ##### a ####
    # a_infection 
    
    newPop[, "a_undiag"] <- entry1[,"a_undiag"] +oldPop[ ,"a_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "a_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_undiag"])) - 
      death[,"a_undiag"] -leave[,"a_undiag"] - death_hcv[, "a_undiag"] -
      transition_dt[, "a_f0"]*oldPop[, "a_undiag"] - 
      spc1_dt[, t]*oldPop[ ,"a_undiag"] -
      tau_ab_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[, "a_undiag"] - 
      tau_ab_sc_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[, "a_undiag"] - 
      tau_poct_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[, "a_undiag"]- 
      tau_poct_sc_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[, "a_undiag"]+ 
      #foi_dt[, ]*I[,]*Ps[,] +
      #reinfP[, t]*foi_dt[, ]*I[,]*Pca[,] + 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf0[,]
      foi_dt[, ]*II[,]*oldPop[, "s"]+ 
      reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "a_cured"] +
      reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "f0_cured"]
    
    
    
    ## a_testing, ab+
    newPop[, "a_diag_ab"] <- entry1[,"a_diag_ab"] + oldPop[,"a_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "a_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_diag_ab"])) -
      death[,"a_diag_ab"] -leave[,"a_diag_ab"] - death_hcv[, "a_diag_ab"] - 
      transition_dt[, "a_f0"]*oldPop[,"a_diag_ab"] + 
      tau_ab_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[,"a_undiag"] - 
      tau_RNA_dt[, "a", t]*oldPop[,"a_diag_ab"] + 
      tau_ab_sc_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[,"a_undiag"]- 
      tau_RNA_sc_dt[, "a", t]*tau_ab_sc_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[,"a_undiag"]
    
    ## a_testing, RNA+
    newPop[, "a_diag_RNA"] <- entry1[,"a_diag_RNA"] + oldPop[ ,"a_diag_RNA"] +
      (colSums(pop_array[, , t]*oldPop[, "a_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_diag_RNA"])) -
      death[,"a_diag_RNA"] -leave[,"a_diag_RNA"] - death_hcv[, "a_diag_RNA"] - 
      transition_dt[, "a_f0"]*oldPop[, "a_diag_RNA"] + 
      tau_RNA_dt[, "a", t]*oldPop[, "a_diag_ab"] +
      tau_poct_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[, "a_undiag"] -
      eta_dt[, "a", t]*oldPop[,"a_diag_RNA"] + 
      tau_RNA_sc_dt[, "a", t]*tau_ab_sc_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[,"a_undiag"]+
      tau_poct_sc_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[, "a_undiag"]-
      eta_sc_dt[, "a", t]*(tau_RNA_sc_dt[, "a", t]*tau_ab_sc_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[,"a_undiag"]+ 
                             tau_poct_sc_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[, "a_undiag"])
    
    ## a_testing, treat
    newPop[, "a_treat"] <- entry1[,"a_treat"] + oldPop[,"a_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "a_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_treat"])) - 
      death[,"a_treat"] -leave[,"a_treat"] - death_hcv[, "a_treat"] - 
      transition_dt[, "a_f0"]*oldPop[,"a_treat"] +
      eta_dt[, "a", t]*oldPop[,"a_diag_RNA"] +
      eta_sc_dt[, "a", t]*(tau_RNA_sc_dt[, "a", t]*tau_ab_sc_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[,"a_undiag"]+ 
                             tau_poct_sc_dt[, "a", t]*(1-spc1_dt[, t])*oldPop[, "a_undiag"]) - 
      lota_dt[, "a", t]*oldPop[,"a_treat"] -
      cure_dt[, "a", t]*oldPop[,"a_treat"] +
      rho_dt[, "a", t]*oldPop[,"a_treat_f"]
    
    # a_treat_failed  
    newPop[, "a_treat_f"] <- entry1[,"a_treat_f"]+ oldPop[,"a_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "a_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_treat_f"])) - 
      death[,"a_treat_f"] -leave[,"a_treat_f"] - death_hcv[, "a_treat_f"] - 
      transition_dt[, "a_f0"]*oldPop[,"a_treat_f"] -
      rho_dt[, "a", t]*oldPop[,"a_treat_f"] +
      lota_dt[, "a", t]*oldPop[,"a_treat"]
    
    # a_cured 
    newPop[, "a_cured"] <- entry1[,"a_cured"] + oldPop[,"a_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "a_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "a_cured"])) - 
      death[,"a_cured"] -leave[,"a_cured"] - death_hcv[, "a_cured"] + 
      cure_dt[, "a", t]*oldPop[, "a_treat"] + 
      spc1_dt[, t]*oldPop[ ,"a_undiag"] -
      #reinfP[, t]*foi_dt[, ]*I[,]*Pca[,]
      reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "a_cured"]
    
    
    ##### f0 #####
    ## f0 infection
    newPop[, "f0_undiag"] <- entry1[,"f0_undiag"] + oldPop[ ,"f0_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_undiag"])) - 
      death[,"f0_undiag"] -leave[,"f0_undiag"] - death_hcv[, "f0_undiag"] +
      transition_dt[, "a_f0"]*oldPop[, "a_undiag"] -
      transition_dt[, "f0_f1"]*oldPop[, "f0_undiag"] -
      tau_ab_dt[, "f0", t]*oldPop[, "f0_undiag"] -  
      tau_poct_dt[, "f0", t]*oldPop[, "f0_undiag"] -
      tau_ab_sc_dt[, "f0", t]*oldPop[, "f0_undiag"] -  
      tau_poct_sc_dt[, "f0", t]*oldPop[, "f0_undiag"] 
    
    
    ## f0 testing, ab+
    newPop[, "f0_diag_ab"] <- entry1[,"f0_diag_ab"] + oldPop[,"f0_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_diag_ab"])) - 
      death[,"f0_diag_ab"] -leave[,"f0_diag_ab"] - death_hcv[, "f0_diag_ab"] + 
      transition_dt[, "a_f0"]*oldPop[, "a_diag_ab"] -
      transition_dt[, "f0_f1"]*oldPop[, "f0_diag_ab"] +
      tau_ab_dt[, "f0", t]*oldPop[, "f0_undiag"] -
      tau_RNA_dt[, "f0", t]*oldPop[, "f0_diag_ab"] +
      tau_ab_sc_dt[, "f0", t]*oldPop[, "f0_undiag"]-
      tau_RNA_sc_dt[, "f0", t]*tau_ab_sc_dt[, "f0", t]*oldPop[, "f0_undiag"]   
    
    
    
    ## f0 testing, ag+/RNA
    newPop[, "f0_diag_RNA"] <- entry1[,"f0_diag_RNA"] + oldPop[ ,"f0_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_diag_RNA"])) -
      death[,"f0_diag_RNA"] -leave[,"f0_diag_RNA"] - death_hcv[, "f0_diag_RNA"] +
      transition_dt[, "a_f0"]*oldPop[, "a_diag_RNA"] -
      transition_dt[, "f0_f1"]*oldPop[, "f0_diag_RNA"] + 
      tau_RNA_dt[, "f0", t]*oldPop[, "f0_diag_ab"] +
      tau_poct_dt[, "f0", t]*oldPop[, "f0_undiag"] -
      eta_dt[, "f0", t]*oldPop[,"f0_diag_RNA"] + 
      tau_RNA_sc_dt[, "f0", t]*tau_ab_sc_dt[, "f0", t]*oldPop[, "f0_undiag"]+
      tau_poct_sc_dt[, "f0", t]*oldPop[, "f0_undiag"]-
      eta_sc_dt[, "f0", t]*(tau_RNA_sc_dt[, "f0", t]*tau_ab_sc_dt[, "f0", t]*oldPop[, "f0_undiag"]+ 
                              tau_poct_sc_dt[, "f0", t]*oldPop[, "f0_undiag"])     
    
    
    ## f0 treat
    newPop[, "f0_treat"] <- entry1[,"f0_treat"] + oldPop[,"f0_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_treat"])) - 
      death[,"f0_treat"] -leave[,"f0_treat"] - death_hcv[, "f0_treat"] +
      transition_dt[, "a_f0"]*oldPop[, "a_treat"] - 
      transition_dt[, "f0_f1"]*oldPop[, "f0_treat"] +  
      eta_dt[, "f0", t]*oldPop[,"f0_diag_RNA"] +  
      eta_sc_dt[, "f0", t]*(tau_RNA_sc_dt[, "f0", t]*tau_ab_sc_dt[, "f0", t]*oldPop[, "f0_undiag"]+ 
                              tau_poct_sc_dt[, "f0", t]*oldPop[, "f0_undiag"])   - 
      lota_dt[, "f0", t]*oldPop[,"f0_treat"] -  
      cure_dt[, "f0", t]*oldPop[,"f0_treat"] +  
      rho_dt[, "f0", t]*oldPop[,"f0_treat_f"]
    
    ## f0 treat failed
    newPop[, "f0_treat_f"] <- entry1[,"f0_treat_f"] + oldPop[,"f0_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_treat_f"])) -
      death[,"f0_treat_f"] -leave[,"f0_treat_f"] - death_hcv[, "f0_treat_f"] + 
      transition_dt[, "a_f0"]*oldPop[, "a_treat_f"] - 
      transition_dt[, "f0_f1"]*oldPop[, "f0_treat_f"] + 
      lota_dt[, "f0", t]*oldPop[,"f0_treat"] - 
      rho_dt[, "f0", t]*oldPop[,"f0_treat_f"]
    
    ## f0 cured
    newPop[, "f0_cured"] <- entry1[,"f0_cured"] + oldPop[,"f0_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "f0_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "f0_cured"]))  - 
      death[,"f0_cured"] -leave[,"f0_cured"] - death_hcv[, "f0_cured"] + 
      cure_dt[, "f0", t]*oldPop[, "f0_treat"] -
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf0[,]
      reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "f0_cured"]
    
    ##### f1 #####
    ## f1 infection
    newPop[, "f1_undiag"] <- entry1[,"f1_undiag"] + oldPop[,"f1_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_undiag"])) - 
      death[,"f1_undiag"] -leave[,"f1_undiag"] - death_hcv[, "f1_undiag"] +                                           
      transition_dt[, "f0_f1"]*oldPop[, "f0_undiag"] - 
      transition_dt[, "f1_f2"]*oldPop[, "f1_undiag"] - 
      tau_ab_dt[, "f1", t]*oldPop[, "f1_undiag"] - 
      tau_poct_dt[, "f1", t]*oldPop[, "f1_undiag"]  - 
      tau_ab_sc_dt[, "f1", t]*oldPop[, "f1_undiag"]- 
      tau_poct_sc_dt[, "f1", t]*oldPop[, "f1_undiag"]+
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf1[,]
      reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "f1_cured"]
    
    
    ## f1 testing, ab+
    newPop[, "f1_diag_ab"] <- entry1[,"f1_diag_ab"] + oldPop[ ,"f1_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_diag_ab"])) - 
      death[,"f1_diag_ab"] -leave[,"f1_diag_ab"] - death_hcv[, "f1_diag_ab"] + 
      transition_dt[, "f0_f1"]*oldPop[, "f0_diag_ab"] - 
      transition_dt[, "f1_f2"]*oldPop[, "f1_diag_ab"] + 
      tau_ab_dt[, "f1", t]*oldPop[, "f1_undiag"] - 
      tau_RNA_dt[, "f1", t]*oldPop[, "f1_diag_ab"] + 
      tau_ab_sc_dt[, "f1", t]*oldPop[, "f1_undiag"]- 
      tau_RNA_sc_dt[, "f1", t]*tau_ab_sc_dt[, "f1", t]*oldPop[, "f1_undiag"]    
    
    
    ## f1 testing, ag+/RNA
    newPop[, "f1_diag_RNA"] <-  entry1[,"f1_diag_RNA"] + oldPop[ ,"f1_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_diag_RNA"])) - 
      death[,"f1_diag_RNA"] -leave[,"f1_diag_RNA"] - death_hcv[, "f1_diag_RNA"] + 
      transition_dt[, "f0_f1"]*oldPop[, "f0_diag_RNA"] - 
      transition_dt[, "f1_f2"]*oldPop[, "f1_diag_RNA"] + 
      tau_RNA_dt[, "f1", t]*oldPop[, "f1_diag_ab"] + 
      tau_poct_dt[, "f1", t]*oldPop[, "f1_undiag"] - 
      eta_dt[, "f1", t]*oldPop[,"f1_diag_RNA"] + 
      tau_RNA_sc_dt[, "f1", t]*tau_ab_sc_dt[, "f1", t]*oldPop[, "f1_undiag"]  + 
      tau_poct_sc_dt[, "f1", t]*oldPop[, "f1_undiag"]- 
      eta_sc_dt[, "f1", t]*(tau_RNA_sc_dt[, "f1", t]*tau_ab_sc_dt[, "f1", t]*oldPop[, "f1_undiag"]  + 
                              tau_poct_sc_dt[, "f1", t]*oldPop[, "f1_undiag"])      
    
    
    ## f1 treat
    newPop[, "f1_treat"] <- entry1[,"f1_treat"] + oldPop[ ,"f1_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_treat"])) - 
      death[,"f1_treat"] -leave[,"f1_treat"] - death_hcv[, "f1_treat"] + 
      transition_dt[, "f0_f1"]*oldPop[, "f0_treat"] - 
      transition_dt[, "f1_f2"]*oldPop[, "f1_treat"] +  
      eta_dt[, "f1", t]*oldPop[,"f1_diag_RNA"] +  
      eta_sc_dt[, "f1", t]*(tau_RNA_sc_dt[, "f1", t]*tau_ab_sc_dt[, "f1", t]*oldPop[, "f1_undiag"]  + 
                              tau_poct_sc_dt[, "f1", t]*oldPop[, "f1_undiag"]) - 
      lota_dt[, "f1", t]*oldPop[,"f1_treat"] -  
      cure_dt[, "f1", t]*oldPop[,"f1_treat"] +  
      rho_dt[, "f1", t]*oldPop[,"f1_treat_f"]
    
    ## f1 treat failed
    newPop[, "f1_treat_f"] <- entry1[,"f1_treat_f"] + oldPop[,"f1_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_treat_f"])) - 
      death[,"f1_treat_f"] -leave[,"f1_treat_f"] - death_hcv[, "f1_treat_f"] + 
      transition_dt[, "f0_f1"]*oldPop[, "f0_treat_f"] - 
      transition_dt[, "f1_f2"]*oldPop[, "f1_treat_f"] + 
      lota_dt[, "f1", t]*oldPop[,"f1_treat"] - 
      rho_dt[, "f1", t]*oldPop[,"f1_treat_f"]
    
    ## f1 cured
    newPop[, "f1_cured"] <- entry1[,"f1_cured"] + oldPop[,"f1_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "f1_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "f1_cured"])) - 
      death[,"f1_cured"] -leave[,"f1_cured"] - death_hcv[, "f1_cured"] + 
      cure_dt[, "f1", t]*oldPop[, "f1_treat"] -
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf1[,] 
      reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "f1_cured"]
    
    ##### f2 #####
    ## f2 infection
    newPop[, "f2_undiag"] <- entry1[,"f2_undiag"] + oldPop[,"f2_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_undiag"])) - 
      death[,"f2_undiag"] -leave[,"f2_undiag"] - death_hcv[, "f2_undiag"] +                                                                                             
      transition_dt[, "f1_f2"]*oldPop[, "f1_undiag"] - 
      transition_dt[, "f2_f3"]*oldPop[, "f2_undiag"] - 
      tau_ab_dt[, "f2", t]*oldPop[, "f2_undiag"] -  
      tau_poct_dt[, "f2", t]*oldPop[, "f2_undiag"] - 
      tau_ab_sc_dt[, "f2", t]*oldPop[, "f2_undiag"]-  
      tau_poct_sc_dt[, "f2", t]*oldPop[, "f2_undiag"]+
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf2[,]
      reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "f2_cured"]
    
    ## f2 testing, ab+
    newPop[, "f2_diag_ab"] <- entry1[,"f2_diag_ab"] + oldPop[,"f2_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_diag_ab"])) - 
      death[,"f2_diag_ab"] -leave[,"f2_diag_ab"] - death_hcv[, "f2_diag_ab"] + 
      transition_dt[, "f1_f2"]*oldPop[, "f1_diag_ab"] - 
      transition_dt[, "f2_f3"]*oldPop[, "f2_diag_ab"] + 
      tau_ab_dt[, "f2", t]*oldPop[, "f2_undiag"] - 
      tau_RNA_dt[, "f2", t]*oldPop[, "f2_diag_ab"] + 
      tau_ab_sc_dt[, "f2", t]*oldPop[, "f2_undiag"]- 
      tau_RNA_sc_dt[, "f2", t]*tau_ab_sc_dt[, "f2", t]*oldPop[, "f2_undiag"]  
    
    
    ## f2 testing, ag+/RNA
    newPop[, "f2_diag_RNA"] <- entry1[,"f2_diag_RNA"] + oldPop[,"f2_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_diag_RNA"])) - 
      death[,"f2_diag_RNA"] -leave[,"f2_diag_RNA"] - death_hcv[, "f2_diag_RNA"] + 
      transition_dt[, "f1_f2"]*oldPop[, "f1_diag_RNA"] - 
      transition_dt[, "f2_f3"]*oldPop[, "f2_diag_RNA"] + 
      tau_RNA_dt[, "f2", t]*oldPop[, "f2_diag_ab"] + 
      tau_poct_dt[, "f2", t]*oldPop[, "f2_undiag"] - 
      eta_dt[, "f2", t]*oldPop[,"f2_diag_RNA"] + 
      tau_RNA_sc_dt[, "f2", t]*tau_ab_sc_dt[, "f2", t]*oldPop[, "f2_undiag"]+ 
      tau_poct_sc_dt[, "f2", t]*oldPop[, "f2_undiag"]- 
      eta_sc_dt[, "f2", t]*(tau_RNA_sc_dt[, "f2", t]*tau_ab_sc_dt[, "f2", t]*oldPop[, "f2_undiag"]+ 
                              tau_poct_sc_dt[, "f2", t]*oldPop[, "f2_undiag"])      
    
    
    ## f2 treat
    newPop[, "f2_treat"] <- entry1[,"f2_treat"] + oldPop[,"f2_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_treat"])) - 
      death[,"f2_treat"] -leave[,"f2_treat"] - death_hcv[, "f2_treat"] + 
      transition_dt[, "f1_f2"]*oldPop[, "f1_treat"] - 
      transition_dt[, "f2_f3"]*oldPop[, "f2_treat"] + 
      eta_dt[, "f2", t]*oldPop[,"f2_diag_RNA"] + 
      eta_sc_dt[, "f2", t]*(tau_RNA_sc_dt[, "f2", t]*tau_ab_sc_dt[, "f2", t]*oldPop[, "f2_undiag"]+ 
                              tau_poct_sc_dt[, "f2", t]*oldPop[, "f2_undiag"]) - 
      lota_dt[, "f2", t]*oldPop[,"f2_treat"] - 
      cure_dt[, "f2", t]*oldPop[,"f2_treat"] + 
      rho_dt[, "f2", t]*oldPop[,"f2_treat_f"]
    
    ## f2 treat failed
    newPop[, "f2_treat_f"] <- entry1[,"f2_treat_f"] + oldPop[,"f2_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_treat_f"]))  - 
      death[,"f2_treat_f"] - leave[,"f2_treat_f"] - death_hcv[, "f2_treat_f"] + 
      transition_dt[, "f1_f2"]*oldPop[, "f1_treat_f"] - 
      transition_dt[, "f2_f3"]*oldPop[, "f2_treat_f"] + 
      lota_dt[, "f2", t]*oldPop[,"f2_treat"] - 
      rho_dt[, "f2", t]*oldPop[,"f2_treat_f"]
    
    ## f2 cured
    newPop[, "f2_cured"] <- entry1[,"f2_cured"] + oldPop[,"f2_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "f2_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "f2_cured"])) - 
      death[,"f2_cured"] -leave[,"f2_cured"] - death_hcv[, "f2_cured"] + 
      cure_dt[, "f2", t]*oldPop[, "f2_treat"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf2[,] 
      reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "f2_cured"]
    
    ##### f3 #####
    ## f3 infection
    newPop[, "f3_undiag"] <- entry1[,"f3_undiag"] + oldPop[,"f3_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_undiag"])) - 
      death[,"f3_undiag"] -leave[,"f3_undiag"] - death_hcv[, "f3_undiag"]  + 
      transition_dt[, "f2_f3"]*oldPop[, "f2_undiag"] - 
      transition_dt[, "f3_f4"]*oldPop[, "f3_undiag"] - 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_undiag"] - 
      tau_ab_dt[, "f3", t]*oldPop[, "f3_undiag"] - 
      tau_poct_dt[, "f3", t]*oldPop[, "f3_undiag"] - 
      tau_ab_sc_dt[, "f3", t]*oldPop[, "f3_undiag"]- 
      tau_poct_sc_dt[, "f3", t]*oldPop[, "f3_undiag"]+ 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf3[,]
      reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "f3_cured"]
    
    
    ## f3 testing, ab+
    newPop[, "f3_diag_ab"] <- entry1[,"f3_diag_ab"] + oldPop[,"f3_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_diag_ab"]))  - 
      death[,"f3_diag_ab"] -leave[,"f3_diag_ab"] - death_hcv[, "f3_diag_ab"]  + 
      transition_dt[, "f2_f3"]*oldPop[, "f2_diag_ab"] - 
      transition_dt[, "f3_f4"]*oldPop[, "f3_diag_ab"] - 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_diag_ab"] + 
      tau_ab_dt[, "f3", t]*oldPop[, "f3_undiag"] - 
      tau_RNA_dt[, "f3", t]*oldPop[, "f3_diag_ab"] + 
      tau_ab_sc_dt[, "f3", t]*oldPop[, "f3_undiag"]- 
      tau_RNA_sc_dt[, "f3", t]*tau_ab_sc_dt[, "f3", t]*oldPop[, "f3_undiag"]    
    
    
    ## f3 testing, ag+/RNA
    newPop[, "f3_diag_RNA"] <- entry1[,"f3_diag_RNA"] + oldPop[,"f3_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_diag_RNA"])) - 
      death[,"f3_diag_RNA"] -leave[,"f3_diag_RNA"] - death_hcv[, "f3_diag_RNA"] +
      transition_dt[, "f2_f3"]*oldPop[, "f2_diag_RNA"] - 
      transition_dt[, "f3_f4"]*oldPop[, "f3_diag_RNA"] - 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_diag_RNA"] + 
      tau_RNA_dt[, "f3", t]*oldPop[, "f3_diag_ab"] + 
      tau_poct_dt[, "f3", t]*oldPop[, "f3_undiag"] - 
      eta_dt[, "f3", t]*oldPop[,"f3_diag_RNA"] + 
      tau_RNA_sc_dt[, "f3", t]* tau_ab_sc_dt[, "f3", t]*oldPop[, "f3_undiag"]+ 
      tau_poct_sc_dt[, "f3", t]*oldPop[, "f3_undiag"]- 
      eta_sc_dt[, "f3", t]*( tau_RNA_sc_dt[, "f3", t]* tau_ab_sc_dt[, "f3", t]*oldPop[, "f3_undiag"]+ 
                               tau_poct_sc_dt[, "f3", t]*oldPop[, "f3_undiag"])     
    
    
    ## f3 treat
    newPop[, "f3_treat"] <- entry1[,"f3_treat"] + oldPop[,"f3_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_treat"])) - 
      death[,"f3_treat"] -leave[,"f3_treat"] - death_hcv[, "f3_treat"] + 
      transition_dt[, "f2_f3"]*oldPop[, "f2_treat"] - 
      transition_dt[, "f3_f4"]*oldPop[, "f3_treat"] - 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_treat"] + 
      eta_dt[, "f3", t]*oldPop[,"f3_diag_RNA"] + 
      eta_sc_dt[, "f3", t]*(tau_RNA_sc_dt[, "f3", t]* tau_ab_sc_dt[, "f3", t]*oldPop[, "f3_undiag"]+ 
                               tau_poct_sc_dt[, "f3", t]*oldPop[, "f3_undiag"]) - 
      lota_dt[, "f3", t]*oldPop[,"f3_treat"] - 
      cure_dt[, "f3", t]*oldPop[,"f3_treat"] + 
      rho_dt[, "f3", t]*oldPop[,"f3_treat_f"]
    
    ## f3 treat failed
    newPop[, "f3_treat_f"] <- entry1[,"f3_treat_f"] + oldPop[,"f3_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_treat_f"])) - 
      death[,"f3_treat_f"] -leave[,"f3_treat_f"] - death_hcv[, "f3_treat_f"] + 
      transition_dt[, "f2_f3"]*oldPop[, "f2_treat_f"] - 
      transition_dt[, "f3_f4"]*oldPop[, "f3_treat_f"] - 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_treat_f"] + 
      lota_dt[, "f3", t]*oldPop[,"f3_treat"] - 
      rho_dt[, "f3", t]*oldPop[,"f3_treat_f"]
    
    ## f3 cured
    newPop[, "f3_cured"] <- entry1[,"f3_cured"] + oldPop[,"f3_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "f3_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "f3_cured"])) - 
      death[,"f3_cured"] -leave[,"f3_cured"] - death_hcv[, "f3_cured"] + 
      cure_dt[, "f3", t]*oldPop[, "f3_treat"]  -
      fibprog_dt[, "f3_cured_f4_cured"]*oldPop[, "f3_cured"] - 
      fibprog_dt[, "f3_cured_hcc_cured"]*oldPop[, "f3_cured"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf3[,]
      reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "f3_cured"]
    
    ##### f4 ##### 
    ## f4 infection
    newPop[, "f4_undiag"] <- entry1[,"f4_undiag"] + oldPop[,"f4_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_undiag"])) - 
      death[,"f4_undiag"] -leave[,"f4_undiag"] - death_hcv[, "f4_undiag"] + 
      transition_dt[, "f3_f4"]*oldPop[, "f3_undiag"] - 
      transition_dt[, "f4_dc"]*oldPop[, "f4_undiag"] - 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_undiag"] - 
      tau_ab_dt[, "f4", t]*oldPop[, "f4_undiag"] - 
      tau_poct_dt[, "f4", t]*oldPop[, "f4_undiag"] - 
      tau_ab_sc_dt[, "f4", t]*oldPop[, "f4_undiag"]- 
      tau_poct_sc_dt[, "f4", t]*oldPop[, "f4_undiag"]+
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf4[,]
      0*reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "f4_cured"]
    
    ## f4 testing, ab+
    newPop[, "f4_diag_ab"] <- entry1[,"f4_diag_ab"] + oldPop[,"f4_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_diag_ab"])) - 
      death[,"f4_diag_ab"] -leave[,"f4_diag_ab"] - death_hcv[, "f4_diag_ab"] + 
      transition_dt[, "f3_f4"]*oldPop[, "f3_diag_ab"] - 
      transition_dt[, "f4_dc"]*oldPop[, "f4_diag_ab"] - 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_diag_ab"] + 
      tau_ab_dt[, "f4", t]*oldPop[, "f4_undiag"] - 
      tau_RNA_dt[, "f4", t]*oldPop[, "f4_diag_ab"] + 
      tau_ab_sc_dt[, "f4", t]*oldPop[, "f4_undiag"]- 
      tau_RNA_sc_dt[, "f4", t]*tau_ab_sc_dt[, "f4", t]*oldPop[, "f4_undiag"]     
    
    
    ## f4 testing, ag+/RNA
    newPop[, "f4_diag_RNA"] <- entry1[,"f4_diag_RNA"] + oldPop[,"f4_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_diag_RNA"])) - 
      death[,"f4_diag_RNA"] - leave[,"f4_diag_RNA"] - death_hcv[, "f4_diag_RNA"] + 
      transition_dt[, "f3_f4"]*oldPop[, "f3_diag_RNA"] - 
      transition_dt[, "f4_dc"]*oldPop[, "f4_diag_RNA"] - 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_diag_RNA"] + 
      tau_RNA_dt[, "f4", t]*oldPop[, "f4_diag_ab"] + 
      tau_poct_dt[, "f4", t]*oldPop[, "f4_undiag"] - 
      eta_dt[, "f4", t]*oldPop[,"f4_diag_RNA"] + 
      tau_RNA_sc_dt[, "f4", t]*tau_ab_sc_dt[, "f4", t]*oldPop[, "f4_undiag"]+ 
      tau_poct_sc_dt[, "f4", t]*oldPop[, "f4_undiag"]- 
      eta_sc_dt[, "f4", t]*(tau_RNA_sc_dt[, "f4", t]*tau_ab_sc_dt[, "f4", t]*oldPop[, "f4_undiag"]+ 
                              tau_poct_sc_dt[, "f4", t]*oldPop[, "f4_undiag"])  
    
    
    ## f4 treat
    newPop[, "f4_treat"] <- entry1[,"f4_treat"] + oldPop[,"f4_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_treat"])) - 
      death[,"f4_treat"] -leave[,"f4_treat"] - death_hcv[, "f4_treat"] + 
      transition_dt[, "f3_f4"]*oldPop[, "f3_treat"] - 
      transition_dt[, "f4_dc"]*oldPop[, "f4_treat"] - 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_treat"] + 
      eta_dt[, "f4", t]*oldPop[,"f4_diag_RNA"] + 
      eta_sc_dt[, "f4", t]*(tau_RNA_sc_dt[, "f4", t]*tau_ab_sc_dt[, "f4", t]*oldPop[, "f4_undiag"]+ 
                              tau_poct_sc_dt[, "f4", t]*oldPop[, "f4_undiag"]) - 
      lota_dt[, "f4", t]*oldPop[,"f4_treat"] - 
      cure_dt[, "f4", t]*oldPop[,"f4_treat"] + 
      rho_dt[, "f4", t]*oldPop[,"f4_treat_f"]
    
    ## f4 treat failed
    newPop[, "f4_treat_f"] <- entry1[,"f4_treat_f"] + oldPop[,"f4_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_treat_f"])) - 
      death[,"f4_treat_f"] -leave[,"f4_treat_f"] - death_hcv[, "f4_treat_f"]  + 
      transition_dt[, "f3_f4"]*oldPop[, "f3_treat_f"] - 
      transition_dt[, "f4_dc"]*oldPop[, "f4_treat_f"] - 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_treat_f"] + 
      lota_dt[, "f4", t]*oldPop[,"f4_treat"] - 
      rho_dt[, "f4", t]*oldPop[,"f4_treat_f"]
    
    ## f4 cured
    newPop[, "f4_cured"] <- entry1[,"f4_cured"] + oldPop[,"f4_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "f4_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "f4_cured"])) - 
      death[,"f4_cured"] -leave[,"f4_cured"] - death_hcv[, "f4_cured"] + 
      cure_dt[, "f4", t]*oldPop[, "f4_treat"]  + 
      fibprog_dt[, "f3_cured_f4_cured"]*oldPop[, "f3_cured"] - 
      fibprog_dt[, "f4_cured_dc_cured"]*oldPop[, "f4_cured"] - 
      fibprog_dt[, "f4_cured_hcc_cured"]*oldPop[, "f4_cured"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcf4[,]
      0*reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "f4_cured"]
    
    
    
    
    ##### decompensated cirrhosis  (dc) #####
    ## DAA immediately (w/wo liver treatment)
    ## dc infection
    ## dc testing, ab+
    newPop[, "dc_undiag"] <- entry1[,"dc_undiag"] + oldPop[,"dc_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_undiag"])) - 
      death[,"dc_undiag"] -leave[,"dc_undiag"]  - 
      death_hcv[,"dc_undiag"]  + 
      transition_dt[, "f4_dc"]*oldPop[, "f4_undiag"] - 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_undiag"] - 
      transition_dt[, "dc_lt"]*oldPop[, "dc_undiag"] - 
      tau_ab_dt[, "dc", t]*oldPop[, "dc_undiag"] - 
      tau_poct_dt[, "dc", t]*oldPop[, "dc_undiag"] - 
      tau_ab_sc_dt[, "dc", t]*oldPop[, "dc_undiag"]- 
      tau_poct_sc_dt[, "dc", t]*oldPop[, "dc_undiag"]+
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcdc[,]
      0*reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "dc_cured"]
    
    newPop[, "dc_diag_ab"] <- entry1[,"dc_diag_ab"] + oldPop[,"dc_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_diag_ab"])) - 
      death[,"dc_diag_ab"] -leave[,"dc_diag_ab"] - 
      death_hcv[,"dc_diag_ab"] + 
      transition_dt[, "f4_dc"]*oldPop[, "f4_diag_ab"] - 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_diag_ab"] - 
      transition_dt[, "dc_lt"]*oldPop[, "dc_diag_ab"] + 
      tau_ab_dt[, "dc", t]*oldPop[, "dc_undiag"] - 
      tau_RNA_dt[, "dc", t]*oldPop[, "dc_diag_ab"] + 
      tau_ab_sc_dt[, "dc", t]*oldPop[, "dc_undiag"]- 
      tau_RNA_sc_dt[, "dc", t]*tau_ab_sc_dt[, "dc", t]*oldPop[, "dc_undiag"]     
    
    
    ## dc testing, ag+/RNA
    newPop[, "dc_diag_RNA"] <- entry1[,"dc_diag_RNA"] + oldPop[,"dc_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_diag_RNA"])) - 
      death[,"dc_diag_RNA"] -leave[,"dc_diag_RNA"] - 
      death_hcv[,"dc_diag_RNA"] + 
      transition_dt[, "f4_dc"]*oldPop[, "f4_diag_RNA"] - 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_diag_RNA"] - 
      transition_dt[, "dc_lt"]*oldPop[, "dc_diag_RNA"] + 
      tau_RNA_dt[, "dc", t]*oldPop[, "dc_diag_ab"] + 
      tau_poct_dt[, "dc", t]*oldPop[, "dc_undiag"] - 
      eta_dt[, "dc", t]*oldPop[,"dc_diag_RNA"] + 
      tau_RNA_sc_dt[, "dc", t]*tau_ab_sc_dt[, "dc", t]*oldPop[, "dc_undiag"]+ 
      tau_poct_sc_dt[, "dc", t]*oldPop[, "dc_undiag"]- 
      eta_sc_dt[, "dc", t]*(tau_RNA_sc_dt[, "dc", t]*tau_ab_sc_dt[, "dc", t]*oldPop[, "dc_undiag"]+ 
                              tau_poct_sc_dt[, "dc", t]*oldPop[, "dc_undiag"])       
    
    
    ## dc treat
    newPop[, "dc_treat"] <- entry1[,"dc_treat"] + oldPop[,"dc_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_treat"])) - 
      death[,"dc_treat"] -leave[,"dc_treat"] - 
      death_hcv[,"dc_treat"] + 
      transition_dt[, "f4_dc"]*oldPop[, "f4_treat"] - 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_treat"] - 
      transition_dt[, "dc_lt"]*oldPop[, "dc_treat"] + 
      eta_dt[, "dc", t]*oldPop[,"dc_diag_RNA"] + 
      eta_sc_dt[, "dc", t]*(tau_RNA_sc_dt[, "dc", t]*tau_ab_sc_dt[, "dc", t]*oldPop[, "dc_undiag"]+ 
      tau_poct_sc_dt[, "dc", t]*oldPop[, "dc_undiag"]) - 
      lota_dt[, "dc", t]*oldPop[,"dc_treat"] - 
      cure_dt[, "dc", t]*oldPop[,"dc_treat"] + 
      rho_dt[, "dc", t]*oldPop[,"dc_treat_f"]
    
    ## dc treat failed
    newPop[, "dc_treat_f"] <- entry1[,"dc_treat_f"] + oldPop[,"dc_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_treat_f"])) - 
      death[,"dc_treat_f"] -leave[,"dc_treat_f"] - 
      death_hcv[,"dc_treat_f"] + 
      transition_dt[, "f4_dc"]*oldPop[, "f4_treat_f"] - 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_treat_f"] - 
      transition_dt[, "dc_lt"]*oldPop[, "dc_treat_f"] + 
      lota_dt[, "dc", t]*oldPop[,"dc_treat"] - 
      rho_dt[, "dc", t]*oldPop[,"dc_treat_f"]
    
    ## dc cured
    newPop[, "dc_cured"] <- entry1[,"dc_cured"] + oldPop[,"dc_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "dc_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "dc_cured"])) - 
      death[,"dc_cured"] -leave[,"dc_cured"] - 
      death_hcv[,"dc_cured"] + 
      cure_dt[, "dc", t]*oldPop[, "dc_treat"]  - 
      fibprog_dt[, "dc_cured_lt_cured"]*oldPop[, "dc_cured"] - 
      fibprog_dt[, "dc_cured_hcc_cured"]*oldPop[, "dc_cured"] + 
      fibprog_dt[, "f4_cured_dc_cured"]*oldPop[, "f4_cured"] -
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcdc[,]
      0*reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "dc_cured"]
    
    
    
    
    ##### hep carcinoma (HCC) #####
    ## guideline recommend treat HCC first then giving DAA 
    
    ## HCC infection
    newPop[, "hcc_undiag"] <- entry1[,"hcc_undiag"] + oldPop[,"hcc_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_undiag"])) - 
      death[,"hcc_undiag"] -leave[,"hcc_undiag"] - 
      death_hcv[,"hcc_undiag"]  + 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_undiag"] + 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_undiag"] + 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_undiag"] - 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_undiag"] - 
      tau_ab_dt[, "hcc", t]*oldPop[, "hcc_undiag"] -  
      tau_poct_dt[, "hcc", t]*oldPop[, "hcc_undiag"] - 
      tau_ab_sc_dt[, "hcc", t]*oldPop[, "hcc_undiag"]-  
      tau_poct_sc_dt[, "hcc", t]*oldPop[, "hcc_undiag"]+ 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pchcc[,] 
      0*reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "hcc_cured"]
    
    
    ## hcc testing, ab+
    newPop[, "hcc_diag_ab"] <- entry1[,"hcc_diag_ab"] + oldPop[,"hcc_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_diag_ab"])) - 
      death[,"hcc_diag_ab"] -leave[,"hcc_diag_ab"] - 
      death_hcv[,"hcc_diag_ab"] + 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_diag_ab"] + 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_diag_ab"] + 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_diag_ab"] - 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_diag_ab"] + 
      tau_ab_dt[, "hcc", t]*oldPop[, "hcc_undiag"] - 
      tau_RNA_dt[, "hcc", t]*oldPop[, "hcc_diag_ab"] + 
      tau_ab_sc_dt[, "hcc", t]*oldPop[, "hcc_undiag"]- 
      tau_RNA_sc_dt[, "hcc", t]* tau_ab_sc_dt[, "hcc", t]*oldPop[, "hcc_undiag"]
    
    
    ## hcc testing, ag+/RNA
    newPop[, "hcc_diag_RNA"] <- entry1[,"hcc_diag_RNA"] + oldPop[,"hcc_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_diag_RNA"])) - 
      death[,"hcc_diag_RNA"] -leave[,"hcc_diag_RNA"] - 
      death_hcv[,"hcc_diag_RNA"] + 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_diag_RNA"] +
      transition_dt[, "dc_hcc"]*oldPop[, "dc_diag_RNA"] + 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_diag_RNA"] - 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_diag_RNA"] + 
      tau_RNA_dt[, "hcc", t]*oldPop[, "hcc_diag_ab"] +
      tau_poct_dt[, "hcc", t]*oldPop[, "hcc_undiag"] - 
      eta_dt[, "hcc", t]*oldPop[,"hcc_diag_RNA"] + 
      tau_RNA_sc_dt[, "hcc", t]* tau_ab_sc_dt[, "hcc", t]*oldPop[, "hcc_undiag"]+
      tau_poct_sc_dt[, "hcc", t]*oldPop[, "hcc_undiag"]- 
      eta_sc_dt[, "hcc", t]*(tau_RNA_sc_dt[, "hcc", t]* tau_ab_sc_dt[, "hcc", t]*oldPop[, "hcc_undiag"]+
                               tau_poct_sc_dt[, "hcc", t]*oldPop[, "hcc_undiag"])
    
    
    ## hcc treat
    newPop[, "hcc_treat"] <- entry1[,"hcc_treat"] + oldPop[,"hcc_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_treat"])) - 
      death[,"hcc_treat"] -leave[,"hcc_treat"] -
      death_hcv[,"hcc_treat"] + 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_treat"] + 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_treat"] + 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_treat"] - 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_treat"] + 
      eta_dt[, "hcc", t]*oldPop[,"hcc_diag_RNA"] + 
      eta_sc_dt[, "hcc", t]*(tau_RNA_sc_dt[, "hcc", t]* tau_ab_sc_dt[, "hcc", t]*oldPop[, "hcc_undiag"]+
                               tau_poct_sc_dt[, "hcc", t]*oldPop[, "hcc_undiag"]) - 
      lota_dt[, "hcc", t]*oldPop[,"hcc_treat"] -  
      cure_dt[, "hcc", t]*oldPop[,"hcc_treat"] +
      rho_dt[, "hcc", t]*oldPop[,"hcc_treat_f"]
    
    ## hcc treat failed
    newPop[, "hcc_treat_f"] <- entry1[,"hcc_treat_f"] + oldPop[,"hcc_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_treat_f"])) - 
      death[,"hcc_treat_f"] -leave[,"hcc_treat_f"] - 
      death_hcv[,"hcc_treat_f"]  + 
      transition_dt[, "f4_hcc"]*oldPop[, "f4_treat_f"] + 
      transition_dt[, "f3_hcc"]*oldPop[, "f3_treat_f"] + 
      transition_dt[, "dc_hcc"]*oldPop[, "dc_treat_f"] - 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_treat_f"] + 
      lota_dt[, "hcc", t]*oldPop[,"hcc_treat"] - 
      rho_dt[, "hcc", t]*oldPop[,"hcc_treat_f"]
    
    ## hcc cured
    newPop[, "hcc_cured"] <- entry1[,"hcc_cured"] + oldPop[,"hcc_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "hcc_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "hcc_cured"])) - 
      death[,"hcc_cured"] -leave[,"hcc_cured"] - 
      death_hcv[,"hcc_cured"] + 
      cure_dt[, "hcc", t]*oldPop[, "hcc_treat"]  + 
      fibprog_dt[, "f4_cured_hcc_cured"]*oldPop[, "f4_cured"] + 
      fibprog_dt[, "dc_cured_hcc_cured"]*oldPop[, "dc_cured"] +
      fibprog_dt[, "f3_cured_hcc_cured"]*oldPop[, "f3_cured"] - 
      fibprog_dt[, "hcc_cured_lt_cured"]*oldPop[, "hcc_cured"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pchcc[,] 
      0*reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "hcc_cured"]
    
    
    
    #####  liver failure needs to liver transplant ##### 
    ### LT infection
    newPop[, "lt_undiag"] <- entry1[,"lt_undiag"] + oldPop[,"lt_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_undiag"]))  - 
      death[,"lt_undiag"] -leave[,"lt_undiag"] - 
      death_hcv[,"lt_undiag"] + 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_undiag"] + 
      transition_dt[, "dc_lt"]*oldPop[, "dc_undiag"]  -
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_undiag"] - 
      tau_ab_dt[, "lt", t]*oldPop[, "lt_undiag"] - 
      tau_poct_dt[, "lt", t]*oldPop[, "lt_undiag"] - 
      tau_ab_sc_dt[, "lt", t]*oldPop[, "lt_undiag"]- 
      tau_poct_sc_dt[, "lt", t]*oldPop[, "lt_undiag"]+
      #reinfP[, t]*foi_dt[, ]*I[,]*Pclt[,] 
      0*reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "lt_cured"]
    
    ### LT testing, ab+
    newPop[, "lt_diag_ab"] <- entry1[,"lt_diag_ab"] + oldPop[,"lt_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_diag_ab"]))- 
      death[,"lt_diag_ab"] -leave[,"lt_diag_ab"] - 
      death_hcv[,"lt_diag_ab"] + 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_diag_ab"] + 
      transition_dt[, "dc_lt"]*oldPop[, "dc_diag_ab"]  - 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_diag_ab"] + 
      tau_ab_dt[, "lt", t]*oldPop[, "lt_undiag"] - 
      tau_RNA_dt[, "lt", t]*oldPop[, "lt_diag_ab"] + 
      tau_ab_sc_dt[, "lt", t]*oldPop[, "lt_undiag"]- 
      tau_RNA_sc_dt[, "lt", t]*tau_ab_sc_dt[, "lt", t]*oldPop[, "lt_undiag"]
    
    
    ### LT testing, ag+/RNA
    newPop[, "lt_diag_RNA"] <- entry1[,"lt_diag_RNA"] + oldPop[,"lt_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_diag_RNA"])) - 
      death[,"lt_diag_RNA"] -leave[,"lt_diag_RNA"] - 
      death_hcv[,"lt_diag_RNA"] + 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_diag_RNA"] + 
      transition_dt[, "dc_lt"]*oldPop[, "dc_diag_RNA"]  - 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_diag_RNA"] +
      tau_RNA_dt[, "lt", t]*oldPop[, "lt_diag_ab"] + 
      tau_poct_dt[, "lt", t]*oldPop[, "lt_undiag"] - 
      eta_dt[, "lt", t]*oldPop[,"lt_diag_RNA"] +
      tau_RNA_sc_dt[, "lt", t]* tau_ab_sc_dt[, "lt", t]*oldPop[, "lt_undiag"]+ 
      tau_poct_sc_dt[, "lt", t]*oldPop[, "lt_undiag"]- 
      eta_sc_dt[, "lt", t]*(tau_RNA_sc_dt[, "lt", t]* tau_ab_sc_dt[, "lt", t]*oldPop[, "lt_undiag"]+ 
                              tau_poct_sc_dt[, "lt", t]*oldPop[, "lt_undiag"])
    
    ### LT treat
    newPop[, "lt_treat"] <- entry1[,"lt_treat"] + oldPop[,"lt_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_treat"]))  -  
      death[,"lt_treat"] -leave[,"lt_treat"] - 
      death_hcv[,"lt_treat"] + 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_treat"] + 
      transition_dt[, "dc_lt"]*oldPop[, "dc_treat"] + 
      eta_dt[, "lt", t]*oldPop[,"lt_diag_RNA"] + 
      eta_sc_dt[, "lt", t]*(tau_RNA_sc_dt[, "lt", t]* tau_ab_sc_dt[, "lt", t]*oldPop[, "lt_undiag"]+ 
                              tau_poct_sc_dt[, "lt", t]*oldPop[, "lt_undiag"]) -  
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_treat"] - 
      lota_dt[, "lt", t]*oldPop[,"lt_treat"] - 
      cure_dt[, "lt", t]*oldPop[,"lt_treat"] + 
      rho_dt[, "lt", t]*oldPop[,"lt_treat_f"]
    
    ### LT treat failed
    newPop[, "lt_treat_f"] <- entry1[,"lt_treat_f"] + oldPop[,"lt_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_treat_f"])) - 
      death[,"lt_treat_f"] -leave[,"lt_treat_f"] - 
      death_hcv[,"lt_treat_f"]  + 
      transition_dt[, "hcc_lt"]*oldPop[, "hcc_treat_f"] + 
      transition_dt[, "dc_lt"]*oldPop[, "dc_treat_f"]  + 
      lota_dt[, "lt", t]*oldPop[,"lt_treat"] - 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_treat_f"] - 
      rho_dt[, "lt", t]*oldPop[,"lt_treat_f"]
    
    ### LT cured
    newPop[, "lt_cured"] <- entry1[,"lt_cured"] + oldPop[,"lt_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "lt_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "lt_cured"]))  - 
      death[,"lt_cured"] -leave[,"lt_cured"] - 
      death_hcv[,"lt_cured"] + 
      cure_dt[, "lt", t]*oldPop[, "lt_treat"]  + 
      fibprog_dt[, "hcc_cured_lt_cured"]*oldPop[, "hcc_cured"] + 
      fibprog_dt[, "dc_cured_lt_cured"]*oldPop[, "dc_cured"] - 
      fibprog_dt[, "lt_cured_plt_cured"]*oldPop[, "lt_cured"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pclt[,]
      0*reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "lt_cured"]
    
    ##### PLT #####
    ### PLT infection
    newPop[, "plt_undiag"] <- entry1[,"plt_undiag"] + oldPop[,"plt_undiag"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_undiag"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_undiag"])) - 
      death[,"plt_undiag"] -leave[,"plt_undiag"] - 
      death_hcv[,"plt_undiag"]  - 
      tau_ab_dt[, "plt", t]*oldPop[, "plt_undiag"] -
      tau_poct_dt[, "plt", t]*oldPop[, "plt_undiag"] - 
      tau_ab_sc_dt[, "plt", t]*oldPop[, "plt_undiag"]-
      tau_poct_sc_dt[, "plt", t]*oldPop[, "plt_undiag"]+ 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_undiag"] +
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcplt[,]
      0*reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "plt_cured"]
    
    ### PLT testing, ab+
    newPop[, "plt_diag_ab"] <- entry1[,"plt_diag_ab"] + oldPop[,"plt_diag_ab"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_diag_ab"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_diag_ab"])) - 
      death[,"plt_diag_ab"] -leave[,"plt_diag_ab"] - 
      death_hcv[,"plt_diag_ab"] + 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_diag_ab"] + 
      tau_ab_dt[, "plt", t]*oldPop[, "plt_undiag"] - 
      tau_RNA_dt[, "plt", t]*oldPop[, "plt_diag_ab"] + 
      tau_ab_sc_dt[, "plt", t]*oldPop[, "plt_undiag"]- 
      tau_RNA_sc_dt[, "plt", t]*tau_ab_sc_dt[, "plt", t]*oldPop[, "plt_undiag"]
    
    ### PLT testing, ag+/RNA
    newPop[, "plt_diag_RNA"] <- entry1[,"plt_diag_RNA"] + oldPop[,"plt_diag_RNA"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_diag_RNA"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_diag_RNA"])) - 
      death[,"plt_diag_RNA"] -leave[,"plt_diag_RNA"] - 
      death_hcv[,"plt_diag_RNA"]  + 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_diag_RNA"] + 
      tau_RNA_dt[, "plt", t]*oldPop[, "plt_diag_ab"] + 
      tau_poct_dt[, "plt", t]*oldPop[, "plt_undiag"] - 
      eta_dt[, "plt", t]*oldPop[,"plt_diag_RNA"] + 
      tau_RNA_sc_dt[, "plt", t]*tau_ab_sc_dt[, "plt", t]*oldPop[, "plt_undiag"]+ 
      tau_poct_sc_dt[, "plt", t]*oldPop[, "plt_undiag"]- 
      eta_sc_dt[, "plt", t]*(tau_RNA_sc_dt[, "plt", t]*tau_ab_sc_dt[, "plt", t]*oldPop[, "plt_undiag"]+ 
                               tau_poct_sc_dt[, "plt", t]*oldPop[, "plt_undiag"])
    
    ### PLT treat
    newPop[, "plt_treat"] <- entry1[,"plt_treat"] + oldPop[,"plt_treat"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_treat"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_treat"]))- 
      death[,"plt_treat"] -leave[,"plt_treat"] - 
      death_hcv[,"plt_treat"] + 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_treat"] + 
      eta_dt[, "plt", t]*oldPop[,"plt_diag_RNA"] + 
      eta_sc_dt[, "plt", t]*(tau_RNA_sc_dt[, "plt", t]*tau_ab_sc_dt[, "plt", t]*oldPop[, "plt_undiag"]+ 
                               tau_poct_sc_dt[, "plt", t]*oldPop[, "plt_undiag"]) - 
      lota_dt[, "plt", t]*oldPop[,"plt_treat"] - 
      cure_dt[, "plt", t]*oldPop[,"plt_treat"] +
      rho_dt[, "plt", t]*oldPop[,"plt_treat_f"]
    
    ### PLT treat failed
    newPop[, "plt_treat_f"] <- entry1[,"plt_treat_f"] + oldPop[,"plt_treat_f"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_treat_f"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_treat_f"])) - 
      death[,"plt_treat_f"] -leave[,"plt_treat_f"] - 
      death_hcv[,"plt_treat_f"] + 
      fibprog_dt[, "lt_plt"]*oldPop[, "lt_treat_f"] + 
      lota_dt[, "plt", t]*oldPop[,"plt_treat"] - 
      rho_dt[, "plt", t]*oldPop[,"plt_treat_f"]
    
    ### PLT cured
    newPop[, "plt_cured"] <- entry1[,"plt_cured"] + oldPop[,"plt_cured"] + 
      (colSums(pop_array[, , t]*oldPop[, "plt_cured"]) -
         rowSums(pop_array[, , t]*oldPop[, "plt_cured"])) - 
      death[,"plt_cured"] -leave[,"plt_cured"] - 
      death_hcv[,"plt_cured"] + 
      cure_dt[, "plt", t]*oldPop[, "plt_treat"] + 
      fibprog_dt[, "lt_cured_plt_cured"]*oldPop[, "lt_cured"] - 
      #reinfP[, t]*foi_dt[, ]*I[,]*Pcplt[,] 
      0*reinfP[, t]*foi_dt[, ]*II[,]*oldPop[, "plt_cured"]
    
    
    ifelse(newPop[,] < 0, 0, newPop[, ]) 
    
    
    
    #### uncertainty 
   
    
    #### result  aggregate ####  
    ##### population #####
    allPops[, , t] <- newPop

    ##### test S #####
    
    newS[, t] <- newPop[ , "s"]
    
    
    ##### new infection #####
    newInfections[, t] <-  foi_dt[, ]*II[,]*oldPop[, "s"]+ 
      reinfP[, t]*foi_dt[, ]*II[,]*(oldPop[, "a_cured"] + 
                                        oldPop[, "f0_cured"] + 
                                        oldPop[, "f1_cured"] + 
                                        oldPop[, "f2_cured"] + 
                                        oldPop[, "f3_cured"] + 
                                        0*oldPop[, "f4_cured"] + 
                                        0*oldPop[, "dc_cured"] + 
                                        0*oldPop[, "hcc_cured"] + 
                                        0*oldPop[, "lt_cured"] + 
                                        0*oldPop[, "plt_cured"])
    
    
    
    
    ##### background death ##### 
    
    for(i in 1: npops){ 
      newEntry[i,t] <- sum(entry1[i, ])
      newDeath[i, t] <- sum(death[i,])
      newDeathState[i, ,t] <- morb_dt[i,t]*oldPop[i,] 
      newLeave[i ,t] <- sum(leave[i, ])
      
      }

    ##### new HCV deaths in this timestep ##### 
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
    
    for(i in 1: npops){
      ##### new HCV death #####
      newHCVdeaths[i, t] <- sum( 0*oldPop[i, "s"],
                                 0*oldPop[i, "a_undiag"],
                                 0*oldPop[i, "a_diag_ab"],
                                 0*oldPop[i, "a_diag_RNA"],
                                 0*oldPop[i, "a_treat"],
                                 0*oldPop[i, "a_treat_f"],
                                 0*oldPop[i, "a_cured"],
                                 0*oldPop[i, "f0_undiag"],
                                 0*oldPop[i, "f0_diag_ab"],
                                 0*oldPop[i, "f0_diag_RNA"],
                                 0*oldPop[i, "f0_treat"],
                                 0*oldPop[i, "f0_treat_f"],
                                 0*oldPop[i, "f0_cured"],
                                 0*oldPop[i, "f1_undiag"],
                                 0*oldPop[i, "f1_diag_ab"],
                                 0*oldPop[i, "f1_diag_RNA"],
                                 0*oldPop[i, "f1_treat"],
                                 0*oldPop[i, "f1_treat_f"],
                                 0*oldPop[i, "f1_cured"],
                                 0*oldPop[i, "f2_undiag"],
                                 0*oldPop[i, "f2_diag_ab"],
                                 0*oldPop[i, "f2_diag_RNA"],
                                 0*oldPop[i, "f2_treat"],
                                 0*oldPop[i, "f2_treat_f"],
                                 0*oldPop[i, "f2_cured"],
                                 0*oldPop[i, "f3_undiag"],
                                 0*oldPop[i, "f3_diag_ab"],
                                 0*oldPop[i, "f3_diag_RNA"],
                                 0*oldPop[i, "f3_treat"],
                                 0*oldPop[i, "f3_treat_f"],
                                 0*oldPop[i, "f3_cured"],
                                 
                                 0*oldPop[i, "f4_undiag"],
                                 0*oldPop[i, "f4_diag_ab"],
                                 0*oldPop[i, "f4_diag_RNA"],
                                 0*oldPop[i, "f4_treat"],
                                 0*oldPop[i, "f4_treat_f"],
                                 0*oldPop[i, "f4_cured"],
                                mordc_dt[i, t]*oldPop[i, "dc_undiag"],
                                mordc_dt[i, t]*oldPop[i, "dc_diag_ab"],
                                mordc_dt[i, t]*oldPop[i, "dc_diag_RNA"],
                                mordc_dt[i,t]*oldPop[i, "dc_treat"],
                                mordc_dt[i,t]*oldPop[i, "dc_treat_f"],
                                mordcCure_dt[i,t]*oldPop[i, "dc_cured"],                      
                                morhcc_dt[i, t]*oldPop[i, "hcc_undiag"],                        
                                morhcc_dt[i, t]*oldPop[i, "hcc_diag_ab"],                       
                                morhcc_dt[i,t]*oldPop[i, "hcc_diag_RNA"],
                                morhcc_dt[i,t]*oldPop[i, "hcc_treat"],
                                morhcc_dt[i,t]*oldPop[i, "hcc_treat_f"],
                                morhccCure_dt[i,t]*oldPop[i, "hcc_cured"],
                                morlt_dt[i, t]*oldPop[i, "lt_undiag"],
                                morlt_dt[i, t]*oldPop[i, "lt_diag_ab"],
                                morlt_dt[i,t]*oldPop[i, "lt_diag_RNA"],
                                morlt_dt[i,t]*oldPop[i, "lt_treat"],
                                morlt_dt[i,t]*oldPop[i, "lt_treat_f"],
                                morlt_dt[i,t]*oldPop[i, "lt_cured"],
                                morplt_dt[i, t]*oldPop[i, "plt_undiag"],
                                morplt_dt[i, t]*oldPop[i, "plt_diag_ab"],
                                morplt_dt[i,t]*oldPop[i, "plt_diag_RNA"],
                                morplt_dt[i,t]*oldPop[i, "plt_treat"],
                                morplt_dt[i,t]*oldPop[i, "plt_treat_f"],
                                morplt_dt[i,t]*oldPop[i, "plt_cured"])
      
      ##### new treat #####
      
      newTreatment[i, t] <- sum(eta_dt[i, "a", t]*oldPop[i,"a_diag_RNA"],
                                eta_dt[i, "f0", t]*oldPop[i,"f0_diag_RNA"],
                                eta_dt[i, "f1", t]*oldPop[i,"f1_diag_RNA"],
                                eta_dt[i, "f2", t]*oldPop[i,"f2_diag_RNA"],
                                eta_dt[i, "f3", t]*oldPop[i,"f3_diag_RNA"],
                                eta_dt[i, "f4", t]*oldPop[i,"f4_diag_RNA"],
                                eta_dt[i, "dc", t]*oldPop[i,"dc_diag_RNA"],
                                eta_dt[i, "hcc", t]*oldPop[i,"hcc_diag_RNA"],
                                eta_dt[i, "lt", t]*oldPop[i,"lt_diag_RNA"],
                                eta_dt[i, "plt", t]*oldPop[i,"plt_diag_RNA"]) 
      
      ##### new retreat #####  
      
      newRetreat[i ,t] <- sum(rho_dt[i, "a", t]*oldPop[i,"a_treat_f"],
                              rho_dt[i, "f0", t]*oldPop[i,"f0_treat_f"],
                              rho_dt[i, "f1", t]*oldPop[i,"f1_treat_f"],
                              rho_dt[i, "f2", t]*oldPop[i,"f2_treat_f"],
                              rho_dt[i, "f3", t]*oldPop[i,"f3_treat_f"],
                              rho_dt[i, "f4", t]*oldPop[i,"f4_treat_f"],
                              rho_dt[i, "dc", t]*oldPop[i,"dc_treat_f"],
                              rho_dt[i, "hcc", t]*oldPop[i,"hcc_treat_f"],
                              rho_dt[i, "lt", t]*oldPop[i,"lt_treat_f"],
                              rho_dt[i, "plt", t]*oldPop[i,"plt_treat_f"])
      
      
      #### newtreatment failed ####
      
      newtreatfailed[i ,t] <- sum(lota_dt[i, "a", t]*oldPop[i,"a_treat"],
                                  lota_dt[i, "f0", t]*oldPop[i,"f0_treat"],
                                  lota_dt[i, "f1", t]*oldPop[i,"f1_treat"],
                                  lota_dt[i, "f2", t]*oldPop[i,"f2_treat"],
                                  lota_dt[i, "f3", t]*oldPop[i,"f3_treat"],
                                  lota_dt[i, "f4", t]*oldPop[i,"f4_treat"],
                                  lota_dt[i, "dc", t]*oldPop[i,"dc_treat"],
                                  lota_dt[i, "hcc", t]*oldPop[i,"hcc_treat"],
                                  lota_dt[i, "lt", t]*oldPop[i,"lt_treat"],
                                  lota_dt[i, "plt", t]*oldPop[i,"plt_treat"])
      
      
      ##### antibody test ##### 
      
      newTestingAb[i, t] <- sum( tau_ab_dt[i,"a", t]*(1-spc1_dt[i, t])*oldPop[i,"a_undiag"],
                                  tau_ab_dt[i, "f0", t]*oldPop[i,"f0_undiag"],
                                  tau_ab_dt[i, "f1", t]*oldPop[i,"f1_undiag"],
                                  tau_ab_dt[i, "f2", t]*oldPop[i,"f2_undiag"],
                                  tau_ab_dt[i, "f3", t]*oldPop[i,"f3_undiag"],
                                  tau_ab_dt[i, "f4", t]*oldPop[i,"f4_undiag"],
                                  tau_ab_dt[i, "dc", t]*oldPop[i,"dc_undiag"],
                                  tau_ab_dt[i, "hcc", t]*oldPop[i,"hcc_undiag"],
                                  tau_ab_dt[i, "lt", t]*oldPop[i,"lt_undiag"],
                                  tau_ab_dt[i, "plt", t]*oldPop[i,"plt_undiag"])
      
      newTestingAb_neg[i, t] <- 0.1*tau_ab_dt[i, "f0", t]*(oldPop[i,"s"] + oldPop[i,"a_cured"] )
                                 
      
      ##### antigen test ##### 
      
      newTestingAg[i, t] <- sum(tau_RNA_dt[i, "a", t]*oldPop[i,"a_diag_ab"],
                                tau_RNA_dt[i, "f0", t]*oldPop[i,"f0_diag_ab"],
                                tau_RNA_dt[i, "f1", t]*oldPop[i,"f1_diag_ab"],
                                tau_RNA_dt[i, "f2", t]*oldPop[i,"f2_diag_ab"],
                                tau_RNA_dt[i, "f3", t]*oldPop[i,"f3_diag_ab"],
                                tau_RNA_dt[i, "f4", t]*oldPop[i,"f4_diag_ab"],
                                tau_RNA_dt[i, "dc", t]*oldPop[i,"dc_diag_ab"],
                                tau_RNA_dt[i, "hcc", t]*oldPop[i,"hcc_diag_ab"],
                                tau_RNA_dt[i, "lt", t]*oldPop[i,"lt_diag_ab"],
                                tau_RNA_dt[i, "plt", t]*oldPop[i,"plt_diag_ab"])
      
      newTestingAg_neg[i, t] <- 0.1*tau_RNA_dt[i,"f0",t]*(oldPop[i , "f0_cured"] +
                                                        oldPop[i , "f1_cured"] +
                                                        oldPop[i , "f2_cured"] +
                                                        oldPop[i , "f3_cured"] +
                                                        oldPop[i , "f4_cured"] +
                                                        oldPop[i , "dc_cured"] +
                                                        oldPop[i , "hcc_cured"] +
                                                        oldPop[i , "lt_cured"] +
                                                        oldPop[i , "plt_cured"])
      
      ##### POCT test #####
      
      newTestingPOCT[i, t] <- sum( tau_poct_dt[i, "a", t]*(1-spc1_dt[i, t])*oldPop[i,"a_undiag"],
                                   tau_poct_dt[i, "f0", t]*oldPop[i,"f0_undiag"],
                                   tau_poct_dt[i, "f1", t]*oldPop[i,"f1_undiag"],
                                   tau_poct_dt[i, "f2", t]*oldPop[i,"f2_undiag"],
                                   tau_poct_dt[i, "f3", t]*oldPop[i,"f3_undiag"],
                                   tau_poct_dt[i, "f4", t]*oldPop[i,"f4_undiag"],
                                   tau_poct_dt[i, "dc", t]*oldPop[i,"dc_undiag"],
                                   tau_poct_dt[i, "hcc", t]*oldPop[i,"hcc_undiag"],
                                   tau_poct_dt[i, "lt", t]*oldPop[i,"lt_undiag"],
                                   tau_poct_dt[i, "plt", t]*oldPop[i,"plt_undiag"])
      
      
      ##### Cured ##### 
      newCured[i, t] <- sum( cure_dt[i, "a", t]*oldPop[i,"a_treat"],
                             cure_dt[i, "f0", t]*oldPop[i,"f0_treat"],
                             cure_dt[i, "f1", t]*oldPop[i,"f1_treat"],
                             cure_dt[i, "f2", t]*oldPop[i,"f2_treat"],
                             cure_dt[i, "f3", t]*oldPop[i,"f3_treat"],
                             cure_dt[i, "f4", t]*oldPop[i,"f4_treat"],
                             cure_dt[i, "dc", t]*oldPop[i,"dc_treat"],
                             cure_dt[i, "hcc", t]*oldPop[i,"hcc_treat"],
                             cure_dt[i, "lt", t]*oldPop[i,"lt_treat"],
                             cure_dt[i, "plt", t]*oldPop[i,"plt_treat"])
      
      
      #### intervention  ####
      
      ##### new treat #####
      
      newTreatment_sc[i, t] <- 
        sum(eta_sc_dt[i, "a", t]*oldPop[i, "a_undiag"]*
              (tau_RNA_sc_dt[i, "a", t]*tau_ab_sc_dt[i, "a", t]*(1-spc1_dt[i, t])+ tau_poct_sc_dt[i, "a", t]*(1-spc1_dt[i, t])), 
            eta_sc_dt[i, "f0", t]*oldPop[i, "f0_undiag"]*
              (tau_RNA_sc_dt[i, "f0", t]*tau_ab_sc_dt[i, "f0", t] + tau_poct_sc_dt[i, "f0", t]),
            eta_sc_dt[i, "f1", t]*oldPop[i, "f1_undiag"]*
              (tau_RNA_sc_dt[i, "f1", t]*tau_ab_sc_dt[i, "f1", t] + tau_poct_sc_dt[i, "f1", t]),
            eta_sc_dt[i, "f2", t]*oldPop[i, "f2_undiag"]*
              (tau_RNA_sc_dt[i, "f2", t]*tau_ab_sc_dt[i, "f2", t] + tau_poct_sc_dt[i, "f2", t]),
            eta_sc_dt[i, "f3", t]*oldPop[i, "f3_undiag"]*
              (tau_RNA_sc_dt[i, "f3", t]*tau_ab_sc_dt[i, "f3", t] + tau_poct_sc_dt[i, "f3", t]),
            eta_sc_dt[i, "f4", t]*oldPop[i, "f4_undiag"]*
              (tau_RNA_sc_dt[i, "f4", t]*tau_ab_sc_dt[i, "f4", t] + tau_poct_sc_dt[i, "f4", t]),
            eta_sc_dt[i, "dc", t]*oldPop[i, "dc_undiag"]*
              (tau_RNA_sc_dt[, "dc", t]*tau_ab_sc_dt[i, "dc", t] + tau_poct_sc_dt[i, "dc", t]),
            eta_sc_dt[i, "hcc", t]*oldPop[i, "hcc_undiag"]*
              (tau_RNA_sc_dt[i, "hcc", t]*tau_ab_sc_dt[i, "hcc", t] + tau_poct_sc_dt[i, "hcc", t]),
            eta_sc_dt[i, "lt", t]*oldPop[i, "lt_undiag"]*
              (tau_RNA_sc_dt[i, "lt", t]*tau_ab_sc_dt[i, "lt", t] + tau_poct_sc_dt[i, "lt", t]),
            eta_sc_dt[i, "plt", t]*oldPop[i, "plt_undiag"]*
              (tau_RNA_sc_dt[i, "plt", t]*tau_ab_sc_dt[i, "plt", t] + tau_poct_sc_dt[i, "plt", t])
            
            ) 
      
      
      ##### antibody test ##### 
      
      newTestingAb_sc[i, t] <- sum(tau_ab_sc_dt[i,"a", t]*(1-spc1_dt[i, t])*oldPop[i,"a_undiag"],
                                 tau_ab_sc_dt[i, "f0", t]*oldPop[i,"f0_undiag"],
                                 tau_ab_sc_dt[i, "f1", t]*oldPop[i,"f1_undiag"],
                                 tau_ab_sc_dt[i, "f2", t]*oldPop[i,"f2_undiag"],
                                 tau_ab_sc_dt[i, "f3", t]*oldPop[i,"f3_undiag"],
                                 tau_ab_sc_dt[i, "f4", t]*oldPop[i,"f4_undiag"],
                                 tau_ab_sc_dt[i, "dc", t]*oldPop[i,"dc_undiag"],
                                 tau_ab_sc_dt[i, "hcc", t]*oldPop[i,"hcc_undiag"],
                                 tau_ab_sc_dt[i, "lt", t]*oldPop[i,"lt_undiag"],
                                 tau_ab_sc_dt[i, "plt", t]*oldPop[i,"plt_undiag"])
      
      newTestingAb_sc_neg[i, t] <- fc[i, t]*tau_ab_sc_dt[i,"f0", t]*(oldPop[i , "f0_cured"] +
                                            oldPop[i , "f1_cured"] +
                                            oldPop[i , "f2_cured"] +
                                            oldPop[i , "f3_cured"] +
                                            oldPop[i , "f4_cured"] +
                                            oldPop[i , "dc_cured"] +
                                            oldPop[i , "hcc_cured"] +
                                            oldPop[i , "lt_cured"] +
                                            oldPop[i , "plt_cured"] + 
                                            oldPop[i,"s"]) + 
        fc[i,t]*tau_ab_sc_dt[i,"a", t]*((oldPop[i,"a_cured"])) 
      
      ##### antigen test ##### 
      
      newTestingAg_sc[i, t] <- sum(tau_RNA_sc_dt[i, "a", t]*tau_ab_sc_dt[i, "a", t]*(1-spc1_dt[i, t])*oldPop[i, "a_undiag"],
                                tau_RNA_sc_dt[i, "f0", t]*tau_ab_sc_dt[i, "f0", t]*oldPop[i,"f0_undiag"],
                                tau_RNA_sc_dt[i, "f1", t]*tau_ab_sc_dt[i, "f1", t]*oldPop[i,"f1_undiag"],
                                tau_RNA_sc_dt[i, "f2", t]*tau_ab_sc_dt[i, "f2", t]*oldPop[i,"f2_undiag"],
                                tau_RNA_sc_dt[i, "f3", t]*tau_ab_sc_dt[i, "f3", t]*oldPop[i,"f3_undiag"],
                                tau_RNA_sc_dt[i, "f4", t]*tau_ab_sc_dt[i, "f4", t]*oldPop[i,"f4_undiag"],
                                tau_RNA_sc_dt[i, "dc", t]*tau_ab_sc_dt[i, "dc", t]*oldPop[i,"dc_undiag"],
                                tau_RNA_sc_dt[i, "hcc", t]*tau_ab_sc_dt[i, "hcc", t]*oldPop[i,"hcc_undiag"],
                                tau_RNA_sc_dt[i, "lt", t]*tau_ab_sc_dt[i, "lt", t]*oldPop[i,"lt_undiag"],
                                tau_RNA_sc_dt[i, "plt", t]*tau_ab_sc_dt[i, "plt", t]*oldPop[i,"plt_undiag"])
      
      newTestingAg_sc_neg[i,t] <- tau_RNA_sc_dt[i, "f0", t]*fc[i, t]*tau_ab_sc_dt[i,"f0", t]*(oldPop[i , "f0_cured"] +
                                                                                                        oldPop[i , "f1_cured"] +
                                                                                                        oldPop[i , "f2_cured"] +
                                                                                                        oldPop[i , "f3_cured"] +
                                                                                                        oldPop[i , "f4_cured"] +
                                                                                                        oldPop[i , "dc_cured"] +
                                                                                                        oldPop[i , "hcc_cured"] +
                                                                                                        oldPop[i , "lt_cured"] +
                                                                                                        oldPop[i , "plt_cured"]) + 
        
        tau_RNA_sc_dt[i, "a", t]*fc[i,t]*tau_ab_sc_dt[i,"a", t]*((oldPop[i,"a_cured"])) 
      
      ##### POCT test #####
      
      newTestingPOCT_sc[i, t] <- sum( tau_poct_sc_dt[i, "a", t]*(1-spc1_dt[i, t])*oldPop[i,"a_undiag"],
                                   tau_poct_sc_dt[i, "f0", t]*oldPop[i,"f0_undiag"],
                                   tau_poct_sc_dt[i, "f1", t]*oldPop[i,"f1_undiag"],
                                   tau_poct_sc_dt[i, "f2", t]*oldPop[i,"f2_undiag"],
                                   tau_poct_sc_dt[i, "f3", t]*oldPop[i,"f3_undiag"],
                                   tau_poct_sc_dt[i, "f4", t]*oldPop[i,"f4_undiag"],
                                   tau_poct_sc_dt[i, "dc", t]*oldPop[i,"dc_undiag"],
                                   tau_poct_sc_dt[i, "hcc", t]*oldPop[i,"hcc_undiag"],
                                   tau_poct_sc_dt[i, "lt", t]*oldPop[i,"lt_undiag"],
                                   tau_poct_sc_dt[i, "plt", t]*oldPop[i,"plt_undiag"])
      
      
      
      
    }
    # POCT negative 
    newTestingPOCT_neg[1, t] <- tau_poct_dt[1, "f0", t]*(oldPop[1 , "f0_cured"] +
                                                             oldPop[1 , "f1_cured"] +
                                                             oldPop[1 , "f2_cured"] +
                                                             oldPop[1 , "f3_cured"] +
                                                             oldPop[1 , "f4_cured"] +
                                                             oldPop[1 , "dc_cured"] +
                                                             oldPop[1 , "hcc_cured"] +
                                                             oldPop[1 , "lt_cured"] +
                                                             oldPop[1 , "plt_cured"]+
                                                           oldPop[1,"s"] ) + 
        tau_poct_dt[1, "a", t]*((oldPop[1,"a_cured"] )) 
    
    newTestingPOCT_neg[2, t] <- tau_poct_dt[2, "f0", t]*(oldPop[2 , "f0_cured"] +
                                                           oldPop[2 , "f1_cured"] +
                                                           oldPop[2 , "f2_cured"] +
                                                           oldPop[2 , "f3_cured"] +
                                                           oldPop[2 , "f4_cured"] +
                                                           oldPop[2 , "dc_cured"] +
                                                           oldPop[2 , "hcc_cured"] +
                                                           oldPop[2 , "lt_cured"] +
                                                           oldPop[2 , "plt_cured"] + 
                                                           oldPop[2,"s"] ) + 
      tau_poct_dt[2, "a", t]*((oldPop[2,"a_cured"] )) 
    
    newTestingPOCT_neg[3, t] <- tau_poct_dt[3, "f0", t]*(oldPop[3 , "f0_cured"] +
                                                           oldPop[3 , "f1_cured"] +
                                                           oldPop[3 , "f2_cured"] +
                                                           oldPop[3 , "f3_cured"] +
                                                           oldPop[3 , "f4_cured"] +
                                                           oldPop[3 , "dc_cured"] +
                                                           oldPop[3 , "hcc_cured"] +
                                                           oldPop[3 , "lt_cured"] +
                                                           oldPop[3 , "plt_cured"]+
                                                           oldPop[3,"s"]) + 
      tau_poct_dt[3, "a", t]*(( oldPop[3,"a_cured"] ))
    
    newTestingPOCT_neg[4, t] <- tau_poct_dt[4, "f0", t]*(oldPop[4 , "f0_cured"] +
                                                           oldPop[4 , "f1_cured"] +
                                                           oldPop[4 , "f2_cured"] +
                                                           oldPop[4 , "f3_cured"] +
                                                           oldPop[4 , "f4_cured"] +
                                                           oldPop[4 , "dc_cured"] +
                                                           oldPop[4 , "hcc_cured"] +
                                                           oldPop[4 , "lt_cured"] +
                                                           oldPop[4 , "plt_cured"] + 
                                                           oldPop[4,"s"]) + 
      tau_poct_dt[4, "a", t]*(( oldPop[4,"a_cured"] )) 
    
    
      newTestingPOCT_neg[5, t] <- tau_poct_dt[5, "f0", t]*(oldPop[5 , "f0_cured"] +
                                                               oldPop[5 , "f1_cured"] +
                                                               oldPop[5 , "f2_cured"] +
                                                               oldPop[5 , "f3_cured"] +
                                                               oldPop[5 , "f4_cured"] +
                                                               oldPop[5 , "dc_cured"] +
                                                               oldPop[5 , "hcc_cured"] +
                                                               oldPop[5 , "lt_cured"] +
                                                               oldPop[5 , "plt_cured"] + 
                                                             oldPop[5,"s"]) + 
       tau_poct_dt[5, "a", t]*((oldPop[5,"a_cured"] )) 
      
      # POCT negative scenario
      newTestingPOCT_sc_neg[1, t] <- fc[1, t]*tau_poct_sc_dt[1, "f0", t]*(oldPop[1 , "f0_cured"] +
                                                             oldPop[1 , "f1_cured"] +
                                                             oldPop[1 , "f2_cured"] +
                                                             oldPop[1 , "f3_cured"] +
                                                             oldPop[1 , "f4_cured"] +
                                                             oldPop[1 , "dc_cured"] +
                                                             oldPop[1 , "hcc_cured"] +
                                                             oldPop[1 , "lt_cured"] +
                                                             oldPop[1 , "plt_cured"] +
                                                             oldPop[1 , "s"]) + 
        fc[1,t]*tau_poct_sc_dt[1, "a", t]*((oldPop[1,"a_cured"]   ))
      
      newTestingPOCT_sc_neg[2, t] <- fc[2, t]*tau_poct_sc_dt[2, "f0", t]*(oldPop[2 , "f0_cured"] +
                                                             oldPop[2 , "f1_cured"] +
                                                             oldPop[2 , "f2_cured"] +
                                                             oldPop[2 , "f3_cured"] +
                                                             oldPop[2 , "f4_cured"] +
                                                             oldPop[2 , "dc_cured"] +
                                                             oldPop[2 , "hcc_cured"] +
                                                             oldPop[2 , "lt_cured"] +
                                                             oldPop[2 , "plt_cured"] + 
                                                               oldPop[2,"s"]) + 
        fc[2,t]*tau_poct_sc_dt[2, "a", t]*((oldPop[2,"a_cured"]   ))
      
      newTestingPOCT_sc_neg[3, t] <- fc[3,t]*tau_poct_sc_dt[3, "f0", t]*(oldPop[3 , "f0_cured"] +
                                                             oldPop[3 , "f1_cured"] +
                                                             oldPop[3 , "f2_cured"] +
                                                             oldPop[3 , "f3_cured"] +
                                                             oldPop[3 , "f4_cured"] +
                                                             oldPop[3 , "dc_cured"] +
                                                             oldPop[3 , "hcc_cured"] +
                                                             oldPop[3 , "lt_cured"] +
                                                             oldPop[3 , "plt_cured"] + 
                                                               oldPop[3,"s"]) + 
        fc[3,t]*tau_poct_sc_dt[3, "a", t]*((oldPop[3,"a_cured"] ))
      
      newTestingPOCT_sc_neg[4, t] <- fc[4,t]*tau_poct_sc_dt[4, "f0", t]*(oldPop[4 , "f0_cured"] +
                                                             oldPop[4 , "f1_cured"] +
                                                             oldPop[4 , "f2_cured"] +
                                                             oldPop[4 , "f3_cured"] +
                                                             oldPop[4 , "f4_cured"] +
                                                             oldPop[4 , "dc_cured"] +
                                                             oldPop[4 , "hcc_cured"] +
                                                             oldPop[4 , "lt_cured"] +
                                                             oldPop[4 , "plt_cured"] + 
                                                               oldPop[4,"s"]) + 
        fc[4,t]*tau_poct_sc_dt[4, "a", t]*((oldPop[4,"a_cured"] )) 
      
      
      newTestingPOCT_sc_neg[5, t] <- fc[5,t]*tau_poct_sc_dt[5, "f0", t]*(oldPop[5 , "f0_cured"] +
                                                               oldPop[5 , "f1_cured"] +
                                                               oldPop[5 , "f2_cured"] +
                                                               oldPop[5 , "f3_cured"] +
                                                               oldPop[5 , "f4_cured"] +
                                                               oldPop[5 , "dc_cured"] +
                                                               oldPop[5 , "hcc_cured"] +
                                                               oldPop[5 , "lt_cured"] +
                                                               oldPop[5 , "plt_cured"]+ 
                                                                 oldPop[5,"s"]) + 
        fc[5,t]*tau_poct_sc_dt[5, "a", t]*((oldPop[5,"a_cured"] ))   
    
    
    
    #####reinfection##### 
    
    
    newreinfection[, t] <- reinfP[, t]*foi_dt[, ]*II[,]*(oldPop[, "a_cured"] + 
                                                             oldPop[, "f0_cured"] + 
                                                             oldPop[, "f1_cured"] + 
                                                             oldPop[, "f2_cured"] + 
                                                             oldPop[, "f3_cured"] + 
                                                             0*oldPop[, "f4_cured"] + 
                                                             0*oldPop[, "dc_cured"] + 
                                                             0*oldPop[, "hcc_cured"] + 
                                                             0*oldPop[, "lt_cured"] + 
                                                             0*oldPop[, "plt_cured"])
      
      newreinfection_chronic[, t] <- reinfP[, t]*foi_dt[, ]*II[,]*(
                                                             oldPop[, "f0_cured"] + 
                                                             oldPop[, "f1_cured"] + 
                                                             oldPop[, "f2_cured"] + 
                                                             oldPop[, "f3_cured"] + 
                                                             0*oldPop[, "f4_cured"] + 
                                                             0*oldPop[, "dc_cured"] + 
                                                             0*oldPop[, "hcc_cured"] + 
                                                             0*oldPop[, "lt_cured"] + 
                                                             0*oldPop[, "plt_cured"])
    #### in and out flow #### 
    # each state in and out flow 
    
    newarray_state["s", , , t] <- pop_array[, , t]*oldPop[,"s"]
    newarray_state["a_undiag", , , t] <- pop_array[, , t]*oldPop[,"a_undiag"]
    newarray_state["f0_undiag", , , t] <- pop_array[, , t]*oldPop[,"f0_undiag"]
    newarray_state["f1_undiag", , , t] <- pop_array[, , t]*oldPop[,"f1_undiag"]
    newarray_state["f2_undiag", , , t] <- pop_array[, , t]*oldPop[,"f2_undiag"]
    newarray_state["f3_undiag", , , t] <- pop_array[, , t]*oldPop[,"f3_undiag"]
    newarray_state["f4_undiag", , , t] <- pop_array[, , t]*oldPop[,"f4_undiag"]
    newarray_state["dc_undiag", , , t] <- pop_array[, , t]*oldPop[,"dc_undiag"]
    newarray_state["hcc_undiag", , , t] <- pop_array[, , t]*oldPop[,"hcc_undiag"]
    newarray_state["lt_undiag", , , t] <- pop_array[, , t]*oldPop[,"lt_undiag"]
    newarray_state["plt_undiag", , , t] <- pop_array[, , t]*oldPop[,"plt_undiag"]
    
    newarray_state["a_diag_ab", , , t] <- pop_array[, , t]*oldPop[,"a_diag_ab"]
    newarray_state["f0_diag_ab", , , t] <- pop_array[, , t]*oldPop[,"f0_diag_ab"]
    newarray_state["f1_diag_ab", , , t] <- pop_array[, , t]*oldPop[,"f1_diag_ab"]
    newarray_state["f2_diag_ab", , , t] <- pop_array[, , t]*oldPop[,"f2_diag_ab"]
    newarray_state["f3_diag_ab", , , t] <- pop_array[, , t]*oldPop[,"f3_diag_ab"]
    newarray_state["f4_diag_ab", , , t] <- pop_array[, , t]*oldPop[,"f4_diag_ab"]
    newarray_state["dc_diag_ab", , , t] <- pop_array[, , t]*oldPop[,"dc_diag_ab"]
    newarray_state["hcc_diag_ab", , , t] <- pop_array[, , t]*oldPop[,"hcc_diag_ab"]
    newarray_state["lt_diag_ab", , , t] <- pop_array[, , t]*oldPop[,"lt_diag_ab"]
    newarray_state["plt_diag_ab", , , t] <- pop_array[, , t]*oldPop[,"plt_diag_ab"]
    
    newarray_state["a_diag_RNA", , , t] <- pop_array[, , t]*oldPop[,"a_diag_RNA"]
    newarray_state["f0_diag_RNA", , , t] <- pop_array[, , t]*oldPop[,"f0_diag_RNA"]
    newarray_state["f1_diag_RNA", , , t] <- pop_array[, , t]*oldPop[,"f1_diag_RNA"]
    newarray_state["f2_diag_RNA", , , t] <- pop_array[, , t]*oldPop[,"f2_diag_RNA"]
    newarray_state["f3_diag_RNA", , , t] <- pop_array[, , t]*oldPop[,"f3_diag_RNA"]
    newarray_state["f4_diag_RNA", , , t] <- pop_array[, , t]*oldPop[,"f4_diag_RNA"]
    newarray_state["dc_diag_RNA", , , t] <- pop_array[, , t]*oldPop[,"dc_diag_RNA"]
    newarray_state["hcc_diag_RNA", , , t] <- pop_array[, , t]*oldPop[,"hcc_diag_RNA"]
    newarray_state["lt_diag_RNA", , , t] <- pop_array[, , t]*oldPop[,"lt_diag_RNA"]
    newarray_state["plt_diag_RNA", , , t] <- pop_array[, , t]*oldPop[,"plt_diag_RNA"]
    
    newarray_state["a_treat", , , t] <- pop_array[, , t]*oldPop[,"a_treat"]
    newarray_state["f0_treat" , , , t] <- pop_array[, , t]*oldPop[,"f0_treat"]
    newarray_state["f1_treat" , , , t] <- pop_array[, , t]*oldPop[,"f1_treat"]
    newarray_state["f2_treat", , , t] <- pop_array[, , t]*oldPop[,"f2_treat"]
    newarray_state["f3_treat", , , t] <- pop_array[, , t]*oldPop[,"f3_treat"]
    newarray_state["f4_treat", , , t] <- pop_array[, , t]*oldPop[,"f4_treat"]
    newarray_state["dc_treat", , , t] <- pop_array[, , t]*oldPop[,"dc_treat"]
    newarray_state["hcc_treat", , , t] <- pop_array[, , t]*oldPop[,"hcc_treat"]
    newarray_state["lt_treat", , , t] <- pop_array[, , t]*oldPop[,"lt_treat"]
    newarray_state["plt_treat", , , t ] <- pop_array[, , t]*oldPop[,"plt_treat"]
    
    newarray_state["a_treat_f", , , t] <- pop_array[, , t]*oldPop[,"a_treat_f"]
    newarray_state["f0_treat_f", , , t] <- pop_array[, , t]*oldPop[,"f0_treat_f"]
    newarray_state["f1_treat_f", , , t] <- pop_array[, , t]*oldPop[,"f1_treat_f"]
    newarray_state["f2_treat_f", , , t] <- pop_array[, , t]*oldPop[,"f2_treat_f"]
    newarray_state["f3_treat_f" , , , t] <- pop_array[, , t]*oldPop[,"f3_treat_f"]
    newarray_state["f4_treat_f", , , t] <- pop_array[, , t]*oldPop[,"f4_treat_f"]
    newarray_state["dc_treat_f", , , t] <- pop_array[, , t]*oldPop[,"dc_treat_f"]
    newarray_state["hcc_treat_f", , , t] <- pop_array[, , t]*oldPop[,"hcc_treat_f"]
    newarray_state["lt_treat_f", , , t] <- pop_array[, , t]*oldPop[,"lt_treat_f"]
    newarray_state["plt_treat_f", , , t] <- pop_array[, , t]*oldPop[,"plt_treat_f"]
    
    newarray_state["a_cured", , , t] <- pop_array[, , t]*oldPop[,"a_cured"]
    newarray_state["f0_cured", , , t] <- pop_array[, , t]*oldPop[,"f0_cured"]
    newarray_state["f1_cured", , , t] <- pop_array[, , t]*oldPop[,"f1_cured"]
    newarray_state["f2_cured", , , t] <- pop_array[, , t]*oldPop[,"f2_cured"]
    newarray_state["f3_cured", , , t] <- pop_array[, , t]*oldPop[,"f3_cured"]
    newarray_state["f4_cured", , , t] <- pop_array[, , t]*oldPop[,"f4_cured"]
    newarray_state["dc_cured", , , t] <- pop_array[, , t]*oldPop[,"dc_cured"]
    newarray_state["hcc_cured", , , t] <- pop_array[, , t]*oldPop[,"hcc_cured"]
    newarray_state["lt_cured", , , t] <- pop_array[, , t]*oldPop[,"lt_cured"]
    newarray_state["plt_cured", , , t] <- pop_array[, , t]*oldPop[,"plt_cured"]
    
    
    
    # pop-array 
    newarray[ , , t] <- 
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
  #### in and out flow bewtween subpops ####  
    
    inflow_hcv[, t] <- 
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
      colSums(pop_array[, , t]*oldPop[, "plt_treat_f"]) 
      
    
    #
    
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
    #
    
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
    
    outflow_hcv[, t]<- 
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
      rowSums(pop_array[, , t]*oldPop[, "plt_treat_f"]) 
    
    
    #### cost ####
    if (!is.null(cost) & proj != "POC_AU"){ 
      costPops[, , t] <- allPops[, , t]*cost$state[,]
      
      costTestingAb[, t] <- (costflow[,"ctau_ab", t]*newTestingAb[, t]) + 
        (costflow_Neg[,"ctau_ab", t]*((oldPop[, "s"] + oldPop[, "a_cured"])*tau_ab_dt[, "f0", t]))
    
    if(sum(tau_poct_dt[, "f0", t])==0){ # if not through poct pathway
      costTestingAg[, t] <- costflow[,"ctau_RNA", t]*(newTestingAg[, t] + 
         ((S[,] - oldPop[,"s"] - oldPop[, "a_cured"])*tau_RNA_dt[, "f0", t])) # all those cured from HCV (except for a_cured) had same probability to receive RNA testing 
    }else if(sum(tau_poct_dt[, "f0", t])!=0){ 
      costTestingAg[, t] <- costflow[,"ctau_RNA", t]*newTestingAg[, t]
      
      }
    
    costnewTestingPOCT[, t] <- costflow[,"ctau_poct",t ]*
      (newTestingPOCT[, t] + S[,]*tau_poct_dt[, "f0", t])
    
    
    
    costTreatment[, t] <- newTreatment[, t]*costflow[,"ceta", t]
    
    
    costCured[, t] <- newCured[, t]*costflow[,"ccured", t]
    
    
    costRetreat[, t] <- newRetreat[, t]*costflow[,"crho", t]
    
    
    QALYPops[, , t] <- allPops[, , t]*cost$QALY[,,t]
    
    }
    
    if (!is.null(cost)& proj == "POC_AU"){ 
      costPops[, , t] <- allPops[, , t]*cost$state[,]
      # costflow is a list and contain the unit cost for the standard of care and intervnetion
      # costflow[[1]]: unit cost for standard of care 
      # costflow[[2]]: unit cost for intervention 
      # same rules apply to costflow_Neg 
      
      costTestingAb[, t] <- (costflow[[1]][,"ctau_ab"]*(newTestingAb[, t])) +
        (costflow[[2]][,"ctau_ab"]*(newTestingAb_sc[, t])) + 
        (costflow_Neg[[1]][,"ctau_ab"]*(newTestingAb_neg[,t])) + 
        (costflow_Neg[[2]][,"ctau_ab"]*newTestingAb_sc_neg[,t])
      
      costTestingAg[, t] <- (costflow[[1]][,"ctau_ag"]*(newTestingAg[, t])) + 
        (costflow[[2]][,"ctau_ag"]*(newTestingAg_sc[, t])) + 
        (costflow_Neg[[1]][,"ctau_ag"]*(newTestingAg_neg[, t])) +  # all those cured from HCV (except for a_cured) had same probability to receive RNA testing 
        costflow_Neg[[2]][,"ctau_ag"]*(newTestingAg_sc_neg[, t])
      
      
      costTestingPOCT[, t] <- costflow[[1]][,"ctau_poct"]*(newTestingPOCT[, t]) + 
        costflow[[2]][,"ctau_poct"]*(newTestingPOCT_sc[, t]) + 
        (costflow_Neg[[1]][,"ctau_poct"]*(newTestingPOCT_neg[, t])) + 
        (costflow_Neg[[2]][,"ctau_poct"]*(newTestingPOCT_sc_neg[, t]))
      
      costTreatment[, t] <- newTreatment[, t]*costflow[[1]][,"ceta"] + 
        newTreatment_sc[, t]*costflow[[1]][,"ceta"]
      
      costCured[, t] <- newCured[, t]*costflow[[1]][,"ccured"]
      
      
      costRetreat[, t] <- newRetreat[, t]*costflow[[1]][,"crho"]
      
      
      QALYPops[, , t] <- allPops[, , t]*cost$QALY[,]
      
      
      costTestingAb_sc[, t] <- 
        (costflow[[2]][,"ctau_ab"]*(newTestingAb_sc[, t])) + 
        (costflow_Neg[[2]][,"ctau_ab"]*newTestingAb_sc_neg[,t])
      
      costTestingAg_sc[, t] <- 
        (costflow[[2]][,"ctau_ag"]*(newTestingAg_sc[, t])) + 
        costflow_Neg[[2]][,"ctau_ag"]*(newTestingAg_sc_neg[, t])
      
      
      costTestingPOCT_sc[, t] <- 
        costflow[[2]][,"ctau_poct"]*(newTestingPOCT_sc[, t]) + 
        (costflow_Neg[[2]][,"ctau_poct"]*(newTestingPOCT_sc_neg[, t]))
      
      costTreatment_sc[, t] <- newTreatment_sc[, t]*costflow[[1]][,"ceta"]
      
      
      
      
    }
    
  }
  
  #------------------------- results output ---------------------------------#
  
  if (!is.null(cost) & proj == "POC_AU"){
    results <- list(allPops = allPops,
                    newS = newS,
                    newEntry = newEntry,
                    newDeath = newDeath,
                    newLeave = newLeave, 
                    newInfections = newInfections, 
                    newHCVdeaths = newHCVdeaths, 
                    newTreatment = newTreatment, 
                    newRetreat = newRetreat, 
                    newTestingAb = newTestingAb, 
                    newTestingAg = newTestingAg,
                    newTestingPOCT = newTestingPOCT,
                    newTestingAb_sc = newTestingAb_sc,
                    newTestingAg_sc = newTestingAg_sc,
                    newTestingPOCT_sc = newTestingPOCT_sc,
                    newTreatment_sc = newTreatment_sc, 
                    newTestingAb_neg = newTestingAb_neg, 
                    newTestingAg_neg = newTestingAg_neg,
                    newTestingPOCT_neg = newTestingPOCT_neg,
                    newTestingAb_sc_neg = newTestingAb_sc_neg,
                    newTestingAg_sc_neg = newTestingAg_sc_neg,
                    newTestingPOCT_sc_neg = newTestingPOCT_sc_neg,
                    newCured = newCured,
                    newtreatfailed = newtreatfailed,
                    newreinfection = newreinfection,
                    newreinfection_chronic = newreinfection_chronic, 
                    newpop_tran = newarray,
                    newpop_tranState = newarray_state,
                    inflow = inflow,
                    inflow_hcv = inflow_hcv,
                    outflow = outflow,
                    outflow_hcv = outflow_hcv, 
                    death_hcv = death_hcv,
                    HCVdeathState = HCVdeathState,
                    newDeathState = newDeathState,
                    costPops = costPops,
                    QALYPops = QALYPops,
                    costTestingAb = costTestingAb,
                    costTestingAg = costTestingAg,
                    costTestingPOCT = costTestingPOCT,
                    costTreatment = costTreatment,
                    costCured = costCured,
                    costRetreat = costRetreat,
                    costTestingAb_sc = costTestingAb_sc,
                    costTestingAg_sc = costTestingAg_sc,
                    costTestingPOCT_sc = costTestingPOCT_sc,
                    costTreatment_sc = costTreatment_sc)
  }else if (!is.null(cost) & proj != "POC_AU"){
    results <- list(allPops = allPops,
                    newS = newS,
                    newEntry = newEntry,
                    newDeath = newDeath,
                    newLeave = newLeave, 
                    newInfections = newInfections, 
                    newHCVdeaths = newHCVdeaths, 
                    newTreatment = newTreatment, 
                    newRetreat = newRetreat, 
                    newTestingAb = newTestingAb, 
                    newTestingAg = newTestingAg,
                    newTestingPOCT = newTestingPOCT,
                    newTestingAb_sc = newTestingAb_sc,
                    newTestingAg_sc = newTestingAg_sc,
                    newTestingPOCT_sc = newTestingPOCT_sc,
                    newTreatment_sc = newTreatment_sc, 
                    newTestingAb_neg = newTestingAb_neg, 
                    newTestingAg_neg = newTestingAg_neg,
                    newTestingPOCT_neg = newTestingPOCT_neg,
                    newTestingAb_sc_neg = newTestingAb_sc_neg,
                    newTestingAg_sc_neg = newTestingAg_sc_neg,
                    newTestingPOCT_sc_neg = newTestingPOCT_sc_neg,
                    newCured = newCured,
                    newreinfection = newreinfection,
                    newreinfection_chronic = newreinfection_chronic,
                    newpop_tran = newarray,
                    newpop_tranState = newarray_state,
                    inflow = inflow,
                    inflow_hcv = inflow_hcv,
                    outflow = outflow,
                    outflow_hcv = outflow_hcv, 
                    death_hcv = death_hcv,
                    HCVdeathState = HCVdeathState,
                    newDeathState = newDeathState,
                    costPops = costPops,
                    QALYPops = QALYPops,
                    costTestingAb = costTestingAb,
                    costTestingAg = costTestingAg,
                    costTestingPOCT = costTestingPOCT,
                    costTreatment = costTreatment,
                    costCured = costCured,
                    costRetreat = costRetreat)
    
    
  }else if (is.null(cost)){ 
    results <- list(allPops = allPops,
                    newS = newS,
                    newEntry = newEntry,
                    newDeath = newDeath,
                    newLeave = newLeave,
                    newInfections = newInfections, 
                    newHCVdeaths = newHCVdeaths, 
                    newTreatment = newTreatment, 
                    newRetreat = newRetreat, 
                    newTestingAb = newTestingAb, 
                    newTestingAg = newTestingAg,
                    newTestingPOCT = newTestingPOCT,
                    newTestingAb_sc = newTestingAb_sc,
                    newTestingAg_sc = newTestingAg_sc,
                    newTestingPOCT_sc = newTestingPOCT_sc, 
                    newTreatment_sc = newTreatment_sc,
                    newTestingAb_neg = newTestingAb_neg, 
                    newTestingAg_neg = newTestingAg_neg,
                    newTestingPOCT_neg = newTestingPOCT_neg,
                    newTestingAb_sc_neg = newTestingAb_sc_neg,
                    newTestingAg_sc_neg = newTestingAg_sc_neg,
                    newTestingPOCT_sc_neg = newTestingPOCT_sc_neg,
                    newCured = newCured,
                    newtreatfailed = newtreatfailed,
                    newreinfection = newreinfection,
                    newreinfection_chronic = newreinfection_chronic,
                    newpop_tran = newarray,
                    newpop_tranState = newarray_state,
                    inflow = inflow,
                    inflow_hcv = inflow_hcv,
                    outflow = outflow,
                    outflow_hcv = outflow_hcv, 
                    death_hcv = death_hcv,
                    HCVdeathState = HCVdeathState,
                    newDeathState = newDeathState)
    }
    
 
  
  
  
}
