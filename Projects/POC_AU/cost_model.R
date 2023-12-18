# this script is a function for cost calculation in the model 
cost_model <- function(pj, costdfList, coverage, dfList, dfList_NP, 
                       Sce_npscale,
                       S_Yint, S_Yfir = NULL, S_Ymid = NULL, S_Yend,
                       censor_Y = NULL){
  # pj: project 
  # costdfList: cost for each state and flows, including QALY! 
  # coverage: the coverage of an intervention
  # dfList: the coverage of each step of cascade in status quo 
  # dfList_NP : the coverage of each step of cascade in a scenario 
  # Sce_npscale: the model outcome of a scenario 
  # S_Yint: the beinging year of a scenario  
  # S_Yfir: the final year of the first stage intervention 
  # S_Ymid: the finished year of a scenario
  # S_Yend: the final year of the model
  
  SYpoint_int <- (S_Yint - pj$cabY)/pj$timestep + 1
  SYpoint_fir <- (S_Yfir - pj$cabY)/pj$timestep + 1
  SYpoint_mid <- (S_Ymid - pj$cabY)/pj$timestep + 1
  SYpoint_end <- (S_Yend - pj$cabY)/pj$timestep + 1
  
  SY_leng <- SYpoint_end - SYpoint_int
  SY_lengmid <- SYpoint_mid - SYpoint_int
  SY_lengfirtomid <- SYpoint_mid - SYpoint_fir
  # applying cost to each timestep 
  
  cost <- list()
  
  cost[["state"]] <- costdfList$costPops*pj$timestep
  
  cost[["QALY"]] <- costdfList$QALYPops*pj$timestep
  
  cost[["flow"]] <- costdfList$costFlow*pj$timestep
  cost[["flow_NEG"]] <- costdfList$costFlow_NEG*pj$timestep
  
  cost[["np"]] <- costdfList$costFlow_POCRNA*pj$timestep
  
  cost[["np_NEG"]] <- costdfList$`costFlow_POCRNA _NEG`*pj$timestep
  
  # set up the format for storing cost outcomes
  costPops <- array(0, c(pj$npops, pj$ncomponent, pj$npts), 
                    dimnames = list(pj$popNames, pj$component_name))
  QALYPops <- array(0, c(pj$npops, pj$ncomponent, pj$npts), 
                    dimnames = list(pj$popNames, pj$component_name))
  
  ResultMatrix <- matrix(0, ncol = pj$npts, nrow = pj$npops)
  
  costTestingAb <- ResultMatrix 
  
  costTestingRNA <- ResultMatrix
  
  costTestingPOCT <- ResultMatrix
  
  costTreatment <- ResultMatrix
  
  costRetreat <- ResultMatrix
  
  costCured <- ResultMatrix 
  
  # soc: status quo 
  # np: scenario
  
  newTestingAb_soc <- ResultMatrix
  
  newTestingAb_np <- ResultMatrix
  
  newTestingAbNEG_soc <- ResultMatrix
  
  newTestingAbNEG_np <- ResultMatrix
  
  newTestingRNA_soc <- ResultMatrix
  
  newTestingRNA_np <- ResultMatrix
  
  newTestingRNANEG_soc <- ResultMatrix
  
  newTestingRNANEG_np <- ResultMatrix
  
  newTestingPOCT_soc <- ResultMatrix
  
  newTestingPOCT_np <- ResultMatrix
  
  newTestingPOCTNEG_soc <- ResultMatrix
  
  newTestingPOCTNEG_np <- ResultMatrix
  
  newTreatment_soc <- ResultMatrix
  
  newTreatment_np <- ResultMatrix
  
  newRetreat_soc <- ResultMatrix
  
  newRetreat_np <- ResultMatrix
  
  newCured_soc <- ResultMatrix
  
  newCured_np <- ResultMatrix
  
  spc1 <-  matrix(0, ncol = pj$npts + 1, nrow = pj$npops)
 
  for(i in 1:pj$npops){
    spc1[i, ] <- best_estimates[1 , paste0("spc", i)]
  }
  
  spc1_dt <- 1-(1-spc1)^pj$timestep
  
  
  soc_tau_ab <- 1-(1-dfList$tau_ab)^pj$timestep
  np_tau_ab <- 1-(1-dfList_NP$tau_ab)^pj$timestep
  
  soc_tau_rna <- 1-(1-dfList$tau_RNA)^pj$timestep
  np_tau_rna <- 1-(1-dfList_NP$tau_RNA)^pj$timestep
  
  soc_tau_poct <- 1-(1-dfList$tau_poct)^pj$timestep
  np_tau_poct <- 1-(1-dfList_NP$tau_poct)^pj$timestep
  
  soc_eta <- 1-(1-dfList$eta)^pj$timestep
  np_eta <- 1-(1-dfList_NP$eta)^pj$timestep
  
  soc_eta <- 1-(1-dfList$eta)^pj$timestep
  np_eta <- 1-(1-dfList_NP$eta)^pj$timestep
  
  soc_rho <- 1-(1-dfList$rho)^pj$timestep
  np_rho <- 1-(1-dfList_NP$rho)^pj$timestep
  
  soc_cured <- 1-(1-dfList$cured)^pj$timestep
  np_cured <- 1-(1-dfList_NP$cured)^pj$timestep
  
  
  
  # set up the additional coverage due to a scneario 
  for(i in 1: pj$npops){ 
    # soc_tau_ab == np_tau_ab
    soc_tau_ab[i,,] <- soc_tau_ab[i,,]*(1-coverage[i,])
    
    np_tau_ab[i,,] <- np_tau_ab[i,,]*coverage[i,]
    
    
    soc_rho[i,,] <- soc_rho[i,,]*(1-coverage[i,])
    
    np_rho[i,,] <- np_rho[i,,]*coverage[i,]
    
    soc_cured[i,,] <- soc_cured[i,,]*(1-coverage[i,])
    
    np_cured[i,,] <- np_cured[i,,]*coverage[i,]
  }
  
  if(!is.null(censor_Y)){ 
    censor_Y_point <- (censor_Y + 1 - pj$cabY)/pj$timestep + 1 
    
    for(i in 1: pj$npops){ 
      for(t in censor_Y_point: dim(soc_tau_ab)[3]){
        # soc_tau_ab == np_tau_ab
        soc_tau_ab[i,,t] <- soc_tau_ab[i,,t]
        
        np_tau_ab[i,,t] <- 0
        
        soc_rho[i,,t] <- soc_rho[i,,t]
        
        np_rho[i,,t] <- 0
        
        soc_cured[i,,t] <- soc_cured[i,,t]
        
        np_cured[i,,t] <- 0
        
        np_tau_rna[i,,t] <- soc_tau_rna[i,,t]
        
        np_tau_poct[i,,t] <- soc_tau_poct[i,,t]
        
        np_eta[i,,t] <- soc_eta[i,,t]
        
      }
      
    }
    
      }
    
  
  
  
  
  
  for(i in 1: pj$npops){
    for(t in SYpoint_int: SYpoint_end){ 
      newTestingAb_soc[i, t] <- 
        sum(soc_tau_ab[i,"a",t]*(1-spc1_dt[i, t])*Sce_npscale$allPops[i,"a_undiag", (t-1)],
            soc_tau_ab[i,"f0",t]*Sce_npscale$allPops[i,"f0_undiag", (t-1)],
            soc_tau_ab[i,"f1",t]*Sce_npscale$allPops[i,"f1_undiag", (t-1)],
            soc_tau_ab[i,"f2",t]*Sce_npscale$allPops[i,"f2_undiag", (t-1)],
            soc_tau_ab[i,"f3",t]*Sce_npscale$allPops[i,"f3_undiag", (t-1)],
            soc_tau_ab[i,"f4",t]*Sce_npscale$allPops[i,"f4_undiag", (t-1)],
            soc_tau_ab[i,"dc",t]*Sce_npscale$allPops[i,"dc_undiag", (t-1)],
            soc_tau_ab[i,"hcc",t]*Sce_npscale$allPops[i,"hcc_undiag", (t-1)],
            soc_tau_ab[i,"lt",t]*Sce_npscale$allPops[i,"lt_undiag", (t-1)],
            soc_tau_ab[i,"plt",t]*Sce_npscale$allPops[i,"plt_undiag", (t-1)])
      
      newTestingAbNEG_soc[i, t] <- 
        (Sce_npscale$allPops[i, "s", (t-1)] + Sce_npscale$allPops[i, "a_cured", (t-1)])*soc_tau_ab[i,"f0",t]
      
      newTestingAb_np[i, t] <- 
        sum(np_tau_ab[i,"a",t]*(1-spc1_dt[i, t])*Sce_npscale$allPops[i,"a_undiag", (t-1)],
            np_tau_ab[i,"f0",t]*Sce_npscale$allPops[i,"f0_undiag", (t-1)],
            np_tau_ab[i,"f1",t]*Sce_npscale$allPops[i,"f1_undiag", (t-1)],
            np_tau_ab[i,"f2",t]*Sce_npscale$allPops[i,"f2_undiag", (t-1)],
            np_tau_ab[i,"f3",t]*Sce_npscale$allPops[i,"f3_undiag", (t-1)],
            np_tau_ab[i,"f4",t]*Sce_npscale$allPops[i,"f4_undiag", (t-1)],
            np_tau_ab[i,"dc",t]*Sce_npscale$allPops[i,"dc_undiag", (t-1)],
            np_tau_ab[i,"hcc",t]*Sce_npscale$allPops[i,"hcc_undiag", (t-1)],
            np_tau_ab[i,"lt",t]*Sce_npscale$allPops[i,"lt_undiag", (t-1)],
            np_tau_ab[i,"plt",t]*Sce_npscale$allPops[i,"plt_undiag", (t-1)])
      
      newTestingAbNEG_np[i, t] <- 
        (Sce_npscale$allPops[i, "s", (t-1)] + Sce_npscale$allPops[i, "a_cured", (t-1)])*np_tau_ab[i,"f0",t]
      
      # RNA testing 
      newTestingRNA_soc[i, t] <- 
        sum(soc_tau_rna[i,"a",t]*(1-spc1_dt[i, t])*Sce_npscale$allPops[i,"a_diag_ab", (t-1)],
            soc_tau_rna[i,"f0",t]*Sce_npscale$allPops[i,"f0_diag_ab", (t-1)],
            soc_tau_rna[i,"f1",t]*Sce_npscale$allPops[i,"f1_diag_ab", (t-1)],
            soc_tau_rna[i,"f2",t]*Sce_npscale$allPops[i,"f2_diag_ab", (t-1)],
            soc_tau_rna[i,"f3",t]*Sce_npscale$allPops[i,"f3_diag_ab", (t-1)],
            soc_tau_rna[i,"f4",t]*Sce_npscale$allPops[i,"f4_diag_ab", (t-1)],
            soc_tau_rna[i,"dc",t]*Sce_npscale$allPops[i,"dc_diag_ab", (t-1)],
            soc_tau_rna[i,"hcc",t]*Sce_npscale$allPops[i,"hcc_diag_ab", (t-1)],
            soc_tau_rna[i,"lt",t]*Sce_npscale$allPops[i,"lt_diag_ab", (t-1)],
            soc_tau_ab[i,"plt",t]*Sce_npscale$allPops[i,"plt_diag_ab", (t-1)])
      
      newTestingRNANEG_soc[i, t] <- 
        (Sce_npscale$allPops[i, "f0_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f1_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f2_cured", (t-1)] +
           Sce_npscale$allPops[i, "f3_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f4_cured", (t-1)] +
           Sce_npscale$allPops[i, "dc_cured", (t-1)] +
           Sce_npscale$allPops[i, "hcc_cured", (t-1)] +
           Sce_npscale$allPops[i, "lt_cured", (t-1)] +
           Sce_npscale$allPops[i, "plt_cured", (t-1)])*soc_tau_rna[i,"f0",t]
      
      newTestingRNA_np[i, t] <- 
        sum((np_tau_rna[i,"a",t] - soc_tau_rna[i,"a",t])*(1-spc1_dt[i, t])*Sce_npscale$allPops[i,"a_diag_ab", (t-1)],
            (np_tau_rna[i,"f0",t] - soc_tau_rna[i,"f0",t])*Sce_npscale$allPops[i,"f0_diag_ab", (t-1)],
            (np_tau_rna[i,"f1",t] - soc_tau_rna[i,"f1",t])*Sce_npscale$allPops[i,"f1_diag_ab", (t-1)],
            (np_tau_rna[i,"f2",t] - soc_tau_rna[i,"f2",t])*Sce_npscale$allPops[i,"f2_diag_ab", (t-1)],
            (np_tau_rna[i,"f3",t] - soc_tau_rna[i,"f3",t])*Sce_npscale$allPops[i,"f3_diag_ab", (t-1)],
            (np_tau_rna[i,"f4",t] - soc_tau_rna[i,"f4",t])*Sce_npscale$allPops[i,"f4_diag_ab", (t-1)],
            (np_tau_rna[i,"dc",t] - soc_tau_rna[i,"dc",t])*Sce_npscale$allPops[i,"dc_diag_ab", (t-1)],
            (np_tau_rna[i,"hcc",t] - soc_tau_rna[i,"hcc",t])*Sce_npscale$allPops[i,"hcc_diag_ab", (t-1)],
            (np_tau_rna[i,"lt",t] - soc_tau_rna[i,"lt",t])*Sce_npscale$allPops[i,"lt_diag_ab", (t-1)],
            (np_tau_rna[i,"plt",t] - soc_tau_rna[i,"plt",t])*Sce_npscale$allPops[i,"plt_diag_ab", (t-1)])
      
      newTestingRNANEG_np[i, t] <- 
        (Sce_npscale$allPops[i, "f0_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f1_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f2_cured", (t-1)] +
           Sce_npscale$allPops[i, "f3_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f4_cured", (t-1)] +
           Sce_npscale$allPops[i, "dc_cured", (t-1)] +
           Sce_npscale$allPops[i, "hcc_cured", (t-1)] +
           Sce_npscale$allPops[i, "lt_cured", (t-1)] +
           Sce_npscale$allPops[i, "plt_cured", (t-1)])*
        (np_tau_rna[i,"f0",t] - soc_tau_rna[i,"f0",t])
      
      
      # POCT 
      newTestingPOCT_soc[i, t] <- 
        sum(soc_tau_poct[i,"a",t]*(1-spc1_dt[i, t])*Sce_npscale$allPops[i,"a_undiag", (t-1)],
            soc_tau_poct[i,"f0",t]*Sce_npscale$allPops[i,"f0_undiag", (t-1)],
            soc_tau_poct[i,"f1",t]*Sce_npscale$allPops[i,"f1_undiag", (t-1)],
            soc_tau_poct[i,"f2",t]*Sce_npscale$allPops[i,"f2_undiag", (t-1)],
            soc_tau_poct[i,"f3",t]*Sce_npscale$allPops[i,"f3_undiag", (t-1)],
            soc_tau_poct[i,"f4",t]*Sce_npscale$allPops[i,"f4_undiag", (t-1)],
            soc_tau_poct[i,"dc",t]*Sce_npscale$allPops[i,"dc_undiag", (t-1)],
            soc_tau_poct[i,"hcc",t]*Sce_npscale$allPops[i,"hcc_undiag", (t-1)],
            soc_tau_poct[i,"lt",t]*Sce_npscale$allPops[i,"lt_undiag", (t-1)],
            soc_tau_poct[i,"plt",t]*Sce_npscale$allPops[i,"plt_undiag", (t-1)])
      
      newTestingPOCTNEG_soc[i, t] <- 
        (Sce_npscale$allPops[i, "s", (t-1)] + 
           Sce_npscale$allPops[i, "a_cured", (t-1)] +
           Sce_npscale$allPops[i, "f0_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f1_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f2_cured", (t-1)] +
           Sce_npscale$allPops[i, "f3_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f4_cured", (t-1)] +
           Sce_npscale$allPops[i, "dc_cured", (t-1)] +
           Sce_npscale$allPops[i, "hcc_cured", (t-1)] +
           Sce_npscale$allPops[i, "lt_cured", (t-1)] +
           Sce_npscale$allPops[i, "plt_cured", (t-1)])*soc_tau_poct[i,"f0",t]
      
      newTestingPOCT_np[i, t] <- 
        sum((np_tau_poct[i,"a",t] - soc_tau_poct[i,"a",t])*(1-spc1_dt[i, t])*Sce_npscale$allPops[i,"a_undiag", (t-1)],
            (np_tau_poct[i,"f0",t] - soc_tau_poct[i,"f0",t])*Sce_npscale$allPops[i,"f0_undiag", (t-1)],
            (np_tau_poct[i,"f1",t] - soc_tau_poct[i,"f1",t])*Sce_npscale$allPops[i,"f1_undiag", (t-1)],
            (np_tau_poct[i,"f2",t] - soc_tau_poct[i,"f2",t])*Sce_npscale$allPops[i,"f2_undiag", (t-1)],
            (np_tau_poct[i,"f3",t] - soc_tau_poct[i,"f3",t])*Sce_npscale$allPops[i,"f3_undiag", (t-1)],
            (np_tau_poct[i,"f4",t] - soc_tau_poct[i,"f4",t])*Sce_npscale$allPops[i,"f4_undiag", (t-1)],
            (np_tau_poct[i,"dc",t] - soc_tau_poct[i,"dc",t])*Sce_npscale$allPops[i,"dc_undiag", (t-1)],
            (np_tau_poct[i,"hcc",t] - soc_tau_poct[i,"hcc",t])*Sce_npscale$allPops[i,"hcc_undiag", (t-1)],
            (np_tau_poct[i,"lt",t] - soc_tau_poct[i,"lt",t])*Sce_npscale$allPops[i,"lt_undiag", (t-1)],
            (np_tau_poct[i,"plt",t] - soc_tau_poct[i,"plt",t])*Sce_npscale$allPops[i,"plt_undiag", (t-1)])
      
      newTestingPOCTNEG_np[i, t] <- 
        (Sce_npscale$allPops[i, "s", (t-1)] + 
           Sce_npscale$allPops[i, "a_cured", (t-1)] +
           Sce_npscale$allPops[i, "f0_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f1_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f2_cured", (t-1)] +
           Sce_npscale$allPops[i, "f3_cured", (t-1)] + 
           Sce_npscale$allPops[i, "f4_cured", (t-1)] +
           Sce_npscale$allPops[i, "dc_cured", (t-1)] +
           Sce_npscale$allPops[i, "hcc_cured", (t-1)] +
           Sce_npscale$allPops[i, "lt_cured", (t-1)] +
           Sce_npscale$allPops[i, "plt_cured", (t-1)])*
        (np_tau_poct[i,"f0",t] - soc_tau_poct[i,"f0",t])
      
      
      
      # treatment 
      newTreatment_soc[i, t] <- 
        sum(soc_eta[i, "a", t]*Sce_npscale$allPops[i,"a_diag_RNA", (t-1)],
            soc_eta[i, "f0", t]*Sce_npscale$allPops[i,"f0_diag_RNA", (t-1)],
            soc_eta[i, "f1", t]*Sce_npscale$allPops[i,"f1_diag_RNA", (t-1)],
            soc_eta[i, "f2", t]*Sce_npscale$allPops[i,"f2_diag_RNA", (t-1)],
            soc_eta[i, "f3", t]*Sce_npscale$allPops[i,"f3_diag_RNA", (t-1)],
            soc_eta[i, "f4", t]*Sce_npscale$allPops[i,"f4_diag_RNA", (t-1)],
            soc_eta[i, "dc", t]*Sce_npscale$allPops[i,"dc_diag_RNA", (t-1)],
            soc_eta[i, "hcc", t]*Sce_npscale$allPops[i,"hcc_diag_RNA", (t-1)],
            soc_eta[i, "lt", t]*Sce_npscale$allPops[i,"lt_diag_RNA", (t-1)],
            soc_eta[i, "plt", t]*Sce_npscale$allPops[i,"plt_diag_RNA", (t-1)])
      
      newTreatment_np[i,t] <- 
        sum((np_eta[i, "a", t] - soc_eta[i, "a", t])*Sce_npscale$allPops[i,"a_diag_RNA", (t-1)],
            (np_eta[i, "f0", t] - soc_eta[i, "f0", t])*Sce_npscale$allPops[i,"f0_diag_RNA", (t-1)],
            (np_eta[i, "f1", t] - soc_eta[i, "f1", t])*Sce_npscale$allPops[i,"f1_diag_RNA", (t-1)],
            (np_eta[i, "f2", t] - soc_eta[i, "f2", t])*Sce_npscale$allPops[i,"f2_diag_RNA", (t-1)],
            (np_eta[i, "f3", t] - soc_eta[i, "f3", t])*Sce_npscale$allPops[i,"f3_diag_RNA", (t-1)],
            (np_eta[i, "f4", t] - soc_eta[i, "f4", t])*Sce_npscale$allPops[i,"f4_diag_RNA", (t-1)],
            (np_eta[i, "dc", t] - soc_eta[i, "dc", t])*Sce_npscale$allPops[i,"dc_diag_RNA", (t-1)],
            (np_eta[i, "hcc", t] - soc_eta[i, "hcc", t])*Sce_npscale$allPops[i,"hcc_diag_RNA", (t-1)],
            (np_eta[i, "lt", t] - soc_eta[i, "lt", t])*Sce_npscale$allPops[i,"lt_diag_RNA", (t-1)],
            (np_eta[i, "plt", t] - soc_eta[i, "plt", t])*Sce_npscale$allPops[i,"plt_diag_RNA", (t-1)])
      
      
      # retreat
      newRetreat_soc[i,t] <- 
        sum(soc_rho[i, "a", t]*Sce_npscale$allPops[i,"a_treat_f",(t-1) ],
            soc_rho[i, "f0", t]*Sce_npscale$allPops[i,"f0_treat_f",(t-1)],
            soc_rho[i, "f1", t]*Sce_npscale$allPops[i,"f1_treat_f",(t-1)],
            soc_rho[i, "f2", t]*Sce_npscale$allPops[i,"f2_treat_f",(t-1)],
            soc_rho[i, "f3", t]*Sce_npscale$allPops[i,"f3_treat_f",(t-1)],
            soc_rho[i, "f4", t]*Sce_npscale$allPops[i,"f4_treat_f",(t-1)],
            soc_rho[i, "dc", t]*Sce_npscale$allPops[i,"dc_treat_f",(t-1)],
            soc_rho[i, "hcc", t]*Sce_npscale$allPops[i,"hcc_treat_f",(t-1)],
            soc_rho[i, "lt", t]*Sce_npscale$allPops[i,"lt_treat_f",(t-1)],
            soc_rho[i, "plt", t]*Sce_npscale$allPops[i,"plt_treat_f",(t-1)])
      
      newRetreat_np[i,t] <- 
        sum(np_rho[i, "a", t]*Sce_npscale$allPops[i,"a_treat_f",(t-1) ],
            np_rho[i, "f0", t]*Sce_npscale$allPops[i,"f0_treat_f",(t-1)],
            np_rho[i, "f1", t]*Sce_npscale$allPops[i,"f1_treat_f",(t-1)],
            np_rho[i, "f2", t]*Sce_npscale$allPops[i,"f2_treat_f",(t-1)],
            np_rho[i, "f3", t]*Sce_npscale$allPops[i,"f3_treat_f",(t-1)],
            np_rho[i, "f4", t]*Sce_npscale$allPops[i,"f4_treat_f",(t-1)],
            np_rho[i, "dc", t]*Sce_npscale$allPops[i,"dc_treat_f",(t-1)],
            np_rho[i, "hcc", t]*Sce_npscale$allPops[i,"hcc_treat_f",(t-1)],
            np_rho[i, "lt", t]*Sce_npscale$allPops[i,"lt_treat_f",(t-1)],
            np_rho[i, "plt", t]*Sce_npscale$allPops[i,"plt_treat_f",(t-1)])
      
      # cured 
      newCured_soc[i, t] <- sum(soc_cured[i, "a", t]*Sce_npscale$allPops[i,"a_treat",(t-1)],
                                soc_cured[i, "f0", t]*Sce_npscale$allPops[i,"f0_treat",(t-1)],
                                soc_cured[i, "f1", t]*Sce_npscale$allPops[i,"f1_treat",(t-1)],
                                soc_cured[i, "f2", t]*Sce_npscale$allPops[i,"f2_treat",(t-1)],
                                soc_cured[i, "f3", t]*Sce_npscale$allPops[i,"f3_treat",(t-1)],
                                soc_cured[i, "f4", t]*Sce_npscale$allPops[i,"f4_treat",(t-1)],
                                soc_cured[i, "dc", t]*Sce_npscale$allPops[i,"dc_treat",(t-1)],
                                soc_cured[i, "hcc", t]*Sce_npscale$allPops[i,"hcc_treat",(t-1)],
                                soc_cured[i, "lt", t]*Sce_npscale$allPops[i,"lt_treat",(t-1)],
                                soc_cured[i, "plt", t]*Sce_npscale$allPops[i,"plt_treat",(t-1)])
      
      newCured_np[i, t] <- sum(np_cured[i, "a", t]*Sce_npscale$allPops[i,"a_treat",(t-1)],
                               np_cured[i, "f0", t]*Sce_npscale$allPops[i,"f0_treat",(t-1)],
                               np_cured[i, "f1", t]*Sce_npscale$allPops[i,"f1_treat",(t-1)],
                               np_cured[i, "f2", t]*Sce_npscale$allPops[i,"f2_treat",(t-1)],
                               np_cured[i, "f3", t]*Sce_npscale$allPops[i,"f3_treat",(t-1)],
                               np_cured[i, "f4", t]*Sce_npscale$allPops[i,"f4_treat",(t-1)],
                               np_cured[i, "dc", t]*Sce_npscale$allPops[i,"dc_treat",(t-1)],
                               np_cured[i, "hcc", t]*Sce_npscale$allPops[i,"hcc_treat",(t-1)],
                               np_cured[i, "lt", t]*Sce_npscale$allPops[i,"lt_treat",(t-1)],
                               np_cured[i, "plt", t]*Sce_npscale$allPops[i,"plt_treat",(t-1)])
      
      
      
    }
    
  }
  
  
  for(t in SYpoint_int: SYpoint_end){ 
    costPops[, ,t] <- Sce_npscale$allPops[, , t]*cost[["state"]]
    QALYPops[, ,t] <- Sce_npscale$allPops[, , t]*0
    costTestingAb[, t] <- 
      ((cost[["flow"]][,"ctau_ab"]*newTestingAb_soc[, t]) +
         (cost[["np"]][,"ctau_ab"]*newTestingAb_np[, t])) +
      (cost[["flow_NEG"]][,"ctau_ab"]*newTestingAbNEG_soc[, t]) +
      (cost[["np_NEG"]][,"ctau_ab"]*newTestingAbNEG_np[, t])
    
    costTestingRNA[, t] <- 
      ((cost[["flow"]][,"ctau_ag"]*newTestingRNA_soc[, t]) +
         (cost[["np"]][,"ctau_ag"]*newTestingRNA_np[, t])) +
      (cost[["flow_NEG"]][,"ctau_ag"]*newTestingRNANEG_soc[, t]) +
      (cost[["np_NEG"]][,"ctau_ag"]*newTestingRNANEG_np[, t])
    
    costTestingPOCT[, t] <- 
      ((cost[["flow"]][,"ctau_poct"]*newTestingPOCT_soc[, t]) +
         (cost[["np"]][,"ctau_poct"]*newTestingPOCT_np[, t])) +
      (cost[["flow_NEG"]][,"ctau_poct"]*newTestingPOCTNEG_soc[, t]) +
      (cost[["np_NEG"]][,"ctau_poct"]*newTestingPOCTNEG_np[, t])
    
    costTreatment[, t] <- 
      (cost[["flow"]][,"ceta"]*newTreatment_soc[, t]) +
      (cost[["np"]][,"ceta"]*newTreatment_np[, t])
    
    costRetreat[, t] <- 
      (cost[["flow"]][,"ceta"]*newRetreat_soc[, t]) +
      (cost[["np"]][,"ceta"]*newRetreat_np[, t])
    
    costCured[, t] <- 
      (cost[["flow"]][,"ccured"]*newCured_soc[, t]) +
      (cost[["np"]][,"ccured"]*newCured_np[, t])
    
  }
  
  Sce_npscale <- append(Sce_npscale,list("costTestingAb" = costTestingAb,
                                          "costTestingRNA" = costTestingRNA,
                                          "costTestingPOCT" = costTestingPOCT,
                                          "costTreatment" = costTreatment,
                                          "costRetreat" = costRetreat,
                                          "costCured" = costCured,
                                          "costPops" = costPops,
                                          "QALYPops" = costPops))
  
  return(Sce_npscale)
}

