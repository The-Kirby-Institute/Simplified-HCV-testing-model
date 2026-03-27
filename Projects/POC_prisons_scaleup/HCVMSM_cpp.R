# ══════════════════════════════════════════════════════════════════════════════
# Compile the Rcpp code
# ══════════════════════════════════════════════════════════════════════════════
library(Rcpp)
library(RcppArmadillo)
sourceCpp("hcvmsm_core.cpp")

# ══════════════════════════════════════════════════════════════════════════════
# R wrapper — same interface as HCVMSM
# ══════════════════════════════════════════════════════════════════════════════
HCVMSM_cpp <- function(HCV, parama, initialPop, disease_progress,
                       pop_array, param_cascade, param_cascade_sc, fib,
                       end_Y = NULL, modelrun = NULL, proj = NULL,
                       fc_sc = NULL, fp = NULL) {
  
  dt <- HCV$timestep
  
  # ── Same preprocessing as original HCVMSM ────────────────────────────────
  if (is.null(end_Y)) {
    npts <- HCV$npts - 1
  } else if (!is.data.frame(parama)) {
    pts   <- seq(HCV$startYear, end_Y, by = dt)
    pts   <- pts[-1]
    npts  <- length(pts)
    parama    <- lapply(parama, function(x) x[1:(npts + 1)])
    pop_array <- pop_array[, , 1:npts]
  } else {
    pts   <- seq(HCV$startYear, end_Y, by = dt)
    pts   <- pts[-1]
    npts  <- length(pts)
    parama    <- parama[1:(npts + 1), ]
    pop_array <- pop_array[, , 1:npts]
  }
  
  npops     <- HCV$npops
  ncomponent <- HCV$ncomponent
  
  # ── fc ───────────────────────────────────────────────────────────────────
  fc <- matrix(0, ncol = npts + 1, nrow = npops)
  if (!is.null(fc_sc)) fc <- fc_sc else fc[, ] <- 1
  if (is.null(fp)) fp <- 0
  
  # ── reinfP ───────────────────────────────────────────────────────────────
  reinfP <- matrix(1, ncol = npts + 1, nrow = npops)
  
  # ── Mortality dt conversions ──────────────────────────────────────────────
  morb   <- matrix(0, npops, npts + 1)
  mordc  <- matrix(0, npops, npts + 1)
  morhcc <- matrix(0, npops, npts + 1)
  morlt  <- matrix(0, npops, npts + 1)
  morplt <- matrix(0, npops, npts + 1)
  l      <- matrix(0, npops, npts + 1)
  
  for (i in 1:npops) {
    morb[i,]  <- parama[, paste0("morb",  i)]
    mordc[i,] <- parama[, paste0("mordc", i)]
    morhcc[i,]<- parama[, paste0("morhcc",i)]
    morlt[i,] <- parama[, paste0("morlt", i)]
    morplt[i,]<- parama[, paste0("morplt",i)]
    l[i,]     <- parama[, paste0("leave", i)]
  }
  
  morb_dt   <- 1-(1-morb)^dt
  mordc_dt  <- 1-(1-mordc)^dt
  morhcc_dt <- 1-(1-morhcc)^dt
  morlt_dt  <- 1-(1-morlt)^dt
  morplt_dt <- 1-(1-morplt)^dt
  leave_dt  <- 1-(1-l)^dt
  
  mordcCure  <- matrix(0, npops, npts + 1)
  morhccCure <- matrix(0, npops, npts + 1)
  for (i in 1:npops) {
    mordcCure[i,]  <- mordc[i,]  * parama$Cure_mordc_Reduction
    morhccCure[i,] <- morhcc[i,] * (1 - parama$Cure_morhcc_Reduction)
  }
  mordcCure_dt  <- 1-(1-mordcCure)^dt
  morhccCure_dt <- 1-(1-morhccCure)^dt
  
  # ── Spontaneous clearance ─────────────────────────────────────────────────
  spc1 <- matrix(0, npops, npts + 1)
  for (i in 1:npops) spc1[i,] <- parama[, paste0("spc", i)]
  spc1_dt <- 1-(1-spc1)^dt
  
  # ── Disease progression ───────────────────────────────────────────────────
  transition    <- as.matrix(disease_progress)
  transition_dt <- 1-(1-transition)^dt
  
  # ── Fibrosis ──────────────────────────────────────────────────────────────
  fibprog    <- as.matrix(fib)
  fibprog_dt <- 1-(1-fibprog)^dt
  
  # ── Cascade dt conversions ────────────────────────────────────────────────
  tau_ab_dt   <- 1-(1-param_cascade$tau_ab)^dt
  tau_RNA_dt  <- 1-(1-param_cascade$tau_RNA)^dt
  tau_poct_dt <- 1-(1-param_cascade$tau_poct)^dt
  eta_dt      <- 1-(1-param_cascade$eta)^dt
  lota_dt     <- 1-(1-param_cascade$lota)^dt
  rho_dt      <- 1-(1-param_cascade$rho)^dt
  cure_dt     <- 1-(1-param_cascade$cured)^dt
  
  # SVR adjustments
  lota_dt[, c("f0","f1","f2","f3"), ]        <- lota_dt[, c("f0","f1","f2","f3"), ] * (1 - parama$SVR[1])
  lota_dt[, c("f4","dc","hcc","lt","plt"), ] <- lota_dt[, c("f4","dc","hcc","lt","plt"), ] * (1 - parama$SVRf4[1])
  cure_dt[, c("f0","f1","f2","f3"), ]        <- cure_dt[, c("f0","f1","f2","f3"), ] * parama$SVR[1]
  cure_dt[, c("f4","dc","hcc","lt","plt"), ] <- cure_dt[, c("f4","dc","hcc","lt","plt"), ] * parama$SVRf4[1]
  
  # Scenario cascade
  tau_ab_sc_dt   <- 1-(1-param_cascade_sc$tau_ab)^dt
  tau_RNA_sc_dt  <- 1-(1-param_cascade_sc$tau_RNA)^dt
  tau_poct_sc_dt <- 1-(1-param_cascade_sc$tau_poct)^dt
  eta_sc_dt      <- 1-(1-param_cascade_sc$eta)^dt
  
  # ── pop_array scaled by dt ────────────────────────────────────────────────
  pop_array_dt <- pop_array * dt
  
  # ── FOI matrix ────────────────────────────────────────────────────────────
  foi_mat <- matrix(0, npops, npts + 1)
  for (i in 1:npops)
    foi_mat[i,] <- parama[, paste0("beta", i)] * dt
  
  # ── Call C++ loop ─────────────────────────────────────────────────────────
  cpp_res <- hcvmsm_loop_cpp(
    init_pop       = initialPop,
    npts           = npts,
    morb_dt        = morb_dt,
    mordc_dt       = mordc_dt,
    morhcc_dt      = morhcc_dt,
    morlt_dt       = morlt_dt,
    morplt_dt      = morplt_dt,
    leave_dt       = leave_dt,
    mordcCure_dt   = mordcCure_dt,
    morhccCure_dt  = morhccCure_dt,
    spc1_dt        = spc1_dt,
    tau_ab_dt      = tau_ab_dt,
    tau_RNA_dt     = tau_RNA_dt,
    tau_poct_dt    = tau_poct_dt,
    eta_dt         = eta_dt,
    lota_dt        = lota_dt,
    rho_dt         = rho_dt,
    cure_dt        = cure_dt,
    tau_ab_sc_dt   = tau_ab_sc_dt,
    tau_RNA_sc_dt  = tau_RNA_sc_dt,
    tau_poct_sc_dt = tau_poct_sc_dt,
    eta_sc_dt      = eta_sc_dt,
    transition_dt  = transition_dt,
    fibprog_dt     = fibprog_dt,
    pop_array      = pop_array_dt,
    foi_dt         = foi_mat,
    fc             = fc,
    reinfP         = reinfP,
    is_POC_AU      = isTRUE(proj == "POC_AU")
  )
  
  dimnames(cpp_res$allPops) <- list(
    POC_AU$popNames,
    POC_AU$component_name,
    NULL
  )
  comp_names <- HCV$component_name
  pop_names  <- HCV$popNames
  
  dimnames(cpp_res$allPops) <- list(pop_names, comp_names, NULL)
  rownames(cpp_res$newS)              <- pop_names
  rownames(cpp_res$newInfections)     <- pop_names
  rownames(cpp_res$newHCVdeaths)      <- pop_names
  rownames(cpp_res$newTreatment)      <- pop_names
  rownames(cpp_res$newRetreat)        <- pop_names
  rownames(cpp_res$newCured)          <- pop_names
  rownames(cpp_res$newtreatfailed)    <- pop_names
  rownames(cpp_res$newreinfection)    <- pop_names
  rownames(cpp_res$newEntry)          <- pop_names
  rownames(cpp_res$newDeath)          <- pop_names
  rownames(cpp_res$newLeave)          <- pop_names
  rownames(cpp_res$newTestingAb_sc)      <- pop_names
  rownames(cpp_res$newTestingAg_sc)      <- pop_names
  rownames(cpp_res$newTestingPOCT_sc)    <- pop_names
  rownames(cpp_res$newTreatment_sc)      <- pop_names
  rownames(cpp_res$newTestingAb_sc_neg)  <- pop_names
  rownames(cpp_res$newTestingAg_sc_neg)  <- pop_names
  rownames(cpp_res$newTestingPOCT_sc_neg)<- pop_names
  # ── Stub zero matrices for unused outputs ─────────────────────────────────
  zero_mat   <- matrix(0, nrow = npops, ncol = npts)
  zero_arr3  <- array(0, c(npops, npops, npts))
  zero_arr3s <- array(0, c(1, 1, 1, 1))
  zero_arr_comp <- array(0, c(npops, ncomponent, npts))
  
  list(
    allPops              = cpp_res$allPops,
    newS                 = cpp_res$newS,
    newEntry             = cpp_res$newEntry,
    newDeath             = cpp_res$newDeath,
    newLeave             = cpp_res$newLeave,
    newInfections        = cpp_res$newInfections,
    newHCVdeaths         = cpp_res$newHCVdeaths,
    newTreatment         = cpp_res$newTreatment,
    newRetreat           = cpp_res$newRetreat,
    newTestingAb         = zero_mat,
    newTestingAg         = zero_mat,
    newTestingPOCT       = zero_mat,
    newTestingAb_sc      = cpp_res$newTestingAb_sc,
    newTestingAg_sc      = cpp_res$newTestingAg_sc,
    newTestingPOCT_sc    = cpp_res$newTestingPOCT_sc,
    newTreatment_sc      = cpp_res$newTreatment_sc,
    newTestingAb_neg     = zero_mat,
    newTestingAg_neg     = zero_mat,
    newTestingPOCT_neg   = zero_mat,
    newTestingAb_sc_neg  = cpp_res$newTestingAb_sc_neg,
    newTestingAg_sc_neg  = cpp_res$newTestingAg_sc_neg,
    newTestingPOCT_sc_neg= cpp_res$newTestingPOCT_sc_neg,
    newCured             = cpp_res$newCured,
    newtreatfailed       = cpp_res$newtreatfailed,
    newreinfection       = cpp_res$newreinfection,
    newreinfection_chronic = zero_mat,
    newpop_tran          = zero_arr3,
    newpop_tranState     = zero_arr3s,
    inflow               = zero_mat,
    inflow_hcv           = zero_mat,
    outflow              = zero_mat,
    outflow_hcv          = zero_mat,
    death_hcv            = zero_mat,
    HCVdeathState        = array(0, c(npops, 20, npts)),
    newDeathState        = zero_arr_comp
  )
}
