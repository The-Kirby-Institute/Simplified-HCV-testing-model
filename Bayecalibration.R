#### This is the script to work on bayesian calibration test #### 
# I will use the TWMSMHCV model to test the bayesian calibration 
# here's some thought to work toward the turn this model to bayesian calibration 


# Load necessary libraries

# 1. Decide which parameters to calibrate
# prevalence and incidence here 
# Define parameter vector
# Let 𝑃= HCV$npops.
# θ = [ beta_1 … beta_P,
#       tau_ab_1 … tau_ab_P,
#       tau_RNA_1 … tau_RNA_P,
#       eta_1 … eta_P ]
nParams <- 4 * P
paramNames <- c(
  paste0("beta_",   1:P),
  paste0("tau_ab_",  1:P),
  paste0("tau_RNA_", 1:P),
  paste0("eta_",     1:P)
)

# Load observation data
# Prevalence
obs_prev_years      <- c(2015, 2018, 2019, 2020, 2021, 2022)
obs_prev_overall    <- c(150, 140, 130, 125, 118, 110)  # positives
obs_prev_overall_N  <- c(1000,1000,1000,1000,1000,1000) # sample sizes
# one list per subpop i:
obs_prev_subpop     <- lapply(1:P, function(i) {
  list(years = obs_prev_years,
       pos   = c(...),   # replace with your subpop‐i positives
       N     = c(...))   # and sample sizes
})

# Incidence
obs_inc_years       <- c(2015, 2018, 2019, 2022)
obs_inc_overall     <- c(20,   18,   16,   13)   # new‐cases
# For the one subpop of interest:
obs_inc_subpop      <- list(
  years = obs_inc_years,
  counts = c(5, 4, 3, 2)          # replace with your data
)


# 3. Choose initial likelihood families
# prev --> Binominal 
# incidence --> Poisson 

# 4. Check dispersion before fitting
# Incidence
disp_inc <- var(obs_inc_overall) / mean(obs_inc_overall)

# Prevalence
props         <- obs_prev_overall / obs_prev_overall_N
var_prop      <- var(props)
p_bar         <- mean(props)
expected_var  <- mean(p_bar*(1-p_bar) / obs_prev_overall_N)

# If disp_inc ≫ 1, consider Negative-Binomial

# If var_prop ≫ expected_var, consider Beta-Binomial

# 5. Write the model wrapper
predict_HCV <- function(params) {
  P <- HCV$npops
  β    <- params[        1:P]
  τ_ab <- params[(P+1):(2*P)]
  τ_RN <- params[(2*P+1):(3*P)]
  η    <- params[(3*P+1):(4*P)]
  
  # Scale parameters in your lists...
  par2 <- parama; cas2 <- param_cascade
  for(i in 1:P) {
    par2$beta[[i]]           <- parama$beta[[i]]   * β[i]
    cas2$tau_ab[i,,]         <- param_cascade$tau_ab[i,,]  * τ_ab[i]
    cas2$tau_RNA[i,,]        <- param_cascade$tau_RNA[i,,] * τ_RN[i]
    cas2$eta[i,,]            <- param_cascade$eta[i,,]     * η[i]
  }
  
  # Run model to max(obs years)
  maxYear <- max(obs_prev_years, obs_inc_years)
  out <- HCVMSM(HCV, par2, initialPop, disease_progress,
                pop_array, cas2, fib, end_Y = maxYear, proj="YOURPROJ")
  
  # Map calendar year → time index
  idx_prev <- obs_prev_years - HCV$startYear + 1
  idx_inc  <- obs_inc_years  - HCV$startYear + 1
  
  # Compute outputs
  inf_states <- setdiff(dimnames(out$allPops)[[2]],
                        c("s", grep("_cured", dimnames(out$allPops)[[2]], value=TRUE)))
  prev_all <- sapply(idx_prev,  function(t) sum(out$allPops[,inf_states,t]) / sum(out$allPops[,,t]))
  prev_sub <- sapply(1:P, function(i)
    sapply(idx_prev, function(t) sum(out$allPops[i,inf_states,t]) / sum(out$allPops[i,,t]))
  )
  inc_all <- sapply(idx_inc,  function(t) sum(out$newInfections[, t]))
  k_sub   <- 1  # your subpopulation index for incidence
  inc_sub <- sapply(idx_inc,  function(t) out$newInfections[k_sub, t])
  
  list(prev_all = prev_all, prev_sub = prev_sub,
       inc_all  = inc_all,  inc_sub  = inc_sub)
}

# 6. Define the likelihood explorer
likelihood <- function(params) {
  pred <- predict_HCV(params)
  logL <- 0
  
  # --- Prevalence: Binomial or (alt) Beta-Binomial --- #
  for(j in seq_along(obs_prev_years)) {
    k <- obs_prev_overall[j]; N <- obs_prev_overall_N[j]
    p <- pred$prev_all[j]
    logL <- logL + dbinom(k, size=N, prob=p, log=TRUE)
    # If extra dispersion → uncomment:
    # alpha <- 5;  beta_param <- alpha*(1-p)/p
    # logL <- logL + dbetabinom.ab(k, size=N, alpha=alpha, beta=beta_param, log=TRUE)
  }
  # Subpop prevalences
  for(i in 1:P) {
    ds <- obs_prev_subpop[[i]]
    for(j in seq_along(ds$years)) {
      k <- ds$pos[j]; N <- ds$N[j]
      p <- pred$prev_sub[i,j]
      logL <- logL + dbinom(k, size=N, prob=p, log=TRUE)
      # alternative:
      # logL <- logL + dbetabinom.ab(k, size=N, alpha=alpha, beta=beta_param, log=TRUE)
    }
  }
  
  # --- Incidence: Poisson or (alt) Neg-Binomial --- #
  for(j in seq_along(obs_inc_years)) {
    y <- obs_inc_overall[j]; λ <- pred$inc_all[j]
    logL <- logL + dpois(y, lambda=λ, log=TRUE)
    # If overdispersed → uncomment:
    # size_nb <- 10
    # logL <- logL + dnbinom(y, size=size_nb, mu=λ, log=TRUE)
  }
  # Subpop incidence
  for(j in seq_along(obs_inc_subpop$years)) {
    y <- obs_inc_subpop$counts[j]; λ <- pred$inc_sub[j]
    logL <- logL + dpois(y, lambda=λ, log=TRUE)
    # alternative:
    # logL <- logL + dnbinom(y, size=size_nb, mu=λ, log=TRUE)
  }
  
  return(logL)
}

# 7. specify priors 
prior <- function(params) {
  βs   <- params[   1:P]
  τabs <- params[(P+1):(2*P)]
  τRNs ← params[(2*P+1):(3*P)]
  ηs   <- params[(3*P+1):(4*P)]
  if ( any(βs   < 0) || any(βs   > 5) ||
       any(τabs < 0) || any(τabs > 2) ||
       any(τRNs < 0) || any(τRNs > 2) ||
       any(ηs   < 0) || any(ηs   > 1) ) {
    return(-Inf)
  }
  return(0)
}

# 8. Run MCMC & diagnose
library(BayesianTools)
setup <- createBayesianSetup(likelihoodFunction=likelihood,
                             priorFunction     =prior,
                             names             =paramNames)
settings <- list(iterations=50000, nrChains=3)
outMCMC <- runMCMC(bayesianSetup=setup, sampler="DEzs", settings=settings)



# 9. Convergence diagnostics
print( getDiagnostics(outMCMC)$gelman.rubin )    # R̂ per parameter
print( getDiagnostics(outMCMC)$ess )              # effective sample size

# 10. Posterior summaries
post_samples <- getSample(outMCMC, coda = FALSE)  # matrix [nSamples × nParams]
theta_mean   <- colMeans(post_samples)
theta_med    <- apply(post_samples, 2, median)
theta_CI     <- apply(post_samples, 2, quantile, probs = c(0.025, 0.975))

# 3. Trace & density plots
plot(outMCMC)    # built-in: trace (top) and marginal density (bottom)

# 4. Pairwise scatter (to spot correlations)
pairs(post_samples[, 1:6], pch = 20, cex = 0.3)   # e.g. first 6 parameters


