# =============================================================================
# Param_cal  —  v3 (chain framework, no lambda)
#
# All NP parameters scale by Cov_np uniformly:
#   tau_ab_sc   = NP_tau_ab   × frac_reflex  × Cov_np
#   tau_poct_sc = NP_tau_poct × frac_immeRNA × Cov_np
#   tau_RNA_sc  = NP_tau_RNA  × Cov_np                  (reflex RNA, no frac)
#   eta_sc      = NP_eta      × Cov_np
#
# fp argument is the per-subpop Cov_np vector (length 5):
#   c(C_PWID, C_fPWID, P_PWID, P_fPWID, P_nPWID)
# Already in annual probability units — no monthly→annual conversion needed
# (the chain framework derives Cov_np from annual Tr, so already annual).
#
# Community (1:2) and Prison (3:5) both use STEP function (not ramp) since
# Cov_np is now an annual program parameter, not a per-month testing rate.
# =============================================================================
Param_cal <- function(pj, dlist, index, S_Yint, S_Yend, r_Yend, NPlst,
                      frac_testing = NULL, fp) {
  
  # ---------------------------------------------------------------------------
  # Time indices
  # ---------------------------------------------------------------------------
  SYpoint_int <- as.integer((S_Yint - pj$cabY) / pj$timestep + 1)
  SYpoint_end <- as.integer((S_Yend - pj$cabY) / pj$timestep)
  SY_leng     <- SYpoint_end - SYpoint_int + 1
  
  rYpoint_end <- as.integer((r_Yend - pj$cabY) / pj$timestep)
  rY_leng     <- rYpoint_end - SYpoint_end
  
  intVal <- dlist[[index]][, 3, (SYpoint_int - 1)]
  
  # fp is annual cumulative coverage under the chain framework
  fp_ann <- pmin(pmax(fp, 0), 0.999999)
  
  # ---------------------------------------------------------------------------
  # Conversion helpers
  # HCV_test appears to store annual probabilities and internally converts to
  # monthly. Therefore we store annual cumulative values here.
  # ---------------------------------------------------------------------------
  annual_to_monthly <- function(x) {
    x <- pmin(pmax(x, 0), 0.999999)
    1 - (1 - x)^pj$timestep
  }
  
  monthly_to_annual <- function(x) {
    x <- pmin(pmax(x, 0), 0.999999)
    1 - (1 - x)^(1 / pj$timestep)
  }
  
  # ---------------------------------------------------------------------------
  # Target annual value per population
  # ---------------------------------------------------------------------------
  scVal_dt <- numeric(pj$npops)
  
  if (index == "eta") {
    
    # Keep eta unchanged.
    scVal_dt[1:2] <- NPlst[["C"]][["eta"]] * fp_ann[1:2]
    scVal_dt[3:5] <- NPlst[["P"]][["eta"]] * fp_ann[3:5]
    
  } else if (index %in% c("tau_ab", "tau_poct", "tau_RNA")) {
    
    # Base annual parameter by population
    base_ann <- numeric(pj$npops)
    base_ann[1:2] <- NPlst[["C"]][[index]]
    base_ann[3:5] <- NPlst[["P"]][[index]]
    base_ann <- pmin(pmax(base_ann, 0), 0.999999)
    
    # First combine annual base parameter and annual coverage.
    # This is the key correction.
    total_ann <- pmin(base_ann * fp_ann, 0.999999)
    
    if (is.null(frac_testing)) {
      
      # No pathway split: preserve original annual total.
      scVal_dt <- total_ann
      
    } else {
      
      # Pathway fraction vector. These are fractions, not probabilities over time.
      frac_vec <- numeric(pj$npops)
      
      if (index == "tau_ab") {
        
        frac_vec[1:2] <- as.numeric(frac_testing[["C"]][["reflex"]])
        frac_vec[3:5] <- as.numeric(frac_testing[["P"]][["reflex"]])
        
      } else if (index == "tau_poct") {
        
        frac_vec[1:2] <- as.numeric(frac_testing[["C"]][["immeRNA"]])
        frac_vec[3:5] <- as.numeric(frac_testing[["P"]][["immeRNA"]])
        
      } else if (index == "tau_RNA") {
        
        # RNA here is reflex-pathway RNA.
        frac_vec[1:2] <- as.numeric(frac_testing[["C"]][["reflex"]])
        frac_vec[3:5] <- as.numeric(frac_testing[["P"]][["reflex"]])
      }
      
      frac_vec <- pmin(pmax(frac_vec, 0), 1)
      
      # Correct scale:
      #   total annual parameter -> monthly total rate
      #   monthly total rate * pathway fraction
      #   monthly pathway rate -> stored annual cumulative parameter
      total_monthly <- annual_to_monthly(total_ann)
      pathway_monthly <- total_monthly * frac_vec
      scVal_dt <- monthly_to_annual(pathway_monthly)
    }
    
  } else {
    
    stop(sprintf("Param_cal: unsupported index '%s'", index))
  }
  
  scVal_dt <- pmin(pmax(scVal_dt, 0), 0.999999)
  
  # ---------------------------------------------------------------------------
  # Fill parameter arrays — step function for all populations
  # ---------------------------------------------------------------------------
  for (i in 2:dim(dlist[[index]])[[2]]) {
    
    if (rYpoint_end > SYpoint_end) {
      
      for (pop in 1:pj$npops) {
        dlist[[index]][pop, i, SYpoint_int:pj$npts] <-
          c(
            rep(scVal_dt[pop], SY_leng + rY_leng),
            rep(intVal[pop],   pj$npts - rYpoint_end)
          )
      }
      
    } else {
      
      for (pop in 1:pj$npops) {
        dlist[[index]][pop, i, SYpoint_int:pj$npts] <-
          c(
            rep(scVal_dt[pop], SY_leng),
            rep(intVal[pop],   pj$npts - SYpoint_end)
          )
      }
    }
  }
  
  return(dlist[[index]])
}


# ─── carry_forward_prison ─────────────────────────────────────────────────────
# Unchanged — community resets to base, prison carries forward last value
carry_forward_prison <- function(dlist_np, dlist_base, index, b_pt, end_dt) {
  
  dim_length <- dim(dlist_np[[index]])[3]
  n_stages   <- dim(dlist_np[[index]])[2]
  
  dlist_np[[index]][1:2, , b_pt:dim_length] <-
    dlist_base[[index]][1:2, , b_pt:dim_length]
  
  for (s in 1:n_stages) {
    for (pop in 3:5) {
      dlist_np[[index]][pop, s, b_pt:dim_length] <-
        dlist_np[[index]][pop, s, end_dt]
    }
  }
  
  return(dlist_np[[index]])
}


# ─── scale_CT_eta / RNA / ab — unchanged from previous version ─────────────────
# All three use C_np vector of length 5; displacement = 1 - alpha * C_np
annual_to_monthly <- function(x, pj) {
  x <- pmin(pmax(x, 0), 0.999999999999999999)
  1 - (1 - x)^pj$timestep
}

monthly_to_annual <- function(x, pj) {
  x <- pmin(pmax(x, 0), 0.99999999999999999)
  1 - (1 - x)^(1 / pj$timestep)
}




scale_CT_eta <- function(dfList_CT, dfList_base, pj,
                         S_Yint, S_Yend, C_np, alpha = 1.0) {
  
  SYpoint_int <- as.integer((S_Yint - pj$cabY) / pj$timestep + 1)
  SYpoint_end <- as.integer((S_Yend - pj$cabY) / pj$timestep)
  
  SY_leng <- SYpoint_end - SYpoint_int + 1
  t_range <- SYpoint_int:SYpoint_end
  
  n_stages <- dim(dfList_base[["eta"]])[2]
  
  for (s in 1:n_stages) {
    
    for (pop in 1:5) {
      
      base_ann <- dfList_base[["eta"]][pop, s, SYpoint_int - 1]
      
      # -----------------------------------------------------------------------
      # Community (1:2)
      # Keep original annual-scale displacement.
      # -----------------------------------------------------------------------
      if (pop %in% 1:2) {
        
        scale <- pmax(0, 1 - alpha * C_np[pop])
        
        dfList_CT[["eta"]][pop, s, t_range] <-
          rep(base_ann * scale, SY_leng)
        
      } else {
        
        # ---------------------------------------------------------------------
        # Prison (3:5)
        # Displacement must happen on monthly scale.
        # ---------------------------------------------------------------------
        
        C_np_monthly <- annual_to_monthly(C_np[pop], pj)
        
        scale_monthly <- pmax(0, 1 - alpha * C_np_monthly)
        
        base_monthly <- annual_to_monthly(base_ann, pj)
        
        scaled_monthly <- base_monthly * scale_monthly
        
        scaled_ann <- monthly_to_annual(scaled_monthly, pj)
        
        dfList_CT[["eta"]][pop, s, t_range] <-
          rep(scaled_ann, SY_leng)
      }
    }
  }
  
  dfList_CT
}

scale_CT_RNA <- function(dfList_CT, dfList_base, pj,
                         S_Yint, S_Yend, C_np, alpha = 1.0) {
  
  SYpoint_int <- as.integer((S_Yint - pj$cabY) / pj$timestep + 1)
  SYpoint_end <- as.integer((S_Yend - pj$cabY) / pj$timestep)
  
  SY_leng <- SYpoint_end - SYpoint_int + 1
  t_range <- SYpoint_int:SYpoint_end
  
  n_stages <- dim(dfList_base[["tau_RNA"]])[2]
  
  for (s in 1:n_stages) {
    
    for (pop in 1:5) {
      
      base_ann <- dfList_base[["tau_RNA"]][pop, s, SYpoint_int - 1]
      
      # -----------------------------------------------------------------------
      # Community
      # -----------------------------------------------------------------------
      if (pop %in% 1:2) {
        
        scale <- pmax(0, 1 - alpha * C_np[pop])
        
        dfList_CT[["tau_RNA"]][pop, s, t_range] <-
          rep(base_ann * scale, SY_leng)
        
      } else {
        
        # ---------------------------------------------------------------------
        # Prison
        # ---------------------------------------------------------------------
        C_np_monthly <- annual_to_monthly(C_np[pop], pj)
        
        scale_monthly <- pmax(0, 1 - alpha * C_np_monthly)
        
        base_monthly <- annual_to_monthly(base_ann, pj)
        
        scaled_monthly <- base_monthly * scale_monthly
        
        scaled_ann <- monthly_to_annual(scaled_monthly, pj)
        
        dfList_CT[["tau_RNA"]][pop, s, t_range] <-
          rep(scaled_ann, SY_leng)
      }
    }
  }
  
  dfList_CT
}

scale_CT_ab <- function(dfList_CT, dfList_base, pj,
                        S_Yint, S_Yend, C_np, alpha = 1.0) {
  
  SYpoint_int <- as.integer((S_Yint - pj$cabY) / pj$timestep + 1)
  SYpoint_end <- as.integer((S_Yend - pj$cabY) / pj$timestep)
  
  SY_leng <- SYpoint_end - SYpoint_int + 1
  t_range <- SYpoint_int:SYpoint_end
  
  n_stages <- dim(dfList_base[["tau_ab"]])[2]
  
  for (s in 1:n_stages) {
    
    for (pop in 1:5) {
      
      base_ann <- dfList_base[["tau_ab"]][pop, s, SYpoint_int - 1]
      
      # -----------------------------------------------------------------------
      # Community
      # -----------------------------------------------------------------------
      if (pop %in% 1:2) {
        
        scale <- pmax(0, 1 - alpha * C_np[pop])
        
        dfList_CT[["tau_ab"]][pop, s, t_range] <-
          rep(base_ann * scale, SY_leng)
        
      } else {
        
        # ---------------------------------------------------------------------
        # Prison
        # ---------------------------------------------------------------------
        C_np_monthly <- annual_to_monthly(C_np[pop], pj)
        
        scale_monthly <- pmax(0, 1 - alpha * C_np_monthly)
        
        base_monthly <- annual_to_monthly(base_ann, pj)
        
        scaled_monthly <- base_monthly * scale_monthly
        
        scaled_ann <- monthly_to_annual(scaled_monthly, pj)
        
        dfList_CT[["tau_ab"]][pop, s, t_range] <-
          rep(scaled_ann, SY_leng)
      }
    }
  }
  
  dfList_CT
}

fs_estimate_v3 <- function(year, fc_list) {
  yr_str <- as.character(year)
  if (is.null(fc_list[[yr_str]])) {
    stop(sprintf("fc_list does not contain year %s — run compute_ccal_chain.R first.",
                 yr_str))
  }
  return(list(fc_list[[yr_str]]))
}


# =============================================================================
# COMPATIBILITY WRAPPER for old fs_estimate signature
#
# The original function had signature:
#   fs_estimate(num_ab, cov_np, frac_ab, fp, year, endY, modsim)
#
# All of these inputs are now redundant — fc is fully determined by data
# and model stocks. This wrapper accepts the old args and ignores them.
# =============================================================================
fs_estimate <- function(num_ab = NULL, cov_np = NULL, frac_ab = NULL,
                        fp = NULL, year, endY = NULL, modsim = NULL) {
  if (!exists("fc_list")) {
    stop("fc_list not found in workspace — run compute_ccal_chain.R before fs_estimate.")
  }
  fs_estimate_v3(year, fc_list)
}