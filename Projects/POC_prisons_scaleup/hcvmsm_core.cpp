#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// BASE_VAL: old pop + entry + net_flow - background_death - leave - hcv_death
#define BASE_VAL(j) (entry1(i,j) + oldPop(i,j) + net_flow(i,j) \
- death(i,j) - lv(i,j) - death_hcv(i,j))

// [[Rcpp::export]]
List hcvmsm_loop_cpp(
    arma::mat  init_pop,
    int        npts,
    
    arma::mat  morb_dt,
    arma::mat  mordc_dt,
    arma::mat  morhcc_dt,
    arma::mat  morlt_dt,
    arma::mat  morplt_dt,
    arma::mat  leave_dt,
    arma::mat  mordcCure_dt,
    arma::mat  morhccCure_dt,
    arma::mat  spc1_dt,
    
    arma::cube tau_ab_dt,
    arma::cube tau_RNA_dt,
    arma::cube tau_poct_dt,
    arma::cube eta_dt,
    arma::cube lota_dt,
    arma::cube rho_dt,
    arma::cube cure_dt,
    
    arma::cube tau_ab_sc_dt,
    arma::cube tau_RNA_sc_dt,
    arma::cube tau_poct_sc_dt,
    arma::cube eta_sc_dt,
    
    arma::mat  transition_dt,
    arma::mat  fibprog_dt,
    arma::cube pop_array,
    arma::mat  foi_dt,
    arma::mat  fc,
    arma::mat  reinfP,
    bool       is_POC_AU
) {
  int npops = init_pop.n_rows;
  int ncomp = init_pop.n_cols;
  
  // ── Compartment layout (0-based) ─────────────────────────────────────────
  // s=0
  // a:   ud=1,  dab=2,  drna=3,  tr=4,  trf=5,  cu=6
  // f0:  ud=7,  dab=8,  drna=9,  tr=10, trf=11, cu=12
  // f1:  ud=13, dab=14, drna=15, tr=16, trf=17, cu=18
  // f2:  ud=19, dab=20, drna=21, tr=22, trf=23, cu=24
  // f3:  ud=25, dab=26, drna=27, tr=28, trf=29, cu=30
  // f4:  ud=31, dab=32, drna=33, tr=34, trf=35, cu=36
  // dc:  ud=37, dab=38, drna=39, tr=40, trf=41, cu=42
  // hcc: ud=43, dab=44, drna=45, tr=46, trf=47, cu=48
  // lt:  ud=49, dab=50, drna=51, tr=52, trf=53, cu=54
  // plt: ud=55, dab=56, drna=57, tr=58, trf=59, cu=60
  
  // Progress index: a=0,f0=1,f1=2,f2=3,f3=4,f4=5,dc=6,hcc=7,lt=8,plt=9
  // Transition index:
  //   a_f0=0, f0_f1=1, f1_f2=2, f2_f3=3, f3_f4=4,
  //   f3_hcc=5, f4_dc=6, f4_hcc=7, dc_hcc=8, dc_lt=9, hcc_lt=10
  // Fibprog index:
  //   lt_plt=0, f3c_f4c=1, f3c_hccc=2, f4c_dcc=3, f4c_hccc=4,
  //   dcc_hccc=5, dcc_ltc=6, hccc_ltc=7, ltc_pltc=8
  
  // Compartment index arrays (progress stage order)
  int ud_idx[10]  = { 1,  7, 13, 19, 25, 31, 37, 43, 49, 55};
  int dab_idx[10] = { 2,  8, 14, 20, 26, 32, 38, 44, 50, 56};
  int drna_idx[10]= { 3,  9, 15, 21, 27, 33, 39, 45, 51, 57};
  int tr_idx[10]  = { 4, 10, 16, 22, 28, 34, 40, 46, 52, 58};
  int trf_idx[10] = { 5, 11, 17, 23, 29, 35, 41, 47, 53, 59};
  int cu_idx[10]  = { 6, 12, 18, 24, 30, 36, 42, 48, 54, 60};
  
  // ── Allocate outputs ─────────────────────────────────────────────────────
  arma::cube allPops(npops, ncomp, npts, fill::zeros);
  allPops.slice(0) = init_pop;
  
  arma::mat newS(npops,           npts, fill::zeros);
  arma::mat newEntry(npops,       npts, fill::zeros);
  arma::mat newDeath(npops,       npts, fill::zeros);
  arma::mat newLeave(npops,       npts, fill::zeros);
  arma::mat newInfections(npops,  npts, fill::zeros);
  arma::mat newHCVdeaths(npops,   npts, fill::zeros);
  arma::mat newTreatment(npops,   npts, fill::zeros);
  arma::mat newRetreat(npops,     npts, fill::zeros);
  arma::mat newCured(npops,       npts, fill::zeros);
  arma::mat newtreatfailed(npops, npts, fill::zeros);
  arma::mat newreinfection(npops, npts, fill::zeros);
  arma::mat newTestingAb_sc(npops,       npts, fill::zeros);
  arma::mat newTestingAg_sc(npops,       npts, fill::zeros);
  arma::mat newTestingPOCT_sc(npops,     npts, fill::zeros);
  arma::mat newTreatment_sc(npops,       npts, fill::zeros);
  arma::mat newTestingAb_sc_neg(npops,   npts, fill::zeros);
  arma::mat newTestingAg_sc_neg(npops,   npts, fill::zeros);
  arma::mat newTestingPOCT_sc_neg(npops, npts, fill::zeros);
  
  // ── Main time loop ───────────────────────────────────────────────────────
  // R loop: for(t in 2:npts) uses 1-based t, reads oldPop[,,t-1], writes newPop[,,t]
  // C++ equivalent: t=1..npts-1, reads slice(t-1), writes slice(t)
  for (int t = 1; t < npts; t++) {
    
    arma::mat oldPop = allPops.slice(t - 1);
    arma::mat newPop(npops, ncomp, fill::zeros);
    
    // ── Population sizes ─────────────────────────────────────────────────
    arma::vec N(npops, fill::zeros);
    arma::vec S(npops, fill::zeros);   // susceptible = s + all cured
    arma::vec I(npops, fill::zeros);
    
    for (int i = 0; i < npops; i++) {
      N(i) = arma::accu(oldPop.row(i));
      S(i) = oldPop(i, 0);                      // s
      for (int p = 0; p < 10; p++)
        S(i) += oldPop(i, cu_idx[p]);          // all cured stages
      I(i) = N(i) - S(i);
    }
    
    // ── Prevalence II ─────────────────────────────────────────────────────
    // Matches R: POC_AU uses grouped prevalence, else global
    arma::vec II(npops, fill::zeros);
    if (is_POC_AU) {
      double N12 = N(0)+N(1), I12 = I(0)+I(1);
      double N34 = N(2)+N(3), I34 = I(2)+I(3);
      II(0) = (N12 > 0) ? I12/N12 : 0.0;
      II(1) = II(0);
      II(2) = (N34 > 0) ? I34/N34 : 0.0;
      II(3) = II(2);
      II(4) = (N(4) > 0) ? I(4)/N(4) : 0.0;
    } else {
      double Nall = arma::accu(N), Iall = arma::accu(I);
      double IIall = (Nall > 0) ? Iall/Nall : 0.0;
      for (int i = 0; i < npops; i++) II(i) = IIall;
    }
    
    // foi_II = foi * II (pre-computed for each pop)
    arma::vec foi_II(npops);
    for (int i = 0; i < npops; i++)
      foi_II(i) = foi_dt(i, t) * II(i);
    
    // ── Background death and leave ────────────────────────────────────────
    arma::mat death(npops, ncomp, fill::zeros);
    arma::mat death_hcv(npops, ncomp, fill::zeros);
    arma::mat lv(npops, ncomp, fill::zeros);
    
    for (int i = 0; i < npops; i++) {
      for (int j = 0; j < ncomp; j++) {
        death(i, j) = morb_dt(i, t)  * oldPop(i, j);
        lv(i, j)    = leave_dt(i, t) * oldPop(i, j);
      }
      // HCV-specific excess deaths by stage
      // dc (37-42)
      death_hcv(i,37) = mordc_dt(i,t)      * oldPop(i,37);
      death_hcv(i,38) = mordc_dt(i,t)      * oldPop(i,38);
      death_hcv(i,39) = mordc_dt(i,t)      * oldPop(i,39);
      death_hcv(i,40) = mordc_dt(i,t)      * oldPop(i,40);
      death_hcv(i,41) = mordc_dt(i,t)      * oldPop(i,41);
      death_hcv(i,42) = mordcCure_dt(i,t)  * oldPop(i,42);
      // hcc (43-48)
      death_hcv(i,43) = morhcc_dt(i,t)     * oldPop(i,43);
      death_hcv(i,44) = morhcc_dt(i,t)     * oldPop(i,44);
      death_hcv(i,45) = morhcc_dt(i,t)     * oldPop(i,45);
      death_hcv(i,46) = morhcc_dt(i,t)     * oldPop(i,46);
      death_hcv(i,47) = morhcc_dt(i,t)     * oldPop(i,47);
      death_hcv(i,48) = morhccCure_dt(i,t) * oldPop(i,48);
      // lt (49-54)
      for (int j = 49; j <= 54; j++)
        death_hcv(i,j) = morlt_dt(i,t)  * oldPop(i,j);
      // plt (55-60)
      for (int j = 55; j <= 60; j++)
        death_hcv(i,j) = morplt_dt(i,t) * oldPop(i,j);
    }
    
    // ── Entry ─────────────────────────────────────────────────────────────
    // R (POC_AU): pops 1-4 deaths/leaves replenish pop1[s]; pop5 self-replenishes
    arma::mat entry1(npops, ncomp, fill::zeros);
    if (is_POC_AU) {
      double entry_s = 0.0;
      for (int i = 0; i < 4; i++) {
        entry_s += arma::accu(death.row(i));
        entry_s += arma::accu(death_hcv.row(i));
        entry_s += arma::accu(lv.row(i));
      }
      entry1(0, 0) = entry_s;
      double sum_pop4 = arma::accu(oldPop.row(4));
      if (sum_pop4 > 0) {
        double lv4_sum = arma::accu(lv.row(4));
        for (int j = 0; j < ncomp; j++)
          entry1(4, j) = lv4_sum * oldPop(4, j) / sum_pop4;
      }
    }
    // Non-POC_AU entry handled externally (entry_dt passed via pop_array logic)
    
    // ── Net population flow between sub-populations ───────────────────────
    // net_flow(i,j) = colSums(pa*oldPop[,j])[i] - oldPop(i,j)*rowSums(pa)[i]
    // Matches R: colSums(pop_array[,,t]*oldPop[,j]) - rowSums(pop_array[,,t]*oldPop[,j])
    arma::mat pa_t = pop_array.slice(t);           // npops x npops, already dt-transformed
    // R uses pop_array[,,t] with 1-based t=2..npts,
    // which equals 0-based slice t (C++ t=1..npts-1)
    arma::vec pa_rowsum = arma::sum(pa_t, 1);      // row sums
    arma::mat net_flow(npops, ncomp, fill::zeros);
    for (int j = 0; j < ncomp; j++) {
      arma::vec col_j    = oldPop.col(j);
      arma::vec inflow_j = pa_t.t() * col_j;    // colSums(pa * oldPop[,j])
      for (int i = 0; i < npops; i++)
        net_flow(i, j) = inflow_j(i) - oldPop(i, j) * pa_rowsum(i);
    }
    
    // ── Compartment equations ─────────────────────────────────────────────
    for (int i = 0; i < npops; i++) {
      
      // ── s (compartment 0) ─────────────────────────────────────────────
      newPop(i,0) = BASE_VAL(0)
      - foi_II(i)*oldPop(i,0);
      
      // ── Stage a (progress index p=0, base compartment b=1) ────────────
      {
        double sp   = spc1_dt(i,t);
        double ta   = tau_ab_dt(i,0,t);
        double tasc = tau_ab_sc_dt(i,0,t);
        double tp   = tau_poct_dt(i,0,t);
        double tpsc = tau_poct_sc_dt(i,0,t);
        double tr   = tau_RNA_dt(i,0,t);
        double trsc = tau_RNA_sc_dt(i,0,t);
        double e    = eta_dt(i,0,t);
        double esc  = eta_sc_dt(i,0,t);
        double lo   = lota_dt(i,0,t);
        double ro   = rho_dt(i,0,t);
        double cu0  = cure_dt(i,0,t);
        double td   = transition_dt(i,0);  // a_f0
        
        // a_undiag (1)
        newPop(i,1) = BASE_VAL(1)
          - td*oldPop(i,1)
          - sp*oldPop(i,1)
          - ta*(1.0-sp)*oldPop(i,1)
          - tasc*(1.0-sp)*oldPop(i,1)
          - tp*(1.0-sp)*oldPop(i,1)
          - tpsc*(1.0-sp)*oldPop(i,1)
          + foi_II(i)*oldPop(i,0)
          + reinfP(i,t)*foi_II(i)*oldPop(i,6)    // a_cured
          + reinfP(i,t)*foi_II(i)*oldPop(i,12);  // f0_cured
          
          // a_diag_ab (2)
          // FIX: sc RNA subtraction is -trsc*oldPop(i,2), NOT -trsc*tasc*undiag
          newPop(i,2) = BASE_VAL(2)
            - td*oldPop(i,2)
            + ta*(1.0-sp)*oldPop(i,1)
            - tr*oldPop(i,2)
            + tasc*(1.0-sp)*oldPop(i,1)
            - trsc*oldPop(i,2);
            
            // a_diag_RNA (3)
            // FIX: sc inflow = trsc*oldPop(diag_ab) + tpsc*(1-sp)*oldPop(undiag)
            //      sc outflow = -esc*oldPop(diag_RNA)
            newPop(i,3) = BASE_VAL(3)
              - td*oldPop(i,3)
              + tr*oldPop(i,2)
              + tp*(1.0-sp)*oldPop(i,1)
              - e*oldPop(i,3)
              + trsc*oldPop(i,2)
              + tpsc*(1.0-sp)*oldPop(i,1)
              - esc*oldPop(i,3);
              
              // a_treat (4)
              newPop(i,4) = BASE_VAL(4)
                - td*oldPop(i,4)
                + e*oldPop(i,3)
                + esc*oldPop(i,3)
                - lo*oldPop(i,4)
                - cu0*oldPop(i,4)
                + ro*oldPop(i,5);
                
                // a_treat_f (5)
                newPop(i,5) = BASE_VAL(5)
                  - td*oldPop(i,5)
                  - ro*oldPop(i,5)
                  + lo*oldPop(i,4);
                  
                  // a_cured (6)
                  newPop(i,6) = BASE_VAL(6)
                    + cu0*oldPop(i,4)
                    + sp*oldPop(i,1)
                    - reinfP(i,t)*foi_II(i)*oldPop(i,6);
      }
      
      // ── Generic fibrosis stage helper ─────────────────────────────────
      // Stages f0-f3 share the same structure. Differences:
      //   - f0: inflow from a, no reinfection into cured (cured[f0] loses to reinfection)
      //   - f1,f2: inflow from prior stage, reinfection into undiag from cured
      //   - f3: extra transition out (f3_hcc), cured has fibprog to f4c+hccc
      //
      // For stages f0..f3 we use generic code with flags for extras.
      
      // f0 (p=1, b=7, tin=a_f0=0, tout=f0_f1=1)
      {
        int b=7, p=1, tin=0, tout=1;
        int prev_ud=ud_idx[0], prev_dab=dab_idx[0],
                                               prev_drna=drna_idx[0], prev_tr=tr_idx[0],
                                                                                    prev_trf=trf_idx[0];
        double ta2   = tau_ab_dt(i,p,t);
        double tasc2 = tau_ab_sc_dt(i,p,t);
        double tp2   = tau_poct_dt(i,p,t);
        double tpsc2 = tau_poct_sc_dt(i,p,t);
        double tr2   = tau_RNA_dt(i,p,t);
        double trsc2 = tau_RNA_sc_dt(i,p,t);
        double e2    = eta_dt(i,p,t);
        double esc2  = eta_sc_dt(i,p,t);
        double lo2   = lota_dt(i,p,t);
        double ro2   = rho_dt(i,p,t);
        double cu2   = cure_dt(i,p,t);
        double td_in = transition_dt(i,tin);
        double td_out= transition_dt(i,tout);
        
        // f0_undiag (7): no reinfection inflow (f0_cured reinfections go to a_undiag)
        newPop(i,b) = BASE_VAL(b)
          + td_in*oldPop(i,prev_ud)
          - td_out*oldPop(i,b)
          - ta2*oldPop(i,b)
          - tp2*oldPop(i,b)
          - tasc2*oldPop(i,b)
          - tpsc2*oldPop(i,b);
          
          // f0_diag_ab (8)
          // FIX: sc RNA subtraction is -trsc2*oldPop(diag_ab)
          newPop(i,b+1) = BASE_VAL(b+1)
            + td_in*oldPop(i,prev_dab)
            - td_out*oldPop(i,b+1)
            + ta2*oldPop(i,b)
            - tr2*oldPop(i,b+1)
            + tasc2*oldPop(i,b)
            - trsc2*oldPop(i,b+1);
            
            // f0_diag_RNA (9)
            // FIX: sc inflow = trsc2*oldPop(diag_ab) + tpsc2*oldPop(undiag)
            newPop(i,b+2) = BASE_VAL(b+2)
              + td_in*oldPop(i,prev_drna)
              - td_out*oldPop(i,b+2)
              + tr2*oldPop(i,b+1)
              + tp2*oldPop(i,b)
              - e2*oldPop(i,b+2)
              + trsc2*oldPop(i,b+1)
              + tpsc2*oldPop(i,b)
              - esc2*oldPop(i,b+2);
              
              // f0_treat (10)
              newPop(i,b+3) = BASE_VAL(b+3)
                + td_in*oldPop(i,prev_tr)
                - td_out*oldPop(i,b+3)
                + e2*oldPop(i,b+2)
                + esc2*oldPop(i,b+2)
                - lo2*oldPop(i,b+3)
                - cu2*oldPop(i,b+3)
                + ro2*oldPop(i,b+4);
                
                // f0_treat_f (11)
                newPop(i,b+4) = BASE_VAL(b+4)
                  + td_in*oldPop(i,prev_trf)
                  - td_out*oldPop(i,b+4)
                  + lo2*oldPop(i,b+3)
                  - ro2*oldPop(i,b+4);
                  
                  // f0_cured (12): reinfection out (goes back to a_undiag)
                  newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    - reinfP(i,t)*foi_II(i)*oldPop(i,b+5);
      }
      
      // f1 (p=2, b=13, tin=f0_f1=1, tout=f1_f2=2)
      {
        int b=13, p=2, tin=1, tout=2;
        int prev_ud=ud_idx[1], prev_dab=dab_idx[1],
                                               prev_drna=drna_idx[1], prev_tr=tr_idx[1],
                                                                                    prev_trf=trf_idx[1];
        double ta2   = tau_ab_dt(i,p,t);
        double tasc2 = tau_ab_sc_dt(i,p,t);
        double tp2   = tau_poct_dt(i,p,t);
        double tpsc2 = tau_poct_sc_dt(i,p,t);
        double tr2   = tau_RNA_dt(i,p,t);
        double trsc2 = tau_RNA_sc_dt(i,p,t);
        double e2    = eta_dt(i,p,t);
        double esc2  = eta_sc_dt(i,p,t);
        double lo2   = lota_dt(i,p,t);
        double ro2   = rho_dt(i,p,t);
        double cu2   = cure_dt(i,p,t);
        double td_in = transition_dt(i,tin);
        double td_out= transition_dt(i,tout);
        
        // f1_undiag (13): reinfection inflow from f1_cured
        newPop(i,b) = BASE_VAL(b)
          + td_in*oldPop(i,prev_ud)
          - td_out*oldPop(i,b)
          - ta2*oldPop(i,b)
          - tp2*oldPop(i,b)
          - tasc2*oldPop(i,b)
          - tpsc2*oldPop(i,b)
          + reinfP(i,t)*foi_II(i)*oldPop(i,cu_idx[2]);  // f1_cured=18
          
          newPop(i,b+1) = BASE_VAL(b+1)
            + td_in*oldPop(i,prev_dab)
            - td_out*oldPop(i,b+1)
            + ta2*oldPop(i,b)
            - tr2*oldPop(i,b+1)
            + tasc2*oldPop(i,b)
            - trsc2*oldPop(i,b+1);
            
            newPop(i,b+2) = BASE_VAL(b+2)
              + td_in*oldPop(i,prev_drna)
              - td_out*oldPop(i,b+2)
              + tr2*oldPop(i,b+1)
              + tp2*oldPop(i,b)
              - e2*oldPop(i,b+2)
              + trsc2*oldPop(i,b+1)
              + tpsc2*oldPop(i,b)
              - esc2*oldPop(i,b+2);
              
              newPop(i,b+3) = BASE_VAL(b+3)
                + td_in*oldPop(i,prev_tr)
                - td_out*oldPop(i,b+3)
                + e2*oldPop(i,b+2)
                + esc2*oldPop(i,b+2)
                - lo2*oldPop(i,b+3)
                - cu2*oldPop(i,b+3)
                + ro2*oldPop(i,b+4);
                
                newPop(i,b+4) = BASE_VAL(b+4)
                  + td_in*oldPop(i,prev_trf)
                  - td_out*oldPop(i,b+4)
                  + lo2*oldPop(i,b+3)
                  - ro2*oldPop(i,b+4);
                  
                  newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    - reinfP(i,t)*foi_II(i)*oldPop(i,b+5);
      }
      
      // f2 (p=3, b=19, tin=f1_f2=2, tout=f2_f3=3)
      {
        int b=19, p=3, tin=2, tout=3;
        int prev_ud=ud_idx[2], prev_dab=dab_idx[2],
                                               prev_drna=drna_idx[2], prev_tr=tr_idx[2],
                                                                                    prev_trf=trf_idx[2];
        double ta2   = tau_ab_dt(i,p,t);
        double tasc2 = tau_ab_sc_dt(i,p,t);
        double tp2   = tau_poct_dt(i,p,t);
        double tpsc2 = tau_poct_sc_dt(i,p,t);
        double tr2   = tau_RNA_dt(i,p,t);
        double trsc2 = tau_RNA_sc_dt(i,p,t);
        double e2    = eta_dt(i,p,t);
        double esc2  = eta_sc_dt(i,p,t);
        double lo2   = lota_dt(i,p,t);
        double ro2   = rho_dt(i,p,t);
        double cu2   = cure_dt(i,p,t);
        double td_in = transition_dt(i,tin);
        double td_out= transition_dt(i,tout);
        
        // f2_undiag (19): reinfection inflow from f2_cured
        newPop(i,b) = BASE_VAL(b)
          + td_in*oldPop(i,prev_ud)
          - td_out*oldPop(i,b)
          - ta2*oldPop(i,b)
          - tp2*oldPop(i,b)
          - tasc2*oldPop(i,b)
          - tpsc2*oldPop(i,b)
          + reinfP(i,t)*foi_II(i)*oldPop(i,cu_idx[3]);  // f2_cured=24
          
          newPop(i,b+1) = BASE_VAL(b+1)
            + td_in*oldPop(i,prev_dab)
            - td_out*oldPop(i,b+1)
            + ta2*oldPop(i,b)
            - tr2*oldPop(i,b+1)
            + tasc2*oldPop(i,b)
            - trsc2*oldPop(i,b+1);
            
            newPop(i,b+2) = BASE_VAL(b+2)
              + td_in*oldPop(i,prev_drna)
              - td_out*oldPop(i,b+2)
              + tr2*oldPop(i,b+1)
              + tp2*oldPop(i,b)
              - e2*oldPop(i,b+2)
              + trsc2*oldPop(i,b+1)
              + tpsc2*oldPop(i,b)
              - esc2*oldPop(i,b+2);
              
              newPop(i,b+3) = BASE_VAL(b+3)
                + td_in*oldPop(i,prev_tr)
                - td_out*oldPop(i,b+3)
                + e2*oldPop(i,b+2)
                + esc2*oldPop(i,b+2)
                - lo2*oldPop(i,b+3)
                - cu2*oldPop(i,b+3)
                + ro2*oldPop(i,b+4);
                
                newPop(i,b+4) = BASE_VAL(b+4)
                  + td_in*oldPop(i,prev_trf)
                  - td_out*oldPop(i,b+4)
                  + lo2*oldPop(i,b+3)
                  - ro2*oldPop(i,b+4);
                  
                  newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    - reinfP(i,t)*foi_II(i)*oldPop(i,b+5);
      }
      
      // f3 (p=4, b=25, tin=f2_f3=3, tout=f3_f4=4, extra out: f3_hcc=5)
      {
        int b=25, p=4, tin=3, tout=4, thcc=5;
        int prev_ud=ud_idx[3], prev_dab=dab_idx[3],
                                               prev_drna=drna_idx[3], prev_tr=tr_idx[3],
                                                                                    prev_trf=trf_idx[3];
        double ta2   = tau_ab_dt(i,p,t);
        double tasc2 = tau_ab_sc_dt(i,p,t);
        double tp2   = tau_poct_dt(i,p,t);
        double tpsc2 = tau_poct_sc_dt(i,p,t);
        double tr2   = tau_RNA_dt(i,p,t);
        double trsc2 = tau_RNA_sc_dt(i,p,t);
        double e2    = eta_dt(i,p,t);
        double esc2  = eta_sc_dt(i,p,t);
        double lo2   = lota_dt(i,p,t);
        double ro2   = rho_dt(i,p,t);
        double cu2   = cure_dt(i,p,t);
        double td_in  = transition_dt(i,tin);
        double td_out = transition_dt(i,tout);
        double td_hcc = transition_dt(i,thcc);
        
        // f3_undiag (25): reinfection inflow from f3_cured
        newPop(i,b) = BASE_VAL(b)
          + td_in*oldPop(i,prev_ud)
          - td_out*oldPop(i,b)
          - td_hcc*oldPop(i,b)
          - ta2*oldPop(i,b)
          - tp2*oldPop(i,b)
          - tasc2*oldPop(i,b)
          - tpsc2*oldPop(i,b)
          + reinfP(i,t)*foi_II(i)*oldPop(i,cu_idx[4]);  // f3_cured=30
          
          newPop(i,b+1) = BASE_VAL(b+1)
            + td_in*oldPop(i,prev_dab)
            - td_out*oldPop(i,b+1)
            - td_hcc*oldPop(i,b+1)
            + ta2*oldPop(i,b)
            - tr2*oldPop(i,b+1)
            + tasc2*oldPop(i,b)
            - trsc2*oldPop(i,b+1);
            
            newPop(i,b+2) = BASE_VAL(b+2)
              + td_in*oldPop(i,prev_drna)
              - td_out*oldPop(i,b+2)
              - td_hcc*oldPop(i,b+2)
              + tr2*oldPop(i,b+1)
              + tp2*oldPop(i,b)
              - e2*oldPop(i,b+2)
              + trsc2*oldPop(i,b+1)
              + tpsc2*oldPop(i,b)
              - esc2*oldPop(i,b+2);
              
              newPop(i,b+3) = BASE_VAL(b+3)
                + td_in*oldPop(i,prev_tr)
                - td_out*oldPop(i,b+3)
                - td_hcc*oldPop(i,b+3)
                + e2*oldPop(i,b+2)
                + esc2*oldPop(i,b+2)
                - lo2*oldPop(i,b+3)
                - cu2*oldPop(i,b+3)
                + ro2*oldPop(i,b+4);
                
                newPop(i,b+4) = BASE_VAL(b+4)
                  + td_in*oldPop(i,prev_trf)
                  - td_out*oldPop(i,b+4)
                  - td_hcc*oldPop(i,b+4)
                  + lo2*oldPop(i,b+3)
                  - ro2*oldPop(i,b+4);
                  
                  // f3_cured (30): fibprog to f4_cured and hcc_cured
                  newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    - fibprog_dt(i,1)*oldPop(i,b+5)  // f3c→f4c
                    - fibprog_dt(i,2)*oldPop(i,b+5)  // f3c→hccc
                    - reinfP(i,t)*foi_II(i)*oldPop(i,b+5);
      }
      
      // f4 (p=5, b=31, tin=f3_f4=4, out: f4_dc=6, f4_hcc=7)
      // No reinfection into f4_cured (0* in R)
      {
        int b=31, p=5, tin=4, tdc=6, thcc=7;
        int prev_ud=ud_idx[4], prev_dab=dab_idx[4],
                                               prev_drna=drna_idx[4], prev_tr=tr_idx[4],
                                                                                    prev_trf=trf_idx[4];
        double ta2   = tau_ab_dt(i,p,t);
        double tasc2 = tau_ab_sc_dt(i,p,t);
        double tp2   = tau_poct_dt(i,p,t);
        double tpsc2 = tau_poct_sc_dt(i,p,t);
        double tr2   = tau_RNA_dt(i,p,t);
        double trsc2 = tau_RNA_sc_dt(i,p,t);
        double e2    = eta_dt(i,p,t);
        double esc2  = eta_sc_dt(i,p,t);
        double lo2   = lota_dt(i,p,t);
        double ro2   = rho_dt(i,p,t);
        double cu2   = cure_dt(i,p,t);
        double td_in = transition_dt(i,tin);
        double td_dc = transition_dt(i,tdc);
        double td_hcc= transition_dt(i,thcc);
        
        newPop(i,b) = BASE_VAL(b)
          + td_in*oldPop(i,prev_ud)
          - td_dc*oldPop(i,b)
          - td_hcc*oldPop(i,b)
          - ta2*oldPop(i,b)
          - tp2*oldPop(i,b)
          - tasc2*oldPop(i,b)
          - tpsc2*oldPop(i,b);
          
          newPop(i,b+1) = BASE_VAL(b+1)
            + td_in*oldPop(i,prev_dab)
            - td_dc*oldPop(i,b+1)
            - td_hcc*oldPop(i,b+1)
            + ta2*oldPop(i,b)
            - tr2*oldPop(i,b+1)
            + tasc2*oldPop(i,b)
            - trsc2*oldPop(i,b+1);
            
            newPop(i,b+2) = BASE_VAL(b+2)
              + td_in*oldPop(i,prev_drna)
              - td_dc*oldPop(i,b+2)
              - td_hcc*oldPop(i,b+2)
              + tr2*oldPop(i,b+1)
              + tp2*oldPop(i,b)
              - e2*oldPop(i,b+2)
              + trsc2*oldPop(i,b+1)
              + tpsc2*oldPop(i,b)
              - esc2*oldPop(i,b+2);
              
              newPop(i,b+3) = BASE_VAL(b+3)
                + td_in*oldPop(i,prev_tr)
                - td_dc*oldPop(i,b+3)
                - td_hcc*oldPop(i,b+3)
                + e2*oldPop(i,b+2)
                + esc2*oldPop(i,b+2)
                - lo2*oldPop(i,b+3)
                - cu2*oldPop(i,b+3)
                + ro2*oldPop(i,b+4);
                
                newPop(i,b+4) = BASE_VAL(b+4)
                  + td_in*oldPop(i,prev_trf)
                  - td_dc*oldPop(i,b+4)
                  - td_hcc*oldPop(i,b+4)
                  + lo2*oldPop(i,b+3)
                  - ro2*oldPop(i,b+4);
                  
                  // f4_cured (36): receives f3c, loses to dcc+hccc; no reinfection
                  newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    + fibprog_dt(i,1)*oldPop(i,30)   // f3c→f4c
                    - fibprog_dt(i,3)*oldPop(i,b+5)  // f4c→dcc
                    - fibprog_dt(i,4)*oldPop(i,b+5); // f4c→hccc
      }
      
      // dc (p=6, b=37, in: f4_dc=6, out: dc_hcc=8, dc_lt=9)
      // No reinfection into dc_cured (0* in R)
      {
        int b=37, p=6, tdc=6, thcc=8, tlt=9;
        int prev_ud=ud_idx[5], prev_dab=dab_idx[5],
                                               prev_drna=drna_idx[5], prev_tr=tr_idx[5],
                                                                                    prev_trf=trf_idx[5];
        double ta2   = tau_ab_dt(i,p,t);
        double tasc2 = tau_ab_sc_dt(i,p,t);
        double tp2   = tau_poct_dt(i,p,t);
        double tpsc2 = tau_poct_sc_dt(i,p,t);
        double tr2   = tau_RNA_dt(i,p,t);
        double trsc2 = tau_RNA_sc_dt(i,p,t);
        double e2    = eta_dt(i,p,t);
        double esc2  = eta_sc_dt(i,p,t);
        double lo2   = lota_dt(i,p,t);
        double ro2   = rho_dt(i,p,t);
        double cu2   = cure_dt(i,p,t);
        double td_dc = transition_dt(i,tdc);
        double td_hcc= transition_dt(i,thcc);
        double td_lt = transition_dt(i,tlt);
        
        newPop(i,b) = BASE_VAL(b)
          + td_dc*oldPop(i,prev_ud)
          - td_hcc*oldPop(i,b)
          - td_lt*oldPop(i,b)
          - ta2*oldPop(i,b)
          - tp2*oldPop(i,b)
          - tasc2*oldPop(i,b)
          - tpsc2*oldPop(i,b);
          
          newPop(i,b+1) = BASE_VAL(b+1)
            + td_dc*oldPop(i,prev_dab)
            - td_hcc*oldPop(i,b+1)
            - td_lt*oldPop(i,b+1)
            + ta2*oldPop(i,b)
            - tr2*oldPop(i,b+1)
            + tasc2*oldPop(i,b)
            - trsc2*oldPop(i,b+1);
            
            newPop(i,b+2) = BASE_VAL(b+2)
              + td_dc*oldPop(i,prev_drna)
              - td_hcc*oldPop(i,b+2)
              - td_lt*oldPop(i,b+2)
              + tr2*oldPop(i,b+1)
              + tp2*oldPop(i,b)
              - e2*oldPop(i,b+2)
              + trsc2*oldPop(i,b+1)
              + tpsc2*oldPop(i,b)
              - esc2*oldPop(i,b+2);
              
              newPop(i,b+3) = BASE_VAL(b+3)
                + td_dc*oldPop(i,prev_tr)
                - td_hcc*oldPop(i,b+3)
                - td_lt*oldPop(i,b+3)
                + e2*oldPop(i,b+2)
                + esc2*oldPop(i,b+2)
                - lo2*oldPop(i,b+3)
                - cu2*oldPop(i,b+3)
                + ro2*oldPop(i,b+4);
                
                newPop(i,b+4) = BASE_VAL(b+4)
                  + td_dc*oldPop(i,prev_trf)
                  - td_hcc*oldPop(i,b+4)
                  - td_lt*oldPop(i,b+4)
                  + lo2*oldPop(i,b+3)
                  - ro2*oldPop(i,b+4);
                  
                  // dc_cured (42): receives f4c, loses to hccc+ltc; no reinfection
                  newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    + fibprog_dt(i,3)*oldPop(i,36)   // f4c→dcc
                    - fibprog_dt(i,5)*oldPop(i,b+5)  // dcc→hccc
                    - fibprog_dt(i,6)*oldPop(i,b+5); // dcc→ltc
      }
      
      // hcc (p=7, b=43, in: f4_hcc=7, dc_hcc=8, f3_hcc=5, out: hcc_lt=10)
      // No reinfection into hcc_cured (0* in R)
      {
        int b=43, p=7, tf4h=7, tdch=8, tf3h=5, thlt=10;
        int prev_f4_ud=ud_idx[5], prev_dc_ud=ud_idx[6], prev_f3_ud=ud_idx[4];
        double ta2   = tau_ab_dt(i,p,t);
        double tasc2 = tau_ab_sc_dt(i,p,t);
        double tp2   = tau_poct_dt(i,p,t);
        double tpsc2 = tau_poct_sc_dt(i,p,t);
        double tr2   = tau_RNA_dt(i,p,t);
        double trsc2 = tau_RNA_sc_dt(i,p,t);
        double e2    = eta_dt(i,p,t);
        double esc2  = eta_sc_dt(i,p,t);
        double lo2   = lota_dt(i,p,t);
        double ro2   = rho_dt(i,p,t);
        double cu2   = cure_dt(i,p,t);
        
        newPop(i,b) = BASE_VAL(b)
          + transition_dt(i,tf4h)*oldPop(i,prev_f4_ud)
          + transition_dt(i,tdch)*oldPop(i,prev_dc_ud)
          + transition_dt(i,tf3h)*oldPop(i,prev_f3_ud)
          - transition_dt(i,thlt)*oldPop(i,b)
          - ta2*oldPop(i,b)
          - tp2*oldPop(i,b)
          - tasc2*oldPop(i,b)
          - tpsc2*oldPop(i,b);
          
          newPop(i,b+1) = BASE_VAL(b+1)
            + transition_dt(i,tf4h)*oldPop(i,dab_idx[5])
            + transition_dt(i,tdch)*oldPop(i,dab_idx[6])
            + transition_dt(i,tf3h)*oldPop(i,dab_idx[4])
            - transition_dt(i,thlt)*oldPop(i,b+1)
            + ta2*oldPop(i,b)
            - tr2*oldPop(i,b+1)
            + tasc2*oldPop(i,b)
            - trsc2*oldPop(i,b+1);
            
            newPop(i,b+2) = BASE_VAL(b+2)
              + transition_dt(i,tf4h)*oldPop(i,drna_idx[5])
              + transition_dt(i,tdch)*oldPop(i,drna_idx[6])
              + transition_dt(i,tf3h)*oldPop(i,drna_idx[4])
              - transition_dt(i,thlt)*oldPop(i,b+2)
              + tr2*oldPop(i,b+1)
              + tp2*oldPop(i,b)
              - e2*oldPop(i,b+2)
              + trsc2*oldPop(i,b+1)
              + tpsc2*oldPop(i,b)
              - esc2*oldPop(i,b+2);
              
              newPop(i,b+3) = BASE_VAL(b+3)
                + transition_dt(i,tf4h)*oldPop(i,tr_idx[5])
                + transition_dt(i,tdch)*oldPop(i,tr_idx[6])
                + transition_dt(i,tf3h)*oldPop(i,tr_idx[4])
                - transition_dt(i,thlt)*oldPop(i,b+3)
                + e2*oldPop(i,b+2)
                + esc2*oldPop(i,b+2)
                - lo2*oldPop(i,b+3)
                - cu2*oldPop(i,b+3)
                + ro2*oldPop(i,b+4);
                
                newPop(i,b+4) = BASE_VAL(b+4)
                  + transition_dt(i,tf4h)*oldPop(i,trf_idx[5])
                  + transition_dt(i,tf3h)*oldPop(i,trf_idx[4])
                  + transition_dt(i,tdch)*oldPop(i,trf_idx[6])
                  - transition_dt(i,thlt)*oldPop(i,b+4)
                  + lo2*oldPop(i,b+3)
                  - ro2*oldPop(i,b+4);
                  
                  // hcc_cured (48): receives f4c+dcc+f3c, loses to ltc; no reinfection
                  newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    + fibprog_dt(i,4)*oldPop(i,36)   // f4c→hccc
                    + fibprog_dt(i,5)*oldPop(i,42)   // dcc→hccc
                    + fibprog_dt(i,2)*oldPop(i,30)   // f3c→hccc
                    - fibprog_dt(i,7)*oldPop(i,b+5); // hccc→ltc
      }
      
      // lt (p=8, b=49, in: hcc_lt=10, dc_lt=9, out: lt_plt via fibprog)
      // No reinfection into lt_cured (0* in R)
      {
        int b=49, p=8, thlt=10, tdlt=9;
        double ta2   = tau_ab_dt(i,p,t);
        double tasc2 = tau_ab_sc_dt(i,p,t);
        double tp2   = tau_poct_dt(i,p,t);
        double tpsc2 = tau_poct_sc_dt(i,p,t);
        double tr2   = tau_RNA_dt(i,p,t);
        double trsc2 = tau_RNA_sc_dt(i,p,t);
        double e2    = eta_dt(i,p,t);
        double esc2  = eta_sc_dt(i,p,t);
        double lo2   = lota_dt(i,p,t);
        double ro2   = rho_dt(i,p,t);
        double cu2   = cure_dt(i,p,t);
        double fp_lt = fibprog_dt(i,0);  // lt_plt
        
        newPop(i,b) = BASE_VAL(b)
          + transition_dt(i,thlt)*oldPop(i,ud_idx[7])   // hcc_ud
          + transition_dt(i,tdlt)*oldPop(i,ud_idx[6])   // dc_ud
          - fp_lt*oldPop(i,b)
          - ta2*oldPop(i,b)
          - tp2*oldPop(i,b)
          - tasc2*oldPop(i,b)
          - tpsc2*oldPop(i,b);
          
          newPop(i,b+1) = BASE_VAL(b+1)
            + transition_dt(i,thlt)*oldPop(i,dab_idx[7])
            + transition_dt(i,tdlt)*oldPop(i,dab_idx[6])
            - fp_lt*oldPop(i,b+1)
            + ta2*oldPop(i,b)
            - tr2*oldPop(i,b+1)
            + tasc2*oldPop(i,b)
            - trsc2*oldPop(i,b+1);
            
            newPop(i,b+2) = BASE_VAL(b+2)
              + transition_dt(i,thlt)*oldPop(i,drna_idx[7])
              + transition_dt(i,tdlt)*oldPop(i,drna_idx[6])
              - fp_lt*oldPop(i,b+2)
              + tr2*oldPop(i,b+1)
              + tp2*oldPop(i,b)
              - e2*oldPop(i,b+2)
              + trsc2*oldPop(i,b+1)
              + tpsc2*oldPop(i,b)
              - esc2*oldPop(i,b+2);
              
              newPop(i,b+3) = BASE_VAL(b+3)
                + transition_dt(i,thlt)*oldPop(i,tr_idx[7])
                + transition_dt(i,tdlt)*oldPop(i,tr_idx[6])
                + e2*oldPop(i,b+2)
                + esc2*oldPop(i,b+2)
                - fp_lt*oldPop(i,b+3)
                - lo2*oldPop(i,b+3)
                - cu2*oldPop(i,b+3)
                + ro2*oldPop(i,b+4);
                
                newPop(i,b+4) = BASE_VAL(b+4)
                  + transition_dt(i,thlt)*oldPop(i,trf_idx[7])
                  + transition_dt(i,tdlt)*oldPop(i,trf_idx[6])
                  + lo2*oldPop(i,b+3)
                  - fp_lt*oldPop(i,b+4)
                  - ro2*oldPop(i,b+4);
                  
                  // lt_cured (54): receives hccc+dcc, loses to pltc; no reinfection
                  newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    + fibprog_dt(i,7)*oldPop(i,48)   // hccc→ltc
                    + fibprog_dt(i,6)*oldPop(i,42)   // dcc→ltc
                    - fibprog_dt(i,8)*oldPop(i,b+5); // ltc→pltc
      }
      
      // plt (p=9, b=55, in: lt via fibprog)
      // No reinfection into plt_cured (0* in R)
      {
        int b=55, p=9;
        double ta2   = tau_ab_dt(i,p,t);
        double tasc2 = tau_ab_sc_dt(i,p,t);
        double tp2   = tau_poct_dt(i,p,t);
        double tpsc2 = tau_poct_sc_dt(i,p,t);
        double tr2   = tau_RNA_dt(i,p,t);
        double trsc2 = tau_RNA_sc_dt(i,p,t);
        double e2    = eta_dt(i,p,t);
        double esc2  = eta_sc_dt(i,p,t);
        double lo2   = lota_dt(i,p,t);
        double ro2   = rho_dt(i,p,t);
        double cu2   = cure_dt(i,p,t);
        double fp_lt = fibprog_dt(i,0);  // lt_plt
        
        newPop(i,b) = BASE_VAL(b)
          + fp_lt*oldPop(i,49)             // lt_ud→plt_ud
          - ta2*oldPop(i,b)
          - tp2*oldPop(i,b)
          - tasc2*oldPop(i,b)
          - tpsc2*oldPop(i,b);
          
          newPop(i,b+1) = BASE_VAL(b+1)
            + fp_lt*oldPop(i,50)
            + ta2*oldPop(i,b)
            - tr2*oldPop(i,b+1)
            + tasc2*oldPop(i,b)
            - trsc2*oldPop(i,b+1);
            
            newPop(i,b+2) = BASE_VAL(b+2)
              + fp_lt*oldPop(i,51)
              + tr2*oldPop(i,b+1)
              + tp2*oldPop(i,b)
              - e2*oldPop(i,b+2)
              + trsc2*oldPop(i,b+1)
              + tpsc2*oldPop(i,b)
              - esc2*oldPop(i,b+2);
              
              newPop(i,b+3) = BASE_VAL(b+3)
                + fp_lt*oldPop(i,52)
                + e2*oldPop(i,b+2)
                + esc2*oldPop(i,b+2)
                - lo2*oldPop(i,b+3)
                - cu2*oldPop(i,b+3)
                + ro2*oldPop(i,b+4);
                
                newPop(i,b+4) = BASE_VAL(b+4)
                  + fp_lt*oldPop(i,53)
                  + lo2*oldPop(i,b+3)
                  - ro2*oldPop(i,b+4);
                  
                  // plt_cured (60): receives ltc; no reinfection
                  newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    + fibprog_dt(i,8)*oldPop(i,54);  // ltc→pltc
      }
      
    } // end for each pop i
    
    // Floor at zero to match R: newPop[newPop < 0] <- 0
    newPop = arma::clamp(newPop, 0.0, arma::datum::inf);
    allPops.slice(t) = newPop;
    
    // ── Result aggregates ─────────────────────────────────────────────────
    for (int i = 0; i < npops; i++) {
      
      newS(i,t) = newPop(i,0);
      
      // Infections: primary + reinfection (a-f3 cured only, 0* for f4+)
      newInfections(i,t) = foi_II(i)*oldPop(i,0)
        + reinfP(i,t)*foi_II(i)*(
            oldPop(i,6) + oldPop(i,12) + oldPop(i,18)
        + oldPop(i,24) + oldPop(i,30));
      
      // HCV deaths (dc+hcc+lt+plt stages only)
      newHCVdeaths(i,t) =
      mordc_dt(i,t)*(oldPop(i,37)+oldPop(i,38)+oldPop(i,39)
                       +oldPop(i,40)+oldPop(i,41))
        + mordcCure_dt(i,t)*oldPop(i,42)
        + morhcc_dt(i,t)*(oldPop(i,43)+oldPop(i,44)+oldPop(i,45)
        +oldPop(i,46)+oldPop(i,47))
        + morhccCure_dt(i,t)*oldPop(i,48)
        + morlt_dt(i,t)*(oldPop(i,49)+oldPop(i,50)+oldPop(i,51)
        +oldPop(i,52)+oldPop(i,53)+oldPop(i,54))
        + morplt_dt(i,t)*(oldPop(i,55)+oldPop(i,56)+oldPop(i,57)
        +oldPop(i,58)+oldPop(i,59)+oldPop(i,60));
        
        newEntry(i,t) = arma::accu(entry1.row(i));
        newDeath(i,t) = arma::accu(death.row(i));
        newLeave(i,t) = arma::accu(lv.row(i));
        
        // Treatment: eta * diag_RNA for all stages
        newTreatment(i,t) =
          eta_dt(i,0,t)*oldPop(i,drna_idx[0]) +
          eta_dt(i,1,t)*oldPop(i,drna_idx[1]) +
          eta_dt(i,2,t)*oldPop(i,drna_idx[2]) +
          eta_dt(i,3,t)*oldPop(i,drna_idx[3]) +
          eta_dt(i,4,t)*oldPop(i,drna_idx[4]) +
          eta_dt(i,5,t)*oldPop(i,drna_idx[5]) +
          eta_dt(i,6,t)*oldPop(i,drna_idx[6]) +
          eta_dt(i,7,t)*oldPop(i,drna_idx[7]) +
          eta_dt(i,8,t)*oldPop(i,drna_idx[8]) +
          eta_dt(i,9,t)*oldPop(i,drna_idx[9]);
        
        // Retreat: rho * treat_f for all stages
        newRetreat(i,t) =
          rho_dt(i,0,t)*oldPop(i,trf_idx[0]) +
          rho_dt(i,1,t)*oldPop(i,trf_idx[1]) +
          rho_dt(i,2,t)*oldPop(i,trf_idx[2]) +
          rho_dt(i,3,t)*oldPop(i,trf_idx[3]) +
          rho_dt(i,4,t)*oldPop(i,trf_idx[4]) +
          rho_dt(i,5,t)*oldPop(i,trf_idx[5]) +
          rho_dt(i,6,t)*oldPop(i,trf_idx[6]) +
          rho_dt(i,7,t)*oldPop(i,trf_idx[7]) +
          rho_dt(i,8,t)*oldPop(i,trf_idx[8]) +
          rho_dt(i,9,t)*oldPop(i,trf_idx[9]);
        
        // Cured: cure * treat for all stages
        newCured(i,t) =
          cure_dt(i,0,t)*oldPop(i,tr_idx[0]) +
          cure_dt(i,1,t)*oldPop(i,tr_idx[1]) +
          cure_dt(i,2,t)*oldPop(i,tr_idx[2]) +
          cure_dt(i,3,t)*oldPop(i,tr_idx[3]) +
          cure_dt(i,4,t)*oldPop(i,tr_idx[4]) +
          cure_dt(i,5,t)*oldPop(i,tr_idx[5]) +
          cure_dt(i,6,t)*oldPop(i,tr_idx[6]) +
          cure_dt(i,7,t)*oldPop(i,tr_idx[7]) +
          cure_dt(i,8,t)*oldPop(i,tr_idx[8]) +
          cure_dt(i,9,t)*oldPop(i,tr_idx[9]);
        
        // Treatment failed: lota * treat for all stages
        newtreatfailed(i,t) =
          lota_dt(i,0,t)*oldPop(i,tr_idx[0]) +
          lota_dt(i,1,t)*oldPop(i,tr_idx[1]) +
          lota_dt(i,2,t)*oldPop(i,tr_idx[2]) +
          lota_dt(i,3,t)*oldPop(i,tr_idx[3]) +
          lota_dt(i,4,t)*oldPop(i,tr_idx[4]) +
          lota_dt(i,5,t)*oldPop(i,tr_idx[5]) +
          lota_dt(i,6,t)*oldPop(i,tr_idx[6]) +
          lota_dt(i,7,t)*oldPop(i,tr_idx[7]) +
          lota_dt(i,8,t)*oldPop(i,tr_idx[8]) +
          lota_dt(i,9,t)*oldPop(i,tr_idx[9]);
        
        // Reinfection: a-f3 cured only (0* for f4+)
        newreinfection(i,t) = reinfP(i,t)*foi_II(i)*(
          oldPop(i,6) + oldPop(i,12) + oldPop(i,18)
          + oldPop(i,24) + oldPop(i,30));
        
        // ── Scenario testing outputs ──────────────────────────────────────
        double sp = spc1_dt(i,t);
        
        // newTestingAb_sc: tau_ab_sc * undiag (with (1-sp) for stage a)
        newTestingAb_sc(i,t) =
          tau_ab_sc_dt(i,0,t)*(1.0-sp)*oldPop(i,ud_idx[0]) +
          tau_ab_sc_dt(i,1,t)*oldPop(i,ud_idx[1]) +
          tau_ab_sc_dt(i,2,t)*oldPop(i,ud_idx[2]) +
          tau_ab_sc_dt(i,3,t)*oldPop(i,ud_idx[3]) +
          tau_ab_sc_dt(i,4,t)*oldPop(i,ud_idx[4]) +
          tau_ab_sc_dt(i,5,t)*oldPop(i,ud_idx[5]) +
          tau_ab_sc_dt(i,6,t)*oldPop(i,ud_idx[6]) +
          tau_ab_sc_dt(i,7,t)*oldPop(i,ud_idx[7]) +
          tau_ab_sc_dt(i,8,t)*oldPop(i,ud_idx[8]) +
          tau_ab_sc_dt(i,9,t)*oldPop(i,ud_idx[9]);
        
        // newTestingAg_sc: tau_RNA_sc * diag_ab (from R: uses diag_ab compartments)
        // FIX: was wrongly tau_RNA_sc * tau_ab_sc * undiag
        newTestingAg_sc(i,t) =
          tau_RNA_sc_dt(i,0,t)*oldPop(i,dab_idx[0]) +
          tau_RNA_sc_dt(i,1,t)*oldPop(i,dab_idx[1]) +
          tau_RNA_sc_dt(i,2,t)*oldPop(i,dab_idx[2]) +
          tau_RNA_sc_dt(i,3,t)*oldPop(i,dab_idx[3]) +
          tau_RNA_sc_dt(i,4,t)*oldPop(i,dab_idx[4]) +
          tau_RNA_sc_dt(i,5,t)*oldPop(i,dab_idx[5]) +
          tau_RNA_sc_dt(i,6,t)*oldPop(i,dab_idx[6]) +
          tau_RNA_sc_dt(i,7,t)*oldPop(i,dab_idx[7]) +
          tau_RNA_sc_dt(i,8,t)*oldPop(i,dab_idx[8]) +
          tau_RNA_sc_dt(i,9,t)*oldPop(i,dab_idx[9]);
        
        // newTestingPOCT_sc: tau_poct_sc * undiag (with (1-sp) for stage a)
        newTestingPOCT_sc(i,t) =
          tau_poct_sc_dt(i,0,t)*(1.0-sp)*oldPop(i,ud_idx[0]) +
          tau_poct_sc_dt(i,1,t)*oldPop(i,ud_idx[1]) +
          tau_poct_sc_dt(i,2,t)*oldPop(i,ud_idx[2]) +
          tau_poct_sc_dt(i,3,t)*oldPop(i,ud_idx[3]) +
          tau_poct_sc_dt(i,4,t)*oldPop(i,ud_idx[4]) +
          tau_poct_sc_dt(i,5,t)*oldPop(i,ud_idx[5]) +
          tau_poct_sc_dt(i,6,t)*oldPop(i,ud_idx[6]) +
          tau_poct_sc_dt(i,7,t)*oldPop(i,ud_idx[7]) +
          tau_poct_sc_dt(i,8,t)*oldPop(i,ud_idx[8]) +
          tau_poct_sc_dt(i,9,t)*oldPop(i,ud_idx[9]);
        
        // newTreatment_sc: eta_sc * diag_RNA for all stages
        // FIX: was complex formula; R is simply eta_sc * diag_RNA
        newTreatment_sc(i,t) =
          eta_sc_dt(i,0,t)*oldPop(i,drna_idx[0]) +
          eta_sc_dt(i,1,t)*oldPop(i,drna_idx[1]) +
          eta_sc_dt(i,2,t)*oldPop(i,drna_idx[2]) +
          eta_sc_dt(i,3,t)*oldPop(i,drna_idx[3]) +
          eta_sc_dt(i,4,t)*oldPop(i,drna_idx[4]) +
          eta_sc_dt(i,5,t)*oldPop(i,drna_idx[5]) +
          eta_sc_dt(i,6,t)*oldPop(i,drna_idx[6]) +
          eta_sc_dt(i,7,t)*oldPop(i,drna_idx[7]) +
          eta_sc_dt(i,8,t)*oldPop(i,drna_idx[8]) +
          eta_sc_dt(i,9,t)*oldPop(i,drna_idx[9]);
        
        // ── Negative test outputs ─────────────────────────────────────────
        // Sum of non-a cured compartments (f0_cured..plt_cured)
        double cured_notA = oldPop(i,12)+oldPop(i,18)+oldPop(i,24)+oldPop(i,30)
          +oldPop(i,36)+oldPop(i,42)+oldPop(i,48)+oldPop(i,54)
          +oldPop(i,60);
          double a_cured    = oldPop(i,6);
          double s_pop      = oldPop(i,0);
          
          // newTestingAb_sc_neg:
          // fc * tau_ab_sc[f0] * (non-a_cured + s)  +  fc * tau_ab_sc[a] * a_cured
          // (R uses "f0" progress index for cured_notA and s, "a" index for a_cured)
          newTestingAb_sc_neg(i,t) =
            fc(i,t)*tau_ab_sc_dt(i,1,t)*(cured_notA + s_pop)
            + fc(i,t)*tau_ab_sc_dt(i,0,t)*a_cured;
          
          // newTestingAg_sc_neg: pop-specific multipliers (matches R exactly)
          //   pop 1 (i=0): in R, fc[i,t] is used where `i` is the surrounding loop variable,
          //                so the final stored value reflects fc(npops-1, t).
          //   pop 2 (i=1): 0.5 * rna_f0 * fc(1,t) * cured_notA
          //   pops 3-5 (i=2,3,4): rna_f0 * fc(i,t) * cured_notA
          //                       + rna_a * fc(i,t) * ab_a * a_cured
            {
              double rna_f0 = tau_RNA_sc_dt(i,1,t);  // "f0" progress index = 1
              double rna_a  = tau_RNA_sc_dt(i,0,t);  // "a"  progress index = 0
              double ab_a   = tau_ab_sc_dt(i,0,t);
              if (i == 0) {
                newTestingAg_sc_neg(i,t) =
                  0.5*rna_f0*fc(npops-1, t)*cured_notA;
              } else if (i == 1) {
                newTestingAg_sc_neg(i,t) =
                  0.5*rna_f0*fc(i,t)*cured_notA;
              } else {
                newTestingAg_sc_neg(i,t) =
                  rna_f0*fc(i,t)*cured_notA
                + rna_a*fc(i,t)*ab_a*a_cured;
              }
            }
            
            // newTestingPOCT_sc_neg:
            // fc * tau_poct_sc[f0] * (cured_notA + s)  +  fc * tau_poct_sc[a] * a_cured
            newTestingPOCT_sc_neg(i,t) =
              fc(i,t)*tau_poct_sc_dt(i,1,t)*(cured_notA + s_pop)
              + fc(i,t)*tau_poct_sc_dt(i,0,t)*a_cured;
            
    } // end aggregate loop
  } // end time loop
  
  return List::create(
    Named("allPops")               = allPops,
    Named("newS")                  = newS,
    Named("newEntry")              = newEntry,
    Named("newDeath")              = newDeath,
    Named("newLeave")              = newLeave,
    Named("newInfections")         = newInfections,
    Named("newHCVdeaths")          = newHCVdeaths,
    Named("newTreatment")          = newTreatment,
    Named("newRetreat")            = newRetreat,
    Named("newCured")              = newCured,
    Named("newtreatfailed")        = newtreatfailed,
    Named("newreinfection")        = newreinfection,
    Named("newTestingAb_sc")       = newTestingAb_sc,
    Named("newTestingAg_sc")       = newTestingAg_sc,
    Named("newTestingPOCT_sc")     = newTestingPOCT_sc,
    Named("newTreatment_sc")       = newTreatment_sc,
    Named("newTestingAb_sc_neg")   = newTestingAb_sc_neg,
    Named("newTestingAg_sc_neg")   = newTestingAg_sc_neg,
    Named("newTestingPOCT_sc_neg") = newTestingPOCT_sc_neg
  );
}