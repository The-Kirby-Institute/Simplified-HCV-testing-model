#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

#define BASE_VAL(j) (entry1(i,j) + oldPop(i,j) + net_flow(i,j) - death(i,j) - lv(i,j) - death_hcv(i,j))
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

    // ── Compartment index constants ──────────────────────────────────────────
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

    // Progress indices: a=0,f0=1,f1=2,f2=3,f3=4,f4=5,dc=6,hcc=7,lt=8,plt=9
    // Transition indices:
    //   a_f0=0, f0_f1=1, f1_f2=2, f2_f3=3, f3_f4=4,
    //   f3_hcc=5, f4_dc=6, f4_hcc=7, dc_hcc=8, dc_lt=9, hcc_lt=10
    // Fibprog indices:
    //   lt_plt=0, f3c_f4c=1, f3c_hccc=2, f4c_dcc=3, f4c_hccc=4,
    //   dcc_hccc=5, dcc_ltc=6, hccc_ltc=7, ltc_pltc=8

    // Undiag compartments by progress stage
    int ud[10] = {1, 7, 13, 19, 25, 31, 37, 43, 49, 55};
    // diag_ab compartments
    int dab[10] = {2, 8, 14, 20, 26, 32, 38, 44, 50, 56};
    // diag_RNA compartments
    int drna[10] = {3, 9, 15, 21, 27, 33, 39, 45, 51, 57};
    // treat compartments
    int idx_tr[10] = {4, 10, 16, 22, 28, 34, 40, 46, 52, 58};
    // treat_f compartments
    int trf[10] = {5, 11, 17, 23, 29, 35, 41, 47, 53, 59};
    // cured compartments
    int cu[10] = {6, 12, 18, 24, 30, 36, 42, 48, 54, 60};

    // ── Allocate output arrays ───────────────────────────────────────────────
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
    for (int t = 1; t < npts; t++) {

        arma::mat oldPop = allPops.slice(t - 1);
        arma::mat newPop(npops, ncomp, fill::zeros);

        // ── Population sizes ─────────────────────────────────────────────────
        arma::vec N(npops, fill::zeros);
        arma::vec S(npops, fill::zeros);
        arma::vec I(npops, fill::zeros);

        for (int i = 0; i < npops; i++) {
            N(i) = arma::accu(oldPop.row(i));
            S(i) = oldPop(i,0);
            for (int p = 0; p < 10; p++) S(i) += oldPop(i, cu[p]);
            I(i) = N(i) - S(i);
        }

        // ── Prevalence (II) ───────────────────────────────────────────────────
        arma::vec II(npops, fill::zeros);
        if (is_POC_AU) {
            double N12 = N(0) + N(1), I12 = I(0) + I(1);
            double N34 = N(2) + N(3), I34 = I(2) + I(3);
            II(0) = (N12 > 0) ? I12 / N12 : 0.0;
            II(1) = II(0);
            II(2) = (N34 > 0) ? I34 / N34 : 0.0;
            II(3) = II(2);
            II(4) = (N(4) > 0) ? I(4) / N(4) : 0.0;
        } else {
            double Nall = arma::accu(N), Iall = arma::accu(I);
            double IIall = (Nall > 0) ? Iall / Nall : 0.0;
            for (int i = 0; i < npops; i++) II(i) = IIall;
        }

        // ── foi_II ────────────────────────────────────────────────────────────
        arma::vec foi_II(npops);
        for (int i = 0; i < npops; i++)
            foi_II(i) = foi_dt(i, t) * II(i);

        // ── Deaths and leave ──────────────────────────────────────────────────
        arma::mat death(npops, ncomp, fill::zeros);
        arma::mat death_hcv(npops, ncomp, fill::zeros);
        arma::mat lv(npops, ncomp, fill::zeros);

        for (int i = 0; i < npops; i++) {
            for (int j = 0; j < ncomp; j++) {
                death(i, j) = morb_dt(i, t)   * oldPop(i, j);
                lv(i, j)    = leave_dt(i, t)  * oldPop(i, j);
            }
            // HCV deaths (dc stage: 37-42)
            death_hcv(i,37) = mordc_dt(i,t)     * oldPop(i,37);
            death_hcv(i,38) = mordc_dt(i,t)     * oldPop(i,38);
            death_hcv(i,39) = mordc_dt(i,t)     * oldPop(i,39);
            death_hcv(i,40) = mordc_dt(i,t)     * oldPop(i,40);
            death_hcv(i,41) = mordc_dt(i,t)     * oldPop(i,41);
            death_hcv(i,42) = mordcCure_dt(i,t) * oldPop(i,42);
            // hcc: 43-48
            death_hcv(i,43) = morhcc_dt(i,t)    * oldPop(i,43);
            death_hcv(i,44) = morhcc_dt(i,t)    * oldPop(i,44);
            death_hcv(i,45) = morhcc_dt(i,t)    * oldPop(i,45);
            death_hcv(i,46) = morhcc_dt(i,t)    * oldPop(i,46);
            death_hcv(i,47) = morhcc_dt(i,t)    * oldPop(i,47);
            death_hcv(i,48) = morhccCure_dt(i,t)* oldPop(i,48);
            // lt: 49-54
            for (int j = 49; j <= 54; j++)
                death_hcv(i,j) = morlt_dt(i,t) * oldPop(i,j);
            // plt: 55-60
            for (int j = 55; j <= 60; j++)
                death_hcv(i,j) = morplt_dt(i,t) * oldPop(i,j);
        }

        // ── Entry (POC_AU specific) ───────────────────────────────────────────
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

        // ── Net population flow (colSums - rowSums of pop_array * oldPop) ──────
        // net_flow(i,j) = sum_k[pa(k,i)*oldPop(k,j)] - oldPop(i,j)*sum_k[pa(i,k)]
        arma::mat pa_t = pop_array.slice(t - 1);  // npops x npops
        arma::vec pa_rowsum = arma::sum(pa_t, 1); // rowSums: npops vector
        arma::mat net_flow(npops, ncomp, fill::zeros);
        for (int j = 0; j < ncomp; j++) {
            arma::vec col_j = oldPop.col(j);
            arma::vec inflow_j = pa_t.t() * col_j;
            for (int i = 0; i < npops; i++)
                net_flow(i, j) = inflow_j(i) - oldPop(i, j) * pa_rowsum(i);
        }

        // ── Compartment equations ─────────────────────────────────────────────

        // Helper macro to avoid repetition
        

        for (int i = 0; i < npops; i++) {

            // s (0)
            newPop(i,0) = BASE_VAL(0) - foi_II(i)*oldPop(i,0);

            // ── Stage a ──────────────────────────────────────────────────────
            double sp  = spc1_dt(i,t);
            double ta  = tau_ab_dt(i,0,t);
            double tasc= tau_ab_sc_dt(i,0,t);
            double tp  = tau_poct_dt(i,0,t);
            double tpsc= tau_poct_sc_dt(i,0,t);
            double tr_rna   = tau_RNA_dt(i,0,t);
            double tr_rna_sc= tau_RNA_sc_dt(i,0,t);
            double e   = eta_dt(i,0,t);
            double esc = eta_sc_dt(i,0,t);
            double lo  = lota_dt(i,0,t);
            double ro  = rho_dt(i,0,t);
            double cu0 = cure_dt(i,0,t);
            double td  = transition_dt(i,0);  // a_f0

            double sc_diag_rna = tr_rna_sc*tasc*(1.0-sp)*oldPop(i,1)
                               + tpsc*(1.0-sp)*oldPop(i,1);

            newPop(i,1) = BASE_VAL(1)
                - td*oldPop(i,1) - sp*oldPop(i,1)
                - ta*(1.0-sp)*oldPop(i,1) - tasc*(1.0-sp)*oldPop(i,1)
                - tp*(1.0-sp)*oldPop(i,1) - tpsc*(1.0-sp)*oldPop(i,1)
                + foi_II(i)*oldPop(i,0)
                + reinfP(i,t)*foi_II(i)*oldPop(i,6)
                + reinfP(i,t)*foi_II(i)*oldPop(i,12);

            newPop(i,2) = BASE_VAL(2) - td*oldPop(i,2)
                + ta*(1.0-sp)*oldPop(i,1) - tr_rna*oldPop(i,2)
                + tasc*(1.0-sp)*oldPop(i,1) - tr_rna_sc*tasc*(1.0-sp)*oldPop(i,1);

            newPop(i,3) = BASE_VAL(3) - td*oldPop(i,3)
                + tr_rna*oldPop(i,2) + tp*(1.0-sp)*oldPop(i,1)
                - e*oldPop(i,3)
                + tr_rna_sc*tasc*(1.0-sp)*oldPop(i,1) + tpsc*(1.0-sp)*oldPop(i,1)
                - esc*sc_diag_rna;

            newPop(i,4) = BASE_VAL(4) - td*oldPop(i,4)
                + e*oldPop(i,3) + esc*sc_diag_rna
                - lo*oldPop(i,4) - cu0*oldPop(i,4) + ro*oldPop(i,5);

            newPop(i,5) = BASE_VAL(5) - td*oldPop(i,5)
                - ro*oldPop(i,5) + lo*oldPop(i,4);

            newPop(i,6) = BASE_VAL(6)
                + cu0*oldPop(i,4) + sp*oldPop(i,1)
                - reinfP(i,t)*foi_II(i)*oldPop(i,6);

            // ── Stages f0-f3: parametric loop ────────────────────────────────
            // Each fibrosis stage follows same pattern with different indices
            // Base: ud_base = {7,13,19,25}, stage transitions = {1,2,3,4}
            // Extra for f3: f3_hcc transition (idx 5)
            // Cured fib: f0→none, f1→none, f2→none, f3→f4+hcc

            // f0 (stage 1, progress_idx=1, transition in=a_f0(0), out=f0_f1(1))
            {
                int b = 7; // base compartment
                int p = 1; // progress index
                int tin = 0; int tout = 1;
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
                double sc2   = trsc2*tasc2*oldPop(i,b) + tpsc2*oldPop(i,b);

                newPop(i,b)   = BASE_VAL(b)
                    + transition_dt(i,tin)*oldPop(i,ud[0])
                    - transition_dt(i,tout)*oldPop(i,b)
                    - ta2*oldPop(i,b) - tp2*oldPop(i,b)
                    - tasc2*oldPop(i,b) - tpsc2*oldPop(i,b);

                newPop(i,b+1) = BASE_VAL(b+1)
                    + transition_dt(i,tin)*oldPop(i,dab[0])
                    - transition_dt(i,tout)*oldPop(i,b+1)
                    + ta2*oldPop(i,b) - tr2*oldPop(i,b+1)
                    + tasc2*oldPop(i,b) - trsc2*tasc2*oldPop(i,b);

                newPop(i,b+2) = BASE_VAL(b+2)
                    + transition_dt(i,tin)*oldPop(i,drna[0])
                    - transition_dt(i,tout)*oldPop(i,b+2)
                    + tr2*oldPop(i,b+1) + tp2*oldPop(i,b) - e2*oldPop(i,b+2)
                    + trsc2*tasc2*oldPop(i,b) + tpsc2*oldPop(i,b) - esc2*sc2;

                newPop(i,b+3) = BASE_VAL(b+3)
                    + transition_dt(i,tin)*oldPop(i,idx_tr[0])
                    - transition_dt(i,tout)*oldPop(i,b+3)
                    + e2*oldPop(i,b+2) + esc2*sc2
                    - lo2*oldPop(i,b+3) - cu2*oldPop(i,b+3) + ro2*oldPop(i,b+4);

                newPop(i,b+4) = BASE_VAL(b+4)
                    + transition_dt(i,tin)*oldPop(i,trf[0])
                    - transition_dt(i,tout)*oldPop(i,b+4)
                    + lo2*oldPop(i,b+3) - ro2*oldPop(i,b+4);

                newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    - reinfP(i,t)*foi_II(i)*oldPop(i,b+5);
            }

            // f1 (stage 2, progress_idx=2, transition in=f0_f1(1), out=f1_f2(2))
            {
                int b = 13;
                int p = 2;
                int tin = 1; int tout = 2;
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
                double sc2   = trsc2*tasc2*oldPop(i,13) + tpsc2*oldPop(i,13);

                newPop(i,b)   = BASE_VAL(b)
                    + transition_dt(i,tin)*oldPop(i,7)
                    - transition_dt(i,tout)*oldPop(i,b)
                    - ta2*oldPop(i,b) - tp2*oldPop(i,b)
                    - tasc2*oldPop(i,b) - tpsc2*oldPop(i,b)
                    + reinfP(i,t)*foi_II(i)*oldPop(i,18);

                newPop(i,b+1) = BASE_VAL(b+1)
                    + transition_dt(i,tin)*oldPop(i,8)
                    - transition_dt(i,tout)*oldPop(i,b+1)
                    + ta2*oldPop(i,b) - tr2*oldPop(i,b+1)
                    + tasc2*oldPop(i,b) - trsc2*tasc2*oldPop(i,b);

                newPop(i,b+2) = BASE_VAL(b+2)
                    + transition_dt(i,tin)*oldPop(i,9)
                    - transition_dt(i,tout)*oldPop(i,b+2)
                    + tr2*oldPop(i,b+1) + tp2*oldPop(i,b) - e2*oldPop(i,b+2)
                    + trsc2*tasc2*oldPop(i,b) + tpsc2*oldPop(i,b) - esc2*sc2;

                newPop(i,b+3) = BASE_VAL(b+3)
                    + transition_dt(i,tin)*oldPop(i,10)
                    - transition_dt(i,tout)*oldPop(i,b+3)
                    + e2*oldPop(i,b+2) + esc2*sc2
                    - lo2*oldPop(i,b+3) - cu2*oldPop(i,b+3) + ro2*oldPop(i,b+4);

                newPop(i,b+4) = BASE_VAL(b+4)
                    + transition_dt(i,tin)*oldPop(i,11)
                    - transition_dt(i,tout)*oldPop(i,b+4)
                    + lo2*oldPop(i,b+3) - ro2*oldPop(i,b+4);

                newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    - reinfP(i,t)*foi_II(i)*oldPop(i,b+5);
            }

            // f2 (stage 3, progress_idx=3, transition in=f1_f2(2), out=f2_f3(3))
            {
                int b = 19;
                int p = 3;
                int tin = 2; int tout = 3;
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
                double sc2   = trsc2*tasc2*oldPop(i,19) + tpsc2*oldPop(i,19);

                newPop(i,b)   = BASE_VAL(b)
                    + transition_dt(i,tin)*oldPop(i,13)
                    - transition_dt(i,tout)*oldPop(i,b)
                    - ta2*oldPop(i,b) - tp2*oldPop(i,b)
                    - tasc2*oldPop(i,b) - tpsc2*oldPop(i,b)
                    + reinfP(i,t)*foi_II(i)*oldPop(i,24);

                newPop(i,b+1) = BASE_VAL(b+1)
                    + transition_dt(i,tin)*oldPop(i,14)
                    - transition_dt(i,tout)*oldPop(i,b+1)
                    + ta2*oldPop(i,b) - tr2*oldPop(i,b+1)
                    + tasc2*oldPop(i,b) - trsc2*tasc2*oldPop(i,b);

                newPop(i,b+2) = BASE_VAL(b+2)
                    + transition_dt(i,tin)*oldPop(i,15)
                    - transition_dt(i,tout)*oldPop(i,b+2)
                    + tr2*oldPop(i,b+1) + tp2*oldPop(i,b) - e2*oldPop(i,b+2)
                    + trsc2*tasc2*oldPop(i,b) + tpsc2*oldPop(i,b) - esc2*sc2;

                newPop(i,b+3) = BASE_VAL(b+3)
                    + transition_dt(i,tin)*oldPop(i,16)
                    - transition_dt(i,tout)*oldPop(i,b+3)
                    + e2*oldPop(i,b+2) + esc2*sc2
                    - lo2*oldPop(i,b+3) - cu2*oldPop(i,b+3) + ro2*oldPop(i,b+4);

                newPop(i,b+4) = BASE_VAL(b+4)
                    + transition_dt(i,tin)*oldPop(i,17)
                    - transition_dt(i,tout)*oldPop(i,b+4)
                    + lo2*oldPop(i,b+3) - ro2*oldPop(i,b+4);

                newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    - reinfP(i,t)*foi_II(i)*oldPop(i,b+5);
            }

            // f3 (stage 4, progress_idx=4, has f3_hcc extra transition)
            {
                int b = 25;
                int p = 4;
                int tin = 3; int tout = 4; int thcc = 5;
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
                double sc2   = trsc2*tasc2*oldPop(i,25) + tpsc2*oldPop(i,25);

                newPop(i,b)   = BASE_VAL(b)
                    + transition_dt(i,tin)*oldPop(i,19)
                    - transition_dt(i,tout)*oldPop(i,b)
                    - transition_dt(i,thcc)*oldPop(i,b)
                    - ta2*oldPop(i,b) - tp2*oldPop(i,b)
                    - tasc2*oldPop(i,b) - tpsc2*oldPop(i,b)
                    + reinfP(i,t)*foi_II(i)*oldPop(i,30);

                newPop(i,b+1) = BASE_VAL(b+1)
                    + transition_dt(i,tin)*oldPop(i,20)
                    - transition_dt(i,tout)*oldPop(i,b+1)
                    - transition_dt(i,thcc)*oldPop(i,b+1)
                    + ta2*oldPop(i,b) - tr2*oldPop(i,b+1)
                    + tasc2*oldPop(i,b) - trsc2*tasc2*oldPop(i,b);

                newPop(i,b+2) = BASE_VAL(b+2)
                    + transition_dt(i,tin)*oldPop(i,21)
                    - transition_dt(i,tout)*oldPop(i,b+2)
                    - transition_dt(i,thcc)*oldPop(i,b+2)
                    + tr2*oldPop(i,b+1) + tp2*oldPop(i,b) - e2*oldPop(i,b+2)
                    + trsc2*tasc2*oldPop(i,b) + tpsc2*oldPop(i,b) - esc2*sc2;

                newPop(i,b+3) = BASE_VAL(b+3)
                    + transition_dt(i,tin)*oldPop(i,22)
                    - transition_dt(i,tout)*oldPop(i,b+3)
                    - transition_dt(i,thcc)*oldPop(i,b+3)
                    + e2*oldPop(i,b+2) + esc2*sc2
                    - lo2*oldPop(i,b+3) - cu2*oldPop(i,b+3) + ro2*oldPop(i,b+4);

                newPop(i,b+4) = BASE_VAL(b+4)
                    + transition_dt(i,tin)*oldPop(i,23)
                    - transition_dt(i,tout)*oldPop(i,b+4)
                    - transition_dt(i,thcc)*oldPop(i,b+4)
                    + lo2*oldPop(i,b+3) - ro2*oldPop(i,b+4);

                newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    - fibprog_dt(i,1)*oldPop(i,b+5)
                    - fibprog_dt(i,2)*oldPop(i,b+5)
                    - reinfP(i,t)*foi_II(i)*oldPop(i,b+5);
            }

            // f4 (stage 5, progress_idx=5, transitions to dc+hcc)
            {
                int b = 31;
                int p = 5;
                int tin = 4; int tdc = 6; int thcc = 7;
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
                double sc2   = trsc2*tasc2*oldPop(i,31) + tpsc2*oldPop(i,31);

                newPop(i,b)   = BASE_VAL(b)
                    + transition_dt(i,tin)*oldPop(i,25)
                    - transition_dt(i,tdc)*oldPop(i,b)
                    - transition_dt(i,thcc)*oldPop(i,b)
                    - ta2*oldPop(i,b) - tp2*oldPop(i,b)
                    - tasc2*oldPop(i,b) - tpsc2*oldPop(i,b);

                newPop(i,b+1) = BASE_VAL(b+1)
                    + transition_dt(i,tin)*oldPop(i,26)
                    - transition_dt(i,tdc)*oldPop(i,b+1)
                    - transition_dt(i,thcc)*oldPop(i,b+1)
                    + ta2*oldPop(i,b) - tr2*oldPop(i,b+1)
                    + tasc2*oldPop(i,b) - trsc2*tasc2*oldPop(i,b);

                newPop(i,b+2) = BASE_VAL(b+2)
                    + transition_dt(i,tin)*oldPop(i,27)
                    - transition_dt(i,tdc)*oldPop(i,b+2)
                    - transition_dt(i,thcc)*oldPop(i,b+2)
                    + tr2*oldPop(i,b+1) + tp2*oldPop(i,b) - e2*oldPop(i,b+2)
                    + trsc2*tasc2*oldPop(i,b) + tpsc2*oldPop(i,b) - esc2*sc2;

                newPop(i,b+3) = BASE_VAL(b+3)
                    + transition_dt(i,tin)*oldPop(i,28)
                    - transition_dt(i,tdc)*oldPop(i,b+3)
                    - transition_dt(i,thcc)*oldPop(i,b+3)
                    + e2*oldPop(i,b+2) + esc2*sc2
                    - lo2*oldPop(i,b+3) - cu2*oldPop(i,b+3) + ro2*oldPop(i,b+4);

                newPop(i,b+4) = BASE_VAL(b+4)
                    + transition_dt(i,tin)*oldPop(i,29)
                    - transition_dt(i,tdc)*oldPop(i,b+4)
                    - transition_dt(i,thcc)*oldPop(i,b+4)
                    + lo2*oldPop(i,b+3) - ro2*oldPop(i,b+4);

                newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    + fibprog_dt(i,1)*oldPop(i,30)   // f3c→f4c
                    - fibprog_dt(i,3)*oldPop(i,b+5)  // f4c→dcc
                    - fibprog_dt(i,4)*oldPop(i,b+5); // f4c→hccc
            }

            // dc (stage 6, progress_idx=6, transitions to hcc+lt)
            {
                int b = 37;
                int p = 6;
                int tin = 6; int thcc = 8; int tlt = 9;
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
                double sc2   = trsc2*tasc2*oldPop(i,37) + tpsc2*oldPop(i,37);

                newPop(i,b)   = BASE_VAL(b)
                    + transition_dt(i,tin)*oldPop(i,31)
                    - transition_dt(i,thcc)*oldPop(i,b)
                    - transition_dt(i,tlt)*oldPop(i,b)
                    - ta2*oldPop(i,b) - tp2*oldPop(i,b)
                    - tasc2*oldPop(i,b) - tpsc2*oldPop(i,b);

                newPop(i,b+1) = BASE_VAL(b+1)
                    + transition_dt(i,tin)*oldPop(i,32)
                    - transition_dt(i,thcc)*oldPop(i,b+1)
                    - transition_dt(i,tlt)*oldPop(i,b+1)
                    + ta2*oldPop(i,b) - tr2*oldPop(i,b+1)
                    + tasc2*oldPop(i,b) - trsc2*tasc2*oldPop(i,b);

                newPop(i,b+2) = BASE_VAL(b+2)
                    + transition_dt(i,tin)*oldPop(i,33)
                    - transition_dt(i,thcc)*oldPop(i,b+2)
                    - transition_dt(i,tlt)*oldPop(i,b+2)
                    + tr2*oldPop(i,b+1) + tp2*oldPop(i,b) - e2*oldPop(i,b+2)
                    + trsc2*tasc2*oldPop(i,b) + tpsc2*oldPop(i,b) - esc2*sc2;

                newPop(i,b+3) = BASE_VAL(b+3)
                    + transition_dt(i,tin)*oldPop(i,34)
                    - transition_dt(i,thcc)*oldPop(i,b+3)
                    - transition_dt(i,tlt)*oldPop(i,b+3)
                    + e2*oldPop(i,b+2) + esc2*sc2
                    - lo2*oldPop(i,b+3) - cu2*oldPop(i,b+3) + ro2*oldPop(i,b+4);

                newPop(i,b+4) = BASE_VAL(b+4)
                    + transition_dt(i,tin)*oldPop(i,35)
                    - transition_dt(i,thcc)*oldPop(i,b+4)
                    - transition_dt(i,tlt)*oldPop(i,b+4)
                    + lo2*oldPop(i,b+3) - ro2*oldPop(i,b+4);

                newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    + fibprog_dt(i,3)*oldPop(i,36)   // f4c→dcc
                    - fibprog_dt(i,5)*oldPop(i,b+5)  // dcc→hccc
                    - fibprog_dt(i,6)*oldPop(i,b+5); // dcc→ltc
            }

            // hcc (stage 7, progress_idx=7, receives from f4+dc+f3, → lt)
            {
                int b = 43;
                int p = 7;
                int thcc_lt = 10;
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
                double sc2   = trsc2*tasc2*oldPop(i,43) + tpsc2*oldPop(i,43);

                newPop(i,b)   = BASE_VAL(b)
                    + transition_dt(i,7)*oldPop(i,31)   // f4_hcc
                    + transition_dt(i,8)*oldPop(i,37)   // dc_hcc
                    + transition_dt(i,5)*oldPop(i,25)   // f3_hcc
                    - transition_dt(i,thcc_lt)*oldPop(i,b)
                    - ta2*oldPop(i,b) - tp2*oldPop(i,b)
                    - tasc2*oldPop(i,b) - tpsc2*oldPop(i,b);

                newPop(i,b+1) = BASE_VAL(b+1)
                    + transition_dt(i,7)*oldPop(i,32)
                    + transition_dt(i,8)*oldPop(i,38)
                    + transition_dt(i,5)*oldPop(i,26)
                    - transition_dt(i,thcc_lt)*oldPop(i,b+1)
                    + ta2*oldPop(i,b) - tr2*oldPop(i,b+1)
                    + tasc2*oldPop(i,b) - trsc2*tasc2*oldPop(i,b);

                newPop(i,b+2) = BASE_VAL(b+2)
                    + transition_dt(i,7)*oldPop(i,33)
                    + transition_dt(i,8)*oldPop(i,39)
                    + transition_dt(i,5)*oldPop(i,27)
                    - transition_dt(i,thcc_lt)*oldPop(i,b+2)
                    + tr2*oldPop(i,b+1) + tp2*oldPop(i,b) - e2*oldPop(i,b+2)
                    + trsc2*tasc2*oldPop(i,b) + tpsc2*oldPop(i,b) - esc2*sc2;

                newPop(i,b+3) = BASE_VAL(b+3)
                    + transition_dt(i,7)*oldPop(i,34)
                    + transition_dt(i,8)*oldPop(i,40)
                    + transition_dt(i,5)*oldPop(i,28)
                    - transition_dt(i,thcc_lt)*oldPop(i,b+3)
                    + e2*oldPop(i,b+2) + esc2*sc2
                    - lo2*oldPop(i,b+3) - cu2*oldPop(i,b+3) + ro2*oldPop(i,b+4);

                newPop(i,b+4) = BASE_VAL(b+4)
                    + transition_dt(i,7)*oldPop(i,35)
                    + transition_dt(i,5)*oldPop(i,29)
                    + transition_dt(i,8)*oldPop(i,41)
                    - transition_dt(i,thcc_lt)*oldPop(i,b+4)
                    + lo2*oldPop(i,b+3) - ro2*oldPop(i,b+4);

                newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    + fibprog_dt(i,4)*oldPop(i,36)   // f4c→hccc
                    + fibprog_dt(i,5)*oldPop(i,42)   // dcc→hccc
                    + fibprog_dt(i,2)*oldPop(i,30)   // f3c→hccc
                    - fibprog_dt(i,7)*oldPop(i,b+5); // hccc→ltc
            }

            // lt (stage 8, progress_idx=8, receives from hcc+dc, → plt via fibprog)
            {
                int b = 49;
                int p = 8;
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
                double sc2   = trsc2*tasc2*oldPop(i,49) + tpsc2*oldPop(i,49);
                double fp_lt = fibprog_dt(i,0); // lt_plt

                newPop(i,b)   = BASE_VAL(b)
                    + transition_dt(i,10)*oldPop(i,43)  // hcc_lt
                    + transition_dt(i,9)*oldPop(i,37)   // dc_lt
                    - fp_lt*oldPop(i,b)
                    - ta2*oldPop(i,b) - tp2*oldPop(i,b)
                    - tasc2*oldPop(i,b) - tpsc2*oldPop(i,b);

                newPop(i,b+1) = BASE_VAL(b+1)
                    + transition_dt(i,10)*oldPop(i,44)
                    + transition_dt(i,9)*oldPop(i,38)
                    - fp_lt*oldPop(i,b+1)
                    + ta2*oldPop(i,b) - tr2*oldPop(i,b+1)
                    + tasc2*oldPop(i,b) - trsc2*tasc2*oldPop(i,b);

                newPop(i,b+2) = BASE_VAL(b+2)
                    + transition_dt(i,10)*oldPop(i,45)
                    + transition_dt(i,9)*oldPop(i,39)
                    - fp_lt*oldPop(i,b+2)
                    + tr2*oldPop(i,b+1) + tp2*oldPop(i,b) - e2*oldPop(i,b+2)
                    + trsc2*tasc2*oldPop(i,b) + tpsc2*oldPop(i,b) - esc2*sc2;

                newPop(i,b+3) = BASE_VAL(b+3)
                    + transition_dt(i,10)*oldPop(i,46)
                    + transition_dt(i,9)*oldPop(i,40)
                    + e2*oldPop(i,b+2) + esc2*sc2
                    - fp_lt*oldPop(i,b+3)
                    - lo2*oldPop(i,b+3) - cu2*oldPop(i,b+3) + ro2*oldPop(i,b+4);

                newPop(i,b+4) = BASE_VAL(b+4)
                    + transition_dt(i,10)*oldPop(i,47)
                    + transition_dt(i,9)*oldPop(i,41)
                    + lo2*oldPop(i,b+3) - fp_lt*oldPop(i,b+4)
                    - ro2*oldPop(i,b+4);

                newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    + fibprog_dt(i,7)*oldPop(i,48)   // hccc→ltc
                    + fibprog_dt(i,6)*oldPop(i,42)   // dcc→ltc
                    - fibprog_dt(i,8)*oldPop(i,b+5); // ltc→pltc
            }

            // plt (stage 9, progress_idx=9, receives from lt via fibprog)
            {
                int b = 55;
                int p = 9;
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
                double sc2   = trsc2*tasc2*oldPop(i,55) + tpsc2*oldPop(i,55);
                double fp_lt = fibprog_dt(i,0);

                newPop(i,b)   = BASE_VAL(b)
                    - ta2*oldPop(i,b) - tp2*oldPop(i,b)
                    - tasc2*oldPop(i,b) - tpsc2*oldPop(i,b)
                    + fp_lt*oldPop(i,49);

                newPop(i,b+1) = BASE_VAL(b+1)
                    + fp_lt*oldPop(i,50)
                    + ta2*oldPop(i,b) - tr2*oldPop(i,b+1)
                    + tasc2*oldPop(i,b) - trsc2*tasc2*oldPop(i,b);

                newPop(i,b+2) = BASE_VAL(b+2)
                    + fp_lt*oldPop(i,51)
                    + tr2*oldPop(i,b+1) + tp2*oldPop(i,b) - e2*oldPop(i,b+2)
                    + trsc2*tasc2*oldPop(i,b) + tpsc2*oldPop(i,b) - esc2*sc2;

                newPop(i,b+3) = BASE_VAL(b+3)
                    + fp_lt*oldPop(i,52)
                    + e2*oldPop(i,b+2) + esc2*sc2
                    - lo2*oldPop(i,b+3) - cu2*oldPop(i,b+3) + ro2*oldPop(i,b+4);

                newPop(i,b+4) = BASE_VAL(b+4)
                    + fp_lt*oldPop(i,53)
                    + lo2*oldPop(i,b+3) - ro2*oldPop(i,b+4);

                newPop(i,b+5) = BASE_VAL(b+5)
                    + cu2*oldPop(i,b+3)
                    + fibprog_dt(i,8)*oldPop(i,54);  // ltc→pltc
            }
        } // end for each pop i

        

        // ── Floor at zero ────────────────────────────────────────────────────
        // clamp removed to match R behavior: newPop = arma::clamp(newPop, 0.0, arma::datum::inf);
        allPops.slice(t) = newPop;

        // ── Result aggregates ─────────────────────────────────────────────────
        for (int i = 0; i < npops; i++) {
            newS(i, t) = newPop(i, 0);

            // Infections
            newInfections(i,t) = foi_II(i)*oldPop(i,0)
                + reinfP(i,t)*foi_II(i)*(
                    oldPop(i,6)+oldPop(i,12)+oldPop(i,18)+
                    oldPop(i,24)+oldPop(i,30));

            // HCV deaths
            newHCVdeaths(i,t) =
                mordc_dt(i,t) *(oldPop(i,37)+oldPop(i,38)+oldPop(i,39)+oldPop(i,40)+oldPop(i,41)) +
                mordcCure_dt(i,t)*oldPop(i,42) +
                morhcc_dt(i,t)*(oldPop(i,43)+oldPop(i,44)+oldPop(i,45)+oldPop(i,46)+oldPop(i,47)) +
                morhccCure_dt(i,t)*oldPop(i,48) +
                morlt_dt(i,t)*(oldPop(i,49)+oldPop(i,50)+oldPop(i,51)+
                               oldPop(i,52)+oldPop(i,53)+oldPop(i,54)) +
                morplt_dt(i,t)*(oldPop(i,55)+oldPop(i,56)+oldPop(i,57)+
                                oldPop(i,58)+oldPop(i,59)+oldPop(i,60));

            newEntry(i,t) = arma::accu(entry1.row(i));
            newDeath(i,t) = arma::accu(death.row(i));
            newLeave(i,t) = arma::accu(lv.row(i));

            // Treatment/retreat/cured/failed
            newTreatment(i,t)   = eta_dt(i,0,t)*oldPop(i,3)  + eta_dt(i,1,t)*oldPop(i,9) +
                                   eta_dt(i,2,t)*oldPop(i,15) + eta_dt(i,3,t)*oldPop(i,21)+
                                   eta_dt(i,4,t)*oldPop(i,27) + eta_dt(i,5,t)*oldPop(i,33)+
                                   eta_dt(i,6,t)*oldPop(i,39) + eta_dt(i,7,t)*oldPop(i,45)+
                                   eta_dt(i,8,t)*oldPop(i,51) + eta_dt(i,9,t)*oldPop(i,57);

            newRetreat(i,t)     = rho_dt(i,0,t)*oldPop(i,5)  + rho_dt(i,1,t)*oldPop(i,11)+
                                   rho_dt(i,2,t)*oldPop(i,17) + rho_dt(i,3,t)*oldPop(i,23)+
                                   rho_dt(i,4,t)*oldPop(i,29) + rho_dt(i,5,t)*oldPop(i,35)+
                                   rho_dt(i,6,t)*oldPop(i,41) + rho_dt(i,7,t)*oldPop(i,47)+
                                   rho_dt(i,8,t)*oldPop(i,53) + rho_dt(i,9,t)*oldPop(i,59);

            newCured(i,t)       = cure_dt(i,0,t)*oldPop(i,4)  + cure_dt(i,1,t)*oldPop(i,10)+
                                   cure_dt(i,2,t)*oldPop(i,16) + cure_dt(i,3,t)*oldPop(i,22)+
                                   cure_dt(i,4,t)*oldPop(i,28) + cure_dt(i,5,t)*oldPop(i,34)+
                                   cure_dt(i,6,t)*oldPop(i,40) + cure_dt(i,7,t)*oldPop(i,46)+
                                   cure_dt(i,8,t)*oldPop(i,52) + cure_dt(i,9,t)*oldPop(i,58);

            newtreatfailed(i,t) = lota_dt(i,0,t)*oldPop(i,4)  + lota_dt(i,1,t)*oldPop(i,10)+
                                   lota_dt(i,2,t)*oldPop(i,16) + lota_dt(i,3,t)*oldPop(i,22)+
                                   lota_dt(i,4,t)*oldPop(i,28) + lota_dt(i,5,t)*oldPop(i,34)+
                                   lota_dt(i,6,t)*oldPop(i,40) + lota_dt(i,7,t)*oldPop(i,46)+
                                   lota_dt(i,8,t)*oldPop(i,52) + lota_dt(i,9,t)*oldPop(i,58);

            newreinfection(i,t) = reinfP(i,t)*foi_II(i)*(
                oldPop(i,6)+oldPop(i,12)+oldPop(i,18)+oldPop(i,24)+oldPop(i,30));

            // Scenario testing
            double sp = spc1_dt(i,t);
            double cured_sum = oldPop(i,12)+oldPop(i,18)+oldPop(i,24)+oldPop(i,30)+
                               oldPop(i,36)+oldPop(i,42)+oldPop(i,48)+oldPop(i,54)+oldPop(i,60);

            newTestingAb_sc(i,t) =
                tau_ab_sc_dt(i,0,t)*(1.0-sp)*oldPop(i,1)  +
                tau_ab_sc_dt(i,1,t)*oldPop(i,7)  + tau_ab_sc_dt(i,2,t)*oldPop(i,13)+
                tau_ab_sc_dt(i,3,t)*oldPop(i,19) + tau_ab_sc_dt(i,4,t)*oldPop(i,25)+
                tau_ab_sc_dt(i,5,t)*oldPop(i,31) + tau_ab_sc_dt(i,6,t)*oldPop(i,37)+
                tau_ab_sc_dt(i,7,t)*oldPop(i,43) + tau_ab_sc_dt(i,8,t)*oldPop(i,49)+
                tau_ab_sc_dt(i,9,t)*oldPop(i,55);

            newTestingAb_sc_neg(i,t) =
                fc(i,t)*tau_ab_sc_dt(i,1,t)*(cured_sum + oldPop(i,0)) +
                fc(i,t)*tau_ab_sc_dt(i,0,t)*oldPop(i,6);

            newTestingAg_sc(i,t) =
                tau_RNA_sc_dt(i,0,t)*tau_ab_sc_dt(i,0,t)*(1.0-sp)*oldPop(i,1) +
                tau_RNA_sc_dt(i,1,t)*tau_ab_sc_dt(i,1,t)*oldPop(i,7)  +
                tau_RNA_sc_dt(i,2,t)*tau_ab_sc_dt(i,2,t)*oldPop(i,13) +
                tau_RNA_sc_dt(i,3,t)*tau_ab_sc_dt(i,3,t)*oldPop(i,19) +
                tau_RNA_sc_dt(i,4,t)*tau_ab_sc_dt(i,4,t)*oldPop(i,25) +
                tau_RNA_sc_dt(i,5,t)*tau_ab_sc_dt(i,5,t)*oldPop(i,31) +
                tau_RNA_sc_dt(i,6,t)*tau_ab_sc_dt(i,6,t)*oldPop(i,37) +
                tau_RNA_sc_dt(i,7,t)*tau_ab_sc_dt(i,7,t)*oldPop(i,43) +
                tau_RNA_sc_dt(i,8,t)*tau_ab_sc_dt(i,8,t)*oldPop(i,49) +
                tau_RNA_sc_dt(i,9,t)*tau_ab_sc_dt(i,9,t)*oldPop(i,55);

            newTestingAg_sc_neg(i,t) =
                tau_RNA_sc_dt(i,1,t)*fc(i,t)*tau_ab_sc_dt(i,1,t)*cured_sum +
                tau_RNA_sc_dt(i,0,t)*fc(i,t)*tau_ab_sc_dt(i,0,t)*oldPop(i,6);

            newTestingPOCT_sc(i,t) =
                tau_poct_sc_dt(i,0,t)*(1.0-sp)*oldPop(i,1) +
                tau_poct_sc_dt(i,1,t)*oldPop(i,7)  + tau_poct_sc_dt(i,2,t)*oldPop(i,13)+
                tau_poct_sc_dt(i,3,t)*oldPop(i,19) + tau_poct_sc_dt(i,4,t)*oldPop(i,25)+
                tau_poct_sc_dt(i,5,t)*oldPop(i,31) + tau_poct_sc_dt(i,6,t)*oldPop(i,37)+
                tau_poct_sc_dt(i,7,t)*oldPop(i,43) + tau_poct_sc_dt(i,8,t)*oldPop(i,49)+
                tau_poct_sc_dt(i,9,t)*oldPop(i,55);

            newTestingPOCT_sc_neg(i,t) =
                fc(i,t)*tau_poct_sc_dt(i,1,t)*(cured_sum + oldPop(i,0)) +
                fc(i,t)*tau_poct_sc_dt(i,0,t)*oldPop(i,6);

            double sc_sum = 0.0;

            // Stage a: (1-sp) factor
            sc_sum += eta_sc_dt(i,0,t)*oldPop(i,ud[0])*(
                tau_RNA_sc_dt(i,0,t)*tau_ab_sc_dt(i,0,t)*(1.0-sp) +
                tau_poct_sc_dt(i,0,t)*(1.0-sp));
            // Stages f0-f3 (p=1-4): standard
            for (int p = 1; p <= 4; p++) {
                sc_sum += eta_sc_dt(i,p,t)*oldPop(i,ud[p])*(
                    tau_RNA_sc_dt(i,p,t)*tau_ab_sc_dt(i,p,t) +
                    tau_poct_sc_dt(i,p,t));
            }
            // Stages f4 (p=5): standard
            sc_sum += eta_sc_dt(i,5,t)*oldPop(i,ud[5])*(
                tau_RNA_sc_dt(i,5,t)*tau_ab_sc_dt(i,5,t) +
                tau_poct_sc_dt(i,5,t));
            // Stage dc (p=6): replicate R bug — tau_RNA_sc uses [,] = sum over all pops
            {
                double tau_rna_dc_sum = 0.0;
                for (int k = 0; k < npops; k++)
                    tau_rna_dc_sum += tau_RNA_sc_dt(k,6,t);
                sc_sum += eta_sc_dt(i,6,t)*oldPop(i,ud[6])*(
                    tau_rna_dc_sum * tau_ab_sc_dt(i,6,t) +
                    (double)npops * tau_poct_sc_dt(i,6,t));
            }
            // Stages hcc, lt, plt (p=7-9): standard
            for (int p = 7; p < 10; p++) {
                sc_sum += eta_sc_dt(i,p,t)*oldPop(i,ud[p])*(
                    tau_RNA_sc_dt(i,p,t)*tau_ab_sc_dt(i,p,t) +
                    tau_poct_sc_dt(i,p,t));
            }
            newTreatment_sc(i,t) = sc_sum;
        }
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
