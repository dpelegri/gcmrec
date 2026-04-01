#include "comp_ak.h"
#include "nsm.h"
#include "at_risk_utils.h"
#include "math_utils.h" // Assuming rho and psi are declared here

using namespace arma;

/**
 * @title CompAK (Fiel al original Fortran)
 * @description Computes Ki's and Ai's to obtain new Z-values (frailty).
 */
void CompAK(double s, int n, int nvar, const ivec& k, int nk, const vec& tau,
            const vec& caltimes, const vec& gaptimes, const vec& censored,
            const vec& intercepts, const vec& slopes, const vec& lastperrep,
            const vec& perrepind, const vec& effagebegin, const vec& effage,
            int ndiseff, const vec& diseff, const mat& cov, double alpha, 
            const vec& beta, const vec& deltalambdafunc, const vec& Z, 
            const vec& offset, int rhoFunc, ivec& KK, vec& AA, vec& BB) {
    
    // Internal variables
    double S0_val, GrS0A1, Gr2S0A1;
    vec GrS0Be(nvar), Gr2S0A1Be(nvar);
    mat Gr2S0Be(nvar, nvar);
    vec Ysubj(n);
    
    // Step 1: Initialize AA and compute KK using nsm
    AA.zeros();
    nsm(s, n, nk, caltimes, k, KK);
    
    // Step 2: Cumulative hazard component (AA)
    for (int i = 0; i < ndiseff; ++i) {
        // We reuse AtRisk to get Ysubj for the current baseline hazard time
        AtRisk(s, KK, diseff[i], n, nvar, nk, k, tau, caltimes,
               gaptimes, censored, intercepts, slopes, lastperrep,
               perrepind, effagebegin, effage, cov, alpha, beta, Z, offset,
               rhoFunc, Ysubj, S0_val, GrS0A1, GrS0Be, Gr2S0A1, Gr2S0A1Be, Gr2S0Be);
        
        // AA(j) = AA(j) + (Ysubj(j) * deltalambdafunc(i))
        // AA += (Ysubj % deltalambdafunc[i]); 
        AA += (Ysubj * deltalambdafunc[i]);
    }
    
    // Step 3: Log-likelihood component (BB)
    BB.zeros();
    int pos = 0; // 0-based index for caltimes/effage/cov
    
    for (int i = 0; i < n; ++i) {
        int ki = k[i];
        vec covariate(nvar);
        vec effageOK(ki); 
        
        // Extract subject's data
        for (int r = 0; r < ki; ++r) {
            covariate = cov.col(pos);
            effageOK[r] = effage[pos];
            pos++;
        }
        
        // Logic for BB if KK(i) > 1
        if (KK[i] > 1) {
            int j = 1; // Corresponds to Fortran j=2 (0-based j=1)
            while (j < KK[i]) {
                // BB(i) = BB(i) + log(rho(j-2)) + log(psi(covariate, beta, offset))
                // Note: Using j-1 because rho seems to expect the event index
                BB[i] += std::log(rho(j - 1, alpha, rhoFunc)) + 
                    std::log(psi(nvar, covariate, beta, offset[i]));
                
                int u = 0;
                int saca = 0;
                while (u < ndiseff && saca == 0) {
                    if (effageOK[j] == diseff[u]) {
                        BB[i] += std::log(deltalambdafunc[u]);
                        saca = 1;
                    }
                    u++;
                }
                j++;
            }
        }
    }
}