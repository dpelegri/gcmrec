#include "max_xi.h"
#include <cmath>

using namespace arma;

/**
 * @title MaxWrtXi
 * @description Maximizes the log-likelihood for the frailty parameter Xi.
 */
void MaxWrtXi(double xiOld, int n, const ivec& K, const vec& A, const vec& B,
              double tol, double& xiNew, double& G0, int maxiter, int& search) {
    
    const double BIG_VALUE = 1.0e60;
    ivec KK = K - 1; // Adjustment as per original logic
    
    int niter = 0;
    search = 0;
    double dist = 100.0;
    
    double etaOld = std::log(xiOld);
    double xiratOld = xiOld / (1.0 + xiOld);
    
    double G1, G2;
    double xiratNew, etaNew;
    
    while ((dist > tol) && (niter < maxiter)) {
        niter++;
        
        // Compute log-likelihood and its derivatives
        LogLikXi(xiOld, n, KK, A, B, G0, G1, G2);
        
        // Newton-Raphson update on the log scale (eta)
        etaNew = etaOld - (G1 / (G1 + xiOld * G2));
        xiNew = std::exp(etaNew);
        
        if (xiNew > BIG_VALUE) {
            xiNew = BIG_VALUE;
            search = 1;
            return; 
        } else {
            xiratNew = xiNew / (1.0 + xiNew);
            dist = std::abs(xiratOld - xiratNew);
            
            // Updates for next iteration
            xiOld = xiNew;
            etaOld = etaNew;
            xiratOld = xiratNew;
        }
    }
    
    if (dist <= tol) {
        search = 1;
    }
}