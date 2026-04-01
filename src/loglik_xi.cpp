#include "max_xi.h"
#include <cmath>

using namespace arma;

/**
 * @title LogLikXi
 * @description Computes the marginal Log-likelihood (G0) and its first (G1) 
 * and second (G2) derivatives with respect to xi.
 */
void LogLikXi(double xi, int n, const ivec& KK, const vec& A, 
              const vec& B, double& G0, double& G1, double& G2) {
    
    G0 = 0.0;
    G1 = 0.0;
    G2 = 0.0;
    
    for (int i = 0; i < n; ++i) {
        double q1 = 0.0;
        double q2 = 0.0;
        double q3 = 0.0;
        
        // Summation loop for components involving xi + (j-1)
        // In Fortran: do j=1, K(i)
        if (KK[i] > 0) {
            for (int j = 1; j <= KK[i]; ++j) {
                double term = xi + static_cast<double>(j - 1);
                q1 += std::log(term);
                q2 += (1.0 / term);
                q3 += (1.0 / std::pow(term, 2));
            }
        }
        
        // G0: Marginal Log-likelihood
        G0 += q1 + (xi * std::log(xi)) - ((xi + KK[i]) * std::log(xi + A[i])) + B[i];
        
        // G1: First derivative
        G1 += q2 + (std::log(xi) + 1.0) - std::log(xi + A[i]) - 
            ((xi + static_cast<double>(KK[i])) / (xi + A[i]));
        
        // G2: Second derivative
        G2 += -q3 + (1.0 / xi) - (1.0 / (xi + A[i])) - 
            ((A[i] - static_cast<double>(KK[i])) / std::pow(xi + A[i], 2));
    }
}