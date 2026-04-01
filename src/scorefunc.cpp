#include "scorefunc.h"
#include "at_risk_utils.h"
#include "math_utils.h"

using namespace arma;


//' @title Likelihood and Score Calculation
//' @description Computes the log-likelihood, score vector, and information matrix 
//' for a given set of parameters in a recurrent event model.
//' @param s Current time point.
//' @param n Number of subjects.
//' @param nvar Number of covariates.
//' @param k Observations per subject.
//' @param nk Total number of observations.
//' @param tau Censoring times.
//' @param caltimes Calendar times.
//' @param gaptimes Gap times.
//' @param censored Censoring indicators.
//' @param intercepts Model intercepts.
//' @param slopes Model slopes.
//' @param lastperrep Periodic report indices.
//' @param perrepind Periodic report indicators.
//' @param effagebegin Effective age starts.
//' @param effage Effective age ends.
//' @param cov Covariate matrix.
//' @param alpha Recurrence parameter.
//' @param beta Regression coefficients.
//' @param Z Fragility/Random effect vector.
//' @param offset Offsets vector.
//' @param rhoFunc Type of rho function.
//' @param ns Vector of event counts per subject.
//' @param loglik [Out] Computed log-likelihood.
//' @param score [Out] Full score vector (alpha + beta).
//' @param info [Out] Full information matrix.
//' @param scorealpha [Out] Score component for alpha.
//' @param scorebeta [Out] Score component for beta.
//' @param infoalpha [Out] Information component for alpha.
//' @param infobeta [Out] Information matrix for beta.
//' @export
 // [[Rcpp::export]]
void scorefunc(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, arma::vec caltimes, 
            arma::vec gaptimes, arma::vec censored, arma::vec intercepts, arma::vec slopes, 
            arma::vec lastperrep, arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
            arma::mat cov, double alpha, arma::vec beta, arma::vec Z, arma::vec offset, int rhoFunc, 
            arma::ivec ns, double& loglik, arma::vec& score, arma::mat& info, 
            double& scorealpha, arma::vec& scorebeta, double& infoalpha, arma::mat& infobeta) 
{
    
    // Global Initializations
    loglik = 0.0;
    scorealpha = 0.0;
    infoalpha = 0.0;
    scorebeta.zeros(nvar);
    infobeta.zeros(nvar, nvar);
    
    vec infoalphabeta(nvar, fill::zeros);
    
    // Auxiliary variables for AtRisk results
    double S0, GrS0A1, Gr2S0A1;
    vec GrS0Be(nvar), Gr2S0A1Be(nvar), Ysubj(n);
    mat Gr2S0Be(nvar, nvar);
    
    // Get number of events per subject
    nsm(s, n, nk, caltimes, k, ns);
    
    int pos = 0; 
    for (int i = 0; i < n; ++i) {
        int ki = k[i];
        
        // --- DYNAMIC FIX: Removing the 200 limit ---
        // We create a view of the effective age for the current subject
        // without copying data, just referencing the memory from 'pos' to 'pos + ki - 1'
        vec effageOK = effage.subvec(pos, pos + ki - 1);
        
        // We pick the covariate vector for the current subject 
        // (Assuming first observation's covariates as per original Fortran logic)
        vec subjectCovariate = cov.col(pos);
        
        // Event-level loop
        if (ns[i] > 1) {
            for (int j = 1; j < ns[i]; ++j) { // j=1 is the 2nd event (0-based)
                double w = effageOK[j];
                
                // Call AtRisk (which is already dynamic)
                AtRisk(s, ns, w, n, nvar, nk, k, tau, caltimes, gaptimes, censored, 
                       intercepts, slopes, lastperrep, perrepind, effagebegin, 
                       effage, cov, alpha, beta, Z, offset, rhoFunc, 
                       Ysubj, S0, GrS0A1, GrS0Be, Gr2S0A1, Gr2S0A1Be, Gr2S0Be);
                
                // Update Log-Likelihood
                loglik += std::log(rho(j - 1, alpha, rhoFunc)) + 
                    std::log(psi(nvar, subjectCovariate, beta, offset[i])) - 
                    std::log(S0);
                
                // Update Alpha components
                double ealpha = GrS0A1 / S0;
                scorealpha += (double(j - 1) / alpha) - ealpha;
                infoalpha += (double(j - 1) / (alpha * alpha)) + (Gr2S0A1 / S0) - std::pow(ealpha, 2);
                
                // Update Beta and Cross components
                for (int t = 0; t < nvar; ++t) {
                    double ebeta_t = GrS0Be[t] / S0;
                    scorebeta[t] += subjectCovariate[t] - ebeta_t;
                    infoalphabeta[t] += (Gr2S0A1Be[t] / S0) - (ealpha * ebeta_t);
                    
                    for (int r = 0; r < nvar; ++r) {
                        infobeta(t, r) += (Gr2S0Be(t, r) / S0) - (ebeta_t * (GrS0Be[r] / S0));
                    }
                }
            }
        }
        pos += ki; // Advance pointer to the next subject's data
    }
    
    // Final matrix assembly (Score vector and Information matrix)
    score[0] = scorealpha;
    info(0, 0) = infoalpha;
    for (int i = 0; i < nvar; ++i) {
        score[i + 1] = scorebeta[i];
        info(0, i + 1) = infoalphabeta[i];
        info(i + 1, 0) = infoalphabeta[i];
        for (int j = 0; j < nvar; ++j) {
            info(i + 1, j + 1) = infobeta(i, j);
        }
    }
}