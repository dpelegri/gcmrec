#include "est_lamb_surv.h" 
#include "nsm.h"           
#include "at_risk_utils.h" 


// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

//' @title Non-parametric Hazard and Survival Estimation
//' @description Computes the baseline cumulative hazard (Nelson-Aalen) and 
//' survival (Kaplan-Meier) functions using the at-risk process.
//' @param s Current calendar time.
//' @param n Number of subjects.
//' @param nvar Number of covariates.
//' @param k Vector with number of observations per subject.
//' @param nk Total number of observations.
//' @param tau Truncation times.
//' @param caltimes Calendar times of events.
//' @param gaptimes Gap times between events.
//' @param censored Censoring indicators.
//' @param intercepts Covariate intercepts.
//' @param slopes Covariate slopes.
//' @param lastperrep Last periodic report time.
//' @param perrepind Periodic report indicator.
//' @param effagebegin Effective age at start.
//' @param effage Effective age at event.
//' @param ndiseff Number of unique event times.
//' @param diseff Sorted unique event times.
//' @param cov Covariate matrix.
//' @param alpha Alpha parameter (intensity scaling).
//' @param beta Beta coefficients vector.
//' @param Z Vector of frailties.
//' @param offset Offset vector.
//' @param rhoFunc Model type indicator.
//' @return A List containing lambdafunc, deltalambdafunc, and survfunc.
//' @export
// [[Rcpp::export]]
Rcpp::List EstLambSurv(double s, int n, int nvar, const arma::ivec& k, int nk, 
                       const arma::vec& tau, const arma::vec& caltimes, 
                       const arma::vec& gaptimes, const arma::vec& censored, 
                       const arma::vec& intercepts, const arma::vec& slopes, 
                       const arma::vec& lastperrep, const arma::vec& perrepind, 
                       const arma::vec& effagebegin, const arma::vec& effage, 
                       int ndiseff, const arma::vec& diseff, const arma::mat& cov, 
                       double alpha, const arma::vec& beta, const arma::vec& Z, 
                       const arma::vec& offset, int rhoFunc,
                       arma::vec& lambda, arma::vec& deltalambda, arma::vec& surv) {
    
    // Initialize vectors that were previously passed as arguments in Fortran
    vec lambdafunc(ndiseff, fill::zeros);
    vec deltalambdafunc(ndiseff, fill::zeros);
    vec survfunc(ndiseff, fill::ones);
    
    vec DeltaNst(ndiseff, fill::zeros);
    vec S0st(ndiseff, fill::zeros);
    
    // Auxiliary variables for AtRisk
    vec Ysubj(n);
    double GrS0A1, Gr2S0A1, S0_val;
    vec GrS0Be(nvar), Gr2S0A1Be(nvar);
    mat Gr2S0Be(nvar, nvar);
    
    ivec ns(n);
    nsm(s, n, nk, caltimes, k, ns);
    
    for (int i = 1; i < ndiseff; ++i) {
        for (int j = 0; j < nk; ++j) {
            if (caltimes[j] <= s && effage[j] == diseff[i]) {
                DeltaNst[i] += 1.0;
            }
        }
        
        AtRisk(s, ns, diseff[i], n, nvar, nk, k, tau, caltimes,
               gaptimes, censored, intercepts, slopes, lastperrep, perrepind,
               effagebegin, effage, cov, alpha, beta, Z, offset, rhoFunc,
               Ysubj, S0_val, GrS0A1, GrS0Be, Gr2S0A1, Gr2S0A1Be, Gr2S0Be);
        
        S0st[i] = S0_val;
        
        if (S0st[i] > 0.0) {
            double dellamb = std::min(DeltaNst[i] / S0st[i], 1.0);
            lambdafunc[i] = lambdafunc[i - 1] + dellamb;
            survfunc[i] = survfunc[i - 1] * (1.0 - dellamb);
        } else {
            lambdafunc[i] = lambdafunc[i - 1];
            survfunc[i] = survfunc[i - 1];
        }
    }
    
    for (int i = 1; i < ndiseff; ++i) {
        deltalambdafunc[i] = lambdafunc[i] - lambdafunc[i - 1];
    }
    
    // Return as a named list for R
    return Rcpp::List::create(
        Rcpp::_["Lamb"]      = lambdafunc,
        Rcpp::_["DeltaLamb"] = deltalambdafunc,
        Rcpp::_["Surv"]      = survfunc
    );
}