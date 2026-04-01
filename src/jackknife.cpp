#include "jackknife.h"
#include "newtraph.h"
#include "frailty_engine.h"

using namespace arma;

//' @title Jacknife estimation (No Frailties)
//' @description Internal C++ function to perform Jackknife leave-one-out cross-validation for models without frailty.
//' @param s Scale parameter.
//' @param n Number of subjects.
//' @param nvar Number of covariates.
//' @param k Integer vector with the number of events per subject.
//' @param nk Total number of observations.
//' @param tau Vector of event/censoring times.
//' @param caltimes Vector of calendar times.
//' @param gaptimes Vector of gap times.
//' @param censored Censoring indicator vector.
//' @param intercepts Vector of intercepts for the effective age.
//' @param slopes Vector of slopes for the effective age.
//' @param lastperrep Vector of last perfect repair times.
//' @param perrepind Indicator of perfect repair.
//' @param effagebegin Effective age at the beginning of the interval.
//' @param effage Effective age at the end of the interval.
//' @param cov Matrix of covariates.
//' @param alphaSeed Initial value for alpha.
//' @param betaSeed Initial vector for betas.
//' @param Z Vector of frailty (fixed to 1 in this case).
//' @param offset Vector of offsets.
//' @param rhoFunc Integer indicating the type of rho function.
//' @param ns Integer vector for subjects (internal use).
//' @param tol Tolerance for convergence.
//' @param maxiter Maximum number of iterations.
//' @return A matrix of size (n x nvar+1) containing the parameter estimates for each jackknife iteration.
//' @export
 // [[Rcpp::export]]
arma::mat Jacknife(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, 
                    arma::vec caltimes, arma::vec gaptimes, arma::vec censored, 
                    arma::vec intercepts, arma::vec slopes, arma::vec lastperrep, 
                    arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
                    arma::mat cov, double alphaSeed, arma::vec betaSeed, 
                    arma::vec Z, arma::vec offset, int rhoFunc, arma::ivec ns, 
                    double tol, int maxiter) {
    
    // Results container: n subjects (rows), nvar+1 params (columns)
    mat estiEndJack(n, nvar + 1, fill::zeros);
    
    // Indices pre-calculation
    uvec stop_idx(n);
    uvec start_idx(n);
    stop_idx(0) = k(0) - 1;
    start_idx(0) = 0;
    for (int i = 1; i < n; ++i) {
        stop_idx(i) = stop_idx(i-1) + k(i);
        start_idx(i) = start_idx(i-1) + k(i-1);
    }
    
    for (int i = 0; i < n; ++i) {
        // "Body" of extract functions: Armadillo's .shed method
        ivec kOK = k; kOK.shed_row(i);
        vec tauOK = tau; tauOK.shed_row(i);
        vec offsetOK = offset; offsetOK.shed_row(i);
        vec censoredOK = censored; censoredOK.shed_row(i);
        
        uvec ids_to_remove = regspace<uvec>(start_idx(i), stop_idx(i));
        
        vec cal_OK = caltimes; cal_OK.shed_rows(ids_to_remove);
        vec gap_OK = gaptimes; gap_OK.shed_rows(ids_to_remove);
        vec int_OK = intercepts; int_OK.shed_rows(ids_to_remove);
        vec slo_OK = slopes; slo_OK.shed_rows(ids_to_remove);
        vec lpr_OK = lastperrep; lpr_OK.shed_rows(ids_to_remove);
        vec pri_OK = perrepind; pri_OK.shed_rows(ids_to_remove);
        vec eab_OK = effagebegin; eab_OK.shed_rows(ids_to_remove);
        vec eff_OK = effage; eff_OK.shed_rows(ids_to_remove);
        mat covOK = cov; covOK.shed_cols(ids_to_remove);
        
        double loglik;
        vec estiEnd(nvar + 1);
        mat info(nvar + 1, nvar + 1);
        int searchBoth, kkiter;
        
        // Re-estimation
        newtraph(s, n - 1, nvar, kOK, nk - k(i), tauOK, cal_OK, gap_OK,
                 censoredOK, int_OK, slo_OK, lpr_OK, pri_OK,
                 eab_OK, eff_OK, covOK, alphaSeed, betaSeed, Z, 
                 offsetOK, rhoFunc, ns, tol, maxiter, loglik, estiEnd, 
                 info, searchBoth, kkiter);
        
        estiEndJack.row(i) = estiEnd.t();
    }
    return estiEndJack;
}



//' @title Jacknife estimation (With Frailties)
//' @description Internal C++ function to perform Jackknife leave-one-out cross-validation for models with frailty.
//' @param s Scale parameter.
//' @param n Number of subjects.
//' @param nvar Number of covariates.
//' @param k Integer vector with the number of events per subject.
//' @param nk Total number of observations.
//' @param tau Vector of event/censoring times.
//' @param caltimes Vector of calendar times.
//' @param gaptimes Vector of gap times.
//' @param censored Censoring indicator vector.
//' @param intercepts Vector of intercepts.
//' @param slopes Vector of slopes.
//' @param lastperrep Vector of last perfect repair times.
//' @param perrepind Indicator of perfect repair.
//' @param effagebegin Effective age start.
//' @param effage Effective age end.
//' @param ndiseff Number of distinct effective ages.
//' @param diseff Vector of distinct effective ages.
//' @param cov Matrix of covariates.
//' @param alphaSeed Initial alpha.
//' @param betaSeed Initial betas.
//' @param xi Initial xi (frailty variance).
//' @param Z Initial frailty vector.
//' @param offset Vector of offsets.
//' @param rhoFunc Type of rho function.
//' @param tol Tolerance.
//' @param maxiter Max iterations.
//' @param maxXi Optimization method for xi.
//' @return A matrix of size (n x nvar+2) containing parameter estimates (alpha, betas, xi) for each iteration.
//' @export
 // [[Rcpp::export]]
 arma::mat Jacknife2(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, 
                     arma::vec caltimes, arma::vec gaptimes, arma::vec censored, 
                     arma::vec intercepts, arma::vec slopes, arma::vec lastperrep, 
                     arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
                     int ndiseff, arma::vec diseff, arma::mat cov, double alphaSeed, 
                     arma::vec betaSeed, double xi, arma::vec Z, arma::vec offset, 
                     int rhoFunc, double tol, int maxiter, int maxXi) {
    
    mat estiEndJack(n, nvar + 2, fill::zeros);
    
    uvec stop_idx(n);
    uvec start_idx(n);
    stop_idx(0) = k(0) - 1;
    start_idx(0) = 0;
    for (int i = 1; i < n; ++i) {
        stop_idx(i) = stop_idx(i-1) + k(i);
        start_idx(i) = start_idx(i-1) + k(i-1);
    }
    
    for (int i = 0; i < n; ++i) {
        ivec kOK = k; kOK.shed_row(i);
        vec tauOK = tau; tauOK.shed_row(i);
        vec offsetOK = offset; offsetOK.shed_row(i);
        vec censoredOK = censored; censoredOK.shed_row(i);
        
        uvec ids_to_remove = regspace<uvec>(start_idx(i), stop_idx(i));
        
        vec cal_OK = caltimes; cal_OK.shed_rows(ids_to_remove);
        vec gap_OK = gaptimes; gap_OK.shed_rows(ids_to_remove);
        vec int_OK = intercepts; int_OK.shed_rows(ids_to_remove);
        vec slo_OK = slopes; slo_OK.shed_rows(ids_to_remove);
        vec lpr_OK = lastperrep; lpr_OK.shed_rows(ids_to_remove);
        vec pri_OK = perrepind; pri_OK.shed_rows(ids_to_remove);
        vec eab_OK = effagebegin; eab_OK.shed_rows(ids_to_remove);
        vec eff_OK = effage; eff_OK.shed_rows(ids_to_remove);
        mat covOK = cov; covOK.shed_cols(ids_to_remove);
        
        // Vector for params + frailties (Z)
        vec estiEnd(nvar + 2 + (n - 1)); 
        // int control_search, control_iter;
        int control_iter;
        arma::ivec control_search(2, arma::fill::zeros);
        double loglikEnd;
        
        // EstimWithFrailty(s, n - 1, nvar, kOK, nk - k(i), tauOK, cal_OK,
        //                  gap_OK, censoredOK, int_OK, slo_OK, 
        //                  lpr_OK, pri_OK, eab_OK, eff_OK, 
        //                  ndiseff, diseff, covOK, alphaSeed, betaSeed, xi, Z, 
        //                  offsetOK, rhoFunc, tol, maxiter, maxXi, 
        //                  estiEnd, control_search, control_iter, loglikEnd);
        
        EstimWithFrailty(s, n - 1, nvar, kOK, nk - k(i), tauOK, cal_OK,
                         gap_OK, censoredOK, int_OK, slo_OK, 
                         lpr_OK, pri_OK, eab_OK, eff_OK, 
                         ndiseff, diseff, covOK, alphaSeed, betaSeed, xi, Z, 
                         offsetOK, rhoFunc, tol, maxiter, maxXi, 
                         estiEnd, control_search, loglikEnd);
        
        // Map back to output row: [alpha, betas..., xi]
        estiEndJack(i, 0) = estiEnd(0);
        for(int j=0; j<nvar; ++j) estiEndJack(i, j+1) = estiEnd(j+1);
        estiEndJack(i, nvar+1) = estiEnd(nvar+1);
    }
    return estiEndJack;
}