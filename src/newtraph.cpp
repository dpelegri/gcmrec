#include "newtraph.h"
#include "newtraph_utils.h"

// #include <RcppArmadillo.h>
// // Incluimos las cabeceras de las funciones que se llaman internamente
// // Asegúrate de que estas funciones estén declaradas en sus respectivos .h
// #include "newtraphAlpha.h"
// #include "newtraphBeta.h"
// #include "newtraphBoth.h"

using namespace arma;
using namespace Rcpp;

//' @title Master Newton-Raphson Optimizer
//' @description Coordinates the optimization process for the intensity parameters by alternating between Alpha and Beta.
//' @param s Current time point.
//' @param n Number of subjects.
//' @param nvar Number of covariates.
//' @param k Number of observations per subject.
//' @param nk Total number of observations across all subjects.
//' @param tau Vector of times.
//' @param caltimes Calendar times.
//' @param gaptimes Gap times.
//' @param censored Censoring indicator.
//' @param intercepts Intercepts for effective age.
//' @param slopes Slopes for effective age.
//' @param lastperrep Last perfect repair.
//' @param perrepind Perfect repair indicator.
//' @param effagebegin Effective age at start.
//' @param effage Effective age at end.
//' @param cov Matrix of covariates.
//' @param alphaSeed Initial value for alpha.
//' @param betaSeed Initial vector for betas.
//' @param Z Frailty vector.
//' @param offset Offset vector.
//' @param rhoFunc Type of rho function.
//' @param ns Internal subject vector.
//' @param tol Convergence tolerance.
//' @param maxiter Maximum number of iterations.
//' @param loglik Output: Final log-likelihood.
//' @param estiEnd Output: Vector of final estimates (alpha + betas).
//' @param info Output: Final information matrix.
//' @param searchBoth Output: Status of the final simultaneous optimization.
//' @param kkiter Output: Number of iterations performed.
//' @export
// [[Rcpp::export]]
void newtraph(double s, int n, int nvar, arma::ivec k, int nk, 
              arma::vec tau, arma::vec caltimes, arma::vec gaptimes,
              arma::vec censored, arma::vec intercepts, arma::vec slopes,
              arma::vec lastperrep, arma::vec perrepind,
              arma::vec effagebegin, arma::vec effage, arma::mat cov,
              double alphaSeed, arma::vec betaSeed, arma::vec Z, 
              arma::vec offset, int rhoFunc, arma::ivec ns, double tol, 
              int maxiter, double& loglik, arma::vec& estiEnd, 
              arma::mat& info, int& searchBoth, int& kkiter) {
    
    double tol1 = 100.0 * tol;
    double distance = 1000.0;
    kkiter = 1;
    
    double alphaOld = alphaSeed;
    vec betaOld = betaSeed; 
    
    int search = 0;
    vec estiOld(nvar + 1);
    vec estiNew(nvar + 1);
    
    double estiNewAlpha;
    vec estiNewBeta(nvar);
    int searchAlpha, searchBeta, kkitera, kkiterb, kkiterab;
    
    // Newton-Raphson procedure alpha and beta 
    while ((distance > tol1) && (kkiter < maxiter)) {
        
        kkiter++;
        
        // Guardar estado previo (equivalente a estiOld(1)=alphaOld...)
        estiOld(0) = alphaOld;
        for(int i = 0; i < nvar; ++i) {
            estiOld(i + 1) = betaOld(i);
        }
        
        // Maximum Alpha 
        newtraphAlpha(s, n, nvar, k, nk, tau, caltimes, gaptimes, censored, 
                      intercepts, slopes, lastperrep, perrepind, effagebegin, 
                      effage, cov, alphaOld, betaOld, Z, offset, rhoFunc, 
                      ns, tol1, maxiter, estiNewAlpha, searchAlpha, kkitera);
        
        alphaOld = estiNewAlpha;
        
        if (searchAlpha != 1) {
            search = 121;
        }
        
        // Maximum Beta
        newtraphBeta(s, n, nvar, k, nk, tau, caltimes, gaptimes, censored, 
                     intercepts, slopes, lastperrep, perrepind, effagebegin, 
                     effage, cov, alphaOld, betaOld, Z, offset, rhoFunc, 
                     ns, tol1, maxiter, estiNewBeta, searchBeta, kkiterb);
        
        for(int i = 0; i < nvar; ++i) {
            betaOld(i) = estiNewBeta(i);
        }
        
        if (searchBeta != 1) {
            search = 122;
        }
        
        // Estado nuevo
        estiNew(0) = alphaOld;
        for(int i = 0; i < nvar; ++i) {
            estiNew(i + 1) = betaOld(i);
        }
        
        // stopping rule (Euclidean distance)
        distance = arma::norm(estiNew - estiOld, 2);
    }
    
    if (distance <= tol1) {
        search = 1;
    }
    
    // Maximum both
    if (search == 1) {
        // En Fortran: alpha = estiNew(1), beta = estiNew(2:nvar+1)
        double finalAlpha = estiNew(0);
        vec finalBeta = estiNew.subvec(1, nvar);
        
        newtraphBoth(s, n, nvar, k, nk, tau, caltimes, gaptimes, censored, 
                     intercepts, slopes, lastperrep, perrepind, effagebegin, 
                     effage, cov, finalAlpha, finalBeta, Z, 
                     offset, rhoFunc, ns, tol, maxiter, loglik, estiEnd, 
                     info, searchBoth, kkiterab);
    }
}