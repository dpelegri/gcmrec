#include <RcppArmadillo.h>
#include "scorefunc.h" 

using namespace arma;
using namespace Rcpp;

//' @title Newton-Raphson for Beta coefficients (Fixed Alpha)
//' @description Optimizes the beta vector using the information matrix.
//' @inheritParams Jacknife
//' @return A list with: loglik, estim (betas), info (infobeta), search, and kiter.
//' @export
// [[Rcpp::export]]
Rcpp::List nrBeta(double s, int n, int nvar, const arma::ivec& k, int nk, const arma::vec& tau,
                  const arma::vec& caltimes, const arma::vec& gaptimes, const arma::vec& censored,
                  const arma::vec& intercepts, const arma::vec& slopes, const arma::vec& lastperrep,
                  const arma::vec& perrepind, const arma::vec& effagebegin, const arma::vec& effage,
                  const arma::mat& cov, double alphaSeed, const arma::vec& betaSeed, 
                  const arma::vec& Z, const arma::vec& offset, int rhoFunc,
                  const arma::ivec& ns, double tol, int maxiter, double& loglik, 
                  arma::vec& estiNew, arma::mat& infobeta, int& search, int& kiter) 
{
    
    double distance = 1000.0;
    // int kiter = 1; // Empezamos en 1 como en Fortran
    // int search = 0;
    // double loglik = 0.0;
    kiter = 1;        // Correcto: ya existe, solo le das valor
    search = 0;       // Correcto
    loglik = 0.0;     // Correcto
    
    
    vec betaOld = betaSeed;
    // vec estiNew = betaSeed;
    estiNew = betaSeed; // Correcto (asumiendo que betaSeed es un argumento)
    
    
    // Variables de trabajo para scorefunc
    vec score(nvar + 1);
    mat info(nvar + 1, nvar + 1);
    double scorealpha, infoalpha;
    vec scorebeta(nvar);
    // mat infobeta(nvar, nvar);
    infobeta.zeros(nvar, nvar); // Correcto: la limpias con ceros
    
    while ((distance > tol) && (kiter < maxiter)) {
        kiter++;
        
        scorefunc(s, n, nvar, k, nk, tau, caltimes, gaptimes, censored,
                  intercepts, slopes, lastperrep, perrepind, effagebegin, effage,
                  cov, alphaSeed, betaOld, Z, offset, rhoFunc, ns,
                  loglik, score, info, scorealpha, scorebeta, infoalpha, infobeta);
        
        vec estiOld = betaOld;
        
        // Equivale a inv(infobeta) * scorebeta pero más estable
        vec DecDirec;
        bool success = arma::solve(DecDirec, infobeta, scorebeta);
        
        if(!success) {
            search = 0; // Error de convergencia por matriz singular
            break;
        }
        
        estiNew = estiOld + DecDirec;
        distance = arma::norm(estiNew - estiOld, 2);
        betaOld = estiNew;
    }
    
    if (distance <= tol) search = 1;
    
    return List::create(
        Named("loglik") = loglik,
        Named("estim")  = estiNew,
        Named("info")   = infobeta,
        Named("search") = search,
        Named("kiter")  = kiter
    );
}

//' @title Newton-Raphson for Alpha (Fixed Betas)
//' @description Optimizes the alpha parameter (scalar) using the information.
//' @inheritParams Jacknife
//' @return A list with: loglik, alphaNew, infoalpha, search, and kiter.
//' @export
// [[Rcpp::export]]
Rcpp::List nrAlpha(double s, int n, int nvar, const arma::ivec& k, int nk, const arma::vec& tau,
             const arma::vec& caltimes, const arma::vec& gaptimes, const arma::vec& censored,
             const arma::vec& intercepts, const arma::vec& slopes, const arma::vec& lastperrep,
             const arma::vec& perrepind, const arma::vec& effagebegin, const arma::vec& effage,
             const arma::mat& cov, double alphaSeed, const arma::vec& betaSeed, 
             const arma::vec& Z, const arma::vec& offset, int rhoFunc,
             const arma::ivec& ns, double tol, int maxiter, double& loglik, 
             double& alphaNew, double& infoalpha, int& search, int& kiter) {
    
    double distance = 1000.0;
    kiter = 1;
    search = 0;
    loglik = 0.0;
    double alphaOld = alphaSeed;
    alphaNew = alphaSeed;
    
    vec score(nvar + 1);
    mat info(nvar + 1, nvar + 1);
    double scorealpha = 0.0;
    infoalpha = 0.0;
    vec scorebeta(nvar);
    mat infobeta(nvar, nvar);
    
    while ((distance > tol) && (kiter < maxiter)) {
        kiter++;
        
        scorefunc(s, n, nvar, k, nk, tau, caltimes, gaptimes, censored,
                  intercepts, slopes, lastperrep, perrepind, effagebegin, effage,
                  cov, alphaOld, betaSeed, Z, offset, rhoFunc, ns,
                  loglik, score, info, scorealpha, scorebeta, infoalpha, infobeta);
        
        double estiOld = alphaOld;
        
        // Paso de Newton escalar: alpha_new = alpha_old + score / info
        alphaNew = estiOld + (scorealpha / infoalpha);
        
        distance = std::abs(alphaNew - estiOld);
        alphaOld = alphaNew;
    }
    
    if (distance <= tol) search = 1;
    
    return List::create(
        Named("loglik") = loglik,
        Named("estim")  = alphaNew, // Lo llamo estim para consistencia en R
        Named("info")   = infoalpha,
        Named("search") = search,
        Named("kiter")  = kiter
    );
}