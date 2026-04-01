#include "newtraph_utils.h"
#include "scorefunc.h"

using namespace arma;


//' @title Newton-Raphson for Alpha and Beta
//' @description Simultaneous optimization for alpha and beta parameters.
//' @keywords internal
//' @noRd
 // [[Rcpp::export]]
void newtraphBoth(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, 
                  arma::vec caltimes, arma::vec gaptimes, arma::vec censored, 
                  arma::vec intercepts, arma::vec slopes, arma::vec lastperrep, 
                  arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
                  arma::mat cov, double alphaSeed, arma::vec betaSeed, arma::vec Z, 
                  arma::vec offset, int rhoFunc, arma::ivec ns, double tol, 
                  int maxiter, double& loglik, arma::vec& estiNew, arma::mat& info, 
                  int& search, int& kiter)
{
    
    double distance = 1000.0;
    kiter = 1;
    double alphaOld = alphaSeed;
    vec betaOld = betaSeed;
    search = 0;
    
    // Buffers for scorefunc
    vec score(nvar + 1);
    double scorealpha, infoalpha;
    vec scorebeta(nvar);
    mat infobeta(nvar, nvar);
    vec estiOld(nvar + 1);
    
    while ((distance > tol) && (kiter < maxiter)) {
        kiter++;
        
        scorefunc(s, n, nvar, k, nk, tau, caltimes, gaptimes, censored, intercepts, 
                  slopes, lastperrep, perrepind, effagebegin, effage, cov, 
                  alphaOld, betaOld, Z, offset, rhoFunc, ns, 
                  loglik, score, info, scorealpha, scorebeta, infoalpha, infobeta);
        
        estiOld(0) = alphaOld;
        for (int i = 0; i < nvar; ++i) estiOld(i + 1) = betaOld(i);
        
        // Compute descent direction: Solve H * x = g
        vec DecDirec = arma::solve(info, score, arma::solve_opts::fast);
        
        estiNew = estiOld + DecDirec;
        
        // Euclidean distance using Armadillo norm
        distance = arma::norm(estiNew - estiOld, 2);
        
        alphaOld = estiNew(0);
        for (int i = 0; i < nvar; ++i) betaOld(i) = estiNew(i + 1);
        
        if (distance <= tol) search = 1;
    }
}

//' @title Newton-Raphson for Alpha
//' @description Optimization restricted to the alpha parameter.
//' @keywords internal
//' @noRd
 // [[Rcpp::export]]
void newtraphAlpha(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, 
                   arma::vec caltimes, arma::vec gaptimes, arma::vec censored, 
                   arma::vec intercepts, arma::vec slopes, arma::vec lastperrep, 
                   arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
                   arma::mat cov, double alphaSeed, arma::vec betaSeed, arma::vec Z, 
                   arma::vec offset, int rhoFunc, arma::ivec ns, double tol, 
                   int maxiter, double& estiNew, int& search, int& kiter)
{
    
    double distance = 1000.0;
    kiter = 1;
    double alphaOld = alphaSeed;
    search = 0;
    
    double loglik, scorealpha, infoalpha;
    vec score(nvar + 1), scorebeta(nvar);
    mat info(nvar + 1, nvar + 1), infobeta(nvar, nvar);
    
    while ((distance > tol) && (kiter < maxiter)) {
        kiter++;
        
        scorefunc(s, n, nvar, k, nk, tau, caltimes, gaptimes, censored, intercepts, 
                  slopes, lastperrep, perrepind, effagebegin, effage, cov, 
                  alphaOld, betaSeed, Z, offset, rhoFunc, ns, 
                  loglik, score, info, scorealpha, scorebeta, infoalpha, infobeta);
        
        double estiOld = alphaOld;
        estiNew = estiOld + (scorealpha / infoalpha); // Scalar update
        
        distance = std::abs(estiNew - estiOld);
        alphaOld = estiNew;
        
        if (distance <= tol) search = 1;
    }
}

//' @title Newton-Raphson for Beta
//' @description Optimization restricted to beta parameters.
//' @keywords internal
//' @noRd
 // [[Rcpp::export]]
void newtraphBeta(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, 
                  arma::vec caltimes, arma::vec gaptimes, arma::vec censored, 
                  arma::vec intercepts, arma::vec slopes, arma::vec lastperrep, 
                  arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
                  arma::mat cov, double alphaSeed, arma::vec betaSeed, arma::vec Z, 
                  arma::vec offset, int rhoFunc, arma::ivec ns, double tol, 
                  int maxiter, arma::vec& estiNew, int& search, int& kiter)
{
    
    double distance = 1000.0;
    kiter = 1;
    vec betaOld = betaSeed;
    search = 0;
    
    double loglik, scorealpha, infoalpha;
    vec score(nvar + 1), scorebeta(nvar);
    mat info(nvar + 1, nvar + 1), infobeta(nvar, nvar);
    
    while ((distance > tol) && (kiter < maxiter)) {
        kiter++;
        
        scorefunc(s, n, nvar, k, nk, tau, caltimes, gaptimes, censored, intercepts, 
                  slopes, lastperrep, perrepind, effagebegin, effage, cov, 
                  alphaSeed, betaOld, Z, offset, rhoFunc, ns, 
                  loglik, score, info, scorealpha, scorebeta, infoalpha, infobeta);
        
        vec estiOld = betaOld;
        
        // Solve linear system for Beta update
        vec DecDirec = arma::solve(infobeta, scorebeta, arma::solve_opts::fast);
        
        estiNew = estiOld + DecDirec;
        distance = arma::norm(estiNew - estiOld, 2);
        
        betaOld = estiNew;
        if (distance <= tol) search = 1;
    }
}