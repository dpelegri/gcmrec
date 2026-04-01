#ifndef NR_OPTIMIZERS_H
#define NR_OPTIMIZERS_H

#include <RcppArmadillo.h>

// Optimization for Beta only
void nrBeta(double s, int n, int nvar, const arma::ivec& k, int nk, const arma::vec& tau,
            const arma::vec& caltimes, const arma::vec& gaptimes, const arma::vec& censored,
            const arma::vec& intercepts, const arma::vec& slopes, const arma::vec& lastperrep,
            const arma::vec& perrepind, const arma::vec& effagebegin, const arma::vec& effage,
            const arma::mat& cov, double alphaSeed, const arma::vec& betaSeed, 
            const arma::vec& Z, const arma::vec& offset, int rhoFunc,
            const arma::ivec& ns, double tol, int maxiter, double& loglik, 
            arma::vec& estiNew, arma::mat& infobeta, int& search, int& kiter);

// Optimization for Alpha only
void nrAlpha(double s, int n, int nvar, const arma::ivec& k, int nk, const arma::vec& tau,
             const arma::vec& caltimes, const arma::vec& gaptimes, const arma::vec& censored,
             const arma::vec& intercepts, const arma::vec& slopes, const arma::vec& lastperrep,
             const arma::vec& perrepind, const arma::vec& effagebegin, const arma::vec& effage,
             const arma::mat& cov, double alphaSeed, const arma::vec& betaSeed, 
             const arma::vec& Z, const arma::vec& offset, int rhoFunc,
             const arma::ivec& ns, double tol, int maxiter, double& loglik, 
             double& alphaNew, double& infoalpha, int& search, int& kiter);

#endif