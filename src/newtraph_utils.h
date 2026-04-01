#ifndef NEWTRAPH_UTILS_H
#define NEWTRAPH_UTILS_H

#include <RcppArmadillo.h>

// Optimization for both parameters simultaneously
void newtraphBoth(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, 
                  arma::vec caltimes, arma::vec gaptimes, arma::vec censored, 
                  arma::vec intercepts, arma::vec slopes, arma::vec lastperrep, 
                  arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
                  arma::mat cov, double alphaSeed, arma::vec betaSeed, arma::vec Z, 
                  arma::vec offset, int rhoFunc, arma::ivec ns, double tol, 
                  int maxiter, double& loglik, arma::vec& estiNew, arma::mat& info, 
                  int& search, int& kiter);

// Optimization restricted to Alpha
void newtraphAlpha(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, 
                   arma::vec caltimes, arma::vec gaptimes, arma::vec censored, 
                   arma::vec intercepts, arma::vec slopes, arma::vec lastperrep, 
                   arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
                   arma::mat cov, double alphaSeed, arma::vec betaSeed, arma::vec Z, 
                   arma::vec offset, int rhoFunc, arma::ivec ns, double tol, 
                   int maxiter, double& estiNew, int& search, int& kiter);

// Optimization restricted to Beta
void newtraphBeta(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, 
                  arma::vec caltimes, arma::vec gaptimes, arma::vec censored, 
                  arma::vec intercepts, arma::vec slopes, arma::vec lastperrep, 
                  arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
                  arma::mat cov, double alphaSeed, arma::vec betaSeed, arma::vec Z, 
                  arma::vec offset, int rhoFunc, arma::ivec ns, double tol, 
                  int maxiter, arma::vec& estiNew, int& search, int& kiter);

#endif