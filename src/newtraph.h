#ifndef NEWTRAPH_H
#define NEWTRAPH_H

#include <RcppArmadillo.h>


void newtraph(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, 
              arma::vec caltimes, arma::vec gaptimes, arma::vec censored, 
              arma::vec intercepts, arma::vec slopes, arma::vec lastperrep, 
              arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
              arma::mat cov, double alphaSeed, arma::vec betaSeed, arma::vec Z, 
              arma::vec offset, int rhoFunc, arma::ivec ns, double tol, 
              int maxiter, double& loglik, arma::vec& estiNew, arma::mat& info, 
              int& search, int& kiter);

#endif