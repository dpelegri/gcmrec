#ifndef SCOREFUNC_H
#define SCOREFUNC_H

#include <RcppArmadillo.h>

// Forward declaration of nsm if it's in another file
void nsm(double s, int n, int nk, const arma::vec& caltimes, const arma::ivec& k, arma::ivec& ns);

// Main score function exported to R
void scorefunc(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, arma::vec caltimes, 
               arma::vec gaptimes, arma::vec censored, arma::vec intercepts, arma::vec slopes, 
               arma::vec lastperrep, arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
               arma::mat cov, double alpha, arma::vec beta, arma::vec Z, arma::vec offset, int rhoFunc, 
               arma::ivec ns, double& loglik, arma::vec& score, arma::mat& info, 
               double& scorealpha, arma::vec& scorebeta, double& infoalpha, arma::mat& infobeta);

#endif