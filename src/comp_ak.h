#ifndef COMP_AK_H
#define COMP_AK_H

#include <RcppArmadillo.h>

/**
 * @title Compute AK, BK and KK coefficients
 * @description Calculates the necessary components for the E-step of the EM algorithm.
 */
void CompAK(double s, int n, int nvar, const arma::ivec& k, int nk, const arma::vec& tau,
            const arma::vec& caltimes, const arma::vec& gaptimes, const arma::vec& censored,
            const arma::vec& intercepts, const arma::vec& slopes, const arma::vec& lastperrep,
            const arma::vec& perrepind, const arma::vec& effagebegin, const arma::vec& effage,
            int ndiseff, const arma::vec& diseff, const arma::mat& cov, double alpha, 
            const arma::vec& beta, const arma::vec& deltalambda, const arma::vec& Z, 
            const arma::vec& offset, int rhoFunc, arma::ivec& KK, arma::vec& AA, arma::vec& BB);

#endif