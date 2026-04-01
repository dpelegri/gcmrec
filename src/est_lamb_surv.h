#ifndef EST_LAMB_SURV_H
#define EST_LAMB_SURV_H

#include <RcppArmadillo.h>

/**
 * @title EstLambSurv Header
 * @description Provides the interface for non-parametric hazard and survival estimation.
 * Note: Returning Rcpp::List for direct access from R.
 */
// Rcpp::List EstLambSurv(double s, int n, int nvar, const arma::ivec& k, int nk, 
//                        const arma::vec& tau, const arma::vec& caltimes, 
//                        const arma::vec& gaptimes, const arma::vec& censored,
//                        const arma::vec& intercepts, const arma::vec& slopes, 
//                        const arma::vec& lastperrep, const arma::vec& perrepind, 
//                        const arma::vec& effagebegin, const arma::vec& effage,
//                        int ndiseff, const arma::vec& diseff, const arma::mat& cov, 
//                        double alpha, const arma::vec& beta, const arma::vec& Z, 
//                        const arma::vec& offset, int rhoFunc);

// Rcpp::List EstLambSurv(double s, int n, int nvar, const arma::ivec& k, int nk, 
//                        const arma::vec& tau, const arma::vec& caltimes, 
//                        const arma::vec& gaptimes, const arma::vec& censored, 
//                        const arma::vec& intercepts, const arma::vec& slopes, 
//                        const arma::vec& lastperrep, const arma::vec& perrepind, 
//                        const arma::vec& effagebegin, const arma::vec& effage, 
//                        int ndiseff, const arma::vec& diseff, const arma::mat& cov, 
//                        double alpha, const arma::vec& beta, const arma::vec& Z, 
//                        const arma::vec& offset, int rhoFunc,
//                        arma::vec& lambda, arma::vec& deltalambda, arma::vec& surv);

Rcpp::List EstLambSurv(double s, int n, int nvar, const arma::ivec& k, int nk, 
                       const arma::vec& tau, const arma::vec& caltimes, 
                       const arma::vec& gaptimes, const arma::vec& censored, 
                       const arma::vec& intercepts, const arma::vec& slopes, 
                       const arma::vec& lastperrep, const arma::vec& perrepind, 
                       const arma::vec& effagebegin, const arma::vec& effage, 
                       int ndiseff, const arma::vec& diseff, const arma::mat& cov, 
                       double alpha, const arma::vec& beta, const arma::vec& Z, 
                       const arma::vec& offset, int rhoFunc,
                       arma::vec& lambda, arma::vec& deltalambda, arma::vec& surv);

#endif