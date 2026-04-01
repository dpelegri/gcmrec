#ifndef JACKNIFE_H
#define JACKNIFE_H

#include <RcppArmadillo.h>

/**
 * @title Jacknife (Non-Frailty)
 */
arma::mat Jacknife(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, 
                   arma::vec caltimes, arma::vec gaptimes, arma::vec censored, 
                   arma::vec intercepts, arma::vec slopes, arma::vec lastperrep, 
                   arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
                   arma::mat cov, double alphaSeed, arma::vec betaSeed, 
                   arma::vec Z, arma::vec offset, int rhoFunc, arma::ivec ns, 
                   double tol, int maxiter);

/**
 * @title Jacknife2 (Frailty)
 */
arma::mat Jacknife2(double s, int n, int nvar, arma::ivec k, int nk, arma::vec tau, 
                    arma::vec caltimes, arma::vec gaptimes, arma::vec censored, 
                    arma::vec intercepts, arma::vec slopes, arma::vec lastperrep, 
                    arma::vec perrepind, arma::vec effagebegin, arma::vec effage, 
                    int ndiseff, arma::vec diseff, arma::mat cov, double alphaSeed, 
                    arma::vec betaSeed, double xi, arma::vec Z, arma::vec offset, 
                    int rhoFunc, double tol, int maxiter, int maxXi);

#endif