#ifndef MAX_XI_H
#define MAX_XI_H

#include <RcppArmadillo.h>

void MaxWrtXi(double xiOld, int n, const arma::ivec& K, const arma::vec& A, const arma::vec& B,
              double tol, double& xiNew, double& G0, int maxiter, int& search);

// Prototype for the LogLikXi function that we'll need next
void LogLikXi(double xi, int n, const arma::ivec& KK, const arma::vec& A, 
              const arma::vec& B, double& G0, double& G1, double& G2);

#endif