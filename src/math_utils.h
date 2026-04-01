#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <RcppArmadillo.h>

double rho(int k, double alpha, int rhoFunc);
double psi(int nvar, const arma::vec& cov, const arma::vec& beta, double offset);

#endif