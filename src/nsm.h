#ifndef NSM_H
#define NSM_H

#include <RcppArmadillo.h>

/**
 * @title Internal Event Counter for Subject i
 * @description Computes n_i^\dagger(s-) logic.
 */
int nism(double s, const arma::vec& caltimes, int k);

/**
 * @title Master Event Counter
 * @description Iterates over all subjects to populate the event count vector.
 */
void nsm(double s, int n, int nk, const arma::vec& caltimes, const arma::ivec& k, arma::ivec& nm);

#endif