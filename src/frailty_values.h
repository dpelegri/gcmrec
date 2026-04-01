#ifndef FRAILTY_VALUES_H
#define FRAILTY_VALUES_H

#include <RcppArmadillo.h>

/**
 * @title Compute Frailty Values (E-Step)
 * @description Updates the frailty estimates (Z) based on the current xi, KK, and AA.
 */
void FrailtyValues(int n, double xi, const arma::ivec& KK, const arma::vec& AA, arma::vec& ZNew);

#endif