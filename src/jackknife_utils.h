#ifndef JACKKNIFE_UTILS_H
#define JACKKNIFE_UTILS_H

#include <RcppArmadillo.h>

/**
 * @title Jackknife Utilities
 * @description Helper functions to perform leave-one-out data extraction.
 */

// Removes one element from an integer vector
arma::ivec remove_unit(const arma::ivec& x, int pos);

// Removes one element from a double vector
arma::vec remove_unit(const arma::vec& x, int pos);

// Removes a range of observations (rows or columns) from a vector
arma::vec remove_range(const arma::vec& x, int start, int stop);

// Removes a range of columns from a covariate matrix
arma::mat remove_range_mat(const arma::mat& x, int start, int stop);

#endif