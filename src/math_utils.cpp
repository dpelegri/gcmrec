#include "math_utils.h"
#include <cmath>

using namespace arma;

/**
 * @title Rho function
 * @description Computes the intensity scaling factor based on the model type.
 */
double rho(int k, double alpha, int rhoFunc) {
    if (rhoFunc == 1) return 1.0;
    if (rhoFunc == 2) return std::pow(alpha, k);
    return 1.0;
}

/**
 * @title Psi function
 * @description Exponential link function for the intensity.
 */
double psi(int nvar, const vec& cov, const vec& beta, double offset) {
    return std::exp(arma::dot(cov, beta) + offset);
}