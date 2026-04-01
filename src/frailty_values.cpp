#include "frailty_values.h"

using namespace arma;

/**
 * @title FrailtyValues
 * @description Estimates individual frailty values (Z). 
 * Adjusts K by -1 as per the original Fortran logic.
 */
void FrailtyValues(int n, double xi, const ivec& KK_in, const vec& AA, vec& ZNew) {
    
    const double BIG_VALUE = 1.0e60;
    
    // According to Fortran: KK = K - 1
    // We can do this in one line with Armadillo
    ivec KK = KK_in - 1;
    
    if (xi >= BIG_VALUE) {
        // If xi is huge, frailties are non-existent (all equal to 1)
        ZNew.fill(1.0);
    } else {
        // Vectorized calculation: Z = (xi + KK) / (xi + AA)
        // conv_to<vec>::from converts integer vector to double for the division
        ZNew = (xi + conv_to<vec>::from(KK)) / (xi + AA);
    }
}