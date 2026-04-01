#include "nsm.h"

using namespace arma;

/**
 * @title nism logic (n_i^\dagger(s-))
 * @description Implementation of the Fortran routine that computes the number of 
 * events strictly before time s for a specific subject.
 * @param s Evaluation time.
 * @param caltimes Vector of calendar times for the subject.
 * @param k Number of observations for this subject.
 * @return Integer count of events.
 */
int nism(double s, const arma::vec& caltimes, int k)
{
    // Índices en C++ empiezan en 0, por eso caltimes(k) es caltimes[k-1]
    if (s > caltimes[k - 1]) {
        return k;
    } else {
        int count = 0;
        int i = 0;
        // Lógica de Fortran: do while (caltimes(i) < s)
        while (i < k && caltimes[i] < s) {
            count++;
            i++;
        }
        return count;
    }
}

//' @title Count Events per Subject (nsm)
//' @description Master routine to fill the event count vector for all subjects.
//' @param s Current evaluation time.
//' @param n Number of subjects.
//' @param nk Total number of observations in caltimes.
//' @param caltimes Global vector of calendar times.
//' @param k Vector with the number of observations for each subject.
//' @param nm [Out] Vector to store the results (integer).
//' @export
 // [[Rcpp::export]]
void nsm(double s, int n, int nk, const arma::vec& caltimes, const arma::ivec& k, arma::ivec& nm)
{
    
    int pos = 0;
    for (int i = 0; i < n; ++i) {
        int ki = k[i];
        
        // Extraemos la vista dinámica de los caltimes del sujeto actual
        vec caltimesOK = caltimes.subvec(pos, pos + ki - 1);
        
        // Aplicamos la lógica nism fiel al código original
        nm[i] = nism(s, caltimesOK, ki);
        
        pos += ki;
    }
}