#ifndef FRAILTY_ENGINE_H
#define FRAILTY_ENGINE_H

#include <RcppArmadillo.h>

    // Función principal exportada a R
    void EstimWithFrailty(double s, int n, int nvar, arma::ivec k, int nk, 
                          arma::vec tau, arma::vec caltimes, arma::vec gaptimes,
                          arma::vec censored, arma::vec intercepts, arma::vec slopes, 
                          arma::vec lastperrep, arma::vec perrepind, arma::vec effagebegin, 
                          arma::vec effage, int ndiseff, arma::vec diseff, arma::mat cov, 
                          double alphaSeed, arma::vec betaSeed, double xi, arma::vec Z, 
                          arma::vec offset, int rhoFunc, double tol, int maxiter, 
                          int maxXi, arma::vec& estimEnd, arma::ivec& control, double& loglikEnd);
    
    // // Prototipos de funciones auxiliares que implementaremos a continuación
    // void EstLambSurv(double s, int n, int nvar, const arma::ivec& k, int nk, const arma::vec& tau,
    //                  const arma::vec& caltimes, const arma::vec& gaptimes, const arma::vec& censored,
    //                  const arma::vec& intercepts, const arma::vec& slopes, const arma::vec& lastperrep,
    //                  const arma::vec& perrepind, const arma::vec& effagebegin, const arma::vec& effage,
    //                  int ndiseff, const arma::vec& diseff, const arma::mat& cov, double alpha, 
    //                  const arma::vec& beta, const arma::vec& Z, const arma::vec& offset, int rhoFunc,
    //                  arma::vec& lambda, arma::vec& deltalambda, arma::vec& surv);
    
    void CompAK(double s, int n, int nvar, const arma::ivec& k, int nk, const arma::vec& tau,
                const arma::vec& caltimes, const arma::vec& gaptimes, const arma::vec& censored,
                const arma::vec& intercepts, const arma::vec& slopes, const arma::vec& lastperrep,
                const arma::vec& perrepind, const arma::vec& effagebegin, const arma::vec& effage,
                int ndiseff, const arma::vec& diseff, const arma::mat& cov, double alpha, 
                const arma::vec& beta, const arma::vec& deltalambda, const arma::vec& Z, 
                const arma::vec& offset, int rhoFunc, arma::ivec& KK, arma::vec& AA, arma::vec& BB);
    
    void FrailtyValues(int n, double xi, const arma::ivec& KK, const arma::vec& AA, arma::vec& ZNew);
    
    void MaxWrtXi(double xiOld, int n, const arma::ivec& KK, const arma::vec& AA, const arma::vec& BB,
                  double tol, double& xiNew, double& loglikMarg, int maxiter, int& searchxi);
    
    // double GammaLogLikOptim(double xiOld, int n, const arma::ivec& KK, const arma::vec& AA, const arma::vec& BB);
    double GammaLogLikOptim(double seedXi, int n, const arma::ivec& K, const arma::vec& A, const arma::vec& B);
    
    double LogLikXi2(double xi, int n, const arma::ivec& K, const arma::vec& A, const arma::vec& B);
    
    void mnbrak(double& ax, double& bx, double& cx, double& fa, double& fb, double& fc,
                double (*func)(double, int, const arma::ivec&, const arma::vec&, const arma::vec&),
                int nn, const arma::ivec& KK, const arma::vec& AA, const arma::vec& BB);
    
    void brent(double ax, double bx, double cx,
               double (*f)(double, int, const arma::ivec&, const arma::vec&, const arma::vec&),
               double tol, double& xmin, int nn, const arma::ivec& KK, const arma::vec& AA, const arma::vec& BB);

#endif