#ifndef AT_RISK_UTILS_H
#define AT_RISK_UTILS_H

#include <RcppArmadillo.h>

void AtRiskSubj(double s, int nsi, double w, int nvar, int ki, double tau_i,
                const arma::vec& caltimes, const arma::vec& gaptimes, double censored,
                const arma::vec& intercepts, const arma::vec& slopes, const arma::vec& lastperrep,
                const arma::vec& perrepind, const arma::vec& effagebegin, const arma::vec& effage,
                const arma::mat& cov, double alpha, const arma::vec& beta, int rhoFunc,
                double offset, double& Ysubj, double& GrYsubjA1, arma::vec& GrYsubjBe,
                double& Gr2YsubjA1, arma::vec& Gr2YsubjA1Be, arma::mat& Gr2YsubjBe);

void AtRisk(double s, const arma::ivec& ns, double w, int n, int nvar, int nk, 
            const arma::ivec& k, const arma::vec& tau, const arma::vec& caltimes, const arma::vec& gaptimes,
            const arma::vec& censored, const arma::vec& intercepts, const arma::vec& slopes, 
            const arma::vec& lastperrep, const arma::vec& perrepind, const arma::vec& effagebegin, 
            const arma::vec& effage, const arma::mat& cov, double alpha, const arma::vec& beta, 
            const arma::vec& Z, const arma::vec& offset, int rhoFunc, arma::vec& Ysubj, 
            double& S0, double& GrS0A1, arma::vec& GrS0Be, double& Gr2S0A1, 
            arma::vec& Gr2S0A1Be, arma::mat& Gr2S0Be);

#endif