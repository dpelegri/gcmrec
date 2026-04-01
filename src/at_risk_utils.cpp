#include "at_risk_utils.h"
#include "math_utils.h"

using namespace arma;



/**
 * @title Internal Subject Risk Calculation
 * @description Computes indicators (Q, R) and their derivatives for a single subject.
 * @param s Evaluation time.
 * @param nsi Number of events for the subject.
 * @param w Effective age.
 * @param nvar Number of covariates.
 * @param ki Number of observations.
 * @param tau_i Individual censoring time.
 * @param caltimes Vector of calendar times.
 * @param gaptimes Vector of gap times.
 * @param censored Censoring indicator.
 * @param intercepts Model intercepts.
 * @param slopes Model slopes.
 * @param lastperrep Index of last periodic report.
 * @param perrepind Periodic report indicator.
 * @param effagebegin Effective age at start.
 * @param effage Effective age at end.
 * @param cov Matrix of covariates.
 * @param alpha Recurrence parameter.
 * @param beta Vector of regression coefficients.
 * @param rhoFunc Type of rho function.
 * @param offset Subject-specific offset.
 * @param Ysubj [Out] Subject's risk indicator.
 * @param GrYsubjA1 [Out] Gradient wrt alpha.
 * @param GrYsubjBe [Out] Gradient wrt beta.
 * @param Gr2YsubjA1 [Out] Second derivative wrt alpha.
 * @param Gr2YsubjA1Be [Out] Cross derivative alpha-beta.
 * @param Gr2YsubjBe [Out] Second derivative matrix wrt beta.
 */
void AtRiskSubj(double s, int nsi, double w, int nvar, int ki, double tau_i,
                const vec& caltimes, const vec& gaptimes, double censored,
                const vec& intercepts, const vec& slopes, const vec& lastperrep,
                const vec& perrepind, const vec& effagebegin, const vec& effage,
                const mat& cov, double alpha, const vec& beta, int rhoFunc,
                double offset, double& Ysubj, double& GrYsubjA1, vec& GrYsubjBe,
                double& Gr2YsubjA1, vec& Gr2YsubjA1Be, mat& Gr2YsubjBe) {
    
    // Dynamic initialization based on nsi (number of events for subject i)
    vec Q(nsi, fill::zeros);
    double R = 0.0;
    vec GrQA1(nsi, fill::zeros);
    mat GrQBe(nvar, nsi, fill::zeros);
    double GrRA1 = 0.0, Gr2RA1 = 0.0;
    vec GrRBe(nvar, fill::zeros), Gr2RA1Be(nvar, fill::zeros);
    mat Gr2RBe(nvar, nvar, fill::zeros);
    
    mat Gr2QBe(nvar, nvar, fill::zeros); // Aggregated internally
    vec Gr2QA1Be_vec(nvar, fill::zeros);
    double Gr2QA1_sum = 0.0;
    
    // Equation (23) Q_ij
    if (nsi > 1) {
        for (int j = 1; j < nsi; ++j) { // 0-based: j=1 is the 2nd event
            if (w > effagebegin[j - 1] && w <= effage[j]) {
                vec covariate = cov.col(j - 1);
                
                Q[j] = (rho(j - 1, alpha, rhoFunc) * psi(nvar, covariate, beta, offset)) / slopes[j - 1];
                
                // Gradients based on rhoFunc logic
                if (rhoFunc == 1) {
                    GrQA1[j] = Q[j];
                    Gr2QA1_sum += Q[j];
                } else {
                    GrQA1[j] = (double(j - 1) / alpha) * Q[j];
                    Gr2QA1_sum += (double((j - 1) * (j - 2)) / std::pow(alpha, 2)) * Q[j];
                }
                
                for (int i = 0; i < nvar; ++i) {
                    GrQBe(i, j) = covariate[i] * Q[j];
                    if (rhoFunc == 1) {
                        Gr2QA1Be_vec[i] += covariate[i] * Q[j];
                    } else {
                        Gr2QA1Be_vec[i] += (double(j - 1) / alpha) * covariate[i] * Q[j];
                    }
                    for (int t = 0; t < nvar; ++t) {
                        Gr2QBe(i, t) += Q[j] * covariate[i] * covariate[t];
                    }
                }
            }
        }
    }
    
    // Equation (24) R_i
    // Note: lastperrep values from Fortran are used as indices (convert to 0-based)
    int last_idx = int(lastperrep[ki - 1]) - 1; 
    double effageats = intercepts[ki - 1] + slopes[ki - 1] * std::min(s, tau_i) - caltimes[last_idx];
    
    if (w > effagebegin[ki - 1] && w <= effageats) {
        vec covariate = cov.col(ki - 1);
        R = (rho(ki - 1, alpha, rhoFunc) * psi(nvar, covariate, beta, offset)) / slopes[ki - 1];
        
        if (rhoFunc == 1) {
            GrRA1 = R;
            Gr2RA1 = R;
        } else {
            GrRA1 = (double(ki - 1) / alpha) * R;
            Gr2RA1 = (double((ki - 1) * (ki - 2)) / std::pow(alpha, 2)) * R;
        }
        
        for (int i = 0; i < nvar; ++i) {
            GrRBe[i] = covariate[i] * R;
            Gr2RA1Be[i] = (rhoFunc == 1) ? (covariate[i] * R) : ((double(ki - 1) / alpha) * covariate[i] * R);
            for (int t = 0; t < nvar; ++t) {
                Gr2RBe(i, t) = R * covariate[i] * covariate[t];
            }
        }
    }
    
    // Output Assembly
    Ysubj = sum(Q) + R;
    GrYsubjA1 = sum(GrQA1) + GrRA1;
    Gr2YsubjA1 = Gr2QA1_sum + Gr2RA1;
    
    GrYsubjBe = sum(GrQBe, 1) + GrRBe;
    Gr2YsubjA1Be = Gr2QA1Be_vec + Gr2RA1Be;
    Gr2YsubjBe = Gr2QBe + Gr2RBe;
}

/**
 * @title Aggregate At-Risk Statistics
 * @description Sums the risk indicators and gradients across all subjects.
 * @param s Evaluation time.
 * @param ns Vector of event counts per subject.
 * @param w Effective age.
 * @param n Number of subjects.
 * @param nvar Number of covariates.
 * @param nk Total observations.
 * @param k Observations per subject.
 * @param tau Censoring times.
 * @param caltimes Calendar times.
 * @param gaptimes Gap times.
 * @param censored Censoring indicators.
 * @param intercepts Model intercepts.
 * @param slopes Model slopes.
 * @param lastperrep Periodic report indices.
 * @param perrepind Periodic report indicators.
 * @param effagebegin Effective age starts.
 * @param effage Effective age ends.
 * @param cov Covariate matrix.
 * @param alpha Recurrence parameter.
 * @param beta Regression coefficients.
 * @param Z Fragility/Random effect vector.
 * @param offset Offsets vector.
 * @param rhoFunc Rho function index.
 * @param Ysubj [Out] Vector of risk indicators per subject.
 * @param S0 [Out] Aggregate risk set value.
 * @param GrS0A1 [Out] Aggregate gradient wrt alpha.
 * @param GrS0Be [Out] Aggregate gradient wrt beta.
 * @param Gr2S0A1 [Out] Aggregate second derivative wrt alpha.
 * @param Gr2S0A1Be [Out] Aggregate cross derivative alpha-beta.
 * @param Gr2S0Be [Out] Aggregate second derivative matrix wrt beta.
 */
void AtRisk(double s, const ivec& ns, double w, int n, int nvar, int nk, 
            const ivec& k, const vec& tau, const vec& caltimes, const vec& gaptimes,
            const vec& censored, const vec& intercepts, const vec& slopes, 
            const vec& lastperrep, const vec& perrepind, const vec& effagebegin, 
            const vec& effage, const mat& cov, double alpha, const vec& beta, 
            const vec& Z, const vec& offset, int rhoFunc, vec& Ysubj, 
            double& S0, double& GrS0A1, vec& GrS0Be, double& Gr2S0A1, 
            vec& Gr2S0A1Be, mat& Gr2S0Be) {
    
    // Initialization
    S0 = 0.0; GrS0A1 = 0.0; Gr2S0A1 = 0.0;
    GrS0Be.zeros(nvar); Gr2S0A1Be.zeros(nvar); Gr2S0Be.zeros(nvar, nvar);
    Ysubj.zeros(n);
    
    int pos = 0;
    for (int i = 0; i < n; ++i) {
        int ki = k[i];
        
        // Dynamic extraction of subject data using sub-views
        vec calOK = caltimes.subvec(pos, pos + ki - 1);
        vec gapOK = gaptimes.subvec(pos, pos + ki - 1);
        vec intOK = intercepts.subvec(pos, pos + ki - 1);
        vec sloOK = slopes.subvec(pos, pos + ki - 1);
        vec lastOK = lastperrep.subvec(pos, pos + ki - 1);
        vec perOK = perrepind.subvec(pos, pos + ki - 1);
        vec effBOK = effagebegin.subvec(pos, pos + ki - 1);
        vec effOK = effage.subvec(pos, pos + ki - 1);
        mat covOK = cov.cols(pos, pos + ki - 1);
        
        double Ys_i, GrYA1_i, Gr2YA1_i;
        vec GrYBe_i(nvar), Gr2YA1Be_i(nvar);
        mat Gr2YBe_i(nvar, nvar);
        
        AtRiskSubj(s, ns[i], w, nvar, ki, tau[i], calOK, gapOK, censored[i],
                   intOK, sloOK, lastOK, perOK, effBOK, effOK, covOK,
                   alpha, beta, rhoFunc, offset[i], Ys_i, GrYA1_i, 
                   GrYBe_i, Gr2YA1_i, Gr2YA1Be_i, Gr2YBe_i);
        
        Ysubj[i] = Ys_i;
        double Zi = Z[i];
        
        // Global aggregation
        S0 += Zi * Ysubj[i];
        GrS0A1 += Zi * GrYA1_i;
        Gr2S0A1 += Zi * Gr2YA1_i;
        
        for (int j = 0; j < nvar; ++j) {
            GrS0Be[j] += Zi * GrYBe_i[j];
            Gr2S0A1Be[j] += Zi * Gr2YA1Be_i[j];
            for (int t = 0; t < nvar; ++t) {
                Gr2S0Be(j, t) += Zi * Gr2YBe_i(j, t);
            }
        }
        pos += ki;
    }
}