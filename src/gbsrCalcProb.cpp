// [[depends(RcppParallel)]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include "gbsrutil.h"
using namespace Rcpp;
using namespace RcppParallel;
using namespace std;

vector<double> calcGenoprob(const double & ref,
                            const double & alt,
                            const double & eseq0,
                            const double & eseq1,
                            const double & w1,
                            const double & w2,
                            const bool & het){
    vector<double> prob(3);

    double logeseq0 = log10_safe_d(eseq0);
    double logeseq1 = log10_safe_d(eseq1);
    double logw1 = log10_safe_d(w1);
    double logw2 = log10_safe_d(w2);
    vector<double> ref_multiplier = {logeseq0, logw1, logeseq1};
    vector<double> alt_multiplier = {logeseq1, logw2, logeseq0};
    for(int g = 0; g < 3; ++g){
        prob[g] = ref * ref_multiplier[g] +
            alt * alt_multiplier[g];
    }

    if(!het){ prob[1] = 0; }

    lognorm_vec(prob);
    for(int g = 0; g < 3; ++g){
        prob[g] = pow(10, prob[g]);
    }
    return prob;
}

// Calculate mismap accounted genotype probabilities
void calcMissmap(vector<double> & prob,
                 const double & mismap1,
                 const double & mismap2,
                 const bool & het){
    vector<double> v1 = {1 - mismap1, mismap1, 0};
    vector<double> v2 = {0, 1, 0};
    vector<double> v3 = {0, mismap2, 1 - mismap2};
    double sum_v1 = 0.0;
    double sum_v2 = 0.0;
    double sum_v3 = 0.0;
    double sum_v = 0.0;
    double prob_lowest = 0.005;

    for(size_t g = 0; g < 3; ++g){
        sum_v1 += v1[g] * prob[g];
        sum_v2 += v2[g] * prob[g];
        sum_v3 += v3[g] * prob[g];
    }

    sum_v = sum_v1 + sum_v2 + sum_v3;
    if(sum_v == 0){
        if(het){
            prob.assign(3, 1/3);
        } else {
            prob = {0.5, 0, 0.5};
        }

    } else {
        if(sum_v1 < prob_lowest){
            sum_v1 += prob_lowest;
            sum_v2 += prob_lowest;
            sum_v3 += prob_lowest;
        }
        if(sum_v2 < prob_lowest){
            sum_v1 += prob_lowest;
            sum_v2 += prob_lowest;
            sum_v3 += prob_lowest;
        }
        if(sum_v3 < prob_lowest){
            sum_v1 += prob_lowest;
            sum_v2 += prob_lowest;
            sum_v3 += prob_lowest;
        }
        prob[0] = sum_v1 / sum_v;
        prob[1] = sum_v2 / sum_v;
        prob[2] = sum_v3 / sum_v;
    }
}


// Function to calculate probabilities of founder genotype patterns.
NumericVector calcPemit(NumericMatrix p_ref,
                        NumericMatrix p_alt,
                        NumericVector eseq,
                        NumericVector w1,
                        NumericVector w2,
                        NumericVector mismap1,
                        NumericVector mismap2,
                        IntegerVector possiblegeno,
                        int & m,
                        IntegerVector n_f,
                        IntegerVector n_p,
                        LogicalVector het,
                        IntegerVector ploidy
){
    vector<double> prob;
    double p_prob;
    int col_i;
    NumericVector p_emit(n_p[0], 1.0);

    for(int i = 0; i < n_f[0]; ++i){
        NumericMatrix::Row ref_i = p_ref.row(i);
        NumericMatrix::Row alt_i = p_alt.row(i);

        prob = calcGenoprob(ref_i[m], alt_i[m],
                            eseq[0], eseq[1],
                                         w1[m], w2[m], het[0]);
        calcMissmap(prob, mismap1[m], mismap2[m], het[0]);
        for(int j = 0; j < n_p[0]; ++j){
            col_i = j * n_f[0] + i;
            p_prob = prob[possiblegeno[col_i]];

            if(p_prob < 0.01){
                p_prob = 0;
            }
            p_emit[j] = p_emit[j] * p_prob;

        }
    }

    for(int j = 0; j < n_p[0]; ++j){
        if(p_emit[j] == 0){
            double neg_inf = -numeric_limits<double>::infinity();
            p_emit[j] = neg_inf;
        } else {
            log10_safe(p_emit[j]);
        }
    }
    p_emit = lognorm(p_emit);

    return p_emit;
}

// Function to calculate probabilities of founder genotype patterns.
vector<double> calcEmit(RMatrix<double> ref,
                        RMatrix<double> alt,
                        RVector<double> eseq,
                        RVector<double> w1,
                        RVector<double> w2,
                        RVector<double> mismap1,
                        RVector<double> mismap2,
                        int m,
                        int & sample_i,
                        bool & het
){
    vector<double> prob;
    RMatrix<double>::Row ref_i = ref.row(sample_i);
    RMatrix<double>::Row alt_i = alt.row(sample_i);

    prob = calcGenoprob(ref_i[m], alt_i[m], eseq[0], eseq[1], w1[m], w2[m], het);
    calcMissmap(prob, mismap1[m], mismap2[m], het);

    for(size_t i = 0; i < prob.size(); ++i){
        log10_safe(prob[i]);
    }

    return prob;
}
