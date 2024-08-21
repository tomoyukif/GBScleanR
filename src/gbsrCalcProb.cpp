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
                            const bool & het,
                            const int & ploidy){
    vector<double> prob(ploidy + 1);
    double d_ploidy = (double)ploidy;
    double interval = eseq0 - eseq1;
    double denomi;

    double ref_multiplier;
    double alt_multiplier;

    for(int g = 0; g < prob.size(); ++g){
        double d_g = (double)g;
        denomi = w1 * (d_ploidy - d_g) + (1 - w1) * d_g;
        ref_multiplier = w1 * (d_ploidy - d_g) / denomi;
        ref_multiplier = ref_multiplier * interval + eseq1;
        alt_multiplier = 1 - ref_multiplier;
        ref_multiplier = log10_safe_d(ref_multiplier);
        alt_multiplier = log10_safe_d(alt_multiplier);
        prob[g] = ref * ref_multiplier + alt * alt_multiplier;
    }

    if(!het){
        for(int g = 0; g < prob.size(); ++g){
            if(g != 0 & g != ploidy){
                prob[g] = -numeric_limits<double>::infinity();
            }
        }
    }

    lognorm_vec(prob);
    return prob;
}

// Calculate mismap accounted genotype probabilities
void calcMissmap(vector<double> & prob,
                 const double & mismap1,
                 const double & mismap2,
                 const bool & het,
                 const int & ploidy){
    vector<double> mismap_prob(prob.size());
    vector<double> sum_prob(prob.size());
    double prob_lowest = -2.30103; //log10(0.005)
    double n_het = prob.size() - 2;
    bool check = false;

    for(size_t g = 0; g < prob.size(); ++g){
        if(g == 0){
            for(size_t h = 0; h < prob.size(); ++h){
                if(h == 0){
                    mismap_prob[h] = 1 - mismap1;

                } else {
                    mismap_prob[h] = 0;
                }
            }

        } else if(g == ploidy){
            for(size_t h = 0; h < prob.size(); ++h){
                if(h == ploidy){
                    mismap_prob[h] = 1 - mismap2;

                } else {
                    mismap_prob[h] = 0;
                }
            }

        } else {
            for(size_t h = 0; h < prob.size(); ++h){
                if(h == 0){
                    mismap_prob[h] = mismap1 / n_het;

                } else if(h == ploidy){
                    mismap_prob[h] = mismap2 / n_het;

                } else if(h == g){
                    mismap_prob[h] = 1;

                } else {
                    mismap_prob[h] = 0;
                }
            }
        }

        if(g == 0){
            for(size_t h = 0; h < prob.size(); ++h){
                log10_safe(mismap_prob[h]);
                mismap_prob[h] = mismap_prob[h] + prob[g];
                sum_prob[h] = mismap_prob[h];
            }

        } else {
            for(size_t h = 0; h < prob.size(); ++h){
                log10_safe(mismap_prob[h]);
                mismap_prob[h] = mismap_prob[h] + prob[g];
                logsum2(sum_prob[h], mismap_prob[h]);
            }
        }
    }

    double all_sum = logsum(sum_prob);
    if(all_sum == 0){
        if(het){
            double even_prob = 1/(double)prob.size();
            log10_safe(even_prob);
            prob.assign(prob.size(), even_prob);

        } else {
            prob.assign(prob.size(), -pow(10, 100));
            for(size_t g = 0; g < prob.size(); ++g){
                if(g == 0){
                    prob[g] = -0.30103;

                } else if(g == ploidy){
                    prob[g] = -0.30103;
                }
            }
        }

    } else {
        lognorm_vec(sum_prob);

        for(size_t g = 0; g < prob.size(); ++g){
            if(!check){
                check = sum_prob[g] < prob_lowest;
            }
        }

        if(check){
            for(size_t g = 0; g < prob.size(); ++g){
                logsum2(sum_prob[g], prob_lowest);
            }
            lognorm_vec(sum_prob);
        }

        for(size_t g = 0; g < prob.size(); ++g){
            prob[g] = sum_prob[g];
        }
    }
}

// Function to calculate probabilities of founder genotype patterns.
NumericVector calcPemit(NumericMatrix p_ref,
                        NumericMatrix p_alt,
                        NumericVector eseq,
                        NumericVector w1,
                        NumericVector mismap1,
                        NumericVector mismap2,
                        IntegerVector possiblegeno,
                        int & m,
                        IntegerVector n_f,
                        IntegerVector n_p,
                        LogicalVector het,
                        int ploidy){
    vector<double> prob;
    double p_prob;
    int col_i;
    NumericVector p_emit(n_p[0]);

    for(int i = 0; i < n_f[0]; ++i){
        NumericMatrix::Row ref_i = p_ref.row(i);
        NumericMatrix::Row alt_i = p_alt.row(i);

        prob = calcGenoprob(ref_i[m], alt_i[m],
                            eseq[0], eseq[1],
                                         w1[m], het[0], ploidy);
        calcMissmap(prob, mismap1[m], mismap2[m], het[0], ploidy);
        for(int j = 0; j < n_p[0]; ++j){
            col_i = j * n_f[0] + i;
            p_prob = prob[possiblegeno[col_i]];

            if(p_prob < -2){ // log10(0.01)
                p_prob = -numeric_limits<double>::infinity();
            }
            p_emit[j] = p_emit[j] + p_prob;
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
                        RVector<double> mismap1,
                        RVector<double> mismap2,
                        int m,
                        int & sample_i,
                        bool & het,
                        int ploidy){
    vector<double> prob(ploidy + 1);
    RMatrix<double>::Row ref_i = ref.row(sample_i);
    RMatrix<double>::Row alt_i = alt.row(sample_i);

    prob = calcGenoprob(ref_i[m], alt_i[m], eseq[0], eseq[1], w1[m], het, ploidy);
    calcMissmap(prob, mismap1[m], mismap2[m], het, ploidy);

    return prob;
}
