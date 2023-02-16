// [[depends(RcppParallel)]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include "gbsrutil.h"
#include "gbsrcalcprob.h"
using namespace Rcpp;
using namespace RcppParallel;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Functions for the forward-backward algorithm

struct ParFB : public Worker {

    RMatrix<double> gamma;
    const RVector<int> iter_sample;
    const RMatrix<double> ref;
    const RMatrix<double> alt;
    const RVector<double> eseq;
    const RVector<double> w1;
    const RVector<double> w2;
    const RVector<double> mismap1;
    const RVector<double> mismap2;
    const RVector<int> possiblehap;
    const RVector<double> init_prob;
    const RMatrix<double> trans_prob;
    const RVector<int> dim;
    const RVector<int> p_geno;
    const RVector<int> ploidy;

    ParFB(NumericMatrix gamma,
          const LogicalVector iter_sample,
          const NumericMatrix ref,
          const NumericMatrix alt,
          const NumericVector eseq,
          const NumericVector w1,
          const NumericVector w2,
          const NumericVector mismap1,
          const NumericVector mismap2,
          const IntegerVector possiblehap,
          const NumericVector init_prob,
          const NumericMatrix trans_prob,
          const IntegerVector dim,
          const IntegerVector p_geno,
          const IntegerVector ploidy)
        : gamma(gamma),
          iter_sample(iter_sample),
          ref(ref),
          alt(alt),
          eseq(eseq),
          w1(w1),
          w2(w2),
          mismap1(mismap1),
          mismap2(mismap2),
          possiblehap(possiblehap),
          init_prob(init_prob),
          trans_prob(trans_prob),
          dim(dim),
          p_geno(p_geno),
          ploidy(ploidy) {}

    void operator()(size_t begin, size_t end) {
        bool het = true;
        for(RVector<int>::const_iterator i=iter_sample.begin() + begin;
            i<iter_sample.begin() + end; ++i){
            size_t sample_i = distance(iter_sample.begin(), i);

            RMatrix<double>::Row gamma_i = gamma.row(sample_i);

            size_t n_m = dim[0];
            size_t n_h = dim[2];

            vector<vector<double>> alpha(n_m,
                                         vector<double>(n_h));
            vector<vector<double>> emit(n_m,
                                        vector<double>(n_h));
            vector<double> beta(n_h, 0);

            vector<double> score_k(n_h);
            double hap_prob;
            double trans_kk;
            double sum_k;
            size_t col_i;
            size_t j;
            double neg_inf = -numeric_limits<double>::infinity();

            for(size_t m=0; m<n_m; ++m){
                vector<double> prob_i = calcEmit(ref,
                                                 alt,
                                                 eseq,
                                                 w1,
                                                 w2,
                                                 mismap1,
                                                 mismap2,
                                                 m,
                                                 sample_i,
                                                 het);
                if(m == 0){
                    j = p_geno[0];
                    for(size_t k=0; k<n_h; ++k){
                        col_i = j * n_h + k;
                        hap_prob = prob_i[possiblehap[col_i]];
                        alpha[m][k] = hap_prob + init_prob[k];
                    }

                } else {
                    size_t j = p_geno[m];

                    for(size_t k2=0; k2<n_h; ++k2){
                        RMatrix<double>::Column trans_prob_k =
                            trans_prob.column((m-1)*n_h + k2);
                        for(size_t k1=0; k1<n_h; ++k1){
                            trans_kk = trans_prob_k[k1];
                            score_k.at(k1) = alpha[m-1][k1] + trans_kk;
                        }
                        sum_k = logsum(score_k);
                        col_i = j * n_h + k2;
                        hap_prob = prob_i[possiblehap[col_i]];
                        emit[m][k2] = hap_prob;
                        alpha[m][k2] = hap_prob + sum_k;
                    }
                }
            }

            vector<double> gamma_i1{neg_inf};
            vector<double> gamma_i2{neg_inf};
            vector<double> gamma_i3{neg_inf};
            double gamma_tmp;
            vector<double> gamma_i_tmp(3);
            for(size_t k=0; k<n_h; ++k){
                size_t j = p_geno[n_m-1];
                col_i = j * n_h + k;
                gamma_tmp = beta[k] + alpha[n_m-1][k];
                if(!isinf(gamma_tmp)){
                    if(possiblehap[col_i] == 0){
                        gamma_i1.push_back(gamma_tmp);
                    }
                    if(possiblehap[col_i] == 1){
                        gamma_i2.push_back(gamma_tmp);
                    }
                    if(possiblehap[col_i] == 2){
                        gamma_i3.push_back(gamma_tmp);
                    }
                }
            }

            gamma_i_tmp[0] = logsum(gamma_i1);
            gamma_i_tmp[1] = logsum(gamma_i2);
            gamma_i_tmp[2] = logsum(gamma_i3);
            lognorm_vec(gamma_i_tmp);
            for(int g=0; g<3;++g){
                gamma_i[(n_m-1) * 3 + g] = gamma_i_tmp[g];
            }

            for(size_t m=n_m-1; m>0; --m){

                vector<double> gamma_i1{neg_inf};
                vector<double> gamma_i2{neg_inf};
                vector<double> gamma_i3{neg_inf};
                double gamma_tmp;
                vector<double> gamma_i_tmp(3);

                for(size_t k1=0; k1<n_h; ++k1){
                    RMatrix<double>::Row trans_prob_k = trans_prob.row(k1);
                    for(size_t k2=0; k2<n_h; ++k2){
                        trans_kk = trans_prob_k[(m-1)*n_h + k2];
                        score_k[k2] = emit[m][k2] + beta[k2] + trans_kk;
                    }
                    sum_k = logsum(score_k);
                    beta[k1] = sum_k;
                    j = p_geno[m-1];
                    col_i = j * n_h + k1;
                    gamma_tmp = sum_k + alpha[m-1][k1];
                    if(!isinf(gamma_tmp)){
                        if(possiblehap[col_i] == 0){
                            gamma_i1.push_back(gamma_tmp);
                        }
                        if(possiblehap[col_i] == 1){
                            gamma_i2.push_back(gamma_tmp);
                        }
                        if(possiblehap[col_i] == 2){
                            gamma_i3.push_back(gamma_tmp);
                        }
                    }
                }

                gamma_i_tmp[0] = logsum(gamma_i1);
                gamma_i_tmp[1] = logsum(gamma_i2);
                gamma_i_tmp[2] = logsum(gamma_i3);
                lognorm_vec(gamma_i_tmp);
                for(int g=0; g<3;++g){
                    gamma_i[(m-1) * 3 + g] = gamma_i_tmp[g];
                }
            }
        }
    }
};

// Solve the HMM for offspring
// [[Rcpp::export]]
NumericMatrix run_fb(NumericMatrix ref,
                     NumericMatrix alt,
                     NumericVector eseq_in,
                     NumericVector bias,
                     NumericMatrix mismap,
                     IntegerVector possiblehap,
                     NumericMatrix trans_prob,
                     NumericVector init_prob,
                     int & n_h,
                     int & n_o,
                     int & n_m,
                     IntegerVector p_geno,
                     IntegerVector ploidy
){

    // Initialize arrays to store output, alpha values,
    // // emittion probs, and beta values.
    NumericMatrix gamma(n_o,  n_m * 3);

    // Convert values to ones used here.
    IntegerVector dim = {n_m, n_o, n_h};
    NumericVector w1(n_m);
    NumericVector w2(n_m);
    NumericVector eseq(2);
    NumericVector mismap1 = mismap( _ , 0 );
    NumericVector mismap2 = mismap( _ , 1 );
    eseq = clone(eseq_in);
    w1 = clone(bias);
    w2 = clone(bias);
    w2 = 1 - w2;

    LogicalVector iter_sample(n_o);

    ParFB calc_fb(gamma,
                  iter_sample,
                  ref,
                  alt,
                  eseq,
                  w1,
                  w2,
                  mismap1,
                  mismap2,
                  possiblehap,
                  init_prob,
                  trans_prob,
                  dim,
                  p_geno,
                  ploidy);

    parallelFor(0, iter_sample.length(), calc_fb);

    gamma.attr("dim") = Dimension(n_o, 3, n_m);
    return gamma;
}
