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
    const RVector<double> mismap1;
    const RVector<double> mismap2;
    const RVector<int> possiblehap;
    const RVector<double> init_prob;
    const RVector<double> trans_prob;
    const RVector<int> nonzero_prob;
    const RVector<int> n_hap;
    const RVector<int> n_offspring;
    const RVector<int> n_marker;
    const RVector<int> n_nonzero_prob;
    const RVector<int> pedigree;
    const RVector<int> hap_offset;
    const RVector<int> init_offset;
    const RVector<int> trans_offset;
    const RVector<int> p_geno;
    const RVector<int> ploidy;

    ParFB(NumericMatrix gamma,
          const LogicalVector iter_sample,
          const NumericMatrix ref,
          const NumericMatrix alt,
          const NumericVector eseq,
          const NumericVector w1,
          const NumericVector mismap1,
          const NumericVector mismap2,
          const IntegerVector possiblehap,
          const NumericVector init_prob,
          const NumericVector trans_prob,
          const LogicalVector nonzero_prob,
          const IntegerVector n_hap,
          const IntegerVector n_offspring,
          const IntegerVector n_marker,
          const IntegerVector n_nonzero_prob,
          const IntegerVector pedigree,
          const IntegerVector hap_offset,
          const IntegerVector init_offset,
          const IntegerVector trans_offset,
          const IntegerVector p_geno,
          const IntegerVector ploidy)
        : gamma(gamma),
          iter_sample(iter_sample),
          ref(ref),
          alt(alt),
          eseq(eseq),
          w1(w1),
          mismap1(mismap1),
          mismap2(mismap2),
          possiblehap(possiblehap),
          init_prob(init_prob),
          trans_prob(trans_prob),
          nonzero_prob(nonzero_prob),
          n_hap(n_hap),
          n_offspring(n_offspring),
          n_marker(n_marker),
          n_nonzero_prob(n_nonzero_prob),
          pedigree(pedigree),
          hap_offset(hap_offset),
          init_offset(init_offset),
          trans_offset(trans_offset),
          p_geno(p_geno),
          ploidy(ploidy) {}

    void operator()(size_t begin, size_t end) {
        bool het = true;
        for(RVector<int>::const_iterator i = iter_sample.begin() + begin;
            i < iter_sample.begin() + end; ++i){
            int sample_i = distance(iter_sample.begin(), i);
            int pedigree_i = pedigree[sample_i];
            int n_hap_i = n_hap[pedigree_i];
            int n_nonzero_prob_i = n_nonzero_prob[pedigree_i];

            RMatrix<double>::Row gamma_i = gamma.row(sample_i);

            vector<vector<double>> alpha(n_marker[0],
                                         vector<double>(n_hap_i));
            vector<vector<double>> emit(n_marker[0],
                                        vector<double>(n_hap_i));
            vector<double> beta(n_hap_i, 0);

            vector<double> score_k(n_hap_i);
            double hap_prob;
            double trans_kk;
            double sum_k;
            int target_i;
            int target_ik2;
            int j;
            int trans_prob_target;
            int hap_offset_i = hap_offset[pedigree_i];
            int hap_offset_ij;
            int trans_offset_i = trans_offset[pedigree_i];
            int trans_offset_im;
            double neg_inf = -numeric_limits<double>::infinity();

            int nonzero_prob_offset = 0;
            double n_jap_j;
            for(int j = 0; j < pedigree_i; ++j){
                n_jap_j = n_hap[j];
                nonzero_prob_offset += n_jap_j * n_jap_j;
            }
            int nonzero_prob_offset_kk;
            bool nonzero_prob_kk;
            int nonzero_prob_k1;
            int nonzero_prob_k2;

            for(int m = 0; m < n_marker[0]; ++m){
                vector<double> prob_i = calcEmit(ref,
                                                 alt,
                                                 eseq,
                                                 w1,
                                                 mismap1,
                                                 mismap2,
                                                 m,
                                                 sample_i,
                                                 het,
                                                 ploidy[0]);
                if(m == 0){
                    j = p_geno[0];
                    hap_offset_ij = hap_offset_i + j * n_hap_i;
                    for(int k = 0; k < n_hap_i; ++k){
                        target_i = hap_offset_ij + k;
                        hap_prob = prob_i[possiblehap[target_i]];
                        alpha[m][k] = hap_prob + init_prob[k];
                    }

                } else {
                    j = p_geno[m];
                    hap_offset_ij = hap_offset_i + j * n_hap_i;
                    trans_offset_im = trans_offset_i + n_nonzero_prob_i * (m - 1);
                    for(int k2 = 0; k2 < n_hap_i; ++k2){
                        target_ik2 = n_nonzero_prob_i * k2;
                        nonzero_prob_k1 = 0;
                        for(int k1 = 0; k1 < n_hap_i; ++k1){
                            nonzero_prob_offset_kk = nonzero_prob_offset + n_hap_i * k2 + k1;
                            nonzero_prob_kk = nonzero_prob[nonzero_prob_offset_kk];
                            if(nonzero_prob_kk){
                                trans_prob_target = trans_offset_im + target_ik2 + nonzero_prob_k1;
                                trans_kk = trans_prob[trans_prob_target];
                                nonzero_prob_k1 += 1;

                            } else {
                                trans_kk = neg_inf;

                            }
                            score_k.at(k1) = alpha[m - 1][k1] + trans_kk;
                        }
                        sum_k = logsum(score_k);
                        target_i = hap_offset_ij + k2;
                        hap_prob = prob_i[possiblehap[target_i]];
                        emit[m][k2] = hap_prob;
                        alpha[m][k2] = hap_prob + sum_k;
                    }
                }
            }

            int n_levels = ploidy[0] + 1;
            double gamma_tmp;
            int target_hap;
            vector<double> gamma_i_sum(n_levels, neg_inf);
            for(int k = 0; k < n_hap_i; ++k){
                j = p_geno[n_marker[0] - 1];
                target_i = hap_offset_ij + k;
                gamma_tmp = beta[k] + alpha[n_marker[0] - 1][k];
                if(!isinf(gamma_tmp)){
                    target_hap = possiblehap[target_i];
                    logsum2(gamma_i_sum[target_hap], gamma_tmp);
                }
            }

            lognorm_vec(gamma_i_sum);
            for(int g = 0; g < n_levels; ++g){
                gamma_i[(n_marker[0] - 1) * n_levels + g] = gamma_i_sum[g];
            }

            int trans_offset_imk1;
            for(int m = n_marker[0] - 1; m > 0; --m){
                double gamma_tmp;
                int target_hap;
                vector<double> gamma_i_sum(n_levels, neg_inf);
                trans_offset_im = trans_offset_i + n_nonzero_prob_i * (m - 1);

                for(int k1 = 0; k1 < n_hap_i; ++k1){
                    trans_offset_imk1 = trans_offset_im + k1;
                    nonzero_prob_k2 = 0;
                    for(int k2 = 0; k2 < n_hap_i; ++k2){

                        nonzero_prob_offset_kk = nonzero_prob_offset + n_hap_i * k2 + k1;
                        nonzero_prob_kk = nonzero_prob[nonzero_prob_offset_kk];
                        if(nonzero_prob_kk){
                            trans_prob_target = trans_offset_imk1 + n_nonzero_prob_i * nonzero_prob_k2;
                            trans_kk = trans_prob[trans_prob_target];
                            nonzero_prob_k1 += 1;

                        } else {
                            trans_kk = neg_inf;

                        }

                        score_k[k2] = emit[m][k2] + beta[k2] + trans_kk;
                    }
                    sum_k = logsum(score_k);
                    beta[k1] = sum_k;
                    j = p_geno[m - 1];
                    target_i = hap_offset_i + j * n_hap_i + k1;
                    gamma_tmp = sum_k + alpha[m - 1][k1];
                    if(!isinf(gamma_tmp)){
                        target_hap = possiblehap[target_i];
                        logsum2(gamma_i_sum[target_hap], gamma_tmp);
                    }
                }

                lognorm_vec(gamma_i_sum);
                for(int g = 0; g < n_levels; ++g){
                    gamma_i[(m-1) * n_levels + g] = gamma_i_sum[g];
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
                     NumericVector trans_prob,
                     NumericVector init_prob,
                     LogicalVector nonzero_prob,
                     IntegerVector n_pgeno,
                     IntegerVector n_hap,
                     IntegerVector n_offspring,
                     IntegerVector n_marker,
                     IntegerVector n_nonzero_prob,
                     IntegerVector pedigree,
                     IntegerVector p_geno,
                     IntegerVector ploidy
){

    // Initialize arrays to store output, alpha values,
    // // emittion probs, and beta values.
    int n_levels = ploidy[0] + 1;
    NumericMatrix gamma(n_offspring[0],  n_levels * n_marker[0]);

    // Convert values to ones used here.
    NumericVector w1(n_marker[0]);
    NumericVector eseq(2);
    NumericVector mismap1 = mismap( _ , 0 );
    NumericVector mismap2 = mismap( _ , 1 );
    eseq = clone(eseq_in);
    w1 = clone(bias);

    // Calculate offsets to access pedigree dependent parameters
    IntegerVector hap_offset(n_hap.size());
    IntegerVector init_offset(n_hap.size());
    IntegerVector trans_offset(n_hap.size());
    for(int i = 1; i < n_hap.size(); ++i){
        hap_offset[i] = hap_offset[i - 1] + n_pgeno[0] * n_hap[i - 1];
        init_offset[i] = init_offset[i - 1] + n_hap[i - 1];
        trans_offset[i] = trans_offset[i - 1] +
            n_nonzero_prob[i - 1] * n_nonzero_prob[i - 1] * (n_marker[0] - 1);
    }

    LogicalVector iter_sample(n_offspring[0]);

    ParFB calc_fb(gamma,
                  iter_sample,
                  ref,
                  alt,
                  eseq,
                  w1,
                  mismap1,
                  mismap2,
                  possiblehap,
                  init_prob,
                  trans_prob,
                  nonzero_prob,
                  n_hap,
                  n_offspring,
                  n_marker,
                  n_nonzero_prob,
                  pedigree,
                  hap_offset,
                  init_offset,
                  trans_offset,
                  p_geno,
                  ploidy);

    parallelFor(0, iter_sample.length(), calc_fb);

    gamma.attr("dim") = Dimension(n_offspring[0], n_levels, n_marker[0]);
    return gamma;
}
