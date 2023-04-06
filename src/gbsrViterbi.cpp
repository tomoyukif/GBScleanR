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
// Functions for the Viterbi algorithm
// Run Viterbi algorithm for founder genotype and offspring genotype separately

// Initialize Viterbi scores
struct ParInitVit : public Worker {

    RMatrix<double> vit_score;
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
    const RVector<int> n_hap;
    const RVector<int> pedigree;
    const RVector<int> hap_offset;
    const RVector<int> init_offset;
    const RVector<int> valid_p_indices;
    const RVector<int> m;
    const RVector<int> ploidy;

    ParInitVit(NumericMatrix vit_score,
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
               const IntegerVector n_hap,
               const IntegerVector pedigree,
               const IntegerVector hap_offset,
               const IntegerVector init_offset,
               const IntegerVector valid_p_indices,
               const IntegerVector m,
               const IntegerVector ploidy)
        : vit_score(vit_score),
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
          n_hap(n_hap),
          pedigree(pedigree),
          hap_offset(hap_offset),
          init_offset(init_offset),
          valid_p_indices(valid_p_indices),
          m(m),
          ploidy(ploidy) {}

    void operator()(size_t begin, size_t end) {
        bool het = true;
        for(RVector<int>::const_iterator i = iter_sample.begin() + begin;
            i < iter_sample.begin() + end; ++i){
            int sample_i = distance(iter_sample.begin(), i);
            int pedigree_i = pedigree[sample_i];
            RMatrix<double>::Row vit_i = vit_score.row(sample_i);

            // Calculate genotype probabilies
            vector<double> prob_i = calcEmit(ref,
                                             alt,
                                             eseq,
                                             w1,
                                             w2,
                                             mismap1,
                                             mismap2,
                                             m[0],
                                             sample_i,
                                             het);
            double hap_prob = 0.0;
            int n_hap_i = n_hap[pedigree_i];
            int target_i;
            int valid_j;
            int target_k;

            for(int j = 0; j < (int)valid_p_indices.size(); ++j){
                for(int k = 0; k < n_hap_i; ++k){
                    valid_j = valid_p_indices[j];
                    target_i = hap_offset[pedigree_i] + valid_j * n_hap_i + k;
                    hap_prob = prob_i[possiblehap[target_i]];
                    target_i = j * n_hap_i + k;
                    target_k = init_offset[pedigree_i] + k;
                    vit_i[target_i] = hap_prob + init_prob[target_k];
                }
            }
        }
    }
};

// Calculate Viterbi scores
struct ParCalcVitFounder : public Worker {

    RMatrix<double> in_score;
    const RVector<int> iter_sample;
    const RVector<double> trans_prob;
    const RVector<int> m;
    const RVector<int> n_hap;
    const RVector<int> pedigree;
    const RVector<int> trans_offset;
    const RVector<int> valid_p_indices;

    ParCalcVitFounder(NumericMatrix in_score,
                      const LogicalVector iter_sample,
                      const NumericVector trans_prob,
                      const IntegerVector m,
                      const IntegerVector n_hap,
                      const IntegerVector pedigree,
                      const IntegerVector trans_offset,
                      const IntegerVector valid_p_indices)
        : in_score(in_score),
          iter_sample(iter_sample),
          trans_prob(trans_prob),
          m(m),
          n_hap(n_hap),
          pedigree(pedigree),
          trans_offset(trans_offset),
          valid_p_indices(valid_p_indices) {}

    void operator()(size_t begin, size_t end) {

        for(RVector<int>::const_iterator i = iter_sample.begin() + begin;
            i < iter_sample.begin() + end; ++i){
            int sample_i = distance(iter_sample.begin(), i);
            int pedigree_i = pedigree[sample_i];
            RMatrix<double>::Row in_i = in_score.row(sample_i);

            int n_hap_i = n_hap[pedigree_i];
            int target_i;
            double trans_kk;
            int max_i;
            vector<double> score_jkk(n_hap_i);
            vector<double> max_scores_jk(n_hap_i);
            int trans_prob_target;

            for(int j = 0; j < (int)valid_p_indices.size(); ++j){
                for(int k2 = 0; k2 < n_hap_i; ++k2){
                    for(int k1 = 0; k1 < n_hap_i; ++k1){
                        trans_prob_target = trans_offset[pedigree_i] +
                            n_hap_i * n_hap_i * (m[0] - 1) +
                            n_hap_i * k2 + k1;
                        trans_kk = trans_prob[trans_prob_target];
                        target_i = j * n_hap_i + k1;
                        score_jkk[k1] = in_i[target_i] + trans_kk;
                    }
                    max_i = get_max_int(score_jkk);
                    max_scores_jk[k2] = score_jkk[max_i];
                }
                for(int k2 = 0; k2 < n_hap_i; ++k2){
                    target_i = j * n_hap_i + k2;
                    in_i[target_i] = max_scores_jk[k2];
                }
            }
        }
    }
};

// Calculate Founder Viterbi scores
struct ParCalcPathFounder : public Worker {

    RMatrix<int> f_path;
    RMatrix<double> vit_score;
    const RMatrix<double> in_score;
    const RVector<int> iter_p_pat;
    const RMatrix<double> ref;
    const RMatrix<double> alt;
    const RVector<double> eseq;
    const RVector<double> w1;
    const RVector<double> w2;
    const RVector<double> mismap1;
    const RVector<double> mismap2;
    const RVector<int> possiblehap;
    const RVector<int> n_offspring;
    const RVector<int> n_pgeno;
    const RVector<int> n_hap;
    const RVector<int> pedigree;
    const RVector<int> hap_offset;
    const RVector<double> p_emit1;
    const RVector<double> p_emit2;
    const RVector<int> valid_p_indices1;
    const RVector<int> valid_p_indices2;
    const RVector<int> m;
    const RVector<int> ploidy;

    ParCalcPathFounder(IntegerMatrix f_path,
                       NumericMatrix vit_score,
                       const NumericMatrix in_score,
                       const LogicalVector iter_p_pat,
                       const NumericMatrix ref,
                       const NumericMatrix alt,
                       const NumericVector eseq,
                       const NumericVector w1,
                       const NumericVector w2,
                       const NumericVector mismap1,
                       const NumericVector mismap2,
                       const IntegerVector possiblehap,
                       const IntegerVector n_offspring,
                       const IntegerVector n_pgeno,
                       const IntegerVector n_hap,
                       const IntegerVector pedigree,
                       const IntegerVector hap_offset,
                       const NumericVector p_emit1,
                       const NumericVector p_emit2,
                       const IntegerVector valid_p_indices1,
                       const IntegerVector valid_p_indices2,
                       const IntegerVector m,
                       const IntegerVector ploidy)
        : f_path(f_path),
          vit_score(vit_score),
          in_score(in_score),
          iter_p_pat(iter_p_pat),
          ref(ref),
          alt(alt),
          eseq(eseq),
          w1(w1),
          w2(w2),
          mismap1(mismap1),
          mismap2(mismap2),
          possiblehap(possiblehap),
          n_offspring(n_offspring),
          n_pgeno(n_pgeno),
          n_hap(n_hap),
          pedigree(pedigree),
          hap_offset(hap_offset),
          p_emit1(p_emit1),
          p_emit2(p_emit2),
          valid_p_indices1(valid_p_indices1),
          valid_p_indices2(valid_p_indices2),
          m(m),
          ploidy(ploidy) {}

    void operator()(size_t begin, size_t end) {
        bool het = true;
        for(RVector<int>::const_iterator i = iter_p_pat.begin() + begin;
            i < iter_p_pat.begin() + end; ++i){
            RMatrix<int>::Row f_path_m = f_path.row(m[0]);

            int j2 = distance(iter_p_pat.begin(), i);
            double neg_inf = -numeric_limits<double>::infinity();
            double hap_prob;
            int target_i;
            vector<double> score_j(n_pgeno[0]);
            vector<double> score_ijk(n_hap[0]);
            int max_j;

            if(isinf(p_emit2[j2])){
                f_path_m[j2] = -1;

            } else {
                for(int sample_i = 0; sample_i < n_offspring[0]; ++sample_i){
                    RMatrix<double>::Row in_i = in_score.row(sample_i);
                    int pedigree_i = pedigree[sample_i];
                    int n_hap_i = n_hap[pedigree_i];

                    // Calculate genotype probabilies
                    vector<double> prob_i = calcEmit(ref,
                                                     alt,
                                                     eseq,
                                                     w1,
                                                     w2,
                                                     mismap1,
                                                     mismap2,
                                                     m[0],
                                                     sample_i,
                                                     het);

                    for(int j1 = 0; j1 < n_pgeno[0]; ++j1){
                        if(isinf(p_emit1[j1])){
                            score_j[j1] = neg_inf;
                        } else {

                            int valid_j1 = 0;
                            for(size_t j = 0; j < valid_p_indices1.size(); ++j){
                                valid_j1 = valid_p_indices1[j];
                                if(j1 == valid_j1){
                                    valid_j1 = j;
                                    break;
                                }
                            }
                            for(int k = 0; k < n_hap_i; ++k){
                                target_i = hap_offset[pedigree_i] + j2 * n_hap_i + k;
                                hap_prob = prob_i[possiblehap[target_i]];
                                target_i = valid_j1 * n_hap_i + k;
                                score_ijk[k] = in_i[target_i] + hap_prob;
                            }
                            score_j[j1] += logsum(score_ijk);
                        }
                    }
                }
                for(int j1 = 0; j1 < n_pgeno[0]; ++j1){
                    score_j[j1] = p_emit1[j1] + score_j[j1];
                }

                max_j = get_max_int(score_j);
                f_path_m[j2] = max_j;

                int max_valid_index = 0;
                int valid_j1 = 0;
                for(int j = 0; j < (int)valid_p_indices1.size(); ++j){
                    valid_j1 = valid_p_indices1[j];
                    if(max_j == valid_j1){
                        max_valid_index = j;
                    }
                }

                int j2_valid_index = 0;
                int valid_j2;
                for(int j = 0; j < (int)valid_p_indices2.size(); ++j){
                    valid_j2 = valid_p_indices2[j];
                    if(j2 == valid_j2){
                        j2_valid_index = j;
                    }
                }

                int out_i;

                for(int sample_i = 0; sample_i < n_offspring[0]; ++sample_i){
                    RMatrix<double>::Row vit_i = vit_score.row(sample_i);
                    RMatrix<double>::Row in_i = in_score.row(sample_i);
                    int pedigree_i = pedigree[sample_i];
                    int n_hap_i = n_hap[pedigree_i];

                    // Calculate genotype probabilies
                    vector<double> prob_i = calcEmit(ref,
                                                     alt,
                                                     eseq,
                                                     w1,
                                                     w2,
                                                     mismap1,
                                                     mismap2,
                                                     m[0],
                                                     sample_i,
                                                     het);

                    for(int k = 0; k < n_hap_i; ++k){
                        target_i = hap_offset[pedigree_i] + j2 * n_hap_i + k;
                        hap_prob = prob_i[possiblehap[target_i]];
                        target_i = max_valid_index * n_hap_i + k;
                        out_i = j2_valid_index * n_hap_i + k;
                        vit_i[out_i] = in_i[target_i] + hap_prob;
                    }
                }
            }
        }
    }
};

// Calculate alpha values in the forward algorithm
void last_vit_founder(IntegerVector f_seq,
                      NumericMatrix vit_score,
                      NumericVector p_emit,
                      IntegerVector n_offspring,
                      IntegerVector n_pgeno,
                      IntegerVector n_hap,
                      IntegerVector pedigree,
                      IntegerVector n_marker,
                      IntegerVector valid_p_indices){
    size_t m = n_marker[0]-1;
    size_t vit_i;
    double neg_inf = -numeric_limits<double>::infinity();

    // Calculate emit, alpha, and gamma
    IntegerMatrix tmp_argmax(n_offspring[0], valid_p_indices.size());
    vector<double> score_ijk(n_hap[0]);
    vector<double> score_j(n_pgeno[0]);

    for(int i = 0; i < n_offspring[0]; ++i){
        for(int j1 = 0; j1 < n_pgeno[0]; ++j1){
            int valid_j1 = 0;
            for(R_xlen_t j = 0; j < valid_p_indices.size(); ++j){
                valid_j1 = valid_p_indices[j];
                if(j1 == valid_j1){
                    valid_j1 = j;
                    break;
                }
            }

            if(isinf(p_emit[j1])){
                score_j[j1] = neg_inf;

            } else {
                for(int k = 0; k < n_hap[0]; ++k){
                    vit_i = valid_j1 * n_hap[0] + k;
                    score_ijk[k] = vit_score(i, vit_i);
                }
                score_j[j1] += logsum(score_ijk);
            }
        }
    }
    for(int j = 0; j < n_pgeno[0]; ++j){
        score_j[j] = p_emit[j] + score_j[j];
    }
    size_t max_j = get_max_int(score_j);
    f_seq[m] = max_j;
}

// Find the best state sequences backwardly
void backtrack(IntegerMatrix f_path,
               IntegerVector f_seq,
               IntegerVector n_marker){
    size_t f_prev;

    for(int m = n_marker[0] - 1; m > 0; --m){
        if(m % 10 == 9){
            Rcpp::Rcout << "\r" <<
                "Backtracking best genotype sequences at marker#: " <<
                    m+1 << string(70, ' ');
        }
        f_prev = f_seq[m];
        f_seq[m-1] = f_path(m , f_prev);
    }
    Rcpp::Rcout << "\r" <<
        "Backtracking best genotype sequences: Done!" <<
            string(70, ' ');
}

// Functions for the Viterbi algorithm For OFFSPRING
struct ParVitOffspring : public Worker {

    RMatrix<int> o_seq;
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
    const RVector<double> trans_prob;
    const RVector<int> n_marker;
    const RVector<int> n_hap;
    const RVector<int> pedigree;
    const RVector<int> hap_offset;
    const RVector<int> init_offset;
    const RVector<int> trans_offset;
    const RVector<int> f_seq;
    const RVector<int> ploidy;

    ParVitOffspring(IntegerMatrix o_seq,
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
                    const NumericVector trans_prob,
                    const IntegerVector n_marker,
                    const IntegerVector n_hap,
                    const IntegerVector pedigree,
                    const IntegerVector hap_offset,
                    const IntegerVector init_offset,
                    const IntegerVector trans_offset,
                    const IntegerVector f_seq,
                    const IntegerVector ploidy)
        : o_seq(o_seq),
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
          n_marker(n_marker),
          n_hap(n_hap),
          pedigree(pedigree),
          hap_offset(hap_offset),
          init_offset(init_offset),
          trans_offset(trans_offset),
          f_seq(f_seq),
          ploidy(ploidy) {}

    void operator()(size_t begin, size_t end) {
        bool het = true;
        for(RVector<int>::const_iterator i = iter_sample.begin() + begin;
            i < iter_sample.begin() + end; ++i){
            int sample_i = distance(iter_sample.begin(), i);
            int pedigree_i = pedigree[sample_i];
            int n_hap_i = n_hap[pedigree_i];

            RMatrix<int>::Column o_seq_i = o_seq.column(sample_i);
            vector<vector<unsigned short>> o_path(n_marker[0],
                                                  vector<unsigned short>(n_hap[0]));
            int target_i;
            double hap_prob = 0.0;
            vector<double> vit(n_hap[0]);
            double trans_kk;
            int max_i;
            vector<double> score_jkk(n_hap[0]);
            vector<double> max_scores_jk(n_hap[0]);
            int o_prev;

            // Viterbi path
            for(int m = 0; m < n_marker[0]; ++m){
                // Calculate genotype probabilies
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

                int f_geno = f_seq[m];
                int trans_prob_target;

                if(m == 0){
                    for(int k = 0; k < n_hap_i; ++k){
                        target_i = hap_offset[pedigree_i] + f_geno * n_hap_i + k;
                        hap_prob = prob_i[possiblehap[target_i]];
                        vit[k] = hap_prob + init_prob[k];
                    }
                } else {

                    for(int k2 = 0; k2 < n_hap_i; ++k2){
                        for(int k1 = 0; k1 < n_hap_i; ++k1){
                            trans_prob_target = trans_offset[pedigree_i] +
                                n_hap_i * n_hap_i * (m - 1) +
                                n_hap_i * k2 + k1;
                            trans_kk = trans_prob[trans_prob_target];
                            score_jkk[k1] = vit[k1] + trans_kk;
                        }
                        max_i = get_max_int(score_jkk);
                        o_path[m][k2] = max_i;
                        max_scores_jk[k2] = score_jkk[max_i];
                    }
                    for(int k2 = 0; k2 < n_hap_i; ++k2){
                        target_i = hap_offset[pedigree_i] + f_geno * n_hap_i + k2;
                        hap_prob = prob_i[possiblehap[target_i]];
                        vit[k2] = max_scores_jk[k2] + hap_prob;
                    }
                }
                if(m == n_marker[0] - 1){
                    o_seq_i[m] = get_max_int(vit);
                }
            }

            // Backtracking
            for(int m = n_marker[0] - 1; m > 0; --m){
                o_prev = o_seq_i[m];
                o_seq_i[m-1] = o_path[m][o_prev];
            }
        }
    }
};


// Solve the HMM
// [[Rcpp::export]]
List run_viterbi(NumericMatrix p_ref,
                 NumericMatrix p_alt,
                 NumericMatrix ref,
                 NumericMatrix alt,
                 NumericVector eseq_in,
                 NumericVector bias,
                 NumericMatrix mismap,
                 NumericVector trans_prob,
                 NumericVector init_prob,
                 IntegerVector n_pgeno,
                 IntegerVector n_hap,
                 IntegerVector n_offspring,
                 IntegerVector n_founder,
                 IntegerVector n_marker,
                 LogicalVector het,
                 IntegerVector pedigree,
                 IntegerVector possiblehap,
                 IntegerVector possiblegeno,
                 IntegerVector p_geno_fix,
                 IntegerVector ploidy
){
    // Initialize arrays to store output, alpha values,
    // emittion probs, and beta values.
    double neg_inf = -numeric_limits<double>::infinity();
    IntegerVector f_seq(n_marker[0]);
    IntegerMatrix f_path(n_marker[0], n_pgeno[0]);
    IntegerMatrix o_seq(n_marker[0], n_offspring[0]);

    // Convert values to ones used here.
    NumericVector w1(n_marker[0]);
    NumericVector w2(n_marker[0]);
    NumericVector eseq(2);
    NumericVector mismap1 = mismap( _ , 0 );
    NumericVector mismap2 = mismap( _ , 1 );
    eseq = clone(eseq_in);
    w1 = clone(bias);
    w2 = clone(bias);
    w2 = 1 - w2;


    // Calculate offsets to access pedigree dependent parameters
    IntegerVector hap_offset(n_hap.size());
    IntegerVector init_offset(n_hap.size());
    IntegerVector trans_offset(n_hap.size());
    for(int i = 1; i < n_hap.size(); ++i){
        hap_offset[i] = hap_offset[i - 1] + n_pgeno[0] * n_hap[i - 1];
        init_offset[i] = init_offset[i - 1] + n_hap[i - 1];
        trans_offset[i] = trans_offset[i - 1] +
            n_hap[i - 1] * n_hap[i - 1] * (n_marker[0] - 1);
    }

    LogicalVector iter_sample(n_offspring[0]);
    LogicalVector iter_p_pat(n_pgeno[0]);
    NumericVector p_emit1;
    NumericVector p_emit2;

    // Initialize Viterbi trellis.
    int m = 0;
    p_emit1 = calcPemit(p_ref,
                        p_alt,
                        eseq,
                        w1,
                        w2,
                        mismap1,
                        mismap2,
                        possiblegeno,
                        m,
                        n_founder,
                        n_pgeno,
                        het,
                        ploidy);

    Rcpp::Rcout << "\r" <<
        "Founder genotype probability calculation ..." <<
            string(70, ' ');
    int int_fix_p = p_geno_fix[0];
    if(int_fix_p >= 0){
        R_xlen_t fix_p = p_geno_fix[0];
        for(R_xlen_t p = 0; p < p_emit1.size(); ++p){
            if(p == fix_p){
                p_emit1[p] = 0;
            } else {
                p_emit1[p] = neg_inf;
            }
        }
    }

    IntegerVector valid_p_indices1;
    for(R_xlen_t j = 0; j < p_emit1.size(); ++j){
        if(!isinf(p_emit1[j])){
            valid_p_indices1.push_back(j);
        }
    }
    int valid_size = valid_p_indices1.size();
    int max_hap = n_hap[0];
    for(int i = 1; i < (int)n_hap.size(); ++i){
        if(max_hap < n_hap[i]){
            max_hap = n_hap[i];
        }
    }
    NumericMatrix vit_score(n_offspring[0], valid_size * max_hap);
    IntegerVector in_m = {m};
    ParInitVit init_vit(vit_score,
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
                        n_hap,
                        pedigree,
                        hap_offset,
                        init_offset,
                        valid_p_indices1,
                        in_m,
                        ploidy);
    parallelFor(0, iter_sample.length(), init_vit);
    NumericMatrix in_score;
    in_score = clone(vit_score);

    for(int m = 1; m < n_marker[0]; ++m){
        if(m % 10 == 9){
            Rcpp::Rcout << "\r" <<
                "Founder genotype probability calculation at marker#: " <<
                    m + 1 << string(70, ' ');
        }
        IntegerVector in_m = {m};
        ParCalcVitFounder calc_vit(in_score,
                                   iter_sample,
                                   trans_prob,
                                   in_m,
                                   n_hap,
                                   pedigree,
                                   trans_offset,
                                   valid_p_indices1);
        parallelFor(0, iter_sample.length(), calc_vit);

        p_emit2 = calcPemit(p_ref,
                            p_alt,
                            eseq,
                            w1,
                            w2,
                            mismap1,
                            mismap2,
                            possiblegeno,
                            m,
                            n_founder,
                            n_pgeno,
                            het,
                            ploidy);

        R_xlen_t fix_p_len = p_geno_fix.size();
        if(fix_p_len > m){
            R_xlen_t fix_p = p_geno_fix[m];
            for(R_xlen_t p = 0; p < p_emit2.size(); ++p){
                if(p == fix_p){
                    p_emit2[p] = 0;
                } else {
                    p_emit2[p] = neg_inf;
                }
            }
        };

        IntegerVector valid_p_indices2;
        for(R_xlen_t j = 0; j < p_emit2.size(); ++j){
            if(!isinf(p_emit2[j])){
                valid_p_indices2.push_back(j);
            }
        }

        size_t valid_size = valid_p_indices2.size();
        NumericMatrix vit_score(n_offspring[0], valid_size * max_hap);
        ParCalcPathFounder calc_path(f_path,
                                     vit_score,
                                     in_score,
                                     iter_p_pat,
                                     ref,
                                     alt,
                                     eseq,
                                     w1,
                                     w2,
                                     mismap1,
                                     mismap2,
                                     possiblehap,
                                     n_offspring,
                                     n_pgeno,
                                     n_hap,
                                     pedigree,
                                     hap_offset,
                                     p_emit1,
                                     p_emit2,
                                     valid_p_indices1,
                                     valid_p_indices2,
                                     in_m,
                                     ploidy);

        parallelFor(0, iter_p_pat.length(), calc_path);
        p_emit1 = clone(p_emit2);
        valid_p_indices1 = clone(valid_p_indices2);
        in_score = clone(vit_score);

        if(m == n_marker[0] - 1){

            last_vit_founder(f_seq,
                             in_score,
                             p_emit1,
                             n_offspring,
                             n_pgeno,
                             n_hap,
                             pedigree,
                             n_marker,
                             valid_p_indices1);
        }
    }

    backtrack(f_path,
              f_seq,
              n_marker);

    Rcpp::Rcout << "\r" <<
        "Offspring genotype probability calculation ..." <<
            string(70, ' ');

    ParVitOffspring vit_offspring(o_seq,
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
                                  n_marker,
                                  n_hap,
                                  pedigree,
                                  hap_offset,
                                  init_offset,
                                  trans_offset,
                                  f_seq,
                                  ploidy);
    parallelFor(0, iter_sample.length(), vit_offspring);

    Rcpp::Rcout << "\r" << string(70, ' ');
    List out_list = List::create(_["p_geno"] = f_seq, _["best_seq"] = o_seq);
    return out_list;
}
