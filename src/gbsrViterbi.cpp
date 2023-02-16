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
    const RVector<int> dim;
    const RVector<int> valid_p_indices;
    const RVector<int> vec_m;
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
               const IntegerVector dim,
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
          dim(dim),
          valid_p_indices(valid_p_indices),
          vec_m(m),
          ploidy(ploidy) {}

    void operator()(size_t begin, size_t end) {
        bool het = true;
        for(RVector<int>::const_iterator i=iter_sample.begin() + begin;
            i<iter_sample.begin() + end; ++i){
            size_t sample_i = distance(iter_sample.begin(), i);
            size_t m = vec_m[0];
            RMatrix<double>::Row vit_i = vit_score.row(sample_i);

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
            double hap_prob = 0.0;
            size_t col_i;
            size_t n_h = dim[4];
            size_t valid_j;

            for(size_t j=0; j<valid_p_indices.size(); ++j){
                for(size_t k=0; k<n_h; ++k){
                    valid_j = valid_p_indices[j];
                    col_i = valid_j * n_h + k;
                    hap_prob = prob_i[possiblehap[col_i]];
                    col_i = j * n_h + k;
                    vit_i[col_i] = hap_prob + init_prob[k];
                }
            }
        }
    }
};

// Calculate Viterbi scores
struct ParCalcVitFounder : public Worker {

    RMatrix<double> in_score;
    const RVector<int> iter_sample;
    const RMatrix<double> trans_prob_m;
    const RVector<int> dim;
    const RVector<int> valid_p_indices;

    ParCalcVitFounder(NumericMatrix in_score,
                      const LogicalVector iter_sample,
                      const NumericMatrix trans_prob_m,
                      const IntegerVector dim,
                      const IntegerVector valid_p_indices)
        : in_score(in_score),
          iter_sample(iter_sample),
          trans_prob_m(trans_prob_m),
          dim(dim),
          valid_p_indices(valid_p_indices) {}

    void operator()(size_t begin, size_t end) {

        for(RVector<int>::const_iterator i=iter_sample.begin() + begin;
            i<iter_sample.begin() + end; ++i){
            size_t sample_i = distance(iter_sample.begin(), i);
            RMatrix<double>::Row in_i = in_score.row(sample_i);

            size_t col_i;
            size_t n_h = dim[4];
            double trans_kk;
            size_t max_i;
            vector<double> score_jkk(n_h);
            vector<double> max_scores_jk(n_h);

            for(size_t j=0; j<valid_p_indices.size(); ++j){
                for(size_t k2=0; k2<n_h; ++k2){
                    RMatrix<double>::Column trans_prob_k = trans_prob_m.column(k2);
                    for(size_t k1=0; k1<n_h; ++k1){
                        col_i = j * n_h + k1;
                        trans_kk = trans_prob_k[k1];
                        score_jkk[k1] = in_i[col_i] + trans_kk;
                    }
                    max_i = get_max_int(score_jkk);
                    max_scores_jk[k2] = score_jkk[max_i];
                }
                for(size_t k2=0; k2<n_h; ++k2){
                    col_i = j * n_h + k2;
                    in_i[col_i] = max_scores_jk[k2];
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
    const RVector<int> dim;
    const RVector<double> p_emit1;
    const RVector<double> p_emit2;
    const RVector<int> valid_p_indices1;
    const RVector<int> valid_p_indices2;
    const RVector<int> vec_m;
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
                       const IntegerVector dim,
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
          dim(dim),
          p_emit1(p_emit1),
          p_emit2(p_emit2),
          valid_p_indices1(valid_p_indices1),
          valid_p_indices2(valid_p_indices2),
          vec_m(m),
          ploidy(ploidy) {}

    void operator()(size_t begin, size_t end) {
        bool het = true;
        for(RVector<int>::const_iterator i=iter_p_pat.begin() + begin;
            i<iter_p_pat.begin() + end; ++i){
            size_t m = vec_m[0];
            RMatrix<int>::Row f_path_m = f_path.row(m);

            size_t j2 = distance(iter_p_pat.begin(), i);
            double neg_inf = -numeric_limits<double>::infinity();
            double hap_prob;
            size_t n_o = dim[2];
            size_t n_p = dim[3];
            size_t n_h = dim[4];
            size_t col_i;
            vector<double> score_j(n_p);
            vector<double> score_ijk(n_h);
            size_t max_j;
            size_t col_in;


            if(isinf(p_emit2[j2])){
                f_path_m[j2] = -1;

            } else {
                for(size_t sample_i=0; sample_i<n_o; ++sample_i){
                    RMatrix<double>::Row in_i = in_score.row(sample_i);

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

                    for(size_t j1=0; j1<n_p; ++j1){
                        if(isinf(p_emit1[j1])){
                            score_j[j1] = neg_inf;
                        } else {

                            size_t valid_j1 = 0;
                            for(size_t j=0;j<valid_p_indices1.size();++j){
                                valid_j1 = valid_p_indices1[j];
                                if(j1 == valid_j1){
                                    valid_j1 = j;
                                    break;
                                }
                            }
                            for(size_t k=0; k<n_h; ++k){
                                hap_prob = prob_i[possiblehap[j2 * n_h + k]];
                                col_in = valid_j1 * n_h + k;
                                score_ijk[k] = in_i[col_in] + hap_prob;
                            }
                            score_j[j1] += logsum(score_ijk);
                        }
                    }
                }
                for(size_t j1=0; j1<n_p; ++j1){
                    score_j[j1] = p_emit1[j1] + score_j[j1];
                }

                max_j = get_max_int(score_j);
                f_path_m[j2] = max_j;

                size_t max_valid_index = 0;
                size_t valid_j1 = 0;
                for(size_t j=0;j<valid_p_indices1.size();++j){
                    valid_j1 = valid_p_indices1[j];
                    if(max_j == valid_j1){
                        max_valid_index = j;
                    }
                }

                size_t j2_valid_index = 0;
                size_t valid_j2;
                for(size_t j=0;j<valid_p_indices2.size();++j){
                    valid_j2 = valid_p_indices2[j];
                    if(j2 == valid_j2){
                        j2_valid_index = j;
                    }
                }

                size_t out_i;

                for(size_t sample_i=0; sample_i<n_o; ++sample_i){
                    RMatrix<double>::Row vit_i = vit_score.row(sample_i);
                    RMatrix<double>::Row in_i = in_score.row(sample_i);

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

                    for(size_t k=0; k<n_h; ++k){
                        col_i = j2 * n_h + k;
                        hap_prob = prob_i[possiblehap[col_i]];
                        col_in = max_valid_index * n_h + k;
                        out_i = j2_valid_index * n_h + k;
                        vit_i[out_i] = in_i[col_in] + hap_prob;
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
                      IntegerVector dim,
                      IntegerVector valid_p_indices){
    size_t n_o = dim[2];
    size_t n_p = dim[3];
    size_t n_h = dim[4];
    size_t m = dim[0]-1;
    size_t vit_i;
    double neg_inf = -numeric_limits<double>::infinity();

    // Calculate emit, alpha, and gamma
    IntegerMatrix tmp_argmax(n_o, valid_p_indices.size());
    vector<double> score_ijk(n_h);
    vector<double> score_j(n_p);

    for(size_t i=0; i<n_o; ++i){
        for(size_t j1=0; j1<n_p; ++j1){
            size_t valid_j1 = 0;
            for(R_xlen_t j=0;j<valid_p_indices.size();++j){
                valid_j1 = valid_p_indices[j];
                if(j1 == valid_j1){
                    valid_j1 = j;
                    break;
                }
            }

            if(isinf(p_emit[j1])){
                score_j[j1] = neg_inf;

            } else {
                for(size_t k=0; k<n_h; ++k){
                    vit_i = valid_j1 * n_h + k;
                    score_ijk[k] = vit_score(i, vit_i);
                }
                score_j[j1] += logsum(score_ijk);
            }
        }
    }
    for(size_t j=0; j<n_p; ++j){
        score_j[j] = p_emit[j] + score_j[j];
    }
    size_t max_j = get_max_int(score_j);
    f_seq[m] = max_j;
}

// Find the best state sequences backwardly
void backtrack(IntegerMatrix f_path,
               IntegerVector f_seq,
               IntegerVector dim){
    size_t f_prev;
    size_t n_m = dim[0];

    for(size_t m=n_m-1; m>0; --m){
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
    const RMatrix<double> trans_prob;
    const RVector<int> dim;
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
                    const NumericMatrix trans_prob,
                    const IntegerVector dim,
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
          dim(dim),
          f_seq(f_seq),
          ploidy(ploidy) {}

    void operator()(size_t begin, size_t end) {
        bool het = true;
        for(RVector<int>::const_iterator i=iter_sample.begin() + begin;
            i<iter_sample.begin() + end; ++i){
            size_t sample_i = distance(iter_sample.begin(), i);
            RMatrix<int>::Column o_seq_i = o_seq.column(sample_i);
            size_t n_m = dim[0];
            size_t n_h = dim[4];
            vector<vector<unsigned short>> o_path(n_m,
                                                  vector<unsigned short>(n_h));
            size_t col_i;
            double hap_prob = 0.0;
            vector<double> vit(n_h);
            double trans_kk;
            size_t max_i;
            vector<double> score_jkk(n_h);
            vector<double> max_scores_jk(n_h);
            size_t o_prev;

            // Viterbi path
            for(size_t m=0; m<n_m; ++m){
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

                size_t f_geno = f_seq[m];
                size_t trans_prob_col;

                if(m == 0){
                    for(size_t k=0; k<n_h; ++k){
                        col_i = f_geno * n_h + k;
                        hap_prob = prob_i[possiblehap[col_i]];
                        vit[k] = hap_prob + init_prob[k];
                    }
                } else {

                    for(size_t k2=0; k2<n_h; ++k2){
                        trans_prob_col = (m-1)*n_h + k2;
                        RMatrix<double>::Column trans_prob_k =
                            trans_prob.column(trans_prob_col);
                        for(size_t k1=0; k1<n_h; ++k1){
                            trans_kk = trans_prob_k[k1];
                            score_jkk[k1] = vit[k1] + trans_kk;
                        }
                        max_i = get_max_int(score_jkk);
                        o_path[m][k2] = max_i;
                        max_scores_jk[k2] = score_jkk[max_i];
                    }
                    for(size_t k2=0; k2<n_h; ++k2){
                        col_i = f_geno * n_h + k2;
                        hap_prob = prob_i[possiblehap[col_i]];
                        vit[k2] = max_scores_jk[k2] + hap_prob;
                    }
                }
                if(m == n_m - 1){
                    o_seq_i[m] = get_max_int(vit);
                }
            }

            // Backtracking
            for(size_t m=n_m-1; m>0; --m){
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
                 NumericMatrix trans_prob,
                 NumericVector init_prob,
                 int & n_p,
                 int & n_h,
                 int & n_o,
                 int & n_f,
                 int & n_m,
                 LogicalVector het,
                 IntegerVector possiblehap,
                 IntegerVector possiblegeno,
                 IntegerVector p_geno_fix,
                 IntegerVector ploidy
){
    // Initialize arrays to store output, alpha values,
    // emittion probs, and beta values.
    double neg_inf = -numeric_limits<double>::infinity();
    IntegerVector f_seq(n_m);
    IntegerMatrix f_path(n_m, n_p);
    IntegerMatrix o_seq(n_m, n_o);

    // Convert values to ones used here.
    IntegerVector dim = {n_m, n_f, n_o, n_p, n_h};
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
    LogicalVector iter_p_pat(n_p);
    NumericVector p_emit1;
    NumericVector p_emit2;

    // Initialize Viterbi trellis.
    bool p_het = het[0];
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
                        n_f,
                        n_p,
                        p_het,
                        ploidy);

    int int_fix_p = p_geno_fix[0];
    if(int_fix_p >= 0){
        R_xlen_t fix_p = p_geno_fix[0];
        for(R_xlen_t p=0; p<p_emit1.size(); ++p){
            if(p == fix_p){
                p_emit1[p] = 0;
            } else {
                p_emit1[p] = neg_inf;
            }
        }
    }

    IntegerVector valid_p_indices1;
    for(R_xlen_t j=0; j<p_emit1.size(); ++j){
        if(!isinf(p_emit1[j])){
            valid_p_indices1.push_back(j);
        }
    }
    size_t valid_size = valid_p_indices1.size();
    NumericMatrix vit_score(n_o, valid_size * n_h);
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
                        dim,
                        valid_p_indices1,
                        in_m,
                        ploidy);
    parallelFor(0, iter_sample.length(), init_vit);
    NumericMatrix in_score;
    in_score = clone(vit_score);

    for(int m=1; m<n_m; ++m){
        if(m % 10 == 9){
            Rcpp::Rcout << "\r" <<
                "Forward founder genotype probability calculation at marker#: " <<
                    m+1 << string(70, ' ');
        }
        NumericMatrix trans_prob_m = trans_prob( _ ,
                                                 Range( (m-1)*n_h ,
                                                        (m-1)*n_h + n_h - 1 ) );

        ParCalcVitFounder calc_vit(in_score,
                                   iter_sample,
                                   trans_prob_m,
                                   dim,
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
                            n_f,
                            n_p,
                            p_het,
                            ploidy);

        R_xlen_t fix_p_len = p_geno_fix.size();
        if(fix_p_len > m){
            R_xlen_t fix_p = p_geno_fix[m];
            for(R_xlen_t p=0; p<p_emit2.size(); ++p){
                if(p == fix_p){
                    p_emit2[p] = 0;
                } else {
                    p_emit2[p] = neg_inf;
                }
            }
        };

        IntegerVector valid_p_indices2;
        for(R_xlen_t j=0; j<p_emit2.size(); ++j){
            if(!isinf(p_emit2[j])){
                valid_p_indices2.push_back(j);
            }
        }

        size_t valid_size = valid_p_indices2.size();
        NumericMatrix vit_score(n_o, valid_size * n_h);
        IntegerVector in_m = {m};
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
                                     dim,
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

        if(m==n_m-1){

            last_vit_founder(f_seq,
                             in_score,
                             p_emit1,
                             dim,
                             valid_p_indices1);
        }
    }
    Rcpp::Rcout << "\r" <<
        "Forward founder genotype probability calculation: Done!" <<
            string(70, ' ');

    backtrack(f_path,
              f_seq,
              dim);

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
                                  dim,
                                  f_seq,
                                  ploidy);
    parallelFor(0, iter_sample.length(), vit_offspring);

    Rcpp::Rcout << "\r" << string(70, ' ');
    List out_list = List::create(_["p_geno"] = f_seq, _["best_seq"] = o_seq);
    return out_list;
}
