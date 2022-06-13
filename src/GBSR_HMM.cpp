// [[depends(RcppParallel)]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include <random>
#include <vector>
using namespace Rcpp;
using namespace RcppParallel;
using namespace std;

// Calculate the log10 value
void log10_safe(double & d){
    if(d == 0){
        d = -pow(10, 100);
    } else {
        d = log10(d);
    }
}

double log10_safe_d(const double & d){
    double out;
    if(d == 0){
        out = -pow(10, 100);
    } else {
        out = log10(d);
    }
    return out;
}

// Calculate the log10 of the sum of probabilities.
double logsum(vector<double> & v){
    if(v.size() == 1){
        return v[0];
    }
    double log_sum = 0;
    vector<double>::size_type count = 0;
    double neg_inf = -numeric_limits<double>::infinity();
    double v_max = *max_element(v.begin(), v.end());
    if(isinf(v_max)){
        return neg_inf;
    };

    for(vector<double>::size_type i=0; i<v.size(); ++i){
        if(!isinf(v[i])){
            count = i + 1;
            log_sum = v[i];
            break;
        }
    }

    for(vector<double>::size_type i=count; i<v.size(); ++i){
        if(log_sum > v[i]){
            log_sum = log_sum + log10(1 + pow(10, v[i] - log_sum));
        } else {
            log_sum = v[i] + log10(1 + pow(10, log_sum - v[i]));
        }
    }
    return log_sum;
}

// Calculate the log10 values of normalized probabilities
// for an Rcpp's NumericVector object.
NumericVector lognorm(NumericVector v){
    double log_sum = 0;
    vector<double>::size_type count = 0;
    double v_max = *max_element(v.begin(), v.end());
    if(isinf(v_max)){
        double even = 1 / (double)v.size();
        log10_safe(even);
        v.fill(even);
        return v;
    };
    for(R_xlen_t i=0; i<v.size(); ++i){
        if(!isinf(v[i])){
            count = i + 1;
            log_sum = v[i];
            break;
        }
    }

    for(R_xlen_t i=count; i<v.size(); ++i){
        if(log_sum > v[i]){
            log_sum = log_sum + log10(1 + pow(10, v[i] - log_sum));
        } else {
            log_sum = v[i] + log10(1 + pow(10, log_sum - v[i]));
        }
    }

    v = v - log_sum;
    return v;
}

// Calculate the log10 values of normalized probabilities
// for a std::vector object.
void lognorm_vec(vector<double> & v){
    double log_sum = 0;
    vector<double>::size_type count = 0;
    double v_max = *max_element(v.begin(), v.end());

    if(isinf(v_max)){
        double even = 1 / (double)v.size();
        log10_safe(even);
        for(vector<double>::size_type i=0; i<v.size(); ++i){
            v[i] = even;
        }
    } else {
        for(vector<double>::size_type i=0; i<v.size(); ++i){
            if(!isinf(v[i])){
                count = i + 1;
                log_sum = v[i];
                break;
            }
        }

        for(vector<double>::size_type i=count; i<v.size(); ++i){
            if(log_sum > v[i]){
                log_sum = log_sum + log10(1 + pow(10, v[i] - log_sum));
            } else {
                log_sum = v[i] + log10(1 + pow(10, log_sum - v[i]));
            }
        }

        for(vector<double>::size_type i=0; i<v.size(); ++i){
            v[i] = v[i] - log_sum;
        }
    }
}

// Power to ten
double pow10(double & d){
    return pow(10, d);
}

// Get the index of maximum value in a vector.
size_t get_max_int(vector<double> & v){
    size_t out_index;
    vector<size_t> max_indices;
    double v_max;
    bool check;
    v_max = *max_element(v.begin(), v.end());

    for(vector<double>::size_type l=0;l<v.size(); ++l){
        check = fabs(v.at(l) - v_max) < 0.0000000001;
        if(check){
            max_indices.push_back(l);
        }
    }

    if(max_indices.size() == 0){
        random_device rnd;
        mt19937 mt(rnd());
        int tmp_len = v.size();
        uniform_int_distribution<> rand1(0, tmp_len - 1);
        out_index = rand1(mt);
        return out_index;
    }

    if(max_indices.size() == 1){
        return max_indices[0];
    }

    random_device rnd;
    mt19937 mt(rnd());
    int tmp_len = max_indices.size();
    uniform_int_distribution<> rand1(0, tmp_len - 1);
    int tmp = rand1(mt);
    out_index = max_indices[tmp];
    return out_index;
}


// Function to calculate the value of the probability density function.
double calcpdf(const double & ratio,
               const double & mean){
    double var = 0.25;
    double inv_sqrt_2pi = 0.3989422804014327;
    double a = (ratio - mean) / var;
    return inv_sqrt_2pi / var * exp(-0.5 * a * a);
}

void setHetZero(vector<double> & prob){
    prob[1] = 0;
}

vector<double> calcGenoprob(const double & ref,
                            const double & alt,
                            const double & eseq0,
                            const double & eseq1,
                            const double & w1,
                            const double & w2,
                            const int & het){
    vector<double> prob(3);
    const double dp = ref + alt;

    if(dp > 5){
        const double ratio = ref / dp;
        prob[0] = calcpdf(ratio, eseq0);
        prob[1] = calcpdf(ratio, w1);
        prob[2] = calcpdf(ratio, eseq1);

        if(het){ setHetZero(prob); }

        double sum_prob;
        for(int g=0; g<3;++g){
            sum_prob += prob[g];
        }

        for(int g=0; g<3;++g){
            prob[g] = prob[g] / sum_prob;
        }

    } else {
        double logeseq0 = log10_safe_d(eseq0);
        double logeseq1 = log10_safe_d(eseq1);
        double logw1 = log10_safe_d(w1);
        double logw2 = log10_safe_d(w2);
        vector<double> ref_multiplier = {logeseq0, logw1, logeseq1};
        vector<double> alt_multiplier = {logeseq1, logw2, logeseq0};
        for(int g=0; g<3;++g){
            prob[g] = ref * ref_multiplier[g] +
                alt * alt_multiplier[g];
        }

        if(het){ setHetZero(prob); }

        lognorm_vec(prob);
        for(int g=0; g<3;++g){
            prob[g] = pow10(prob[g]);
        }
    }
    return prob;
}

// Calculate mismap accounted genotype probabilities
void calcMissmap(vector<double> & prob,
                 const double & mismap1,
                 const double & mismap2){
    vector<double> v1 = {1 - mismap1, mismap1, 0};
    vector<double> v2 = {0, 1, 0};
    vector<double> v3 = {0, mismap2, 1 - mismap2};
    double sum_v1 = 0.0;
    double sum_v2 = 0.0;
    double sum_v3 = 0.0;
    double sum_v = 0.0;

    for(size_t g=0; g<3;++g){
        sum_v1 += v1[g] * prob[g];
        sum_v2 += v2[g] * prob[g];
        sum_v3 += v3[g] * prob[g];
    }
    sum_v = sum_v1 + sum_v2 + sum_v3;
    if(sum_v == 0){
        prob.assign(3, 1/3);
    } else {
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
                        int & n_f,
                        int & n_p,
                        LogicalVector het,
                        IntegerVector ploidy
){
    vector<double> prob;
    double p_prob;
    int col_i;
    NumericVector p_emit(n_p, 1.0);

    for(int i=0; i<n_f; ++i){
        NumericMatrix::Row ref_i = p_ref.row(i);
        NumericMatrix::Row alt_i = p_alt.row(i);

        prob = calcGenoprob(ref_i[m], alt_i[m],
                            eseq[0], eseq[1],
                                         w1[m], w2[m], het[0]);
        for(int j=0; j<n_p; ++j){
            col_i = j * n_f + i;
            p_prob = prob[possiblegeno[col_i]];

            if(p_prob < 0.01){
                p_prob = 0;
            }
            p_emit[j] = p_emit[j] * p_prob;

        }
    }

    for(int j=0; j<n_p; ++j){
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
                        size_t & m,
                        size_t & sample_i,
                        int & het
){
    vector<double> prob;
    RMatrix<double>::Row ref_i = ref.row(sample_i);
    RMatrix<double>::Row alt_i = alt.row(sample_i);

    prob = calcGenoprob(ref_i[m], alt_i[m], eseq[0], eseq[1], w1[m], w2[m], het);
    calcMissmap(prob, mismap1[m], mismap2[m]);

    return prob;
}

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
        int het = 0;
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
                    log10_safe(hap_prob);
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
        int het = 0;
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
                                log10_safe(hap_prob);
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
                        log10_safe(hap_prob);
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
        int het = 0;
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
                        log10_safe(hap_prob);
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
                        log10_safe(hap_prob);
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
                        het,
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
                            het,
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
        int het = 0;
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
            vector<double> beta(n_h, 1);

            vector<double> score_k(n_h);
            double hap_prob;
            double trans_kk;
            double sum_k;
            size_t col_i;
            size_t j;

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
                        log10_safe(hap_prob);
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
                        log10_safe(hap_prob);
                        emit[m][k2] = hap_prob;
                        alpha[m][k2] = hap_prob + sum_k;
                    }
                }
            }

            vector<double> gamma_i1;
            vector<double> gamma_i2;
            vector<double> gamma_i3;
            vector<double> gamma_i_tmp(3);
            for(size_t k=0; k<n_h; ++k){
                size_t j = p_geno[n_m-1];
                col_i = j * n_h + k;
                if(possiblehap[col_i] == 0){
                    gamma_i1.push_back(beta[k] + alpha[n_m-1][k]);
                }
                if(possiblehap[col_i] == 1){
                    gamma_i2.push_back(beta[k] + alpha[n_m-1][k]);
                }
                if(possiblehap[col_i] == 2){
                    gamma_i3.push_back(beta[k] + alpha[n_m-1][k]);
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
                vector<double> gamma_i1;
                vector<double> gamma_i2;
                vector<double> gamma_i3;
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
                    if(possiblehap[col_i] == 0){
                        gamma_i1.push_back(sum_k + alpha[m-1][k1]);
                    }
                    if(possiblehap[col_i] == 1){
                        gamma_i2.push_back(sum_k + alpha[m-1][k1]);
                    }
                    if(possiblehap[col_i] == 2){
                        gamma_i3.push_back(sum_k + alpha[m-1][k1]);
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
    // emittion probs, and beta values.
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

////////////////////////////////////////////////////////////////////////////////

// Function to calculate probabilities of founder genotype patterns.
vector<double> calcDosProb(RMatrix<double> ref,
                           RMatrix<double> alt,
                           RMatrix<double> ratio,
                           RVector<double> eseq,
                           size_t & m,
                           size_t & sample_i,
                           RVector<int> ploidy,
                           RVector<int> mindp
){

    RMatrix<double>::Row ref_i = ref.row(sample_i);
    RMatrix<double>::Row alt_i = alt.row(sample_i);
    RMatrix<double>::Row ratio_i = ratio.row(sample_i);

    size_t plex = ploidy[0] + 1;
    vector<double> prob(plex);
    vector<double> w(plex);
    const double dp = ref_i[m] + alt_i[m];

    double p;
    for(size_t i = 0; i < plex; ++i){
        p = (1 / (double)ploidy[0]) * (double)i;
        if(p > eseq[0]){
            w[i] = eseq[0];
        } else if(p < eseq[1]){
            w[i] = eseq[1];
        } else {
            w[i] = p;
        }
    }

    if(ref_i[m] == -1){
        if(ratio_i[m] < 0){
            double even_p = 1 / (ploidy[0] + 1);
            for(size_t i = 0; i < plex; ++i){
                prob[i] = even_p;
            }
        } else {
            for(size_t i = 0; i < plex; ++i){
                prob[i] = calcpdf(ratio_i[m], w[i]);
            }
        }

        double sum_prob;
        for(size_t g = 0; g < plex; ++g){
            sum_prob += prob[g];
        }
        for(size_t g = 0; g < plex; ++g){
            prob[g] = prob[g] / sum_prob;
        }

    } else {
        if(dp > mindp[0]){
            if(ratio_i[m] < 0){
                double even_p = 1 / plex;
                for(size_t i = 0; i < plex; ++i){
                    prob[i] = even_p;
                }
            } else {
                for(size_t i = 0; i < plex; ++i){
                    prob[i] = calcpdf(ratio_i[m], w[i]);
                }
            }

            double sum_prob;
            for(size_t g = 0; g < plex; ++g){
                sum_prob += prob[g];
            }
            for(size_t g = 0; g < plex; ++g){
                prob[g] = prob[g] / sum_prob;
            }

        } else {

            for(size_t i = 0; i < plex; ++i){
                w[i] = log10_safe_d(w[i]);
            }
            for(size_t g = 0; g < plex; ++g){
                prob[g] = ref_i[m] * w[g] + alt_i[m] * w[ploidy[0] - g];
            }

            lognorm_vec(prob);
            for(size_t g = 0; g < plex; ++g){
                prob[g] = pow10(prob[g]);
            }
        }
    }
    double sum_prob;
    for(size_t i = 0; i < plex; ++i){
        sum_prob += prob[i];
    }
    if(sum_prob == 0){
        prob.assign(plex, 1/plex);
    }

    return prob;
}

// Functions for the Viterbi algorithm For OFFSPRING
struct ParVitDosage : public Worker {

    RMatrix<int> o_dos;
    const RVector<int> iter_sample;
    const RMatrix<double> ref;
    const RMatrix<double> alt;
    const RMatrix<double> ratio;
    const RVector<double> eseq;
    const RVector<double> init_prob;
    const RMatrix<double> trans_prob;
    const RVector<int> dim;
    const RVector<int> ploidy;
    const RVector<int> mindp;

    ParVitDosage(IntegerMatrix o_dos,
                 const LogicalVector iter_sample,
                 const NumericMatrix ref,
                 const NumericMatrix alt,
                 const NumericMatrix ratio,
                 const NumericVector eseq,
                 const NumericVector init_prob,
                 const NumericMatrix trans_prob,
                 const IntegerVector dim,
                 const IntegerVector ploidy,
                 const IntegerVector mindp)
        : o_dos(o_dos),
          iter_sample(iter_sample),
          ref(ref),
          alt(alt),
          ratio(ratio),
          eseq(eseq),
          init_prob(init_prob),
          trans_prob(trans_prob),
          dim(dim),
          ploidy(ploidy),
          mindp(mindp) {}

    void operator()(size_t begin, size_t end) {
        for(RVector<int>::const_iterator i=iter_sample.begin() + begin;
            i<iter_sample.begin() + end; ++i){
            size_t sample_i = distance(iter_sample.begin(), i);
            RMatrix<int>::Column o_dos_i = o_dos.column(sample_i);
            size_t n_m = dim[0];
            size_t plex = ploidy[0] + 1;
            vector<vector<unsigned short>> o_path(n_m,
                                                  vector<unsigned short>(plex));
            vector<double> vit(plex);
            double trans_kk;
            size_t max_i;
            vector<double> score_jkk(plex);
            vector<double> max_scores_jk(plex);
            size_t o_prev;

            for(size_t m=0; m<n_m; ++m){
                vector<double> prob_i = calcDosProb(ref,
                                                    alt,
                                                    ratio,
                                                    eseq,
                                                    m,
                                                    sample_i,
                                                    ploidy,
                                                    mindp);

                for(size_t i = 0; i < prob_i.size(); ++i){
                    log10_safe(prob_i[i]);
                }

                size_t trans_prob_col;

                if(m == 0){
                    for(size_t k=0; k<plex; ++k){
                        vit[k] = prob_i[k] + init_prob[k];
                    }
                } else {
                    for(size_t k2 = 0; k2 < plex; ++k2){
                        trans_prob_col = (m-1) * plex + k2;
                        RMatrix<double>::Column trans_prob_k =
                            trans_prob.column(trans_prob_col);

                        for(size_t k1 = 0; k1 < plex; ++k1){
                            trans_kk = trans_prob_k[k1];
                            score_jkk[k1] = vit[k1] + trans_kk;
                        }
                        max_i = get_max_int(score_jkk);
                        o_path[m][k2] = max_i;
                        max_scores_jk[k2] = score_jkk[max_i];
                    }
                    for(size_t k2 = 0; k2 < plex; ++k2){
                        vit[k2] = max_scores_jk[k2] + prob_i[k2];
                    }
                }
                if(m == n_m - 1){
                    o_dos_i[m] = get_max_int(vit);
                }
            }

            // Backtracking
            for(size_t m=n_m-1; m>0; --m){
                o_prev = o_dos_i[m];
                o_dos_i[m-1] = o_path[m][o_prev];
            }
        }
    }
};


// Solve the HMM
// [[Rcpp::export]]
IntegerMatrix dosage_viterbi(NumericMatrix ref,
                             NumericMatrix alt,
                             NumericMatrix ratio,
                             NumericVector eseq_in,
                             NumericMatrix trans_prob,
                             NumericVector init_prob,
                             int & n_o,
                             int & n_m,
                             IntegerVector & ploidy,
                             IntegerVector & mindp
){
    // Initialize arrays to store output, alpha values,
    // emittion probs, and beta values.
    IntegerMatrix o_dos(n_m, n_o);

    // Convert values to ones used here.
    IntegerVector dim = {n_m, n_o};
    NumericVector eseq(2);
    eseq = clone(eseq_in);

    LogicalVector iter_sample(n_o);

    ParVitDosage vit_dosage(o_dos,
                            iter_sample,
                            ref,
                            alt,
                            ratio,
                            eseq,
                            init_prob,
                            trans_prob,
                            dim,
                            ploidy,
                            mindp);

    parallelFor(0, iter_sample.length(), vit_dosage);

    Rcpp::Rcout << "\r" << string(70, ' ');
    return o_dos;
}
