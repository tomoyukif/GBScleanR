// [[depends(RcppParallel)]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include <random>
using namespace Rcpp;
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
        if(!isinf(v[i])){
            if(log_sum > v[i]){
                log_sum = log_sum + log10(1 + pow(10, v[i] - log_sum));
            } else {
                log_sum = v[i] + log10(1 + pow(10, log_sum - v[i]));
            }
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
