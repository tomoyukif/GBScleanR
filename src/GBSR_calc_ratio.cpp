// [[depends(RcppParallel)]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace RcppParallel;
using namespace std;


struct ParSR : public Worker {

    const RMatrix<double> ref;
    const RMatrix<double> alt;
    const RVector<double> pos;
    RMatrix<double> ratio;
    const RVector<double> window;
    const RVector<int> iter_sample;

    ParSR(const NumericMatrix ref,
          const NumericMatrix alt,
          const NumericVector pos,
          NumericMatrix ratio,
          const NumericVector window,
          const LogicalVector iter_sample)
        : ref(ref),
          alt(alt),
          pos(pos),
          ratio(ratio),
          window(window),
          iter_sample(iter_sample) {}

    void operator()(size_t begin, size_t end) {
        for(RVector<int>::const_iterator i = iter_sample.begin() + begin;
            i < iter_sample.begin() + end; ++i){
            size_t sample_i = distance(iter_sample.begin(), i);
            RMatrix<double>::Row ratio_i = ratio.row(sample_i);
            RMatrix<double>::Row ref_i = ref.row(sample_i);
            RMatrix<double>::Row alt_i = alt.row(sample_i);
            double w = window[0];

            // Calculate ratio
            double dp;
            vector<double> rr(pos.length());
            for(size_t i = 0; i < pos.length(); ++i){
                dp = ref_i[i] + alt_i[i];
                if(dp != 0){
                    rr[i] = ref_i[i] / dp;
                } else {
                    rr[i] = -1;
                }
            }

            double i_pos;
            size_t w_start;
            size_t w_end;
            double len;
            double w_ratio;
            double non_zero;
            for(size_t i = 0; i < rr.size(); ++i){
                w_ratio = 0;
                non_zero = 0;
                i_pos = pos[i];
                w_start = 0;
                w_end = rr.size() - 1;
                for(size_t m = i; m > 0; --m){
                    len = i_pos - pos[m];
                    if(len > w){
                        w_start = m;
                        break;
                    }
                }
                for(size_t m = i; m < rr.size(); ++m){
                    len = pos[m] - i_pos;
                    if(len > w){
                        w_end = m;
                        break;
                    }
                }
                for(size_t m = w_start; m < w_end; ++m){
                    if(rr[m] >= 0){
                        non_zero += 1;
                        w_ratio += rr[m];
                    }
                }
                ratio_i[i] = w_ratio / non_zero;
            }
        }
    }
};

// Calculate mean read ratio in the sliding window for each sample.
// [[Rcpp::export]]
NumericMatrix calc_ratio(NumericMatrix ref,
                         NumericMatrix alt,
                         NumericVector pos,
                         NumericVector window
){
    // Initialize arrays to store output read ratio.
    NumericMatrix ratio(ref.nrow(), ref.ncol());

    LogicalVector iter_sample(ref.nrow());

    ParSR smooth_ratio(ref,
                       alt,
                       pos,
                       ratio,
                       window,
                       iter_sample);

    parallelFor(0, iter_sample.length(), smooth_ratio);

    return ratio;
}
