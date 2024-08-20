// [[depends(RcppParallel)]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include "gbsrutil.h"
#include "gbsrcalcprob.h"
using namespace Rcpp;
using namespace RcppParallel;
using namespace std;

struct ParGenoProb : public Worker {

    RMatrix<int> genocall;
    const RVector<int> iter_sample;
    const RMatrix<double> ref;
    const RMatrix<double> alt;
    const RVector<double> eseq;
    const RVector<double> w1;
    const RVector<double> mismap1;
    const RVector<double> mismap2;
    const RVector<int> ploidy;

    ParGenoProb(LogicalMatrix genocall,
                const LogicalVector iter_sample,
                const NumericMatrix ref,
                const NumericMatrix alt,
                const NumericVector eseq,
                const NumericVector w1,
                const NumericVector mismap1,
                const NumericVector mismap2,
                const IntegerVector ploidy)
        : genocall(genocall),
          iter_sample(iter_sample),
          ref(ref),
          alt(alt),
          eseq(eseq),
          w1(w1),
          mismap1(mismap1),
          mismap2(mismap2),
          ploidy(ploidy) {}

    void operator()(size_t begin, size_t end) {
        bool het = true;
        for(RVector<int>::const_iterator i=iter_sample.begin() + begin;
            i<iter_sample.begin() + end; ++i){
            size_t sample_i = distance(iter_sample.begin(), i);

            RMatrix<int>::Row genocall_i = genocall.row(sample_i);
            RMatrix<double>::Row ref_i = ref.row(sample_i);
            RMatrix<double>::Row alt_i = alt.row(sample_i);
            vector<double> prob(ploidy[0]);
            int sel;
            double threshold = 0.99;

            for(size_t m=0; m<ref_i.size(); ++m){
                prob = calcGenoprob(ref_i[m], alt_i[m], eseq[0], eseq[1], w1[m], het, ploidy[0]);
                sel = get_max_int(prob);
                if(prob[sel] > threshold){
                    genocall_i[m] = 0;
                } else {
                    genocall_i[m] = 1;
                }
            }
        }
    }
};

// [[Rcpp::export]]
LogicalMatrix get_genocall(NumericMatrix ref,
                           NumericMatrix alt,
                           NumericVector eseq_in,
                           NumericVector bias,
                           NumericMatrix mismap,
                           int & n_o,
                           int & n_m,
                           IntegerVector ploidy

){
    // Initialize arrays to store output
    LogicalMatrix genocall(n_o,  n_m);

    // Convert values to ones used here.
    NumericVector w1(n_m);
    NumericVector w2(n_m);
    NumericVector eseq(2);
    NumericVector mismap1 = mismap( _ , 0 );
    NumericVector mismap2 = mismap( _ , 1 );
    eseq = clone(eseq_in);
    w1 = clone(bias);

    LogicalVector iter_sample(n_o);

    ParGenoProb calc_gp(genocall,
                        iter_sample,
                        ref,
                        alt,
                        eseq,
                        w1,
                        mismap1,
                        mismap2,
                        ploidy);

    parallelFor(0, iter_sample.length(), calc_gp);

    return genocall;
}
