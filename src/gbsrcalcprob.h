#ifndef GBSRCALCPROB_H
#define GBSRCALCPROB_H

#include <Rcpp.h>
#include <RcppParallel.h>

std::vector<double> calcGenoprob(const double & ref,
                                 const double & alt,
                                 const double & eseq0,
                                 const double & eseq1,
                                 const double & w1,
                                 const double & w2,
                                 const bool & het);

void calcMissmap(std::vector<double> & prob,
                 const double & mismap1,
                 const double & mismap2);

Rcpp::NumericVector calcPemit(Rcpp::NumericMatrix p_ref,
                              Rcpp::NumericMatrix p_alt,
                              Rcpp::NumericVector eseq,
                              Rcpp::NumericVector w1,
                              Rcpp::NumericVector w2,
                              Rcpp::NumericVector mismap1,
                              Rcpp::NumericVector mismap2,
                              Rcpp::IntegerVector possiblegeno,
                              int & m,
                              int & n_f,
                              int & n_p,
                              bool & het,
                              Rcpp::IntegerVector ploidy);

std::vector<double> calcEmit(RcppParallel::RMatrix<double> ref,
                             RcppParallel::RMatrix<double> alt,
                             RcppParallel::RVector<double> eseq,
                             RcppParallel::RVector<double> w1,
                             RcppParallel::RVector<double> w2,
                             RcppParallel::RVector<double> mismap1,
                             RcppParallel::RVector<double> mismap2,
                             size_t & m,
                             size_t & sample_i,
                             bool & het);

#endif
