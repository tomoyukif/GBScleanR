#ifndef GBSRUTIL_H
#define GBSRUTIL_H

#include <Rcpp.h>

void log10_safe(double & d);
double log10_safe_d(const double & d);
double logsum(std::vector<double> & v);
Rcpp::NumericVector lognorm(Rcpp::NumericVector v);
void lognorm_vec(std::vector<double> & v);
size_t get_max_int(std::vector<double> & v);

#endif
