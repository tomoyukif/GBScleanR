// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerVector count_geno(IntegerMatrix geno){
  IntegerVector out(7);
  int g;
  for(R_xlen_t i = 0; i < geno.ncol(); i ++){
    if(geno(0, i) == 3){
      g = 3;
    } else {
      g = geno(0, i) + geno(1, i);
    }
    out[g] += 1;
  }

  out[4] = out[1] + out[0] * 2;
  out[5] = out[1] + out[2] * 2;
  out[6] = out[3] * 2;

  return out;
}

// [[Rcpp::export]]
NumericVector count_read(NumericVector read,
                         NumericVector tot_read){
  NumericVector out(8);
  NumericVector ref(read.size() / 2);
  NumericVector alt(read.size() / 2);
  R_xlen_t ref_i;
  R_xlen_t alt_i;

  if(tot_read.size() == 1){
    for(R_xlen_t i = 0; i < ref.size(); i ++){
      ref_i = i * 2;
      alt_i = i * 2 + 1;

      if(read[ref_i] == 0){
        ref[i] = NA_REAL;
      } else {
        ref[i] = read[ref_i];
      }
      if(read[alt_i] == 0){
        alt[i] = NA_REAL;
      } else {
        alt[i] = read[alt_i];
      }
    }
    ref = na_omit(ref);
    alt = na_omit(alt);
    out[0] = sum(ref);
    out[1] = sum(alt);
    double dp = out[0] + out[1];
    ref = ref / dp * 1e+06;
    alt = alt / dp * 1e+06;
    out[2] = mean(ref);
    out[3] = mean(alt);
    out[4] = sd(ref);
    out[5] = sd(alt);
    out[6] = median(ref);
    out[7] = median(alt);

  } else {
    for(R_xlen_t i = 0; i < read.size(); i ++){
      if(i < ref.size()){
        if(read[i] == 0){
          ref[i] = NA_REAL;
        } else {
          ref[i] = read[i];
        }

      } else {
        alt_i = i - ref.size();
        if(read[i] == 0){
          alt[alt_i] = NA_REAL;
        } else {
          alt[alt_i] = read[i];
        }
      }
    }

    out[0] = sum(na_omit(ref));
    out[1] = sum(na_omit(alt));
    ref = ref / tot_read * 1e+06;
    alt = alt / tot_read * 1e+06;
    ref = na_omit(ref);
    alt = na_omit(alt);
    out[2] = mean(ref);
    out[3] = mean(alt);
    out[4] = sd(ref);
    out[5] = sd(alt);
    out[6] = median(ref);
    out[7] = median(alt);
  }

  return out;
}

// [[Rcpp::export]]
LogicalVector thinout_marker(IntegerVector chr,
                             IntegerVector pos,
                             IntegerVector missing_count,
                             int range){
  LogicalVector valid(pos.size(), true);
  R_xlen_t i = 0;
  R_xlen_t j = 1;
  double mar1;
  double mar2;
  while(true){
    if(chr[i] == chr[j]){
      mar1 = pos[i];
      mar2 = pos[j];
      if(mar2 - mar1 <= range){
        if(missing_count[i] == missing_count[j]){
          valid[j] = false;
          j = j + 1;
        } else {
          i = j;
          j = j + 1;
        }
      } else {
        i = j;
        j = j + 1;
      }
    } else {
      i = j;
      j = j + 1;
    }
    if(j >= pos.size()){
      break;
    }
  }
  return valid;
}

