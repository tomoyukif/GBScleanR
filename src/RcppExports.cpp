// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// run_viterbi
List run_viterbi(NumericMatrix p_ref, NumericMatrix p_alt, NumericMatrix ref, NumericMatrix alt, NumericVector eseq_in, NumericVector bias, NumericMatrix mismap, NumericMatrix trans_prob, NumericVector init_prob, int& n_p, int& n_h, int& n_o, int& n_f, int& n_m, LogicalVector nohet, IntegerVector possiblehap, IntegerVector possiblegeno, IntegerVector p_geno_fix);
RcppExport SEXP _GBScleanR_run_viterbi(SEXP p_refSEXP, SEXP p_altSEXP, SEXP refSEXP, SEXP altSEXP, SEXP eseq_inSEXP, SEXP biasSEXP, SEXP mismapSEXP, SEXP trans_probSEXP, SEXP init_probSEXP, SEXP n_pSEXP, SEXP n_hSEXP, SEXP n_oSEXP, SEXP n_fSEXP, SEXP n_mSEXP, SEXP nohetSEXP, SEXP possiblehapSEXP, SEXP possiblegenoSEXP, SEXP p_geno_fixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p_ref(p_refSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p_alt(p_altSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ref(refSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alt(altSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eseq_in(eseq_inSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bias(biasSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mismap(mismapSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type trans_prob(trans_probSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type init_prob(init_probSEXP);
    Rcpp::traits::input_parameter< int& >::type n_p(n_pSEXP);
    Rcpp::traits::input_parameter< int& >::type n_h(n_hSEXP);
    Rcpp::traits::input_parameter< int& >::type n_o(n_oSEXP);
    Rcpp::traits::input_parameter< int& >::type n_f(n_fSEXP);
    Rcpp::traits::input_parameter< int& >::type n_m(n_mSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type nohet(nohetSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type possiblehap(possiblehapSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type possiblegeno(possiblegenoSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p_geno_fix(p_geno_fixSEXP);
    rcpp_result_gen = Rcpp::wrap(run_viterbi(p_ref, p_alt, ref, alt, eseq_in, bias, mismap, trans_prob, init_prob, n_p, n_h, n_o, n_f, n_m, nohet, possiblehap, possiblegeno, p_geno_fix));
    return rcpp_result_gen;
END_RCPP
}
// run_fb
NumericMatrix run_fb(NumericMatrix ref, NumericMatrix alt, NumericVector eseq_in, NumericVector bias, NumericMatrix mismap, IntegerVector possiblehap, NumericMatrix trans_prob, NumericVector init_prob, int& n_h, int& n_o, int& n_m, IntegerVector p_geno);
RcppExport SEXP _GBScleanR_run_fb(SEXP refSEXP, SEXP altSEXP, SEXP eseq_inSEXP, SEXP biasSEXP, SEXP mismapSEXP, SEXP possiblehapSEXP, SEXP trans_probSEXP, SEXP init_probSEXP, SEXP n_hSEXP, SEXP n_oSEXP, SEXP n_mSEXP, SEXP p_genoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ref(refSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alt(altSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eseq_in(eseq_inSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bias(biasSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mismap(mismapSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type possiblehap(possiblehapSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type trans_prob(trans_probSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type init_prob(init_probSEXP);
    Rcpp::traits::input_parameter< int& >::type n_h(n_hSEXP);
    Rcpp::traits::input_parameter< int& >::type n_o(n_oSEXP);
    Rcpp::traits::input_parameter< int& >::type n_m(n_mSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p_geno(p_genoSEXP);
    rcpp_result_gen = Rcpp::wrap(run_fb(ref, alt, eseq_in, bias, mismap, possiblehap, trans_prob, init_prob, n_h, n_o, n_m, p_geno));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GBScleanR_run_viterbi", (DL_FUNC) &_GBScleanR_run_viterbi, 18},
    {"_GBScleanR_run_fb", (DL_FUNC) &_GBScleanR_run_fb, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_GBScleanR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}