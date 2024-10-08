# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

run_fb <- function(ref, alt, eseq_in, bias, mismap, possiblehap, trans_prob, init_prob, nonzero_prob, n_pgeno, n_hap, n_offspring, n_marker, n_nonzero_prob, pedigree, p_geno, ploidy) {
    .Call(`_GBScleanR_run_fb`, ref, alt, eseq_in, bias, mismap, possiblehap, trans_prob, init_prob, nonzero_prob, n_pgeno, n_hap, n_offspring, n_marker, n_nonzero_prob, pedigree, p_geno, ploidy)
}

get_genocall <- function(ref, alt, eseq_in, bias, mismap, n_o, n_m, ploidy) {
    .Call(`_GBScleanR_get_genocall`, ref, alt, eseq_in, bias, mismap, n_o, n_m, ploidy)
}

count_geno <- function(geno) {
    .Call(`_GBScleanR_count_geno`, geno)
}

count_read <- function(read, tot_read) {
    .Call(`_GBScleanR_count_read`, read, tot_read)
}

thinout_marker <- function(chr, pos, missing_count, range) {
    .Call(`_GBScleanR_thinout_marker`, chr, pos, missing_count, range)
}

run_viterbi <- function(p_ref, p_alt, ref, alt, eseq_in, bias, mismap, trans_prob, init_prob, nonzero_prob, n_pgeno, n_hap, n_offspring, n_founder, n_marker, n_nonzero_prob, het, pedigree, possiblehap, possiblegeno, p_geno_fix, ploidy) {
    .Call(`_GBScleanR_run_viterbi`, p_ref, p_alt, ref, alt, eseq_in, bias, mismap, trans_prob, init_prob, nonzero_prob, n_pgeno, n_hap, n_offspring, n_founder, n_marker, n_nonzero_prob, het, pedigree, possiblehap, possiblegeno, p_geno_fix, ploidy)
}

