// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// beta_update
arma::vec beta_update(arma::mat x, arma::mat z, double sigma2_beta, arma::vec w, arma::vec gamma, arma::vec theta_old);
RcppExport SEXP _GPCW_beta_update(SEXP xSEXP, SEXP zSEXP, SEXP sigma2_betaSEXP, SEXP wSEXP, SEXP gammaSEXP, SEXP theta_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_beta(sigma2_betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_old(theta_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(beta_update(x, z, sigma2_beta, w, gamma, theta_old));
    return rcpp_result_gen;
END_RCPP
}
// GPCW
Rcpp::List GPCW(int mcmc_samples, arma::vec y, arma::mat x, arma::mat z, double sigma2_beta, double alpha_phi0, double beta_phi0, double a_phi1, double b_phi1, double mhvar_phi1_trans, arma::vec beta_init, arma::vec theta_init, double phi0_init, double phi1_init);
RcppExport SEXP _GPCW_GPCW(SEXP mcmc_samplesSEXP, SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP sigma2_betaSEXP, SEXP alpha_phi0SEXP, SEXP beta_phi0SEXP, SEXP a_phi1SEXP, SEXP b_phi1SEXP, SEXP mhvar_phi1_transSEXP, SEXP beta_initSEXP, SEXP theta_initSEXP, SEXP phi0_initSEXP, SEXP phi1_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type mcmc_samples(mcmc_samplesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_beta(sigma2_betaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_phi0(alpha_phi0SEXP);
    Rcpp::traits::input_parameter< double >::type beta_phi0(beta_phi0SEXP);
    Rcpp::traits::input_parameter< double >::type a_phi1(a_phi1SEXP);
    Rcpp::traits::input_parameter< double >::type b_phi1(b_phi1SEXP);
    Rcpp::traits::input_parameter< double >::type mhvar_phi1_trans(mhvar_phi1_transSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< double >::type phi0_init(phi0_initSEXP);
    Rcpp::traits::input_parameter< double >::type phi1_init(phi1_initSEXP);
    rcpp_result_gen = Rcpp::wrap(GPCW(mcmc_samples, y, x, z, sigma2_beta, alpha_phi0, beta_phi0, a_phi1, b_phi1, mhvar_phi1_trans, beta_init, theta_init, phi0_init, phi1_init));
    return rcpp_result_gen;
END_RCPP
}
// neg_two_loglike_update
double neg_two_loglike_update(arma::vec y, arma::mat x, arma::mat z, arma::vec beta, arma::vec theta);
RcppExport SEXP _GPCW_neg_two_loglike_update(SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP betaSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(neg_two_loglike_update(y, x, z, beta, theta));
    return rcpp_result_gen;
END_RCPP
}
// phi0_update
double phi0_update(double alpha_phi0, double beta_phi0, arma::vec theta, arma::mat corr_inv);
RcppExport SEXP _GPCW_phi0_update(SEXP alpha_phi0SEXP, SEXP beta_phi0SEXP, SEXP thetaSEXP, SEXP corr_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha_phi0(alpha_phi0SEXP);
    Rcpp::traits::input_parameter< double >::type beta_phi0(beta_phi0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type corr_inv(corr_invSEXP);
    rcpp_result_gen = Rcpp::wrap(phi0_update(alpha_phi0, beta_phi0, theta, corr_inv));
    return rcpp_result_gen;
END_RCPP
}
// phi1_update
Rcpp::List phi1_update(double phi1_old, double phi0, arma::vec theta, Rcpp::List temporal_corr_info, double a_phi1, double b_phi1, double mhvar_phi1_trans, int acctot_phi1_trans);
RcppExport SEXP _GPCW_phi1_update(SEXP phi1_oldSEXP, SEXP phi0SEXP, SEXP thetaSEXP, SEXP temporal_corr_infoSEXP, SEXP a_phi1SEXP, SEXP b_phi1SEXP, SEXP mhvar_phi1_transSEXP, SEXP acctot_phi1_transSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type phi1_old(phi1_oldSEXP);
    Rcpp::traits::input_parameter< double >::type phi0(phi0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type temporal_corr_info(temporal_corr_infoSEXP);
    Rcpp::traits::input_parameter< double >::type a_phi1(a_phi1SEXP);
    Rcpp::traits::input_parameter< double >::type b_phi1(b_phi1SEXP);
    Rcpp::traits::input_parameter< double >::type mhvar_phi1_trans(mhvar_phi1_transSEXP);
    Rcpp::traits::input_parameter< int >::type acctot_phi1_trans(acctot_phi1_transSEXP);
    rcpp_result_gen = Rcpp::wrap(phi1_update(phi1_old, phi0, theta, temporal_corr_info, a_phi1, b_phi1, mhvar_phi1_trans, acctot_phi1_trans));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_pgdraw
arma::vec rcpp_pgdraw(double b, arma::vec c);
RcppExport SEXP _GPCW_rcpp_pgdraw(SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_pgdraw(b, c));
    return rcpp_result_gen;
END_RCPP
}
// temporal_corr_fun
Rcpp::List temporal_corr_fun(int p_z, double phi1);
RcppExport SEXP _GPCW_temporal_corr_fun(SEXP p_zSEXP, SEXP phi1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p_z(p_zSEXP);
    Rcpp::traits::input_parameter< double >::type phi1(phi1SEXP);
    rcpp_result_gen = Rcpp::wrap(temporal_corr_fun(p_z, phi1));
    return rcpp_result_gen;
END_RCPP
}
// theta_update
arma::vec theta_update(arma::mat x, arma::mat z, arma::vec w, arma::vec gamma, arma::vec beta, double phi0_old, arma::mat corr_inv);
RcppExport SEXP _GPCW_theta_update(SEXP xSEXP, SEXP zSEXP, SEXP wSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP phi0_oldSEXP, SEXP corr_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type phi0_old(phi0_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type corr_inv(corr_invSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_update(x, z, w, gamma, beta, phi0_old, corr_inv));
    return rcpp_result_gen;
END_RCPP
}
// w_update
Rcpp::List w_update(arma::vec y, arma::mat x, arma::mat z, arma::vec beta_old, arma::vec theta_old);
RcppExport SEXP _GPCW_w_update(SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP beta_oldSEXP, SEXP theta_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_old(theta_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(w_update(y, x, z, beta_old, theta_old));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GPCW_beta_update", (DL_FUNC) &_GPCW_beta_update, 6},
    {"_GPCW_GPCW", (DL_FUNC) &_GPCW_GPCW, 14},
    {"_GPCW_neg_two_loglike_update", (DL_FUNC) &_GPCW_neg_two_loglike_update, 5},
    {"_GPCW_phi0_update", (DL_FUNC) &_GPCW_phi0_update, 4},
    {"_GPCW_phi1_update", (DL_FUNC) &_GPCW_phi1_update, 8},
    {"_GPCW_rcpp_pgdraw", (DL_FUNC) &_GPCW_rcpp_pgdraw, 2},
    {"_GPCW_temporal_corr_fun", (DL_FUNC) &_GPCW_temporal_corr_fun, 2},
    {"_GPCW_theta_update", (DL_FUNC) &_GPCW_theta_update, 7},
    {"_GPCW_w_update", (DL_FUNC) &_GPCW_w_update, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_GPCW(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
