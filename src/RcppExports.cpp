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
Rcpp::List GPCW(int mcmc_samples, arma::vec y, arma::mat x, arma::mat z, double mhvar_phi_trans, Rcpp::Nullable<double> sigma2_beta_prior, Rcpp::Nullable<double> a_sigma2_theta_prior, Rcpp::Nullable<double> b_sigma2_theta_prior, Rcpp::Nullable<double> a_phi_prior, Rcpp::Nullable<double> b_phi_prior, Rcpp::Nullable<Rcpp::NumericVector> beta_init, Rcpp::Nullable<Rcpp::NumericVector> theta_init, Rcpp::Nullable<double> sigma2_theta_init, Rcpp::Nullable<double> phi_init);
RcppExport SEXP _GPCW_GPCW(SEXP mcmc_samplesSEXP, SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP mhvar_phi_transSEXP, SEXP sigma2_beta_priorSEXP, SEXP a_sigma2_theta_priorSEXP, SEXP b_sigma2_theta_priorSEXP, SEXP a_phi_priorSEXP, SEXP b_phi_priorSEXP, SEXP beta_initSEXP, SEXP theta_initSEXP, SEXP sigma2_theta_initSEXP, SEXP phi_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type mcmc_samples(mcmc_samplesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type mhvar_phi_trans(mhvar_phi_transSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type sigma2_beta_prior(sigma2_beta_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type a_sigma2_theta_prior(a_sigma2_theta_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type b_sigma2_theta_prior(b_sigma2_theta_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type a_phi_prior(a_phi_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type b_phi_prior(b_phi_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type sigma2_theta_init(sigma2_theta_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type phi_init(phi_initSEXP);
    rcpp_result_gen = Rcpp::wrap(GPCW(mcmc_samples, y, x, z, mhvar_phi_trans, sigma2_beta_prior, a_sigma2_theta_prior, b_sigma2_theta_prior, a_phi_prior, b_phi_prior, beta_init, theta_init, sigma2_theta_init, phi_init));
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
// phi_update
Rcpp::List phi_update(double phi_old, double sigma2_theta, arma::vec theta, Rcpp::List temporal_corr_info, double a_phi, double b_phi, double mhvar_phi_trans, double acctot_phi_trans);
RcppExport SEXP _GPCW_phi_update(SEXP phi_oldSEXP, SEXP sigma2_thetaSEXP, SEXP thetaSEXP, SEXP temporal_corr_infoSEXP, SEXP a_phiSEXP, SEXP b_phiSEXP, SEXP mhvar_phi_transSEXP, SEXP acctot_phi_transSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type phi_old(phi_oldSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_theta(sigma2_thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type temporal_corr_info(temporal_corr_infoSEXP);
    Rcpp::traits::input_parameter< double >::type a_phi(a_phiSEXP);
    Rcpp::traits::input_parameter< double >::type b_phi(b_phiSEXP);
    Rcpp::traits::input_parameter< double >::type mhvar_phi_trans(mhvar_phi_transSEXP);
    Rcpp::traits::input_parameter< double >::type acctot_phi_trans(acctot_phi_transSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_update(phi_old, sigma2_theta, theta, temporal_corr_info, a_phi, b_phi, mhvar_phi_trans, acctot_phi_trans));
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
// sigma2_theta_update
double sigma2_theta_update(double a_sigma2_theta, double b_sigma2_theta, arma::vec theta, arma::mat corr_inv);
RcppExport SEXP _GPCW_sigma2_theta_update(SEXP a_sigma2_thetaSEXP, SEXP b_sigma2_thetaSEXP, SEXP thetaSEXP, SEXP corr_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a_sigma2_theta(a_sigma2_thetaSEXP);
    Rcpp::traits::input_parameter< double >::type b_sigma2_theta(b_sigma2_thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type corr_inv(corr_invSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma2_theta_update(a_sigma2_theta, b_sigma2_theta, theta, corr_inv));
    return rcpp_result_gen;
END_RCPP
}
// temporal_corr_fun
Rcpp::List temporal_corr_fun(int p_z, double phi);
RcppExport SEXP _GPCW_temporal_corr_fun(SEXP p_zSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p_z(p_zSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    rcpp_result_gen = Rcpp::wrap(temporal_corr_fun(p_z, phi));
    return rcpp_result_gen;
END_RCPP
}
// theta_update
arma::vec theta_update(arma::mat x, arma::mat z, arma::vec w, arma::vec gamma, arma::vec beta, double sigma2_theta_old, arma::mat corr_inv);
RcppExport SEXP _GPCW_theta_update(SEXP xSEXP, SEXP zSEXP, SEXP wSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP sigma2_theta_oldSEXP, SEXP corr_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_theta_old(sigma2_theta_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type corr_inv(corr_invSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_update(x, z, w, gamma, beta, sigma2_theta_old, corr_inv));
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
    {"_GPCW_phi_update", (DL_FUNC) &_GPCW_phi_update, 8},
    {"_GPCW_rcpp_pgdraw", (DL_FUNC) &_GPCW_rcpp_pgdraw, 2},
    {"_GPCW_sigma2_theta_update", (DL_FUNC) &_GPCW_sigma2_theta_update, 4},
    {"_GPCW_temporal_corr_fun", (DL_FUNC) &_GPCW_temporal_corr_fun, 2},
    {"_GPCW_theta_update", (DL_FUNC) &_GPCW_theta_update, 7},
    {"_GPCW_w_update", (DL_FUNC) &_GPCW_w_update, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_GPCW(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
