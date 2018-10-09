#include "RcppArmadillo.h"
#include "GPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List phi1_update(double phi1_old,
                       double phi0,
                       arma::vec theta,
                       Rcpp::List temporal_corr_info,
                       double a_phi1,
                       double b_phi1,
                       double mhvar_phi1_trans,
                       int acctot_phi1_trans){

/*Second*/
Rcpp::List temporal_corr_info_old = temporal_corr_info;
arma::mat corr_inv_old = temporal_corr_info_old[0];
double log_deter_old = temporal_corr_info_old[1];
double phi1_trans_old = log((phi1_old - a_phi1)/(b_phi1 - phi1_old));

double second = -0.50*log_deter_old - 
                (1/phi0)*0.50*dot(theta, (corr_inv_old*theta)) + 
                phi1_trans_old -
                2.0*log(1 + exp(phi1_trans_old));

/*First*/
double phi1_trans = R::rnorm(phi1_trans_old, 
                             sqrt(mhvar_phi1_trans));
double phi1 = (b_phi1*exp(phi1_trans) + a_phi1)/(exp(phi1_trans) + 1);
temporal_corr_info = temporal_corr_fun(theta.size(), phi1);
arma::mat corr_inv = temporal_corr_info[0];
double log_deter = temporal_corr_info[1];

double first = -0.50*log_deter - 
               (1/phi0)*0.50*dot(theta, (corr_inv*theta)) + 
               phi1_trans -
               2.0*log(1 + exp(phi1_trans));

/*Decision*/
double ratio = exp(first - second);   
double acc = 1;
if(ratio < R::runif(0, 1)){
  phi1 = phi1_old;
  temporal_corr_info = temporal_corr_info_old;
  acc = 0;
  }
acctot_phi1_trans = acctot_phi1_trans + 
                    acc;

return Rcpp::List::create(Rcpp::Named("phi1") = phi1,
                          Rcpp::Named("acctot_phi1_trans") = acctot_phi1_trans,
                          Rcpp::Named("temporal_corr_info") = temporal_corr_info);

}



