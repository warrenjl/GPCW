#include "RcppArmadillo.h"
#include "GPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List GPCW(int mcmc_samples,
                arma::vec y,
                arma::mat x,
                arma::mat z,
                double mhvar_phi1_trans,
                Rcpp::Nullable<double> sigma2_beta_prior = R_NilValue,
                Rcpp::Nullable<double> a_phi0_prior = R_NilValue,
                Rcpp::Nullable<double> b_phi0_prior = R_NilValue,
                Rcpp::Nullable<double> a_phi1_prior = R_NilValue,
                Rcpp::Nullable<double> b_phi1_prior = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> theta_init = R_NilValue,
                Rcpp::Nullable<double> phi0_init = R_NilValue,
                Rcpp::Nullable<double> phi1_init = R_NilValue){

//Defining Parameters and Quantities of Interest
arma::mat beta(x.n_cols, mcmc_samples); beta.fill(0);
arma::mat theta(z.n_cols, mcmc_samples); theta.fill(0);
arma::vec phi0(mcmc_samples); phi0.fill(0);
arma::vec phi1(mcmc_samples); phi1.fill(0);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0);

//Prior Information
double sigma2_beta = 10000;
if(sigma2_beta_prior.isNotNull()){
  sigma2_beta = Rcpp::as<double>(sigma2_beta_prior);
  }

double a_phi0 = 3.00;
if(a_phi0_prior.isNotNull()){
  a_phi0 = Rcpp::as<double>(a_phi0_prior);
  }
  
double b_phi0 = 2.00;
if(b_phi0_prior.isNotNull()){
  b_phi0 = Rcpp::as<double>(b_phi0_prior);
  }

double a_phi1 = log(0.9999)/(-(z.n_cols - 1));
if(a_phi1_prior.isNotNull()){
  a_phi1 = Rcpp::as<double>(a_phi1_prior);
  }
  
double b_phi1 = log(0.0001)/(-1);
if(b_phi1_prior.isNotNull()){
  b_phi1 = Rcpp::as<double>(b_phi1_prior);
  }

//Initial Values
beta.col(0).fill(0);
if(beta_init.isNotNull()){
  beta.col(0) = Rcpp::as<arma::vec>(beta_init);
  }

theta.col(0).fill(0);
if(theta_init.isNotNull()){
  theta.col(0) = Rcpp::as<arma::vec>(theta_init);
  }

phi0(0) = 0.25;
if(phi0_init.isNotNull()){
  phi0(0) = Rcpp::as<double>(phi0_init);
  }

phi1(0) = 0.01;
if(phi1_init.isNotNull()){
  phi1(0) = Rcpp::as<double>(phi1_init);
  }

Rcpp::List temporal_corr_info = temporal_corr_fun(z.n_cols, phi1(0));
neg_two_loglike(0) = neg_two_loglike_update(y,
                                            x,
                                            z, 
                                            beta.col(0),
                                            theta.col(0));

//Metropolis Settings
double acctot_phi1_trans = 0;

//Main Sampling Loop
for(int j = 1; j < mcmc_samples; ++j){
  
  //w Update
  Rcpp::List w_output = w_update(y,
                                 x,
                                 z,
                                 beta.col(j-1),
                                 theta.col(j-1));
  arma::vec w = w_output[0];
  arma::vec gamma = w_output[1];
  
  //beta Update
  beta.col(j) = beta_update(x, 
                            z,
                            sigma2_beta,
                            w,
                            gamma,
                            theta.col(j-1));
  
  //theta Update
  theta.col(j) = theta_update(x, 
                              z,
                              w,
                              gamma,
                              beta.col(j),
                              phi0(j-1),
                              temporal_corr_info(0));

  //phi0 Update
  phi0(j) = phi0_update(a_phi0,
                        b_phi0,
                        theta.col(j),
                        temporal_corr_info(0));
  
  //phi1 Update
  Rcpp::List phi1_output = phi1_update(phi1(j-1),
                                       phi0(j),
                                       theta.col(j),
                                       temporal_corr_info,
                                       a_phi1,
                                       b_phi1,
                                       mhvar_phi1_trans,
                                       acctot_phi1_trans);

  phi1(j) = phi1_output(0);
  acctot_phi1_trans = phi1_output(1);
  temporal_corr_info = phi1_output(2);

  //neg_two_loglike Update
  neg_two_loglike(j) = neg_two_loglike_update(y,
                                              x,
                                              z, 
                                              beta.col(j),
                                              theta.col(j));
  
  //Progress
  if(j % 10 == 0){ 
    Rcpp::checkUserInterrupt();
    }
  
  if(j % int(round(mcmc_samples*0.02)) == 0){
    double completion = round(100*(j/(double)mcmc_samples));
    Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
    double accrate_phi1_trans = round(100*(acctot_phi1_trans/j));
    Rcpp::Rcout << "phi1 Acceptance: " << accrate_phi1_trans << "%" << std::endl;
    Rcpp::Rcout << "********************" << std::endl;
    }
  
  }
                                  
return Rcpp::List::create(Rcpp::Named("beta") = beta,
                          Rcpp::Named("theta") = theta,
                          Rcpp::Named("phi0") = phi0,
                          Rcpp::Named("phi1") = phi1,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_phi1_trans") = acctot_phi1_trans);

}

