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
                double metrop_var_phi_trans,
                Rcpp::Nullable<double> sigma2_beta_prior = R_NilValue,
                Rcpp::Nullable<double> a_sigma2_theta_prior = R_NilValue,
                Rcpp::Nullable<double> b_sigma2_theta_prior = R_NilValue,
                Rcpp::Nullable<double> a_phi_prior = R_NilValue,
                Rcpp::Nullable<double> b_phi_prior = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> theta_init = R_NilValue,
                Rcpp::Nullable<double> sigma2_theta_init = R_NilValue,
                Rcpp::Nullable<double> phi_init = R_NilValue){

//Defining Parameters and Quantities of Interest
int p_x = x.n_cols;
int p_z = z.n_cols;
arma::mat beta(p_x, mcmc_samples); beta.fill(0.00);
arma::mat theta(p_z, mcmc_samples); theta.fill(0.00);
arma::vec sigma2_theta(mcmc_samples); sigma2_theta.fill(0.00);
arma::vec phi(mcmc_samples); phi.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

//Prior Information
double sigma2_beta = 10000.00;
if(sigma2_beta_prior.isNotNull()){
  sigma2_beta = Rcpp::as<double>(sigma2_beta_prior);
  }

double a_sigma2_theta = 3.00;
if(a_sigma2_theta_prior.isNotNull()){
  a_sigma2_theta = Rcpp::as<double>(a_sigma2_theta_prior);
  }
  
double b_sigma2_theta = 2.00;
if(b_sigma2_theta_prior.isNotNull()){
  b_sigma2_theta = Rcpp::as<double>(b_sigma2_theta_prior);
  }

double a_phi = log(0.9999)/(-(p_z - 1.00));
if(a_phi_prior.isNotNull()){
  a_phi = Rcpp::as<double>(a_phi_prior);
  }
  
double b_phi = log(0.0001)/(-1.00);
if(b_phi_prior.isNotNull()){
  b_phi = Rcpp::as<double>(b_phi_prior);
  }

//Initial Values
beta.col(0).fill(0.00);
if(beta_init.isNotNull()){
  beta.col(0) = Rcpp::as<arma::vec>(beta_init);
  }

theta.col(0).fill(0.00);
if(theta_init.isNotNull()){
  theta.col(0) = Rcpp::as<arma::vec>(theta_init);
  }

sigma2_theta(0) = 1.00;
if(sigma2_theta_init.isNotNull()){
  sigma2_theta(0) = Rcpp::as<double>(sigma2_theta_init);
  }

phi(0) = (b_phi - a_phi)*0.01;
if(phi_init.isNotNull()){
  phi(0) = Rcpp::as<double>(phi_init);
  }

Rcpp::List temporal_corr_info = temporal_corr_fun(p_z, phi(0));
neg_two_loglike(0) = neg_two_loglike_update(y,
                                            x,
                                            z, 
                                            beta.col(0),
                                            theta.col(0));

//Metropolis Settings
int acctot_phi_trans = 0;

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
                              sigma2_theta(j-1),
                              temporal_corr_info(0));

  //sigma2_theta Update
  sigma2_theta(j) = sigma2_theta_update(a_sigma2_theta,
                                        b_sigma2_theta,
                                        theta.col(j),
                                        temporal_corr_info(0));
  
  //phi Update
  Rcpp::List phi_output = phi_update(phi(j-1),
                                     sigma2_theta(j),
                                     theta.col(j),
                                     temporal_corr_info,
                                     a_phi,
                                     b_phi,
                                     metrop_var_phi_trans,
                                     acctot_phi_trans);

  phi(j) = Rcpp::as<double>(phi_output[0]);
  acctot_phi_trans = phi_output[1];
  temporal_corr_info = phi_output[2];

  //neg_two_loglike Update
  neg_two_loglike(j) = neg_two_loglike_update(y,
                                              x,
                                              z, 
                                              beta.col(j),
                                              theta.col(j));
  
  //Progress
  if((j + 1) % 10 == 0){ 
    Rcpp::checkUserInterrupt();
    }
  
  if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
    double completion = round(100*((j + 1)/(double)mcmc_samples));
    Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
    double accrate_phi_trans = round(100*(acctot_phi_trans/(double)j));
    Rcpp::Rcout << "phi Acceptance: " << accrate_phi_trans << "%" << std::endl;
    Rcpp::Rcout << "*******************" << std::endl;
    }
  
  }
                                  
return Rcpp::List::create(Rcpp::Named("beta") = beta,
                          Rcpp::Named("theta") = theta,
                          Rcpp::Named("sigma2_theta") = sigma2_theta,
                          Rcpp::Named("phi") = phi,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_phi_trans") = acctot_phi_trans);

}

