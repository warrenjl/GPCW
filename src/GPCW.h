#ifndef __GPCW__
#define __GPCW__

arma::vec rcpp_pgdraw(double b, 
                      arma::vec c);

Rcpp::List temporal_corr_fun(int p_z,
                             double phi);

Rcpp::List w_update(arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec beta_old,
                    arma::vec theta_old);

arma::vec beta_update(arma::mat x, 
                      arma::mat z,
                      double sigma2_beta,
                      arma::vec w,
                      arma::vec gamma,
                      arma::vec theta_old);

arma::vec theta_update(arma::mat x, 
                       arma::mat z,
                       arma::vec w,
                       arma::vec gamma,
                       arma::vec beta,
                       double sigma2_theta_old,
                       arma::mat corr_inv);

double sigma2_theta_update(double a_sigma2_theta,
                           double b_sigma2_theta,
                           arma::vec theta,
                           arma::mat corr_inv);

Rcpp::List phi_update(double phi_old,
                      double sigma2_theta,
                      arma::vec theta,
                      Rcpp::List temporal_corr_info,
                      double a_phi,
                      double b_phi,
                      double mhvar_phi_trans,
                      double acctot_phi_trans);

double neg_two_loglike_update(arma::vec y,
                              arma::mat x,
                              arma::mat z, 
                              arma::vec beta,
                              arma::vec theta);

Rcpp::List GPCW(int mcmc_samples,
                arma::vec y,
                arma::mat x,
                arma::mat z,
                double mhvar_phi_trans,
                Rcpp::Nullable<double> sigma2_beta_prior,
                Rcpp::Nullable<double> a_sigma2_theta_prior,
                Rcpp::Nullable<double> b_sigma2_theta_prior,
                Rcpp::Nullable<double> a_phi_prior,
                Rcpp::Nullable<double> b_phi_prior,
                Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                Rcpp::Nullable<Rcpp::NumericVector> theta_init,
                Rcpp::Nullable<double> sigma2_theta_init,
                Rcpp::Nullable<double> phi_init); 

#endif // __GPCW__
