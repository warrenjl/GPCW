#ifndef __GPCW__
#define __GPCW__

arma::vec rcpp_pgdraw(double b, 
                      arma::vec c);

Rcpp::List temporal_corr_fun(int p_z,
                             double phi1);

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
                       double phi0_old,
                       arma::mat corr_inv);

double phi0_update(double alpha_phi0,
                   double beta_phi0,
                   arma::vec theta,
                   arma::mat corr_inv);

Rcpp::List phi1_update(double phi1_old,
                       double phi0,
                       arma::vec theta,
                       Rcpp::List temporal_corr_info,
                       double a_phi1,
                       double b_phi1,
                       double mhvar_phi1_trans,
                       int acctot_phi1_trans);

double neg_two_loglike_update(arma::vec y,
                              arma::mat x,
                              arma::mat z, 
                              arma::vec beta,
                              arma::vec theta);

Rcpp::List GPCW(int mcmc_samples,
                arma::vec y,
                arma::mat x,
                arma::mat z,
                double sigma2_beta,
                double alpha_phi0,
                double beta_phi0,
                double a_phi1,
                double b_phi1,
                double mhvar_phi1_trans,
                arma::vec beta_init,
                arma::vec theta_init,
                double phi0_init,
                double phi1_init,
                double neg_two_loglike_init); 

#endif // __GPCW__
