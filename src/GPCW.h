#ifndef __GPCW__
#define __GPCW__

arma::vec rcpp_pgdraw(arma::vec b, 
                      arma::vec c);

Rcpp::NumericVector sampleRcpp(Rcpp::NumericVector x,
                               int size,
                               bool replace,
                               Rcpp::NumericVector prob = Rcpp::NumericVector::create());

Rcpp::List temporal_corr_fun(int p_z,
                             double phi);

int r_update(arma::vec y,
             arma::mat x,
             arma::mat z,
             arma::vec off_set,
             int a_r,
             int b_r,
             arma::vec beta,
             arma::vec theta);

double sigma2_epsilon_update(arma::vec y,
                             arma::mat x,
                             arma::mat z,
                             double a_sigma2_epsilon,
                             double b_sigma2_epsilon,
                             arma::vec beta_old,
                             arma::vec theta_old);

Rcpp::List w_update(arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec off_set,
                    arma::vec tri_als,
                    int likelihood_indicator,
                    int r,
                    arma::vec beta_old,
                    arma::vec theta_old);

arma::vec beta_update(arma::mat x, 
                      arma::mat z,
                      arma::vec off_set,
                      double sigma2_beta,
                      arma::vec w,
                      arma::vec gamma,
                      arma::vec theta_old);

arma::vec theta_update(arma::mat x, 
                       arma::mat z,
                       arma::vec off_set,
                       arma::vec w,
                       arma::vec gamma,
                       arma::vec beta,
                       double sigma2_theta_old,
                       arma::mat corr_inv);

double sigma2_theta_update(double a_sigma2_theta,
                           double b_sigma2_theta,
                           arma::vec theta,
                           arma::mat corr_inv);

Rcpp::List phi_update(double a_phi,
                      double b_phi,
                      Rcpp::List temporal_corr_info,
                      arma::vec theta,
                      double sigma2_theta,
                      double phi_old,
                      double metrop_var_phi_trans,
                      int acctot_phi_trans);

double neg_two_loglike_update(arma::vec y,
                              arma::mat x,
                              arma::mat z,
                              arma::vec off_set,
                              int likelihood_indicator,
                              int r,
                              double sigma2_epsilon,
                              arma::vec beta,
                              arma::vec theta);

Rcpp::List GPCW(int mcmc_samples,
                arma::vec y,
                arma::mat x,
                arma::mat z,
                double metrop_var_phi_trans,
                int likelihood_indicator,
                Rcpp::Nullable<Rcpp::NumericVector> offset,
                Rcpp::Nullable<Rcpp::NumericVector> trials,
                Rcpp::Nullable<double> a_r_prior,
                Rcpp::Nullable<double> b_r_prior,
                Rcpp::Nullable<double> a_sigma2_epsilon_prior,
                Rcpp::Nullable<double> b_sigma2_epsilon_prior,
                Rcpp::Nullable<double> sigma2_beta_prior,
                Rcpp::Nullable<double> a_sigma2_theta_prior,
                Rcpp::Nullable<double> b_sigma2_theta_prior,
                Rcpp::Nullable<double> a_phi_prior,
                Rcpp::Nullable<double> b_phi_prior,
                Rcpp::Nullable<double> r_init,
                Rcpp::Nullable<double> sigma2_epsilon_init,
                Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                Rcpp::Nullable<Rcpp::NumericVector> theta_init,
                Rcpp::Nullable<double> sigma2_theta_init,
                Rcpp::Nullable<double> phi_init); 

#endif // __GPCW__
