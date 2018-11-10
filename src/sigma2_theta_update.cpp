#include "RcppArmadillo.h"
#include "GPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_theta_update(double a_sigma2_theta,
                           double b_sigma2_theta,
                           arma::vec theta,
                           arma::mat corr_inv){

int p_z = theta.size();
double a_sigma2_theta_update = 0.50*p_z + 
                               a_sigma2_theta;

double b_sigma2_theta_update = 0.50*dot(theta, ((corr_inv)*theta)) + 
                               b_sigma2_theta;

double sigma2_theta = 1/R::rgamma(a_sigma2_theta_update,
                                  (1.00/b_sigma2_theta_update));

return(sigma2_theta);

}





