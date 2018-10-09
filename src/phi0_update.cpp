#include "RcppArmadillo.h"
#include "GPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double phi0_update(double alpha_phi0,
                   double beta_phi0,
                   arma::vec theta,
                   arma::mat corr_inv){

int p_z = theta.size();
double alpha_phi0_update = 0.50*p_z + 
                           alpha_phi0;

double beta_phi0_update = 0.50*dot(theta, ((corr_inv)*theta)) + 
                          beta_phi0;

double phi0 = 1/R::rgamma(alpha_phi0_update,
                          (1/beta_phi0_update));

return(phi0);

}





