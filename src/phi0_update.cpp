#include "RcppArmadillo.h"
#include "GPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double phi0_update(double a_phi0,
                   double b_phi0,
                   arma::vec theta,
                   arma::mat corr_inv){

int p_z = theta.size();
double a_phi0_update = 0.50*p_z + 
                       a_phi0;

double b_phi0_update = 0.50*dot(theta, ((corr_inv)*theta)) + 
                       b_phi0;

double phi0 = 1/R::rgamma(a_phi0_update,
                          (1/b_phi0_update));

return(phi0);

}





