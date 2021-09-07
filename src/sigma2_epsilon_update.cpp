#include "RcppArmadillo.h"
#include "GPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_epsilon_update(arma::vec y,
                             arma::mat x,
                             arma::mat z,
                             double a_sigma2_epsilon,
                             double b_sigma2_epsilon,
                             arma::vec beta_old,
                             arma::vec theta_old){

int n = y.size();
double a_sigma2_epsilon_update = 0.50*n + 
                                 a_sigma2_epsilon;

double b_sigma2_epsilon_update = 0.50*dot((y - x*beta_old - z*theta_old), (y - x*beta_old - z*theta_old)) + 
                                 b_sigma2_epsilon;

double sigma2_epsilon = 1.00/R::rgamma(a_sigma2_epsilon_update,
                                       (1.00/b_sigma2_epsilon_update));

return(sigma2_epsilon);

}





