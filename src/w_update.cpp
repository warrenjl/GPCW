#include "RcppArmadillo.h"
#include "GPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w_update(arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec beta_old,
                    arma::vec theta_old){

arma::vec mean_w = x*beta_old + 
                   z*theta_old;

arma::vec w = rcpp_pgdraw(1.00,
                          mean_w);

arma::vec gamma = (y - 0.50)/w;

return Rcpp::List::create(Rcpp::Named("w") = w,
                          Rcpp::Named("gamma") = gamma);

}
































































