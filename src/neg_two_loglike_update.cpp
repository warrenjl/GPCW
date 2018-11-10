#include "RcppArmadillo.h"
#include "GPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double neg_two_loglike_update(arma::vec y,
                              arma::mat x,
                              arma::mat z, 
                              arma::vec beta,
                              arma::vec theta){

int n = y.size();
arma::vec dens(n); dens.fill(0.00);

arma::vec logit_probs = x*beta + 
                        z*theta;

arma::vec probs = exp(logit_probs)/(1.00 + exp(logit_probs));

for(int j = 0; j < n; ++j){
   dens(j) = R::dbinom(y(j),
                       1,
                       probs(j),
                       TRUE);
   }

double neg_two_loglike = -2.00*sum(dens);

return neg_two_loglike;

}

























































