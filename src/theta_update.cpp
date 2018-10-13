#include "RcppArmadillo.h"
#include "GPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec theta_update(arma::mat x, 
                       arma::mat z,
                       arma::vec w,
                       arma::vec gamma,
                       arma::vec beta,
                       double sigma2_theta_old,
                       arma::mat corr_inv){

int p_z = z.n_cols;
int n = w.size();
arma::mat w_mat(n, p_z);
for(int j = 0; j < p_z; ++j){
   w_mat.col(j) = w;
   }

arma::mat z_trans = trans(z);

arma::mat cov_theta = inv_sympd(z_trans*(w_mat%z) + 
                                (1/sigma2_theta_old)*corr_inv);

arma::vec mean_theta = cov_theta*(z_trans*(w%(gamma - x*beta)));

arma::mat ind_norms = arma::randn(1, p_z);
arma::vec theta = mean_theta + 
                  trans(ind_norms*arma::chol(cov_theta));

return(theta);

}






