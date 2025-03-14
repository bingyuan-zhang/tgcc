// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
// [[Rcpp::export]]
List qpp_arma(const arma::mat xmat, const int K, const double r){
  double tmp = 0, tmp_dist = 0;
  int n = xmat.n_rows, d = xmat.n_cols;
  arma::mat tmp_cmat = arma::zeros(K,d);
  arma::vec tmp_dvec = arma::ones(n), index_vec(K), all_index = arma::linspace(0, n-1, n), tmp_index(1);
  for(int k = 0; k < K; ++k){//main iteration
    tmp_index = RcppArmadillo::sample(all_index, 1, false, tmp_dvec);
    index_vec(k) = tmp_index(0);
    tmp_cmat(k,arma::span::all) = xmat(tmp_index(0),arma::span::all);
    for(int i = 0; i < n; ++i){
      tmp = norm( xmat(i,arma::span::all) - tmp_cmat(k,arma::span::all) );
      if(r == 1){
        tmp_dist = tmp;
      }else{
        tmp_dist = pow(tmp,r);
      }
      if(k == 0){
        tmp_dvec(i) = tmp_dist;
      }else if(tmp_dvec(i) > tmp_dist){
        tmp_dvec(i) = tmp_dist;
      }
    }//END FOR i
  }//main iteration
  List res;
  res["centers"] = tmp_cmat;
  res["loss"] = mean(tmp_dvec);
  res["index"] = index_vec;
  return res;
}