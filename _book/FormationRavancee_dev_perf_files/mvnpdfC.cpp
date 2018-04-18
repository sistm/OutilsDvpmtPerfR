#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
const double log2pi2 = log(2.0 * M_PI)/2.0;

//'Based on the implementation from Nino Hardt and Dicko Ahmadou
//'http://gallery.rcpp.org/articles/dmvnorm_arma/
//'(accessed in August 2014)
//'
//'@rdname mvnpdf
//'@export
//'
// [[Rcpp::export]]
NumericVector mvnpdfC(arma::mat x,
                      arma::colvec mean,
                      arma::mat varcovM,
                      bool Log=true){

  int p = x.n_rows;
  int n = x.n_cols;
  NumericVector y(n);

  mat Rinv = inv(trimatu(chol(varcovM)));
  //mat R = chol(as<arma::mat>(varcovM));
  double logSqrtDetvarcovM = sum(log(Rinv.diag()));
  double constant = - p*log2pi2;

  for (int i=0; i < n; i++) {
    colvec x_i = x.col(i) - mean;
    rowvec xRinv = trans(x_i)*Rinv;
    //vec xRinv = solve(trimatl(R.t()), x_i);
    double quadform = sum(xRinv%xRinv);
    if (!Log) {
      y(i) = exp(-0.5*quadform + logSqrtDetvarcovM + constant);
    } else{
      y(i) = -0.5*quadform + logSqrtDetvarcovM + constant;
    }
  }

  return y;

}

