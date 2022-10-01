#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* CPP ROUTINES RE-WRITTEN IN RCPPARMADILLO
 * 1. getconfusion : compute the confusion matrix
 *
 */


// 1. getconfusion -------------------------------------------------------------
// [[Rcpp::export("get.confusion")]]
arma::mat genconfusion(arma::vec &x, arma::vec &y, arma::vec &ux, arma::vec &uy){
  // parameters
  int k = ux.n_elem;
  int l = uy.n_elem;
  int n = x.n_elem;
  arma::mat confmat(k,l,fill::zeros);

  // main iteration
  int tgtx = 0;
  int tgty = 0;

  for (int it1=0; it1<k; it1++){
    tgtx = ux(it1);
    for (int it2=0; it2<l; it2++){
      tgty = uy(it2);

      // for all elements, add if they match
      for (int it3=0; it3<n; it3++){
        if ((x(it3)==tgtx)&&(y(it3)==tgty)){
          confmat(it1,it2) += 1;
        }
      }

    }
  }

  // return
  return(confmat);
}
