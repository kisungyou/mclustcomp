// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

/*
 * Aux 1. compute confusion matrix
 */
//' Compute confusion matrix
//'
//' @keywords internal
// [[Rcpp::export("get.confusion")]]
NumericMatrix getconfusion(NumericVector x, NumericVector y, NumericVector ux, NumericVector uy){
  // 1. preprocessing
  const int k = ux.length();
  const int l = uy.length();
  const int n = x.length();
  NumericMatrix confmat(k,l);

  // 2. main iteration
  for (int it1=0;it1<k;it1++){
    int tgtx = ux[it1];
    for (int it2=0;it2<l;it2++){
      int tgty = uy[it2];
      for (int it3=0;it3<n;it3++){
        if ((x[it3]==tgtx)&&(y[it3]==tgty)){
          confmat(it1,it2) += 1;
        }
      }
    }
  }

  // 3. return output
  return confmat;
}
/*
 * Aux 2. compute community size
 */
//' Compute community size of a clustering
//'
//' @keywords internal
// [[Rcpp::export("get.commsize")]]
NumericVector getcommsize(NumericVector x, NumericVector ux){
  // 1. preliminary
  const int n = x.length();
  const int lux = ux.length();
  NumericVector sizeinfo(lux);
  // 2. main iteration
  for (int i=0;i<lux;i++){
    int tgt = ux[i];
    for (int j=0;j<n;j++){
      if (x[j]==tgt){
        sizeinfo[i]+=1;
      }
    }
  }
  // 3. return
  return sizeinfo;
}
/*
 * Aux 3. compute comembership
 */
//' Comembership matrix of size \code{(2-by-2)}
//'
//' @keywords internal
// [[Rcpp::export("get.pair")]]
NumericMatrix getpair(NumericVector x, NumericVector y){
  // 0. length of vectors
  const int n = x.length();
  // 1. prepare an empty confusion matrix
  NumericMatrix pairmat(2,2);
  // 2. iterate over nC2 case
  for (int i=0;i<n;i++){
    for (int j=(i+1);j<n;j++){
      bool flagx;
      bool flagy;
      // 2-1. check on x
      if (x[i]==x[j]){
        flagx = true;
      } else{
        flagx = false;
      }
      // 2-2. check on y
      if (y[i]==y[j]){
        flagy = true;
      } else{
        flagy = false;
      }
      // 2-3. add on pairmat
      if (flagx){
        if (flagy){
          pairmat(0,0) += 1;
        } else {
          pairmat(0,1) += 1;
        }
      } else {
        if (flagy){
          pairmat(1,0) += 1;
        } else {
          pairmat(1,1) += 1;
        }
      }
    }
  }
  // 3. return
  return pairmat;
}

/*
 * Aux 4. compute probability related materials
 */
//' Compute confusion matrix
//'
//' @keywords internal
// [[Rcpp::export("get.probs")]]
List getprobs(NumericMatrix confmat, NumericVector scx, NumericVector scy, const int n, double threps){
  // 1. preliminary
  const int nk = scx.length();
  const int nl = scy.length();
  NumericVector altthr(2);
  altthr[0] = threps;
  altthr[1] = 1e-7;
  double maxthr = max(altthr);
  int warningint = 0; // 1<-Px, 2<-Py, 3<-Pxy

  // 2. compute::basics
  NumericVector Px  = scx/n;
  NumericVector Py  = scy/n;
  NumericMatrix Pxy = confmat/n;

  // 3. compute::3 measures
  double valx, valy, valxy;
  double log2 = log(2.0);
  double Hx = 0;
  for (int i=0;i<nk;i++){
    valx = Px[i];
    if (valx>0){
      Hx -= valx*log(valx)/log2;
    }
  }
  double Hy = 0;
  for (int j=0;j<nl;j++){
    valy = Py[j];
    if (valy>0){
      Hy -= valy*log(valy)/log2;
    }
  }
  double Ixy = 0;
  for (int i=0;i<nk;i++){
    valx = Px[i];
    for (int j=0;j<nl;j++){
      valy = Py[j];
      valxy = Pxy(i,j);
      if ((valx>0)&&(valy>0)&&(valxy>0)){
        Ixy += valxy*((log(valxy)-log(valx)-log(valy))/log2);
      }
    }
  }

  // 4. return
  List output;
  output["Hx"] = Hx;
  output["Hy"] = Hy;
  output["Ixy"] = Ixy;
  return output;
}
//
//
// List getprobs_old(NumericMatrix confmat, NumericVector scx, NumericVector scy, const int n, double threps){
//   // 1. preliminary
//   const int nk = scx.length();
//   const int nl = scy.length();
//   NumericVector altthr(2);
//   altthr[0] = threps;
//   altthr[1] = 1e-7;
//   double maxthr = max(altthr);
//   int warningint = 0; // 1<-Px, 2<-Py, 3<-Pxy
//
//   // 2. compute::basics
//   NumericVector Px  = scx/n;
//   NumericVector Py  = scy/n;
//   NumericMatrix Pxy = confmat/n;
//   for (int i=0;i<nk;i++){
//     if (Px[i]<threps){
//       Px[i] = maxthr;
//       warningint = 1;
//     } else if (Px[i]>1-threps){
//       Px[i] = 1-maxthr;
//       warningint = 1;
//     }
//   }
//   for (int j=0;j<nl;j++){
//     if (Py[j]<threps){
//       Py[j] = maxthr;
//       warningint = 2;
//     } else if (Py[j]>1-threps){
//       Py[j] = 1-maxthr;
//       warningint = 2;
//     }
//   }
//   for (int i=0;i<nk;i++){
//     for (int j=0;j<nl;j++){
//       if (Pxy(i,j)<threps){
//         Pxy(i,j) = maxthr;
//         warningint = 3;
//       } else if (Pxy(i,j)>1-threps){
//         Pxy(i,j) = 1-maxthr;
//         warningint = 3;
//       }
//     }
//   }
//
//   // 3. compute::3 measures
//   double log2 = log(2.0);
//   double Hx = 0;
//   for (int i=0;i<nk;i++){
//     Hx -= Px[i]*log(Px[i])/log2;
//   }
//   double Hy = 0;
//   for (int j=0;j<nl;j++){
//     Hy -= Py[j]*log(Py[j])/log2;
//   }
//   double Ixy = 0;
//   for (int i=0;i<nk;i++){
//     for (int j=0;j<nl;j++){
//       Ixy += Pxy(i,j)*((log(Pxy(i,j))-log(Px[i])-log(Py[j]))/log2);
//     }
//   }
//
//   // 4. return
//   List output;
//   output["Hx"] = Hx;
//   output["Hy"] = Hy;
//   output["Ixy"] = Ixy;
//   output["warningint"] = warningint;
//   return output;
// }
