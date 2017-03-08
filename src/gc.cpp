#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>

#include "gc.h"

// [[Rcpp::export]]
double gc(int n, double a){if(n>0)return(gc(n-1,a))+log(n+a-1); else return(0);}

// [[Rcpp::export]]
double gc_all(IntegerVector cc, double a){
 	int m=cc.size();
 	if(m==1)return(gc((int)cc(0),a)); else return(gc_all(tail(cc,m-1),a)+gc((int)cc(0),a));
}
