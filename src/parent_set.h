#ifndef PARENT_SET_H
#define PARENT_SET_H

#include <Rcpp.h>

Rcpp::DataFrame parent(Rcpp::NumericMatrix, int, int, int);
Rcpp::IntegerMatrix fftable(Rcpp::NumericMatrix, int);
double Bayes_score(Rcpp::IntegerMatrix, int, int);
double Jeffreys_score(Rcpp::IntegerMatrix, int);
double MDL_score(Rcpp::IntegerMatrix, int);
double BDeu_score(Rcpp::IntegerMatrix, int);
double bound(Rcpp::IntegerMatrix, int, int);
double Jeffreys_bound(Rcpp::IntegerMatrix, int);
double MDL_bound(Rcpp::IntegerMatrix, int);
double BDeu_bound(Rcpp::IntegerMatrix, int);

#endif // PARENT_SET_H
