// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// aster_cpp
NumericVector aster_cpp(NumericMatrix matrix, int tree_width, int proc);
RcppExport SEXP BNSL_aster_cpp(SEXP matrixSEXP, SEXP tree_widthSEXP, SEXP procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< int >::type tree_width(tree_widthSEXP);
    Rcpp::traits::input_parameter< int >::type proc(procSEXP);
    rcpp_result_gen = Rcpp::wrap(aster_cpp(matrix, tree_width, proc));
    return rcpp_result_gen;
END_RCPP
}
// gc
double gc(int n, double a);
RcppExport SEXP BNSL_gc(SEXP nSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(gc(n, a));
    return rcpp_result_gen;
END_RCPP
}
// gc_all
double gc_all(IntegerVector cc, double a);
RcppExport SEXP BNSL_gc_all(SEXP ccSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(gc_all(cc, a));
    return rcpp_result_gen;
END_RCPP
}
// kruskal
IntegerMatrix kruskal(NumericMatrix W);
RcppExport SEXP BNSL_kruskal(SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(kruskal(W));
    return rcpp_result_gen;
END_RCPP
}
// empirical_mi
double empirical_mi(NumericVector x, NumericVector y);
RcppExport SEXP BNSL_empirical_mi(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(empirical_mi(x, y));
    return rcpp_result_gen;
END_RCPP
}
// mi
double mi(NumericVector x, NumericVector y, int proc);
RcppExport SEXP BNSL_mi(SEXP xSEXP, SEXP ySEXP, SEXP procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type proc(procSEXP);
    rcpp_result_gen = Rcpp::wrap(mi(x, y, proc));
    return rcpp_result_gen;
END_RCPP
}
// MDL_mi
double MDL_mi(NumericVector x, NumericVector y, int m_x, int m_y);
RcppExport SEXP BNSL_MDL_mi(SEXP xSEXP, SEXP ySEXP, SEXP m_xSEXP, SEXP m_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type m_x(m_xSEXP);
    Rcpp::traits::input_parameter< int >::type m_y(m_ySEXP);
    rcpp_result_gen = Rcpp::wrap(MDL_mi(x, y, m_x, m_y));
    return rcpp_result_gen;
END_RCPP
}
// Jeffreys_mi
double Jeffreys_mi(NumericVector x, NumericVector y, int m_x, int m_y);
RcppExport SEXP BNSL_Jeffreys_mi(SEXP xSEXP, SEXP ySEXP, SEXP m_xSEXP, SEXP m_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type m_x(m_xSEXP);
    Rcpp::traits::input_parameter< int >::type m_y(m_ySEXP);
    rcpp_result_gen = Rcpp::wrap(Jeffreys_mi(x, y, m_x, m_y));
    return rcpp_result_gen;
END_RCPP
}
// BDeu_mi
double BDeu_mi(NumericVector x, NumericVector y, int m_x, int m_y, double d);
RcppExport SEXP BNSL_BDeu_mi(SEXP xSEXP, SEXP ySEXP, SEXP m_xSEXP, SEXP m_ySEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type m_x(m_xSEXP);
    Rcpp::traits::input_parameter< int >::type m_y(m_ySEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(BDeu_mi(x, y, m_x, m_y, d));
    return rcpp_result_gen;
END_RCPP
}
// empirical_cmi
double empirical_cmi(NumericVector x, NumericVector y, NumericVector z);
RcppExport SEXP BNSL_empirical_cmi(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(empirical_cmi(x, y, z));
    return rcpp_result_gen;
END_RCPP
}
// cmi
double cmi(NumericVector x, NumericVector y, NumericVector z, int proc);
RcppExport SEXP BNSL_cmi(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type proc(procSEXP);
    rcpp_result_gen = Rcpp::wrap(cmi(x, y, z, proc));
    return rcpp_result_gen;
END_RCPP
}
// MDL_cmi
double MDL_cmi(NumericVector x, NumericVector y, NumericVector z, int m_x, int m_y, int m_z);
RcppExport SEXP BNSL_MDL_cmi(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP m_xSEXP, SEXP m_ySEXP, SEXP m_zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type m_x(m_xSEXP);
    Rcpp::traits::input_parameter< int >::type m_y(m_ySEXP);
    Rcpp::traits::input_parameter< int >::type m_z(m_zSEXP);
    rcpp_result_gen = Rcpp::wrap(MDL_cmi(x, y, z, m_x, m_y, m_z));
    return rcpp_result_gen;
END_RCPP
}
// Jeffreys_cmi
double Jeffreys_cmi(NumericVector x, NumericVector y, NumericVector z, int m_x, int m_y, int m_z);
RcppExport SEXP BNSL_Jeffreys_cmi(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP m_xSEXP, SEXP m_ySEXP, SEXP m_zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type m_x(m_xSEXP);
    Rcpp::traits::input_parameter< int >::type m_y(m_ySEXP);
    Rcpp::traits::input_parameter< int >::type m_z(m_zSEXP);
    rcpp_result_gen = Rcpp::wrap(Jeffreys_cmi(x, y, z, m_x, m_y, m_z));
    return rcpp_result_gen;
END_RCPP
}
// BDeu_cmi
double BDeu_cmi(NumericVector x, NumericVector y, NumericVector z, int m_x, int m_y, int m_z, double d);
RcppExport SEXP BNSL_BDeu_cmi(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP m_xSEXP, SEXP m_ySEXP, SEXP m_zSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type m_x(m_xSEXP);
    Rcpp::traits::input_parameter< int >::type m_y(m_ySEXP);
    Rcpp::traits::input_parameter< int >::type m_z(m_zSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(BDeu_cmi(x, y, z, m_x, m_y, m_z, d));
    return rcpp_result_gen;
END_RCPP
}
// mi_matrix
NumericMatrix mi_matrix(DataFrame df, int proc);
RcppExport SEXP BNSL_mi_matrix(SEXP dfSEXP, SEXP procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type proc(procSEXP);
    rcpp_result_gen = Rcpp::wrap(mi_matrix(df, proc));
    return rcpp_result_gen;
END_RCPP
}
// cont_mi
double cont_mi(NumericVector x, NumericVector y);
RcppExport SEXP BNSL_cont_mi(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(cont_mi(x, y));
    return rcpp_result_gen;
END_RCPP
}
// intervals
NumericVector intervals(int level, NumericVector array);
RcppExport SEXP BNSL_intervals(SEXP levelSEXP, SEXP arraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type level(levelSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type array(arraySEXP);
    rcpp_result_gen = Rcpp::wrap(intervals(level, array));
    return rcpp_result_gen;
END_RCPP
}
// binary_search
int binary_search(NumericVector array, double pattern);
RcppExport SEXP BNSL_binary_search(SEXP arraySEXP, SEXP patternSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type array(arraySEXP);
    Rcpp::traits::input_parameter< double >::type pattern(patternSEXP);
    rcpp_result_gen = Rcpp::wrap(binary_search(array, pattern));
    return rcpp_result_gen;
END_RCPP
}
// parent
DataFrame parent(NumericMatrix df0, int h, int tw, int proc);
RcppExport SEXP BNSL_parent(SEXP df0SEXP, SEXP hSEXP, SEXP twSEXP, SEXP procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type df0(df0SEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type tw(twSEXP);
    Rcpp::traits::input_parameter< int >::type proc(procSEXP);
    rcpp_result_gen = Rcpp::wrap(parent(df0, h, tw, proc));
    return rcpp_result_gen;
END_RCPP
}
// fftable
IntegerMatrix fftable(NumericMatrix df, int w);
RcppExport SEXP BNSL_fftable(SEXP dfSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(fftable(df, w));
    return rcpp_result_gen;
END_RCPP
}
// Bayes_score
double Bayes_score(IntegerMatrix T, int m, int proc);
RcppExport SEXP BNSL_Bayes_score(SEXP TSEXP, SEXP mSEXP, SEXP procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type proc(procSEXP);
    rcpp_result_gen = Rcpp::wrap(Bayes_score(T, m, proc));
    return rcpp_result_gen;
END_RCPP
}
// bound
double bound(IntegerMatrix T, int m, int proc);
RcppExport SEXP BNSL_bound(SEXP TSEXP, SEXP mSEXP, SEXP procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type proc(procSEXP);
    rcpp_result_gen = Rcpp::wrap(bound(T, m, proc));
    return rcpp_result_gen;
END_RCPP
}
// Jeffreys_score
double Jeffreys_score(IntegerMatrix T, int m);
RcppExport SEXP BNSL_Jeffreys_score(SEXP TSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(Jeffreys_score(T, m));
    return rcpp_result_gen;
END_RCPP
}
// BDeu_score
double BDeu_score(IntegerMatrix T, int m);
RcppExport SEXP BNSL_BDeu_score(SEXP TSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(BDeu_score(T, m));
    return rcpp_result_gen;
END_RCPP
}
// MDL_score
double MDL_score(IntegerMatrix T, int m);
RcppExport SEXP BNSL_MDL_score(SEXP TSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(MDL_score(T, m));
    return rcpp_result_gen;
END_RCPP
}
// Jeffreys_bound
double Jeffreys_bound(IntegerMatrix T, int m);
RcppExport SEXP BNSL_Jeffreys_bound(SEXP TSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(Jeffreys_bound(T, m));
    return rcpp_result_gen;
END_RCPP
}
// MDL_bound
double MDL_bound(IntegerMatrix T, int m);
RcppExport SEXP BNSL_MDL_bound(SEXP TSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(MDL_bound(T, m));
    return rcpp_result_gen;
END_RCPP
}
// BDeu_bound
double BDeu_bound(IntegerMatrix T, int m);
RcppExport SEXP BNSL_BDeu_bound(SEXP TSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(BDeu_bound(T, m));
    return rcpp_result_gen;
END_RCPP
}