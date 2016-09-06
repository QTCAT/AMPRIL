// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// chiSqTables
List chiSqTables(RawVector x1, RawVector x2, CharacterVector labels);
RcppExport SEXP AMPRIL_chiSqTables(SEXP x1SEXP, SEXP x2SEXP, SEXP labelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< RawVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< RawVector >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type labels(labelsSEXP);
    __result = Rcpp::wrap(chiSqTables(x1, x2, labels));
    return __result;
END_RCPP
}
// read_ibdData
Rcpp::List read_ibdData(Rcpp::CharacterVector file, bool header, char sep, char quote, bool rowNames, int nPos, Rcpp::CharacterVector na_str, int nrows);
RcppExport SEXP AMPRIL_read_ibdData(SEXP fileSEXP, SEXP headerSEXP, SEXP sepSEXP, SEXP quoteSEXP, SEXP rowNamesSEXP, SEXP nPosSEXP, SEXP na_strSEXP, SEXP nrowsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type file(fileSEXP);
    Rcpp::traits::input_parameter< bool >::type header(headerSEXP);
    Rcpp::traits::input_parameter< char >::type sep(sepSEXP);
    Rcpp::traits::input_parameter< char >::type quote(quoteSEXP);
    Rcpp::traits::input_parameter< bool >::type rowNames(rowNamesSEXP);
    Rcpp::traits::input_parameter< int >::type nPos(nPosSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type na_str(na_strSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    __result = Rcpp::wrap(read_ibdData(file, header, sep, quote, rowNames, nPos, na_str, nrows));
    return __result;
END_RCPP
}
// similarityIBD
NumericMatrix similarityIBD(RawMatrix x, NumericVector pWeight, CharacterVector labels);
RcppExport SEXP AMPRIL_similarityIBD(SEXP xSEXP, SEXP pWeightSEXP, SEXP labelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pWeight(pWeightSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type labels(labelsSEXP);
    __result = Rcpp::wrap(similarityIBD(x, pWeight, labels));
    return __result;
END_RCPP
}
