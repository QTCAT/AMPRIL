#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix similarityIBD(RawMatrix x, NumericVector pWeight, 
                            CharacterVector labels) {
  size_t n = x.nrow();
  std::set<char> allelesSet;
  for (size_t i = 0; i < labels.size(); ++i) {
    allelesSet.insert(labels[i].begin(), labels[i].end());
  }
  std::vector<char> alleles(allelesSet.begin(), allelesSet.end());
  double p = sum( pWeight);
  RawVector rowi, rowj;
  std::string a, b;
  double an, bn, cn;
  NumericMatrix out(n , n);
  for (size_t i = 0; i < (n-1); i ++) {
    rowi = x(i, _);
    for (size_t j = (i+1); j < n; j ++) {
      rowj = x(j, _);
      for (size_t k = 0; k < x.ncol(); k ++) {
        a = labels[rowi[k]-1];
        b = labels[rowj[k]-1];
        for (size_t l = 0; l < alleles.size(); l ++) {
          an = std::count(a.begin(), a.end(), alleles[l]);
          bn = std::count(b.begin(), b.end(), alleles[l]);
          cn += std::min(an, bn)* pWeight[k];
        }
      }
      out(i, j) = cn/(p*2);
      out(j, i) = cn/(p*2);
      cn = 0;
    }
  }
  return out;
}
