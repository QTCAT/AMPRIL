#include <Rcpp.h>
#include <string>
#include <map>

using std::map;
using std::string;
using std::pair;
using namespace Rcpp;

// Estimating allele counts ...
// [[Rcpp::export]]
List chiSqTables(RawVector x1, RawVector x2, CharacterVector labels) {
    int labeSize = labels.size();
    map<string, int> countMap;
    string lab;
    for (int i = 0; i < labeSize; ++i) {
        lab = labels[i];
        countMap.insert(pair<string, int>(lab, 0));
    }
    string A12;
    string B12;
    string A1B1;
    string A2B2;
    int k = 0;
    int l = 0;
    for (int i = 0; i < x1.size(); ++i) {
        k = x1[i] - 1;
        l = x2[i] - 1;
        if (k >= 0 && l >= 0) {
            A12 = labels[k];
            B12 = labels[l];
            A1B1 = A12.substr(0, 1) + B12.substr(0, 1);
            A2B2 = A12.substr(1, 1) + B12.substr(1, 1);
            countMap[A1B1] += 1;
            countMap[A2B2] += 1;
        }
    }
    // observation table
    IntegerVector observed = wrap(countMap);
    // row and col frequencies
    int dimSize = sqrt(labeSize);
    IntegerVector rowCounts(dimSize);
    IntegerVector colCounts(dimSize);
    for (int i = 0; i < dimSize; ++i) {
        for (int j = 0; j < dimSize; ++j) {
            rowCounts[i] += observed[i * dimSize + j];
            colCounts[i] += observed[i + j * dimSize];
        }
    }
    // number of observation
    double countAll = Rcpp::sum(rowCounts);
    // expected table and chi squre value
    NumericVector expected(observed.size());
    int count = 0;
    for (int i = 0; i < dimSize; ++i) {
        for (int j = 0; j < dimSize; ++j) {
            expected[count] = rowCounts[i] * colCounts[j] / countAll;
            count ++;
        }
    }
    // Result
    return List::create(observed, expected);
}  // chiSqTables
