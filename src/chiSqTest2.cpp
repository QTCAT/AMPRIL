// #include <Rcpp.h>
// #include <string>
// #include <map>
//   
// using std::map;
// using std::string;
// using std::pair;
// using std::set;
// using std::vector;
// using namespace Rcpp;
// 
// // Estimating allele counts ...
// // [[Rcpp::export]]
// List tables3V(RawVector x1, RawVector x2, RawVector x3, CharacterVector labels) {
//     int labeSize = labels.size();
//     map<string, int> countMap;
//     string lab;
//     set<char> allelesSet;
//     for (int i = 0; i < labels.size(); ++i) {
//       allelesSet.insert(labels[i].begin(), labels[i].end());
//     }
//     vector<char> allele(allelesSet.begin(), allelesSet.end());
//     int dimSize = allele.size();
//     string ABC;
//     for (int i = 0; i < dimSize; ++i) {
//         for (int j = 0; j < dimSize; ++j) {
//             for (int k = 0; k < dimSize; ++k) {
//                 ABC = string() + allele[i] + allele[j] + allele[k];
//                 countMap.insert(pair<string, int>(ABC, 0));
//             }
//         }
//     }
//     string A12;
//     string B12;
//     string C12;
//     string A1B1C1;
//     string A2B2C2;
//     int k = 0;
//     int l = 0;
//     int m = 0;
//     for (int i = 0; i < x1.size(); ++i) {
//         k = x1[i] - 1;
//         l = x2[i] - 1;
//         m = x3[i] - 1;
//         if (k >= 0 && l >= 0 && m >= 0) {
//             A12 = labels[k];
//             B12 = labels[l];
//             C12 = labels[m];
//             A1B1C1 = A12.substr(0, 1) + B12.substr(0, 1) + C12.substr(0, 1);
//             A2B2C2 = A12.substr(1, 1) + B12.substr(1, 1) + C12.substr(1, 1);
//             countMap[A1B1C1] += 1;
//             countMap[A2B2C2] += 1;
//         }
//     }
//     // observation "table"
//     IntegerVector observed = wrap(countMap);
//     // dim frequencies
//     IntegerVector ACounts(dimSize);
//     IntegerVector BCounts(dimSize);
//     IntegerVector CCounts(dimSize);
//     int obsA, obsB, obsC;
//     for (int i = 0; i < dimSize; ++i) {
//         for (int j = 0; j < dimSize; ++j) {
//             for (int k = 0; k < dimSize; ++k) {
//                 obsA = observed[i * dimSize * dimSize + j * dimSize + k];
//                 obsB = observed[i + j * dimSize * dimSize + k * dimSize];
//                 obsC = observed[i + j * dimSize + k * dimSize * dimSize];
//                 ACounts[k] += obsA;
//                 BCounts[k] += obsB;
//                 CCounts[k] += obsC;
//             }
//         }
//     }
//     // number of observation
//     double countAll = Rcpp::sum(ACounts);
//     // expected "table"
//     NumericVector expected(observed.size());
//     double exp;
//     int counter = 0;
//     for (int i = 0; i < dimSize; ++i) {
//         for (int j = 0; j < dimSize; ++j) {
//             for (int k = 0; k < dimSize; ++k) {
//                 exp = (ACounts[i] / countAll) * 
//                       (BCounts[j] / countAll) * 
//                       (CCounts[k] / countAll) * countAll;
//                 expected[counter] = exp;
//                 counter ++;
//             }
//         }
//     }
//     // Result
//     return List::create(observed, expected, ACounts, BCounts, BCounts);
// }
