#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

// split string in to vector of strings
void split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

// [[Rcpp::export]]
Rcpp::List read_ibdData(Rcpp::CharacterVector file, bool header, char sep, 
                         char quote, bool rowNames, int nPos, 
                         Rcpp::CharacterVector na_str, int nrows) {
    string oneLine;
    vector<string> lineElements;
    const string fname = Rcpp::as<string>(file);
    const string naStr = Rcpp::as<string>(na_str);
    int dataStart = 0;
    int posStart = 0;
    if (rowNames) {
        ++dataStart;
        ++posStart;
    }
    dataStart += nPos;
    // Input file stream
    ifstream fileIn(fname.c_str());
    set<char> allelesSet;
    int inx = 0;
    while(inx < 20 && getline(fileIn, oneLine)) {
        lineElements.clear();
        // check if sep is part of first line
        if ((inx == 0) && (oneLine.find(sep) == string::npos)) {
            Rcpp::stop("The separator character 'sep' doesn't exist in line 1");
        }
        if (header && (inx == 0)) {
            ++inx;
            continue;
        }
        // Delete qoutes if exist
        if (quote != ' ') {
            oneLine.erase(remove(oneLine.begin(), oneLine.end(), quote), 
                          oneLine.end());
        }
        // split string in to vector of strings
        split(oneLine, sep, lineElements);
        // detect Nucleotides in this line (SNP)
        for (int i = dataStart; i < lineElements.size(); ++i) {
            if (lineElements[i] != naStr) {
                allelesSet.insert(lineElements[i].begin(), 
                                  lineElements[i].end());
            }
        }
        ++inx;
    }
    vector<char> alleles(allelesSet.begin(), allelesSet.end());
    set<string> labelsSet;
    for (int i = 0; i < alleles.size(); ++i) {
        for (int j = 0; j < alleles.size(); ++j) {
            string comAllele(1, alleles[i]);
            comAllele += alleles[j];
            labelsSet.insert(comAllele);
        }
    }
    // file stream back to start
    fileIn.clear();
    fileIn.seekg(0, ios::beg);
    // If header == true -> first line are individual names
    Rcpp::CharacterVector indivNames;
    if (header) {
        getline(fileIn, oneLine);
        lineElements.clear();
        // Delete qoutes if exist
        if (quote != ' ') {
            oneLine.erase(remove(oneLine.begin(), oneLine.end(), quote ), 
                          oneLine.end());
        }
        // split string in to vector of strings
        split(oneLine, sep, lineElements);
        // Individual names
        indivNames = Rcpp::wrap(lineElements);
    }
    vector<string> lociNames;
    vector<int> pos;
    vector<int> ibdData;
    set<string>::iterator it;
    int parNumber;
    int col = 0;
    int row = 0;
    while(getline(fileIn, oneLine)) {
        lineElements.clear();
        // Delete qoutes if exist
        if (quote != ' ') {
            oneLine.erase(remove(oneLine.begin(), oneLine.end(), quote ), 
                          oneLine.end());
        }
        // split string in to vector of strings
        split(oneLine, sep, lineElements);
        // length of line
        if (col == 0) {
            row = lineElements.size();
        } else if (row != lineElements.size()){
            Rcpp::Rcerr << "Error: Length of line " << col +1 << " is " << 
                lineElements.size() << " instead of " << row <<  "\n" << endl;
            Rcpp::stop("Reading stopt");
        }
        // raw coding of line (IBD) 
        for (int i = dataStart; i < lineElements.size(); ++i) {
            it = labelsSet.find(lineElements[i]);
            if (it == labelsSet.end()) {
                ibdData.push_back(0x00);
            } else {
                parNumber = distance(labelsSet.begin(), it) +1;
                ibdData.push_back(parNumber);
            }
        }
        // SNP names
        if (rowNames) {
            lociNames.push_back(lineElements[0]);
        }
        // genetic positions
        for (int i = 0; i < nPos; ++i) {
            pos.push_back(atoi(lineElements[posStart+i].c_str()));;
        }
        ++col;
        if(nrows > 0 && nrows == col) {
            break;
        }
    }
    // Labels of ibd alleles
    vector<string> labels(labelsSet.begin(), labelsSet.end());
    // Rcpp conversioan and vector as matrix
    Rcpp::IntegerVector position(pos.size());
    if (nPos > 0) {
        if (header) {
            indivNames.erase(indivNames.begin(), indivNames.begin()+nPos);
        }
        position = Rcpp::wrap(pos);
        position.attr("dim") = Rcpp::Dimension(nPos, col);
    }
    Rcpp::RawVector ibdOutData(ibdData.size());
    copy(ibdData.begin(), ibdData.end(), ibdOutData.begin());
    row  = ibdOutData.size() / col;
    ibdOutData.attr("dim") = Rcpp::Dimension(row, col);
    // Results as List
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("indivNames", indivNames),
                                        Rcpp::Named("lociNames", lociNames),
                                        Rcpp::Named("ibdData", ibdOutData),
                                        Rcpp::Named("alleleLabels", labels),
                                        Rcpp::Named("position", position));
    return out;
}
