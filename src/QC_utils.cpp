// Utilities to support qcsim()
// Author: S Master
// masters@chop.edu

#include <Rcpp.h>
using namespace Rcpp;

// Return a vector that contains the max value for each column in the input matrix
// [[Rcpp::export]]
NumericVector colMax(NumericMatrix X) {
  int ncol = X.ncol();
  int nrow = X.nrow();
  Rcpp::NumericVector out(ncol);
  for (int col = 0; col < ncol; col++){
    double newmax = X(0, col);
    for (int row = 0; row < nrow; row++) {
      newmax = (X(row, col) > newmax) ? X(row, col) : newmax;
    }
    out[col]=newmax;
  }
  return wrap(out);
}

// Return a vector that contains the min value for each column in the input matrix
// [[Rcpp::export]]
NumericVector colMin(NumericMatrix X) {
  int ncol = X.ncol();
  int nrow = X.nrow();
  Rcpp::NumericVector out(ncol);
  for (int col = 0; col < ncol; col++){
    double newmin = X(0, col);
    for (int row = 0; row < nrow; row++) {
      newmin = (X(row, col) < newmin) ? X(row, col) : newmin;
    }
    out[col]=newmin;
  }
  return wrap(out);
}

// return a logical vector that indicates if at least contig sequential values in a column exceed cutoff
// [[Rcpp::export]]
LogicalVector cutoffContig(NumericMatrix X, double cutoff, int contig) {
  int ncol = X.ncol();
  int nrow = X.nrow();
  Rcpp::LogicalVector out(ncol);
  for (int col = 0; col < ncol; col++) {
    int streak = 0;
    int pn = 0;
    for (int row = 0; row < nrow; row++) {
      double val = X(row, col);
      if (val > cutoff) {
        if (pn == 1) {
          streak++;
          if (streak == contig) break;
        } else {
          pn = 1;
          streak = 1;
        }
      } else if (std::fabs(val) > cutoff) {
        if (pn == -1) {
          streak++;
          if (streak == contig) break;
        } else {
          pn = -1;
          streak = 1;
        }
      } else pn = 0;
    }
    out[col] = (streak == contig) ? TRUE : FALSE;
  }
  return wrap(out);
}


// return a logical vector that indicates if T test or F test fails in each column
// [[Rcpp::export]]
LogicalVector tftestMat(NumericMatrix X1, NumericMatrix X2,
                        int useT, int useF,
                        double tcutoff,
                        double fcutofflow, double fcutoffhigh) {
  if (X1.ncol() != X2.ncol()) {
    Rcpp::LogicalVector out(1);
    out[0] = NA_LOGICAL;
    return(wrap(out));
  }
  int ncol = X1.ncol();
  int n1 = X1.nrow();
  int n2 = X2.nrow();
  double npool = sqrt((1/(double)n1) + (1/(double)n2));
  Rcpp::LogicalVector out(ncol);
  bool tout = false, fout = false;

  for (int col = 0; col < ncol; col++) {
    double mean1 = 0.0;
    double mean2 = 0.0;
    for (int row = 0; row < n1; row++) {
      mean1 += X1(row, col);
    }
    for (int row = 0; row < n2; row++) {
      mean2 += X2(row, col);
    }
    mean1 /= n1;
    mean2 /= n2;
    double sx1 = 0.0;
    double sx2 = 0.0;
    double x = 0.0;
    for (int row = 0; row < n1; row++) {
      x = X1(row, col);
      sx1 += (x - mean1) * (x - mean1);
    }
    sx1 /= (n1 - 1);
    for (int row = 0; row < n2; row++) {
      x = X2(row, col);
      sx2 += (x - mean2) * (x - mean2);
    }
    sx2 /= (n2 - 1);

    if (useT) {
      double sp = sqrt((((n1 - 1 ) * sx1) + ((n2 - 1) * sx2)) / (n1 + n2 - 2));
      double t = fabs((mean1 - mean2) / (sp * npool));
      tout = (t > tcutoff) ? true : false;
    }

    if (useF) {
      fout = ((sx1/sx2 > fcutoffhigh) || (sx1/sx2 < fcutofflow)) ? true : false;
    }
    out[col] = tout || fout;
  }
  return wrap(out);
}
