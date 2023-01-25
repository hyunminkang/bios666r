#include <Rcpp.h>
using namespace Rcpp;

//' Calculate IBS matrix efficiently
//' @param genoMat - (n.variants) x (n.sample) matrix of 0,1,2-encoded genotypes
//' @return A matrix containing IBS matrix (symmetric)
// [[Rcpp::export]]
NumericMatrix ibsMatrix(NumericMatrix genoMat) {
    int32_t nvars = genoMat.nrow();
    int32_t nsamps = genoMat.ncol();

    NumericMatrix ibs(nsamps, nsamps);
    for(int i=0; i < nsamps; ++i) {
      for(int j=0; j < i; ++j) {
        int sumDiff = 0;
        for(int k=0; k < nvars; ++k) 
          sumDiff += abs(genoMat(k,i)-genoMat(k,j));
        ibs(i,j) = 1.0 - sumDiff / (2.0 * nvars);
        ibs(j,i) = ibs(i,j);
      }
      ibs(i,i) = 1;
    }
    return ibs;
}

//' Estimate kinship coefficient in homogeneous population
//' @param genoMat - (n.variants) x (n.sample) matrix of 0,1,2-encoded genotypes
//' @return A matrix containing pairwise kinship estimates matrix (symmetric)
// [[Rcpp::export]]
NumericMatrix kingHom(NumericMatrix genoMat) {
    int32_t nvars = genoMat.nrow();
    int32_t nsamps = genoMat.ncol();

    double expH = 0;
    for(int i=0; i < nvars; ++i) {
      int ac = 0;
      for(int j=0; j < nsamps; ++j) {
        ac += (int)genoMat(i,j);
      }
      double af = ac / (2.0 * nsamps);
      expH += (2.0*af*(1.0-af));
    }

    NumericMatrix kin(nsamps, nsamps);
    for(int i=0; i < nsamps; ++i) {
      for(int j=0; j < i; ++j) {
        int sumSqDiff = 0;
        for(int k=0; k < nvars; ++k) 
          sumSqDiff += ((genoMat(k,i)-genoMat(k,j))*(genoMat(k,i)-genoMat(k,j)));
        kin(i,j) = 0.5 - 0.5 * sumSqDiff / expH;
        kin(j,i) = kin(i,j);
      }
      kin(i,i) = 0.5;
    }
    return kin;
}

//' Estimate kinship coefficient in heterogeneous populations
//' @param genoMat - (n.variants) x (n.sample) matrix of 0,1,2-encoded genotypes
//' @return A matrix containing pairwise kinship estimates matrix (symmetric)
// [[Rcpp::export]]
NumericMatrix kingRobust(NumericMatrix genoMat) {
    int32_t nvars = genoMat.nrow();
    int32_t nsamps = genoMat.ncol();

    
    NumericVector nHets(nsamps);
    for(int i=0; i < nsamps; ++i) {
      int h = 0;
      for(int j=0; j < nvars; ++j) {
        if ( genoMat(j,i) == 1 ) ++h;
      }
      nHets[i] = h;
    }

    NumericMatrix kin(nsamps, nsamps);
    for(int i=0; i < nsamps; ++i) {
      for(int j=0; j < i; ++j) {
        int sumSqDiff = 0;
        for(int k=0; k < nvars; ++k) 
          sumSqDiff += ((genoMat(k,i)-genoMat(k,j))*(genoMat(k,i)-genoMat(k,j)));
        kin(i,j) = 0.5 - 0.5 * sumSqDiff / (nHets[j] + nHets[i]);
        kin(j,i) = kin(i,j);
      }
      kin(i,i) = 0.5;
    }
    return kin;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
/*** R
geno_mat = t(matrix(rep(c(0,0,1,1,1,1,1,1,1,1,2,2),2), 6, 4, byrow=TRUE))
ibsMatrix(geno_mat)
kingHom(geno_mat)
kingRobust(geno_mat)
*/
