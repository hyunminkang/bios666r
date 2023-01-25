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
        ibs(i,j) = 2.0 - (double) sumDiff / (double) nvars;
        ibs(j,i) = ibs(i,j);
      }
      ibs(i,i) = 1;
    }
    rownames(ibs) = colnames(genoMat);
    colnames(ibs) = colnames(genoMat);
    return ibs;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
/*** R
geno_mat = t(matrix(rep(c(0,0,1,1,1,1,1,1,1,1,2,2),2), 6, 4, byrow=TRUE))
ibsMatrix(geno_mat)
*/
