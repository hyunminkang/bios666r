#include <Rcpp.h>
using namespace Rcpp;

#include <map>
#include <climits>
#include <utility>
#include <cmath>

// internal function to count unique haplotypes
//   each row in genoMat consists of 0,1,2-encoded genotype
//   this genotype vector is converted into uint64_t key (max of 40 variants)
//   and the count of all unique genotypes are maintained in geno2cnt  
void count_genotypes(NumericMatrix& genoMat, 
                     std::map<uint64_t,int32_t>& geno2cnt) {
  int32_t nrow = genoMat.nrow();
  int32_t ncol = genoMat.ncol();
  geno2cnt.clear(); // clear the genotype counter

  for(int32_t i=0; i < nrow; ++i) {
    // convert genotype vector into a 64-bit integer
    uint64_t key = 0;
    for(int32_t j=0; j < ncol; ++j) 
      key = key * 3 + (uint64_t)genoMat(i,j);
    // count each unique genotype vector
    geno2cnt[key]++;
  }  
}

// type to define a list of haplotype pairs
typedef std::vector<std::pair<int32_t,int32_t> > v_hap_pairs_t;

// internal function to map genotype vector with haplotype pairs
//   construct all possible haplotype pairs and associate them with genotypes
//   possible mapping is maintained in geno2pairs  
void geno_to_hap_pairs(NumericMatrix& possibleHaps, 
                       std::map<uint64_t,v_hap_pairs_t>& geno2pairs) {
  int32_t nhaps = possibleHaps.nrow();
  int32_t nvars = possibleHaps.ncol();
  
  // consider all possible haplotype pairs
  for(int32_t i=0; i < nhaps; ++i) {
    for(int32_t j=0; j <= i; ++j) {
      // construct genotype vector corresponding to the haplotype pair
      uint64_t key = 0;
      for(int32_t k=0; k < nvars; ++k) {
        int32_t h1 = (int32_t)possibleHaps(i, k);
        int32_t h2 = (int32_t)possibleHaps(j, k);
        key = key * 3 + (uint64_t)(h1+h2);
      }
      // associate the genotype vector with possible haplotype pairs
      geno2pairs[key].push_back(std::make_pair(i,j));
    }
  }
}

//' E-M haplotyping with a limited number of markers
//' @param genoMat - (n.samples) x (n.variants) matrix of 0,1,2-encoded genotypes
//' @param possibleHaps - (n.haplotypes) x (n.variants) matrix of 0,1-encoded possible haplotypes
//' @param max_iter - maximum number of E-M iterations
//' @param tol - tolerance of relative log-likelihood changes
//' @return A list containing
//'    * hap_freqs : estimated frequencies of possible haplotypes
//'    * records : matrix keeping track of iteration, log-likelihood, and haplotype frequencies
// [[Rcpp::export]]
List emHaplotyping(NumericMatrix genoMat, NumericMatrix possibleHaps, 
                   int32_t max_iter = 100, double tol = 1e-6) {
    std::map<uint64_t,v_hap_pairs_t> geno2pairs;
    std::map<uint64_t,int32_t> geno2cnt;
    count_genotypes(genoMat, geno2cnt);
    geno_to_hap_pairs(possibleHaps, geno2pairs);
    
    int32_t nsamps = genoMat.nrow();
    int32_t nhaps = possibleHaps.nrow();
    
    std::vector<double> hap_freqs(nhaps);
    std::vector<double> hap_counts(nhaps);
    
    // start with equal haplotype frequencies
    std::fill(hap_freqs.begin(), hap_freqs.end(), 1.0/nhaps);
    
    NumericMatrix records(max_iter, nhaps+2);

    double prev_loglik = -DBL_MAX;    
    int32_t iter;
    for(iter=0; iter < max_iter; ++iter) {
      // E-step, estimate haplotype counts
      std::fill(hap_counts.begin(), hap_counts.end(), 0);
      double loglik = 0;
      for(std::map<uint64_t,int32_t>::iterator it = geno2cnt.begin();
          it != geno2cnt.end(); ++it) {
        uint64_t geno = it->first;
        int32_t cnt = it->second;
        std::map<uint64_t,v_hap_pairs_t>::iterator jt = geno2pairs.find(geno);
        if ( jt == geno2pairs.end() ) {
          stop("Possible haplotypes do not cover all observed genotypes");
        }
        v_hap_pairs_t& pairs = jt->second;
        
        // calculate probabilities of possible haplotype pairs
        std::vector<double> prob_pairs(pairs.size());
        double sum_prob_pairs = 0;
        for(int32_t i=0; i < (int32_t)pairs.size(); ++i) {
          prob_pairs[i] = hap_freqs[pairs[i].first] * hap_freqs[pairs[i].second];
          sum_prob_pairs += prob_pairs[i];
        }
        
        loglik += (cnt * log(sum_prob_pairs));
        
        // normalize the probabilities
        for(int32_t i=0; i < (int32_t)pairs.size(); ++i) {
          prob_pairs[i] /= sum_prob_pairs;
          hap_counts[pairs[i].first]  += (double)(cnt * prob_pairs[i]);
          hap_counts[pairs[i].second] += (double)(cnt * prob_pairs[i]);
        }
      }
      
      // store the intermediate results
      records(iter,0) = iter;
      records(iter,1) = loglik;
      for(int32_t i=0; i < nhaps; ++i) 
        records(iter,2+i) = hap_freqs[i];
 
      // check the convergence
      if ( (loglik - prev_loglik)/fabs(loglik) < tol ) {
        break;
      }
      prev_loglik = loglik;
      
      // M-step, update the haplotype frequencies
      for(int32_t i=0; i < nhaps; ++i) {
        hap_freqs[i] = hap_counts[i] / (nsamps * 2);
      }
    }
    
    // return haplotype frequencies and records at individual iterations
    return List::create(Named("hap_freqs") = hap_freqs,
                        Named("records") = records(Range(0,iter),_));
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
/*** R
geno_mat = matrix(rep(c(0,0,1,1,2,2),10), 30, 2, byrow=TRUE)
emHaplotyping(geno_mat, matrix(c(0,0,0,1,1,0,1,1),4,2,byrow=TRUE))
*/
