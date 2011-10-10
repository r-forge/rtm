#include <vector>

#include <Rcpp.h>

class Sampler {
public:
  Sampler(int K, Rcpp::List posteriors, SEXP corpus) :
    K_(K), posteriors_(posteriors), p_(K),
    corpus_(*wrappedObject<Corpus>(corpus)) {
    p_.resize(K);
  }

  void sampleDocument(const Document& document) {
    int state_size = 0;
    for (int ii = 0; ii < state_size; ++ii) {
      double p_sum = 0.0;
      for (int kk = 0; kk < K_; ++kk) {
        p_[kk] = 1.0;
        for (int jj = 0; jj < posteriors_.size(); ++jj) {
          p_[kk] *=
            posteriors_.getMutableData(jj)->probability(document, ii, kk);
        }
        p_sum += p_[kk];
      }
      for (int kk = 0; kk < K_; ++kk) {
        p_[kk] / p_sum;
      }
    }
  }

  void sampleCorpus() {
    for (int ii = 0; ii < corpus_.numDocuments(); ++ii) {
      sampleDocument(corpus_.getDocument(ii));
    }
  }

private:
  int K_;
  WrappedList<Posterior> posteriors_;
  const Corpus& corpus_;
  std::vector<double> p_;
};
